library(lixoftConnectors)
library(mvtnorm)
library(nlme)
library(MASS)
library(saemix)
library(tidyverse)
library(testit)
library(here)
library(docopt)
library(ggplot2)
initializeLixoftConnectors(software="monolix")
setwd(here("Desktop","Analysis"))

######################## main function: covariates on beta1 in rebound model
# filename is a string to indicate the raw data file name. 
# Raw data includes "PATIENT","CD4_V","RNA_V","RNA_L","treatment", and "days_from_seroco"
# num_bootstrap is the number of bootstrap want to repeat
# outliers is a vector of PATIENT ID which are outliers, e.g., outliers = c(1,5,10)
# num_months is a number to indicate the time interval for viral rebound analysis. 
analysis_beta1 <- function(filename, num_bootstrap, outliers=NULL, num_months=NULL){
  ######################### Clean raw data
  # Read raw data
  raw_data<-read.csv(filename,header = T)
  # If RNA_V = -1, indicate the observation as censored and replace RNA_V with RNA_L.
  # Take log10 to RNA_V. Convert time from days to months. 
  # Set lower boundary (limit) for the censoring interval to 0.
  # Remove rows with unobserved RNA_V and time
  raw_data <- raw_data %>% 
    mutate(censor = ifelse(RNA_V==-1, 1, 0)) %>% 
    mutate(RNA_V = ifelse(RNA_V==-1, RNA_L, RNA_V)) %>% 
    mutate(log10 = log10(RNA_V)) %>% 
    mutate(time = days_from_seroco/30) %>% 
    mutate(limit = 0) %>% 
    filter(!is.na(RNA_V) & !is.na(time))
  # Remove outliers if there is any
  if(!is.null(outliers)){
    raw_data <- raw_data %>% 
      filter(PATIENT != outliers)
  }
  # Find lowest detection limit for each patient
  raw_data <- raw_data %>%
    group_by(PATIENT) %>%
    mutate(RNA_L = min(RNA_L,na.rm = TRUE)) 
  # Split the dataset into before_ART (treatment=1) and after_ART(treatment=0 after treatment=1). 
  # Note: for after_ART, ignore all "treatment=0" before "treatment=1".
  # Change time after ART to the time since ART interruption 
  # For example: if before_ART time is c(1,2,3), after_ART time is c(3.5,4,5), then 
  #             the time since ART interruption for after_ART is c(0.5,1,2).
  before_ART <- raw_data %>% 
    filter(treatment==1) %>% 
    dplyr::select("PATIENT","RNA_L","CD4_V","censor","log10","time","limit")
  after_ART <- raw_data %>% 
    group_by(PATIENT) %>% 
    mutate(ind = case_when(treatment == 1 ~ time)) %>%  
    mutate(ind = max(ind, na.rm=TRUE)) %>% 
    filter(time > ind) %>% 
    mutate(time = time-ind) %>% 
    dplyr::select("PATIENT","RNA_L","censor","log10","time","limit")
  # Only leave first num_months months data for after_ART, only leave useful columns
  if(!is.null(num_months)){
    after_ART <- after_ART %>% 
      group_by(PATIENT) %>% 
      mutate(min = min(time)) %>% 
      filter(time <= num_months+min) %>% 
      dplyr::select("PATIENT","RNA_L","censor","log10","time","limit")
  }
  
  # Save the two datasets locally. This is required when run Monolix.
  write.csv(before_ART,"before_ART.csv",row.names = FALSE)
  write.csv(after_ART,"after_ART.csv",row.names = FALSE)
  
  ######################## Find population parameter estimates before ART interruption in Monolix
  # model1: y_ij = log10(exp(P_1i-lambda_1i*t_ij)+exp(P_2i-lambda_2i*t_ij))+e_ij
  # P_1i = P_1 + b_1i, P_2i = P_2 + b_2i, lambda_1i = lambda_1+b_3i, lambda_2i = lambda_2+b_4i
  # b_i ~ N(0,B), e_ij ~ N(0, sigma^2)
  data = list(dataFile = paste0('before_ART.csv'),
              headerTypes =c("id","ignore","ignore","cens","observation","time","limit"),
              observationTypes = list(RE="continuous"))
  modelFile = paste0('before_ART_model.txt')
  # create a new project by setting a data set and a structural model
  newProject(data = data, modelFile =modelFile)
  setErrorModel(log10_ = "constant")
  # set tasks in scenario
  scenario <- getScenario()
  scenario$tasks = c(populationParameterEstimation = T,
                     conditionalModeEstimation = T,
                     conditionalDistributionSampling = T,
                     standardErrorEstimation=T,
                     logLikelihoodEstimation=T)
  scenario$linearization = FALSE
  setScenario(scenario)
  setIndividualParameterDistribution(b1="normal",p1="normal",p2="normal",b2="normal")
  # run the estimation
  runScenario()
  # store the estimates
  popestimates_before <- getEstimatedPopulationParameters()
  # store the random effects 
  bi <- as.data.frame(getEstimatedRandomEffects()$conditionalMean)
  # store the s.e.
  se_before <- unlist(getEstimatedStandardErrors())
  # if P1 < P2, switch P1,P2 and lambda1,lambda2
  popestimates_tmp <- popestimates_before
  bi_tmp <- bi
  se_tmp <- se_before
  if (popestimates_before[1]<popestimates_before[3]){
    popestimates_before[1]<-popestimates_tmp[3]
    popestimates_before[3]<-popestimates_tmp[1]
    popestimates_before[2]<-popestimates_tmp[4]
    popestimates_before[4]<-popestimates_tmp[2]
    bi[,4]<-bi_tmp[,2]
    bi[,2]<-bi_tmp[,4]
    bi[,3]<-bi_tmp[,5]
    bi[,5]<-bi_tmp[,3]
    se_before[1]<-se_tmp[3]
    se_before[3]<-se_tmp[1]
    se_before[2]<-se_tmp[4]
    se_before[4]<-se_tmp[2]
  }
  # Find population parameter estimates
  P1 <- popestimates_before[1]
  lambda1 <- popestimates_before[2]
  P2 <- popestimates_before[3]
  lambda2 <- popestimates_before[4]
  # Find sigma
  sigma <- popestimates_before[9]
  # Find random effects
  colnames(bi)<-c("PATIENT","b1_est","b3_est","b2_est","b4_est")
  # Find B
  B<-as.matrix(cov(bi[2:5]))
  
  
  ######################## Find population parameter estimates CD4
  # Read before_ART data
  before_ART <- read.csv("before_ART.csv")
  before_ART$CD4_V<-as.numeric(before_ART$CD4_V)
  long.data<-groupedData(CD4_V~time|PATIENT,data=before_ART)
  long.data$time<-as.numeric(long.data$time)
  long.data$PATIENT<-as.factor(long.data$PATIENT)
  long.data <- long.data %>% 
    filter(!is.na(CD4_V))
  # model2: sqrt(z_ij) = alpha_1i + alpha_2i*t_ij + epsilon_ij
  # alpha_1i = alpha_1 + a_1i, alpha_2i = alpha_2 + a_2i
  # a_i ~ N(0, A), epsilon_ij ~ N(0, delta^2)
  lme <- lme(sqrt(CD4_V)~time, data=long.data,
             random = ~time|PATIENT,method="ML")
  # Find population parameter estimates
  alpha1 <- summary(lme)$coefficients$fixed[1]
  alpha2 <- summary(lme)$coefficients$fixed[2]
  # Find random effects
  ai <- ranef(lme)
  ai <- rownames_to_column(ai, "PATIENT")
  colnames(ai)=c("PATIENT","a1_est","a2_est")
  # Find A
  A <- as.matrix(cov(ai[,2:3]))
  # Find delta
  delta <- summary(lme)$sigma
  # Find s.e.
  se_cd4 <- coef(summary(lme))[,2]
  
  
  ######################## Find population parameter estimates after ART interruption
  # Read after ART data
  after_ART <- read.csv("after_ART.csv")
  # Add random effects from model1 and model2 to after ART data
  after_ART <- merge(after_ART, bi, by="PATIENT")
  after_ART <- merge(after_ART, ai, by="PATIENT")
  # Save after ART data with covariates
  write.csv(after_ART,"after_ART_cov.csv",row.names = FALSE)
  # w_ij = beta_1i + t_ij/(t_ij + exp(beta_2i - beta_3i * t_ij)) + beta_4i + xi_ij
  # beta_1i = beta_1+gamma_1*b_1i+gamma_2*b_2i+gamma_3*b_3i+gamma_4*b_4i+gamma_5*a_1i+gamma_6*a_2i+tau_1i
  # beta_ji = beta_j + tau_ji, j=2,3,4
  # tau_i ~ N(0, G), xi_ij ~ N(0, omega^2)
  data = list(dataFile = paste0('after_ART_cov.csv'),
              headerTypes =c("id","ignore","cens","observation","time","limit","contcov","contcov","contcov","contcov","contcov","contcov"),
              observationTypes = list(RE="continuous"))
  modelFile = paste0('after_ART_model.txt')
  newProject(data = data, modelFile =modelFile)
  setIndividualParameterDistribution(beta1="normal",beta2="normal",beta3="normal",beta4="normal")
  # covariates on beta1
  setCovariateModel(beta1 = c(a1_est=TRUE,a2_est=TRUE,b1_est=TRUE,b2_est=TRUE,b3_est=TRUE,b4_est = TRUE))
  setErrorModel(log10_ = "constant")
  # run the estimation
  scenario <- getScenario()
  scenario$tasks = c(populationParameterEstimation = T,
                     conditionalModeEstimation = T,
                     conditionalDistributionSampling = T,
                     standardErrorEstimation=T,
                     logLikelihoodEstimation=T)
  scenario$linearization = FALSE
  setScenario(scenario)
  runScenario()
  # Find population parameter estimates
  popestimates_after <- getEstimatedPopulationParameters()
  beta1 <- popestimates_after[1]
  gamma_5 <- popestimates_after[2]
  gamma_6 <- popestimates_after[3]
  gamma_1 <- popestimates_after[4]
  gamma_2 <- popestimates_after[5]
  gamma_3 <- popestimates_after[6]
  gamma_4 <- popestimates_after[7]
  beta2 <- popestimates_after[8]
  beta3 <- popestimates_after[9]
  beta4 <- popestimates_after[10]
  # Find omega
  omega <- popestimates_after[15]
  # Find G
  taui <- as.data.frame(getEstimatedRandomEffects()$conditionalMean)
  G <- as.matrix(cov(taui[,2:5]))
  # Find s.e.
  se_after <- unlist(getEstimatedStandardErrors())
  
  
  ######################## Bootstrap
  # Create a empty dataset to store all bootstrap values
  names<-c("P1","b1","P2","b2","alpha1","alpha2","beta1","g5","g6","g1","g2","g3","g4","beta2","beta3","beta4")
  bootstrap<-as.data.frame(matrix(NA,nrow=length(names),ncol=num_bootstrap))
  rownames(bootstrap)<-names
  # Create a summary table, first column is population estimates and second column is original s.e.
  popestimates<-c(popestimates_before[1:4],alpha1,alpha2,popestimates_after[1:10])
  se <- c(se_before[1:4], se_cd4, se_after[1:10])
  se[which(se=="NaN")]<-NA
  se <- as.numeric(se)
  summary <- cbind(popestimates,se)
  colnames(summary)=c("Estimate","Original SE")
  rownames(summary)<-names
  # Find unique patient and number of patients
  n_before<-length(unique(before_ART$PATIENT))
  PATIENT_before<-unique(before_ART$PATIENT)
  n_after<-length(unique(after_ART$PATIENT))
  PATIENT_after<-unique(after_ART$PATIENT)
  # repeat bootstrap num_bootstrap times.
  for (i in 1:num_bootstrap){
    print(glue::glue("i = ",i))
    ############# Before ART 
    # simulate b_i
    b_sim<-rmvnorm(n_before, c(0,0,0,0), B)
    colnames(b_sim)=c("b1_sim","b3_sim","b2_sim","b4_sim")
    b_sim<-cbind(PATIENT_before,b_sim)
    colnames(b_sim)[1] <- "PATIENT"
    # add simulated b_i to before_ART.
    before_ART<-merge(before_ART,b_sim,by="PATIENT")
    # simulate e and y_ij, if y_ij < detection limit, then the y_ij = detection limit
    before_ART <- before_ART %>% 
      mutate(e = rnorm(nrow(before_ART),0,sigma)) %>% 
      mutate(simulated = log10(exp(P1+b1_sim-(P2+b3_sim)*time) + exp ((lambda1+b2_sim) -(lambda2+b4_sim)*time))+e) %>% 
      mutate(simulated = ifelse(simulated<log10(RNA_L), log10(RNA_L), simulated))
    # simulate a_i
    a_sim<-rmvnorm(n_before, c(0,0), A)
    colnames(a_sim)=c("a1_sim","a2_sim")
    a_sim<-cbind(PATIENT_before,a_sim)
    colnames(a_sim)[1] <- "PATIENT"
    # add simulated a_i to before_ART.
    before_ART<-merge(before_ART,a_sim,by="PATIENT")
    # simulate epsilon and z_ij, if z_ij < 1, then z_ij = 1
    before_ART <- before_ART %>% 
      mutate(epsilon = rnorm(nrow(before_ART),0, delta)) %>% 
      mutate(simulatedcd4 = alpha1+a1_sim+time*(alpha2+a2_sim)+epsilon) %>% 
      mutate(simulatedcd4 = ifelse(simulatedcd4<1, 1, simulatedcd4))
    # leave useful columns
    before_ART <- before_ART %>% 
      dplyr::select(PATIENT,RNA_L,censor,time,limit,simulated,simulatedcd4)
    # save before ART simulated data
    write.table(before_ART,"before_ART_sim.txt", sep = "," ,quote = FALSE,row.names = FALSE)
    # Find new estimates of P1,P2,lambda1, and lambda2 for simulated viral loads
    data = list(dataFile = paste0('before_ART_sim.txt'),
                headerTypes =c("id","ignore","cens","time","limit", "observation","ignore"),
                observationTypes = list(RE="continuous"))
    modelFile = paste0('before_ART_model.txt')
    newProject(data = data, modelFile =modelFile)
    setErrorModel(simulated = "constant")
    # set tasks in scenario
    scenario <- getScenario()
    scenario$tasks = c(populationParameterEstimation = T,
                       conditionalModeEstimation = T,
                       conditionalDistributionSampling = T,
                       standardErrorEstimation=T,
                       logLikelihoodEstimation=T)
    scenario$linearization = FALSE
    setScenario(scenario)
    setIndividualParameterDistribution(b1="normal",p1="normal",p2="normal",b2="normal")
    # run the estimation
    runScenario()
    # find population parameter estimates
    popestimates_before <- getEstimatedPopulationParameters()
    # find random effects
    bi <- as.data.frame(getEstimatedRandomEffects()$conditionalMean)
    # if P1 < P2, switch P1,P2 and lambda1,lambda2
    popestimates_tmp <- popestimates_before
    bi_tmp <- bi
    if (popestimates_before[1]<popestimates_before[3]){
      popestimates_before[1]<-popestimates_tmp[3]
      popestimates_before[3]<-popestimates_tmp[1]
      popestimates_before[2]<-popestimates_tmp[4]
      popestimates_before[4]<-popestimates_tmp[2]
      bi[,4]<-bi_tmp[,2]
      bi[,2]<-bi_tmp[,4]
      bi[,3]<-bi_tmp[,5]
      bi[,5]<-bi_tmp[,3]
    }
    # store the estimates in bootstrap dataframe
    bootstrap[1:4,i]=popestimates_before[1:4]
    # find the random effects bi
    colnames(bi)<-c("PATIENT","b1_est","b3_est","b2_est","b4_est")
    
    ############# CD4 Data
    long.data<-groupedData(simulatedcd4~time|PATIENT,data=before_ART)
    long.data$time<-as.numeric(long.data$time)
    long.data$PATIENT<-as.factor(long.data$PATIENT)
    # If the lme does not converge, then re-simulate cd4 data
    while (has_error (lme <- lme(simulatedcd4~time, data=long.data, 
                                 random = ~time|PATIENT,method="ML"))){
      a_sim<-rmvnorm(n_after, c(0,0), A)
      colnames(a_sim)=c("a1_sim","a2_sim")
      a_sim<-cbind(PATIENT_before,a_sim)
      colnames(a_sim)[1] <- "PATIENT"
      before_ART<-merge(before_ART,a_sim,by="PATIENT")
      before_ART <- before_ART %>% 
        mutate(epsilon = rnorm(nrow(before_ART),0,delta)) %>% 
        mutate(simulatedcd4 = alpha1+a1_sim+time*(alpha2+a2_sim)+epsilon) %>% 
        mutate(simulatedcd4 = ifelse(simulatedcd4<1, 1, simulatedcd4))
      before_ART <- before_ART %>% 
        dplyr::select(PATIENT,RNA_L,censor,time,limit,simulated,simulatedcd4)
      long.data<-groupedData(simulatedcd4~time|PATIENT,data=before_ART)
      long.data$time<-as.numeric(long.data$time)
      long.data$PATIENT<-as.factor(long.data$PATIENT)
    }
    # store the new estimates in bootstrap dataframe
    bootstrap[5:6,i]<-summary(lme)$coefficients$fixed
    # find the random effects ai
    ai<-ranef(lme)
    ai <- rownames_to_column(ai,"PATIENT")
    colnames(ai)[2:3]<-c("a1_est","a2_est")
    
    ############# Bootstrap before ART 
    # leave useful columns
    after_ART <- after_ART %>% 
      dplyr::select("PATIENT","RNA_L","censor","time","limit")
    # simulate tau_i
    taui<-rmvnorm(n_after, c(0,0,0,0), G)
    colnames(taui)=c("tau1","tau2","tau3","tau4")
    taui<-cbind(PATIENT_after,taui)
    colnames(taui)[1]<-"PATIENT"
    # Add random effects from previous two models and simulated taui to after ART data
    after_ART<-merge(after_ART,taui,by="PATIENT")
    after_ART <- merge(after_ART,bi,by="PATIENT")
    after_ART<-merge(after_ART,ai,by="PATIENT")
    # simulate xi_ij 
    after_ART <- after_ART %>% 
      mutate(xi = rnorm(nrow(after_ART),0,omega))
    # simulate w_ij, if w_ij < detection limit, then w_ij = detection limit.
    after_ART <- after_ART %>% 
      mutate(simulated = (beta1+a1_est*gamma_5+a2_est*gamma_6+b1_est*gamma_1+b2_est*gamma_2+b3_est*gamma_3+b4_est*gamma_4+tau1) * time / ( time + exp( beta2+tau2 - (beta3+tau3) * time )) + (beta4+tau4)+xi)
    # if w_ij < detection limit, then w_ij = detection limit.
    after_ART <- after_ART %>% 
      mutate(simulated = ifelse(simulated<log10(RNA_L),log10(RNA_L),simulated)) 
    # leave useful columns
    after_ART<-after_ART %>% 
      dplyr::select(PATIENT,RNA_L,censor,time,limit,simulated,a1_est,a2_est,b1_est,b2_est,b3_est,b4_est)
    # write simulated after ART data
    write.table(after_ART,"after_ART_sim.txt", sep = "," ,quote = FALSE,row.names = FALSE)
    # find new estimates of beta and gamma for simulated viral rebound data
    data = list(dataFile = paste0('after_ART_sim.txt'),
                headerTypes =c("id","ignore","cens","time","limit","observation","contcov","contcov","contcov","contcov","contcov","contcov"),
                observationTypes = list(RE="continuous"))
    modelFile = paste0('after_ART_model.txt')
    newProject(data = data, modelFile =modelFile)
    setIndividualParameterDistribution(beta1="normal",beta2="normal",beta3="normal",beta4="normal")
    # covariates on beta1
    setCovariateModel(beta1 = c(a1_est=TRUE,a2_est=TRUE,b1_est=TRUE,b2_est=TRUE,b3_est=TRUE,b4_est = TRUE))
    setErrorModel(simulated = "constant")
    # run the estimation
    scenario <- getScenario()
    scenario$tasks = c(populationParameterEstimation = T,
                       conditionalModeEstimation = T,
                       conditionalDistributionSampling = T,
                       standardErrorEstimation=T,
                       logLikelihoodEstimation=T)
    scenario$linearization = FALSE
    setScenario(scenario)
    runScenario()
    # store new estimates in the bootstrap dataframe
    popestimates_after <- getEstimatedPopulationParameters()
    bootstrap[7:16,i]<-popestimates_after[1:10]
  }
  # Save the bootstrap estimates in file "bootstrap_beta1.txt"
  write.table(t(bootstrap),glue::glue("bootstrap_beta1.txt"), sep = "," ,quote = FALSE,row.names =  FALSE)
  bootstrap<-read.table(glue::glue("bootstrap_beta1.txt"), sep = ",", header = TRUE)
  # find boostrap se
  bootstrap_se<-sapply(bootstrap, sd)
  # add bootstrap se to the summary table
  summary<-as.data.frame(cbind(summary,bootstrap_se))
  colnames(summary)[3]="Bootstrap SE"
  # find z-value
  summary$zvalue = abs(summary$Estimate/summary$`Bootstrap SE`)
  # find p-value
  summary$`p-value` = pnorm(summary$`zvalue`, lower.tail = FALSE)*2
  # reorder the rows in the summary table
  summary <- summary[c(1:7,10:13, 8:9, 14:nrow(summary)), ]
  summary <- rownames_to_column(summary,"parameters")
  # save summary table in file "summary_beta1.txt")
  write.table(summary, glue::glue("summary_beta1.txt"), sep = "," ,quote = FALSE,row.names =  FALSE)
  View(summary)
}

######################## main function: covariates on beta2 in rebound model
# filename is a string to indicate the raw data file name. 
# Raw data includes "PATIENT","CD4_V","RNA_V","RNA_L","treatment", and "days_from_seroco"
# num_bootstrap is the number of bootstrap want to repeat
# outliers is a vector of PATIENT ID which are outliers, e.g., outliers = c(1,5,10)
# num_months is a number to indicate the time interval for viral rebound analysis. 
analysis_beta2 <- function(filename, num_bootstrap, outliers=NULL, num_months=NULL){
  ######################### Clean raw data
  # Read raw data
  raw_data<-read.csv(filename,header = T)
  # If RNA_V = -1, indicate the observation as censored and replace RNA_V with RNA_L.
  # Take log10 to RNA_V. Convert time from days to months. 
  # Set lower boundary (limit) for the censoring interval to 0.
  # Remove rows with unobserved RNA_V and time
  raw_data <- raw_data %>% 
    mutate(censor = ifelse(RNA_V==-1, 1, 0)) %>% 
    mutate(RNA_V = ifelse(RNA_V==-1, RNA_L, RNA_V)) %>% 
    mutate(log10 = log10(RNA_V)) %>% 
    mutate(time = days_from_seroco/30) %>% 
    mutate(limit = 0) %>% 
    filter(!is.na(RNA_V) & !is.na(time))
  # Remove outliers if there is any
  if(!is.null(outliers)){
    raw_data <- raw_data %>% 
      filter(PATIENT != outliers)
  }
  # Find lowest detection limit for each patient
  raw_data <- raw_data %>%
    group_by(PATIENT) %>%
    mutate(RNA_L = min(RNA_L,na.rm = TRUE)) 
  # Split the dataset into before_ART (treatment=1) and after_ART(treatment=0 after treatment=1). 
  # Note: for after_ART, ignore all "treatment=0" before "treatment=1".
  # Change time after ART to the time since ART interruption 
  # For example: if before_ART time is c(1,2,3), after_ART time is c(3.5,4,5), then 
  #             the time since ART interruption for after_ART is c(0.5,1,2).
  before_ART <- raw_data %>% 
    filter(treatment==1) %>% 
    dplyr::select("PATIENT","RNA_L","CD4_V","censor","log10","time","limit")
  after_ART <- raw_data %>% 
    group_by(PATIENT) %>% 
    mutate(ind = case_when(treatment == 1 ~ time)) %>%  
    mutate(ind = max(ind, na.rm=TRUE)) %>% 
    filter(time > ind) %>% 
    mutate(time = time-ind) %>% 
    dplyr::select("PATIENT","RNA_L","censor","log10","time","limit")
  # Only leave first num_months months data for after_ART, only leave useful columns
  if(!is.null(num_months)){
    after_ART <- after_ART %>% 
      group_by(PATIENT) %>% 
      mutate(min = min(time)) %>% 
      filter(time <= num_months+min) %>% 
      dplyr::select("PATIENT","RNA_L","censor","log10","time","limit")
  }
  
  # Save the two datasets locally. This is required when run Monolix.
  write.csv(before_ART,"before_ART.csv",row.names = FALSE)
  write.csv(after_ART,"after_ART.csv",row.names = FALSE)
  
  ######################## Find population parameter estimates before ART interruption in Monolix
  # model1: y_ij = log10(exp(P_1i-lambda_1i*t_ij)+exp(P_2i-lambda_2i*t_ij))+e_ij
  # P_1i = P_1 + b_1i, P_2i = P_2 + b_2i, lambda_1i = lambda_1+b_3i, lambda_2i = lambda_2+b_4i
  # b_i ~ N(0,B), e_ij ~ N(0, sigma^2)
  data = list(dataFile = paste0('before_ART.csv'),
              headerTypes =c("id","ignore","ignore","cens","observation","time","limit"),
              observationTypes = list(RE="continuous"))
  modelFile = paste0('before_ART_model.txt')
  # create a new project by setting a data set and a structural model
  newProject(data = data, modelFile =modelFile)
  setErrorModel(log10_ = "constant")
  # set tasks in scenario
  scenario <- getScenario()
  scenario$tasks = c(populationParameterEstimation = T,
                     conditionalModeEstimation = T,
                     conditionalDistributionSampling = T,
                     standardErrorEstimation=T,
                     logLikelihoodEstimation=T)
  scenario$linearization = FALSE
  setScenario(scenario)
  setIndividualParameterDistribution(b1="normal",p1="normal",p2="normal",b2="normal")
  # run the estimation
  runScenario()
  # store the estimates
  popestimates_before <- getEstimatedPopulationParameters()
  # store the random effects 
  bi <- as.data.frame(getEstimatedRandomEffects()$conditionalMean)
  # store the s.e.
  se_before <- unlist(getEstimatedStandardErrors())
  # if P1 < P2, switch P1,P2 and lambda1,lambda2
  popestimates_tmp <- popestimates_before
  bi_tmp <- bi
  se_tmp <- se_before
  if (popestimates_before[1]<popestimates_before[3]){
    popestimates_before[1]<-popestimates_tmp[3]
    popestimates_before[3]<-popestimates_tmp[1]
    popestimates_before[2]<-popestimates_tmp[4]
    popestimates_before[4]<-popestimates_tmp[2]
    bi[,4]<-bi_tmp[,2]
    bi[,2]<-bi_tmp[,4]
    bi[,3]<-bi_tmp[,5]
    bi[,5]<-bi_tmp[,3]
    se_before[1]<-se_tmp[3]
    se_before[3]<-se_tmp[1]
    se_before[2]<-se_tmp[4]
    se_before[4]<-se_tmp[2]
  }
  # Find population parameter estimates
  P1 <- popestimates_before[1]
  lambda1 <- popestimates_before[2]
  P2 <- popestimates_before[3]
  lambda2 <- popestimates_before[4]
  # Find sigma
  sigma <- popestimates_before[9]
  # Find random effects
  colnames(bi)<-c("PATIENT","b1_est","b3_est","b2_est","b4_est")
  # Find B
  B<-as.matrix(cov(bi[2:5]))
  
  
  ######################## Find population parameter estimates CD4
  # Read before_ART data
  before_ART <- read.csv("before_ART.csv")
  before_ART$CD4_V<-as.numeric(before_ART$CD4_V)
  long.data<-groupedData(CD4_V~time|PATIENT,data=before_ART)
  long.data$time<-as.numeric(long.data$time)
  long.data$PATIENT<-as.factor(long.data$PATIENT)
  long.data <- long.data %>% 
    filter(!is.na(CD4_V))
  # model2: sqrt(z_ij) = alpha_1i + alpha_2i*t_ij + epsilon_ij
  # alpha_1i = alpha_1 + a_1i, alpha_2i = alpha_2 + a_2i
  # a_i ~ N(0, A), epsilon_ij ~ N(0, delta^2)
  lme <- lme(sqrt(CD4_V)~time, data=long.data,
             random = ~time|PATIENT,method="ML")
  # Find population parameter estimates
  alpha1 <- summary(lme)$coefficients$fixed[1]
  alpha2 <- summary(lme)$coefficients$fixed[2]
  # Find random effects
  ai <- ranef(lme)
  ai <- rownames_to_column(ai, "PATIENT")
  colnames(ai)=c("PATIENT","a1_est","a2_est")
  # Find A
  A <- as.matrix(cov(ai[,2:3]))
  # Find delta
  delta <- summary(lme)$sigma
  # Find s.e.
  se_cd4 <- coef(summary(lme))[,2]
  
  
  ######################## Find population parameter estimates after ART interruption
  # Read after ART data
  after_ART <- read.csv("after_ART.csv")
  # Add random effects from model1 and model2 to after ART data
  after_ART <- merge(after_ART, bi, by="PATIENT")
  after_ART <- merge(after_ART, ai, by="PATIENT")
  # Save after ART data with covariates
  write.csv(after_ART,"after_ART_cov.csv",row.names = FALSE)
  # w_ij = beta_1i + t_ij/(t_ij + exp(beta_2i - beta_3i * t_ij)) + beta_4i + xi_ij
  # beta_2i = beta_2+gamma_1*b_1i+gamma_2*b_2i+gamma_3*b_3i+gamma_4*b_4i+gamma_5*a_1i+gamma_6*a_2i+tau_2i
  # beta_ji = beta_j + tau_ji, j=1,3,4
  # tau_i ~ N(0, G), xi_ij ~ N(0, omega^2)
  data = list(dataFile = paste0('after_ART_cov.csv'),
              headerTypes =c("id","ignore","cens","observation","time","limit","contcov","contcov","contcov","contcov","contcov","contcov"),
              observationTypes = list(RE="continuous"))
  modelFile = paste0('after_ART_model.txt')
  newProject(data = data, modelFile =modelFile)
  setIndividualParameterDistribution(beta1="normal",beta2="normal",beta3="normal",beta4="normal")
  # covariates on beta2
  setCovariateModel(beta2 = c(a1_est=TRUE,a2_est=TRUE,b1_est=TRUE,b2_est=TRUE,b3_est=TRUE,b4_est = TRUE))
  setErrorModel(log10_ = "constant")
  # run the estimation
  scenario <- getScenario()
  scenario$tasks = c(populationParameterEstimation = T,
                     conditionalModeEstimation = T,
                     conditionalDistributionSampling = T,
                     standardErrorEstimation=T,
                     logLikelihoodEstimation=T)
  scenario$linearization = FALSE
  setScenario(scenario)
  runScenario()
  # Find population parameter estimates
  popestimates_after <- getEstimatedPopulationParameters()
  beta1 <- popestimates_after[1]
  beta2 <- popestimates_after[2]
  gamma_5 <- popestimates_after[3]
  gamma_6 <- popestimates_after[4]
  gamma_1 <- popestimates_after[5]
  gamma_2 <- popestimates_after[6]
  gamma_3 <- popestimates_after[7]
  gamma_4 <- popestimates_after[8]
  beta3 <- popestimates_after[9]
  beta4 <- popestimates_after[10]
  # Find omega
  omega <- popestimates_after[15]
  # Find G
  taui <- as.data.frame(getEstimatedRandomEffects()$conditionalMean)
  G <- as.matrix(cov(taui[,2:5]))
  # Find s.e.
  se_after <- unlist(getEstimatedStandardErrors())
  
  
  ######################## Bootstrap
  # Create a empty dataset to store all bootstrap values
  names<-c("P1","b1","P2","b2","alpha1","alpha2","beta1","beta2","g5","g6","g1","g2","g3","g4","beta3","beta4")
  bootstrap<-as.data.frame(matrix(NA,nrow=length(names),ncol=num_bootstrap))
  rownames(bootstrap)<-names
  # Create a summary table, first column is population estimates and second column is original s.e.
  popestimates<-c(popestimates_before[1:4],alpha1,alpha2,popestimates_after[1:10])
  se <- c(se_before[1:4], se_cd4, se_after[1:10])
  se[which(se=="NaN")]<-NA
  se <- as.numeric(se)
  summary <- cbind(popestimates,se)
  colnames(summary)=c("Estimate","Original SE")
  rownames(summary)<-names
  # Find unique patient and number of patients
  n_before<-length(unique(before_ART$PATIENT))
  PATIENT_before<-unique(before_ART$PATIENT)
  n_after<-length(unique(after_ART$PATIENT))
  PATIENT_after<-unique(after_ART$PATIENT)
  # repeat bootstrap num_bootstrap times.
  for (i in 1:num_bootstrap){
    print(glue::glue("i = ",i))
    ############# Before ART 
    # simulate b_i
    b_sim<-rmvnorm(n_before, c(0,0,0,0), B)
    colnames(b_sim)=c("b1_sim","b3_sim","b2_sim","b4_sim")
    b_sim<-cbind(PATIENT_before,b_sim)
    colnames(b_sim)[1] <- "PATIENT"
    # add simulated b_i to before_ART.
    before_ART<-merge(before_ART,b_sim,by="PATIENT")
    # simulate e and y_ij, if y_ij < detection limit, then the y_ij = detection limit
    before_ART <- before_ART %>% 
      mutate(e = rnorm(nrow(before_ART),0,sigma)) %>% 
      mutate(simulated = log10(exp(P1+b1_sim-(P2+b3_sim)*time) + exp ((lambda1+b2_sim) -(lambda2+b4_sim)*time))+e) %>% 
      mutate(simulated = ifelse(simulated<log10(RNA_L), log10(RNA_L), simulated))
    # simulate a_i
    a_sim<-rmvnorm(n_before, c(0,0), A)
    colnames(a_sim)=c("a1_sim","a2_sim")
    a_sim<-cbind(PATIENT_before,a_sim)
    colnames(a_sim)[1] <- "PATIENT"
    # add simulated a_i to before_ART.
    before_ART<-merge(before_ART,a_sim,by="PATIENT")
    # simulate epsilon and z_ij, if z_ij < 1, then z_ij = 1
    before_ART <- before_ART %>% 
      mutate(epsilon = rnorm(nrow(before_ART),0, delta)) %>% 
      mutate(simulatedcd4 = alpha1+a1_sim+time*(alpha2+a2_sim)+epsilon) %>% 
      mutate(simulatedcd4 = ifelse(simulatedcd4<1, 1, simulatedcd4))
    # leave useful columns
    before_ART <- before_ART %>% 
      dplyr::select(PATIENT,RNA_L,censor,time,limit,simulated,simulatedcd4)
    # save before ART simulated data
    write.table(before_ART,"before_ART_sim.txt", sep = "," ,quote = FALSE,row.names = FALSE)
    # Find new estimates of P1,P2,lambda1, and lambda2 for simulated viral loads
    data = list(dataFile = paste0('before_ART_sim.txt'),
                headerTypes =c("id","ignore","cens","time","limit", "observation","ignore"),
                observationTypes = list(RE="continuous"))
    modelFile = paste0('before_ART_model.txt')
    newProject(data = data, modelFile =modelFile)
    setErrorModel(simulated = "constant")
    # set tasks in scenario
    scenario <- getScenario()
    scenario$tasks = c(populationParameterEstimation = T,
                       conditionalModeEstimation = T,
                       conditionalDistributionSampling = T,
                       standardErrorEstimation=T,
                       logLikelihoodEstimation=T)
    scenario$linearization = FALSE
    setScenario(scenario)
    setIndividualParameterDistribution(b1="normal",p1="normal",p2="normal",b2="normal")
    # run the estimation
    runScenario()
    # find population parameter estimates
    popestimates_before <- getEstimatedPopulationParameters()
    # find random effects
    bi <- as.data.frame(getEstimatedRandomEffects()$conditionalMean)
    # if P1 < P2, switch P1,P2 and lambda1,lambda2
    popestimates_tmp <- popestimates_before
    bi_tmp <- bi
    if (popestimates_before[1]<popestimates_before[3]){
      popestimates_before[1]<-popestimates_tmp[3]
      popestimates_before[3]<-popestimates_tmp[1]
      popestimates_before[2]<-popestimates_tmp[4]
      popestimates_before[4]<-popestimates_tmp[2]
      bi[,4]<-bi_tmp[,2]
      bi[,2]<-bi_tmp[,4]
      bi[,3]<-bi_tmp[,5]
      bi[,5]<-bi_tmp[,3]
    }
    # store the estimates in bootstrap dataframe
    bootstrap[1:4,i]=popestimates_before[1:4]
    # find the random effects bi
    colnames(bi)<-c("PATIENT","b1_est","b3_est","b2_est","b4_est")
    
    ############# CD4 Data
    long.data<-groupedData(simulatedcd4~time|PATIENT,data=before_ART)
    long.data$time<-as.numeric(long.data$time)
    long.data$PATIENT<-as.factor(long.data$PATIENT)
    # If the lme does not converge, then re-simulate cd4 data
    while (has_error (lme <- lme(simulatedcd4~time, data=long.data, 
                                 random = ~time|PATIENT,method="ML"))){
      a_sim<-rmvnorm(n_after, c(0,0), A)
      colnames(a_sim)=c("a1_sim","a2_sim")
      a_sim<-cbind(PATIENT_before,a_sim)
      colnames(a_sim)[1] <- "PATIENT"
      before_ART<-merge(before_ART,a_sim,by="PATIENT")
      before_ART <- before_ART %>% 
        mutate(epsilon = rnorm(nrow(before_ART),0,delta)) %>% 
        mutate(simulatedcd4 = alpha1+a1_sim+time*(alpha2+a2_sim)+epsilon) %>% 
        mutate(simulatedcd4 = ifelse(simulatedcd4<1, 1, simulatedcd4))
      before_ART <- before_ART %>% 
        dplyr::select(PATIENT,RNA_L,censor,time,limit,simulated,simulatedcd4)
      long.data<-groupedData(simulatedcd4~time|PATIENT,data=before_ART)
      long.data$time<-as.numeric(long.data$time)
      long.data$PATIENT<-as.factor(long.data$PATIENT)
    }
    # store the new estimates in bootstrap dataframe
    bootstrap[5:6,i]<-summary(lme)$coefficients$fixed
    # find the random effects ai
    ai<-ranef(lme)
    ai <- rownames_to_column(ai,"PATIENT")
    colnames(ai)[2:3]<-c("a1_est","a2_est")
    
    ############# Bootstrap before ART 
    # leave useful columns
    after_ART <- after_ART %>% 
      dplyr::select("PATIENT","RNA_L","censor","time","limit")
    # simulate tau_i
    taui<-rmvnorm(n_after, c(0,0,0,0), G)
    colnames(taui)=c("tau1","tau2","tau3","tau4")
    taui<-cbind(PATIENT_after,taui)
    colnames(taui)[1]<-"PATIENT"
    # Add random effects from previous two models and simulated taui to after ART data
    after_ART<-merge(after_ART,taui,by="PATIENT")
    after_ART <- merge(after_ART,bi,by="PATIENT")
    after_ART<-merge(after_ART,ai,by="PATIENT")
    # simulate xi_ij 
    after_ART <- after_ART %>% 
      mutate(xi = rnorm(nrow(after_ART),0,omega))
    # simulate w_ij, if w_ij < detection limit, then w_ij = detection limit.
    after_ART <- after_ART %>% 
      mutate(simulated = (beta1+tau1) * time / ( time + exp( beta2+a1_est*gamma_5+a2_est*gamma_6+b1_est*gamma_1+b2_est*gamma_2+b3_est*gamma_3+b4_est*gamma_4+tau2 - (beta3+tau3) * time )) + (beta4+tau4)+xi)
    # if w_ij < detection limit, then w_ij = detection limit.
    after_ART <- after_ART %>% 
      mutate(simulated = ifelse(simulated<log10(RNA_L),log10(RNA_L),simulated)) 
    # leave useful columns
    after_ART<-after_ART %>% 
      dplyr::select(PATIENT,RNA_L,censor,time,limit,simulated,a1_est,a2_est,b1_est,b2_est,b3_est,b4_est)
    # write simulated after ART data
    write.table(after_ART,"after_ART_sim.txt", sep = "," ,quote = FALSE,row.names = FALSE)
    # find new estimates of beta and gamma for simulated viral rebound data
    data = list(dataFile = paste0('after_ART_sim.txt'),
                headerTypes =c("id","ignore","cens","time","limit","observation","contcov","contcov","contcov","contcov","contcov","contcov"),
                observationTypes = list(RE="continuous"))
    modelFile = paste0('after_ART_model.txt')
    newProject(data = data, modelFile =modelFile)
    setIndividualParameterDistribution(beta1="normal",beta2="normal",beta3="normal",beta4="normal")
    # covariates on beta2
    setCovariateModel(beta2 = c(a1_est=TRUE,a2_est=TRUE,b1_est=TRUE,b2_est=TRUE,b3_est=TRUE,b4_est = TRUE))
    setErrorModel(simulated = "constant")
    # run the estimation
    scenario <- getScenario()
    scenario$tasks = c(populationParameterEstimation = T,
                       conditionalModeEstimation = T,
                       conditionalDistributionSampling = T,
                       standardErrorEstimation=T,
                       logLikelihoodEstimation=T)
    scenario$linearization = FALSE
    setScenario(scenario)
    runScenario()
    # store new estimates in the bootstrap dataframe
    popestimates_after <- getEstimatedPopulationParameters()
    bootstrap[7:16,i]<-popestimates_after[1:10]
  }
  # Save the bootstrap estimates in file "bootstrap_beta2.txt"
  write.table(t(bootstrap),glue::glue("bootstrap_beta2.txt"), sep = "," ,quote = FALSE,row.names =  FALSE)
  bootstrap<-read.table(glue::glue("bootstrap_beta2.txt"), sep = ",", header = TRUE)
  # find boostrap se
  bootstrap_se<-sapply(bootstrap, sd)
  # add bootstrap se to the summary table
  summary<-as.data.frame(cbind(summary,bootstrap_se))
  colnames(summary)[3]="Bootstrap SE"
  # find z-value
  summary$zvalue = abs(summary$Estimate/summary$`Bootstrap SE`)
  # find p-value
  summary$`p-value` = pnorm(summary$`zvalue`, lower.tail = FALSE)*2
  # reorder the rows in the summary table
  summary <- summary[c(1:8,11:14, 9:10,15:nrow(summary)), ]
  summary <- rownames_to_column(summary,"parameters")
  # save summary table in file "summary_beta2.txt")
  write.table(summary, glue::glue("summary_beta2.txt"), sep = "," ,quote = FALSE,row.names =  FALSE)
  View(summary)
}


######################## main function: covariates on beta3 in rebound model
# filename is a string to indicate the raw data file name. 
# Raw data includes "PATIENT","CD4_V","RNA_V","RNA_L","treatment", and "days_from_seroco"
# num_bootstrap is the number of bootstrap want to repeat
# outliers is a vector of PATIENT ID which are outliers, e.g., outliers = c(1,5,10)
# num_months is a number to indicate the time interval for viral rebound analysis. 
analysis_beta3 <- function(filename, num_bootstrap, outliers=NULL, num_months=NULL){
  ######################### Clean raw data
  # Read raw data
  raw_data<-read.csv(filename,header = T)
  # If RNA_V = -1, indicate the observation as censored and replace RNA_V with RNA_L.
  # Take log10 to RNA_V. Convert time from days to months. 
  # Set lower boundary (limit) for the censoring interval to 0.
  # Remove rows with unobserved RNA_V and time
  raw_data <- raw_data %>% 
    mutate(censor = ifelse(RNA_V==-1, 1, 0)) %>% 
    mutate(RNA_V = ifelse(RNA_V==-1, RNA_L, RNA_V)) %>% 
    mutate(log10 = log10(RNA_V)) %>% 
    mutate(time = days_from_seroco/30) %>% 
    mutate(limit = 0) %>% 
    filter(!is.na(RNA_V) & !is.na(time))
  # Remove outliers if there is any
  if(!is.null(outliers)){
    raw_data <- raw_data %>% 
      filter(PATIENT != outliers)
  }
  # Find lowest detection limit for each patient
  raw_data <- raw_data %>%
    group_by(PATIENT) %>%
    mutate(RNA_L = min(RNA_L,na.rm = TRUE)) 
  # Split the dataset into before_ART (treatment=1) and after_ART(treatment=0 after treatment=1). 
  # Note: for after_ART, ignore all "treatment=0" before "treatment=1".
  # Change time after ART to the time since ART interruption 
  # For example: if before_ART time is c(1,2,3), after_ART time is c(3.5,4,5), then 
  #             the time since ART interruption for after_ART is c(0.5,1,2).
  before_ART <- raw_data %>% 
    filter(treatment==1) %>% 
    dplyr::select("PATIENT","RNA_L","CD4_V","censor","log10","time","limit")
  after_ART <- raw_data %>% 
    group_by(PATIENT) %>% 
    mutate(ind = case_when(treatment == 1 ~ time)) %>%  
    mutate(ind = max(ind, na.rm=TRUE)) %>% 
    filter(time > ind) %>% 
    mutate(time = time-ind) %>% 
    dplyr::select("PATIENT","RNA_L","censor","log10","time","limit")
  # Only leave first num_months months data for after_ART, only leave useful columns
  if(!is.null(num_months)){
    after_ART <- after_ART %>% 
      group_by(PATIENT) %>% 
      mutate(min = min(time)) %>% 
      filter(time <= num_months+min) %>% 
      dplyr::select("PATIENT","RNA_L","censor","log10","time","limit")
  }
  
  # Save the two datasets locally. This is required when run Monolix.
  write.csv(before_ART,"before_ART.csv",row.names = FALSE)
  write.csv(after_ART,"after_ART.csv",row.names = FALSE)
  
  ######################## Find population parameter estimates before ART interruption in Monolix
  # model1: y_ij = log10(exp(P_1i-lambda_1i*t_ij)+exp(P_2i-lambda_2i*t_ij))+e_ij
  # P_1i = P_1 + b_1i, P_2i = P_2 + b_2i, lambda_1i = lambda_1+b_3i, lambda_2i = lambda_2+b_4i
  # b_i ~ N(0,B), e_ij ~ N(0, sigma^2)
  data = list(dataFile = paste0('before_ART.csv'),
              headerTypes =c("id","ignore","ignore","cens","observation","time","limit"),
              observationTypes = list(RE="continuous"))
  modelFile = paste0('before_ART_model.txt')
  # create a new project by setting a data set and a structural model
  newProject(data = data, modelFile =modelFile)
  setErrorModel(log10_ = "constant")
  # set tasks in scenario
  scenario <- getScenario()
  scenario$tasks = c(populationParameterEstimation = T,
                     conditionalModeEstimation = T,
                     conditionalDistributionSampling = T,
                     standardErrorEstimation=T,
                     logLikelihoodEstimation=T)
  scenario$linearization = FALSE
  setScenario(scenario)
  setIndividualParameterDistribution(b1="normal",p1="normal",p2="normal",b2="normal")
  # run the estimation
  runScenario()
  # store the estimates
  popestimates_before <- getEstimatedPopulationParameters()
  # store the random effects 
  bi <- as.data.frame(getEstimatedRandomEffects()$conditionalMean)
  # store the s.e.
  se_before <- unlist(getEstimatedStandardErrors())
  # if P1 < P2, switch P1,P2 and lambda1,lambda2
  popestimates_tmp <- popestimates_before
  bi_tmp <- bi
  se_tmp <- se_before
  if (popestimates_before[1]<popestimates_before[3]){
    popestimates_before[1]<-popestimates_tmp[3]
    popestimates_before[3]<-popestimates_tmp[1]
    popestimates_before[2]<-popestimates_tmp[4]
    popestimates_before[4]<-popestimates_tmp[2]
    bi[,4]<-bi_tmp[,2]
    bi[,2]<-bi_tmp[,4]
    bi[,3]<-bi_tmp[,5]
    bi[,5]<-bi_tmp[,3]
    se_before[1]<-se_tmp[3]
    se_before[3]<-se_tmp[1]
    se_before[2]<-se_tmp[4]
    se_before[4]<-se_tmp[2]
  }
  # Find population parameter estimates
  P1 <- popestimates_before[1]
  lambda1 <- popestimates_before[2]
  P2 <- popestimates_before[3]
  lambda2 <- popestimates_before[4]
  # Find sigma
  sigma <- popestimates_before[9]
  # Find random effects
  colnames(bi)<-c("PATIENT","b1_est","b3_est","b2_est","b4_est")
  # Find B
  B<-as.matrix(cov(bi[2:5]))
  
  
  ######################## Find population parameter estimates CD4
  # Read before_ART data
  before_ART <- read.csv("before_ART.csv")
  before_ART$CD4_V<-as.numeric(before_ART$CD4_V)
  long.data<-groupedData(CD4_V~time|PATIENT,data=before_ART)
  long.data$time<-as.numeric(long.data$time)
  long.data$PATIENT<-as.factor(long.data$PATIENT)
  long.data <- long.data %>% 
    filter(!is.na(CD4_V))
  # model2: sqrt(z_ij) = alpha_1i + alpha_2i*t_ij + epsilon_ij
  # alpha_1i = alpha_1 + a_1i, alpha_2i = alpha_2 + a_2i
  # a_i ~ N(0, A), epsilon_ij ~ N(0, delta^2)
  lme <- lme(sqrt(CD4_V)~time, data=long.data,
             random = ~time|PATIENT,method="ML")
  # Find population parameter estimates
  alpha1 <- summary(lme)$coefficients$fixed[1]
  alpha2 <- summary(lme)$coefficients$fixed[2]
  # Find random effects
  ai <- ranef(lme)
  ai <- rownames_to_column(ai, "PATIENT")
  colnames(ai)=c("PATIENT","a1_est","a2_est")
  # Find A
  A <- as.matrix(cov(ai[,2:3]))
  # Find delta
  delta <- summary(lme)$sigma
  # Find s.e.
  se_cd4 <- coef(summary(lme))[,2]
  
  
  ######################## Find population parameter estimates after ART interruption
  # Read after ART data
  after_ART <- read.csv("after_ART.csv")
  # Add random effects from model1 and model2 to after ART data
  after_ART <- merge(after_ART, bi, by="PATIENT")
  after_ART <- merge(after_ART, ai, by="PATIENT")
  # Save after ART data with covariates
  write.csv(after_ART,"after_ART_cov.csv",row.names = FALSE)
  # w_ij = beta_1i + t_ij/(t_ij + exp(beta_2i - beta_3i * t_ij)) + beta_4i + xi_ij
  # beta_3i = beta_3+gamma_1*b_1i+gamma_2*b_2i+gamma_3*b_3i+gamma_4*b_4i+gamma_5*a_1i+gamma_6*a_2i+tau_3i
  # beta_ji = beta_j + tau_ji, j=1,2,4
  # tau_i ~ N(0, G), xi_ij ~ N(0, omega^2)
  data = list(dataFile = paste0('after_ART_cov.csv'),
              headerTypes =c("id","ignore","cens","observation","time","limit","contcov","contcov","contcov","contcov","contcov","contcov"),
              observationTypes = list(RE="continuous"))
  modelFile = paste0('after_ART_model.txt')
  newProject(data = data, modelFile =modelFile)
  setIndividualParameterDistribution(beta1="normal",beta2="normal",beta3="normal",beta4="normal")
  # covariates on beta3
  setCovariateModel(beta3 = c(a1_est=TRUE,a2_est=TRUE,b1_est=TRUE,b2_est=TRUE,b3_est=TRUE,b4_est = TRUE))
  setErrorModel(log10_ = "constant")
  # run the estimation
  scenario <- getScenario()
  scenario$tasks = c(populationParameterEstimation = T,
                     conditionalModeEstimation = T,
                     conditionalDistributionSampling = T,
                     standardErrorEstimation=T,
                     logLikelihoodEstimation=T)
  scenario$linearization = FALSE
  setScenario(scenario)
  runScenario()
  # Find population parameter estimates
  popestimates_after <- getEstimatedPopulationParameters()
  beta1 <- popestimates_after[1]
  beta2 <- popestimates_after[2]
  beta3 <- popestimates_after[3]
  gamma_5 <- popestimates_after[4]
  gamma_6 <- popestimates_after[5]
  gamma_1 <- popestimates_after[6]
  gamma_2 <- popestimates_after[7]
  gamma_3 <- popestimates_after[8]
  gamma_4 <- popestimates_after[9]
  beta4 <- popestimates_after[10]
  # Find omega
  omega <- popestimates_after[15]
  # Find G
  taui <- as.data.frame(getEstimatedRandomEffects()$conditionalMean)
  G <- as.matrix(cov(taui[,2:5]))
  # Find s.e.
  se_after <- unlist(getEstimatedStandardErrors())
  
  
  ######################## Bootstrap
  # Create a empty dataset to store all bootstrap values
  names<-c("P1","b1","P2","b2","alpha1","alpha2","beta1","beta2","beta3","g5","g6","g1","g2","g3","g4","beta4")
  bootstrap<-as.data.frame(matrix(NA,nrow=length(names),ncol=num_bootstrap))
  rownames(bootstrap)<-names
  # Create a summary table, first column is population estimates and second column is original s.e.
  popestimates<-c(popestimates_before[1:4],alpha1,alpha2,popestimates_after[1:10])
  se <- c(se_before[1:4], se_cd4, se_after[1:10])
  se[which(se=="NaN")]<-NA
  se <- as.numeric(se)
  summary <- cbind(popestimates,se)
  colnames(summary)=c("Estimate","Original SE")
  rownames(summary)<-names
  # Find unique patient and number of patients
  n_before<-length(unique(before_ART$PATIENT))
  PATIENT_before<-unique(before_ART$PATIENT)
  n_after<-length(unique(after_ART$PATIENT))
  PATIENT_after<-unique(after_ART$PATIENT)
  # repeat bootstrap num_bootstrap times.
  for (i in 1:num_bootstrap){
    print(glue::glue("i = ",i))
    ############# Before ART 
    # simulate b_i
    b_sim<-rmvnorm(n_before, c(0,0,0,0), B)
    colnames(b_sim)=c("b1_sim","b3_sim","b2_sim","b4_sim")
    b_sim<-cbind(PATIENT_before,b_sim)
    colnames(b_sim)[1] <- "PATIENT"
    # add simulated b_i to before_ART.
    before_ART<-merge(before_ART,b_sim,by="PATIENT")
    # simulate e and y_ij, if y_ij < detection limit, then the y_ij = detection limit
    before_ART <- before_ART %>% 
      mutate(e = rnorm(nrow(before_ART),0,sigma)) %>% 
      mutate(simulated = log10(exp(P1+b1_sim-(P2+b3_sim)*time) + exp ((lambda1+b2_sim) -(lambda2+b4_sim)*time))+e) %>% 
      mutate(simulated = ifelse(simulated<log10(RNA_L), log10(RNA_L), simulated))
    # simulate a_i
    a_sim<-rmvnorm(n_before, c(0,0), A)
    colnames(a_sim)=c("a1_sim","a2_sim")
    a_sim<-cbind(PATIENT_before,a_sim)
    colnames(a_sim)[1] <- "PATIENT"
    # add simulated a_i to before_ART.
    before_ART<-merge(before_ART,a_sim,by="PATIENT")
    # simulate epsilon and z_ij, if z_ij < 1, then z_ij = 1
    before_ART <- before_ART %>% 
      mutate(epsilon = rnorm(nrow(before_ART),0, delta)) %>% 
      mutate(simulatedcd4 = alpha1+a1_sim+time*(alpha2+a2_sim)+epsilon) %>% 
      mutate(simulatedcd4 = ifelse(simulatedcd4<1, 1, simulatedcd4))
    # leave useful columns
    before_ART <- before_ART %>% 
      dplyr::select(PATIENT,RNA_L,censor,time,limit,simulated,simulatedcd4)
    # save before ART simulated data
    write.table(before_ART,"before_ART_sim.txt", sep = "," ,quote = FALSE,row.names = FALSE)
    # Find new estimates of P1,P2,lambda1, and lambda2 for simulated viral loads
    data = list(dataFile = paste0('before_ART_sim.txt'),
                headerTypes =c("id","ignore","cens","time","limit", "observation","ignore"),
                observationTypes = list(RE="continuous"))
    modelFile = paste0('before_ART_model.txt')
    newProject(data = data, modelFile =modelFile)
    setErrorModel(simulated = "constant")
    # set tasks in scenario
    scenario <- getScenario()
    scenario$tasks = c(populationParameterEstimation = T,
                       conditionalModeEstimation = T,
                       conditionalDistributionSampling = T,
                       standardErrorEstimation=T,
                       logLikelihoodEstimation=T)
    scenario$linearization = FALSE
    setScenario(scenario)
    setIndividualParameterDistribution(b1="normal",p1="normal",p2="normal",b2="normal")
    # run the estimation
    runScenario()
    # find population parameter estimates
    popestimates_before <- getEstimatedPopulationParameters()
    # find random effects
    bi <- as.data.frame(getEstimatedRandomEffects()$conditionalMean)
    # if P1 < P2, switch P1,P2 and lambda1,lambda2
    popestimates_tmp <- popestimates_before
    bi_tmp <- bi
    if (popestimates_before[1]<popestimates_before[3]){
      popestimates_before[1]<-popestimates_tmp[3]
      popestimates_before[3]<-popestimates_tmp[1]
      popestimates_before[2]<-popestimates_tmp[4]
      popestimates_before[4]<-popestimates_tmp[2]
      bi[,4]<-bi_tmp[,2]
      bi[,2]<-bi_tmp[,4]
      bi[,3]<-bi_tmp[,5]
      bi[,5]<-bi_tmp[,3]
    }
    # store the estimates in bootstrap dataframe
    bootstrap[1:4,i]=popestimates_before[1:4]
    # find the random effects bi
    colnames(bi)<-c("PATIENT","b1_est","b3_est","b2_est","b4_est")
    
    ############# CD4 Data
    long.data<-groupedData(simulatedcd4~time|PATIENT,data=before_ART)
    long.data$time<-as.numeric(long.data$time)
    long.data$PATIENT<-as.factor(long.data$PATIENT)
    # If the lme does not converge, then re-simulate cd4 data
    while (has_error (lme <- lme(simulatedcd4~time, data=long.data, 
                                 random = ~time|PATIENT,method="ML"))){
      a_sim<-rmvnorm(n_after, c(0,0), A)
      colnames(a_sim)=c("a1_sim","a2_sim")
      a_sim<-cbind(PATIENT_before,a_sim)
      colnames(a_sim)[1] <- "PATIENT"
      before_ART<-merge(before_ART,a_sim,by="PATIENT")
      before_ART <- before_ART %>% 
        mutate(epsilon = rnorm(nrow(before_ART),0,delta)) %>% 
        mutate(simulatedcd4 = alpha1+a1_sim+time*(alpha2+a2_sim)+epsilon) %>% 
        mutate(simulatedcd4 = ifelse(simulatedcd4<1, 1, simulatedcd4))
      before_ART <- before_ART %>% 
        dplyr::select(PATIENT,RNA_L,censor,time,limit,simulated,simulatedcd4)
      long.data<-groupedData(simulatedcd4~time|PATIENT,data=before_ART)
      long.data$time<-as.numeric(long.data$time)
      long.data$PATIENT<-as.factor(long.data$PATIENT)
    }
    # store the new estimates in bootstrap dataframe
    bootstrap[5:6,i]<-summary(lme)$coefficients$fixed
    # find the random effects ai
    ai<-ranef(lme)
    ai <- rownames_to_column(ai,"PATIENT")
    colnames(ai)[2:3]<-c("a1_est","a2_est")
    
    ############# Bootstrap before ART 
    # leave useful columns
    after_ART <- after_ART %>% 
      dplyr::select("PATIENT","RNA_L","censor","time","limit")
    # simulate tau_i
    taui<-rmvnorm(n_after, c(0,0,0,0), G)
    colnames(taui)=c("tau1","tau2","tau3","tau4")
    taui<-cbind(PATIENT_after,taui)
    colnames(taui)[1]<-"PATIENT"
    # Add random effects from previous two models and simulated taui to after ART data
    after_ART<-merge(after_ART,taui,by="PATIENT")
    after_ART <- merge(after_ART,bi,by="PATIENT")
    after_ART<-merge(after_ART,ai,by="PATIENT")
    # simulate xi_ij 
    after_ART <- after_ART %>% 
      mutate(xi = rnorm(nrow(after_ART),0,omega))
    # simulate w_ij, if w_ij < detection limit, then w_ij = detection limit.
    after_ART <- after_ART %>% 
      mutate(simulated = (beta1+tau1) * time / ( time + exp( beta2+tau2 - (beta3+a1_est*gamma_5+a2_est*gamma_6+b1_est*gamma_1+b2_est*gamma_2+b3_est*gamma_3+b4_est*gamma_4+tau3) * time )) + (beta4+tau4)+xi)
    # if w_ij < detection limit, then w_ij = detection limit.
    after_ART <- after_ART %>% 
      mutate(simulated = ifelse(simulated<log10(RNA_L),log10(RNA_L),simulated)) 
    # leave useful columns
    after_ART<-after_ART %>% 
      dplyr::select(PATIENT,RNA_L,censor,time,limit,simulated,a1_est,a2_est,b1_est,b2_est,b3_est,b4_est)
    # write simulated after ART data
    write.table(after_ART,"after_ART_sim.txt", sep = "," ,quote = FALSE,row.names = FALSE)
    # find new estimates of beta and gamma for simulated viral rebound data
    data = list(dataFile = paste0('after_ART_sim.txt'),
                headerTypes =c("id","ignore","cens","time","limit","observation","contcov","contcov","contcov","contcov","contcov","contcov"),
                observationTypes = list(RE="continuous"))
    modelFile = paste0('after_ART_model.txt')
    newProject(data = data, modelFile =modelFile)
    setIndividualParameterDistribution(beta1="normal",beta2="normal",beta3="normal",beta4="normal")
    # covariates on beta3
    setCovariateModel(beta3 = c(a1_est=TRUE,a2_est=TRUE,b1_est=TRUE,b2_est=TRUE,b3_est=TRUE,b4_est = TRUE))
    setErrorModel(simulated = "constant")
    # run the estimation
    scenario <- getScenario()
    scenario$tasks = c(populationParameterEstimation = T,
                       conditionalModeEstimation = T,
                       conditionalDistributionSampling = T,
                       standardErrorEstimation=T,
                       logLikelihoodEstimation=T)
    scenario$linearization = FALSE
    setScenario(scenario)
    runScenario()
    # store new estimates in the bootstrap dataframe
    popestimates_after <- getEstimatedPopulationParameters()
    bootstrap[7:16,i]<-popestimates_after[1:10]
  }
  # Save the bootstrap estimates in file "bootstrap_beta3.txt"
  write.table(t(bootstrap),glue::glue("bootstrap_beta3.txt"), sep = "," ,quote = FALSE,row.names =  FALSE)
  bootstrap<-read.table(glue::glue("bootstrap_beta3.txt"), sep = ",", header = TRUE)
  # find boostrap se
  bootstrap_se<-sapply(bootstrap, sd)
  # add bootstrap se to the summary table
  summary<-as.data.frame(cbind(summary,bootstrap_se))
  colnames(summary)[3]="Bootstrap SE"
  # find z-value
  summary$zvalue = abs(summary$Estimate/summary$`Bootstrap SE`)
  # find p-value
  summary$`p-value` = pnorm(summary$`zvalue`, lower.tail = FALSE)*2
  # reorder the rows in the summary table
  summary <- summary[c(1:9,12:15, 10:11,16:nrow(summary)), ]
  summary <- rownames_to_column(summary,"parameters")
  # save summary table in file "summary_beta3.txt")
  write.table(summary, glue::glue("summary_beta3.txt"), sep = "," ,quote = FALSE,row.names =  FALSE)
  View(summary)
}


######################## main function: covariates on beta4 in rebound model
# filename is a string to indicate the raw data file name. 
# Raw data includes "PATIENT","CD4_V","RNA_V","RNA_L","treatment", and "days_from_seroco"
# num_bootstrap is the number of bootstrap want to repeat
# outliers is a vector of PATIENT ID which are outliers, e.g., outliers = c(1,5,10)
# num_months is a number to indicate the time interval for viral rebound analysis. 
analysis_beta4 <- function(filename, num_bootstrap, outliers=NULL, num_months=NULL){ 
  ######################### Clean raw data
  # Read raw data
  raw_data<-read.csv(filename,header = T)
  # If RNA_V = -1, indicate the observation as censored and replace RNA_V with RNA_L.
  # Take log10 to RNA_V. Convert time from days to months. 
  # Set lower boundary (limit) for the censoring interval to 0.
  # Remove rows with unobserved RNA_V and time
  raw_data <- raw_data %>% 
    mutate(censor = ifelse(RNA_V==-1, 1, 0)) %>% 
    mutate(RNA_V = ifelse(RNA_V==-1, RNA_L, RNA_V)) %>% 
    mutate(log10 = log10(RNA_V)) %>% 
    mutate(time = days_from_seroco/30) %>% 
    mutate(limit = 0) %>% 
    filter(!is.na(RNA_V) & !is.na(time))
  # Remove outliers if there is any
  if(!is.null(outliers)){
    raw_data <- raw_data %>% 
      filter(PATIENT != outliers)
  }
  # Find lowest detection limit for each patient
  raw_data <- raw_data %>%
    group_by(PATIENT) %>%
    mutate(RNA_L = min(RNA_L,na.rm = TRUE)) 
  # Split the dataset into before_ART (treatment=1) and after_ART(treatment=0 after treatment=1). 
  # Note: for after_ART, ignore all "treatment=0" before "treatment=1".
  # Change time after ART to the time since ART interruption 
  # For example: if before_ART time is c(1,2,3), after_ART time is c(3.5,4,5), then 
  #             the time since ART interruption for after_ART is c(0.5,1,2).
  before_ART <- raw_data %>% 
    filter(treatment==1) %>% 
    dplyr::select("PATIENT","RNA_L","CD4_V","censor","log10","time","limit")
  after_ART <- raw_data %>% 
    group_by(PATIENT) %>% 
    mutate(ind = case_when(treatment == 1 ~ time)) %>%  
    mutate(ind = max(ind, na.rm=TRUE)) %>% 
    filter(time > ind) %>% 
    mutate(time = time-ind) %>% 
    dplyr::select("PATIENT","RNA_L","censor","log10","time","limit")
  # Only leave first num_months months data for after_ART, only leave useful columns
  if(!is.null(num_months)){
    after_ART <- after_ART %>% 
      group_by(PATIENT) %>% 
      mutate(min = min(time)) %>% 
      filter(time <= num_months+min) %>% 
      dplyr::select("PATIENT","RNA_L","censor","log10","time","limit")
  }
  
  # Save the two datasets locally. This is required when run Monolix.
  write.csv(before_ART,"before_ART.csv",row.names = FALSE)
  write.csv(after_ART,"after_ART.csv",row.names = FALSE)
  
  ######################## Find population parameter estimates before ART interruption in Monolix
  # model1: y_ij = log10(exp(P_1i-lambda_1i*t_ij)+exp(P_2i-lambda_2i*t_ij))+e_ij
  # P_1i = P_1 + b_1i, P_2i = P_2 + b_2i, lambda_1i = lambda_1+b_3i, lambda_2i = lambda_2+b_4i
  # b_i ~ N(0,B), e_ij ~ N(0, sigma^2)
  data = list(dataFile = paste0('before_ART.csv'),
              headerTypes =c("id","ignore","ignore","cens","observation","time","limit"),
              observationTypes = list(RE="continuous"))
  modelFile = paste0('before_ART_model.txt')
  # create a new project by setting a data set and a structural model
  newProject(data = data, modelFile =modelFile)
  setErrorModel(log10_ = "constant")
  # set tasks in scenario
  scenario <- getScenario()
  scenario$tasks = c(populationParameterEstimation = T,
                     conditionalModeEstimation = T,
                     conditionalDistributionSampling = T,
                     standardErrorEstimation=T,
                     logLikelihoodEstimation=T)
  scenario$linearization = FALSE
  setScenario(scenario)
  setIndividualParameterDistribution(b1="normal",p1="normal",p2="normal",b2="normal")
  # run the estimation
  runScenario()
  # store the estimates
  popestimates_before <- getEstimatedPopulationParameters()
  # store the random effects 
  bi <- as.data.frame(getEstimatedRandomEffects()$conditionalMean)
  # store the s.e.
  se_before <- unlist(getEstimatedStandardErrors())
  # if P1 < P2, switch P1,P2 and lambda1,lambda2
  popestimates_tmp <- popestimates_before
  bi_tmp <- bi
  se_tmp <- se_before
  if (popestimates_before[1]<popestimates_before[3]){
    popestimates_before[1]<-popestimates_tmp[3]
    popestimates_before[3]<-popestimates_tmp[1]
    popestimates_before[2]<-popestimates_tmp[4]
    popestimates_before[4]<-popestimates_tmp[2]
    bi[,4]<-bi_tmp[,2]
    bi[,2]<-bi_tmp[,4]
    bi[,3]<-bi_tmp[,5]
    bi[,5]<-bi_tmp[,3]
    se_before[1]<-se_tmp[3]
    se_before[3]<-se_tmp[1]
    se_before[2]<-se_tmp[4]
    se_before[4]<-se_tmp[2]
  }
  # Find population parameter estimates
  P1 <- popestimates_before[1]
  lambda1 <- popestimates_before[2]
  P2 <- popestimates_before[3]
  lambda2 <- popestimates_before[4]
  # Find sigma
  sigma <- popestimates_before[9]
  # Find random effects
  colnames(bi)<-c("PATIENT","b1_est","b3_est","b2_est","b4_est")
  # Find B
  B<-as.matrix(cov(bi[2:5]))
  
  
  ######################## Find population parameter estimates CD4
  # Read before_ART data
  before_ART <- read.csv("before_ART.csv")
  before_ART$CD4_V<-as.numeric(before_ART$CD4_V)
  long.data<-groupedData(CD4_V~time|PATIENT,data=before_ART)
  long.data$time<-as.numeric(long.data$time)
  long.data$PATIENT<-as.factor(long.data$PATIENT)
  long.data <- long.data %>% 
    filter(!is.na(CD4_V))
  # model2: sqrt(z_ij) = alpha_1i + alpha_2i*t_ij + epsilon_ij
  # alpha_1i = alpha_1 + a_1i, alpha_2i = alpha_2 + a_2i
  # a_i ~ N(0, A), epsilon_ij ~ N(0, delta^2)
  lme <- lme(sqrt(CD4_V)~time, data=long.data,
             random = ~time|PATIENT,method="ML")
  # Find population parameter estimates
  alpha1 <- summary(lme)$coefficients$fixed[1]
  alpha2 <- summary(lme)$coefficients$fixed[2]
  # Find random effects
  ai <- ranef(lme)
  ai <- rownames_to_column(ai, "PATIENT")
  colnames(ai)=c("PATIENT","a1_est","a2_est")
  # Find A
  A <- as.matrix(cov(ai[,2:3]))
  # Find delta
  delta <- summary(lme)$sigma
  # Find s.e.
  se_cd4 <- coef(summary(lme))[,2]
  
  
  ######################## Find population parameter estimates after ART interruption
  # Read after ART data
  after_ART <- read.csv("after_ART.csv")
  # Add random effects from model1 and model2 to after ART data
  after_ART <- merge(after_ART, bi, by="PATIENT")
  after_ART <- merge(after_ART, ai, by="PATIENT")
  # Save after ART data with covariates
  write.csv(after_ART,"after_ART_cov.csv",row.names = FALSE)
  # w_ij = beta_1i + t_ij/(t_ij + exp(beta_2i - beta_3i * t_ij)) + beta_4i + xi_ij
  # beta_4i = beta_4+gamma_1*b_1i+gamma_2*b_2i+gamma_3*b_3i+gamma_4*b_4i+gamma_5*a_1i+gamma_6*a_2i+tau_4i
  # beta_ji = beta_j + tau_ji, j=1,2,3
  # tau_i ~ N(0, G), xi_ij ~ N(0, omega^2)
  data = list(dataFile = paste0('after_ART_cov.csv'),
              headerTypes =c("id","ignore","cens","observation","time","limit","contcov","contcov","contcov","contcov","contcov","contcov"),
              observationTypes = list(RE="continuous"))
  modelFile = paste0('after_ART_model.txt')
  newProject(data = data, modelFile =modelFile)
  setIndividualParameterDistribution(beta1="normal",beta2="normal",beta3="normal",beta4="normal")
  # covariates on beta4
  setCovariateModel(beta4 = c(a1_est=TRUE,a2_est=TRUE,b1_est=TRUE,b2_est=TRUE,b3_est=TRUE,b4_est = TRUE))
  setErrorModel(log10_ = "constant")
  # run the estimation
  scenario <- getScenario()
  scenario$tasks = c(populationParameterEstimation = T,
                     conditionalModeEstimation = T,
                     conditionalDistributionSampling = T,
                     standardErrorEstimation=T,
                     logLikelihoodEstimation=T)
  scenario$linearization = FALSE
  setScenario(scenario)
  runScenario()
  # Find population parameter estimates
  popestimates_after <- getEstimatedPopulationParameters()
  beta1 <- popestimates_after[1]
  beta2 <- popestimates_after[2]
  beta3 <- popestimates_after[3]
  beta4 <- popestimates_after[4]
  gamma_5 <- popestimates_after[5]
  gamma_6 <- popestimates_after[6]
  gamma_1 <- popestimates_after[7]
  gamma_2 <- popestimates_after[8]
  gamma_3 <- popestimates_after[9]
  gamma_4 <- popestimates_after[10]
  # Find omega
  omega <- popestimates_after[15]
  # Find G
  taui <- as.data.frame(getEstimatedRandomEffects()$conditionalMean)
  G <- as.matrix(cov(taui[,2:5]))
  # Find s.e.
  se_after <- unlist(getEstimatedStandardErrors())
  
  
  ######################## Bootstrap
  # Create a empty dataset to store all bootstrap values
  names<-c("P1","b1","P2","b2","alpha1","alpha2","beta1","beta2","beta3","beta4","g5","g6","g1","g2","g3","g4")
  bootstrap<-as.data.frame(matrix(NA,nrow=length(names),ncol=num_bootstrap))
  rownames(bootstrap)<-names
  # Create a summary table, first column is population estimates and second column is original s.e.
  popestimates<-c(popestimates_before[1:4],alpha1,alpha2,popestimates_after[1:10])
  se <- c(se_before[1:4], se_cd4, se_after[1:10])
  se[which(se=="NaN")]<-NA
  se <- as.numeric(se)
  summary <- cbind(popestimates,se)
  colnames(summary)=c("Estimate","Original SE")
  rownames(summary)<-names
  # Find unique patient and number of patients
  n_before<-length(unique(before_ART$PATIENT))
  PATIENT_before<-unique(before_ART$PATIENT)
  n_after<-length(unique(after_ART$PATIENT))
  PATIENT_after<-unique(after_ART$PATIENT)
  # repeat bootstrap num_bootstrap times.
  for (i in 1:num_bootstrap){
    print(glue::glue("i = ",i))
    ############# Before ART 
    # simulate b_i
    b_sim<-rmvnorm(n_before, c(0,0,0,0), B)
    colnames(b_sim)=c("b1_sim","b3_sim","b2_sim","b4_sim")
    b_sim<-cbind(PATIENT_before,b_sim)
    colnames(b_sim)[1] <- "PATIENT"
    # add simulated b_i to before_ART.
    before_ART<-merge(before_ART,b_sim,by="PATIENT")
    # simulate e and y_ij, if y_ij < detection limit, then the y_ij = detection limit
    before_ART <- before_ART %>% 
      mutate(e = rnorm(nrow(before_ART),0,sigma)) %>% 
      mutate(simulated = log10(exp(P1+b1_sim-(P2+b3_sim)*time) + exp ((lambda1+b2_sim) -(lambda2+b4_sim)*time))+e) %>% 
      mutate(simulated = ifelse(simulated<log10(RNA_L), log10(RNA_L), simulated))
    # simulate a_i
    a_sim<-rmvnorm(n_before, c(0,0), A)
    colnames(a_sim)=c("a1_sim","a2_sim")
    a_sim<-cbind(PATIENT_before,a_sim)
    colnames(a_sim)[1] <- "PATIENT"
    # add simulated a_i to before_ART.
    before_ART<-merge(before_ART,a_sim,by="PATIENT")
    # simulate epsilon and z_ij, if z_ij < 1, then z_ij = 1
    before_ART <- before_ART %>% 
      mutate(epsilon = rnorm(nrow(before_ART),0, delta)) %>% 
      mutate(simulatedcd4 = alpha1+a1_sim+time*(alpha2+a2_sim)+epsilon) %>% 
      mutate(simulatedcd4 = ifelse(simulatedcd4<1, 1, simulatedcd4))
    # leave useful columns
    before_ART <- before_ART %>% 
      dplyr::select(PATIENT,RNA_L,censor,time,limit,simulated,simulatedcd4)
    # save before ART simulated data
    write.table(before_ART,"before_ART_sim.txt", sep = "," ,quote = FALSE,row.names = FALSE)
    # Find new estimates of P1,P2,lambda1, and lambda2 for simulated viral loads
    data = list(dataFile = paste0('before_ART_sim.txt'),
                headerTypes =c("id","ignore","cens","time","limit", "observation","ignore"),
                observationTypes = list(RE="continuous"))
    modelFile = paste0('before_ART_model.txt')
    newProject(data = data, modelFile =modelFile)
    setErrorModel(simulated = "constant")
    # set tasks in scenario
    scenario <- getScenario()
    scenario$tasks = c(populationParameterEstimation = T,
                       conditionalModeEstimation = T,
                       conditionalDistributionSampling = T,
                       standardErrorEstimation=T,
                       logLikelihoodEstimation=T)
    scenario$linearization = FALSE
    setScenario(scenario)
    setIndividualParameterDistribution(b1="normal",p1="normal",p2="normal",b2="normal")
    # run the estimation
    runScenario()
    # find population parameter estimates
    popestimates_before <- getEstimatedPopulationParameters()
    # find random effects
    bi <- as.data.frame(getEstimatedRandomEffects()$conditionalMean)
    # if P1 < P2, switch P1,P2 and lambda1,lambda2
    popestimates_tmp <- popestimates_before
    bi_tmp <- bi
    if (popestimates_before[1]<popestimates_before[3]){
      popestimates_before[1]<-popestimates_tmp[3]
      popestimates_before[3]<-popestimates_tmp[1]
      popestimates_before[2]<-popestimates_tmp[4]
      popestimates_before[4]<-popestimates_tmp[2]
      bi[,4]<-bi_tmp[,2]
      bi[,2]<-bi_tmp[,4]
      bi[,3]<-bi_tmp[,5]
      bi[,5]<-bi_tmp[,3]
    }
    # store the estimates in bootstrap dataframe
    bootstrap[1:4,i]=popestimates_before[1:4]
    # find the random effects bi
    colnames(bi)<-c("PATIENT","b1_est","b3_est","b2_est","b4_est")
    
    ############# CD4 Data
    long.data<-groupedData(simulatedcd4~time|PATIENT,data=before_ART)
    long.data$time<-as.numeric(long.data$time)
    long.data$PATIENT<-as.factor(long.data$PATIENT)
    # If the lme does not converge, then re-simulate cd4 data
    while (has_error (lme <- lme(simulatedcd4~time, data=long.data, 
                                 random = ~time|PATIENT,method="ML"))){
      a_sim<-rmvnorm(n_after, c(0,0), A)
      colnames(a_sim)=c("a1_sim","a2_sim")
      a_sim<-cbind(PATIENT_before,a_sim)
      colnames(a_sim)[1] <- "PATIENT"
      before_ART<-merge(before_ART,a_sim,by="PATIENT")
      before_ART <- before_ART %>% 
        mutate(epsilon = rnorm(nrow(before_ART),0,delta)) %>% 
        mutate(simulatedcd4 = alpha1+a1_sim+time*(alpha2+a2_sim)+epsilon) %>% 
        mutate(simulatedcd4 = ifelse(simulatedcd4<1, 1, simulatedcd4))
      before_ART <- before_ART %>% 
        dplyr::select(PATIENT,RNA_L,censor,time,limit,simulated,simulatedcd4)
      long.data<-groupedData(simulatedcd4~time|PATIENT,data=before_ART)
      long.data$time<-as.numeric(long.data$time)
      long.data$PATIENT<-as.factor(long.data$PATIENT)
    }
    # store the new estimates in bootstrap dataframe
    bootstrap[5:6,i]<-summary(lme)$coefficients$fixed
    # find the random effects ai
    ai<-ranef(lme)
    ai <- rownames_to_column(ai,"PATIENT")
    colnames(ai)[2:3]<-c("a1_est","a2_est")
    
    ############# Bootstrap before ART 
    # leave useful columns
    after_ART <- after_ART %>% 
      dplyr::select("PATIENT","RNA_L","censor","time","limit")
    # simulate tau_i
    taui<-rmvnorm(n_after, c(0,0,0,0), G)
    colnames(taui)=c("tau1","tau2","tau3","tau4")
    taui<-cbind(PATIENT_after,taui)
    colnames(taui)[1]<-"PATIENT"
    # Add random effects from previous two models and simulated taui to after ART data
    after_ART<-merge(after_ART,taui,by="PATIENT")
    after_ART <- merge(after_ART,bi,by="PATIENT")
    after_ART<-merge(after_ART,ai,by="PATIENT")
    # simulate xi_ij 
    after_ART <- after_ART %>% 
      mutate(xi = rnorm(nrow(after_ART),0,omega))
    # simulate w_ij, if w_ij < detection limit, then w_ij = detection limit.
    after_ART <- after_ART %>% 
      mutate(simulated = (beta1+tau1) * time / ( time + exp( beta2+tau2 - (beta3+tau3) * time )) + (beta4+a1_est*gamma_5+a2_est*gamma_6+b1_est*gamma_1+b2_est*gamma_2+b3_est*gamma_3+b4_est*gamma_4+tau4)+xi)
    # if w_ij < detection limit, then w_ij = detection limit.
    after_ART <- after_ART %>% 
      mutate(simulated = ifelse(simulated<log10(RNA_L),log10(RNA_L),simulated)) 
    # leave useful columns
    after_ART<-after_ART %>% 
      dplyr::select(PATIENT,RNA_L,censor,time,limit,simulated,a1_est,a2_est,b1_est,b2_est,b3_est,b4_est)
    # write simulated after ART data
    write.table(after_ART,"after_ART_sim.txt", sep = "," ,quote = FALSE,row.names = FALSE)
    # find new estimates of beta and gamma for simulated viral rebound data
    data = list(dataFile = paste0('after_ART_sim.txt'),
                headerTypes =c("id","ignore","cens","time","limit","observation","contcov","contcov","contcov","contcov","contcov","contcov"),
                observationTypes = list(RE="continuous"))
    modelFile = paste0('after_ART_model.txt')
    newProject(data = data, modelFile =modelFile)
    setIndividualParameterDistribution(beta1="normal",beta2="normal",beta3="normal",beta4="normal")
    # covariates on beta4
    setCovariateModel(beta4 = c(a1_est=TRUE,a2_est=TRUE,b1_est=TRUE,b2_est=TRUE,b3_est=TRUE,b4_est = TRUE))
    setErrorModel(simulated = "constant")
    # run the estimation
    scenario <- getScenario()
    scenario$tasks = c(populationParameterEstimation = T,
                       conditionalModeEstimation = T,
                       conditionalDistributionSampling = T,
                       standardErrorEstimation=T,
                       logLikelihoodEstimation=T)
    scenario$linearization = FALSE
    setScenario(scenario)
    runScenario()
    # store new estimates in the bootstrap dataframe
    popestimates_after <- getEstimatedPopulationParameters()
    bootstrap[7:16,i]<-popestimates_after[1:10]
  }
  
  # Save the bootstrap estimates in file "bootstrap_beta4.txt"
  write.table(t(bootstrap),glue::glue("bootstrap_beta4.txt"), sep = "," ,quote = FALSE,row.names =  FALSE)
  bootstrap<-read.table(glue::glue("bootstrap_beta4.txt"), sep = ",", header = TRUE)
  # find boostrap se
  bootstrap_se<-sapply(bootstrap, sd)
  # add bootstrap se to the summary table
  summary<-as.data.frame(cbind(summary,bootstrap_se))
  colnames(summary)[3]="Bootstrap SE"
  # find z-value
  summary$zvalue = abs(summary$Estimate/summary$`Bootstrap SE`)
  # find p-value
  summary$`p-value` = pnorm(summary$`zvalue`, lower.tail = FALSE)*2
  # reorder the rows in the summary table
  summary <- summary[c(1:10,13:16, 11:12), ]
  # save summary table in file "summary_beta4.txt")
  summary <- rownames_to_column(summary,"parameters")
  write.table(summary, glue::glue("summary_beta4.txt"), sep = "," ,quote = FALSE,row.names =  FALSE)
  View(summary)
}

######################## plot function
# filename, outliers, and num_months same as before
make_plots <- function(filename = "raw_data.csv", outliers=NULL, num_months=NULL){
  ######################### Clean raw data
  # Read raw data
  raw_data<-read.csv(filename,header = T)
  # If RNA_V = -1, indicate the observation as censored and replace RNA_V with RNA_L.
  # Take log10 to RNA_V. Convert time from days to months. 
  # Set lower boundary (limit) for the censoring interval to 0.
  # Remove rows with unobserved RNA_V and time
  raw_data <- raw_data %>% 
    mutate(censor = ifelse(RNA_V==-1, 1, 0)) %>% 
    mutate(RNA_V = ifelse(RNA_V==-1, RNA_L, RNA_V)) %>% 
    mutate(log10 = log10(RNA_V)) %>% 
    mutate(time = days_from_seroco/30) %>% 
    mutate(limit = 0) %>% 
    filter(!is.na(RNA_V) & !is.na(time))
  # Remove outliers if there is any
  if(!is.null(outliers)){
    raw_data <- raw_data %>% 
      filter(PATIENT != outliers)
  }
  # Find lowest detection limit for each patient
  raw_data <- raw_data %>%
    group_by(PATIENT) %>%
    mutate(RNA_L = min(RNA_L,na.rm = TRUE)) 
  # Split the dataset into before_ART (treatment=1) and after_ART(treatment=0 after treatment=1). 
  # Note: for after_ART, ignore all "treatment=0" before "treatment=1".
  # Change time after ART to the time since ART interruption 
  # For example: if before_ART time is c(1,2,3), after_ART time is c(3.5,4,5), then 
  #             the time since ART interruption for after_ART is c(0.5,1,2).
  before_ART <- raw_data %>% 
    filter(treatment==1) 
  after_ART <- raw_data %>% 
    group_by(PATIENT) %>% 
    mutate(ind = case_when(treatment == 1 ~ time)) %>%  
    mutate(ind = max(ind, na.rm=TRUE)) %>% 
    filter(time > ind) %>% 
    mutate(time = time-ind) %>% 
    dplyr::select(-ind)
  # Only leave first num_months months data for after_ART
  if(!is.null(num_months)){
    after_ART <- after_ART %>% 
      group_by(PATIENT) %>% 
      mutate(min = min(time)) %>% 
      filter(time <= num_months+min) %>% 
      dplyr::select(-min)
  }
  # plot raw_data (do not plot any "treatment=0" before "treatment=1")
  combined <- rbind(before_ART,after_ART)
  combined<-combined%>% 
    arrange(PATIENT,days_from_seroco)
  combined$PATIENT<-as.factor(combined$PATIENT)
  combined$censor<-as.factor(combined$censor)
  combined$treatment<-as.factor(combined$treatment)
  plot_all <- ggplot(combined, aes(x=days_from_seroco, y=log10)) + 
    geom_point(aes(fill=factor(censor)),size=2, shape=21, stroke=0) +
    geom_line(aes(group=PATIENT,color=treatment)) +
    scale_x_continuous("Day") + 
    scale_y_continuous(bquote("Viral load (in" ~ log[10]~"-scale)"))+
    scale_fill_manual(values=c("black","red"),labels = c("Observed data", "Censored data"))+
    scale_color_manual(values=c("steelblue","black"),labels = c("After ART interruption", "Before ART interruption"))+
    labs(color = "ART status", fill = "Data type")+ggtitle("Plot for all observations")+
    theme_classic()+theme(panel.grid.major=element_line(colour = "grey90"), panel.grid.minor =element_line(colour = "grey90"))+
      ggsave("plot_all.png")
  before_ART<-before_ART%>% 
    arrange(PATIENT,days_from_seroco)
  plot_before <- ggplot(before_ART, aes(x=days_from_seroco, y=log10)) + 
    geom_point(aes(fill=factor(censor)),size=2, shape=21, stroke=0) +
    geom_line(aes(group=PATIENT)) +
    scale_x_continuous("Day") + 
    scale_y_continuous(bquote("Viral load (in" ~ log[10]~"-scale)"))+
    scale_fill_manual(values=c("black","red"),labels = c("Observed data", "Censored data"))+
    labs(fill = "Data type")+ggtitle("before ART interruption")+
    theme_classic()+theme(panel.grid.major =element_line(colour = "grey90"), panel.grid.minor =element_line(colour = "grey90"))+
      suppressMessages(ggsave("plot_before.png"))
  after_ART<-after_ART%>% 
    arrange(PATIENT,days_from_seroco)
  plot_after <- after_ART %>% 
    mutate(days = time*30) %>% 
    ggplot(aes(x=days, y=log10)) + 
    geom_point(aes(fill=factor(censor)),size=2, shape=21, stroke=0) +
    geom_line(aes(group=PATIENT)) +
    scale_x_continuous("Day") + 
    scale_y_continuous(bquote("Viral load (in" ~ log[10]~"-scale)"))+
    scale_fill_manual(values=c("black","red"),labels = c("Observed data", "Censored data"))+
    labs(fill = "Data type")+ggtitle("after ART interruption")+
    theme_classic()+theme(panel.grid.major =element_line(colour = "grey90"), panel.grid.minor =element_line(colour = "grey90"))+
    suppressMessages(ggsave("plot_after.png"))
}



######################## analysis example
make_plots(filename = "raw_data.csv", outliers=c(4,71), num_months=8)
analysis_beta1(filename = "raw_data.csv", num_bootstrap=100, outliers=NULL, num_months=NULL)
analysis_beta2(filename = "raw_data.csv", num_bootstrap=100, outliers=NULL, num_months=NULL)
analysis_beta3(filename = "raw_data.csv", num_bootstrap=100, outliers=NULL, num_months=NULL)
analysis_beta4(filename = "raw_data.csv", num_bootstrap=100, outliers=NULL, num_months=NULL)


