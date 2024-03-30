#core code
####machine learning model####
library(gbm)
library(xgboost)
library(randomForest)
library(kernlab)
library(pROC)
library(tidyverse)
library(sampling)
library(reshape2)
library(caret)
###parameter tuning
#gbm
gbmGrid <-  expand.grid(
  interaction.depth = 3*(3:7),
  n.trees = 1000*(3:6),
  shrinkage = 0.005,
  n.minobsinnode = 3*(2:6)
)
set.seed(seed)
gbmFit <- train(
  x = x.train,
  y = y.train,
  method = 'gbm',
  trControl = fitControl, 
  verbose = FALSE,
  tuneGrid = gbmGrid,
  metric = "ROC"
)
bestTune <- gbmFit$bestTune
#xgb
xgbGrid <-  expand.grid(
  nrounds=2000*(1:3), 
  max_depth=5*(1:3), 
  eta=0.005, 
  gamma=0, 
  subsample=1, 
  colsample_bytree=0.1*(8:10), 
  min_child_weight=2^(0:4)
)
set.seed(seed)
xgbFit <- train(
  x = x.train,
  y = y.train,
  method = "xgbTree",
  trControl = fitControl, 
  verbose = FALSE,
  tuneGrid = xgbGrid,
  metric = "ROC"
)
bestTune <- xgbFit$bestTune
#svm
svmGrid <-  expand.grid(
  sigma = 2^(-10:-1),
  C = 1:10
)
set.seed(seed)
svmFit <- train(
  x = x.train,
  y = y.train,
  method = "svmRadial",
  trControl = fitControl, 
  verbose = FALSE,
  tuneGrid = svmGrid,
  metric = "ROC"
)
bestTune <- svmFit$bestTune
#rf
rfGrid <-  expand.grid(
  mtry = 2:20
)
set.seed(seed)
rfFit <- train(
  x = x.train,
  y = y.train,
  method = "rf",
  trControl = fitControl, 
  verbose = FALSE,
  tuneGrid = rfGrid,
  metric = "ROC"
)
bestTune <- rfFit$bestTune
###sampling
#train:test(4:1)
n <- count(train_data_all,train_data_all$dis==1) %>% as.data.frame()
set.seed(seed)  
train_test_id <- strata(train_data_all,stratanames = 'dis',size = c(0.2*n[2,2],0.2*n[1,2]),method = 'srswor')
train_test <- train_data_all[train_test_id$ID_unit,]
train_train <- train_data_all[-train_test_id$ID_unit,]
#different situations in spatial modelling
if (length(positivepoint$DDD)<376 & length(positivepoint$DDD)>124)  {
  set.seed(seed)
  sample0 <- sample(negativepoint$DDD,500-length(positivepoint$DDD),replace = F,prob = negativepoint$weightindex)%>%as.data.frame()
  names(sample0) <- 'DDD'
  sample1 <- left_join(sample0,allpoint,by='DDD') 
  sample1 <- rbind(sample1,positivepoint)
}
if (length(positivepoint$DDD)<125) {
  set.seed(seed)
  sample0 <- sample(negativepoint$DDD,length(positivepoint$DDD)*3,replace = F,prob = negativepoint$weightindex)%>%as.data.frame()
  names(sample0) <- 'DDD'
  sample1 <- left_join(sample0,allpoint,by='DDD') 
  sample1 <- rbind(sample1,positivepoint)
}
if (length(positivepoint$DDD)>375) {
  set.seed(seed)
  sample0=strata(allpoint,stratanames="dis",size=c(375,125),method="systematic",pik=allpoint$weightindex)
  sample1 <- getdata(allpoint,sample0)
}
###modelling
#train
set.seed(seed+i)
ML.train <- train(
  x = x.train,
  y = y.train,
  method = 'gbm',#"xgbTree","svmRadial","rf"
  trControl = fitControl,
  verbose = FALSE,
  tuneGrid = bestTune,
  metric = "ROC"
)
#prediction
predict(ML.train, newdata = pred_data,type='prob')[2]
#contribution
varImp(ML.train,scale = F)$importance[1]
#assessment
rocobj1 <- roc(dataset1$case,dataset1$pred,auc=T)
assessment<- coords(rocobj1, "best",
                    best.method ="youden",
                    ret=c("threshold", "specificity", "sensitivity", "accuracy",
                          "precision", "recall","youden","tn","tp","fn","fp"), transpose = FALSE)
assessment$auc <- rocobj1$auc
assessment$f1score <- 2*assessment$precision*assessment$recall/(assessment$precision+assessment$recall)
#cutoff
cutoff <- assessment$threshold
####entropy weight####
Rescale = function(x) {
  rng = range(x, na.rm = TRUE)
  (x - rng[1]) / (rng[2] - rng[1]) + 1/10000
}
Entropy_Weight = function(X) {
  X = lapply(X, Rescale)
  P = data.frame(lapply(X, function(x) x / sum(x)))
  e = sapply(P, function(x) sum(x * log(x)) *(-1/log(nrow(P))))
  d = 1 - e
  w = d / sum(d)
  return(w)
}
####durbin####
library(spdep)
library(spatialreg)
#data prepare
nbProvince <- read.gal('map/provincefinal.gal',override.id = T)
OBJECTID <- attr(nbProvince,'region.id')
nbProvince1 <- read.gwt2nb('map/provincefinal1.gwt',region.id = OBJECTID)
ndists <- attr(nbProvince1,'GeoDa')$dist
invdist.province <- lapply(ndists,function(x) 1/x)
lwProvince <- nb2listw(nbProvince1,glist = invdist.province)
#moran test
morancase <- moran.test(seq1$cases,lwProvince)
#GNS
gns <- sacsarlm(fm,data = seq1,lwProvince,Durbin = T)
gns_sum <- summary(gns)
#impact
imps <- impacts(gns,listw=lwProvince,R=2000)
imps_sum <- summary(imps,zstats=T,short=T)
####sequence analysis####
##distance
library(ape)
library(phangorn)
NS_fasta <- read.dna("NS.fasta",format = 'fasta')
NS_fasta_1 <- phyDat(NS_fasta)
NS_dist <- dist.ml(NS_fasta_1)
##similarity
library(Biostrings)
pair_wise <- pairwiseAlignment(sequence1, sequence2)
similarity <- pid(pair_wise)
####two-stage time series analysis####
rm(list=ls())
library(tidyverse)
library(tsModel)
library(dlnm)
library(splines)
library(mvmeta)
library(FluMoDL)
library(lubridate)
library(this.path)
setwd(dirname(this.path()))
## Data preparation
## The holidays
## Spring festival holiday
sf1 = ymd("2015-02-18"):ymd("2015-02-24")
sf2 = ymd("2016-02-07"):ymd("2016-02-13")
sf3 = ymd("2017-01-27"):ymd("2017-02-02")
sf4 = ymd("2018-02-15"):ymd("2018-02-21")
sf5 = ymd("2019-02-04"):ymd("2019-02-10")
sf = as.Date(c(sf1,sf2,sf3,sf4,sf5),origin = "1970-01-01")
## National day holiday
nd1 = ymd("2015-10-01"):ymd("2015-10-07")
nd2 = ymd("2016-10-01"):ymd("2016-10-07")
nd3 = ymd("2017-10-01"):ymd("2017-10-07")
nd4 = ymd("2018-10-01"):ymd("2018-10-07")
nd5 = ymd("2019-10-01"):ymd("2019-10-07")
nd = as.Date(c(nd1,nd2,nd3,nd4,nd5),origin = "1970-01-01")
holy = as.Date(c(nd,sf),origin = "1970-01-01")
runSum <- function(v,
                   lags = 0,
                   na.rm = F) {
  lagMat <- Lag(v, lags)
  rowSums(lagMat, na.rm = na.rm)
}
#dow represents the day of week; 
#tem represents the daily average temperature; 
#rhu represents the daily average relative humidity;
#total represents the daily reported influenza case number in this city
d0 <- rawdata %>%
  left_join(city_info %>% select(Citycode, Population)) %>% 
  group_by(Citycode) %>% 
  mutate(tem1 = runSum(tem,1:7)/7,
         rhu1 = runSum(rhu,1:7)/7) %>% 
  mutate(tt = as.numeric(as.factor(Date))) %>% 
  ungroup() %>% 
  mutate(holiday = case_when(
    Date %in% sf ~ "spring festival holiday",
    Date %in% nd ~ "national day holiday",
    !(Date %in% holy) ~ "common days"
  ))
xdata <- d0 %>%
  group_by(Citycode) %>%
  mutate(ac = runSum(total,1:3)+0.00000000000001) %>% 
  ungroup() %>% 
  mutate(CO = 1000*CO) %>% 
  mutate_at(vars(`PM2.5`:`CO`),function(x)x/10) 
## Choose the very air pollutant
xdata <- xdata %>%
  mutate(plu = `CO`) 
func1 = function(){
  subdata <- lapply(city, function(x)
    xdata[xdata$Citycode == x,])
  coef.plu <-
    matrix(NA, m, 1, dimnames = list(city, paste0("x", seq(1))))
  vcov.plu <- vector("list", m)
  lag = c(1,7)
  arglag = list(fun = "ns", df = 3)
  argvar <- list(fun = "lin")
  fx = as.formula(
    "total~offset(log(Population))+ns(tt,df=35)+log(ac)+factor(dow)+factor(holiday)+ns(tem1,3)+ns(rhu1,3)+cb.plu"
  )
  for (i in 1:m) {
    data <- subdata[[i]]
    cb.plu <-
      crossbasis(data$plu,
                 lag = lag,
                 argvar = argvar,
                 arglag = arglag)
    model = glm(formula = fx,
                family = quasipoisson(link = "log"),
                data = data)
    red.plu <- crossreduce(cb.plu, model)
    coef.plu[i, ] <- coef(red.plu)
    vcov.plu[[i]] <- vcov(red.plu)
  }
  ## calculate the pooled estimate cumulative rr
  mv.plu <- mvmeta(coef.plu, vcov.plu, method = "reml")
  rr_est = exp(coef(summary(mv.plu))[1]) %>% round(., 4) 
  rr_lci = exp(coef(summary(mv.plu))[5]) %>% round(., 4)  
  rr_uci = exp(coef(summary(mv.plu))[6]) %>% round(., 4)  
  cum_rr = paste0(rr_est, " (",rr_lci, "¡ª", rr_uci, ")")
  ## Calculation for the attributable fractions with pooled effects
  af <- matrix(NA, m, 3, dimnames = list(city, c("mean","lci","uci")))
  for (j in 1:m) {
    data <- subdata[[j]]
    cb.plu <-
      crossbasis(data$plu,
                 lag = lag,
                 argvar = argvar,
                 arglag = arglag)
    af[j,"mean"] = attrdl(data$plu, cb.plu, data$total, coef = coef(mv.plu),type = "af",dir = "forw", cen=0) 
    af[j,"lci"] = attrdl(data$plu, cb.plu, data$total, coef = coef(mv.plu)-1.96*sqrt(vcov(mv.plu)[[1]]),type = "af",dir = "forw", cen=0) 
    af[j,"uci"] = attrdl(data$plu, cb.plu, data$total, coef = coef(mv.plu)+1.96*sqrt(vcov(mv.plu)[[1]]),type = "af",dir = "forw", cen=0)
  }
  ## Lag-specific effect estimate
  rr_lag = data.frame()
  for(kkk in 1:7){
    for(i in 1:m){
      data <- subdata[[i]]
      cb.plu <- crossbasis(data$plu,lag=lag,argvar=argvar,arglag=arglag)
      model = glm(
        formula = fx,
        family = quasipoisson(link = "log"),
        data = data
      )
      red.plu <- crossreduce(cb.plu,model)
      red.plu =  crossreduce(cb.plu,model,type = "lag",value = kkk,cen=0)
      coef.plu[i,] <- coef(red.plu)
      vcov.plu[[i]] <- vcov(red.plu)
    }
    mv.plu <- mvmeta(coef.plu,vcov.plu,method = "reml")
    rr_est = exp(coef(summary(mv.plu))[1])
    rr_lci = exp(coef(summary(mv.plu))[5])
    rr_uci = exp(coef(summary(mv.plu))[6])
    ff = c("rr_est" = rr_est, "rr_lci" = rr_lci,"rr_uci" = rr_uci,"lag" = kkk)
    rr_lag = rbind(rr_lag,ff)
  }
  colnames(rr_lag) = c("rr_est","rr_lci","rr_uci","lag")
  plot_lag <- rr_lag %>% ggplot()+
    geom_hline(aes(yintercept=1),size=1)+
    geom_segment(aes(x = factor(lag),y = rr_lci,xend = factor(lag),yend = rr_uci))+
    geom_line(aes(x = lag,y = rr_est)) +
    geom_point(aes(x = factor(lag),y = rr_est)) +
    theme_classic()+
    xlab("Lag days")+
    ylab("RR")
  result = list(
    "plot_lag" = plot_lag,
    "rr_lag" = rr_lag,
    "cum_rr" = cum_rr,
    "af" = af)
  return(result)
}
ff = func1()
a6 = ff$rr_lag %>% mutate(plu = "CO")
xa=xxxx %>%
  mutate(rr_est = round(rr_est,4),
         rr_lci = round(rr_lci,4),
         rr_uci = round(rr_uci,4)) %>% 
  mutate(rr = paste0(rr_est," (",rr_lci,"-",rr_uci,")")) %>% 
  select(rr,plu,lag) %>% spread(plu,rr)
## Multi pollutant model
func1 = function(xitem = "pm25"){
  subdata <- lapply(city, function(x)
    xdata[xdata$Citycode == x,])
  coef.plu <-
    matrix(NA, m, 1, dimnames = list(city, paste0("x", seq(1))))
  vcov.plu <- vector("list", m)
  lag = c(1,7)
  arglag = list(fun = "ns", df = 3)
  argvar <- list(fun = "lin")
  item = paste(c("cb.pm25","cb.pm10","cb.no2","cb.so2","cb.co","cb.o3"),collapse = "+")
  fx = paste0(
    "total~offset(log(Population))+ns(tt,df=35)+log(ac)+factor(dow)+factor(holiday)+ns(tem1,3)+ns(rhu1,3)+",item
  )
  for (i in 1:m) {
    data <- subdata[[i]]
    cb.pm25 <-
      crossbasis(data$PM2.5,
                 lag = lag,
                 argvar = argvar,
                 arglag = arglag)
    cb.pm10 <-
      crossbasis(data$PM10,
                 lag = lag,
                 argvar = argvar,
                 arglag = arglag)
    cb.no2 <-
      crossbasis(data$NO2,
                 lag = lag,
                 argvar = argvar,
                 arglag = arglag)
    cb.so2 <-
      crossbasis(data$SO2,
                 lag = lag,
                 argvar = argvar,
                 arglag = arglag)
    cb.o3 <-
      crossbasis(data$O3,
                 lag = lag,
                 argvar = argvar,
                 arglag = arglag)
    cb.co <-
      crossbasis(data$CO,
                 lag = lag,
                 argvar = argvar,
                 arglag = arglag)
    model = glm(formula = as.formula(fx),
                family = quasipoisson(link = "log"),
                data = data)
    red.pm25 <- crossreduce(cb.pm25, model)
    red.pm10 <- crossreduce(cb.pm10, model)
    red.no2 <- crossreduce(cb.no2, model)
    red.so2 <- crossreduce(cb.so2, model)
    red.co <- crossreduce(cb.co, model)
    red.o3 <- crossreduce(cb.o3, model)
    red.plu = get(paste0("red.",xitem))
    coef.plu[i, ] <- coef(red.plu)
    vcov.plu[[i]] <- vcov(red.plu)
  }
  mv.plu <- mvmeta(coef.plu, vcov.plu, method = "reml")
  rr_est = exp(coef(summary(mv.plu))[1]) %>% round(., 4) 
  rr_lci = exp(coef(summary(mv.plu))[5]) %>% round(., 4)
  rr_uci = exp(coef(summary(mv.plu))[6]) %>% round(., 4) 
  cum_rr = paste0(xitem,": ",rr_est, " (",rr_lci, "¡ª", rr_uci, ")")
  return(cum_rr)
}
func1("pm25")
func1("pm10")
func1("no2")
func1("so2")
func1("co")
func1("o3")