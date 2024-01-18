# -------------------
#   Title: Longevity & protein
#   Author: xjliu@stanford.edu 
#   Date: 8-17-2022
# ---------------------

library(dplyr)
library(tidyr)
library(haven)
library(survival)
library(openxlsx)
library(qqman)
library(pROC)

setwd("D:/RESEARCH/Proteomics")


# read data ---------------------------------------------------------------


# read proteins data
proteins <- read.csv("./dataset/CHS_Values_5K_ANML-Norm_2022_04.csv",as.is = T,check.names = F) # 3171 4986

# read marker information
marker <- read.csv("./dataset/CHS_prot_30Mar2022_5K7K_ANML_fullinfo.csv")
marker <- marker[,c(1,2,5,8,7,6)]
colnames(marker) <- c("SeqId","AptName","SomaId","UniProt","Target","TargetFullName")

# read CHS data
pheno <-  read.csv('./longevity/data/primary_model1x.csv',header=T,as.is=T) # 3067
pheno <- pheno[is.element(pheno$idno,proteins$idno),] # 3067

covariates <- c('age5','male','white','stht13','weight5','waist5','curSmok','pkyrs')


# define function ---------------------------------------------------------


doglm <- function(covar,data,aptdata,model,ttevent,outcome,...){
  aptlist <- marker$SeqId
  result <- c()
  # loop for each aptamer
  for (x in aptlist){
    cdf = data %>% inner_join(aptdata %>% dplyr::select(all_of(x),idno),by="idno")
    cdf = cdf[!is.na(cdf[x]),]
    NM = length(unique(cdf$idno))
    
    # log2 + standardize
    cdf$x1 <- scale(log2(cdf[x]),center=T,scale=T)
    
    if(model=="binary"){
      # glm model
      compformula <-  as.formula(paste(outcome, paste(c(covar,"x1"), collapse=" + "), sep=" ~ "))
      fit <- glm(compformula, data=cdf,na.action = na.omit,family = "binomial")
      tbl <- summary(fit)$coef
      result <- rbind(result,c(x, tbl[nrow(tbl),1], tbl[nrow(tbl),2], tbl[nrow(tbl),4],NM))
    }
    if(model=="survival"){
      survformula <-  paste('Surv(',ttevent,',',outcome,')')
      compformula <-  as.formula(paste(survformula, paste(c(covariates,"x1"), collapse=" + "), sep=" ~ "))
      tbl <- coef(summary(coxph(compformula,data=cdf)))
      result <- rbind(result,c(x, tbl[nrow(tbl),1], tbl[nrow(tbl),3], tbl[nrow(tbl),5],NM))
    }
    cat("Currently in pheno #", which(aptlist==x), "\n")
  }
  # organize result
  result <- as.data.frame(result)
  colnames(result) <- c("SeqId","Beta","SE","P_value","N")
  result$threshold <-  as.factor(as.numeric(result$P_value)*2505 < 0.05) # 2505 is the #PC explained 99% variation
  result <-  marker %>% right_join(result,by="SeqId") %>% 
    mutate(id = row_number()) %>% relocate(id,.before=SeqId)
  return(result)
}


# run function ------------------------------------------------------------
# logistic regression for survival to 90 outcome
binary1 <- doglm(covariates,pheno,proteins,model="binary",outcome = 'age90')
# Cox regression for overall survival
survival1 <- doglm(covariates,pheno,proteins,model="survival",ttevent = 'ttodth',outcome = 'dead' )

write.csv(binary1,file='./longevity/Routput/CHS_Norm_age90.csv',row.names=F,quote=T) 
write.csv(survival1,file='./longevity/Routput/CHS_Norm_surv_all.csv',row.names=F,quote=T)



# LASSO AUC ----------------------------------------------------------------
# read in lasso results
lassoRes <- read.csv('./longevity/Routput/lassoRes.csv')
sigProtlist <- lassoRes$SeqId # aptamer selected by LASSO 

proteins <- proteins[match(pheno$idno,proteins$idno),]
stopifnot(all(is.element(sigProtlist,colnames(proteins))))
sigProt <- proteins[,is.element(colnames(proteins),sigProtlist)] # only keep LASSO selected proteins

ln.std <- function(x) {
  x <- log(x, 2)
  x <- (x - mean(x))/sd(x)
  return(x)
}

sigProtStd <- apply(sigProt,2,ln.std) # log + standardize
colnames(sigProtStd) <- lassoRes$AptName

# test the predictive performance of LASSO selected aptamers for survival to 90
pheno <- cbind(pheno,sigProtMg) 
dim(pheno)

set.seed(5302)
default_idx = sample(nrow(pheno), nrow(pheno)/2) # 50:50 split
trn = pheno[default_idx, ] 
tst = pheno[-default_idx, ]

# protein model
form <- as.formula(paste('age90 ~',paste(c(covariates,colnames(sigProtStd)),collapse='+')))
modelFit0 <- glm(form, data=trn, family='binomial')
test_prob0 <- predict(modelFit0, newdata=tst, type='response')
test_roc0 <- plot(roc(tst$age90 ~ test_prob0,ci=T),col='darkorange',print.auc=T)

# basic model (only covariates)
form <- as.formula(paste('age90 ~',paste(covariates,collapse='+')))
modelFit1 <- glm(form, data=trn, family='binomial')
test_prob1 <- predict(modelFit1, newdata=tst, type='response')
test_roc1 <- plot(roc(tst$age90 ~ test_prob1,ci=T), plot=T, col='darkblue',print.auc=T)

# output AUC plot
tiff('longevity/Routput/fig3.auc.jpg',width = 7, height = 7,units='in',res=300, compression = c("lzw"))
test_roc0 <- plot(roc(tst$age90 ~ test_prob0,ci=T),col='darkorange')#print.auc=T ,
test_roc1 <- plot(roc(tst$age90 ~ test_prob1,ci=T), plot=T, add=T,col='darkblue')
legend('bottomright', legend=c( "Basic model (AUC=0.658)","Protein model (AUC=0.748)"),
       col=c( "darkblue","darkorange"), lty=1,cex=1,bty='n')
dev.off()



