######################################################################################
#####                                                           #####
#####                                                           #####
##### Author: Necla Kochan                                      ##### 
##### Date created:                                             #####
##### Date modified:                                            #####
######################################################################################
rm(list=ls())


# Load libraries and source codes

# doParallel 
library(doParallel)
#registerDoParallel(5)
library(doSNOW)
cl <- makeCluster(6, outfile = "")
registerDoSNOW(cl)

# NBinomNBC, QTBC, QTNBC, qNBC
library(limma)
library(edgeR)


setwd("C:/Users/Necla/Desktop/NeclaKochan/rectangularcodes") #change accordingly
source("QTBC.R")
source("AX_y.R")
source("AY_x.R")
source("H.R")

# Correlation regularization
library(corpcor)

setwd("C:/Users/Necla/Desktop/NeclaKochan/hapmap") #change accordingly
# Load data. Define default and fixed parameters
load("GEUV1_rawCounts_featurecount_Filtered_Normalised.RData")
d$samples <- d$samples[d$samples$group %in% c("CEU","YRI"), ]
levels(d$samples) <- factor(c("CEU", "YRI"))
#levels(d$samples$group) <- factor(c("CEU", "YRI"))
d$counts <- d$counts[, rownames(d$samples)]
colnames(d$counts) <- c(paste0("CEU", c(1:91)), paste0("YRI", c(1:89)))



Y <- d$counts
labels <- rep(c("CEU","YRI"), c(91,89)); types <- unique(labels)

# Define DGEList object
Y <- DGEList(counts=Y, group=labels, genes = rownames(Y))

n = c(20, 50, 100, 200, 300, 500, 1000) # number of selected genes
n = c(20)
# Default values
# Classification Error Rate for different number of genes (CER_n)
CER_n.QTBC <- c()
CER_n.localQTBC <- c() 

# Construct the top DE gene list from the data
design <- model.matrix(~labels)
Y <- estimateDisp(Y,design)

ptm <- proc.time()



for(i in 1:length(n)) {
        
        # Default values  
        C.QTBC <- list()
        covariance.QTBC <- list()
        sigma.QTBC <- list()
        
        # Default values for errors
        errors.localQTBC <- c()
        errors.QTBC <- c()
        
        fit <- glmFit(Y,design)
        fit <- glmLRT(fit,coef=2)
        toptable <- topTags(fit, n[i])
        DE <- toptable$table$genes # DE gene names
        d <- Y$counts[rownames(Y$counts) %in% DE, rownames(Y$samples)]
        y <- DGEList(d, group = labels, genes = rownames(d))
        
        nbs <- 10 # number of bootstrap
        
        for (h in 1:nbs){
                print(h)
                # Set the training set
                training.index <- c(sample(91), sample(92:180))
                tr_index <- c(training.index[28:91], training.index[119:180])
                training <- y$counts[, tr_index]
                
                # DGE List for training set
                training <- DGEList(training)
                train.labels <- y$samples[colnames(training),]$group
                type.train <- unique(train.labels)
                
                # Set validation set
                ts_index <- training.index[-tr_index]
                test <- y$counts[, ts_index]
                
                # DGE List for test set
                test <- DGEList(test)
                test.labels <- y$samples[colnames(test),]$group
                
                # number of observations in validation set
                m <- length(test.labels)
                
                # prior probabilities of each class
                prior.prob <- sapply( 1:length(type.train), function(q) {sum(train.labels==type.train[q])/ length(train.labels) } )
                
                
                # construct design matrix for proposed classifiers
                X <- model.matrix(~0+factor(train.labels))
                training <- estimateDisp(training, X)
                fit <- glmQLFit(training, X)
                
                # ZScoreNBinom transformation for QTBC and QTNBC
                zscoreddata <- edgeR::zscoreNBinom(training$counts, size = 1/training$trended.dispersion, mu = fit$fitted.values)
                
                # data sets displaying each subtype 
                C.QTBC <- lapply( 1:length(type.train), function(q){zscoreddata[ ,grepl( type.train[q], colnames( zscoreddata))] })
                
                # take the transpose of C.QTBC 
                C.transpose <- lapply(1:length(type.train), function(q){t(C.QTBC[[q]])})
                
                #covariance matrix for each class
                covariance.QTBC <- lapply( 1:length(type.train), function(q) { cov(C.transpose[[q]]) })
                
                # shrinking the covariance matrix for each class
                sigma.QTBC <- lapply( 1:length(type.train), function(q) { cov.shrink(covariance.QTBC [[q]]) })    
                
                # default values for decision making  
                
                decision.QTBC <- rep(0, m)
                
                #------ Classification using QTBC ------------- 
                for(l in 1:m){
                        # new observation to be classified 
                        Y.new <- test[,l]$counts
                        
                        result.QTBC <- QTBC(Y.new, prior.prob, fit$coefficients, training$trended.dispersion, sigma.QTBC)
                        decision.QTBC[l] <- types[which(result.QTBC==max(result.QTBC))]
                        
                }
                
                error.QTBC <- sum(!decision.QTBC==test.labels)/m
                errors.QTBC <- c(errors.QTBC, error.QTBC)
                
                
                S <- list() # cov matrices for each class
                
                #------ Classification using localQTBC ------------- 
                Rslt <- foreach(j = 1:m, .combine = cbind, .packages = c("corpcor") ) %dopar% {
                        #for(j in 1:m){
                        # new observation to be classified 
                        Y.new <-test[,j]$counts
                        mean <- exp(log(sum(Y.new)) + fit$coefficients)
                        
                        y1 <- as.numeric(edgeR::zscoreNBinom(Y.new, size = 1/training$trended.dispersion, mu = mean[,1]))
                        y2 <- as.numeric(edgeR::zscoreNBinom(Y.new, size = 1/training$trended.dispersion, mu= mean[,2]))
                        
                        S[[1]] <- S[[2]] <- matrix(0, n[i], n[i])
                        
                        for (k in 1:n[i]){
                                S[[1]][k,] <- as.numeric(lapply(1:n[i], function(s) {H(C.transpose[[1]][,k], C.transpose[[1]][,s], x=y1[k], y=y1[s])}))
                                S[[2]][k,] <- as.numeric(lapply(1:n[i], function(s) {H(C.transpose[[2]][,k], C.transpose[[2]][,s], x=y2[k], y=y2[s])}))
                                
                                
                        }
                        S[[1]][is.nan(S[[1]])] <- 1e-06; S[[2]][is.nan(S[[2]])] <- 1e-06
                        S[[1]][is.na(S[[1]])] <- 0; S[[2]][is.na(S[[2]])] <- 0
                        
                        # shrinking the covariance matrix for each class
                        # For lambda=0 the empirical correlations are recovered
                        
                        sigma.localQTBC <- lapply( 1:length(type.train), function(q) {cov.shrink(S[[q]])})    
                        
                        result.localQTBC <- QTBC(Y.new, prior.prob, fit$coefficients, training$trended.dispersion, sigma.localQTBC)
                        
                        return(types[which(result.localQTBC==max(result.localQTBC))])
                        
                }
                
                
                error.localQTBC <- sum(!Rslt==test.labels)/m
                errors.localQTBC <- c(errors.localQTBC, error.localQTBC)
                
        }
        
        CER_n.QTBC <-c(CER_n.QTBC, mean(errors.QTBC))
        CER_n.localQTBC <- c(CER_n.localQTBC, mean(errors.localQTBC)) 
        
}


proc.time() - ptm


CER_n.localQTBC
CER_n.QTBC
