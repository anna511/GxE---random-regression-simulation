### NRRM simulation - modified ###
### Ani A. Elias - July, 2017 ###

library(lme4) # to test the models on simulated dataset

# loc is the vector of location numbers
      loc <- c(5,10,20,50,100,150,170,190,200)
      # special case minimum is 5 location and 2 hybrids OR 3 loc and 4 hyb
# hyb is the hybrid numbers
    # element of c(2,5,10,20,40,80,160,320,640,700,...)
      # checked up to 640 hybrids which is more than three times the maximum number of hybrids available in real datasets
      # when checking 3 location keep the minimum # of hybrids as 4 
# varL is the vector of loc variance
     varL <- seq(1,200,15)
# varH is the hybrid variance
    # element of seq(1,200,15)
# corHI is the vector of correlation between loc and index
     corHI <- c(-1,-0.5,0,0.5,1)
# varI is the vector of index variance
     varI <- seq(1,200,15)
# repeating the simulation for 20 times

simulate.NRRM <- function(loc, hyb, varL, varH, corHI, varI){

for (i in 1:length(loc)){
    for(j in 1:length(varL)){
        for(k in 1:length(corHI)){
            for(l in 1:length(varI)){

                # defining output matrix
                pred.mb <- matrix(NA, 20, 23)
                colnames(pred.mb) <- c("Field.n", "Hybrid.n","Field.v", "Hybrid.v", "covHI","Index.v", "herit", "Fieldv.m", "Hybridv.m",
                                        "covHI.m", "Indexv.m", "Fieldv.b","Hybridv.b", "pCOR.Hm", "pCOR.Im", "pCOR.Fm", "pCOR.Hb", 
                                            "pCOR.Fb", "pRMSE.Hm", "pRMSE.Im", "pRMSE.Fm", "pRMSE.Hb", "pRMSE.Fb")

                for(s in 1:20){

                    L <- loc[i]
                    H <- hyb
                    var.L <- varL[j]
                    df1 <- data.frame(Field = sort(rep(c(1:L),H)), Hybrid = rep(c(1:H),L))

                    # generating location effect
                    alphaL <- rnorm(L, sd = sqrt(var.L))
                    df1$alphaL <- alphaL[df1$Field]

                    # generating hybrid effects and influence across location
                    covariance <- corHI[k] * (sqrt(varH * varI[l]))
                    cov.matrix <- matrix(c(varH, covariance, covariance, varI[l]), nrow=2, byrow=TRUE)
                    diag(cov.matrix)<-diag(cov.matrix)+1e-6
                    sqrtCov <- chol(cov.matrix)
                    random.effects <- data.frame(t(t(sqrtCov) %*% matrix(rnorm(2*H), nrow=2, ncol = H)))
                    random.effects[,3] <- as.matrix(c(1:H))
                    colnames(random.effects) <- c("alphaH", "beta", "Hybrid")

                    # generating location index values
                    Vk <- as.matrix(runif(L,-2,2))
                    Indices <-  data.frame(Hybrid = sort(rep(c(1:H), L)), Field = rep(c(1:L),H), Index =rep(Vk, H))

                    # merging random effects, adding error, and calculating phentoype
                    random.effects2 <- merge(random.effects, Indices, by="Hybrid")
                    random.effects2$betaVH <- as.matrix(random.effects2$beta * random.effects2$Index)
                    random.effects3 <- merge(random.effects2, df1, by=c("Field", "Hybrid"))

                    #heritability to start 
                    h2 <- varH/(varH + varI[l] + varL[j] + 1)

                    df3 <- within(random.effects3, phenotype <- alphaL + alphaH + betaVH + rnorm(n=L*H))

                    #testing the simulated dataset 
                    model <- lmer(phenotype ~ 1+ (1|Field) + (1+ Index|Hybrid), df3) 
                    base <- lmer(phenotype ~ 1 + (1|Field) + (1|Hybrid), df3)

                    # extracting random effects
                    fittednorms.m <-ranef(model)$Hybrid[,1:2]
                    colnames(fittednorms.m)<-c("intercept.m","slope.m")
                    fittednorms.m$Hybrid <- rownames(fittednorms.m)

                    fittednorms.b <- ranef(base)$Hybrid
                    colnames(fittednorms.b) <- c("intercept.b")
                    fittednorms.b$Hybrid <- rownames(fittednorms.b)

                    fittednorms.Lm <- ranef(model)$Field
                    colnames(fittednorms.Lm) <- c("loc.m")
                    fittednorms.Lm$Field <- rownames(fittednorms.Lm)

                    fittednorms.Lb <- ranef(base)$Field
                    colnames(fittednorms.Lb) <- c("loc.b")
                    fittednorms.Lb$Field <- rownames(fittednorms.Lb)

                    datalist.a <- list(df3, fittednorms.m, fittednorms.b)
                    df4 <- Reduce(function(x, y) { merge(x, y, by=c('Hybrid')) }, datalist.a)

                    datalist.b <- list(df4, fittednorms.Lm, fittednorms.Lb)
                    df5 <- Reduce(function(x, y) { merge(x, y, by=c('Field')) }, datalist.b)

                    # estimated sd, pCOR, and pRMSE
                    pred.mb[s,1] <- L
                    pred.mb[s,2] <- H
                    pred.mb[s,3] <- var.L
                    pred.mb[s,4] <- varH
                    pred.mb[s,5] <- covariance
                    pred.mb[s,6] <- varI[l]
                    pred.mb[s,7] <- h2
                    pred.mb[s,8] <- unlist(VarCorr(model))[5]
                    pred.mb[s,9] <- unlist(VarCorr(model))[1]
                    pred.mb[s,10] <- unlist(VarCorr(model))[2]
                    pred.mb[s,11] <- unlist(VarCorr(model))[4]
                    pred.mb[s,12] <- unlist(VarCorr(base))[2]
                    pred.mb[s,13] <- unlist(VarCorr(base))[1]
                    pred.mb[s,14] <- cor(df5$alphaH, df5$intercept.m)
                    pred.mb[s,15] <- cor(df5$beta, df5$slope.m)
                    pred.mb[s,16] <- cor(df5$alphaL, df5$loc.m)
                    pred.mb[s,17] <- cor(df5$alphaH, df5$intercept.b)
                    pred.mb[s,18] <- cor(df5$alphaL, df5$loc.b)
                    pred.mb[s,19] <- sqrt(mean((df5$alphaH - df5$intercept.m)^2))
                    pred.mb[s,20] <- sqrt(mean((df5$beta - df5$slope.m)^2))
                    pred.mb[s,21] <- sqrt(mean((df5$alphaL - df5$loc.m)^2))
                    pred.mb[s,22] <- sqrt(mean((df5$alphaH - df5$intercept.b)^2))
                    pred.mb[s,23] <- sqrt(mean((df5$alphaL - df5$loc.b)^2))

                } # end of replications 
                write.csv(pred.mb, paste0("loc-",L,"_hyb-",H,"_varL-",var.L, "_varH-",varH,"_varI-",varI[l],"_corHI-",corHI[k],"_NRRM_simout.csv"))  
            } # end of varI        
        } # end of corHI    
    } # end of varL
} # end of loc

} # end of function
