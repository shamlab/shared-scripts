#' burden test
#' @param pheno vector of phenotypes
#' @param genotypeMatrix matrix of genotypes 
#' @return test statistic
burden <- function (pheno, genotypeMatrix) {
        if(0%in%pheno) {
                pheno[which(pheno==0)]  <- -1
        }
        return(sum(t(pheno)%*%genotypeMatrix)^2)
}

#' cmc test
#' @param pheno vector of phenotypes
#' @param genotypeMatrix matrix of genotypes 
#' @return test statistic
cmc <- function (pheno, genotypeMatrix) {
        if(0%in%pheno) {
                pheno[which(pheno==0)]  <- -1
        }
        temp <- rowSums(genotypeMatrix)

        nControls  <- length(which(pheno < 0))
        nCases  <- length(which(pheno > 0))

        temp <- pheno * temp
        cases <- which(temp > 0)
        controls <- which(temp < 0)
        test.stat = (length(cases) / nCases) - (length(controls) / nControls)
        return(test.stat^2)
}

#' permutation test
#' @param pheno vector of phenotypes [0,1]
#' @param genotypeMatrix matrix of genotypes 
#' @param number of iterations 
#' @param FUN rare variant test function (wither CMC or burden or another one)
#' @return list with test statistic and pvalue
permute <- function(pheno, genotypeMatrix, iter, FUN) {
        burdenStat <- FUN(pheno, genotypeMatrix)
        dis <- vector("numeric", length=iter)
        for(i in 1:iter) {
                rand <- sample(1:length(pheno), length(pheno))
                dis[i] <- FUN(pheno[rand], genotypeMatrix)
        }
        pval <- mean(dis >=burdenStat) 
        return(pval)
}


#' liability model (without LD)
#' @param stand.genotypeMatrix standardized genotype matrix
#' @param liability lifetime risk, or liability threshold 
#' @param hera effect size 
#' @param no.causal number of causal variants in %
#' @param n.cases number of cases
#' @param n.controls number of controls
#' @return a vector of row numbers of cases and controls
liability.model <- function (stand.genotypeMatrix, liability, hera, no.causal, n.cases, n.controls) {
        no.var <- ncol(stand.genotypeMatrix)
        no.causal <- floor(no.var*no.causal)
        effect <- c(rep(1, no.causal), rep(0, (no.var-no.causal)))
        if (hera==0) {
                effect <- rep(0, length(effect))
        } else {
                effect <- sqrt(effect^2/(sum(effect^2)/hera))
        }
        risk <- stand.genotypeMatrix %*% effect
        qtil <- qnorm(1-liability)
        nsim <- nrow(stand.genotypeMatrix)
        cases <- vector("numeric", n.cases)
        controls <- vector("numeric", n.controls)

        c.cases <- 1
        c.controls <- 1
        while(c.cases<=n.cases | c.controls<=n.controls) {
                liab.dist <- risk + rnorm(nsim, 0, sqrt(1-hera))
                ln <- sample(1:nsim, nsim)
                for (s in ln) {
                        if(liab.dist[s] >=qtil && c.cases<= n.cases) {
                                cases[c.cases] <- s
                                c.cases <- c.cases+1
                        }
                        if(liab.dist[s] <qtil && c.controls<= n.controls) {
                                controls[c.controls] <- s
                                c.controls <- c.controls+1
                        }
                }
        }
        genoIndex <- c(cases, controls)
        return(genoIndex)
}

#' computes the genomic infaltion factor
#' @param z a numeric vector of type stat_type
#' @param stat_type a character of either PVAL, CHISQ or Z
#' @return lambda
lambda <- function (z, stat_type="PVAL") {
if (stat_type == "Z") {
   message("using Z scores")
} else if (stat_type == "CHISQ") {
   message("using CHISQ")
   z <- sqrt(z)
} else if (stat_type == "PVAL") {
   message("using PVAL")
   z <- qnorm(z/2)
} else {
  stop("unknown type")
}
lambda  <-  round(median(z^2)/.454,3)
return(lambda)
}

