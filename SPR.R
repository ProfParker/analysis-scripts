# R script for analyzing Linger self-paced reading data
# Dan Parker Fall 2012

library(languageR)
library(lme4)
library(plotrix)

MAX.SD <- 2.5

round.to.nearest <- function(x, p) {
  d <- as.integer(x/p);
  if (x %% p != 0) {
  	s <- sign(d);
		a <- abs(d) + 1;
		d <- s*a;
	}
	return(d*p);
}

get.ci <- function(v) {
  return(t.test(v)$conf.int);
}

error.bar <- function(x, upper, lower=upper, length=0.03, lwd=2,...){
  arrows(x,upper, x, lower, angle=90, code=3, length=length, lwd=lwd, ...)
}

trim.outliers <- function(data, by) {
  f <- rep("", nrow(data))
  for (v in by) {
    f <- paste(f, as.character(v))
  }
  f <- as.factor(f)
  means <- tapply(data$RT, f, mean);
  stds <- sqrt(tapply(data$RT, f, var));
  mins <- means - MAX.SD*(stds);
  mins.full <- mins[as.character(f)];
  maxs <- means + MAX.SD*(stds);
  maxs.full <- maxs[as.character(f)];
  data.cor <- data[data$RT >= mins.full & data$RT <= maxs.full,];
  return(data.cor);
}

compute.stats <- function(data) {
  means <- tapply(data$RT, data$Region [drop=TRUE], mean);
  mins <- tapply(data$RT, data$Region [drop=TRUE], min);
  maxs <- tapply(data$RT, data$Region [drop=TRUE], max);
  stds <- sqrt(tapply(data$RT, data$Region [drop=TRUE], var));
  rtns <- sqrt(tapply(data$RT, data$Region [drop=TRUE], length));
  sems <- stds/rtns;
  seml <- means - sems;
  semu <- means + sems;
  cis <- tapply(data$RT, data$Region [drop=TRUE], get.ci)
  cil <- sapply(cis, function(v) {v[1]});
  ciu <- sapply(cis, function(v) {v[2]});
  stats <- data.frame(means, mins, maxs, cil, ciu, seml, semu);
  names(stats) <- c("Mean", "Min", "Max", "95%tlow", "95%tupp", "-SEM", "+SEM");
  return(stats);
}

pvals.fnc <- function (object, nsim = 10000, ndigits = 4, withMCMC = FALSE, 
    addPlot = TRUE, ...) 
{
    summary <- getMethod("summary", "mer")
    require("lme4", quietly = TRUE, character = TRUE)
    if (is(object, "mer")) {
        coefs = summary(object)@coefs
        ncoef = length(coefs[, 1])
        sgma = summary(object)@sigma
        if (nsim > 0) {
            if (colnames(coefs)[3] == "z value") {
                stop("mcmc sampling is not yet implemented for generalized mixed models\n")
            }
            mcmc = try(lme4::mcmcsamp(object, n = nsim), silent = TRUE)
            if (is(mcmc, "try-error")) {
                stop("MCMC sampling is not yet implemented in lme4_0.999375\n  for models with random correlation parameters\n")
            }
            hpd = lme4::HPDinterval(mcmc)
            mcmcfixef = t(mcmc@fixef)
            nr <- nrow(mcmcfixef)
            prop <- colSums(mcmcfixef > 0)/nr
            ans <- 2 * pmax(0.5/nr, pmin(prop, 1 - prop))
            fixed = data.frame(Estimate = round(as.numeric(coefs[, 
                1]), ndigits), MCMCmean = round(apply(t(mcmc@fixef), 
                2, mean), ndigits), HPD95lower = round(hpd$fixef[, 
                1], ndigits), HPD95upper = round(hpd$fixef[, 
                2], ndigits), pMCMC = round(ans, ndigits), pT = round(2 * 
                (1 - pt(abs(coefs[, 3]), nrow(object@frame) - 
                  ncoef)), ndigits), row.names = names(coefs[, 
                1]))
            colnames(fixed)[ncol(fixed)] = "Pr(>|t|)"
            ranefNames = names(object@flist)
            assigned = attr(object@flist, "assign")
            n = length(assigned) + 1
            dfr = data.frame(Groups = rep("", n), Name = rep("", 
                n), Std.Dev. = rep(0, n), MCMCmedian = rep(0, 
                n), MCMCmean = rep(0, n), HPD95lower = rep(0, 
                n), HPD95upper = rep(0, n))
            dfr$Groups = as.character(dfr$Groups)
            dfr$Name = as.character(dfr$Name)
            for (i in 1:length(object@ST)) {
                dfr$Groups[i] = ranefNames[assigned[i]]
                dfr$Name[i] = colnames(object@ST[[i]])
                dfr$Std.Dev.[i] = round(object@ST[[i]] * sgma, 
                  ndigits)
                dfr$MCMCmedian[i] = round(median(mcmc@ST[i, ] * 
                  mcmc@sigma), ndigits)
                dfr$MCMCmean[i] = round(mean(mcmc@ST[i, ] * mcmc@sigma), 
                  ndigits)
                hpdint = as.numeric(lme4::HPDinterval(mcmc@ST[i, 
                  ] * mcmc@sigma))
                dfr$HPD95lower[i] = round(hpdint[1], ndigits)
                dfr$HPD95upper[i] = round(hpdint[2], ndigits)
            }
            dfr[n, 1] = "Residual"
            dfr[n, 2] = " "
            dfr[n, 3] = round(sgma, ndigits)
            dfr[n, 4] = round(median(mcmc@sigma), ndigits)
            dfr[n, 5] = round(mean(mcmc@sigma), ndigits)
            hpdint = as.numeric(lme4::HPDinterval(mcmc@sigma))
            dfr[n, 6] = round(hpdint[1], ndigits)
            dfr[n, 7] = round(hpdint[2], ndigits)
            mcmcM = as.matrix(mcmc)
            k = 0
            for (j in (ncol(mcmcM) - n + 1):(ncol(mcmcM) - 1)) {
                k = k + 1
                mcmcM[, j] = mcmcM[, j] * mcmcM[, "sigma"]
                colnames(mcmcM)[j] = paste(dfr$Group[k], dfr$Name[k], 
                  sep = " ")
            }
            if (addPlot) {
                m = data.frame(Value = mcmcM[, 1], Pblueictor = rep(colnames(mcmcM)[1], 
                  nrow(mcmcM)))
                for (i in 2:ncol(mcmcM)) {
                  mtmp = data.frame(Value = mcmcM[, i], Pblueictor = rep(colnames(mcmcM)[i], 
                    nrow(mcmcM)))
                  m = rbind(m, mtmp)
                }
                print(densityplot(~Value | Pblueictor, data = m, 
                  scales = list(relation = "free"), par.strip.text = list(cex = 0.75), 
                  xlab = "Posterior Values", ylab = "Density", 
                  pch = "."))
            }
            if (withMCMC) {
                return(list(fixed = format(fixed, digits = ndigits, 
                  sci = FALSE), random = dfr, mcmc = as.data.frame(mcmcM)))
            }
            else {
                return(list(fixed = format(fixed, digits = ndigits, 
                  sci = FALSE), random = dfr))
            }
        }
        else {
            coefs = summary(object)@coefs
            ncoef = length(coefs[, 1])
            fixed = data.frame(Estimate = round(as.numeric(coefs[, 
                1]), ndigits), pT = round(2 * (1 - pt(abs(coefs[, 
                3]), nrow(object@frame) - ncoef)), ndigits), 
                row.names = names(coefs[, 1]))
            colnames(fixed)[ncol(fixed)] = "Pr(>|t|)"
            return(list(fixed = format(fixed, digits = ndigits, 
                sci = FALSE)))
        }
    }
    else {
        cat("the input model is not a mer object\n")
        return()
    }
}


##############################################################################
# Begin processing
##############################################################################

data <- read.table("npi.data.txt", quote=" ")
names(data) <- c("Subject", "Trial_type", "Item_number", "Condition",
                 "Position", "Word", "Region", "RT")



#Select my items
data <- data[data$Trial_type == "<your experiment name>", ]
data$Trial_type <- factor(data$Trial_type)


##############################################################################
# Accuracy
##############################################################################

data.q <- data[data$Position == "?",]   
data.q$Accuracy <- as.numeric(as.character(data.q$Region))
accuracy.subj <- tapply(data.q$Accuracy, data.q$Subject, mean)									

include <- accuracy.subj >= 0.8

subjects <- names(accuracy.subj)
include.subj <- subjects[include]
data <- data[data$Subject %in% include.subj, ]
data$Subject <- factor(data$Subject)

print( "Overall Accuracy")
print (mean(as.numeric(as.character(data[data$Position == "?", "Region"]))))

data <- data[data$Position != "?", ]

# By condition:
# mean(data.q$Accuracy[data.q$Condition=="a"])


##############################################################################
# Regions of interest
##############################################################################


reg.ever <- c("01", "02", "03", "04", "05", "07", "08",
             "09", "10", "12","13","14", "15", "16")


##############################################################################
# Trimming the data
##############################################################################

data <- data[data$RT < 2000, ]     
data <- trim.outliers(data, list(data$Region, data$Condition))

##############################################################################
# Region averages
##############################################################################

data.a <- data[data$Condition == "a", ]  
data.b <- data[data$Condition == "b", ]
data.c <- data[data$Condition == "c", ]

data.a.stat <- compute.stats(data.a)
data.b.stat <- compute.stats(data.b)
data.c.stat <- compute.stats(data.c)


##############################################################################
# Plot: data
##############################################################################
pdf ("data.pdf", width=16, height = 4)

min.y <- min(c(data.a.stat[["-SEM"]], data.b.stat[["-SEM"]], data.c.stat[["-SEM"]]));
max.y <- max(c(data.a.stat[["+SEM"]], data.b.stat[["+SEM"]], data.c.stat[["+SEM"]]));

xlimit <- c(1,length(reg.data));
ylimit <- c(min.y,max.y);

par(mar=c(5,5,1.75,5))
plot(NA,xlim=xlimit,ylim=ylimit, xlab="Region", ylab="Reading time (ms)", 
     xaxt="n", yaxt="n", main="", type="n", cex.lab=2);            

legend(.5, 720, c("A","B","C"), 
       pch=c(15, 17, 0), col=c("steelblue", "olivedrab", "maroon"), cex=1.25, pt.cex=1.75, bty="n")

axis(1, at=seq(1, 14, 1), labels=c("1", "2", "3", "4", "5", "6", "7", "8", "9", 
                                   "10", "11", "12", "13", "14"), las=1, cex.axis=1.5);    
axis(2, at=seq(350, 650, by=100), cex.axis=1.45);    

text(9,680,expression(paste(italic("data"))), cex=1.75, col="dodgerblue4" , font=2)

error.bar (1:14, data.a.stat[reg.data,"+SEM"], data.a.stat[reg.data,"-SEM"],col= "gray40", lwd=1)
points (data.a.stat[reg.data,"Mean"], type="b", pch=15, col= "#2EB8E9", lwd=2.25, cex=1.75)

error.bar (1:14, data.b.stat[reg.data,"+SEM"], data.b.stat[reg.data,"-SEM"],col= "gray40", lwd=1)
points (data.b.stat[reg.data,"Mean"], type="b", pch=17, col= "#99BC5F", lwd=2.25, cex=1.75)

error.bar (1:14, data.c.stat[reg.data,"+SEM"], data.c.stat[reg.data,"-SEM"],col= "gray40", lwd=1)
points (data.c.stat[reg.data,"Mean"], type="b", pch=22, bg="white", col= "maroon", lwd=2.25, cex=1.75)

dev.off()



