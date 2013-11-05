# R script to analyze speeded acceptability data.
# Dan Parker Fall 2012

library(languageR)
library(lme4)
library(Hmisc)


# Function for error bars
sem <- function(x,...){
  sqrt(var(x,...)/length(x))
}

error.bar <- function(x, upper, lower=upper, length=0,...){
  arrows(x,upper, x, lower, angle=90, code=3, length=length, ...)
}

get.ci <- function(v) {
  return(t.test(v)$conf.int);
}

compute.stats <- function(data) {
  means <- mean(data);
  mins <- min(data);
  maxs <- max(data);
  stds <- sqrt(var(data));
  rtns <- sqrt(length(data));
  sems <- stds/rtns;
  seml <- means - sems;
  semu <- means + sems;
  cis <- get.ci(data)
  cil <- cis[1];
  ciu <- cis[2];
  stats <- data.frame(means, mins, maxs, cil, ciu, seml, semu);
  names(stats) <- c("Mean", "Min", "Max", "CIlow", "CIup", "-SEM", "+SEM");
  return(stats);
}


##############################################################################
# Process data
##############################################################################

data <- read.table("data.txt");

data$Rating <- as.character(data$Rating)
data$Rating[data$Rating %in% "Yes"] <- "1"
data$Rating[data$Rating %in% "No"] <- "0"

data$Rating <- as.numeric(data$Rating)

conditions <- data[data$Cond=="a" | data$Cond=="b" | data$Cond=="c",];

# compute means and SE
means.conditions <- with(conditions, tapply(conditions$Rating, list(conditions$Cond), mean, na.rm=T))
means.conditions <- means.conditions[1:3]*100
se_means.conditions <- with(conditions,aggregate(Rating,list(Cond),sem))
se_means.conditions <- se_means.conditions[1:3,2]*100


##############################################################################
# PLOTS 
##############################################################################

pdf('data.pdf', width=3.75, height=5.5)
means_conditions <- barplot(means.conditions, ylim=c(0,100), col=c("steelblue", "olivedrab", "maroon", alpha=50), 
                      ylab="Mean Rating", xaxt="n", main="Off-line Ratings: conditions \n n=24",
                      angle = 10+10*1, density = 50);

arrows(means_conditions, means.conditions + se_means.conditions, means_conditions, means.conditions - se_means.conditions,
       angle=90, lwd=2, length=0, col="dimgray")

dev.off()

