# R script for analyzing acceptability judgment data
# Dan Parker Fall 2012

library(languageR)
library(lme4)
library(Hmisc)


# Function for error bars
sem <- function(x,...){
  sqrt(var(x,...)/length(x))
}

error.bar <- function(x, upper, lower=upper, length=0.025,...){
  arrows(x,upper, x, lower, angle=90, code=3, length=length, ...)
}

data <- read.table("data.txt");

conditions <- data[data$Cond=="a" | data$Cond=="b" | data$Cond=="c",];

# Calculate means and SE
means.conditions <- with(conditions, tapply(conditions$Rating, list(conditions$Cond), mean, na.rm=T))
means.conditions<- means.conditions[1:3]
se_means.conditions <- with(conditions,aggregate(Rating,list(Cond),sem))
se_means.conditions <- se_means.conditions[1:3,2]


# Plotting

pdf('data.pdf', width=3.75, height=5.5)
means_conditions <- barplot(means.conditions, ylim=c(0,7), col=c("steelblue", "olivedrab", "maroon", alpha=50), 
                      ylab="Mean Rating", xaxt="n", main="Off-line Ratings: conditions \n n=24",
                      angle = 10+10*1, density = 50);

arrows(means_conditions, means.conditions + se_means.conditions, means_conditions, means.conditions - se_means.conditions,
       angle=90, lwd=2, length=0, col="dimgray")

dev.off()




