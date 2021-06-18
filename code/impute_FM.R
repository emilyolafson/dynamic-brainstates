library(ggplot2)
library(tidyverse)
library(sjPlot)
library(latex2exp)
install.packages("gridExtra")
library(gridExtra)


#theme_set(theme_sjplot(base_size = 12, base_family = ""))

mydata = read.csv('/Users/emilyolafson/GIT/dynamic-brainstates/data/datatable_k5.csv')
summary(mydata)

#set dominant affected variable to categorical
mydata$Dom_affected_factor <- as.factor(mydata$Dom_Affected_Factor)

fit <- lm(data=mydata,ChangeFM ~ ChangeDT + Dom_affected_factor + ChangeDT*Dom_affected_factor)

# summarize the results
summary(fit)

n<-plot_model(fit, type = "eff", terms = "ChangeDT", title=TeX("$\\Delta Fugl-Meyer ~ \\Delta Dwell Time + \\Delta Dwell Time * DominantAffected$$"), axis.title=c('Change in Dwell Time in State 2 (6 months -1 week)', 'Change in Fugl-Meyer scores (6 months -1 week)'), legend.title="", show.data=T)
m<-plot_model(fit, type = "pred", terms = c("ChangeDT", "Dom_affected_factor"), title=TeX("$\\Delta Fugl-Meyer ~ \\Delta Dwell Time + \\Delta Dwell Time * DominantAffected$$"), axis.title=c('Change in Dwell Time in State 2 (6 months -1 week)', 'Change in Fugl-Meyer scores (6 months -1 week)'), legend.title="", show.data=T)
m


grid.arrange(m, n, ncol = 2, heights = c(1, 1))