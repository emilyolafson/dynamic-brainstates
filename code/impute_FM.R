library(ggplot2)
library(tidyverse)
library(sjPlot)
library(latex2exp)
install.packages("gridExtra")
library(gridExtra)


theme_set(theme_sjplot(base_size = 15, base_family = ""))

# Hypothesis 1: Increase in dwell time in state 2 associated with better motor recovery, but only in subjects whose dominant CST is disrupted.
mydata = read.csv('/Users/emilyolafson/GIT/dynamic-brainstates/data/datatable_state2_DT.csv')
summary(mydata)

#set dominant affected variable to categorical
mydata$Dom_affected_factor <- as.factor(mydata$Dom_Affected_Factor)

fit <- lm(data=mydata,ChangeFM ~ ChangeDT + Dom_affected_factor + ChangeDT*Dom_affected_factor)

# summarize the results
summary(fit)

plot_model(fit, type = "pred", terms = c("ChangeDT", "Dom_affected_factor"), colors=c( "firebrick2","darkseagreen"), title=TeX("$\\Delta Fugl-Meyer ~ \\Delta Dwell Time + \\Delta Dwell Time * DominantAffected$$"), axis.title=c('Change in Dwell Time in State 2 (6 months -1 week)', 'Change in Fugl-Meyer scores (6 months -1 week)'), legend.title="", show.data=T)



# Hypothesis 2: Increase in appearance rate in state 4 is associated with better motor recovery, but only in subjects whose non-dominant CST is disrupted
mydata = read.csv('/Users/emilyolafson/GIT/dynamic-brainstates/data/datatable_state4_FO.csv')
#set dominant affected variable to categorical
mydata$Dom_affected_factor <- as.factor(mydata$Dom_Affected_Factor)

fit <- lm(data=mydata,ChangeFM ~ ChangeAR + Dom_affected_factor + ChangeAR*Dom_affected_factor)

# summarize the results
summary(fit)
scale_color_sjplot(discrete=T,reverse = T)
plot_model(fit, type = "pred", terms = c("ChangeAR", "Dom_affected_factor"), colors = c("darkseagreen", "firebrick2"), title=TeX("$\\Delta Fugl-Meyer ~ \\Delta Fractional Occupancy + \\Delta Fractional Occupancy * DominantAffected$$"), axis.title=c('Change in Dwell Time in State 4 (3 months -1 week)', 'Change in Fugl-Meyer scores (3 months -1 week)'), legend.title="", show.data=T)



