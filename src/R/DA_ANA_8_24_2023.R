# Analysis of DA dat
#install.packages("tidyverse")
D = read.csv('C:/Users/cowen/Documents/GitHub/DANA/src/R/DA_DATA_8_24_2023.csv')
D$SUBJECT = factor(D$SUBJECT);
D$LVHZ = paste(D$LV,D$HZ);
D$LVHZ = factor(D$LVHZ);
D$LV = factor(D$LV);
D$HZ = factor(D$HZ);
levels(D$LV)
levels(D$HZ)
# Two within variables
mod <- aov(DA ~ LV*HZ + Error(SUBJECT/(LV*HZ)), data=D)
summary(mod)
#install.packages("ez")
# As a sanity check, let's do it using another package...
library("ez")
ezANOVA(data=D, dv=.(DA), within=.(HZ,LV), wid=.(SUBJECT), detailed=TRUE)
###############################
# OK- I am thinking that this is working as I get the same answer in both.
# Only Hz is a significant main effect, not LV.
# As a result, we are officially only allowed to do post-hoc tests on differences in Hz.
# The benefit of this is that our correction will not be too severe.
# Thus, there are only 4 tests - one for each LV.
t1 = t.test(D[D$LV == "0" & D$HZ == "10","DA"],D[D$LV == "0" & D$HZ == "20","DA"],paired = T)
t2 = t.test(D[D$LV == "0.4" & D$HZ == "10","DA"],D[D$LV == "0.4" & D$HZ == "20","DA"],paired = T)
t3 = t.test(D[D$LV == "1" & D$HZ == "10","DA"],D[D$LV == "1" & D$HZ == "20","DA"],paired = T)
t4 = t.test(D[D$LV == "1.26" & D$HZ == "10","DA"],D[D$LV == "1.26" & D$HZ == "20","DA"],paired = T)
p.adjust(c(t1$p.value, t2$p.value, t3$p.value, t4$p.value), method = 'holm')

library('tidyverse')
ggplot(D, aes(x=LV, y=DA, fill=HZ, color = HZ)) + geom_boxplot(position=position_dodge(1)) + geom_jitter(aes(color = HZ),width = 0.01)


# indicates a main effect of Hz. No effect of LV - oh well.
pairwise.t.test(D$DA,D$LVHZ, p.adjust.method = "holm",  paired = TRUE)

boxplot(DA ~ LV*HZ, data = D)

# The following just confuses me.

#install.packages("multcomp")
library("multcomp")
#install.packages("nlme")
library("nlme")

lme1 = lme(DA ~ LV*HZ, data=D, random = ~1|SUBJECT)

anova(lme1)

summary(glht(lme1), test = adjusted(type = "holm"))
summary(glht(lme1, linfct=mcp(LV="Tukey"), test = adjusted(type = "holm")))
summary(glht(lme1, linfct=mcp(HZ="Tukey"), test = adjusted(type = "holm")))
