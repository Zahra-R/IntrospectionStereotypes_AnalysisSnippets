## Script name: RandomSplitHalfReliabilities.R
##
## Purpose of script: Calculate the split half reliability of an IAT with 100 randomly generated splits of the IAT trials
##
## Author: Zahra Rahmani Azad
##
## Date Created: 2020-July-5
##
## Copyright (c) Zahra Rahmani Azad 2020
## Email: zahra.rahmani@hotmail.de



####################################################################################################################################
## introspection into themselves? or random by cultural norms?
## Simulation of random permutations of predictions with d-scores
####################################################################################################################################

library(PerMallows)
library(lme4)
library(ggplot2)
library(dplyr)


#####preamble
# Add columns (reverse: if the IATscore is reverse coded, IDrank: make the ID numbers run consecutively, and accuracy) 
# that are missing to the dataframe

main_dat$contrastStimOrder<-ifelse(main_dat$type=="own" | main_dat$type=="ratio", -1, 1)

main_dat$reverse<-as.factor(main_dat$contrastStimOrder)
levels(main_dat$reverse)<-c("yes", "no")


func <- function(data)
{
  return(data.frame(Accuracy = cor(data$ZIATscore_n, data$imp)))
}

accuracy_per_ID<-ddply(main_dat, .(ID), func)
accuracy_per_ID$Accuracy<-ifelse(is.na(accuracy_per_ID$Accuracy), 0, accuracy_per_ID$Accuracy)
main_dat<-merge(main_dat, accuracy_per_ID, by='ID')

main_dat$IDrank<-rank(main_dat$ID, ties.method="min")%/%5+1


slim<-main_dat[,c( "IDrank", "reverse","Zpred_grp", "ZIATscore_n","type", "contrastStimOrder" )]
slim$ID<-slim$IDrank

## replace 86 with the number of subjects in your sample and prepare result matrices
IDs<-1:86
results<- matrix( nrow=1000, ncol=2)
res2<-vector(mode="numeric", length=1000)


set.seed(42)
for(i in 1:1000){
  #rdist is a permutation function, the method 'h' makes true derangements without fixpoints!
  shuffIDs<-rep(rdist(n = 1, perm.length = 86, dist.value = 86, "h"), each=5)
  helpframe<-data.frame(ID=shuffIDs, predictions=main_dat$Zpred_grp, type=main_dat$type)
  new<-base::merge(slim, helpframe)
  model1 <- lmer(ZIATscore_n ~-1 + predictions +(-1+Zpred_grp|ID), new)
  results[i,]<-confint(model1)[3,]
  res2[i]<-fixef(model1)
}

nrow(results)
length(res2)

#now visualize the results!

#m1 is the original, true/non-permutated model
m1 <- lme4::lmer(ZIATscore_n ~ -1 + Zpred_grp + (-1+Zpred_grp|ID), main_dat)
df <- data.frame(x =1:1001,
                 valu =c(res2, fixef(m1)),
                 L =c(results[,1], confint(m1)[3,1]),
                 U =c(results[,2], confint(m1)[3,2])
)



df$highlight <- ifelse(df$x == 1001, "highlight", "normal")
textdf <- df[df$x == 1001, ]
mycolours <- c("highlight" = "indianred3", "normal" = "grey50")

ggplot(data = df[850:1001,], aes(x = x, y = sort(valu))) +
  geom_point(size = 3, aes(colour = highlight)) +
  geom_errorbar(aes(ymax = sort(U), ymin = sort(L), color=highlight))+
  scale_color_manual("Status", values = mycolours) +
  geom_text(data = textdf, aes(x = x * 1.01, y = valu, label = "my label")) +
  theme(legend.position = "none") +
  theme()


dfsmall<-df[c(200:400,1001),]
dfsmall$x<-rank(dfsmall$x)
textdf <- dfsmall[dfsmall$x ==202, ]



########## Final line plot!!!!

ggplot(data = dfsmall, aes(x = x, y = sort(valu))) +
  geom_point(size = 3, aes(colour = highlight, shape=highlight)) +
  geom_errorbar(aes(ymax = sort(U), ymin = sort(L), color=highlight))+
  scale_shape_manual(values = c(16, 5))+
  scale_color_manual("Status", values = mycolours) +
  geom_text(data = textdf, aes(x = x * 1.04, y = valu, label = "Original\n dataset"), size=4, color="indianred3") +
  theme_bw()+ theme(legend.position = "none", axis.text.x=element_blank(),
                    axis.ticks.x=element_blank()) + xlab("") + ylab("Regression Weight for Prediction Accuracy")


quantile(res2, c(0.005,0.995, 1)) 


boxplot(df$valu)

df %>%
  ggplot( aes( y=valu, color="red" )) +
  geom_boxplot(outlier.colour = "black") + stat_boxplot(coef=3.1, out=1001)+
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("Boxplot chart") +
  xlab("")


######
df$highlight

df[1:800,] %>%
  ggplot( aes( x=valu)) +
  geom_density()+
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  geom_point(data = textdf, aes(x = valu, y=2)) +
  ggtitle("Violin chart") +
  xlab("")
####################

is_outlier <- function(x) {
  return(x < quantile(x, 0.25) - 3 * IQR(x) | x > quantile(x, 0.75) + 3* IQR(x))
}


######### Final Violin Plot!!!

df[1:1000,] %>%
  mutate(outlier = ifelse(is_outlier(valu), "Original dataset", as.numeric(NA))) %>%
  ggplot(., aes( y = valu)) +
  geom_violin(aes(x=0), color="lightblue4", fill="lightblue") + xlab("")+ theme_classic()+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(), legend.position = "none")+
  ylab("Regression weight for prediction accuracy")+   geom_text(data=textdf, aes(x=-0.01, y=0.635, label = "Original Dataset"),colour="palegreen4", na.rm = TRUE, hjust = -0.3)+
  geom_point(data = textdf, aes(x = 0, y=valu), size=4, color="limegreen")
 