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


rm(list = ls())
gc()


#Load Libraries and prepare 
library(tidyverse)
library(plyr)
library(dplyr)
library(sjlabelled)
library(CTT)
mutate<-dplyr::mutate

# Load Merged Data (Caveat: Merged Data not Clean Data as you will need the order of the IAT for the reliability by order
# you will also need trial by trial data in long format)
setwd("")
load("2-MergedDataEnviron.RData")


# Load clean Dataset as well 
setwd("")
load("CleanDataFinal.RData")

# Merge new Column with the IAT order to the main dataset
main_dat$IAT_order<-main_dat$order



#Function to calculate split half reliability
func.corr <- function(data)
{
  return(data.frame(COR = spearman.brown(cor(data$D1, data$D2))))
}



# use inquisit raw data to calculate mean reaction times for compatible vs. incompatible trials 
helpdata<-inqraw%>% plyr::rename(c("subject"= "ID", "group" = "cond","trialnum" = "order", "latency" = "RT", "stimulusnumber1" = "Stim")) %>%
  select(ID,cond, order,RT,Stim,blockcode,blocknum,trialcode,correct) %>% 
  filter(grepl('pictures_words_incompatible|pictures_words_compatible',blockcode))%>%
  filter(!grepl('Cat1' ,blockcode))%>%
  mutate(type = recode(blockcode, "pictures_words_incompatible_A" = "fam",
                       "pictures_words_compatible_A" = "fam",
                       "pictures_words_incompatible_B" = "ratio",
                       "pictures_words_compatible_B" = "ratio",
                       "pictures_words_incompatible_C" = "lang",
                       "pictures_words_compatible_C" = "lang",
                       "pictures_words_incompatible_D" = "own",
                       "pictures_words_compatible_D" = "own",
                       "pictures_words_incompatible_E" = "art",
                       "pictures_words_compatible_E" = "art"),
         blk1 = recode(blockcode, "pictures_words_incompatible_A" = 2,
                       "pictures_words_compatible_A" = 1,
                       "pictures_words_incompatible_B" = 2,
                       "pictures_words_compatible_B" = 1,
                       "pictures_words_incompatible_C" = 2,
                       "pictures_words_compatible_C" = 1,
                       "pictures_words_incompatible_D" = 2,
                       "pictures_words_compatible_D" = 1,
                       "pictures_words_incompatible_E" = 2,
                       "pictures_words_compatible_E" = 1)) %>%
#  filter(ID != 11 & ID!=24 & ID!=25) %>%                 #### IDs that were excluded posthoc...
  filter(RT<1e4)%>% 
  group_by(ID,type,blk1)%>%
  mutate(meanRT_blk1 = mean(RT))%>% 
  group_by(ID,type)%>% 
  mutate(SD = sd(RT))


helpdata<-merge(helpdata, main_dat[,c("ID", "type", "cond", "IAT_order")], by=c("ID", "type", "cond"))


#Prepare result matrices
rel_dat1<-data.frame(matrix(, nrow = 5, ncol = 101))
rel_dat1[,1]<- c("art",   "fam"   ,"lang"  ,"own"   ,"ratio")

rel_dat2<-data.frame(matrix(, nrow = 5, ncol = 101))
rel_dat2[,1]<-1:5

#generate random splits and calculate mean correlation
set.seed(1)
for(i in 2:101){
  set<-sample(1:40, 20)
    newhelpdata<- helpdata%>% group_by(ID,blk1) %>%
      mutate(part      = ifelse(order %in%  set,"1","2"))%>%
      group_by(ID,type,blk1,part)  %>%
      mutate(meanRT_blk2 = mean(RT)) %>%
      group_by(ID,type,part)  %>%
      mutate(SD_2 = sd(RT))%>%
      group_by(ID,type,blk1)%>%
      select(ID,cond,IAT_order,type,blk1,part, 
             meanRT_blk1,SD,
             meanRT_blk2,SD_2) %>%
      distinct()%>%  pivot_wider(id_cols     = c("ID", "cond", "type", "IAT_order"),
                                           names_from  = c("blk1","part"),
                                           values_from = c("meanRT_blk1", "SD","meanRT_blk2","SD_2"))%>%
      group_by(ID,type)%>%
      mutate(IATscore = (meanRT_blk1_2_1 - meanRT_blk1_1_1)/SD_1_1,
             D1 = (meanRT_blk2_2_1 - meanRT_blk2_1_1)/SD_2_2_1,
             D2 = (meanRT_blk2_2_2 - meanRT_blk2_1_2)/SD_2_2_2)
    
    rel_dat1[,i]<-ddply(newhelpdata, .(type), func.corr)[,2]
    
    rel_dat2[,i]<-ddply(newhelpdata, .(IAT_order), func.corr)[,2]
  
}


rel_dat1
rel_dat2
meanType<-data.frame(type=rel_dat1$X1, reliability=rel_dat1%>%select(X2:X101)%>%rowMeans())
meanOrder<-data.frame(order=rel_dat2$X1, reliabiity=rel_dat2%>%select(X2:X101)%>%rowMeans())
meanType
meanOrder
  