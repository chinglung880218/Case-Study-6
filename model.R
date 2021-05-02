library(tidyverse)
library(dplyr)
library(ggplot2)
library(patchwork)
library(splines)
library(rstanarm)
library(brms)
library(rstan)
library(lme4)
library(data.table)



# Data

setwd("~/STA723/Case Study 6/Case-Study-6")
ToxCast <- readRDS("ToxCast.RDS")

obese_data <- ToxCast$obese_data
hitc_data <- ToxCast$hitc_data
aenm_data <- ToxCast$aenm_data


# EDA

plot_heat <- obese_data %>% 
  group_by(aenm, code) %>% 
  count() %>% 
  # pivot_wider(names_from = code, values_from = n) %>% 
  ggplot(mapping = aes(x = aenm, y = code, fill = n)) +
  labs(x = "Assay Endpoint", y = "Chemical", fill = "Number of \nDoses") +
  ggtitle("Heat Map of Doses for Pairwise Cell")+
  geom_raster() +
  scale_fill_gradientn(colours = c("#7fcdbb", "#1d91c0", "#2c7fb8", "#225ea8","#253494", "#081d58")) +
  theme_minimal(base_size = 16) +
  theme(axis.text.x=element_blank()) +
  theme(axis.text.y=element_blank())


ttt <- obese_data %>% 
  group_by(aenm, code) %>% 
  count() %>% 
  pivot_wider(names_from = code, values_from = n)

plot_unidose <- obese_data %>% 
  group_by(aenm, code) %>% 
  summarise(unique_dose = n_distinct(logc)) %>% 
  ggplot(mapping = aes(x = unique_dose)) +
  geom_histogram(fill = "#238443") +
  labs(x = "Number of Unique Doses", y = "Number of Chemical-Assay Pairs") +
  ggtitle("Histogram of Number of Unique Doses") +
  theme_bw(base_size = 16)


obese_data %>% 
  group_by(aenm, code) %>% 
  count() %>% 
  group_by(aenm) %>%
  count() %>% 
  ggplot(mapping = aes(x = n)) +
  geom_histogram(fill = "#0570b0", binwidth = 1) +
  labs(x = "Number of Unique Chemicals", y = "Number of Assay Endpoints") +
  ggtitle("Histogram of Number of Unique Chemicals") +
  theme_bw(base_size = 16)


# dose-response model for single cell
tmp <- obese_data %>% 
  filter(aenm == "TOX21_AR_LUC_MDAKB2_Antagonist_Specificity", code == "C80057")

fit <- lm(resp ~ logc, tmp)
fit$coefficients[2]

# dose-response model for single cell
tmp2 <- obese_data %>% 
  filter(aenm == "ATG_ERE_CIS_up", code == "C80057")

lm(resp ~ logc, tmp2)

test2 <- obese_data

aenm <- unique(obese_data$aenm)
code <- unique(obese_data$code)
test <- crossing(aenm, code) %>% 
  mutate(id = row_number())

test2 <- left_join(test2,test, by = c('aenm'='aenm', 'code'='code'))


test$activity <- 0

for (i in 1:length(test2$id)){
  
  tmp <- test2 %>% 
    filter(id %in% i)
  
  if (length(tmp$id) != 0){
    
    fit <- lm(resp ~ logc, tmp)
    test$activity[i] <- fit$coefficients[2]
    
  }
  
}

data <- test %>% 
  mutate(status = case_when(
    activity < -1 ~ 1,
    activity <= 1 ~ 0,
    activity > 1 ~ 1
  ))


test3 <- left_join(obese_data,data, by = c('aenm'='aenm', 'code'='code'))



write.csv(data,"activity2.csv", row.names = FALSE)

obese_matrix <- data %>% 
  select(aenm, code, status) %>% 
  pivot_wider(names_from = code, values_from = status)

colors <- c("pink", "red")

data %>% 
  ggplot(mapping = aes(x = aenm, y = code, fill = factor(status))) +
  labs(x = "Assay Endpoint", y = "Chemical", fill = "Activities") +
  ggtitle("Heat Map of Binary Activity Classification")+
  geom_raster() +
  scale_fill_manual(values=colors) +
  theme_minimal(base_size = 16) +
  theme(axis.text.x=element_blank()) +
  theme(axis.text.y=element_blank())



 obese_matrix[is.na(ttt)] <- NA


 new_obese <- obese_matrix %>% 
   pivot_longer(cols = C106467:C95954,
                names_to = "code", values_to = "status")

new_obese %>% 
  ggplot(mapping = aes(x = aenm, y = code, fill = factor(status))) +
  labs(x = "Assay Endpoint", y = "Chemical", fill = "Activities") +
  ggtitle("Heat Map of Binary Activity Classification")+
  geom_raster() +
  scale_fill_manual(values=colors) +
  theme_minimal(base_size = 16) +
  theme(axis.text.x=element_blank()) +
  theme(axis.text.y=element_blank())



hitc_matrix <- hitc_data %>% 
  pivot_wider(names_from = code, values_from = hitc)


tttt<-merge(hitc_data, new_obese, by=c("aenm","code"))

tttt  %>% 
  ggplot(mapping = aes(x = aenm, y = code, fill = factor(hitc))) +
  labs(x = "Assay Endpoint", y = "Chemical", fill = "Activities") +
  ggtitle("Heat Map of Activity Classification")+
  geom_raster() +
  scale_fill_manual(values=colors) +
  theme_minimal(base_size = 16) +
  theme(axis.text.x=element_blank()) +
  theme(axis.text.y=element_blank())

tttt  %>% 
  ggplot(mapping = aes(x = aenm, y = code, fill = factor(status))) +
  labs(x = "Assay Endpoint", y = "Chemical", fill = "Activities") +
  ggtitle("Heat Map of Activity Classification by Linear Regression")+
  geom_raster() +
  scale_fill_manual(values=colors) +
  theme_minimal(base_size = 16) +
  theme(axis.text.x=element_blank()) +
  theme(axis.text.y=element_blank())
