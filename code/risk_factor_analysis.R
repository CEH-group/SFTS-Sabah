########################### Risk factor analysis ############################## 


# load packages
{
  lib_list = c("dplyr", "glmmTMB", "verification")
  inst_pkg = lib_list[!lib_list %in% installed.packages()]
  lapply(inst_pkg, function(x) install.packages(x, dependencies = TRUE))
  sapply(lib_list, require, character = TRUE)
}


# load data 
results <- read.csv("SFTS_data.csv")


# Make the data factor
n_col_cate = length(colnames(results))
factor_data_mod_cate = results %>% dplyr::select(-c(houseID, results_bi)) %>% apply(MARGIN = 2, FUN = factor)

# bind data
data_mod_cate_fac = cbind(factor_data_mod_cate, results %>% dplyr::select(c(houseID, results_bi)))

# Factor level
data_mod_cate_fac$Age_group  = data_mod_cate_fac$Age_group %>% factor(levels = c("Baby (0-4)" ,"Child (5-14)", "Adult (15-59)", "Elderly (60+)"))


data_mod_cate_fac = data_mod_cate_fac %>% mutate(Age_group_2 = 
                                         case_when(Age_group %in% c("Child (5-14)", "Adult (15-59)", "Elderly (60+)") ~ "over_5", 
                                                   Age_group == "Baby (0-4)" ~ "under_5") %>% factor(levels = c("under_5", "over_5")))

n_col_cate = length(colnames(data_mod_cate_fac))

# Null model
null <- glmmTMB(results_bi ~ +(1|houseID), data=data_mod_cate_fac, family="binomial")
summary(null)


# Univariable analysis
list_mod_cate = list(NULL)
list_anova_cate = list(NULL)
for(i in c(1:18, 21)){
  
  
  formula = paste("results_bi ~", colnames(data_mod_cate_fac)[i], "+ (1|houseID)")
  f1 = as.formula(formula)
  
  mod <- glmmTMB(f1, data=data_mod_cate_fac, family="binomial")
  
  comp_results = anova(null, mod, test = "LRT")
  list_mod_cate[[i]] = mod
  list_anova_cate[[i]] = comp_results
  
  print(i)
  
}
names(list_mod_cate) = colnames(data_mod_cate_fac)[c(1:18, 21)]

list_mod_cate$Age_group %>% summary()
list_mod_cate$Gender_factor %>% confint %>% exp() %>% round(3)
list_mod_cate$outsideKmpg %>% confint %>% exp() %>% round(3)


sig_anova_cate = list_anova_cate %>% 
  lapply(FUN = function(x){
    x$`Pr(>Chisq)`[[2]]<0.30
  }) %>% unlist


names(list_anova_cate) = colnames(data_mod_cate_fac)

list_pvalue = list_anova_cate %>% 
  lapply(FUN = function(x){
    x$`Pr(>Chisq)`[[2]]
 }) %>% unlist(use.names = TRUE) %>% round(3) 


list_anova_cate[list_pvalue<0.30]




# Forward selection with Age_group
colnames(data_mod_cate_fac)

mod1 <- glmmTMB(results_bi ~ factor(Age_group_2) +(1|houseID), data=data_mod_cate_fac, family="binomial")
summary(mod1)


mod1 %>% confint() %>% exp() %>% round(3)

n_col = length(colnames(data_mod_cate_fac))

list_mod_cate2 = list(NULL)
list_anova_cate2 = list(NULL)
for(i in 1:c(n_col_cate-3)){
  
  
  formula = paste0("results_bi ~", colnames(data_mod_cate_fac)[i], "+ factor(Age_group_2) + (1|houseID)")
  f1 = as.formula(formula)
  
  mod <- glmmTMB(f1, data=data_mod_cate_fac, family="binomial")
  
  comp_results = anova( mod, mod1, test = "LRT")
  list_mod_cate2[[i]] = mod
  list_anova_cate2[[i]] = comp_results
  
  print(i)
  
}




list_pvalue_forward1 = list_anova_cate2 %>% 
  lapply(FUN = function(x){
    x$`Pr(>Chisq)`[[2]]
  }) %>% unlist(use.names = TRUE) %>% round(3) 




list_pvalue_forward1 = list_anova_cate2 %>% 
  lapply(FUN = function(x){
    x$`Pr(>Chisq)`[[2]]
  }) %>% unlist(use.names = TRUE) %>% round(3) 

list_pvalue_forward1[list_pvalue_forward1<0.05]

