# Load required libraries ####
library(chisq.posthoc.test)
library(tidyr)
require("rstatix")
library(purrr)
library(ineq) #for calculating gini coefficients
library(vegan)
# Figure 2 ----------------------------------------------------------------
# Panel I/J
var_list_green<- list("COMB_maximum_y","COMB_startPoint_x2","COMB_midPoint1_x","COMB_slope1", "COMB_incrementTime")
var_list_red<- list("COMB_maximum_y_red","COMB_startPoint_x_red2","COMB_midPoint1_x_red","COMB_slope1_red", "COMB_incrementTime_red")
var_list_all <- list("COMB_maximum_y","COMB_startPoint_x_red2","COMB_startPoint_x2","COMB_midPoint1_x","COMB_slope1", "COMB_incrementTime")
conditions <- c("Infected_sigmoidal", "Control", "QVD")
#Reconfigure data to give 3 levels of Productive (YES, NO, COMBINED)
plot_data_combined <- plot_data_with_sicegar_2_temp %>%
  filter(Condition_red=="MOI1") %>% 
  dplyr:::filter(DEC_decision_red!="ambiguous"&DEC_decision_red!="no_signal") %>% 
  group_by(unique_ID) %>% 
  slice(1) %>%
  mutate(productive = "COMBINED") 
# Combine the original and copied data
plot_data_all <- rbind(plot_data_with_sicegar_2_temp %>%
                         filter(Condition_red=="MOI1") %>% 
                         dplyr:::filter(DEC_decision_red!="ambiguous"&DEC_decision_red!="no_signal") %>% 
                         group_by(unique_ID) %>% 
                         slice(1), plot_data_combined)
#variables to t test
var_list_all_productive_vs_nonProductive <- c("COMB_startPoint_x_red2", "COMB_startDeclinePoint_x_red", "COMB_incrementTime_red")
#t test impact of Productive (YES vs NO) on var_list_all_productive_vs_nonProductive variables
perform_t_tests_P_vs_NP <- function(df, vars) {
  results_df <- data.frame()
  
  for (var in vars) {
    result <- pairwise.t.test(df[[var]], df$productive, p.adjust.method = "none")
    tidy_result <- broom::tidy(result)
    tidy_result$variable <- var
    
    results_df <- rbind(results_df, tidy_result)
  }
  results_df$p.value <- format(results_df$p.value, nsmall = 3)
  return(results_df)
}
signif_list_allVars_P_vs_NP <- purrr::map_dfr(var_list_all_productive_vs_nonProductive, .f=perform_t_tests_P_vs_NP, df= plot_data_all)

temp_lysis <- plot_data_all %>% 
  dplyr::group_by(productive) %>% 
  filter(productive=="YES"|productive=="NO") %>% 
  dplyr::select(unique_ID, HPI, COMB_startDeclinePoint_x_red,COMB_incrementTime_red) %>% 
  dplyr::summarise(mean=mean(COMB_startDeclinePoint_x_red, na.rm=TRUE),
                   mean=mean(COMB_incrementTime_red, na.rm=TRUE))

kn<- plot_data_with_sicegar_2_temp %>% 
  filter(Condition_red=="MOI1") %>% 
  dplyr:::filter(DEC_decision_red!="ambiguous"&DEC_decision_red!="no_signal") %>% 
  group_by(unique_ID) %>% 
  slice(1) %>% 
  dplyr::summarise(median_start=median(COMB_startPoint_x_red2))

#Panel K
# Create a subset of the data
df_subset <- plot_data_all %>%
  group_by(unique_ID,productive) %>% 
  slice(1) %>% 
  filter(DEC_decision_red %in% c("sigmoidal", "double_sigmoidal"))
# Create a contingency table of productive and lytic/non-lytic
contingency_table <- table(df_subset$productive, df_subset$DEC_decision_red)
# Perform chi-squared test
test <- chisq.test(contingency_table)

# Figure 3 ----------------------------------------------------------------

#Prep data for t test 
temp_tTest_lytic_nonLytic<- plot_data_with_sicegar_2_temp %>%   
  filter(Condition_red=="MOI1") %>% 
  group_by(unique_ID) %>% 
  slice(1) %>%
  filter(DEC_decision_red!="ambiguous"&DEC_decision_red!="no_signal") 

#temporary use for assessing difference between 12, 24, 36HPI cutoffs
temp_start<- plot_data %>% 
  group_by(unique_ID) %>% 
  slice(1) %>% 
  dplyr::select(start_time, start_time_green, unique_ID)
temp_tTest_lytic_nonLytic<- sicegar_results_df_merged %>%   
  filter(DEC_decision_red!="ambiguous"&DEC_decision_red!="no_signal") %>% 
  left_join(temp_start, by="unique_ID") %>% 
  dplyr::mutate(COMB_startPoint_x_red2=(start_time*10)/60,
                COMB_startPoint_x2=(start_time_green*10)/60) %>% 
  filter(Condition_red=="MOI1") %>% 
  filter(DEC_decision_red!="ambiguous"&DEC_decision_red!="no_signal") 

#List of variables to t test
var_list_all <- list("COMB_maximum_y","COMB_startPoint_x_red2","COMB_startPoint_x2","COMB_midPoint1_x","COMB_slope1", "COMB_incrementTime")
#t test impact of model class (lytic vs non-lytic) on var_list_all variables
perform_t_tests_l_vs_nl <- function(df, vars) {
  results_df <- data.frame()
  
  for (var in vars) {
    result <- pairwise.t.test(df[[var]], df$DEC_decision_red, p.adjust.method = "bonferroni")
    tidy_result <- broom::tidy(result)
    tidy_result$variable <- var
    
    results_df <- rbind(results_df, tidy_result)
  }
  results_df$p.value <- format(results_df$p.value, nsmall = 3)
  return(results_df)
}
signif_list_allVars_l_vs_nl <- purrr::map_dfr(var_list_all, .f=perform_t_tests_l_vs_nl, df= temp_tTest_lytic_nonLytic)


# Figure 4 ----------------------------------------------------------------

#Panel A
infection_proportion <- worked_filtered_infection_Scored %>% 
  filter(Condition!="uninfected") %>% 
  group_by(unique_ID) %>% 
  dplyr::summarise(Infected=first(Infected),
            productive=first(productive),
            Condition=first(Condition)) %>% 
  dplyr:::select(Condition, Infected, productive)
#Infection frequency
proportion_table_infected<- table(infection_proportion$Condition, infection_proportion$Infected)
chisq.test(table(infection_proportion$Infected, infection_proportion$Condition))
#Post-hoc
infection_proportion$Condition <- as.factor(infection_proportion$Condition) #combn() requires factors
# Create a list of all pairs of groups
group_pairs <- combn(levels(infection_proportion$Condition), 2, simplify = FALSE)
# Apply chi-square test to each pair
posthoc_tests <- map_dfr(group_pairs, function(pair) {
  sub_table <- proportion_table_infected[pair, ]
  test <- chisq.test(sub_table)
  tibble(pair1 = pair[1], pair2 = pair[2], p.value = test$p.value)
})
# Adjust p-values for multiple comparisons
posthoc_tests$p.adjusted <- p.adjust(posthoc_tests$p.value, method = "bonferroni")
# Print results
print(posthoc_tests)

#Productive infection frequency
proportion_table_productive<- table(infection_proportion$Condition, infection_proportion$productive)
chisq.test(table(infection_proportion$productive, infection_proportion$Condition))
#Post-hoc
infection_proportion$Condition <- as.factor(infection_proportion$Condition) #combn() requires factors
# Create a list of all pairs of groups
group_pairs <- combn(levels(infection_proportion$Condition), 2, simplify = FALSE)
# Apply chi-square test to each pair
posthoc_tests <- map_dfr(group_pairs, function(pair) {
  sub_table <- proportion_table_productive[pair, ]
  test <- chisq.test(sub_table)
  tibble(pair1 = pair[1], pair2 = pair[2], p.value = test$p.value)
})
# Adjust p-values for multiple comparisons
posthoc_tests$p.adjusted <- p.adjust(posthoc_tests$p.value, method = "bonferroni")
# Print results
print(posthoc_tests)

#Lytic frequency
# Create a subset of the data
df_subset_MOI <- plot_data_with_sicegar_2_temp %>%
  group_by(unique_ID,Condition_red) %>% 
  slice(1) %>% 
  filter(DEC_decision_red %in% c("sigmoidal","double_sigmoidal")) 
# Create a contingency table of productive and lytic/non-lytic
contingency_table_MOI <- table(df_subset_MOI$Condition_red, df_subset_MOI$DEC_decision_red)
# Perform chi-squared test
test_MOI <- chisq.test(contingency_table_MOI)

#Post-hoc
df_subset_MOI$Condition_red <- as.factor(df_subset_MOI$Condition_red) #combn() requires factors
# Create a list of all pairs of groups
group_pairs <- combn(levels(df_subset_MOI$Condition_red), 2, simplify = FALSE)
# Apply chi-square test to each pair
posthoc_tests <- map_dfr(group_pairs, function(pair) {
  sub_table <- contingency_table_MOI[pair, ]
  test <- chisq.test(sub_table)
  tibble(pair1 = pair[1], pair2 = pair[2], p.value = test$p.value)
})
# Adjust p-values for multiple comparisons
posthoc_tests$p.adjusted <- p.adjust(posthoc_tests$p.value, method = "bonferroni")
# Print results
print(posthoc_tests)

#Panel B
#Calcualte E-to-L delay
EL_t_test<- plot_data_with_sicegar_2_temp %>%
  group_by(unique_ID) %>% 
  slice(1) %>% 
  group_by(Condition_red) %>% 
  dplyr::mutate(EL_delay=COMB_startPoint_x2-COMB_startPoint_x_red2) 

perform_t_tests_E_to_L_delay <- function(df) {
  results_df <- data.frame()
  
  result <- pairwise.t.test(df[["EL_delay"]], df$Condition_red, p.adjust.method = "bonferroni")
  
  tidy_result <- broom::tidy(result)
  tidy_result$variable <- "EL_delay" 
  
  results_df <- rbind(results_df, tidy_result)
  
  return(results_df)
}

k<-perform_t_tests_E_to_L_delay(EL_t_test)

EL_t_test<- plot_data_with_sicegar_2_temp %>%
  group_by(unique_ID) %>% 
  slice(1) %>% 
  filter(Condition_red=="MOI1") %>% 
  filter(DEC_decision_red!="ambiguous"&DEC_decision_red!="no_signal") %>% 
  group_by(DEC_decision_red) %>% 
  dplyr::mutate(EL_delay=COMB_startPoint_x2-COMB_startPoint_x_red2) 

perform_t_tests_E_to_L_delay <- function(df) {
  results_df <- data.frame()
  
  result <- pairwise.t.test(df[["EL_delay"]], df$DEC_decision_red, p.adjust.method = "bonferroni")
  
  tidy_result <- broom::tidy(result)
  tidy_result$variable <- "EL_delay" 
  
  results_df <- rbind(results_df, tidy_result)
  
  return(results_df)
}

k<-perform_t_tests_E_to_L_delay(EL_t_test)

#Panel D
#Prep data for t test 
temp_tTest<-  plot_data_with_sicegar_2_temp %>% 
  dplyr:::filter(DEC_decision_red!="ambiguous"&DEC_decision_red!="no_signal") %>% 
  group_by(unique_ID) %>% 
  #filter(productive=="NO") %>% 
  slice(1)
#List of variables to t test
var_list_all <- list("COMB_maximum_y_red","COMB_midPoint1_x_red","COMB_slope1_red", "COMB_incrementTime_red")
#t test impact of model class (lytic vs non-lytic) on var_list_all variables
perform_t_tests_MOI <- function(df, vars) {
  results_df <- data.frame()
  
  for (var in vars) {
    result <- pairwise.t.test(df[[var]], df$Condition_red, p.adjust.method = "bonferroni")
    tidy_result <- broom::tidy(result)
    tidy_result$variable <- var
    
    results_df <- rbind(results_df, tidy_result)
  }
  
  return(results_df)
}
signif_list_allVars_MOI_all_red <- purrr::map_dfr(var_list_all, .f=perform_t_tests_MOI, df= temp_tTest)


#Panel C
#Prep data for t test 
temp_tTest<- plot_data_with_sicegar_2_temp %>% 
  group_by(unique_ID) %>% 
  slice(1) %>% 
  dplyr:::filter(DEC_decision_red!="ambiguous"&DEC_decision_red!="no_signal")

#List of variables to t test
var_list_all <- list("COMB_startPoint_x_red2","COMB_startPoint_x2","COMB_midPoint1_x","COMB_slope1", "COMB_incrementTime", "COMB_maximum_y")
#t test impact of model class (lytic vs non-lytic) on var_list_all variables
perform_t_tests_MOI <- function(df, vars) {
  results_df <- data.frame()
  
  for (var in vars) {
    # Perform the t-test
    result <- pairwise.t.test(df[[var]], df$Condition_red, p.adjust.method = "bonferroni")
    
    # Get sample sizes
    sample_sizes <- table(df$Condition_red)
    sample_size_str <- paste(names(sample_sizes), ":", as.vector(sample_sizes), collapse=", ")
    
    # Tidy the results
    tidy_result <- broom::tidy(result)
    tidy_result$variable <- var
    tidy_result$sample_sizes <- sample_size_str
    
    results_df <- rbind(results_df, tidy_result)
  }
  
  return(results_df)
}
signif_list_allVars_MOI_all <- purrr::map_dfr(var_list_all, .f=perform_t_tests_MOI, df= temp_tTest)

#Panel D/E/F
#requires MDS_matrix generated from the figures_cleanVersion_final script
#Modify MDS_matrix there to calculate dist for the correct set of variables

MDS_matrix <- data.frame(MDS_matrix, Condition_red = group_labels_temp$Condition.x)
# Create a data frame containing the distance matrix
dist_df <- as.matrix(dist(MDS_matrix[, -ncol(MDS_matrix)], method = "manhattan"))
# Perform PERMANOVA
permanova_result_host <- adonis2(dist_df ~ group_labels_temp$Condition.x)
# Print the result
print(permanova_result_host)  

# Figure 5 ----------------------------------------------------------------
#List of variables to t test
#Panel A
infection_proportion <- worked_filtered_infection_Scored_fulltimecourse %>% 
  filter(Condition=="Infected"|Condition=="QVD"|Condition=="QVD_Nec") %>% 
  group_by(unique_ID) %>% 
  dplyr::summarise(Infected=first(Infected),
                   productive=first(productive),
                   Condition=first(Condition)) %>% 
  dplyr:::select(Condition, Infected, productive)
#Infection frequency
proportion_table_infected<- table(infection_proportion$Condition, infection_proportion$Infected)
chisq.test(table(infection_proportion$Infected, infection_proportion$Condition))
#Post-hoc
infection_proportion$Condition <- as.factor(infection_proportion$Condition) #combn() requires factors
# Create a list of all pairs of groups
group_pairs <- combn(levels(infection_proportion$Condition), 2, simplify = FALSE)
# Apply chi-square test to each pair
posthoc_tests <- map_dfr(group_pairs, function(pair) {
  sub_table <- proportion_table_infected[pair, ]
  test <- chisq.test(sub_table)
  tibble(pair1 = pair[1], pair2 = pair[2], p.value = test$p.value)
})
# Adjust p-values for multiple comparisons
posthoc_tests$p.adjusted <- p.adjust(posthoc_tests$p.value, method = "bonferroni")
# Print results
print(posthoc_tests)

#Productive infection frequency
proportion_table_productive<- table(infection_proportion$Condition, infection_proportion$productive)
chisq.test(table(infection_proportion$productive, infection_proportion$Condition))
#Post-hoc
infection_proportion$Condition <- as.factor(infection_proportion$Condition) #combn() requires factors
# Create a list of all pairs of groups
group_pairs <- combn(levels(infection_proportion$Condition), 2, simplify = FALSE)
# Apply chi-square test to each pair
posthoc_tests <- map_dfr(group_pairs, function(pair) {
  sub_table <- proportion_table_productive[pair, ]
  test <- chisq.test(sub_table)
  tibble(pair1 = pair[1], pair2 = pair[2], p.value = test$p.value)
})
# Adjust p-values for multiple comparisons
posthoc_tests$p.adjusted <- p.adjust(posthoc_tests$p.value, method = "bonferroni")
# Print results
print(posthoc_tests)

#Lytic frequency
# Create a subset of the data
df_subset_drugs <- plot_data_with_sicegar_drugs_2%>%
  filter(Condition_red=="Infected"|Condition_red=="QVD"|Condition_red=="QVD_Nec") %>%
  group_by(unique_ID,Condition_red) %>% 
  slice(1) %>% 
  filter(DEC_decision_red %in% c("sigmoidal","double_sigmoidal"))
# Create a contingency table of productive and lytic/non-lytic
contingency_table_drugs <- table(df_subset_MOI$Condition_red, df_subset_MOI$DEC_decision_red)
# Perform chi-squared test
test_drugs <- chisq.test(contingency_table)
#pairing combinations
group_pairs <- combn(levels(df_subset_drugs$Condition_red), 3, simplify = FALSE)
# Apply chi-square test to each pair
posthoc_tests <- map_dfr(group_pairs, function(pair) {
  sub_table <- proportion_table_productive[pair, ]
  test <- chisq.test(sub_table)
  tibble(pair1 = pair[1], pair2 = pair[2], p.value = test$p.value)
})
# Adjust p-values for multiple comparisons
posthoc_tests$p.adjusted <- p.adjust(posthoc_tests$p.value, method = "bonferroni")
# Print results
print(posthoc_tests)

#Panel B
#Calcualte E-to-L delay
EL_t_test<- plot_data_with_sicegar_drugs_2 %>% 
  group_by(unique_ID) %>% 
  slice(1) %>% 
  filter(Condition_red=="Infected"|Condition_red=="QVD"|Condition_red=="QVD_Nec") %>% 
  group_by(unique_ID, Condition_red) %>% 
  dplyr::summarise(EL_delay=COMB_startPoint_x-COMB_startPoint_x_red)


perform_t_tests_E_to_L_delay <- function(df) {
  results_df <- data.frame()
  
  # Removed the loop, as it is not necessary if you only want to test against 'mean'
  result <- pairwise.t.test(df[["EL_delay"]], df$Condition_red, p.adjust.method = "bonferroni")
  
  tidy_result <- broom::tidy(result)
  tidy_result$variable <- "EL_delay" # Corrected the variable name assignment
  
  results_df <- rbind(results_df, tidy_result)
  
  return(results_df)
}

drugs_EL_delay_results<-perform_t_tests_E_to_L_delay(EL_t_test)

#Panel C
var_list_all <- c("COMB_maximum_y","COMB_startPoint_x_red2","COMB_startPoint_x2","COMB_midPoint1_x","COMB_slope1", "COMB_incrementTime")
var_list_all_character <- c("COMB_startPoint_x_red2", "COMB_startPoint_x2",
                            "COMB_midPoint1_x", "COMB_slope1", "COMB_incrementTime","COMB_maximum_y")
var_list_all_character <- c("COMB_startPoint_x_red2",
                            "COMB_midPoint1_x_red", "COMB_slope1_red", "COMB_incrementTime_red","COMB_maximum_y_red")
#filter data
temp_tTest_drugs<- plot_data_with_sicegar_drugs_2_fig6_7%>% 
  #filter(unique_ID!="20221010_2_0_14") %>%  #extreme outlier
  filter(DEC_decision=="sigmoidal"|DEC_decision=="double_sigmoidal") %>% 
  filter(decision_Condition=="Infected_sigmoidal"|decision_Condition=="Control"|decision_Condition=="QVD") 
#t test impact of model class (lytic vs non-lytic) on var_list_all variables
perform_t_tests_drugs <- function(df, vars) {
  results_df <- data.frame()
  
  for (var in vars) {
    # Perform the t-test
    result <- pairwise.t.test(df[[var]], df$decision_Condition, p.adjust.method = "no")
    
    # Get sample sizes
    sample_sizes <- table(df$decision_Condition)
    sample_size_str <- paste(names(sample_sizes), ":", as.vector(sample_sizes), collapse=", ")
    
    # Tidy the results
    tidy_result <- broom::tidy(result)
    tidy_result$variable <- var
    tidy_result$sample_sizes <- sample_size_str
    
    results_df <- rbind(results_df, tidy_result)
  }
  
  return(results_df)
}  #No multiple comparison correction for QVD>Control, QVD>Control non-lytic
#Requires the df plot_data_with_sicegar_drugs_2_fig6_7 from figures_cleanVersion_final script
signif_list_allVars_drugs <- purrr::map_dfr(var_list_all, .f=perform_t_tests_drugs, df= temp_tTest_drugs)

#Panel E
#Stats
temp_QVD_diameter <- data_size %>% 
  filter(Condition=="QVD") %>% 
  pull(Feret) 
temp_control_diameter <- data_size %>% 
  filter(Condition=="Control") %>% 
  pull(Feret) 
t.test(temp_QVD_diameter, temp_control_diameter)

#Panel F
temp_QVD_count <- data_count %>% 
  dplyr::rename(Condition=ï..Condition) %>% 
  filter(Condition=="QVD") %>% 
  pull(Count) 
temp_control_count <- data_count %>% 
  dplyr::rename(Condition=ï..Condition) %>% 
  filter(Condition=="Control") %>% 
  pull(Count) 
t.test(temp_QVD_count, temp_control_count)

#Figure 6 #####
#Calculating the PFU range for 95% of cells at various MOIs
MOI <- 1

# Find the 2.5th and 97.5th percentiles
lower_bound <- qpois(0.025, MOI)
upper_bound <- qpois(0.975, MOI)

# Print the range
cat("The range of particle numbers that 95% of cells would receive is:", lower_bound, "to", upper_bound, "\n")
#Calculate means and variance of single cell infection parameters ####

# Calculate the mean, standard deviation, and coefficient of variation for each variable in var_list_all
stats_sd_productive_vs_nonProductive <- plot_data_with_sicegar_2_temp %>%
  filter(Condition_red=="MOI1") %>% 
  dplyr:::filter(DEC_decision_red!="ambiguous"&DEC_decision_red!="no_signal") %>% 
  group_by(productive) %>%
  dplyr::summarise(across(all_of(var_list_all_productive_vs_nonProductive),
                          list(mean = ~mean(., na.rm = TRUE), 
                               std_dev = ~sd(., na.rm = TRUE),
                               coeff_of_variation = ~sd(., na.rm = TRUE) / mean(., na.rm = TRUE))
  )) %>%
  # Reshape the data using a regular expression
  pivot_longer(-productive, 
               names_to = c("variable", ".value"),
               names_pattern = "(.+)_(mean|std_dev|coeff_of_variation)",
               values_to = c("mean", "std_dev", "coeff_of_variation"))
#Figure 3 related
stats_sd_lytic_vs_nonLytic <- plot_data_with_sicegar_2_temp_fig3 %>%
  filter(DEC_decision_red=="double_sigmoidal"|DEC_decision_red=="sigmoidal") %>% 
  group_by(DEC_decision_red_2) %>%
  dplyr::summarise(across(all_of(var_list_all_character),
                          list(mean = ~mean(., na.rm = TRUE), 
                               std_dev = ~sd(., na.rm = TRUE),
                               coeff_of_variation = ~sd(., na.rm = TRUE) / mean(., na.rm = TRUE))
  )) %>%
  # Reshape the data using a regular expression
  pivot_longer(-DEC_decision_red_2, 
               names_to = c("variable", ".value"),
               names_pattern = "(.+)_(mean|std_dev|coeff_of_variation)",
               values_to = c("mean", "std_dev", "coeff_of_variation"))

#Figure 4 related
var_list_all_character <- c("COMB_startPoint_x_red2", "COMB_startPoint_x2",
                            "COMB_midPoint1_x", "COMB_slope1", "COMB_incrementTime","COMB_maximum_y")
stats_sd_MOI <- plot_data_with_sicegar_2_temp %>%
  group_by(unique_ID) %>% 
  slice(1) %>% 
  dplyr:::filter(DEC_decision_red!="ambiguous"&DEC_decision_red!="no_signal") %>% 
  group_by(Condition_red) %>%
  dplyr::summarise(across(all_of(var_list_all_character),
                          list(mean = ~mean(., na.rm = TRUE), 
                               std_dev = ~sd(., na.rm = TRUE),
                               coeff_of_variation = ~sd(., na.rm = TRUE) / mean(., na.rm = TRUE))
  )) %>%
  # Reshape the data using a regular expression
  pivot_longer(-Condition_red, 
               names_to = c("variable", ".value"),
               names_pattern = "(.+)_(mean|std_dev|coeff_of_variation)",
               values_to = c("mean", "std_dev", "coeff_of_variation"))


#Figure 5 related
stats_sd_drugs <- plot_data_with_sicegar_drugs_2_fig6_7 %>%
  filter(decision_Condition=="Infected_sigmoidal"|decision_Condition=="Control"|decision_Condition=="QVD") %>% 
  group_by(decision_Condition) %>%
  dplyr::summarise(across(all_of(var_list_all_character),
                          list(mean = ~mean(., na.rm = TRUE), 
                               std_dev = ~sd(., na.rm = TRUE),
                               coeff_of_variation = ~sd(., na.rm = TRUE) / mean(., na.rm = TRUE))
  )) %>%
  # Reshape the data using a regular expression
  pivot_longer(-decision_Condition, 
               names_to = c("variable", ".value"),
               names_pattern = "(.+)_(mean|std_dev|coeff_of_variation)",
               values_to = c("mean", "std_dev", "coeff_of_variation"))



#Extended data ####
#Extended data 1 ####
#Prep data for t test 
temp_start<- plot_data %>% 
  filter(HPI<=24) %>% 
  group_by(unique_ID) %>% 
  slice(1) %>% 
  dplyr::select(start_time, start_time_green, unique_ID)

temp_tTest_lytic_nonLytic_24HPI<- sicegar_results_df_merged_24HPI %>%  
  filter(Condition_red=="MOI1") %>% 
  filter(DEC_decision_red!="ambiguous"&DEC_decision_red!="no_signal") %>% 
  left_join(temp_start, by="unique_ID") %>% 
  mutate(COMB_startPoint_x_red2=((start_time*10)/60),
         COMB_startPoint_x2=((start_time_green*10)/60)) 
  
#List of variables to t test
var_list_all <- list("COMB_startPoint_x_red2", "COMB_startPoint_x2", "COMB_midPoint1_x", 
                  "COMB_slope1", "COMB_incrementTime", "COMB_maximum_y")

#t test impact of model class (lytic vs non-lytic) on var_list_all variables
perform_t_tests_l_vs_nl <- function(df, vars) {
  results_df <- data.frame()
  
  for (var in vars) {
    result <- pairwise.t.test(df[[var]], df$DEC_decision_red, p.adjust.method = "bonferroni")
    tidy_result <- broom::tidy(result)
    tidy_result$variable <- var
    
    results_df <- rbind(results_df, tidy_result)
  }
  results_df$p.value <- format(results_df$p.value, nsmall = 3)
  return(results_df)
}
signif_list_allVars_l_vs_nl_12HPI <- purrr::map_dfr(var_list_all, .f=perform_t_tests_l_vs_nl, df= temp_tTest_lytic_nonLytic_12HPI)
perform_t_tests_l_vs_nl <- function(df, vars) {
  results_df <- data.frame()
  
  for (var in vars) {
    result <- pairwise.t.test(df[[var]], df$DEC_decision_red, p.adjust.method = "bonferroni")
    tidy_result <- broom::tidy(result)
    tidy_result$variable <- var
    
    results_df <- rbind(results_df, tidy_result)
  }
  results_df$p.value <- format(results_df$p.value, nsmall = 3)
  return(results_df)
}
signif_list_allVars_l_vs_nl_24HPI <- purrr::map_dfr(var_list_all, .f=perform_t_tests_l_vs_nl, df= temp_tTest_lytic_nonLytic_24HPI)

#Extended data 2 ####
QVD_oneStep <- read.csv("QVD_oneStep_growth.csv")
QVD_oneStep_t_test <- QVD_oneStep %>% 
  filter(Source == "Total") %>% 
  pivot_longer(cols=4:5, names_to = "Condition", values_to = "Titre") %>% 
  group_by(HPI) %>% 
  pairwise_t_test(Titre ~ Condition, p.adjust.method = "bonferroni")

#Extended data table 1 ####
#MOI data Early
MOI_data_model_frequencies_Early<- plot_data_with_sicegar_2_temp %>% 
  group_by(unique_ID) %>% 
  slice(1) %>% 
  group_by(Condition_red, DEC_decision_red) %>% 
  dplyr::summarise(n=n())
#MOI data Late
MOI_data_model_frequencies_Late<- plot_data_with_sicegar_2_temp %>% 
  group_by(unique_ID) %>% 
  slice(1) %>% 
  group_by(Condition_red, DEC_decision) %>% 
  dplyr::summarise(n=n())
#Drugs data Early
Drugs_data_model_frequencies_Early<- plot_data_with_sicegar_drugs_2 %>% 
  filter(Condition_red=="Infected"|Condition_red=="QVD") %>% 
  group_by(unique_ID) %>% 
  slice(1) %>% 
  group_by(Condition_red, DEC_decision_red) %>% 
  dplyr::summarise(n=n())
#Drugs data Late
Drugs_data_model_frequencies_Late<- plot_data_with_sicegar_drugs_2 %>% 
  filter(Condition_red=="Infected"|Condition_red=="QVD") %>% 
  group_by(unique_ID) %>% 
  slice(1) %>% 
  group_by(Condition_red, DEC_decision) %>% 
  dplyr::summarise(n=n())

#Extended data table 2 ####
gini_df<- plot_data_with_sicegar_2_temp %>%
  #filter(DEC_decision_red=="double_sigmoidal") %>% 
  group_by(unique_ID) %>% 
  slice(1) 

# List of variables to calculate Gini coefficient
var_list <- c("COMB_startPoint_x_red2","COMB_startPoint_x2",
              "COMB_midPoint1_x","COMB_slope1", "COMB_incrementTime", "COMB_maximum_y")

# Function to calculate Gini coefficient
gini_calc <- function(x) {
  return(ineq::Gini(x))
}

# Calculate Gini coefficient for each variable grouped by Condition_red
result <- gini_df %>%
  group_by(Condition_red) %>%
  dplyr::summarise(across(all_of(var_list), gini_calc, .names = "Gini_{.col}")) 

print(result)
#Extended data table 3 ####
MOI_data_model_frequencies_P_vs_NP_Early<- plot_data_with_sicegar_2_temp %>% 
  group_by(unique_ID) %>% 
  slice(1) %>% 
  filter(DEC_decision_red=="double_sigmoidal"|DEC_decision_red=="sigmoidal") %>% 
  group_by(Condition_red, productive, DEC_decision_red) %>% 
  dplyr::summarise(n=n())
#Extended data table 4 ####
uninfected_dead_count<- read.csv("uninfected_live_dead_count_drugs_data.csv")
uninfected_dead_count_long<- uninfected_dead_count %>% 
  dplyr::select(!Frequency) %>% 
  pivot_longer(cols=3:4, names_to = "State", values_to = "Count") %>% 
  mutate(Well_ID = str_split(FoV_ID, "_", simplify = TRUE)[, 1]) %>% 
  group_by(Well_ID, Condition, State) %>% 
  dplyr::summarise(sum_count=sum(Count)) %>% 
  pivot_wider(names_from = State, values_from = sum_count) %>% 
  mutate(Proportion=Dead/(Dead+Live))




#Supplementary data ####
#Supplementary Figure 1 ####
temp_tTest<- sicegar_results_df_merged_2%>% 
  filter(DEC_decision_red=="double_sigmoidal"|DEC_decision_red=="sigmoidal")
#List of variables to t test
var_list_all <- list("COMB_startPoint_x_red","COMB_startPoint_x","COMB_midPoint1_x","COMB_slope1", "COMB_incrementTime", "COMB_maximum_y")
#t test impact of model class (lytic vs non-lytic) on var_list_all variables
perform_t_tests_MOI <- function(df, vars) {
  results_df <- data.frame()
  
  for (var in vars) {
    # Perform the t-test
    result <- pairwise.t.test(df[[var]], df$Condition_red, p.adjust.method = "bonferroni")
    
    # Get sample sizes
    sample_sizes <- table(df$Condition_red)
    sample_size_str <- paste(names(sample_sizes), ":", as.vector(sample_sizes), collapse=", ")
    
    # Tidy the results
    tidy_result <- broom::tidy(result)
    tidy_result$variable <- var
    tidy_result$sample_sizes <- sample_size_str
    
    results_df <- rbind(results_df, tidy_result)
  }
  
  return(results_df)
}
signif_list_allVars_MOI_all <- purrr::map_dfr(var_list_all, .f=perform_t_tests_MOI, df= temp_tTest)
#Supplementary Figure 3 ####
#read in data
plaque_characterisation <- read.csv("A3_mCh_vacv.csv") %>% 
  dplyr::rename(Virus=ï..Virus) 
plaque_characterisation$Virus <- factor(plaque_characterisation$Virus, 
                                        levels = c("WR", "A3GFP", "mCh", "A3GFP_mCh"))
# Compute the analysis of variance
res.aov <- aov(Diameter ~ Virus, data = plaque_characterisation)
# Summary of the analysis
summary(res.aov)
TukeyHSD(res.aov)
#Supplementary Figure 4 ####
data_AraC <- read.csv("data_AraC.csv")
data_AraC_normalised<- data_AraC %>% 
  mutate(normalised_green=green_IntDen-(Background_green*Area),
         normalised_red=red_IntDen-(Background_red*Area))
data_AraC_normalised$Condition <- factor(data_AraC_normalised$Condition, 
                                         levels = c("Control", "AraC"))
data_AraC_normalised_filtered<- data_AraC_normalised %>% 
  filter(Time=="24")
# Compute the analysis of variance
res.aov <- aov(normalised_red ~ Condition, data = data_AraC_normalised_filtered)
# Summary of the analysis
summary(res.aov)
TukeyHSD(res.aov)
#Supplementary table 1 ####
#Prep data for t test 
temp_tTest<- sicegar_results_df_merged_2%>% 
  filter(DEC_decision_red=="double_sigmoidal")
#List of variables to t test
var_list_all <- list("COMB_startPoint_x_red","COMB_startPoint_x","COMB_midPoint1_x","COMB_slope1", "COMB_incrementTime", "COMB_maximum_y")
#t test impact of model class (lytic vs non-lytic) on var_list_all variables
perform_t_tests_MOI <- function(df, vars) {
  results_df <- data.frame()
  
  for (var in vars) {
    # Perform the t-test
    result <- pairwise.t.test(df[[var]], df$Condition_red, p.adjust.method = "bonferroni")
    
    # Get sample sizes
    sample_sizes <- table(df$Condition_red)
    sample_size_str <- paste(names(sample_sizes), ":", as.vector(sample_sizes), collapse=", ")
    
    # Tidy the results
    tidy_result <- broom::tidy(result)
    tidy_result$variable <- var
    tidy_result$sample_sizes <- sample_size_str
    
    results_df <- rbind(results_df, tidy_result)
  }
  
  return(results_df)
}
signif_list_allVars_MOI_all <- purrr::map_dfr(var_list_all, .f=perform_t_tests_MOI, df= temp_tTest)

#filter data
temp_tTest_drugs<- plot_data_with_sicegar_drugs_2_fig6_7%>% 
  #filter(unique_ID!="20221010_2_0_14") %>%  #extreme outlier
  filter(DEC_decision=="sigmoidal"|DEC_decision=="double_sigmoidal") %>% 
  filter(decision_Condition=="Infected_sigmoidal"|decision_Condition=="Control"|decision_Condition=="QVD") 
#t test impact of model class (lytic vs non-lytic) on var_list_all variables
perform_t_tests_drugs <- function(df, vars) {
  results_df <- data.frame()
  
  for (var in vars) {
    # Perform the t-test
    result <- pairwise.t.test(df[[var]], df$decision_Condition, p.adjust.method = "no")
    
    # Get sample sizes
    sample_sizes <- table(df$decision_Condition)
    sample_size_str <- paste(names(sample_sizes), ":", as.vector(sample_sizes), collapse=", ")
    
    # Tidy the results
    tidy_result <- broom::tidy(result)
    tidy_result$variable <- var
    tidy_result$sample_sizes <- sample_size_str
    
    results_df <- rbind(results_df, tidy_result)
  }
  
  return(results_df)
}  #No multiple comparison correction for QVD>Control, QVD>Control non-lytic
#Requires the df plot_data_with_sicegar_drugs_2_fig6_7 from figures_cleanVersion_final script
signif_list_allVars_drugs <- purrr::map_dfr(var_list_all, .f=perform_t_tests_drugs, df= temp_tTest_drugs)


temp1 <- plot_data_with_sicegar_2_temp %>% 
  group_by(unique_ID) %>% 
  slice(1) %>% 
  filter(Condition_red=="MOI100")%>%
  dplyr::select(contains("COMB_"))
