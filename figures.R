# Load required libraries ####
library(RColorBrewer)
library(ggplot2)
library(GGally)
library(patchwork)
library(gridExtra)
library(dplyr)
require("cowplot")
require("stringr")
require("rlist")
require("rstatix")
library("eply")
library(ggplot2)
library(tidyr)
library(plyr)
library(gganimate)
library(transformr)
library(dplyr)
library(tidyverse)
library(broom)
library(magick)
library(patchwork)
library(zoo)
library(ggsignif)
library(qpcR)
library(car) 
library(FSA)
library(ggrepel)
library(gridExtra) # merge plots
library(gplots) # heatmap
library(tseries) # bootstrap
library(tibble)
library(data.table)
library(chisq.posthoc.test)
library(rlang)
library(extrafont)
library(showtext)
library(scales)

#Font loading
showtext_auto()
font_add("Arial", "C:/Windows/Fonts/arial.ttf")


# Figure 1 ----------------------------------------------------------------
showtext_opts(dpi = 300) #set text DPI
sample_groups <- function(df, unique_ID, num_sample_groups) {
  # Get unique groups
  unique_groups <- df %>% distinct({{unique_ID}})
  
  # Check number of unique groups
  if(nrow(unique_groups) < num_sample_groups){
    warning("Number of unique groups is less than the number of samples required.")
    num_sample_groups <- nrow(unique_groups)
  }
  
  # Sample groups
  sampled_groups <- unique_groups %>% ungroup() %>% 
    dplyr::sample_n(num_sample_groups)
  
  # Filter original data for sampled groups
  sampled_data <- df %>% 
    semi_join(sampled_groups, by = deparse(substitute(unique_ID)))
  
  return(sampled_data)
} #randomly sample n cells

# Select a subset of data where the condition (Multiplicity of Infection, MOI) 
# is "MOI1", and sample 50 unique cells.
temp_sample_keep <- plot_data %>%
  filter(Condition=="MOI1") %>% 
  sample_groups(unique_ID, 50)

# Define the color palette for the plots.
color_palette <- c("grey","blue","magenta","red", "darkorange1")


# Create a line plot where the green fluorescence intensity (A3-GFP) colors the lines that plot time post-infection (HPI) against the normalized red fluorescence intensity (pEL-mCh).
temp_sample_keep <- temp_sample_keep %>%
  arrange(unique_ID, HPI)


p_combined <- temp_sample_keep  %>% 
  ggplot(aes(x=HPI, y=abs(normalised_sum_red), group = unique_ID, colour = normalised_sum_green)) +
  geom_path(size=0.25,lineend="round",linejoin="mitre") +
  scale_color_gradientn(colors = color_palette)+
  scale_x_continuous()+
  labs(x = "Time (HPI)", y = "Early (NFI)",
       title = "", colour = "Late \n(NFI)")+
  theme_minimal()+
  theme(legend.position = c(0.5, 1.05),
        legend.direction = "horizontal",
        legend.key.height = unit(0.2, 'cm'),
        legend.key.width = unit(1, 'cm'),
        axis.text.y = element_blank(),
        axis.text = element_text(size = 5),
        axis.title = element_text(size = 7),
        text=element_text(family ="Arial"),
        legend.title = element_text(size = 7, vjust=3),
        legend.text = element_blank())+
  scale_y_sqrt()

# Save the plot as a TIFF file.
ggsave("single_MOI_combined_linePlot.pdf", width=60, height=60, units="mm", p_combined)

# Figure 2 ----------------------------------------------------------------
showtext_opts(dpi = 300) #Set font dpi
# Calculate the mean and standard error of the green and red fluorescence intensities for each time point (HPI) 
#and condition, then plot the mean green fluorescence intensity (A3-GFP) with error bars for each condition over time.
temp_mean<-plot_data_with_sicegar_2_temp %>% 
  filter(Condition_red=="MOI1") %>% 
  group_by(Condition_red,HPI) %>% 
  dplyr::summarise(meanGreen=mean(normalised_sum_green),
                   sdGreen=sd(normalised_sum_green),
                   seGreen=(sd(normalised_sum_green)/sqrt(217)),
                   meanRed=mean(normalised_sum_red),
                   sdRed=sd(normalised_sum_red),
                   seRed=(sd(normalised_sum_red)/sqrt(217)),
                   start_red=mean(COMB_startPoint_x_red2, na.rm=TRUE),
                   start_green=mean(COMB_startPoint_x2, na.rm=TRUE))

#Panel A
p <- temp_mean  %>% 
  ggplot(aes(x=HPI, y=meanRed+1)) +
  geom_path(size=0.5, alpha=1, color="red2") +
  scale_x_continuous()+
  geom_ribbon(aes(ymin=meanRed+1-seRed, ymax=meanRed+1+seRed, alpha=0.5),
              position=position_dodge(0.05), alpha=0.5)+
  labs(x = "Time (HPI)", y = "Mean Early (NFI)",
       title = "", colour = "")+
  theme_minimal()+
  theme(legend.position = "none",
        axis.text = element_text(size = 5),
        axis.title = element_text(size = 7),
        text=element_text(family ="Arial"))+
  scale_y_sqrt()
# Save the plot as a TIFF file.
ggsave("mean_red_allMOI_linePlot.pdf",width=60, height=40, units="mm", p)

#Panel B
p <- temp_mean  %>% 
  ggplot(aes(x=HPI, y=meanGreen+1)) +
  geom_line(size=1, alpha=1, color="green2") +
  scale_x_continuous()+
  geom_ribbon(aes(ymin=meanGreen+1-seGreen, ymax=meanGreen+1+seGreen, alpha=0.5),
              position=position_dodge(0.05), alpha=0.5)+
  labs(x = "Time (HPI)", y = "Mean Late (NFI)",
       title = "", colour = "")+
  theme_minimal()+
  theme(legend.position = "none",
        axis.text = element_text(size = 5),
        axis.title = element_text(size = 7),
        text=element_text(family ="Arial"))+
  scale_y_sqrt()
# Save the plot as a TIFF file.
ggsave("mean_green_allMOI_linePlot.pdf",width=60, height=40, units="mm", p)

#Panel C
p_single_channel <- temp_sample_keep  %>% 
  ggplot(aes(x=HPI, y=abs(normalised_sum_red), group = unique_ID)) +
  geom_line(size=0.25, color="black") +
  scale_color_gradientn(colors = color_palette)+
  scale_x_continuous()+
  labs(x = "Time (HPI)", y = "Early (NFI)")+
  theme_minimal()+
  theme(legend.position = "none",
        axis.text = element_text(size = 5),
        axis.title = element_text(size = 7),
        text=element_text(family ="Arial"))+
  scale_y_sqrt()
# Save the plot as a TIFF file.
ggsave("single_MOI_red_linePlot.pdf",width=60, height=40, units="mm", p_single_channel)

#Panel D
p_single_channel <- temp_sample_keep  %>% 
  ggplot(aes(x=HPI, y=abs(normalised_sum_green), group = unique_ID)) +
  geom_line(size=0.25, color="black") +
  scale_color_gradientn(colors = color_palette)+
  scale_x_continuous()+
  labs(x = "Time (HPI)", y = "Late (NFI)")+
  theme_minimal()+
  theme(legend.position = "none",
        axis.text = element_text(size = 5),
        axis.title = element_text(size = 7),
        text=element_text(family ="Arial"))+
  scale_y_sqrt()
# Save the plot as a TIFF file.
ggsave("single_MOI_green_linePlot.pdf",width=60, height=40, units="mm", p_single_channel)

#Panel C,D
# Filter the raw data for a single cell, and select the columns corresponding to time (HPI) and red (or green) fluorescence intensity.
# Rename these columns to "time" and "intensity".
temp_sm <- plot_data %>% 
  dplyr::filter(unique_ID=="20230402_10_4_13") %>% 
  ungroup() %>%
  dplyr::select(HPI, normalised_sum_red) %>% 
  dplyr::rename("time"=HPI, "intensity"=normalised_sum_red)
temp_dsm <- plot_data %>% 
  dplyr::filter(unique_ID=="20230402_10_0_15") %>%
  ungroup() %>% 
  dplyr::select(HPI, normalised_sum_red) %>% 
  dplyr::rename("time"=HPI, "intensity"=normalised_sum_red)
# Normalize the data using the sicegar package.
temp1_sm<- sicegar::normalizeData(temp_sm)
temp1_dsm<- sicegar::normalizeData(temp_dsm)
# Fit the normalized data to a sigmoidal model using sicegar.
sigmoidalModel <- sicegar::multipleFitFunction(dataInput = temp1_sm,
                                               model = "sigmoidal",
                                               n_runs_min = 20,
                                               n_runs_max = 500,
                                               showDetails = FALSE)
# Calculate the parameters of the sigmoidal model.
sigmoidalModel <- sicegar::parameterCalculation(sigmoidalModel)
# Fit the normalized data to a double sigmoidal model using sicegar.
doubleSigmoidalModel <- multipleFitFunction(dataInput = temp1_dsm,
                                            model = "doublesigmoidal",
                                            n_runs_min = 20,
                                            n_runs_max = 500,
                                            showDetails = FALSE)
# Calculate the parameters of the double sigmoidal model.
doubleSigmoidalModel <- parameterCalculation(doubleSigmoidalModel)
# Generate plots of the data along with the fitted sigmoidal and double sigmoidal models.
fig01 <- sicegar::figureModelCurves(dataInput = temp1_sm,
                                    sigmoidalFitVector = sigmoidalModel,
                                    showParameterRelatedLines = FALSE)
fig02 <- figureModelCurves(dataInput = temp1_dsm,
                           doubleSigmoidalFitVector = doubleSigmoidalModel,
                           showParameterRelatedLines = FALSE)
# Extract the environments from fig02 (for the double sigmoidal model) and fig01 (for the sigmoidal model).
env <- fig02[["plot_env"]] 
env_2 <- fig01[["plot_env"]] 
# Access the list from the environment
fitted_model <- get("intensityTheoreticalDoubleSigmoidalDf", envir = env) 
raw_data <- get("dataFrameInput", envir = env)

fitted_model_2 <- get("intensityTheoreticalSigmoidalDf", envir = env_2)%>% 
  filter(time<16)
raw_data_2 <- get("dataFrameInput", envir = env_2)%>% 
  filter(time<16)
# Generate a plot for the double sigmoidal infection.

p_DSM <- ggplot()+
  geom_point(data = raw_data, aes_(x = ~time, y = ~intensity), size=0.25) + 
  geom_line(data = fitted_model, 
            ggplot2::aes_(x = ~time, y = ~intensityTheoreticalDoubleSigmoidal), 
            color = "red2", size = 0.5)+
  labs(title= "", y="Early (NFI)",  x="Time (HPI)")+
  theme_void()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "",
        axis.text = element_blank(),
        axis.title.y = element_text(size = 5, angle=90, vjust = -3),
        axis.title.x = element_text(size = 5, vjust=1),
        text=element_text(family ="Arial"))+
  coord_cartesian(ylim=c(-20,600), xlim=c(-2, 37))+
  annotate("segment", x = 0, xend = 37, y = -0.05, yend = -0.05, arrow = arrow(length = unit(0.125, "cm"),type="closed"), lwd=0.5) +
  annotate("segment", x = 0, xend = 0, y = 0, yend = 600, arrow = arrow(length = unit(0.125, "cm"),type="closed"), lwd=0.5)
ggsave("DSM_linePlot.pdf",width=30, height=20, units="mm", p_DSM)

# Generate a plot for the sigmoidal infection.
p_SM <- ggplot()+
  geom_point(data = raw_data_2, aes_(x = ~time, y = ~intensity), size=0.25) + 
  geom_line(data = fitted_model_2, 
            ggplot2::aes_(x = ~time, y = ~intensityTheoreticalSigmoidal), 
            color = "red2", size = 0.5)+
  labs(title= "", y="Early (NFI)",  x="Time (HPI)")+
  theme_void()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "",
        axis.text = element_blank(),
        axis.title.y = element_text(size = 5, angle=90, vjust = -3),
        axis.title.x = element_text(size = 5, vjust=1),
        text=element_text(family ="Arial"))+
  coord_cartesian(ylim=c(-20,600), xlim=c(-1, 17))+ #xlim for plot v2 only.
  annotate("segment", x = 0, xend = 17, y = -0.05, yend = -0.05, arrow = arrow(length = unit(0.125, "cm"),type="closed"), lwd=0.5) +
  annotate("segment", x = 0, xend = 0, y = 0, yend = 600, arrow = arrow(length = unit(0.125, "cm"),type="closed"), lwd=0.5)
ggsave("SM_red_linePlot.pdf",width=30, height=20, units="mm", p_SM)

temp_ns_green <- plot_data %>% 
  dplyr::filter(unique_ID=="20230402_1_0_29") %>% 
  ungroup() %>%
  dplyr::select(HPI, normalised_sum_green) %>% 
  dplyr::rename("time"=HPI, "intensity"=normalised_sum_green)
# Normalize the data using the sicegar package.
temp1_ns_green<- sicegar::normalizeData(temp_ns_green)

NS_green <- sicegar::multipleFitFunction(dataInput = temp1_ns_green,
                                         model = "sigmoidal",
                                         n_runs_min = 20,
                                         n_runs_max = 500,
                                         showDetails = FALSE)
# Calculate the parameters of the sigmoidal model.
NS_green <- sicegar::parameterCalculation(NS_green)
# Generate plots of the data along with the fitted sigmoidal and double sigmoidal models.
fig02_green <- sicegar::figureModelCurves(dataInput = temp1_ns_green,
                                          sigmoidalFitVector = NS_green,
                                          showParameterRelatedLines = FALSE)
env_3_green <- fig02_green[["plot_env"]] 
# Access the list from the environment
fitted_model_2_green_ns <- get("intensityTheoreticalSigmoidalDf", envir = env_3_green) %>% 
  filter(time<15)
raw_data_3_green <- get("dataFrameInput", envir = env_3_green)%>% 
  filter(time<15)


p_NS_green <- ggplot()+
  geom_point(data = raw_data_3_green, aes_(x = ~time, y = ~intensity), size=0.25) + 
  geom_line(data = fitted_model_2_green_ns, 
            ggplot2::aes_(x = ~time, y = ~intensityTheoreticalSigmoidal), 
            color = "green2", size = 0.5)+
  labs(title= "", y="Late (NFI)",  x="Time (HPI)")+
  theme_void()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "",
        axis.text = element_blank(),
        axis.title.y = element_text(size = 5, angle=90, vjust = -3),
        axis.title.x = element_text(size = 5, vjust=1),
        text=element_text(family ="Arial"))+
  coord_cartesian(ylim=c(-10,180), xlim=c(-1, 16))+ #xlim for plot v2 only.
  annotate("segment", x = 0, xend = 16, y = -5, yend = -5, arrow = arrow(length = unit(0.125, "cm"),type="closed"), lwd=0.5) +
  annotate("segment", x = 0, xend = 0, y = -5.2, yend = 180, arrow = arrow(length = unit(0.125, "cm"),type="closed"), lwd=0.5)
# Save the plot
ggsave("No_signal_linePlot.pdf",width=30, height=20, units="mm", p_NS_green)

temp_sm_green <- plot_data %>% 
  dplyr::filter(unique_ID=="20230402_9_5_32") %>% 
  ungroup() %>%
  dplyr::select(HPI, normalised_sum_green) %>% 
  dplyr::rename("time"=HPI, "intensity"=normalised_sum_green)
# Normalize the data using the sicegar package.
temp1_sm_green<- sicegar::normalizeData(temp_sm_green)

sigmoidalModel_green <- sicegar::multipleFitFunction(dataInput = temp1_sm_green,
                                                     model = "sigmoidal",
                                                     n_runs_min = 20,
                                                     n_runs_max = 500,
                                                     showDetails = FALSE)
# Calculate the parameters of the sigmoidal model.
sigmoidalModel_green <- sicegar::parameterCalculation(sigmoidalModel_green)
# Generate plots of the data along with the fitted sigmoidal and double sigmoidal models.
fig01_green <- sicegar::figureModelCurves(dataInput = temp1_sm_green,
                                          sigmoidalFitVector = sigmoidalModel_green,
                                          showParameterRelatedLines = FALSE)
env_2_green <- fig01_green[["plot_env"]] 
# Access the list from the environment
fitted_model_2_green <- get("intensityTheoreticalSigmoidalDf", envir = env_2_green)%>% 
  filter(time<15)
raw_data_2_green <- get("dataFrameInput", envir = env_2_green)%>% 
  filter(time<15)

p_SM_green <- ggplot()+
  geom_point(data = raw_data_2_green, aes_(x = ~time, y = ~intensity), size=0.25) + 
  geom_line(data = fitted_model_2_green, 
            ggplot2::aes_(x = ~time, y = ~intensityTheoreticalSigmoidal), 
            color = "green2", size = 0.5)+
  labs(title= "", y="Late (NFI)",  x="Time (HPI)")+
  theme_void()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "",
        axis.text = element_blank(),
        axis.title.y = element_text(size = 5, angle=90, vjust = -3),
        axis.title.x = element_text(size = 5, vjust=1),
        text=element_text(family ="Arial"))+
  coord_cartesian(ylim=c(-10,180), xlim=c(-1, 16))+ #xlim for plot v2 only.
  annotate("segment", x = 0, xend = 16, y = -0.05, yend = -0.05, arrow = arrow(length = unit(0.125, "cm"),type="closed"), lwd=0.5) +
  annotate("segment", x = 0, xend = 0, y = 0, yend = 180, arrow = arrow(length = unit(0.125, "cm"),type="closed"), lwd=0.5)
ggsave("SM_green_linePlot.pdf",width=30, height=20, units="mm", p_SM_green)

#Panel E


p_SM <- ggplot()+
  geom_point(data = raw_data_2, aes_(x = ~time, y = ~intensity), size=0.5) + 
  geom_line(data = fitted_model_2, 
            ggplot2::aes_(x = ~time, y = ~intensityTheoreticalSigmoidal), 
            color = "red2", size = 1)+
  labs(title= "", y="Early (NFI)",  x="Time (HPI)")+
  theme_void()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "",
        axis.text = element_blank(),
        axis.title.y = element_text(size = 7, angle=90, vjust = -6),
        axis.title.x = element_text(size = 7, vjust=2),
        text=element_text(family ="Arial"))+
  coord_cartesian(ylim=c(-20,600), xlim=c(-1, 17))+ #xlim for plot v2 only.
  annotate("segment", x = 0, xend = 17, y = -0.05, yend = -0.05, arrow = arrow(length = unit(0.25, "cm"),type="closed"), lwd=1) +
  annotate("segment", x = 0, xend = 0, y = 0, yend = 600, arrow = arrow(length = unit(0.25, "cm"),type="closed"), lwd=1)
ggsave("SM_red_linePlot_2.pdf",width=60, height=40, units="mm", p_SM)

#Panel F
p_SM_green <- ggplot()+
  geom_point(data = raw_data_2_green, aes_(x = ~time, y = ~intensity), size=0.5) + 
  geom_line(data = fitted_model_2_green, 
            ggplot2::aes_(x = ~time, y = ~intensityTheoreticalSigmoidal), 
            color = "green2", size = 1)+
  labs(title= "", y="Late (NFI)",  x="Time (HPI)")+
  theme_void()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "",
        axis.text = element_blank(),
        axis.title.y = element_text(size = 7, angle=90, vjust = -6),
        axis.title.x = element_text(size = 7, vjust=2),
        text=element_text(family ="Arial"))+
  coord_cartesian(ylim=c(-10,180), xlim=c(-1, 16))+ #xlim for plot v2 only.
  annotate("segment", x = 0, xend = 16, y = -0.2, yend = -0.2, arrow = arrow(length = unit(0.25, "cm"),type="closed"), lwd=1) +
  annotate("segment", x = 0, xend = 0, y = -0.2, yend = 180, arrow = arrow(length = unit(0.25, "cm"),type="closed"), lwd=1)
ggsave("SM_green_linePlot_2.pdf",width=60, height=40, units="mm", p_SM_green)

#Panel G  
p_sicegar_nonProductive_max <- plot_data_with_sicegar_2_temp %>% 
  filter(Condition_red=="MOI1") %>% 
  dplyr:::filter(DEC_decision_red!="ambiguous"&DEC_decision_red!="no_signal") %>% 
  group_by(unique_ID) %>% 
  slice(1) %>% 
  group_by(productive) %>% 
  mutate(mean=mean(COMB_startPoint_x_red2), na.rm = TRUE) %>%
  mutate(sd=sd(COMB_startPoint_x_red2), na.rm = TRUE) %>% 
  ggplot(aes(x=productive, y=COMB_startPoint_x_red2))+
  geom_violin(aes(fill=productive),width=0.8, size=0.1)+
  #geom_signif(y_position = c(22), xmin = c(1), xmax = c(2),
   #           annotation = c("****"), size=0.15, textsize = 2, tip_length = 0.025)+
  scale_color_manual(values=c("grey","grey"))+
  scale_fill_manual(values=c("grey","grey"))+
  scale_x_discrete(limit = c("NO", "YES"),
                   labels = c("NP", "P"))+ 
  labs(title= "", y="Start Early (hr)",  x="")+
  theme_minimal()+  
  theme(legend.position = "none",
        axis.text = element_text(size = 5),
        axis.title = element_text(size = 7),
        text=element_text(family ="Arial"))+
  coord_cartesian(ylim=c(0,23))
ggsave("COMB_startPoint_x_red_productive_box_whisker.pdf",width=25, height=60, units="mm", p_sicegar_nonProductive_max)
#Panel H
# Make a copy of the filtered data
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


# Create the plot
p_sicegar_nonProductive_max <- plot_data_all %>%
  mutate(mean=mean(COMB_startDeclinePoint_x_red), na.rm = TRUE) %>%
  mutate(sd=sd(COMB_startDeclinePoint_x_red), na.rm = TRUE) %>% 
  ggplot(aes(x=productive, y=COMB_startDeclinePoint_x_red))+
  geom_violin(aes(fill=productive), width=0.8, size=0.1)+
  scale_color_manual(values=c("grey","grey","grey"))+
  scale_fill_manual(values=c("grey","grey","grey"))+
  scale_x_discrete(limit = c("NO", "YES", "COMBINED"),
                   labels = c("NP", "P", "All"))+ 
  labs(title= "", y="Lysis time (hr)",  x="")+
  theme_minimal()+  
  theme(legend.position = "none",
        axis.text = element_text(size = 5),
        axis.title = element_text(size = 7),
        text=element_text(family ="Arial"))+
  geom_signif(comparisons = list(c("NO", "YES")),
              size=0.15, textsize = 2,
              y_position = c(36),
              annotation = c("*"), 
              tip_length = 0.025)+
  coord_cartesian(ylim=c(0, 38))+
  scale_y_continuous(breaks=c(0, 10, 20, 30, 36))
ggsave("COMB_startDeclinePoint_x_redproductive_box_whisker.pdf",width=30, height=60, units="mm", p_sicegar_nonProductive_max)
#Panel I
unique_wells <- plot_data_all %>%
  dplyr:::filter(DEC_decision_red!="ambiguous"&DEC_decision_red!="no_signal") %>% 
  group_by(Metadata_Well_ID) %>% 
  dplyr::slice_head(n = 1) %>% 
  dplyr::select(Condition_red, Metadata_Well_ID) %>%
  distinct() %>%
  group_by(Condition_red) %>%
  dplyr::mutate(Well_Position_Within_Condition = paste("Well", row_number())) %>%
  ungroup()
# Shapes 21-25 can take a fill color
shapes <- c(21, 22, 23, 24, 25, 3)
p_lytic_frequency <- plot_data_all %>% 
  dplyr:::filter(DEC_decision_red!="ambiguous"&DEC_decision_red!="no_signal") %>% 
  left_join(unique_wells, by = c("Condition_red", "Metadata_Well_ID")) %>% 
  filter(Condition_red=="MOI1") %>% 
  dplyr:::group_by(productive, DEC_decision_red, FoV_ID) %>% 
  dplyr:::summarise(count=n_distinct(unique_ID),
                    Well_Position_Within_Condition=first(Well_Position_Within_Condition)) %>% 
  tidyr:::pivot_wider(names_from = DEC_decision_red, values_from = count) %>% 
  mutate(sigmoidal = coalesce(sigmoidal, 0)) %>%    #replace NAs with zeros
  group_by(productive, FoV_ID) %>% 
  dplyr:::summarise(frequency=double_sigmoidal/(sigmoidal+double_sigmoidal),
                    Well_Position_Within_Condition=first(Well_Position_Within_Condition)) %>% 
  ggplot(aes(y = frequency*100, x = productive, shape=Well_Position_Within_Condition)) +
  geom_point(color = "black", fill = "grey", position = position_jitterdodge(jitter.width = 0.7, dodge.width = 0), size = 1.2, stroke = 0.1) +
  scale_shape_manual(values = shapes) +
  labs(x = "", y = "E,L-to-lysis frequency (%)", shape = NULL) +
  scale_x_discrete(limit = c("NO", "YES", "COMBINED"),
                   labels = c("NP", "P", "All"))+ 
  theme_minimal() +
  theme(legend.position = "none",
        axis.text = element_text(size = 5),
        axis.title = element_text(size = 7),
        text=element_text(family ="Arial"))+
  guides(shape = guide_legend(override.aes = list(color = "black", fill = "black"))) +
  scale_y_continuous(breaks=c(0, 20, 40, 60, 80,100))+
  coord_cartesian(ylim = c(0, 110))
ggsave("lytic_frequency_MOI1.pdf", width=45, height=60, units="mm",  p_lytic_frequency)



# Figure 3 ----------------------------------------------------------------
# Data preparation
showtext_opts(dpi = 72) #Set font dpi
plot_data_combined <- plot_data_with_sicegar_2_temp %>%
  filter(Condition_red=="MOI1") %>% 
  dplyr:::filter(DEC_decision_red!="ambiguous"&DEC_decision_red!="no_signal") %>% 
  group_by(unique_ID) %>% 
  slice(1) %>%
  mutate(
    DEC_decision_red_2 = "All",
    unique_ID = sub("^[^_]*_", "20221007_", unique_ID) #Avoids duplication of unique IDs by replacing the date component with an arbitrary string
  ) 
# Combine the original and copied data
plot_data_with_sicegar_2_temp_fig3 <- rbind(plot_data_with_sicegar_2_temp %>%
                                                 dplyr:::filter(DEC_decision_red!="ambiguous"&DEC_decision_red!="no_signal") %>% 
                                                 mutate(DEC_decision_red_2=DEC_decision_red) %>% 
                                                 filter(Condition_red=="MOI1") %>%
                                                 group_by(unique_ID, DEC_decision_red_2) %>% 
                                                 slice(1), plot_data_combined)
# Select relevant columns and rename some of them
sicegar_scatterplot_filtered_green <- plot_data_with_sicegar_2_temp_fig3 %>% 
  # Filter the data where DEC_decision and DEC_decision_red are not "ambiguous" or "no_signal"
  dplyr:::filter(DEC_decision!="ambiguous"&DEC_decision!="no_signal") %>%
  dplyr:::filter(DEC_decision_red!="ambiguous"&DEC_decision_red!="no_signal") %>% 
  dplyr::select(unique_ID, Condition_red, DEC_decision_red, DEC_decision, COMB_maximum_y, DEC_decision, COMB_startPoint_x2, COMB_startPoint_x_red2,
                COMB_incrementTime, COMB_midPoint1_x, COMB_slope1) %>% 
  # Rename the selected columns to more easily interpretable names
  dplyr::rename("Maximum A3-GFP (NFI)"=COMB_maximum_y,
                "Productive Start (h)"=COMB_startPoint_x2,
                "Infection Time (h)"=COMB_incrementTime,
                "Mid Point (h)"=COMB_midPoint1_x,
                "Slope (NFI/h)"=COMB_slope1,
                "Infection Start (h)"=COMB_startPoint_x_red2) %>% 
  arrange(desc(DEC_decision_red_2))
#Filter to MOI 1 data
temp3<- sicegar_scatterplot_filtered_green %>% 
  filter(Condition_red=="MOI1")  
temp3$DEC_decision_red_2 <- fct_rev(factor(temp3$DEC_decision_red_2))
temp4<- sicegar_scatterplot_filtered_green %>% 
  filter(DEC_decision_red_2=="sigmoidal"|DEC_decision_red_2=="double_sigmoidal")

#Create a matrix combining pairwise correlation scatterplots and density distributions
# Define a function to create grouped scatter plots with correlation displayed on the plot
scatter_with_corr_smDSM <- function(data, mapping, colour = c("sigmoidal" ="darkorange1", "double_sigmoidal" = "steelblue2"), sizes = c("sigmoidal" = 1.3, "double_sigmoidal" = 1.1), alphas = c("sigmoidal" = 1, "double_sigmoidal" = 0.5), ...) {
  # Split data by Condition
  data_split <- split(data, data$DEC_decision_red_2)
  
  # Add size and alpha columns
  data$size <- sizes[data$DEC_decision_red_2]
  data$alpha <- alphas[data$DEC_decision_red_2]
  
  # Create an empty ggplot object
  p <- ggplot(data = data, mapping = mapping) + 
    geom_point(aes(color = DEC_decision_red_2, size = size/100, alpha = 1)) +
    scale_color_manual(values = colour) +
    scale_size_identity() +  # Identity scale to use the sizes directly
    scale_alpha_identity() +  # Identity scale to use the alphas directly
    theme_minimal() +
    theme(axis.text = element_blank(),
          axis.title = element_blank(),
          legend.position = "",
          panel.grid.minor = element_line(size = 0.125), panel.grid.major = element_line(size = 0.25))
  
  # Calculate correlation for the entire dataset
  x_all <- eval_data_col(data, mapping$x)
  y_all <- eval_data_col(data, mapping$y)
  correlation_all <- round(cor(x_all, y_all, use = "pairwise.complete.obs"), 2)
  correlation_3dp_all <- sprintf("%.2f", correlation_all)
  correlation_padded_all <- str_pad(correlation_3dp_all, width = 6, side = "left")
  
  label_all <- paste("r = ", correlation_padded_all)
  p <- p + annotate("text", x = Inf, y = Inf, label = label_all, 
                    hjust = 1, vjust = 1.2, size = 2, colour = "black", parse = FALSE)
  
  # Calculate correlation for each group and add to the plot
  for (i in seq_along(data_split)) {
    x <- eval_data_col(data_split[[i]], mapping$x)
    y <- eval_data_col(data_split[[i]], mapping$y)
    correlation <- round(cor(x, y, use = "pairwise.complete.obs"), 2)
    correlation_3dp <- sprintf("%.2f", correlation)
    # Pad the correlation value with spaces on the left
    correlation_padded <- str_pad(correlation_3dp, width = 6, side = "left")
    
    label <- paste("r = ", correlation_padded)
    p <- p + annotate("text", x = Inf, y = Inf, label = label, 
                      hjust = 1, vjust = 1.2 * (i + 1), size = 2, colour = colour[names(data_split)[[i]]], parse = FALSE)
  }
  
  p
}
#Define a function to create density plots
density_with_fake_axes <- function(data, mapping, colour = c("All"="grey40","sigmoidal" ="darkorange1", "double_sigmoidal" = "steelblue2"), ...) {
  
  # Create an empty ggplot object
  p <- ggplot(data = data, mapping = mapping) + 
    geom_density(aes(color = DEC_decision_red_2, fill=DEC_decision_red_2), adjust=4/5, alpha = 0.8) +
    scale_color_manual(values = colour) +
    scale_fill_manual(values = colour) +
    theme_minimal() +
    labs(y="Density")+
    theme(axis.title.x = element_blank(),
          legend.position = "",
          axis.text = element_text(size=5),
          axis.title = element_text(size=5),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(size = 0.25))
  
  # Calculate the max value for the x variable
  max_x <- max(eval_data_col(data, mapping$x), na.rm = TRUE)
  
  # Set the x-axis limits
  p <- p + scale_x_continuous(limits = c(0, max_x*1.4), labels = label_number(accuracy = 1))
  
  return(p)
}
# Define the variables of interest
variables <- c("Infection Start (h)", "Productive Start (h)", "Mid Point (h)", "Slope (NFI/h)", "Infection Time (h)","Maximum A3-GFP (NFI)")
# Create a data frame to store the plot objects
plot_matrix <- data.frame(matrix(vector((length(variables)^2), mode = "list"), ncol = length(variables)))
colnames(plot_matrix) <- variables
rownames(plot_matrix) <- variables
# Generate the upper triangle plots
for(i in 1:(length(variables)-1)) {
  for(j in (i+1):length(variables)) {
    plot_matrix[[variables[i], variables[j]]] <- ggplot() + theme_void()  # Empty plot
  }
}
# Generate the lower triangle plots
for(i in 2:length(variables)) {
  for(j in 1:(i-1)) {
    plot_matrix[[variables[i], variables[j]]] <- scatter_with_corr_smDSM(temp4[, c(variables, "DEC_decision_red_2")], ggplot2::aes_(x = as.name(variables[i]), y = as.name(variables[j])))
  }
}
# Generate the diagonal plots
for(i in 1:length(variables)) {
  plot_matrix[[variables[i], variables[i]]] <- density_with_fake_axes(temp3[, c(variables, "DEC_decision_red_2")], ggplot2::aes_(x = as.name(variables[i])))
}
# Create a list to store the grid of plots
plot_list <- list()
for(i in 1:nrow(plot_matrix)) {
  plot_row <- do.call("arrangeGrob", c(plot_matrix[i, ], ncol = ncol(plot_matrix)))
  plot_list[[i]] <- plot_row
}
labels_txt <- c("Start Early (hr)", "Start Late (hr)",  "Midpoint Late (hr)", " Slope Late (NFI/hr)", " Period Late (hr)","   Max Late (NFI)")
top_labels <- lapply(labels_txt, function(label) grid::textGrob(label, gp = grid::gpar(fontsize = 7), vjust=7, hjust=0.4))
right_labels <- lapply(labels_txt, function(label) grid::textGrob(label, gp = grid::gpar(fontsize = 7), rot = -90, vjust=7, hjust=0.6))
# Convert ggplot objects to grobs
grob_list <- unlist(lapply(plot_list, function(x) lapply(x$grobs, function(y) y[[1]])), recursive = FALSE) 
all_grobs <- c(top_labels, right_labels, grob_list)
layout_matrix <- matrix(NA, nrow = length(variables) + 1, ncol = length(variables) + 1)
layout_matrix[1, -ncol(layout_matrix)] <- 1:length(variables)  # top labels
layout_matrix[-1, ncol(layout_matrix)] <- (length(variables) + 1):(2 * length(variables))  # right labels
layout_matrix[-1, -ncol(layout_matrix)] <- (2 * length(variables) + 1):length(all_grobs)  # plots
# Arrange plots with grid.arrange()
p_lytic_matrix <- grid.arrange(grobs = all_grobs, layout_matrix = layout_matrix)
# Save the plot matrix to a tiff file
ggsave("lytic_vs_nonLytic_matrix.pdf", width=190, height=190, units="mm", plot = p_lytic_matrix)


# Figure 4 ----------------------------------------------------------------
#Panel A
showtext_opts(dpi = 300) #set font size
shapes <- c(21, 22, 23, 24, 25)
#Infection frequency
# Identify the unique wells within each condition
unique_wells <- infected_score %>%
  dplyr::select(Condition, Metadata_Well_ID) %>%
  distinct() %>%
  group_by(Condition) %>%
  dplyr::mutate(Well_Position_Within_Condition = as.character(row_number())) %>%
  ungroup()
#Plot the proportion of all cells that are infected
p_infection_proportion <- infected_score %>%
  filter(Condition != "uninfected") %>%
  left_join(unique_wells, by = c("Condition", "Metadata_Well_ID")) %>% 
  mutate(Condition = fct_relevel(Condition, "MOI100", "MOI50", "MOI10", "MOI1")) %>%
  ggplot(aes(y = proportion * 100, x = Condition, fill = Condition, shape = Well_Position_Within_Condition)) +
  geom_point(colour="black",position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), alpha=0.8, size = 1.2, stroke = 0.1) +
  # Add annotations for statistical significance of differences between Conditions
  geom_signif(y_position = c(105, 112, 119), xmin = c(3, 2, 1), xmax = c(4, 4, 4),
              annotation = c("****", "****", "****"), size=0.15, textsize = 2,tip_length = 0.0125, color="black") +
  scale_fill_manual(values = c("blue","magenta3", "cyan3", "darkgoldenrod2")) +
  labs(x = "", y = "Infected cells (%)", shape = NULL) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text = element_blank(),
        axis.title = element_blank(),
        text=element_text(family ="Arial"))+
  guides(shape = guide_legend(override.aes = list(color = "black", fill = "black"))) +
  coord_flip(clip = "off") +  # Turn off clipping
  scale_shape_manual(values = shapes) +
  scale_y_continuous(breaks = c(0,20,40,60,80,100), minor_breaks = NULL, limits = c(0, 120), expand = expansion(mult = c(0, .1)))
ggsave("infected_proportion.pdf", width=80, height=16, units="mm",  p_infection_proportion)
#Productive infection frequency
unique_wells <- productive_score %>%
  dplyr::select(Condition, Metadata_Well_ID) %>%
  distinct() %>%
  group_by(Condition) %>%
  dplyr::mutate(Well_Position_Within_Condition = paste("Well", row_number())) %>%
  ungroup()
#Plot the proportion of all infections that are productive
p_productive_proportion <- productive_score %>%
  filter(Condition != "uninfected") %>%
  left_join(unique_wells, by = c("Condition", "Metadata_Well_ID")) %>% 
  mutate(Condition = fct_relevel(Condition, "MOI100", "MOI50", "MOI10", "MOI1")) %>%
  ggplot(aes(y = proportion * 100, x = Condition, fill = Condition, shape = Well_Position_Within_Condition)) +
  geom_point(colour="black",position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), alpha=0.8, size = 1.2, stroke = 0.1) +
  # Add annotations for statistical significance of differences between Conditions
  geom_signif(y_position = c(105, 112, 119), xmin = c(3, 2, 1), xmax = c(4, 4, 4),
              annotation = c("****", "****", "****"), size=0.15, textsize = 2,tip_length = 0.0125, color="black") +
  scale_fill_manual(values = c("blue","magenta3", "cyan3", "darkgoldenrod2")) +
  labs(x = "", y = "Productive infections (%)", shape = NULL) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text = element_blank(),
        axis.title = element_blank(),
        text=element_text(family ="Arial"))+
  guides(shape = guide_legend(override.aes = list(color = "black", fill = "black"))) +
  coord_flip(clip = "off") +  # Turn off clipping
  scale_shape_manual(values = shapes) +
  scale_y_continuous(breaks = c(0,20,40,60,80,100), minor_breaks = NULL, limits = c(0, 120), expand = expansion(mult = c(0, .1))) 
ggsave("productive_proportion.pdf", width=80, height=16, units="mm",  p_productive_proportion)
#Frequency of lysis
# Identify the unique wells within each condition
unique_wells <- plot_data_with_sicegar_2_temp %>%
  dplyr:::filter(DEC_decision_red!="ambiguous"&DEC_decision_red!="no_signal") %>% 
  group_by(Metadata_Well_ID) %>% 
  slice_head(n = 1) %>% 
  dplyr::select(Condition_red, Metadata_Well_ID) %>%
  distinct() %>%
  group_by(Condition_red) %>%
  dplyr::mutate(Well_Position_Within_Condition = paste("Well", row_number())) %>%
  ungroup()

p_lytic_frequency <- plot_data_with_sicegar_2_temp %>% 
  dplyr:::filter(DEC_decision_red!="ambiguous"&DEC_decision_red!="no_signal") %>% 
  left_join(unique_wells, by = c("Condition_red", "Metadata_Well_ID")) %>% 
  dplyr:::group_by(Condition_red, DEC_decision_red, FoV_ID) %>% 
  dplyr:::summarise(count=n_distinct(unique_ID),
                    Well_Position_Within_Condition=first(Well_Position_Within_Condition)) %>% 
  tidyr:::pivot_wider(names_from = DEC_decision_red, values_from = count) %>% 
  mutate(sigmoidal = coalesce(sigmoidal, 0)) %>%    #replace NAs with zeros
  group_by(Condition_red) %>% 
  dplyr:::summarise(frequency=double_sigmoidal/(sigmoidal+double_sigmoidal),
                    Well_Position_Within_Condition=Well_Position_Within_Condition) %>% 
  mutate(Condition_red = fct_relevel(Condition_red, "MOI100", "MOI50", "MOI10", "MOI1")) %>%
  ggplot(aes(y = frequency*100, x = Condition_red, fill = Condition_red, shape = Well_Position_Within_Condition)) +
  geom_point(colour="black",position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), alpha=0.8, size = 1.2, stroke = 0.1) +
  # Add annotations for statistical significance of differences between Conditions
  #geom_signif(y_position = c(105, 112, 119), xmin = c(3, 2, 1), xmax = c(4, 4, 4),
   #           annotation = c("****", "****", "****"), size=0.15, textsize = 2,tip_length = 0.0125, color="black") +
  scale_fill_manual(values = c("blue","magenta3", "cyan3", "darkgoldenrod2")) +
  labs(x = "", y = "Frequency (%)", shape = NULL) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 5),
        axis.title.x = element_text(size = 7),
        axis.text.y = element_blank(),
        text=element_text(family ="Arial")) +
  guides(shape = guide_legend(override.aes = list(color = "black", fill = "black"))) +
  coord_flip(clip = "off") +  # Turn off clipping
  scale_shape_manual(values = shapes) +
  scale_y_continuous(breaks = c(0,20,40,60,80,100), minor_breaks = NULL, limits = c(0, 120), expand = expansion(mult = c(0, .1)))
ggsave("lytic_proportion.pdf", width=80, height=24, units="mm", p_lytic_frequency)

#Panel B
p_red_green_lag <- RED_GREEN_LAG_stats_plots_MOI %>% 
  mutate(Condition = fct_relevel(Condition, "MOI1", "MOI10", "MOI50", "MOI100")) %>%
  ggplot(aes(x=Condition, y=mean_lag, fill=Condition)) + 
  geom_violin(size=0.1)+
  geom_signif(y_position = c(24, 27, 30, 33), xmin = c(1, 1, 2, 2), xmax = c(2, 3, 3, 4),
              annotation = c("****", "****", "****", "*"), size=0.15, textsize = 2, tip_length = 0.0125, color="black") +
  labs(title= "", y="E-to-L delay (hr)",  x="")+
  theme_minimal()+
  scale_color_manual(values=c("darkgoldenrod2","cyan3","magenta3", "blue"))+
  scale_fill_manual(values=c("darkgoldenrod2","cyan3","magenta3", "blue"))+
  theme(legend.position = "none",
        axis.text = element_text(size = 5),
        axis.title = element_text(size = 7),
        axis.text.x = element_blank(),
        text=element_text(family ="Arial"))+
  scale_y_continuous(breaks = c(0, 5, 10, 15, 20, 25, 30))+
  coord_cartesian(ylim=c(0,34))
ggsave("red_green_lag.pdf", width=25, height=60, units="mm", p_red_green_lag)

#Panel D
p_red_green_lag <- plot_data_with_sicegar_2_temp %>% 
  dplyr:::filter(DEC_decision_red!="ambiguous"&DEC_decision_red!="no_signal") %>% 
  group_by(unique_ID) %>% 
  #filter(productive=="NO") %>% 
  slice(1) %>% 
  filter(COMB_slope1_red<348) %>% 
  mutate(Condition.x = fct_relevel(Condition.x, "MOI1", "MOI10", "MOI50", "MOI100")) %>%
  ggplot(aes(y=COMB_slope1_red, x=Condition_red, fill=Condition_red))+
  geom_violin(size=0.1)+
  geom_signif(y_position = c(290, 360, 400, 440, 480), xmin = c(1, 1, 1, 2, 2), xmax = c(2, 3, 4, 3, 4),
              annotation = c("****", "****", "****", "*", "*"), size=0.15, textsize = 2, tip_length = 0.025, color="black") +
  labs(title= "", y="Slope Early (NFI/hr)",  x="")+
  theme_minimal()+
  scale_color_manual(values=c("darkgoldenrod2","cyan3","magenta3", "blue"))+
  scale_fill_manual(values=c("darkgoldenrod2","cyan3","magenta3", "blue"))+
  theme(legend.position = "none",
        axis.text = element_text(size = 5),
        axis.title = element_text(size = 7),
        axis.text.x = element_blank(),
        text=element_text(family ="Arial"))+
  coord_cartesian(ylim=c(0,490))
ggsave("Slope_Early_violin.pdf", width=25, height=60, units="mm", p_red_green_lag)

#Panel C
showtext_opts(dpi = 72) #set font size
#Filter data to lytici nfections
sicegar_scatterplot_filtered_green <- plot_data_with_sicegar_2_temp %>% 
  group_by(unique_ID) %>% 
  slice(1) %>% 
  # Filter the data where DEC_decision and DEC_decision_red are not "ambiguous" or "no_signal"
  dplyr:::filter(DEC_decision!="ambiguous"&DEC_decision!="no_signal") %>%
  dplyr:::filter(DEC_decision_red!="ambiguous"&DEC_decision_red!="no_signal") %>% 
  dplyr::select(unique_ID, Condition.x, DEC_decision_red, DEC_decision, COMB_maximum_y, DEC_decision, COMB_startPoint_x2, COMB_startPoint_x_red2,
                COMB_incrementTime, COMB_midPoint1_x, COMB_slope1) %>% 
  # Rename the selected columns to more easily interpretable names
  dplyr::rename("Maximum A3-GFP (NFI)"=COMB_maximum_y,
                "Productive Start (h)"=COMB_startPoint_x2,
                "Infection Time (h)"=COMB_incrementTime,
                "Mid Point (h)"=COMB_midPoint1_x,
                "Slope (NFI/h)"=COMB_slope1,
                "Infection Start (h)"=COMB_startPoint_x_red2) %>% 
  arrange(desc(Condition.x))
temp2<- sicegar_scatterplot_filtered_green %>% 
  mutate(Condition.x = fct_relevel(Condition.x, "MOI1", "MOI10", "MOI50", "MOI100")) 

#Scatterplot function
scatter_with_corr_MOI <- function(data, mapping, colour = c("MOI1" ="darkgoldenrod2", "MOI10" = "cyan3", "MOI50" = "magenta3", "MOI100" = "blue"), ...) {
  # Split data by Condition
  data_split <- split(data, data$Condition.x)
  
  # Create an empty ggplot object
  p <- ggplot(data = data, mapping = mapping) + 
    geom_point(aes(color = Condition.x), alpha = 0.5, size = 0.15, stroke=0.2) +
    scale_color_manual(values = colour) +
    scale_fill_manual(values = colour) +
    theme_minimal() +
    theme(axis.text = element_blank(),
          axis.title = element_blank(),
          legend.position = "",
          panel.grid.minor = element_line(size = 0.125), panel.grid.major = element_line(size = 0.25))
  
  # Calculate correlation for each group and add to the plot
  for (i in seq_along(data_split)) {
    x <- eval_data_col(data_split[[i]], mapping$x)
    y <- eval_data_col(data_split[[i]], mapping$y)
    correlation <- round(cor(x, y, use = "pairwise.complete.obs", method="pearson"), 2)
    correlation_3dp <- sprintf("%.2f", correlation)
    # Pad the correlation value with spaces on the left
    correlation_padded <- str_pad(correlation_3dp, width = 6, side = "left")
    
    label <- paste("r = ", correlation_padded)
    p <- p + annotate("text", x = Inf, y = Inf, label = label, 
                      hjust = 1, vjust = 1.2 * i, size = 2, colour = colour[names(data_split)[[i]]], parse = FALSE)
  }
  
  p
}
# Define a function to create density plots
density_with_fake_axes <- function(data, mapping, colour = c("MOI1" = "darkgoldenrod2", "MOI10" = "cyan3", "MOI50" = "magenta3", "MOI100" = "blue"), ...) {
  # Convert condition.x to a factor and reverse the levels
  data$condition.x <- fct_rev(factor(data$Condition.x))
  # Create an empty ggplot object
  p <- ggplot(data = data, mapping = mapping) + 
    geom_density(aes(color = condition.x, fill=Condition.x), size=0.2, adjust=4/5, alpha = 0.7) +
    scale_color_manual(values = colour) +
    scale_fill_manual(values = colour) +
    labs(y="Density")+
    theme_minimal() +
    theme(axis.title.x = element_blank(),
          legend.position = "",
          axis.text = element_text(size=5),
          axis.title.y = element_text(size=5),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(size = 0.25))
  
  # Calculate the max value for the x variable
  max_x <- max(eval_data_col(data, mapping$x), na.rm = TRUE)
  
  # Set the x-axis limits
  p <- p + scale_x_continuous(limits = c(0, max_x*1.2), labels = label_number(accuracy = 1))
  
  return(p)
}
# Define the variables of interest
variables <- c("Infection Start (h)", "Productive Start (h)", "Mid Point (h)", "Slope (NFI/h)", "Infection Time (h)","Maximum A3-GFP (NFI)")
# Create a data frame to store the plot objects
plot_matrix <- data.frame(matrix(vector((length(variables)^2), mode = "list"), ncol = length(variables)))
colnames(plot_matrix) <- variables
rownames(plot_matrix) <- variables
# Generate the upper triangle plots
for(i in 1:(length(variables)-1)) {
  for(j in (i+1):length(variables)) {
    plot_matrix[[variables[i], variables[j]]] <- ggplot() + theme_void()  # Empty plot
  }
}
# Generate the lower triangle plots
for(i in 2:length(variables)) {
  for(j in 1:(i-1)) {
    plot_matrix[[variables[i], variables[j]]] <- scatter_with_corr_MOI(temp2[, c(variables, "Condition.x")], ggplot2::aes_(x = as.name(variables[i]), y = as.name(variables[j])))
  }
}
# Generate the diagonal plots
for(i in 1:length(variables)) {
  plot_matrix[[variables[i], variables[i]]] <- density_with_fake_axes(temp2[, c(variables, "Condition.x")], ggplot2::aes_(x = as.name(variables[i])))
}
# Create a list to store the grid of plots
plot_list <- list()
for(i in 1:nrow(plot_matrix)) {
  plot_row <- do.call("arrangeGrob", c(plot_matrix[i, ], ncol = ncol(plot_matrix)))
  plot_list[[i]] <- plot_row
}
labels_txt <- c("Start Early (hr)", "Start Late (hr)",  "Midpoint Late (hr)", " Slope Late (NFI/hr)", " Period Late (hr)"," Max Late (NFI)")
top_labels <- lapply(labels_txt, function(label) grid::textGrob(label, gp = grid::gpar(fontsize = 7), vjust=7, hjust=0.4))
right_labels <- lapply(labels_txt, function(label) grid::textGrob(label, gp = grid::gpar(fontsize = 7), rot = -90, vjust=7, hjust=0.6))
# Convert ggplot objects to grobs
grob_list <- unlist(lapply(plot_list, function(x) lapply(x$grobs, function(y) y[[1]])), recursive = FALSE) 
all_grobs <- c(top_labels, right_labels, grob_list)
layout_matrix <- matrix(NA, nrow = length(variables) + 1, ncol = length(variables) + 1)
layout_matrix[1, -ncol(layout_matrix)] <- 1:length(variables)  # top labels
layout_matrix[-1, ncol(layout_matrix)] <- (length(variables) + 1):(2 * length(variables))  # right labels
layout_matrix[-1, -ncol(layout_matrix)] <- (2 * length(variables) + 1):length(all_grobs)  # plots
# Arrange plots with grid.arrange()
p_multipleMOI_matrix <- grid.arrange(grobs = all_grobs, layout_matrix = layout_matrix)
# Save the plot matrix to a tiff file
ggsave("multipleMOI_matrix.pdf", width=170, height=170, units="mm", plot = p_multipleMOI_matrix)

#Panel E
showtext_opts(dpi = 300) #set font size
#Specify df

df2<-plot_data_with_sicegar_2_temp %>% 
  group_by(unique_ID) %>%
  slice(1) %>% 
  ungroup() %>% 
  filter(unique_ID!="20230829_3_2_37") #outlier
set.seed(10)

# Replace NA in numeric columns with 0
df2 <- df2 %>%
  mutate_if(is.numeric, ~replace_na(., 0))

# Replace NA in character columns with ""
df2 <- df2 %>%
  mutate_if(is.character, ~replace_na(., ""))

MDS_matrix<-df2 %>% 
  #filter(DEC_decision_red=="double_sigmoidal")  %>% 
  filter(DEC_decision=="sigmoidal")  %>% 
  #dplyr:::select(!c(DEC_decision, unique_ID, Condition.x, DEC_decision_red, Condition_red)) %>% 
  #dplyr::select(COMB_slope1, COMB_maximum_y, COMB_incrementTime, 
   #         COMB_startPoint_x2, COMB_startPoint_x_red2, COMB_midPoint1_x) %>% #all parameters
  dplyr::select(COMB_startPoint_x2, COMB_startPoint_x_red2, COMB_midPoint1_x) %>% #Timing parameters
  #dplyr::select(COMB_slope1, COMB_maximum_y, COMB_incrementTime) %>% #Production parameters
  scale()
group_labels <- df2 %>% 
  #filter(DEC_decision_red=="double_sigmoidal")  %>% 
  dplyr::filter(DEC_decision=="sigmoidal")  %>% 
  mutate(factor_unique_ID=factor(unique_ID)) %>% 
  dplyr::select(unique_ID)
# Calculate distance matrix using Euclidean distance
dist_matrix <- dist(MDS_matrix, method = "manhattan")
# Perform classical multidimensional scaling (MDS)
mds_coords <- cmdscale(dist_matrix, k = 2) # k is the number of dimensions
group_labels_temp<-df2 %>%
  #filter(DEC_decision_red=="double_sigmoidal")  %>% 
  filter(DEC_decision=="sigmoidal")  %>% 
  dplyr::select(unique_ID,Condition.x,DEC_decision,DEC_decision_red) 
mds_df <- data.frame(mds_coords, unique_ID = group_labels) %>% 
  left_join(., group_labels_temp, by="unique_ID")

# Create the MDS plot using ggplot2
p<-mds_df %>% 
  ggplot(aes(x = X1, y = X2)) +
  geom_point(aes(colour = Condition.x),size=0.5, alpha=0.5, stroke = 0.1) +
  labs(title = "", x = "Dimension 2", y = "Dimension 1", color="") +
  theme_minimal()+
  theme(legend.position = "",
        axis.text = element_blank(),
        axis.title = element_text(size = 5)) +
  guides(colour = guide_legend(override.aes = list(size=10)))+
  scale_color_manual(values=c("darkgoldenrod2","cyan3","magenta3","blue"))
ggsave("MDS_all_factors.pdf", width=40, height=40, unit="mm", p)
#Panel E
ggsave("MDS_host_factors.pdf", width=40, height=40, unit="mm", p)
#Panel F
ggsave("MDS_viral_factors.pdf", width=40, height=40, unit="mm", p)

# Figure 5 ----------------------------------------------------------------

#Panel A
showtext_opts(dpi = 300) #set font size
# Shapes 21-25 can take a fill color
shapes <- c(21, 22, 23, 24, 25, 3)
#Infection frequency
# Determine the unique wells within each condition
unique_wells <- infected_score_drugs %>%
  filter(Condition=="Infected"|Condition=="QVD"|Condition=="QVD_Nec") %>% 
  group_by(Date_Well_ID) %>% 
  dplyr::select(Condition, Date_Well_ID) %>%
  distinct() %>%
  group_by(Condition) %>%
  dplyr::mutate(Well_Position_Within_Condition = paste("Well", row_number())) %>%
  ungroup()

#Plot frequency of exposed cells that become infected
p_infection_proportion <- infected_score_drugs %>% 
  filter(Condition=="Infected"|Condition=="QVD"|Condition=="QVD_Nec") %>% 
  left_join(unique_wells, by = c("Condition", "Date_Well_ID")) %>% 
  mutate(Condition = fct_rev(Condition)) %>%
  ggplot(aes(y = proportion * 100, x = Condition, fill = Condition, shape = Well_Position_Within_Condition)) +
  geom_point(colour="black",position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), size = 1.2, stroke = 0.1) +
  # Add annotations for statistical significance of differences between Conditions
  geom_signif(y_position = c(105, 110), xmin = c(1,1), xmax = c(2,3),
              annotation = c("***", "****"), size=0.15, textsize = 2, tip_length = 0.025)+
  scale_fill_manual(values=c("red3","blue2","green2"))+
  labs(x = "", y = "Infected cells (%)", shape = NULL) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text = element_blank(),
        axis.title = element_blank(),
        text=element_text(family ="Arial"))+
  guides(shape = guide_legend(override.aes = list(color = "black", fill = "black"))) +
  coord_flip(clip = "off") +  # Turn off clipping
  scale_shape_manual(values = shapes) +
  scale_y_continuous(breaks = c(0,20,40,60,80,100), minor_breaks = NULL, limits = c(0, 110), expand = expansion(mult = c(0, .1))) 
ggsave("infected_proportion_drugs.pdf", width=80, height=16, units="mm",  p_infection_proportion)
# Frequency of infections that become productive
p_productive_proportion <- productive_score_drugs %>% 
  filter(Condition=="Infected"|Condition=="QVD"|Condition=="QVD_Nec") %>% 
  left_join(unique_wells, by = c("Condition", "Date_Well_ID")) %>% 
  mutate(Condition = fct_rev(Condition)) %>%
  ggplot(aes(y = proportion * 100, x = Condition, fill = Condition, shape = Well_Position_Within_Condition)) +
  geom_point(colour="black",position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), size = 1.2, stroke = 0.1) +
  # Add annotations for statistical significance of differences between Conditions
  geom_signif(y_position = c(105), xmin = c(1), xmax = c(3),
              annotation = c("*"), size=0.15, textsize = 2, tip_length = 0.025)+
  scale_fill_manual(values=c("red3","blue2","green2"))+
  labs(x = "", y = "Productive cells (%)", shape = NULL) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text = element_blank(),
        axis.title = element_blank(),
        text=element_text(family ="Arial"))+
  guides(shape = guide_legend(override.aes = list(color = "black", fill = "black"))) +
  coord_flip(clip = "off") +  # Turn off clipping
  scale_shape_manual(values = shapes) +
  scale_y_continuous(breaks = c(0,20,40,60,80,100), minor_breaks = NULL, limits = c(0, 110), expand = expansion(mult = c(0, .1)))
ggsave("productive_proportion_drugs.pdf", width=80, height=16, units="mm",  p_productive_proportion)
#Lysis frequency
unique_wells <- plot_data_with_sicegar_drugs %>%
  dplyr:::filter(DEC_decision_red!="ambiguous"&DEC_decision_red!="no_signal") %>% 
  group_by(Metadata_Well_ID) %>% 
  slice_head(n = 1) %>% 
  dplyr::select(Condition_red, Metadata_Well_ID) %>%
  distinct() %>%
  group_by(Condition_red) %>%
  dplyr::mutate(Well_Position_Within_Condition = paste("Well", row_number())) %>%
  ungroup()
p_lytic_frequency_drugs <- plot_data_with_sicegar_drugs %>% 
  dplyr:::filter(DEC_decision_red!="ambiguous"&DEC_decision_red!="no_signal") %>% 
  left_join(unique_wells, by = c("Condition_red", "Metadata_Well_ID")) %>% 
  filter(Condition_red=="Infected"|Condition_red=="QVD"|Condition_red=="QVD_Nec") %>% 
  dplyr:::group_by(Condition_red, DEC_decision_red, Metadata_Well_ID) %>% 
  dplyr:::summarise(count=n_distinct(unique_ID),
                    Well_Position_Within_Condition=first(Well_Position_Within_Condition)) %>% 
  tidyr:::pivot_wider(names_from = DEC_decision_red, values_from = count) %>% 
  mutate(sigmoidal = coalesce(sigmoidal, 0)) %>%    #replace NAs with zeros
  group_by(Condition_red, Metadata_Well_ID) %>% 
  dplyr:::summarise(frequency=double_sigmoidal/(sigmoidal+double_sigmoidal),
                    Well_Position_Within_Condition=first(Well_Position_Within_Condition)) %>% 
  mutate(frequency = coalesce(frequency, 0)) %>%
  mutate(Condition_red = fct_rev(Condition_red)) %>%
  ggplot(aes(y = frequency * 100, x = Condition_red, fill = Condition_red, shape = Well_Position_Within_Condition)) +
  geom_point(colour="black",position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), size = 1.2, stroke = 0.1) +
  # Add annotations for statistical significance of differences between Conditions
  geom_signif(y_position = c(105, 110), xmin = c(1, 1), xmax = c(2,3 ),
              annotation = c("****","****"), size=0.15, textsize = 2,tip_length = 0.025)+
  scale_fill_manual(values=c("red3","blue2","green2"))+
  labs(x = "", y = "Frequency (%)", shape = NULL) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text = element_text(size = 5),
        axis.title = element_text(size = 7),
        axis.text.y = element_blank(),
        text=element_text(family ="Arial"))+
  guides(shape = guide_legend(override.aes = list(color = "black", fill = "black"))) +
  coord_flip(clip = "off") +  # Turn off clipping
  scale_shape_manual(values = shapes) +
  scale_y_continuous(breaks = c(0,20,40,60,80,100), minor_breaks = NULL, limits = c(0, 110), expand = expansion(mult = c(0, .1)))
ggsave("lytic_proportion_drugs.pdf", width=80, height=24, units="mm",  p_lytic_frequency_drugs)

#Panel B
showtext_opts(dpi = 300) #set font size
p_red_green_lag_drugs <- RED_GREEN_LAG_stats_plots_drugs %>% 
  filter(Condition=="Infected"|Condition=="QVD") %>% 
  ggplot(aes(x=Condition, y=mean_lag, fill=Condition)) + 
  geom_violin(size=0.3)+
  #geom_signif(y_position = c(20), xmin = c(1), xmax = c(2),
   #           annotation = c("ns"), tip_length = 0.025)+
  labs(title= "", y="E-to-L delay (hr)",  x="")+
  theme_minimal()+
  scale_color_manual(values=c("green2","blue2"))+
  scale_fill_manual(values=c("green2","blue2"))+
  theme(legend.position = "none",
        axis.text = element_text(size = 5),
        axis.title = element_text(size = 7),
        axis.text.x = element_blank(),
        text=element_text(family ="Arial"))
ggsave("red_green_lag_drugs.pdf", width=40, height=60, units="mm",  p_red_green_lag_drugs)

#Panel C
showtext_opts(dpi = 300) #set font size
# Make a copy of the filtered data to generate QVD, all Control and Control non-lytic groups
# Creating the grouping variable
plot_data_combined_drugs <- plot_data_with_sicegar_drugs_2 %>%
  filter(Condition_red=="Infected") %>% 
  dplyr:::filter(DEC_decision_red!="ambiguous"&DEC_decision_red!="no_signal") %>% 
  group_by(unique_ID) %>% 
  slice(1) %>%
  mutate(
    decision_Condition = "Control",
    unique_ID = sub("^[^_]*_", "20221007_", unique_ID)
  ) 
# Combine the original and copied data
plot_data_with_sicegar_drugs_2_fig6_7 <- rbind(plot_data_with_sicegar_drugs_2 %>%
                                                 dplyr:::filter(DEC_decision_red!="ambiguous"&DEC_decision_red!="no_signal") %>% 
                                                 group_by(unique_ID, decision_Condition) %>% 
                                                 slice(1), plot_data_combined_drugs)
sicegar_scatterplot_filtered_green_drugs <- plot_data_with_sicegar_drugs_2_fig6_7 %>% 
  # Filter the data where DEC_decision and DEC_decision_red are not "ambiguous" or "no_signal"
  dplyr:::filter(DEC_decision!="ambiguous"&DEC_decision!="no_signal") %>%
  dplyr:::filter(DEC_decision_red!="ambiguous"&DEC_decision_red!="no_signal") %>% 
  filter(unique_ID!="20221010_2_0_14") %>%  #extreme outlier
  filter(decision_Condition=="Infected_sigmoidal"|decision_Condition=="Control"|decision_Condition=="QVD") %>% 
  dplyr::select(unique_ID, Condition_red, DEC_decision_red, DEC_decision, decision_Condition, COMB_maximum_y, DEC_decision, COMB_startPoint_x, COMB_startPoint_x_red,
                COMB_incrementTime, COMB_midPoint1_x, COMB_slope1)

#re-order levels of decision_Condition
levels_order <- c("QVD", "Control", "Infected_sigmoidal")
#Re-order levels in grouping variable for aesthetic purposes
sicegar_scatterplot_filtered_green_drugs_2$decision_Condition <- factor(sicegar_scatterplot_filtered_green_drugs_2$decision_Condition, levels = levels_order)
# Create a list of variables and corresponding y-labels
variable_list <- c("COMB_startPoint_x_red", "COMB_startPoint_x", "COMB_midPoint1_x", 
                   "COMB_slope1", "COMB_incrementTime", "COMB_maximum_y")
label_list <- c("Start Early (hr)", "Start Late (hr)",  "Midpoint Late (hr)", 
                "Slope Late (NFI/hr)", "Period Late (hr)","Max Late (NFI)")
# Define the function that creates and saves the plots
create_save_plot <- function(data, variable, ylabel) {
  p <- ggplot(data = data) + 
    geom_density(aes_string(x=variable, color = "decision_Condition", fill="decision_Condition"), alpha = 0.6) +
    scale_color_manual(values = c("Control" = "green2", "Infected_sigmoidal" ="darkorange1","QVD" = "blue2")) +
    scale_fill_manual(values = c("Control" = "green2", "Infected_sigmoidal" ="darkorange1","QVD" = "blue2")) +
    theme_minimal() +
    labs(x = ylabel, y = "") + 
    theme(legend.position = "none",
          axis.text = element_text(size = 5),
          axis.title.y = element_blank(),
          axis.title.x=element_text(size = 7),
          text=element_text(family ="Arial"))
  
  # Calculate the max value for the x variable
  max_x <- max(data[[variable]], na.rm = TRUE)
  
  # Check that max_x is finite
  if(is.finite(max_x)) {
    # Set the x-axis limits
    p <- p + scale_x_continuous(limits = c(0, max_x*1.4), labels = label_number(accuracy = 1))
  }
  
  # Save the plot
  ggsave(paste0(variable, "_drugs.pdf"), width=29, height=40, units="mm",plot=p)
}
# Loop through the variable list and generate the plots
for (i in 1:length(variable_list)) {
  create_save_plot(sicegar_scatterplot_filtered_green_drugs_2, variable_list[i], label_list[i])
}

#Panel E
showtext_opts(dpi = 300) #set font size
#Read in data
data_size <- read.csv(file="data_size.csv")
#Set scale in micron/px
scale <- 1.62
#Change measurements to microns
data_size_worked<- data_size %>% 
  unite("unique_plaque_ID",?..OG_Object_ID, Well_ID, Condition, remove=FALSE) %>% 
  mutate(Area_micron=Area*1.62,
         Feret_micron=Feret*1.62)

#Plot plaque diameter
size_plot <- data_size_worked %>%
  mutate(Condition = fct_rev(Condition)) %>%
  ggplot(aes(y = Feret_micron/1000, x = Condition, fill = Condition, shape = factor(Well_ID))) +
  geom_point(color="black",position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.5), alpha=0.8, size = 1.2, stroke = 0.1) +
  # Add annotations for statistical significance of differences between Conditions
  geom_signif(y_position = c(0.95), xmin = c(1), xmax = c(2),
              annotation = c("****"), size=0.15, textsize = 2,tip_length = 0.025, color="black") +
  scale_fill_manual(values = c("blue2","green2")) +
  labs(x = "", y = "Plaque diameter (mm)") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text = element_text(size = 5),
        axis.title = element_text(size = 7),
        axis.text.y = element_blank(),
        text=element_text(family ="Arial")) +
  guides(shape = guide_legend(override.aes = list(color = "black", fill = "black"))) +
  coord_flip(clip = "off") +  # Turn off clipping
  scale_shape_manual(values = shapes) +
  scale_y_continuous(breaks = c(0,0.25,0.50,0.75), minor_breaks = NULL, limits = c(0, 1), expand = expansion(mult = c(0, .1)))
ggsave("plaque_diameter.pdf", width=80, height=40, units="mm",  size_plot)

#Panel F
showtext_opts(dpi = 300) #set font size
#Read in data
data_count <- read.csv(file="data.csv")
#Plot plaque count
count_plot <- data_count_2 %>%
  mutate(Condition = fct_rev(Condition)) %>%
  ggplot(aes(y = Count_2, x = Condition, fill = Condition, shape = factor(Well_ID))) +
  geom_point(color="black",position = position_jitterdodge(jitter.width = 0, dodge.width = 0), alpha=1, size = 3, stroke = 0.1) +
  # Add annotations for statistical significance of differences between Conditions
  geom_signif(y_position = c(57), xmin = c(1), xmax = c(2),
              annotation = c("*"),size=0.15, textsize = 2, tip_length = 0.025, color="black") +
  scale_fill_manual(values = c("blue2","green2")) +
  labs(x = "", y = "Plaque count") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text = element_text(size = 5),
        axis.title = element_text(size = 7),
        axis.text.y = element_blank(),
        text=element_text(family ="Arial")) +
  guides(shape = guide_legend(override.aes = list(color = "black", fill = "black"))) +
  coord_flip(clip = "off") +  # Turn off clipping
  scale_shape_manual(values = shapes) +
  scale_y_continuous(breaks = c(30,40,50,60), minor_breaks = NULL, limits = c(30, 60), expand = expansion(mult = c(0, .1))) 
ggsave("plaque_count.pdf", width=50, height=40, units="mm",  count_plot)



#Figure 6 ####
#Panel A
VACV_means <- data.frame(
  MOI = c(1, 10, 50, 100),
  start_early = c(1,
                  0.424922127,
                  0.265584663,
                  0.231102808),
  start_late = c(1,
                 0.506102113,
                 0.333038441,
                 0.287729501),
  mid = c(1,
          0.600033248,
          0.49122536,
          0.450739164),
  slope = c(0.792215167,
            0.893619352,
            1,
            0.993461507),
  period = c(0.902239486,
             0.951869415,
             1,
             0.928563115),
  max = c(0.784180873,
          0.853597772,
          1,
          0.93735775) 
) %>% 
  pivot_longer(cols=2:7, values_to = "relative_value", names_to = "Variable")
# Define the levels of Variable and their corresponding colors

custom_colors_2_2 <- c("red","green2", "yellow3","cyan3", "magenta3","blue")
levels_custom <- c("start_early","start_late", "mid", "slope", "period","max")
colors_custom <- c("start_early" = custom_colors_2_2[1], 
                   "start_late" = custom_colors_2_2[2],
                   "mid" = custom_colors_2_2[3], 
                   "slope" = custom_colors_2_2[4],
                   "period" = custom_colors_2_2[5],
                   "max" = custom_colors_2_2[6])
labels_displayed <- c("Start Early", "Start Late","Midpoint Late", "Slope Late", "Period Late", "Max Late")
VACV_means$Variable <- factor(VACV_means$Variable, levels = levels_custom)

# Plotting
plot <- ggplot() +
  geom_line(data = VACV_means, aes(x = MOI, y = relative_value, group=Variable, color = Variable), size=1) +
  geom_point(data=VACV_means, aes(x = MOI, y = relative_value, fill = Variable, color=Variable), size=2, stroke=0.1) +
  scale_color_manual(values = colors_custom, guide = "none") +  
  scale_fill_manual(values = colors_custom, breaks = levels_custom, labels = labels_displayed, name = "Parameter")+
  labs(title = "",
       y = "Relative Value",
       x = "MOI") +
  theme_minimal() +
  theme(legend.position = "right",
        legend.key.size = unit(2, "mm"),          
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 5),
        axis.text = element_text(size = 5),
        axis.title.y = element_text(size = 7, angle=90),
        axis.title.x = element_text(size = 7),
        text=element_text(family ="Arial")) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) # Override shape in legend

ggsave("VACV_MEANS_V2.pdf", width=70, height=55, units="mm", plot)

#Panel B
PV_means <- data.frame(
  MOI = c(0.045,
          0.45,
          4.5),
  start = c(1,
            0.860674157,
            0.637078652),
  mid = c(1,
          0.835304823,
          0.621474067),
  slope = c(0.511111111,
            0.7,
            1),
  period = c(1,
             0.727923628,
             0.556085919),
  max = c(0.867088608,
          0.905063291,
          1) 
) %>% 
  pivot_longer(cols=2:6, values_to = "relative_value", names_to = "Variable")

custom_colors_2_2 <- c("red", "yellow3","cyan3", "magenta3","blue")
levels_custom <- c("start", "mid", "slope", "period","max")
colors_custom <- c("start" = custom_colors_2_2[1], 
                   "mid" = custom_colors_2_2[2], 
                   "slope" = custom_colors_2_2[3],
                   "period" = custom_colors_2_2[4],
                   "max" = custom_colors_2_2[5])
labels_displayed <- c("Start", "Midpoint", "Slope", "Period", "Max")
PV_means$Variable <- factor(PV_means$Variable, levels = levels_custom)

# Plotting
plot <- ggplot() +
  geom_line(data = PV_means, aes(x = MOI, y = relative_value, group=Variable, color = Variable), size=1) +
  geom_point(data=PV_means, aes(x = MOI, y = relative_value, fill = Variable, color=Variable), size=2, stroke=0.1) +
  scale_color_manual(values = colors_custom, guide = "none") +  
  scale_fill_manual(values = colors_custom, breaks=levels_custom, labels = labels_displayed, name = "Parameter") +
  labs(title = "",
       y = "Relative Value",
       x = "MOI") +
  theme_minimal() +
  theme(legend.position = "right",
        legend.key.size = unit(2, "mm"),          
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 5),
        axis.text = element_text(size = 5),
        axis.title.y = element_text(size = 7, angle=90),
        axis.title.x = element_text(size = 7),
        text=element_text(family ="Arial")) +
  guides(fill = guide_legend(override.aes = list(shape = 21)))+
  coord_cartesian(ylim=c(0.2,1))

ggsave("PV_MEANS_V2.pdf", width=70, height=55, units="mm", plot)

#Panel D
# Set lambda values
lambdas <- c(10, 10, 50, 100)

# Generate data frame with x values and lambda parameters
df <- expand.grid(x = seq(0, 150, 1), lambda = lambdas)

# Calculate corresponding Poisson probabilities
df$probability <- with(df, dpois(x, lambda))

# Plot
custom_colors_2_2 <- c("darkgoldenrod2","cyan3","magenta3", "blue")
levels_custom <- c("10", "10", "50", "100")
colors_custom <- c("10" = custom_colors_2_2[1], 
                   "10" = custom_colors_2_2[2], 
                   "50" = custom_colors_2_2[3], 
                   "100" = custom_colors_2_2[4])

plot <- ggplot(df, aes(x = x, y = probability, color = as.factor(lambda))) +
  geom_line(size = 1) +
  scale_color_manual(values = colors_custom, breaks = levels_custom, name = "MOI") +
  theme_minimal() +
  labs(x = "PFU/Cell", y = "Probability", title = "") +
  theme(legend.key = element_blank(), 
        legend.position = "right",
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 7),
        axis.text = element_text(size = 5),
        axis.title.y = element_text(size = 7, angle = 90),
        axis.title.x = element_text(size = 7),
        text = element_text(family = "Arial"))
ggsave("poisson_variability_model.pdf", width=70, height=55, units="mm",  plot)

#Panel E
# Define MOI values
moi_values <- c(1, 10, 50, 100)

# For PV
pv_start_cov <- c(0.569662921,
                  0.567885117,
                  0.55026455)
pv_slope_cov <- c(0.97826087,
                  0.80952381,
                  1.077777778)

pv <- data.frame(MOI = c(1, 10, 100),
                 virus = "PV",
                 variable = c(rep("start", 3), rep("slope", 3)),
                 CoV = c(pv_start_cov, pv_slope_cov))

# For VACV
vacv_start_early_cov <- c(0.4766415,
                          0.3299413,
                          0.2119404,
                          0.1996086)
# vacv_start_late_cov <- c(0.3253476,
#                          0.2710901,
#                          0.298482,
#                          0.2949246)
vacv_slope_late_cov <- c(0.4006489,
                         0.4991867,
                         0.5725533,
                         0.5497348)

vacv <- data.frame(MOI = rep(moi_values, 3),
                   virus = "VACV",
                   variable = c(rep("start early", 4), rep("start late", 4), rep("slope late", 4)),
                   CoV = c(vacv_start_early_cov, vacv_start_late_cov, vacv_slope_late_cov))

# Bind the dataframes together
CoVs <- rbind(pv, vacv)



# Define the levels of Variable and their corresponding colors

custom_colors_2_2 <- c("red", "yellow3","cyan3", "magenta3","blue")
levels_custom <- c("start","start early","start_late", "slope late", "slope")
colors_custom <- c("start" = custom_colors_2_2[1], 
                   "start early" = custom_colors_2_2[2], 
                   "start_late" = custom_colors_2_2[3],
                   "slope late" = custom_colors_2_2[4],
                   "slope" = custom_colors_2_2[5])
labels_displayed <- c("Start","Start Early", "Start Late", "Slope Late", "Slope")
CoVs$variable <- factor(CoVs$variable, levels = levels_custom)

# Plotting
plot <- ggplot() +
  geom_line(data = CoVs, aes(x = MOI, y = CoV, group=variable, color = variable), size=0.5) +
  geom_point(data=CoVs, aes(x = MOI, y = CoV, fill = variable, color=variable, shape=virus), size=2, stroke=0.1) +
  scale_color_manual(values = colors_custom, guide = "none") +  
  scale_fill_manual(values = colors_custom, breaks = levels_custom, labels = labels_displayed, name = "Parameter")+
  labs(title = "",
       y = "Relative value",
       x = "Relative virus dose") +
  theme_minimal() +
  theme(legend.position = "right",
        legend.key.size = unit(2, "mm"),          
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 5),
        axis.text = element_text(size = 5),
        axis.title.y = element_text(size = 7, angle=90),
        axis.title.x = element_text(size = 7),
        text=element_text(family ="Arial")) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) # Override shape in legend

ggsave("variability_V2.pdf", width=70, height=55, units="mm", plot)


# Supplementary data -----------------------------------------------------------

#Supplementary Figure 1 ####
#Lytic vs non-lytic parameters for MOI data with time trimmed to 12 and 24 HPI.
sicegar_results_df_merged_MOI24HPI <-read.csv("sicegar_results_MOI_24HPI.csv")
plot_data <- read.csv("plot_data.csv")

temp_start<- plot_data %>% 
  filter(HPI<=24) %>% 
  group_by(unique_ID) %>% 
  slice(1) %>% 
  dplyr::select(start_time, start_time_green, unique_ID)

temp_tTest_lytic_nonLytic_24HPI<- sicegar_results_df_merged_24HPI %>%   
  filter(DEC_decision_red!="ambiguous"&DEC_decision_red!="no_signal") %>% 
  left_join(temp_start, by="unique_ID") %>% 
  mutate(COMB_startPoint_x_red2=((start_time*10)/60),
         COMB_startPoint_x2=((start_time_green*10)/60))

var_list_all_24 <- c("COMB_startPoint_x_red2", "COMB_startPoint_x2", "COMB_midPoint1_x", 
                     "COMB_slope1", "COMB_incrementTime", "COMB_maximum_y")


# Define the function that creates and saves the plots
create_save_plot <- function(data, variable, ylabel) {
  p <- ggplot(data = data) + 
    geom_density(aes_string(x=variable, color = "DEC_decision_red", fill="DEC_decision_red"), alpha = 0.6) +
    scale_color_manual(values = c("sigmoidal" ="darkorange1", "double_sigmoidal" = "lightblue")) +
    scale_fill_manual(values = c("sigmoidal" ="darkorange1", "double_sigmoidal" = "lightblue")) +
    theme_minimal() +
    labs(x = ylabel, y = "") + 
    theme(legend.position = "none",
          axis.text = element_text(size = 5),
          axis.title.y = element_blank(),
          axis.title.x=element_text(size = 7),
          text=element_text(family ="Arial"))
  
  # Calculate the max value for the x variable
  max_x <- max(data[[variable]], na.rm = TRUE)
  
  # Check that max_x is finite
  if(is.finite(max_x)) {
    # Set the x-axis limits
    p <- p + scale_x_continuous(limits = c(0, max_x*1.2), labels = label_number(accuracy = 1))
  }
  
  # Save the plot
  ggsave(paste0(variable, "_L_NL_24HPI.pdf"), width=29, height=40, units="mm",plot=p)
}
# Loop through the variable list and generate the plots
for (i in 1:length(variable_list)) {
  create_save_plot(temp_tTest_lytic_nonLytic_24HPI, var_list_all_24[i], label_list[i])
}


#Supplementary Figure 3 ####
QVD_oneStep <- read.csv("QVD_oneStep_growth.csv")
QVD_oneStep_worked <- QVD_oneStep %>% 
  filter(Source=="Total") %>% 
  pivot_longer(cols=4:5, names_to = "Condition", values_to = "Titre") %>% 
  group_by(Condition, HPI) %>% 
  dplyr::summarise(mean=mean(Titre),
                   sd=sd(Titre)/sqrt(n()))

p_combined <- ggplot() +
  geom_path(data=QVD_oneStep_worked, aes(x=HPI, y=mean, group = Condition, color=Condition, linetype="solid"), size=0.5, lineend="round", linejoin="mitre") +
  geom_point(data=QVD_oneStep_worked, aes(x=HPI, y=mean, group = Condition, shape="21"), size=0.125) +
  geom_errorbar(data=QVD_oneStep_worked, aes(x=HPI, ymin=mean-sd, ymax=mean+sd, group=Condition), width=1, color="black", size=0.125) +
  scale_color_manual(values=c("green2","blue2"))+
  scale_x_continuous() +
  labs(x = "Time (HPI)", y = "Titre (PFU/mL)", title = "", colour = "") +
  theme_minimal() +
  theme(legend.position = "right",
        axis.text = element_text(size = 5),
        axis.title = element_text(size = 7),
        text=element_text(family ="Arial"),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7)) +
  scale_y_log10()
ggsave("QVD_oneStep_growth.pdf",width=90, height=60, units="mm", p_combined)

#Supplementary figure 2 ####
sicegar_scatterplot_filtered_green_drugs <- plot_data_with_sicegar_drugs_2_fig6_7 %>% 
  # Filter the data where DEC_decision and DEC_decision_red are not "ambiguous" or "no_signal"
  dplyr:::filter(DEC_decision!="ambiguous"&DEC_decision!="no_signal") %>%
  dplyr:::filter(DEC_decision_red!="ambiguous"&DEC_decision_red!="no_signal") %>% 
  filter(unique_ID!="20221010_2_0_14") %>%  #extreme outlier
  filter(decision_Condition=="Infected_sigmoidal"|decision_Condition=="Control"|decision_Condition=="QVD") %>% 
  dplyr::select(unique_ID, Condition_red, DEC_decision_red, DEC_decision, decision_Condition, COMB_maximum_y, DEC_decision, COMB_startPoint_x, COMB_startPoint_x_red,
                COMB_incrementTime, COMB_midPoint1_x, COMB_slope1) %>% 
  # Rename the selected columns to more easily interpretable names
  dplyr::rename("Maximum A3-GFP (NFI)"=COMB_maximum_y,
                "Productive Start (h)"=COMB_startPoint_x,
                "Infection Time (h)"=COMB_incrementTime,
                "Mid Point (h)"=COMB_midPoint1_x,
                "Slope (NFI/h)"=COMB_slope1,
                "Infection Start (h)"=COMB_startPoint_x_red)

temp_drugs <- sicegar_scatterplot_filtered_green_drugs 
showtext_opts(dpi = 72)
scatter_with_corr_drugs <- function(data, mapping, colour = c("Control" = "green2", "Infected_sigmoidal" ="darkorange1","QVD" = "blue2"), ...) {
  # Split data by Condition
  data_split <- split(data, data$decision_Condition)
  
  # Create an empty ggplot object
  p <- ggplot(data = data, mapping = mapping) + 
    geom_point(aes(color = decision_Condition),alpha = 0.5, size = 0.15, stroke=0.2) +
    scale_color_manual(values = colour) +
    theme_minimal() +
    theme(axis.text = element_blank(),
          axis.title = element_blank(),
          legend.position = "")
  
  # Calculate correlation for each group and add to the plot
  for (i in seq_along(data_split)) {
    x <- eval_data_col(data_split[[i]], mapping$x)
    y <- eval_data_col(data_split[[i]], mapping$y)
    correlation <- round(cor(x, y, use = "pairwise.complete.obs"), 3)
    correlation_3dp <- sprintf("%.3f", correlation)
    # Pad the correlation value with spaces on the left
    correlation_padded <- str_pad(correlation_3dp, width = 6, side = "left")
    
    label <- paste("r = ", correlation_padded)
    p <- p + annotate("text", x = Inf, y = Inf, label = label, 
                      hjust = 1, vjust = 1.2 * i, size = 2, colour = colour[names(data_split)[[i]]], parse = FALSE)
  }
  
  p
}

#Define a function to create density plots
density_with_fake_axes_drugs <- function(data, mapping, colour = c("Control" = "green2", "Infected_sigmoidal" ="darkorange1","QVD" = "blue2"), ...) {
  
  data_split <- split(data, data$decision_Condition)
  
  # Create an empty ggplot object
  p <- ggplot(data = data, mapping = mapping) + 
    geom_density(aes(color = decision_Condition, fill=decision_Condition), alpha = 0.6) +
    scale_color_manual(values = colour) +
    scale_fill_manual(values = colour) +
    theme_minimal() +
    theme(axis.title = element_blank(),
          legend.position = "",
          axis.text = element_text(size=15))
  
  # Calculate the max value for the x variable
  max_x <- max(eval_data_col(data, mapping$x), na.rm = TRUE)
  
  # Set the x-axis limits
  p <- p + scale_x_continuous(limits = c(0, max_x*1.4), labels = label_number(accuracy = 1))
  
  return(p)
}


# Define the variables of interest
variables <- c("COMB_startPoint_x_red", "COMB_startPoint_x", "COMB_midPoint1_x", 
               "COMB_slope1", "COMB_incrementTime", "COMB_maximum_y")

# Create a data frame to store the plot objects
plot_matrix <- data.frame(matrix(vector((length(variables)^2), mode = "list"), ncol = length(variables)))
colnames(plot_matrix) <- variables
rownames(plot_matrix) <- variables

# Generate the upper triangle plots
for(i in 1:(length(variables)-1)) {
  for(j in (i+1):length(variables)) {
    plot_matrix[[variables[i], variables[j]]] <- ggplot() + theme_void()  # Empty plot
  }
}

# Generate the lower triangle plots
for(i in 2:length(variables)) {
  for(j in 1:(i-1)) {
    plot_matrix[[variables[i], variables[j]]] <- scatter_with_corr_drugs(temp_drugs[, c(variables, "decision_Condition")], ggplot2::aes_(x = as.name(variables[i]), y = as.name(variables[j])))
  }
}

# Generate the diagonal plots
for(i in 1:length(variables)) {
  plot_matrix[[variables[i], variables[i]]] <- ggplot() + theme_void()  # Empty plot
}

# Create a list to store the grid of plots
plot_list <- list()
for(i in 1:nrow(plot_matrix)) {
  plot_row <- do.call("arrangeGrob", c(plot_matrix[i, ], ncol = ncol(plot_matrix)))
  plot_list[[i]] <- plot_row
}

labels_txt <- c("Start Early (hr)", "Start Late (hr)",  "Midpoint Late (hr)", "  Slope Late (NFI/hr)", "   Period Late (hr)","  Max Late (NFI)")
top_labels <- lapply(labels_txt, function(label) grid::textGrob(label, gp = grid::gpar(fontsize = 7), vjust=7, hjust=0.4))
right_labels <- lapply(labels_txt, function(label) grid::textGrob(label, gp = grid::gpar(fontsize = 7), rot = -90, vjust=7, hjust=0.6))

# Convert ggplot objects to grobs
grob_list <- unlist(lapply(plot_list, function(x) lapply(x$grobs, function(y) y[[1]])), recursive = FALSE) 

all_grobs <- c(top_labels, right_labels, grob_list)


layout_matrix <- matrix(NA, nrow = length(variables) + 1, ncol = length(variables) + 1)
layout_matrix[1, -ncol(layout_matrix)] <- 1:length(variables)  # top labels
layout_matrix[-1, ncol(layout_matrix)] <- (length(variables) + 1):(2 * length(variables))  # right labels
layout_matrix[-1, -ncol(layout_matrix)] <- (2 * length(variables) + 1):length(all_grobs)  # plots

# Arrange plots with grid.arrange()
p_multipleCondition_matrix <- grid.arrange(grobs = all_grobs, layout_matrix = layout_matrix)

# Save the plot matrix to a tiff file
ggsave("lytic_nonLytic_Drugs_matrix.pdf", width=170, height=170, units="mm", plot = p_multipleCondition_matrix)


#Supplementary figure 5 ####
#read in data
plaque_characterisation <- read.csv("A3_mCh_vacv.csv") %>% 
  dplyr::rename(Virus=?..Virus) 
plaque_characterisation$Virus <- factor(plaque_characterisation$Virus, 
                                        levels = c("WR", "A3GFP", "mCh", "A3GFP_mCh"))
p_plaque_characterisation <- plaque_characterisation %>% 
  filter(Virus!="A3GFP") %>% 
  ggplot(aes(y=Diameter, x=Virus, fill="grey")) +
  geom_signif(y_position = c(0.95, 1, 0.9), xmin = c(1, 1, 2), xmax = c(2, 3, 3),
              annotation = c("****", "****", "ns"), tip_length = 0.025, color="black")+
  geom_boxplot()+
  theme_minimal()+
  #scale_x_discrete(labels = c("VACV WR", "VACV WR \npE/L-mCh", parse("VACV~'\nReporter'^'*'"))) + 
  scale_color_manual(values=c("grey"))+
  scale_fill_manual(values=c("grey"))+
  labs(x = "", y = "Plaque Diameter (mm)")+
  theme(legend.position = "none",
        axis.text = element_text(size = 5),
        axis.title = element_text(size = 7),
        axis.text.x = element_blank(),
        text=element_text(family ="Arial"))
ggsave("VACV_reporter_characterisation.pdf", width=100, height=100, units="mm", p_plaque_characterisation)


#Supplementary figure 6 ####
# Reshape data from wide to long format
data_AraC <- read.csv("data_AraC.csv")
data_AraC_normalised<- data_AraC %>% 
  mutate(normalised_green=green_IntDen-(Background_green*Area),
         normalised_red=red_IntDen-(Background_red*Area))
data_AraC_normalised$Condition <- factor(data_AraC_normalised$Condition, 
                                         levels = c("Control", "AraC"))
data_AraC_normalised_long_8 <- data_AraC_normalised %>%
  pivot_longer(cols = c(normalised_green, normalised_red),
               names_to = "Channel",
               values_to = "Fluorescence") %>% 
  filter(Time=="8")%>% 
  mutate(Fluorescence=abs(Fluorescence))
data_AraC_normalised_long_24 <- data_AraC_normalised %>%
  pivot_longer(cols = c(normalised_green, normalised_red),
               names_to = "Channel",
               values_to = "Fluorescence") %>% 
  filter(Time=="24") %>% 
  mutate(Fluorescence=abs(Fluorescence))
p_AraC <- data_AraC_normalised_long_8 %>% 
  mutate(Condition = fct_rev(Condition)) %>% 
  unite("Channel_condition",Channel, Condition, remove=FALSE) %>% 
  mutate(Channel_condition = fct_rev(Channel_condition)) %>% 
  ggplot(aes(x = Channel, y = Fluorescence, fill = Condition, group=Channel_condition)) +
  geom_boxplot(size=0.2,outlier.size = 0.5)+
  scale_fill_manual(values = c("Control" = "grey50", "AraC" = "blue2"),
                    labels = c("Control" = "Control", "AraC" = "Ara C")) +
  scale_x_discrete(labels = c("normalised_green" = "YFP-A3", "normalised_red" = "pE/L-mCherry")) +
  labs(x = "", y = "Mean NFI", fill = "") +
  theme_minimal() +
  theme(legend.position = "right",
        axis.text = element_text(size = 5),
        axis.title = element_text(size = 7),
        axis.text.x = element_text(size = 5))+
  scale_y_log10(limits = c(1, 4500000))   # Add 'expand' argument

ggsave("AraC_8HPI.pdf", width=70, height=40, units="mm", p_AraC)

#Supplementary figure 8 ####
# Set the input directory
input_directory <- "D:/Data/single cell imaging/viability_testing/cp_output_2"

# Check if the input directory exists
if (!dir.exists(input_directory)) {
  stop("Input directory does not exist.")
}
setwd(input_directory) 
dir.create(file.path("R_output/"))  #creates a directory for outputs
output_directory <- paste0(file.path(input_directory, "R_output/"))
setwd(input_directory)
temp <- list.files(path = input_directory, pattern = "\\.csv$", recursive = TRUE)

#Prompts the user to enter two user inputted strings separated by a comma (with no spaces). These elements should be strings contained in the filenames  
#that are unique to that specific condition. This format has been chosen because it's the native naming convention of wells generated by a JOB in NIS-elements
condition_1 <- "01,05,09"
condition_2 <- "02,06,10"
condition_3 <- "03,07,11"
condition_4 <- "04,08,12"

condition_1_split <- str_split_fixed(condition_1, ",", 3)
condition_2_split <- str_split_fixed(condition_2, ",", 3)
condition_3_split <- str_split_fixed(condition_3, ",", 3)
condition_4_split <- str_split_fixed(condition_4, ",", 3)

condition_1_1 <- paste0(condition_1_split[1],"_")
condition_1_2 <- paste0(condition_1_split[2],"_")
condition_1_3 <- paste0(condition_1_split[3],"_")
condition_2_1 <- paste0(condition_2_split[1],"_")
condition_2_2 <- paste0(condition_2_split[2],"_")
condition_2_3 <- paste0(condition_2_split[3],"_")
condition_3_1 <- paste0(condition_3_split[1],"_")
condition_3_2 <- paste0(condition_3_split[2],"_")
condition_3_3 <- paste0(condition_3_split[3],"_")
condition_4_1 <- paste0(condition_4_split[1],"_")
condition_4_2 <- paste0(condition_4_split[2],"_")
condition_4_3 <- paste0(condition_4_split[3],"_")
#SEPARATE FILES BASED ON CONDITION/TREATMENT
#breaks up the last of input directory files based on the user inputted strings above
filenames_condition_1_nuc <- Filter(function(x) 
  (grepl(paste0("^", condition_1_1), x) & grepl("_nuclei", x)) |
    (grepl(paste0("^", condition_1_2), x) & grepl("_nuclei", x)) |
    (grepl(paste0("^", condition_1_3), x) & grepl("_nuclei", x)), 
  temp)
filenames_condition_2_nuc <- Filter(function(x) 
  (grepl(paste0("^", condition_2_1), x) & grepl("_nuclei", x)) |
    (grepl(paste0("^", condition_2_2), x) & grepl("_nuclei", x)) |
    (grepl(paste0("^", condition_2_3), x) & grepl("_nuclei", x)), 
  temp)
filenames_condition_3_nuc <- Filter(function(x) 
  (grepl(paste0("^", condition_3_1), x) & grepl("_nuclei", x)) |
    (grepl(paste0("^", condition_3_2), x) & grepl("_nuclei", x)) |
    (grepl(paste0("^", condition_3_3), x) & grepl("_nuclei", x)), 
  temp)
filenames_condition_4_nuc <- Filter(function(x) 
  (grepl(paste0("^", condition_4_1), x) & grepl("_nuclei", x)) |
    (grepl(paste0("^", condition_4_2), x) & grepl("_nuclei", x)) |
    (grepl(paste0("^", condition_4_3), x) & grepl("_nuclei", x)), 
  temp)

filenames_condition_1_cell <- Filter(function(x) 
  (grepl(paste0("^", condition_1_1), x) & grepl("_cell", x)) |
    (grepl(paste0("^", condition_1_2), x) & grepl("_cell", x)) |
    (grepl(paste0("^", condition_1_3), x) & grepl("_cell", x)), 
  temp)
filenames_condition_2_cell <- Filter(function(x) 
  (grepl(paste0("^", condition_2_1), x) & grepl("_cell", x)) |
    (grepl(paste0("^", condition_2_2), x) & grepl("_cell", x)) |
    (grepl(paste0("^", condition_2_3), x) & grepl("_cell", x)), 
  temp)
filenames_condition_3_cell <- Filter(function(x) 
  (grepl(paste0("^", condition_3_1), x) & grepl("_cell", x)) |
    (grepl(paste0("^", condition_3_2), x) & grepl("_cell", x)) |
    (grepl(paste0("^", condition_3_3), x) & grepl("_cell", x)), 
  temp)
filenames_condition_4_cell <- Filter(function(x) 
  (grepl(paste0("^", condition_4_1), x) & grepl("_cell", x)) |
    (grepl(paste0("^", condition_4_2), x) & grepl("_cell", x)) |
    (grepl(paste0("^", condition_4_3), x) & grepl("_cell", x)), 
  temp)

#READ DATA INTO TABLES
mydata_assembled_condition_1 <- data.table:::rbindlist(lapply(filenames_condition_1_nuc, read.csv), use.names=FALSE) 
mydata_assembled_condition_2 <- data.table:::rbindlist(lapply(filenames_condition_2_nuc, read.csv), use.names=FALSE)
mydata_assembled_condition_3 <- data.table:::rbindlist(lapply(filenames_condition_3_nuc, read.csv), use.names=FALSE)
mydata_assembled_condition_4 <- data.table:::rbindlist(lapply(filenames_condition_4_nuc, read.csv), use.names=FALSE)

df_list <- list(mydata_assembled_condition_1,mydata_assembled_condition_2,mydata_assembled_condition_3,
                mydata_assembled_condition_4) #put above dfs into a list for use in loop functions
condition_list <- list("DMEM_10","CMC_5","CMC_10","CMC_15")#create a list of condition names for use in loop functions


condition_data_prepper = function(x, var) {
  x %>% 
    dplyr::mutate(TrackObjects_Label = ifelse(is.nan(TrackObjects_Label), NA, TrackObjects_Label)) %>%#change cellprofiler's 'NaN' representation to R recognised 'NA'
    filter(!is.na(TrackObjects_Label)) %>% #remove NAs from the dataset
    unite("FoV_ID", Metadata_Well_ID, Metadata_FoV_ID, remove=FALSE)%>% #constructing unique_ID columns that can be used to group rows for various purposes (e.g. by cell, by FoV, etc)
    unite("unique_ID", FoV_ID, TrackObjects_Label, remove=FALSE)%>%
    unite("FoV_Frame_ID", FoV_ID, Metadata_Frame, remove=FALSE) %>% 
    mutate(Condition = var) #set the Condition name
} 

condition_data_nuc_viability <- map2(.x=df_list, .y=condition_list, .f=condition_data_prepper) #Loops .f over each pairwise combination of .x and .y. Output is a list of dfs


mydata_assembled_condition_1 <- data.table:::rbindlist(lapply(filenames_condition_1_cell, read.csv), use.names=FALSE) 
mydata_assembled_condition_2 <- data.table:::rbindlist(lapply(filenames_condition_2_cell, read.csv), use.names=FALSE)
mydata_assembled_condition_3 <- data.table:::rbindlist(lapply(filenames_condition_3_cell, read.csv), use.names=FALSE)
mydata_assembled_condition_4 <- data.table:::rbindlist(lapply(filenames_condition_4_cell, read.csv), use.names=FALSE)

df_list <- list(mydata_assembled_condition_1,mydata_assembled_condition_2,mydata_assembled_condition_3,
                mydata_assembled_condition_4) #put above dfs into a list for use in loop functions
condition_list <- list("DMEM_10","CMC_5","CMC_10","CMC_15")#create a list of condition names for use in loop functions


condition_data_prepper_cell = function(x, var) {
  x %>% 
    unite("FoV_ID", Metadata_Well_ID, Metadata_FoV_ID, remove=FALSE)%>% #constructing unique_ID columns that can be used to group rows for various purposes (e.g. by cell, by FoV, etc)
    unite("FoV_Frame_ID", FoV_ID, Metadata_Frame, remove=FALSE) %>% 
    mutate(Condition = var) #set the Condition name
} 

condition_data_cell_viability <- map2(.x=df_list, .y=condition_list, .f=condition_data_prepper_cell) 

data_merged_cell_viability <- condition_data_cell_viability %>%
  bind_rows() %>%
  mutate(time_hours = Metadata_Frame) %>%
  group_by(Condition,FoV_ID, time_hours) %>% 
  dplyr::summarise(cell_area=sum(AreaShape_Area),
                   confluency=100*(cell_area/(2048*2048))) %>% 
  filter(!c(FoV_ID=="4_1"&time_hours=="26"|FoV_ID=="4_1"&time_hours=="31"|
              FoV_ID=="4_2"&time_hours=="26"|FoV_ID=="4_2"&time_hours=="31"|
              FoV_ID=="4_4"&time_hours=="26"|FoV_ID=="4_4"&time_hours=="31")) %>% 
  group_by(Condition, time_hours) %>% 
  dplyr::summarise(mean_confluency=mean(confluency))%>% 
  filter(time_hours>9) %>% 
  group_by(Condition) %>% 
  dplyr::mutate(normalised_mean_confluency=mean_confluency-min(mean_confluency))

total_cell_count_df <-  condition_data_cell_viability %>%
  bind_rows() %>%
  mutate(time_hours = Metadata_Frame) %>%
  group_by(Condition,FoV_ID, time_hours) %>% 
  dplyr::summarise(cell_area=sum(AreaShape_Area),
                   cell_count=cell_area/4629.32) %>% #4024.32 empirically determined by measuring area of 50 random cells
  dplyr::select(cell_count, FoV_ID, time_hours)

data_merged_nuc_viability <- condition_data_nuc_viability %>%
  bind_rows() %>%
  mutate(time_hours = Metadata_Frame) %>%
  group_by(Condition,FoV_ID, time_hours) %>%
  dplyr::summarise(dead_cell_count=n_distinct(unique_ID)) %>% 
  ungroup() %>% 
  left_join(total_cell_count_df, by=c("FoV_ID", "time_hours")) %>% 
  group_by(Condition.x, time_hours, FoV_ID) %>% 
  dplyr::summarise(viability=(1-dead_cell_count/cell_count)*100) %>% 
  group_by(Condition.x, time_hours) %>% 
  dplyr::summarise(mean_viability=mean(viability, na.rm=TRUE)) %>% 
  filter(time_hours>9)

scale <- 5
p_combined <- ggplot() +
  geom_path(data=data_merged_cell_viability, aes(x=time_hours, y=normalised_mean_confluency, group = Condition, color=Condition, linetype="solid"),size=0.5,lineend="round",linejoin="mitre") +
  geom_point(data=data_merged_cell_viability, aes(x=time_hours, y=normalised_mean_confluency, group = Condition, shape="21"),size=0.125) +
  geom_path(data=data_merged_nuc_viability, aes(x=time_hours, y=mean_viability/scale, group = Condition.x, colour=Condition.x, linetype="dashed"),size=0.5,lineend="round",linejoin="mitre") +
  geom_point(data=data_merged_nuc_viability, aes(x=time_hours, y=mean_viability/scale, group = Condition.x, shape="21"),size=0.125) +
  scale_x_continuous()+
  labs(x = "Time (HPI)", y = "Confluency (%)",
       title = "", colour = "")+
  theme_minimal()+
  theme(legend.position = "right",
        axis.text = element_text(size = 5),
        axis.title = element_text(size = 7),
        text=element_text(family ="Arial"),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7)) + 
  scale_y_continuous(sec.axis = sec_axis(~.*scale, name="Viability (%)"))
ggsave("viability_growth.pdf",width=95, height=60, units="mm", p_combined)