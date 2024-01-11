##-------------- This code is related to epigenetic clocks 
## The package used is 'methylclock' 


# TCGA dataset
# August, 2023
# Analysis of epigenetic clocks

# ---------------------------------
# Data Loading and Preprocessing
# --------------------------------

#Set working directory
getwd()
setwd("~/data")


##Call libraries
library(RColorBrewer)
library(methylclock)
library(ggplot2)
library(dplyr)
library(gtsummary)
library(ChAMP)


#-------------------------------------------------------------------------------------------
#                 P R E P A R E   T H E    D A T A S E T
#--------------------------------------------------------------------------------------------

##This is not needed if you are already working with your dataset and performed all these steps!
# To perform this analysis, you will need a phenotype data and a methylation data (beta-values)

### Load TCGA methylation data already processed by ChAMP package
load("~/ChAMP_input_TCGA.RData")

# Set another directory to save files
setwd("~/results_methylclock/")


#Let's start!

##Discover missing probes in each clock
missing <-checkClocks(met$beta)
#       clock     Cpgs_in_clock  missing_CpG     percentage
#1     Horvath           353           33        9.3
#2      Hannum            71           12       16.9
#3      Levine           513           50        9.7
#4 SkinHorvath           391           75       19.2
#5       PedBE            94           15       16.0
#6          Wu           111            6        5.4
#7          TL           140           49       35.0
#8        BLUP        319607        52871       16.5
#9          EN           514           77       15.0



#-------------------------------------------------------------------------------------------
#                T U M O R     S A M P L E S
#--------------------------------------------------------------------------------------------
#Filtering to keep only tumor samples
tumor_samples_tcga <- met$pd[met$pd$shortLetterCode == "TP", ]

# Filtering the object 'met$beta' to keep only tumor samples
beta_tumor_tcga <- met$beta[, colnames(met$beta) %in% rownames(tumor_samples_tcga)]


# Perform DNAmAge analysis - 'Levine' clock 
Age_tcga <- DNAmAge(beta_tumor_tcga, 
                     #clocks = "Levine", 
                     toBetas = FALSE, 
                     cell.count = FALSE, 
                     age = tumor_samples_tcga$age_at_index)




summary(Age_tcga$age)
##Age_tcga$Levine
###Min.  1st Qu.  Median    Mean    3rd Qu.    Max. 
##40.85   77.93   88.70    91.59    105.35    150.61

##Age_tcga$age
##Min.    1st Qu.  Median    Mean    3rd Qu.    Max. 
##38.00   56.00    62.00    61.74     67.00     83.00


#Remove other epigenetic clocks that I won't use (for now I am only interested in Levine's clock!)
#rm(coefBLUP, coefEN, coefHorvath, coefHannum, coefSkin, coefTL, coefWu, coefPedBE)


#------------
# OPTIONAL!!
#------------
##If I want all clocks
Age_tcga_all_clocks <- DNAmAge(beta_tumor_tcga,
                              age = tumor_samples_tcga$age_at_index,
                              toBetas = FALSE,
                              cell.count = FALSE)



#####################################################################################################################
#--------------------------------------------------------------------------
#                             G R A P H I C    V I S U A L I Z A T I O N
#                                   T U M O R A L   S A M P L E S 
#--------------------------------------------------------------------------

# Create a df with chronological age and tissue epigenetic age
data_tcga <- data.frame(Real_Age = Age_tcga$age, 
                        Methylation_Age = Age_tcga$EN)


# Common boxplot (not for presentations) 
boxplot(data_tcga, 
        names = c("Real_Age", "Methylation_Age"),
        ylab = "Years", 
        col = c("blue", "red"))


# Statistical test -- paired samples!
data_tcga_test <- wilcox.test(data_tcga$Real_Age, 
                                data_tcga$Methylation_Age, 
                                exact = T,
                                paired = T,
                                correct = T)

# Statistics 
data_tcga_test$p.value
##Exact and correct p-value:  4.75343e-27




# Create a new df for better plotting
data_tcga_plotdf <- data.frame(Group = c(rep("Chronological Age",nrow(data_tcga)), 
                                         rep("DNAm Age",nrow(data_tcga))),
                                 Age = c(data_tcga$Real_Age, 
                                         data_tcga$Methylation_Age))


# Visualize plot
ggplot(data_tcga_plotdf, aes(x = Group, y = Age, fill = Group, color = Group)) +
  geom_boxplot(alpha = 0.2, outlier.shape = NA) +
  geom_jitter(width = 0.4, height = 0, size = 2, alpha = 1, shape = 25) +
  labs(x = "", y = "Years") +
  ggtitle("Chronological age versus tumor DNAm age - EN") +
  geom_text(aes(label = data_tcga_test$p.value),
            x = 1.5, y = max(data_tcga_plotdf$Age), vjust = -1.5,
            size = 3, color = "black", fontface = "plain") +
  scale_fill_viridis(discrete = TRUE, option = "turbo") +
  scale_color_viridis(discrete = TRUE, option = "turbo") +
  theme(
    panel.grid.major = element_line(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 13, face = "plain"),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 10),
    text = element_text(size = 12, color = "black", face = "plain")
  )


## Correlation plot
plotDNAmAge(Age_tcga$EN, 
            Age_tcga$age, 
            tit="Bayesian Neural Network to predict DNAm age based on EN’s CpGs")

#rm(coefBLUP, coefEN, coefHannum, coefHorvath, coefPedBE, coefTL, coefWu, coefSkin)




####################################################################################################################
#----------------------------------
#  N O R M A L     S A M P L E S
#----------------------------------

# Filtering to keep only normal samples
normal_samples_tcga <- met$pd[met$pd$shortLetterCode == "NT", ]

# Filtering object 'met$beta' to keep only normal samples
beta_normal_tcga <- met$beta[, colnames(met$beta) %in% rownames(normal_samples_tcga)]

# Perform DNAmAge analysis - Levine's clock - in normal samples
Age_normal_tcga <- DNAmAge(beta_normal_tcga, 
                      clocks = "Levine", 
                      toBetas = FALSE, 
                      cell.count = FALSE, 
                      age = normal_samples_tcga$age_at_index)

rm(coefBLUP, coefEN, coefHannum, coefHorvath, coefPedBE, coefTL, coefWu, coefSkin)


# !!!!!!!!!!!!!!!!!  THERE IS OUTLIERS IN THIS DF!

# Calculate quantiles from Levine's column
Q1 <- quantile(Age_normal_tcga$Levine, 0.25)
Q3 <- quantile(Age_normal_tcga$Levine, 0.75)

# Calculate interquartil interval (IQR)
IQR <- Q3 - Q1

# Define inferior and superior limits to identify outliers
limite_inferior <- Q1 - 1.5 * IQR
limite_superior <- Q3 + 1.5 * IQR

# Identify which samples are outliers
outliers <- Age_normal_tcga[Age_normal_tcga$Levine < limite_inferior | Age_normal_tcga$Levine > limite_superior, ]

# Print it
print(outliers)

##Samples considered outliers
# A tibble: 2 × 6
#id                              Levine   ageAcc.Levine    ageAcc2.Levine   age type  
#<chr>                           <dbl>         <dbl>          <dbl>       <int> <chr> 
#1 TCGA-CV-5978-11A-01D-1684-05   14.4         -38.6          -49.8         53 Normal
#2 TCGA-CV-7089-11A-01D-2014-05   20.0         -54.0          -58.9         74 Normal

# Remove outliers
linhas_para_remover <- c(10, 13)  # Check samples rows (lines)

# Remove it from the df
Age_normal_tcga_f <- Age_normal_tcga[-linhas_para_remover, ]


rm(outliers, IQR, limite_inferior, limite_superior, linhas_para_remover, Q1, Q3)



#--------------------------------------------------------------------------
#            G R A P H I C    V I S U A L I Z A T I O N
#                 N O R M A L   S A M P L E S 
#---------------------------------------------------------------------------

# Create a df with chronological age and tissue epigenetic age
data_normal_f <- data.frame(Real_Age = Age_normal_tcga_f$age, 
                          Methylation_Age = Age_normal_tcga_f$Levine)


# Statistical test (paired!)
test_result_normal_tcga_f <- wilcox.test(data_normal_f$Real_Age,
                                       data_normal_f$Methylation_Age, 
                                       paired = T,
                                       exact = T,
                                       correct = T)



test_result_normal_tcga_f
#p-value = 0.00024414


# Correlation plot
plotDNAmAge(Age_normal_tcga_f$Levine, 
            Age_normal_tcga_f$age, 
            tit="Bayesian Neural Network to predict DNAm age based on Levine’s CpGs")



# Create a new df for better plotting
result_data_normal_tcga_f <- data.frame(Group = c(rep("Chronological Age", nrow(data_normal_f)), 
                                                rep("DNAm Age", nrow(data_normal_f))),
                                      Age = c(data_normal_f$Real_Age, 
                                              data_normal_f$Methylation_Age))



### Visualize plotting
ggplot(result_data_normal_tcga_f, aes(x = Group, y = Age, fill = Group, color = Group)) +
  geom_boxplot(alpha = 0.2, outlier.shape = NA) +
  geom_jitter(width = 0.4, height = 0, size = 2, alpha = 1, shape = 25) +
  labs(x = "", y = "Years") +
  ggtitle("Chronological age versus normal tissue epigenetic age") +
  geom_text(aes(label = test_result_normal_tcga_f$p.value),
            x = 1.5, y = max(result_data_normal_tcga_f$Age), vjust = -1.5,
            size = 3, color = "black", fontface = "plain") +
  scale_fill_viridis(discrete = TRUE, option = "turbo") +
  scale_color_viridis(discrete = TRUE, option = "turbo") +
  theme(
    panel.grid.major = element_line(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 13, face = "plain"),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 10),
    text = element_text(size = 12, color = "black", face = "plain")
  )

#print
result_data_normal_tcga_f




#-----------------------------------------------------------------------------------------
#                                     G R A P H I C    V I S U A L I Z A T I O N

#                                  T U M O R A L  VS  N O R M A L      S A M P L E S 
#-----------------------------------------------------------------------------------------

# Create a new df binding chronological and epigenetic age for both tumor and normal samples! 
# These is for ALL samples!

# Add a column named 'type' in tumor df 
Age_tcga <- Age_tcga %>%
  mutate(type = "Tumor")

# Add a column named 'type' in normal dfs
Age_normal_tcga <- Age_normal_tcga %>%
  mutate(type = "Normal")

# Combine both df with 'rbind'
combined_data_tcga <- rbind(Age_tcga, Age_normal_tcga)
summary(combined_data_tcga)


# Statistical test between epigenetic age of all samples according to their type
diff_ep.age_tcga <- wilcox.test(combined_data_tcga$Levine ~combined_data_tcga$type)
diff_ep.age_tcga
#p-value = 0.00779


# Visualize this result 
ggplot(combined_data_tcga, aes(x = type, y = Levine, fill = type, color = type)) +
  geom_boxplot(alpha = 0.2, outlier.shape = NA) +
  geom_jitter(width = 0.4, height = 0, size = 2, alpha = 1, shape = 25) +
  labs(x = "", y = "Age in Years") +
  ggtitle("DNAm age in tumor and normal tissue - TCGA") +
  geom_text(aes(label = diff_ep.age_tcga$p.value),
            x = 1.5, y = max(combined_data_tcga$Levine), vjust = -1.5,
            size = 3, color = "black", fontface = "plain") +
  scale_fill_viridis(discrete = TRUE, option = "turbo") +
  scale_color_viridis(discrete = TRUE, option = "turbo") +
  theme(
    panel.grid.major = element_line(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 13, face = "plain"),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 10),
    text = element_text(size = 12, color = "black", face = "plain")
  )





#-----------------------------------------------------------------------------------------
#                                     G R A P H I C    V I S U A L I Z A T I O N

#                                           P A I R E D   S A M P L E S 
#-----------------------------------------------------------------------------------------


##First, we need to keep only the first 12 digits of 'id' column in both df; it is the sample 'barcode'
Age_tcga$id <- substring(Age_tcga$id, 1, 12) 
Age_normal_tcga_f$id <- substring(Age_normal_tcga_f$id, 1, 12)


# Then, we intersect normal and tumor samples with the same barcode number
paired_samples_tcga_f <- intersect(Age_tcga$id, Age_normal_tcga_f$id)


# We keep the epigenetic age for these paired samples 
methylation_age_normal_f <- Age_normal_tcga_f$Levine[Age_normal_tcga_f$id %in% paired_samples_tcga_f]
methylation_age_tumor_f <- Age_tcga$Levine[Age_tcga$id %in% paired_samples_tcga_f]


# Create a new df with both chronological and epigenetic age of tumor and normal samples
combined_data_tcga2_f <- data.frame(
  id = c(Age_normal_tcga_f$id[Age_normal_tcga_f$id %in% paired_samples_tcga_f],
         Age_tcga$id[Age_tcga$id %in% paired_samples_tcga_f]),
  Chronological_Age = c(Age_normal_tcga_f$age[Age_normal_tcga_f$id %in% paired_samples_tcga_f],
                        Age_tcga$age[Age_tcga$id %in% paired_samples_tcga_f]),
  Methylation_Age = c(methylation_age_normal_f, methylation_age_tumor_f),
  type = rep(c("Normal", "Tumor"), each = length(paired_samples_tcga_f))
)


# Perform statistical test
paired_tcga_methyl.age_f <- wilcox.test(methylation_age_normal_f, methylation_age_tumor_f, paired = TRUE)
paired_tcga_methyl.age_f

#Check medians
median(methylation_age_normal_f)
#[1] 80.01018
median(methylation_age_tumor_f)
#[1] 91.7769

# Define color palette for plotting
cores.paired <- c("#998ec3", "#f1a340")


#Plot - this plot connects paired samples!
ggpaired(combined_data_tcga2_f, 
         x = 'type', 
         y = 'Methylation_Age', 
         id = "id",
         color = 'type',
         palette = cores.paired,
         line.color = "gray", 
         line.size = 0.4) +
  stat_compare_means(label = "p.format", 
                     method = 'wilcox.test',
                     paired = T,
                     label.y = max(combined_data_tcga2_f$Methylation_Age) + 2) +
  labs(title = "DNAm in paired samples",
       x = "",
       y = "DNAm Age") +
  theme_pubr()



median(Age_normal_tcga_f$Levine)
median(Age_normal_tcga_f$age)




#-----------------------------------------------------------------------------------------
#                                     S U R V I V A L     A N A L Y S I S  
#                                         T U M O R    S A M P L E S 
#-----------------------------------------------------------------------------------------

#   For this analysis, we use the ageAcc2 column.
#   ageAcc2: this measure is the residuals obtained after regressing chronological age and DNAmAge (similar to IEAA). 



#  ATTENTION: always check the column name in "Age_tcga", due to other all other data manipulation, sometimes this column name can be changed. 
#  We want to combine these dataframes according to the 'patient' column in clinical.tcga

# 1) We have to rename this column to combine dataframes - column 'id' is the same as the column 'patient' - 12 digits.
Age_tcga_surv <- Age_tcga %>%
  rename("patient" = "id")
colnames(Age_tcga_surv)[1] = "barcode"

# 2) Combine dataframes
surv_tumor_tcga <- merge(pd_filtered, Age_tcga_surv, by = "barcode", all.x = TRUE)

# 3) Remove normal samples from this df. We want to perform survival analysis based on tumor samples profile. 
surv_tumor_tcga <- surv_tumor_tcga[surv_tumor_tcga$shortLetterCode != "NT", ]


# 4) Create a input df for the survival analysis. It contains selected clinicopathological data. 
surv_tumor_tcga <- surv_tumor_tcga %>%
  select(patient, 
         barcode, 
         age_at_index, 
         Levine, 
         ageAcc2.Levine, 
         ageAcc.Levine,
         ajcc_clinical_stage, 
         ajcc_pathologic_n,
         ajcc_clinical_n, 
         alcohol_history,
         cigarettes_per_day,
         ethnicity, 
         gender, 
         race, 
         gender, 
         days_to_last_follow_up, 
         days_to_death, 
         vital_status)


# 5) Create a column "days_to_event".
# This columns contains data regarding days to death or days to last followup (patient still alive or censored)
surv_tumor_tcga$days_to_event <- ifelse(is.na(surv_tumor_tcga$days_to_death),
                                        surv_tumor_tcga$days_to_last_follow_up,
                                        surv_tumor_tcga$days_to_death)



#-----------------------------------------------------------------------------------------
#      D E F I N E    C A T E G O R I C A L    V A R I A B L E S !
#-----------------------------------------------------------------------------------------

# Epigenetic Age Acceleration: we defined groups "ACC" and "DEC" based on their 'ageAcc2' measures. 
# If positive = 'ACC'; If negative = 'DEC'
surv_tumor_tcga$AgeAcceleration <- ifelse(surv_tumor_tcga$ageAcc2.EN > 0, "ACC", "DEC")


# Vital Status categorization: dead (1) or alive (0)
surv_tumor_tcga$vital_status <- ifelse(surv_tumor_tcga$vital_status == "Dead", 1, 0)


# Nodal Invasion: Positive or Negative; based on pathological evaluation!
surv_tumor_tcga$nodal_invasion_pathologic <- ifelse(
  surv_tumor_tcga$ajcc_pathologic_n == "N0", "Negative",
  ifelse(surv_tumor_tcga$ajcc_pathologic_n == "NX", "NA", "Positive"))


# Clinical Disease Stage: early and advanced
surv_tumor_tcga$stage <- ifelse(
  surv_tumor_tcga$ajcc_clinical_stage %in% c("Stage I", "Stage II"), "early",
  ifelse(is.na(surv_tumor_tcga$ajcc_clinical_stage), "NA", "advanced"))


# Skin Color: white or black
surv_tumor_tcga$skincolor <- ifelse(surv_tumor_tcga$race == "white", "white",
                             ifelse(surv_tumor_tcga$race == "not reported", "NA", "black"))

# Alcohol consumation: yes or no
surv_tumor_tcga$alcohol <- ifelse(surv_tumor_tcga$alcohol_history == "Yes", "Yes",
                                    ifelse(surv_tumor_tcga$alcohol_history == "Not Reported", "NA", "No"))

# Tobacco use: smoker and no-smoker
surv_tumor_tcga$smokingstatus <- ifelse(is.na(surv_tumor_tcga$cigarettes_per_day), "no-smoker", "smoker")




#-----------------------------------------------------------------------------------------
#      P E R F O R M     T H E     A N A L Y S I S
#-----------------------------------------------------------------------------------------

# 1) Create a survival object
surv_obj <- Surv(surv_tumor_tcga$days_to_event, surv_tumor_tcga$vital_status)


# 2) Adjust survival model considering the variable I want to observe (in this case, Age Acceleration)
surv_model_acceleration <- survfit(surv_obj ~ AgeAcceleration, data = surv_tumor_tcga)

# 3) Statistical test
surv_test <- survdiff(surv_obj ~ AgeAcceleration, data = surv_tumor_tcga)

# 4) Print
surv_test

##                   N     Observed    Expected    (O-E)^2/E    (O-E)^2/V
#AgeAcceleration=ACC 48       15         21.7         2.07         3.91
#AgeAcceleration=DEC 66       33         26.3         1.71         3.91
#Chisq= 3.9  on 1 degrees of freedom, p= 0.05


# 5) Create a sfit object!
# This computes an estimate of a survival curve for censored data using Kaplan-Meier method
sfit <- survfit(Surv(surv_tumor_tcga$days_to_event, surv_tumor_tcga$vital_status)~ AgeAcceleration, 
                data=surv_tumor_tcga)
sfit
#                     n   events   median   0.95LCL  0.95UCL
#AgeAcceleration=ACC 48     15     2319      1972      NA
#AgeAcceleration=DEC 66     33      988       663      NA


#----------------------------------------------------------------------------------------------
#      G R A P H I C A L     V I S U A L I Z A T I O N: 'ACC' AND 'DEC' GROUPS OVERALL SURVIVAL
#----------------------------------------------------------------------------------------------

## Plot the survival analysis (Kaplan-Meier method)
ggsurvplot(sfit,
           conf.int=TRUE, 
           pval=TRUE, 
           risk.table=TRUE, 
           legend.labs=c("ACC", "DEC"), 
           legend.title="Group %", 
           xscale = "d_y",
           break.x.by=365,##diaa
           xlim = c(0,1825),##plotar apenas 5 anos
           palette=c("dodgerblue2", "orchid2"), 
           surv.median.line = "hv",
           title="Survival Curve - Age Acceleration - EN", 
           risk.table.height=.23)

summary(sfit, times=seq(0,365*5,365))


#save the dataframe used as input!
write.csv(surv_tumor_tcga, file = "overall_survival_tcga_according_to_ageacc.csv")




#-----------------------------------------------------------------------------------------------------------------------
# p l o t t i n g     o v e r a l l   s u r v i v a l    c u r v e s    a c c o r d i n g    t o    p a r a m e t e r s.
#------------------------------------------------------------------------------------------------------------------------


## for this plot, we need to subset our samples groups 'ACC' and 'DEC' according to the desired parameter
# that we want to observe the effect. 
#Example: we want to see the survival curves of disease stage classification according
# to their age acceleration classification (ACC or DEC). Samples with advanced disease and DEC groups perform worse than
# samples with ACC-advanced disease group?



#For this, we need to subset samples according to the parameter we want to observe and their ageacc categorization

# Filtering samples in 'DEC' group for 'x' parameter (insert any wanted parameter)
dec_early <- subset(surv_tumor_tcga, AgeAcceleration == "DEC" & stage == "early")

dec_adv  <- subset(surv_tumor_tcga, AgeAcceleration == "DEC" & stage == "advanced")


# Filtering samples in 'ACC' group for x' parameter (insert any wanted parameter):
acc_early  <- subset(surv_tumor_tcga, AgeAcceleration == "ACC" & stage == "early")

acc_adv <- subset(surv_tumor_tcga, AgeAcceleration == "ACC" & stage == "advanced")



# Create a df combining these subsets 
combined_data_stage <- rbind(
  transform(dec_early, Group = "DEC - Early Stage"),
  transform(dec_adv, Group = "DEC - Advanced Stage"),
  transform(acc_early, Group = "ACC - Early Stage"),
  transform(acc_adv, Group = "ACC - Advanced Stage")
)

# Adjust the survival model with this dataset
surv_model_stage <- survfit(Surv(days_to_event, vital_status) ~ Group, data = combined_data_stage)


# Plotting the survival curve 
ggsurvplot(surv_model_stage, 
           data = combined_data_stage,
           xlab = "Time (years)", 
           ylab = "Overall Survival",
           break.x.by=365, ##dias
           xscale = "d_y", ##dias para anos
           xlim = c(0,1825),##plotar apenas 5 anos
           main = "Survival Curve by Disease Stage",
           palette=c("dodgerblue2", "orchid2", "orangered2","thistle"),
           risk.table = T, 
           pval = T, 
           pval.method = T,
           conf.int = F,
           surv.median.line = "hv",
           risk.table.height=.34, 
           legend.title = "AJCC clinical stage classification")


summary(survfit(Surv(days_to_event, vital_status) ~ Group, data = combined_data_stage))


## Linear regression
coxph(Surv(days_to_event, vital_status) ~ Group, data = combined_data_stage)%>% 
  tbl_regression(exp = TRUE) 




#----------------
# Smoking Status
#----------------

#Filtering samples in 'DEC' group for 'x' parameter (insert any wanted parameter)

dec_yes <- subset(surv_tumor_tcga, AgeAcceleration == "DEC" & smokingstatus == "smoker")

dec_no  <- subset(surv_tumor_tcga, AgeAcceleration == "DEC" & smokingstatus == "no-smoker")


# # Filtering samples in 'ACC' group for x' parameter (insert any wanted parameter):
acc_yes  <- subset(surv_tumor_tcga, AgeAcceleration == "ACC" & smokingstatus == "smoker")

acc_no <- subset(surv_tumor_tcga, AgeAcceleration == "ACC" & smokingstatus == "no-smoker")



# Create a df combining these subsets 
combined_data_smoking <- rbind(
  transform(dec_yes, Group = "DEC - Smoker"),
  transform(dec_no, Group = "DEC - Non-smoker"),
  transform(acc_yes, Group = "ACC - Smoker"),
  transform(acc_no, Group = "ACC - Non-smoker")
)

# Adjust the survival model with this dataset
surv_model_smoking <- survfit(Surv(days_to_event, vital_status) ~ Group, data = combined_data_smoking)


# Plotting the survival curve 
ggsurvplot(surv_model_smoking, 
           data = combined_data_smoking,
           xlab = "Time (years)", 
           ylab = "Overall Survival",
           break.x.by=365, ##dias
           xscale = "d_y", ##dias para anos
           xlim = c(0,1825),##plotar apenas 5 anos
           main = "Survival Curve by Smoking Status",
           palette=c("dodgerblue2", "orchid2", "orangered2","thistle"),
           risk.table = T, 
           pval = T, 
           pval.method = T,
           conf.int = F,
           surv.median.line = "hv",
           risk.table.height=.34, 
           legend.title = "Smoking Status")


summary(survfit(Surv(days_to_event, vital_status) ~ Group, data = combined_data_smoking))


coxph(Surv(days_to_event, vital_status) ~ Group, data = combined_data_smoking)%>% 
  tbl_regression(exp = TRUE)



#----------------
# Nodal Status
#----------------

#Filtering samples in 'DEC' group for 'x' parameter (insert any wanted parameter)

dec_pos <- subset(surv_tumor_tcga, AgeAcceleration == "DEC" & nodal_invasion_pathologic == "Positive")

dec_neg  <- subset(surv_tumor_tcga, AgeAcceleration == "DEC" & nodal_invasion_pathologic == "Negative")


# # Filtering samples in 'ACC' group for x' parameter (insert any wanted parameter):
acc_pos  <- subset(surv_tumor_tcga, AgeAcceleration == "ACC" & nodal_invasion_pathologic == "Positive")

acc_neg <- subset(surv_tumor_tcga, AgeAcceleration == "ACC" & nodal_invasion_pathologic == "Negative")



# Create a df combining these subsets 
combined_data_nodal <- rbind(
  transform(dec_pos, Group = "DEC - Positive"),
  transform(dec_neg, Group = "DEC - Negative"),
  transform(acc_pos, Group = "ACC - Positive"),
  transform(acc_neg, Group = "ACC - Negative")
)

# Adjust the survival model with this dataset
surv_model_nodal <- survfit(Surv(days_to_event, vital_status) ~ Group, data = combined_data_nodal)


# Plotting the survival curve 
ggsurvplot(surv_model_nodal, 
           data = combined_data_nodal,
           xlab = "Time (years)", 
           ylab = "Overall Survival",
           break.x.by=365, ##dias
           xscale = "d_y", ##dias para anos
           xlim = c(0,1825),##plotar apenas 5 anos
           main = "Survival Curve by Nodal Invasion",
           palette=c("dodgerblue2", "orchid2", "orangered2","thistle"),
           risk.table = T, 
           pval = T, 
           pval.method = T,
           conf.int = F,
           surv.median.line = "hv",
           risk.table.height=.34, 
           legend.title = "")


summary(survfit(Surv(days_to_event, vital_status) ~ Group, data = combined_data_nodal))


library('gtsummary')
coxph(Surv(days_to_event, vital_status) ~ Group, data = combined_data_nodal)%>% 
  tbl_regression(exp = TRUE)





#---------------------
# Alcohol consumption
#---------------------

#Filtering samples in 'DEC' group for 'x' parameter (insert any wanted parameter)

dec_drink <- subset(surv_tumor_tcga, AgeAcceleration == "DEC" & alcohol == "Yes")

dec_nodrink  <- subset(surv_tumor_tcga, AgeAcceleration == "DEC" & alcohol == "No")


# # Filtering samples in 'ACC' group for x' parameter (insert any wanted parameter):
acc_drink  <- subset(surv_tumor_tcga, AgeAcceleration == "ACC" & alcohol == "Yes")

acc_nodrink <- subset(surv_tumor_tcga, AgeAcceleration == "ACC" & alcohol == "No")



# Create a df combining these subsets 
combined_data_alcohol <- rbind(
  transform(dec_drink, Group = "DEC - Yes"),
  transform(dec_nodrink, Group = "DEC - No"),
  transform(acc_drink, Group = "ACC - Yes"),
  transform(acc_nodrink, Group = "ACC - No")
)

# Adjust the survival model with this dataset
surv_model_alcohol <- survfit(Surv(days_to_event, vital_status) ~ Group, data = combined_data_alcohol)


# Plotting the survival curve 
ggsurvplot(surv_model_alcohol, 
           data = combined_data_alcohol,
           xlab = "Time (years)", 
           ylab = "Overall Survival",
           break.x.by=365, ##dias
           xscale = "d_y", ##dias para anos
           xlim = c(0,1825),##plotar apenas 5 anos
           main = "Survival Curve by Alcohol Use",
           palette=c("dodgerblue2", "orchid2", "orangered2","thistle"),
           risk.table = T, 
           pval = T, 
           pval.method = T,
           conf.int = F,
           surv.median.line = "hv",
           risk.table.height=.34, 
           legend.title = "")


summary(survfit(Surv(days_to_event, vital_status) ~ Group, data = combined_data_alcohol))


coxph(Surv(days_to_event, vital_status) ~ Group, data = combined_data_alcohol)%>% 
  tbl_regression(exp = TRUE)





# --------------------------------------------------------------------------------------------

### Chi-square test
str(surv_tumor_tcga)


# Create a contingency table for the varible alcohol
data1 <-surv_tumor_tcga[surv_tumor_tcga$alcohol != "NA", ] ##Exclude NA values
table_alcohol <- table(data1$AgeAcceleration, data1$alcohol)

# Create a contingency table for the varible "tobacco"
table_tobacco <- table(surv_tumor_tcga$AgeAcceleration, surv_tumor_tcga$smokingstatus)


# Create a contingency table for the varible "nodal_invasion"
data2 <-surv_tumor_tcga[surv_tumor_tcga$nodal_invasion_pathologic != "NA", ]  ##Exclude NA values
data2 <- data2[!is.na(data2$nodal_invasion_pathologic), ]
table_nodal_invasion <- table(data2$AgeAcceleration, data2$nodal_invasion_pathologic)
table_nodal_invasion

# Create a contingency table for the varible "stage"
data3 <-surv_tumor_tcga[surv_tumor_tcga$stage != "NA", ] ##Exclude NA values
table_stage <- table(data3$AgeAcceleration, data3$stage)


# Visualize contingency tables
print(table_alcohol)
print(table_tobacco)
print(table_nodal_invasion)
print(table_stage)


plot(table_nodal_invasion)

#Chi-square test
alcohol <- chisq.test(table_alcohol)
alcohol
#X-squared = 0.30354, df = 1, p-value = 0.581

#Chi-square test
tobacco <- chisq.test(table_tobacco)
tobacco
#X-squared = 0.019492, df = 1, p-value = 0.889

#Chi-square test
nodal <- chisq.test(table_nodal_invasion)
nodal
#X-squared = 3.3223, df = 1, p-value = 0.06835

#Chi-square test
stage <- chisq.test(table_stage)
stage
#X-squared = 0.042285, df = 1, p-value = 0.8371



#####################################################################################################################

#                                         S U R V I V A L     A N A L Y S I S  
#                                           N O R M A L     S A M P L E S  
#------------------------------------------------------------------------------------------------------------------
#We repeat the same steps we did in tumor samples, bot with a focus on normal samples!

#Filtering according to patient number
Age_normaltcga_surv <- Age_normal_tcga %>%
  rename(patient = id)


#Remove tumor samples from clinical data
clinical.tcga.normals <- clinical.tcga[clinical.tcga$shortLetterCode != "TP", ]


# Combining dataframes
surv_normal_tcga <- merge(clinical.tcga.normals, Age_normaltcga_surv, by = "patient", all.x = TRUE)

# Create a input df for the survival analysis. It contains selected clinicopathological data. 
surv_normal_tcga <- surv_normal_tcga %>%
  select(barcode, age_at_index, Levine, ageAcc2.Levine, ageAcc.Levine,
         ajcc_clinical_stage, ajcc_pathologic_n,ajcc_clinical_n, 
         alcohol_history,cigarettes_per_day,
         ethnicity, gender, race, gender, 
         days_to_last_follow_up, days_to_death, vital_status)


##Create a column "days_to_event".
# This columns contains data regarding days to death or days to last followup (patient still alive or censored)
surv_normal_tcga$days_to_event <- ifelse(is.na(surv_normal_tcga$days_to_death),
                                         surv_normal_tcga$days_to_last_follow_up,
                                         surv_normal_tcga$days_to_death)

#Define categorical variables

#Epigenetic Age Acceleration groups
surv_normal_tcga$AgeAcceleration <- ifelse(surv_normal_tcga$ageAcc2.Levine > 0, "ACC", "DEC")

#Vital Status
surv_normal_tcga$vital_status <- ifelse(surv_normal_tcga$vital_status == "Dead", 1, 0)


###########################################################################################################################

## Combining 'normal' and 'tumoral' survival datasets into a single one

data <- rbind(
  data.frame(grupo = "Normal", tipo = "ACC", 
             taxa = surv_normal_tcga$ageAcc2.Levine[surv_normal_tcga$AgeAcceleration == "ACC"]),
  data.frame(grupo = "Tumoral", tipo = "ACC", 
             taxa = surv_tumor_tcga$ageAcc2.Levine[surv_tumor_tcga$AgeAcceleration == "ACC"]),
  data.frame(grupo = "Normal", tipo = "DEC", 
             taxa = surv_normal_tcga$ageAcc2.Levine[surv_normal_tcga$AgeAcceleration == "DEC"]),
  data.frame(grupo = "Tumoral", tipo = "DEC",
             taxa = surv_tumor_tcga$ageAcc2.Levine[surv_tumor_tcga$AgeAcceleration == "DEC"])
)


# Define color palette 
cores <- c("Normal" = "#2962FF", 
           "Tumoral" = "#FF6F00")

# Define jittering color palette
cores_jitter <- c("Normal" = "#2962FF", 
                  "Tumoral" = "#FF6F00")


# Create the plot
p <- ggplot(data, aes(x = tipo, y = taxa, fill = grupo)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, alpha = 0.4) +
  geom_jitter(aes(shape = grupo, color = grupo), 
              position = position_jitterdodge(jitter.width = 0.2), size = 3) +
  labs(x = "Group", y = "Levine's Age Acceleration") +
  ggtitle("DNAm Age Acceleration") +
  scale_fill_manual(values = alpha(cores, 0.7), labels = c("Normal", "Tumoral")) +
  scale_shape_manual(values = c("Normal" = 16, "Tumoral" = 17)) +
  scale_color_manual(values = cores_jitter) +
  theme_bw() +
    theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    text = element_text(size = 12, color = "black", face = "plain"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

p

# Add the statistical test result in the plot
p <- p + stat_compare_means(
  aes(group = grupo),
  method = "wilcox.test",
  label = "p.format",
  label.x.npc = "left",
  label.y.npc = c(0.9, 0.7)
)

# Graphic visualization
p





#####################################################################################################################
#                                         S U R V I V A L     A N A L Y S I S  

#                          D I F F E R E N C E     B E T W E E N     N O R M A L  vs  T U M O R
#---------------------------------------------------------------------------------------------------------------------


# Create a new dataset with columns 'type' and 'ageacc2'
acc_comparison <- data.frame(
  Type = c(rep("Normal", length(surv_normal_tcga$ageAcc2.Levine)),
           rep("Tumoral", length(surv_tumor_tcga$ageAcc2.Levine))),
  AgeAcc2 = c(surv_normal_tcga$ageAcc2.Levine, surv_tumor_tcga$ageAcc2.Levine)
)


# Statistical test
test_result <- wilcox.test(AgeAcc2 ~ Type, data = acc_comparison)



# Plotting
ggplot(acc_comparison, aes(x = Type, y = AgeAcc2, fill = Type)) +
  geom_boxplot(alpha = 0.2, outlier.shape = NA) +
  geom_jitter(width = 0.4, height = 0, size = 2, alpha = 1, shape = 25) +
  labs(x = "", y = "Age Acceleration") +
  ggtitle("Age Acceleration: Normal vs Tumoral Tissue") +
  geom_text(aes(x = 1.5, y = max(AgeAcc2), 
                label = paste("p-value =", format(test_result$p.value, digits = 4))),
            vjust = -1.5, size = 3, color = "black", fontface = "plain") +
  theme_minimal() +
  scale_fill_viridis(discrete = TRUE, option = "turbo") +
  theme(
    plot.title = element_text(size = 13, face = "plain"),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 10),
    text = element_text(size = 12, color = "black", face = "plain")
  )



####################################################################################################################

# Subset of 'surv_tumor_tcga': I want to create a df for each group (ACC and DEC) to facilitate
# the comparison between these two

data_acc <- subset(surv_tumor_tcga, AgeAcceleration == "ACC")
data_dec <- subset(surv_tumor_tcga, AgeAcceleration == "DEC")

# Apply the statistyical test (Student's t-test)
t_test_result <- t.test(data_acc$Levine, data_dec$Levine)

print(t_test_result)


# Plot according to a variable of interest! Ex. smokingstatus, stage, gender...
# Boxplot
boxplot_plot <- ggplot(surv_tumor_tcga, aes(x = AgeAcceleration, y = Levine, fill = smokingstatus)) +
  geom_boxplot() +
  geom_boxplot(alpha = 0.2, outlier.shape = NA) +
  geom_jitter(width = 0.4, height = 0, size = 2, alpha = 1, shape = 25) +
  labs(x = "Group", y = "DNAm Age", fill = "Smoking Status", title = "Epigenetic Age Acceleration and Smoking Status") +
  scale_fill_manual(values = c("no-smoker" = "gray", "smoker" = "pink"), labels = c("Non-Smoker", "Smoker")) +
  theme_minimal()

boxplot_plot

#Statistical test: epigenetic age acceleration value according to the variable of interest
t <- t.test(Levine ~ smokingstatus, data = surv_tumor_tcga)
t
##mean in group no-smoker    mean in group smoker 
######92.58324                91.12818

getwd()




########################################################################################################

#                         F I N D I N G     D M P ' s     B E T W E E N      D N A m A G E    G R O U P S 
# EXPLORE DIFFERENCES BETWEEN 'ACC' AND 'DEC' GROUPS!
#--------------------------------------------------------------------------------------------------------

# First I need to go back to my original file 
met

# Filter this object to keep only tumor samples
pd_filtered <- met$pd %>% filter(shortLetterCode == "TP")


# Create a new column with the classification of each sample into "ACC" or "DEC"
pd_filtered$AgeAcc <- surv_tumor_tcga$AgeAcceleration[match(pd_filtered$barcode, surv_tumor_tcga$barcode)]

# As 'surv_tumor_tcga' have columns "barcode" and "AgeAcceleration"
# And 'pd_filtered' have columns  "barcode" and "AgeAcc"

# check for discrepancies: samples are correctly classified in pd_filtered?
discrepancies <- pd_filtered$AgeAcc != surv_tumor_tcga$AgeAcceleration[pd_filtered$barcode %in% surv_tumor_tcga$barcode]

# Show discrepeancies
samples_with_discrepancies <- pd_filtered$barcode[discrepancies]

# If 'samples_with_discrepancies' is empty, it means that all samples are correct
if (length(samples_with_discrepancies) == 0) {
  print("As amostras estão corretas.")
} else {
  print("Amostras incorretas:")
  print(samples_with_discrepancies)
}
##Ok!


# Store rownames of pd_filtered
col_names <- rownames(pd_filtered)

#Filter met$beta object to keep only samples contained in 'col_names' - only tumoral samples beta-values
beta_filtered <- met$beta[, col_names]

#create another 'met' object: i need both beta matrix and pd_filtered in a single list (same input as met)
met_filtered <- list(beta = beta_filtered, pd = pd_filtered)



#Saving files
#save(pd_filtered, beta_filtered, met_filtered, file = "TCGA_MethylationTumorData_AgeAcc.RData")

#If I was to plot something, I do not need to perform all the pipeline again, I can just load the saved data:
#load("~/TCGA_MethylationTumorData_AgeAcc.RData")


#############################################################################################################

##---------------   Run champ.DMP to get DMPs between groups 'ACC' and 'DEC'
dmp_tumor_ageacc <- champ.DMP(beta = met_filtered$beta, 
                              pheno = met_filtered$pd$AgeAcc, 
                              adjPVal = 0.001,
                              adjust.method = "BH")
#90897 obs --  p= 0.05


#39969 obs -- p= 0.001

dmp_dnamage <- as.data.frame(dmp_tumor_ageacc$DEC_to_ACC)


#convert rownames into a column named 'probes' 
dmp_dnam.age <- tibble::rownames_to_column(dmp_dnam.age, "probes")




DMP.GUI(DMP=dmp_tumor_ageacc[[1]],
        beta=met_filtered$beta,
        pheno=met_filtered$pd$AgeAcc)


##Filter to keep probes with a absolute deltabeta value of >= 0.3:
db_ageacc = with(dmp_tumor_ageacc$DEC_to_ACC, which(abs(dmp_tumor_ageacc$DEC_to_ACC$deltaBeta) >= 0.3))
db_ageacc_f = dmp_tumor_ageacc$DEC_to_ACC[db_ageacc,]
## 20 probes remained.


