# LOADING THE LIBRARY -----------------------------------------------------

library(tidyverse)
library(data.table)
library(scales)

# GETTING THE NUMBERS -----------------------------------------------------
SCAN_SFS_df = tbl_df(read.csv("supplementary_file_1.csv"))
screened_swabs = nrow(SCAN_SFS_df)
handle_type = filter(SCAN_SFS_df, side_used == "handle")
handle_type_count = nrow(handle_type)
path_detect = filter(SCAN_SFS_df, !is.na(detected))
path_detect_count = nrow(path_detect)
rnasep_detect = filter(SCAN_SFS_df, !is.na(rnasep_detected))
rnasep_detect_count = nrow(rnasep_detect)

# FAILURE RATES -----------------------------------------------------------
# new dataframe for failure rate
failure_rate_df <- SCAN_SFS_df 
failure_rate_df <- select(SCAN_SFS_df, sample_id, side_used, rnasep_detected)
failure_rate_df  <- filter(failure_rate_df, !is.na(rnasep_detected)) 
# create table for failure rates
handle_swab_table <- table(failure_rate_df$side_used, failure_rate_df$rnasep_detected)
handle_swab_table_w_sum <- addmargins(table(failure_rate_df$side_used, 
                                            failure_rate_df$rnasep_detected), 2)
handle_swab_table
handle_swab_table_w_sum
# calculate p-value for failure rate
failure_rate_p_value <- fisher.test(handle_swab_table)
failure_rate_p_value

# REDCAP ANALYSIS ---------------------------------------------------------
# new dataframe for participants' data analysis
participant_df <- SCAN_SFS_df

# Confidence and Discomfort:
# create a tables of confidence
confidence_table <- table(participant_df$confidence, participant_df$side_used)
confidence_prop <- prop.table(table(participant_df$confidence, participant_df$side_used))
# create a tables of discomfort 
discomfort_table <- table(participant_df$discomfort, participant_df$side_used)
discomfort_prop <- prop.table(table(participant_df$discomfort, participant_df$side_used))
# calculate p-value for confidence of handle versus swab
fisher.test(participant_df$confidence, participant_df$side_used)
# calculate p-value for discomfort for handle versus swab
fisher.test(participant_df$discomfort, participant_df$side_used)
# binding data of confidence and discomfort together 
confidence_prop_table <- data.table(confidence_prop)
colnames(confidence_prop_table) <- c("confidence","side_used","freq")
discomfort_prop_table <- data.table(discomfort_prop)
colnames(discomfort_prop_table) <- c("discomfort","side_used","freq")
nasal_swab_usability_table <- bind_rows(confidence_prop_table,discomfort_prop_table)
# final table is saved as a csv
write_csv(nasal_swab_usability_table,"nasal_swab_usability_table.csv")
# create table with relative frequency 
rel_confidence_prop_table <- confidence_prop_table %>% group_by(side_used) %>% 
  mutate(relfreq = freq / sum(freq))
# rename column names in prep of bar plot generation
library(plyr)
rel_confidence_prop_table$confidence <- revalue(rel_confidence_prop_table$confidence, 
                                                c("not_con" = "Not \n confident", "some_con" = "Somewhat \n confident", 
                                                  "v_con" = "Very \n confident"))
detach("package:plyr", unload = TRUE)
rel_confidence_prop_table$side_used <- factor(rel_confidence_prop_table$side_used, levels = c("swab","handle"))
# create barplot for discomfort 
ggplot(data = rel_confidence_prop_table,aes(x=confidence,y=relfreq,fill=side_used))+
  geom_bar(stat="identity", position = "dodge") +
  scale_y_continuous(labels = scales::percent)+
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major.y = element_line(size = 0.5, linetype = 'solid', color = "grey"),
        panel.grid.major.x = element_blank(),
        legend.position = "none",
        plot.title = element_text(face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_text(face = "bold",vjust= -1),
        axis.ticks = element_blank(),
        axis.title.y = element_blank())+
  geom_text(aes(label=scales::percent(relfreq,accuracy =1),y=relfreq), vjust=-0.5, size=3.5,position = position_dodge(width=1))+
  ggtitle("Participant Confidence that Sample was \n Taken Correctly by Handle and Swab")+
  scale_fill_manual(values = c("swab" = "lightgrey",
                               "handle" = "#095DA8"))+
  labs(x="Confidence
       ", fill = "Side Used")

# removing unneeded dataframe
rm(confidence_prop_table)
# create table with relative frequency 
rel_discomfort_prop_table <- discomfort_prop_table %>% group_by(side_used) %>% 
  mutate(relfreq = freq / sum(freq))
# rename column names in prep of bar plot generation
library(plyr)
rel_discomfort_prop_table$discomfort <- revalue(rel_discomfort_prop_table$discomfort, 
                                                c("no_dis" = "No \n discomfort", 
                                                  "mild_dis" = "Mild \n discomfort", 
                                                  "strong_dis" = "Strong \n discomfort"))
detach("package:plyr", unload = TRUE)
rel_discomfort_prop_table$discomfort <- factor(rel_discomfort_prop_table$discomfort, 
                                               levels = c("No \n discomfort","Mild \n discomfort",
                                                          "Strong \n discomfort"))
rel_discomfort_prop_table$side_used <- factor (rel_discomfort_prop_table$side_used, levels = c("swab","handle"))
# create barplot for discomfort 
ggplot(data = rel_discomfort_prop_table,aes(x=discomfort,y=relfreq,fill=side_used))+
  geom_bar(stat="identity", position = "dodge")+
  scale_y_continuous(labels = scales::percent)+
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major.y = element_line(size = 0.5, linetype = 'solid', color = "grey"),
        panel.grid.major.x = element_blank(),
        legend.position = "none",
        plot.title = element_text(face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_text(face = "bold",vjust= -1),
        axis.ticks = element_blank(),
        axis.title.y = element_blank())+
  geom_text(aes(label=scales::percent(relfreq,accuracy =1),y=relfreq), vjust=-0.5, size=3.5,position = position_dodge(width=1))+
  ggtitle("Discomfort Experienced by Participants by \n Handle and Swab")+
  scale_fill_manual(values = c("swab" = "lightgrey",
                               "handle" = "#095DA8"))+
  labs(x="Discomfort", fill = "Side Used")
# removing unneeded dataframe
rm(discomfort_prop_table)

# Step 2 Analysis and Creation of Plot for Age
participant_df %>%
  dplyr::group_by(side_used) %>%
  dplyr::summarise(mean = mean(age,na.rm=TRUE),median = median(age,na.rm=TRUE), n = n()) 
# t test for difference in age
dfAge <- data.frame(side_used = participant_df$side_used, age = participant_df$age)
dfAge <- dfAge %>% mutate(id = row_number())
age_t_test <- dfAge %>% spread(side_used,age)
t.test(age_t_test$handle,age_t_test$swab)

length(!is.na(participant_df$age))
#Separate into chosen bins
agedata <- participant_df %>% filter (!is.na(age))
age_group <- agedata %>% mutate(age_group = case_when(
  age >= 65 ~ '65+',
  age >= 50 ~ '50-64',
  age >= 35 ~ '35-49',
  age >= 18 ~ '18-34',
  age < 18 ~ 'Under 18'))
age_group2 <- data.frame(age = age_group$age_group, side_used = age_group$side_used)
age_hist <- age_group2 %>%
  group_by(side_used,age,.drop=FALSE) %>%
  dplyr::summarise(n = n()) %>%
  complete(side_used, fill=list(count_a =0))%>%
  mutate(freq = n / sum(n))%>%
  ungroup()

age_hist$age <- factor(age_hist$age, levels = c("Under 18","18-34","35-49","50-64", "65+"))
age_hist$side_used <- factor (age_hist$side_used, levels = c("swab","handle"))
age_hist$percent <- scales::percent(age_hist$freq)

ggplot(age_hist, aes(x=age,y=freq,fill=side_used))+
  geom_bar(stat="identity",position = position_dodge(preserve="single"))+
  scale_y_continuous(labels = c("0%","10%","20%","30%","40%"), breaks = c(0,0.1,0.2,0.3,0.4))+
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major.y = element_line(size = 0.5, linetype = 'solid', color = "grey"),
        panel.grid.major.x = element_blank(),
        legend.position = "right",
        plot.title = element_text(face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_text(face = "bold",vjust= -1),
        axis.ticks = element_blank(),
        axis.title.y = element_blank())+
  geom_text(aes(label=scales::percent(freq,accuracy =1),y=freq), vjust=-0.5, size=3.5,position = position_dodge(width=1))+
  ggtitle("Age of Participants by Handle and Swab")+
  scale_fill_manual(values = c("swab" = "lightgrey",
                               "handle" = "#095DA8"))+
  labs(x="Age
       ", fill = "Side Used")


# Step 3: Analysis and Creation of Plot for Sex
sex_group <- participant_df %>% mutate(sex_group = case_when(
  sex == "male" ~ 'Male',
  sex == "female" ~ 'Female',
  sex == "other" ~ 'NA',
  sex == "dont_say" ~ 'NA'))

# NAs are removed 
sex_group <- sex_group[!is.na(sex_group$sex_group), ]
sex_group <- sex_group %>% filter(!sex_group %in% "NA")

# dataframe created and modified 
sex_hist <- sex_group %>%
  group_by(side_used,sex_group,.drop=FALSE) %>%
  dplyr::summarise(n = n()) %>%
  complete(side_used, fill=list(count_a =0))%>%
  mutate(freq = n / sum(n))%>%
  ungroup()
sex_hist$side_used <- factor(sex_hist$side_used, levels = c("swab","handle"))
ggplot(sex_hist, aes(x=sex_group,y=freq,fill=side_used))+
  geom_bar(stat="identity",position = position_dodge(preserve="single"))+
  scale_y_continuous(labels = scales::percent)+
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major.y = element_line(size = 0.5, linetype = 'solid', color = "grey"),
        panel.grid.major.x = element_blank(),
        legend.position = "right",
        plot.title = element_text(face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_text(face = "bold",vjust= -1),
        axis.ticks = element_blank(),
        axis.title.y = element_blank())+
  geom_text(aes(label=scales::percent(freq,accuracy =1),y=freq), vjust=-0.5, size=3.5,position = position_dodge(width=1))+
  ggtitle("Sex of Participants by \n Handle and Swab")+
  scale_fill_manual(values = c("swab" = "lightgrey",
                               "handle" = "#095DA8"))+
  labs(x="Sex", fill = "Side Used")

# create a 2x2 dataframe in prep of chi-square test 
dfSex <- data.frame(side_used = sex_group$side_used, sex = sex_group$sex_group)
# chi-square test
chisq.test(dfSex$sex,dfSex$side_used)
# display table
table(dfSex)

# Step 4: Analysis and Creation of Plot for Income
# new dataframe created 
dfIncome <- participant_df
dfIncome <- dfIncome %>% filter (!is.na(income))
dfIncome %>% group_by(income) %>% summarise()


# income binnings are created 
income_group <- dfIncome %>% mutate(income_group = case_when(
  income == "100k_125k" ~ '100-125k',
  income == "125k_150k" ~ '125-150k',
  income == "25k_50k" ~ '25-50k',
  income == "50k_75k" ~ '50-75k',
  income == "75k_100k" ~ '75-100k',
  income == "dont_know" ~ "NA",
  income == "dont_say" ~ "NA",
  income == "less_25k" ~ '<25k',
  income == "more_150k" ~ '>150k',))

# p-value calculated for significance
chisq.test(dfIncome$income,dfIncome$side_used)

income_hist <- income_group %>%
  group_by(side_used,income_group,.drop=FALSE) %>%
  dplyr::summarise(n = n()) %>%
  complete(side_used, fill=list(count_a =0))%>%
  mutate(freq = n / sum(n))%>%
  ungroup()
sum(income_hist$freq)
# QC-ing of the dataframe (factor)
income_group %>% group_by(income_group) %>% summarise()
income_hist$income_group <- factor(income_hist$income_group, levels =
                                     c("<25k","25-50k","50-75k","75-100k","100-125k","125-150k",">150k", "NA"))
income_hist$side_used <- factor (income_hist$side_used, levels = c("swab","handle"))
ggplot(income_hist, aes(x=income_group,y=freq,fill=side_used, width = 0.8))+
  geom_bar(stat="identity",position = position_dodge(preserve="single"))+
  scale_y_continuous(labels = scales::percent)+
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major.y = element_line(size = 0.5, linetype = 'solid', color = "grey"),
        panel.grid.major.x = element_blank(),
        legend.position = "right",
        plot.title = element_text(face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 10, angle = 35, vjust=1),
        axis.title.x = element_text(face = "bold",vjust= -1),
        axis.ticks = element_blank(),
        axis.title.y = element_blank())+
  geom_text(aes(label=scales::percent(freq,accuracy =1),y=freq), vjust=-0.5, size=3.5,
            position = position_dodge(width=.8))+
  ggtitle("Household Income of Participants by Handle and Swab")+
  scale_fill_manual(values = c("swab" = "lightgrey",
                               "handle" = "#095DA8"))+
  labs(x="Income", fill = "Side Used")


# MEAN RNASEP CT VALUES ---------------------------------------------------
crt_data <- SCAN_SFS_df

# create and modify dataframe to include only 4 columns
crt_data <- crt_data %>% subset(select = c("side_used","crt","sample_id","rack_id"))
# frop the NAs from the crt column
crt_data <- drop_na(crt_data, crt)
# convert to a data table
crt_data <- as.data.table(crt_data)
# new dataframe created where crt values are assigned, dependent on type
crt_cast <- dcast(crt_data, sample_id ~ side_used, value.var = "crt", fun.aggregate = sum)
# Values of 0 casted as NA
crt_cast$handle <- na_if(crt_cast$handle,0)
crt_cast$swab <- na_if(crt_cast$swab,0)
# calculate mean RNAse P CRTs 
handle_x <- mean(crt_cast$handle,na.rm = TRUE)
swab_x <- mean(crt_cast$swab, na.rm = TRUE)

# calculate 95% confidence intervals 
handle_ci = qnorm(0.975)*sd(na.omit(crt_cast$handle))/sqrt(length(na.omit(crt_cast$handle)))#[[1,1]]
swab_ci <- qnorm(0.975)*sd(na.omit(crt_cast$swab))/sqrt(length(na.omit(crt_cast$swab)))#[[1,1]]

# to find the upper and lower ranges for the mean based on the 95% CI
handle_lower = handle_x - handle_ci
handle_upper = handle_x + handle_ci
swab_lower = swab_x - swab_ci
swab_upper = swab_x + swab_ci
t.test(crt_cast$handle,crt_cast$swab)

# calculate N's
crt_data %>% group_by(side_used) %>% summarise(n=n())

# box plot 
boxplot <- crt_data %>%
  ggplot(aes(x=side_used,y=crt,fill=side_used))+
  geom_boxplot()+
  scale_fill_manual(values = c("swab" = "lightgrey",
                               "handle" = "#095DA8"))+
  theme_light()+
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major.y = element_line(size = 0.5, linetype = 'solid', color = "grey"),
        legend.position = "none",
        plot.title = element_text(face = "bold", hjust = 0.5),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        axis.text = element_text(size = 11))+
  ggtitle("Human RNAse P CRT Values by Handle and Swab")+
  scale_y_continuous("CRT", breaks = c(10,15,20,25,30), labels = c("10","15","20","25","30"))+
  geom_hline(yintercept=28, linetype = "longdash",color="#4D4D4D")+
  #  annotate("rect", xmin= 1.26,xmax=1.74,ymin=27.5,ymax=28.5,fill="white")+
  #  annotate("text", y=28,x=1.5, label="Limit of Detection", color = "black")+
  annotate("text",x=1.3,y=21.7, label = "N=99")+
  annotate("text",x=2.27,y=18.8, label = "N=11904")+
  xlab("Side used")+
  ylab("CRT")
print(boxplot)

# Dot Plot
# rack id of the samples that was collected via handle
handle_racks <- data.table(rack_id = na.omit(unique(ifelse(crt_data$side_used=="handle",crt_data$rack_id,NA))))

# rack id and crt of samples collected via swab
swabs_within_handle_racks <- left_join(handle_racks,crt_data,by="rack_id")
swabs_within_handle_racks <- swabs_within_handle_racks[grepl("swab", swabs_within_handle_racks$side_used),]

# rack id and crt of samples collected via handle
handles_within_handle_racks <- crt_data[grepl("handle", crt_data$side_used),]

# dot plot generation
p <- ggplot()+
  geom_point(
    data=swabs_within_handle_racks,
    aes(x=as.character(rack_id),y=crt,color="swab"),
    size=1)+
  geom_point(
    data=handles_within_handle_racks,
    aes(x=as.character(rack_id),y=crt,color="handle"),
    size=1)+
  theme_bw()+
  theme(panel.grid.minor.y = element_blank(),
        panel.border = element_rect(color="lightgrey"),
        plot.title = element_text(face = "bold", hjust = 0.5),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        axis.text.x = element_text(face = "bold",size = 7),
        axis.ticks = element_blank(),
        #        plot.margin = unit(c(0.2,1,0.2,0.2), "in"),
        #        legend.position = c(1.07,0.5),
        #        legend.background = element_rect(color = "lightgrey"),
        axis.line.x = element_blank())+
  ylim(10,30)+
  scale_x_discrete(breaks = c("021","054","100","133","206"), labels = c("Dec","Jan","Feb","Mar","Apr"))+
  ggtitle("Human RNAse P CRT Values by Handle and Swab")+
  #  geom_text(aes(x="21",y=9.5, label="Dec 2019"))+
  scale_color_manual(name = "Side Used", values = c("swab" = "lightgrey",
                                                    "handle" = "#095DA8"))+
  xlab("Date")+
  ylab("CRT")
print(p) #Print PDF as 3" by 6.5"

# PATHOGEN DETECTION ------------------------------------------------------
path_detect_df  <- SCAN_SFS_df
# Step 1b: Abridge dataframe
path_detect_df <- select(path_detect_df, -c(rack_id, crt, sex, confidence, discomfort, age, income, aliquoted_date))
# Step 2a: Filter data to remove samples with Rnase P failures
ctrldet <- path_detect_df %>% filter(!is.na(rnasep_detected))
ctrldet <- ctrldet %>% filter(rnasep_detected == "Detected")
# Step 3a: Create table to display number of detected/not detected by pathogen
detcol <- colnames(ctrldet[4:19])
TFdetected <- data.frame("Pathogen" = detcol, "Handle_Present" = 0, 
                         "Handle_Not_Present" = 0, "Swab_Present" = 0, 
                         "Swab_Not_Present" = 0, stringsAsFactors = FALSE)
# Step 3b: Group samples by handles vs. swabs
ctrldetH <- ctrldet %>% filter(side_used=="handle")
ctrldetS <- ctrldet %>% filter(side_used=="swab")
# Step 3c: Determine number of samples tested for SARS-CoV-2 and on Open Array
sars_set <- ctrldet %>% filter(!is.na(ctrldet$SARSCOV2))
sars_set_H <- sars_set %>% filter(side_used=="handle")
sars_set_S <- sars_set %>% filter(side_used=="swab")
oa_set <- ctrldet %>% filter(!is.na(ctrldet$Adenovirus))
oa_set_H <- oa_set %>% filter(side_used=="handle")
oa_set_S <- oa_set %>% filter(side_used=="swab")
# Step 3d: Fill table with numbers of detection for each pathogen for handles vs. swabs
for(row in 1:nrow(TFdetected)){
  TFdetected$Handle_Present[row] <- nrow(oa_set_H[oa_set_H[,TFdetected$Pathogen[row]]==TRUE,])
  TFdetected$Handle_Not_Present[row] <- nrow(oa_set_H[oa_set_H[,TFdetected$Pathogen[row]]==FALSE,])
  TFdetected$Swab_Present[row] <- nrow(oa_set_S[oa_set_S[,TFdetected$Pathogen[row]]==TRUE,])
  TFdetected$Swab_Not_Present[row] <- nrow(oa_set_S[oa_set_S[,TFdetected$Pathogen[row]]==FALSE,])
}
# Step 3e: Adjust SARS-CoV-2 numbers using totals found in Step 3c
TFdetected$Handle_Present[TFdetected$Pathogen=="SARSCOV2"] <- nrow(sars_set_H[sars_set_H$SARSCOV2==TRUE,])
TFdetected$Handle_Not_Present[TFdetected$Pathogen=="SARSCOV2"] <- nrow(sars_set_H[sars_set_H$SARSCOV2==FALSE,])
TFdetected$Swab_Present[TFdetected$Pathogen=="SARSCOV2"] <- nrow(sars_set_S[sars_set_S$SARSCOV2==TRUE,])
TFdetected$Swab_Not_Present[TFdetected$Pathogen=="SARSCOV2"] <- nrow(sars_set_S[sars_set_S$SARSCOV2==FALSE,])
# Step 3f: Add decimal form of rates of detection for each pathogen for handles vs. swabs
TFdetected <- TFdetected %>% mutate(dec_rate_H = Handle_Present/nrow(oa_set_H))
TFdetected <- TFdetected %>% mutate(dec_rate_S = Swab_Present/nrow(oa_set_S))
# Step 3g: Adjust SARS-CoV-2 rates based on totals found in Step 3c
TFdetected$dec_rate_H[TFdetected$Pathogen=="SARSCOV2"] <- TFdetected$Handle_Present[TFdetected$Pathogen=="SARSCOV2"]/nrow(sars_set_H)
TFdetected$dec_rate_S[TFdetected$Pathogen=="SARSCOV2"] <- TFdetected$Swab_Present[TFdetected$Pathogen=="SARSCOV2"]/nrow(sars_set_S)
# Step 4a: Run paired t-test on rates of detection
t.test(as.numeric(TFdetected$dec_rate_H),as.numeric(TFdetected$dec_rate_S),paired=TRUE)
# Step 5a: Rename pathogens and rate columns for new table
TFdetected2 <- TFdetected
TFdetected2$Pathogen[TFdetected2$Pathogen=="Flu_A"] <- "Influenza A"
TFdetected2$Pathogen[TFdetected2$Pathogen=="Flu_B"] <- "Influenza B"
TFdetected2$Pathogen[TFdetected2$Pathogen=="Flu_C"] <- "Influenza C"
TFdetected2$Pathogen[TFdetected2$Pathogen=="RSV"] <- "Respiratory Syncytial Virus"
TFdetected2$Pathogen[TFdetected2$Pathogen=="SARSCOV2"] <- "SARS-CoV-2*"
TFdetected2$Pathogen[TFdetected2$Pathogen=="Seasonal_Coronavirus"] <- "Seasonal Coronavirus"
TFdetected2$Pathogen[TFdetected2$Pathogen=="C.Pneumoniae"] <- "C. pneumoniae"
TFdetected2$Pathogen[TFdetected2$Pathogen=="S.Pneumoniae"] <- "S. pneumoniae"
TFdetected2$Pathogen[TFdetected2$Pathogen=="M.Pneumoniae"] <- "M. pneumoniae"
names(TFdetected2)[names(TFdetected2) == "dec_rate_H"] <- "pct_rate_H"
names(TFdetected2)[names(TFdetected2) == "dec_rate_S"] <- "pct_rate_S"
# Step 5b: Create percent function
percent <- function(x, digits = 1, format = "f", ...) {
  paste0(formatC(100 * x, format = format, digits = digits, ...), "%")
}
# Step 5c: Fill in percentages on TFdetected2 table
for (row in 1:nrow(TFdetected)){
  TFdetected2$pct_rate_H[row] <- percent(TFdetected2$Handle_Present[row]/nrow(oa_set_H))
  TFdetected2$pct_rate_S[row] <- percent(TFdetected2$Swab_Present[row]/nrow(oa_set_S))
}
# Step 5d: Recalculate SARS-CoV-2 rates based on totals found in Step 3c
TFdetected2$pct_rate_H[TFdetected2$Pathogen=="SARS-CoV-2*"] <- percent(TFdetected2$Handle_Present[TFdetected2$Pathogen=="SARS-CoV-2*"]/nrow(sars_set_H))
TFdetected2$pct_rate_S[TFdetected2$Pathogen=="SARS-CoV-2*"] <- percent(TFdetected2$Swab_Present[TFdetected2$Pathogen=="SARS-CoV-2*"]/nrow(sars_set_S))
# Step 5e: Create new table for percentages
detcol2 <- TFdetected2$Pathogen
pcttable <- data.frame("Pathogen"= detcol2)
# Step 5f: Fill in values from detection rate table and format for publication
pcttable <- pcttable %>% mutate(Handle_Present = paste(paste(TFdetected2$Handle_Present, TFdetected2$pct_rate_H,sep=" ("), "",sep=")"))
pcttable <- pcttable %>% mutate(Swab_Present = paste(paste(TFdetected2$Swab_Present, TFdetected2$pct_rate_S,sep=" ("), "",sep=")"))
# Step 5g: Sort table by pathogen (alphabetically) and move bacteria to the bottom
pcttable <- pcttable[order(pcttable$Pathogen),]
table_top <- pcttable[-grep("pneumoniae",pcttable$Pathogen),]
table_bottom <- pcttable[grep("pneumoniae",pcttable$Pathogen),]
final_pcttable <- rbind(table_top, table_bottom)
# Step 5h: Display table from Step 5g and write .csv with data from table
view(final_pcttable)
write.csv(pcttable, 'pcttable.csv', row.names = FALSE)