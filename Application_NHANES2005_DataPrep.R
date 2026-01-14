
library(dplyr); library(tidyr)
library(Hmisc)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Data Preparation for 2005-2006
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Inclusion criteria:
## 1. 20 years or older
## 2. Have physical activity data for 24 hours every day
## 3. Have physical activity data for at least one day of weekday and one day of weekend


####### Physical Activity ######
pa.0506.raw = sasxport.get("Raw NHANES Data/paxraw_d.xpt")

# Sum to hour level
pa.0506.hour = pa.0506.raw %>%
  group_by(seqn,paxday,paxhour) %>%
  summarise(paxinten = sum(paxinten)) %>%
  ungroup()

# Remove outliers
# Larger than the 3rd quantile + $3*$IQR
sum.inten = summary(pa.0506.hour$paxinten)
out.inten = sum.inten[5] + 3*(sum.inten[5] - sum.inten[2])

pa.0506.hour = pa.0506.hour %>% filter(paxinten < out.inten)

# Keep participants with physical activity data for 24 hours and with at least one day of weekday and one day of weekend
pa.0506.hour = pa.0506.hour %>% 
  group_by(seqn,paxday) %>% 
  filter(length(unique(paxhour)) == 24) %>% 
  ungroup() %>% 
  mutate(paxday_group = case_when(paxday %in% c(2:6) ~ "Weekday", paxday %in% c(1,7) ~ "Weekend")) %>%
  group_by(seqn) %>% 
  filter(length(unique(paxday_group)) == 2) %>% 
  ungroup()

# Average across days
pa.0506.hour = pa.0506.hour %>%
  group_by(seqn,paxday_group,paxhour) %>%
  summarise(paxinten = mean(paxinten, na.rm = TRUE))

# Convert to wide format
pa.0506 = pa.0506.hour %>% 
  pivot_wider(names_from = c(paxday_group), values_from = c(paxinten)) %>%
  mutate(Weekday_rescale = Weekday/1000,
         Weekend_rescale = Weekend/1000)
# length(unique(pa.0506$seqn))


####### Demographic #######
demo.0506 = sasxport.get("Raw NHANES Data/DEMO_D.XPT") %>%
  filter(seqn %in% unique(pa.0506$seqn)) %>%
  filter(ridageyr >= 20) %>%
  select(seqn,ridageyr,riagendr,ridreth1,dmdeduc2,wtint2yr)


####### Diabetes #######
db.0506 = sasxport.get("Raw NHANES Data/DIQ_D.XPT") %>%
  filter(seqn %in% unique(demo.0506$seqn)) %>% 
  filter(diq010 == 1 | diq010 == 2) %>%  ##only keep yes or no results 
  mutate(db = ifelse(diq010 == 2, 0, 1)) %>%
  select(seqn,db)


####### Body Weight #######
bw.0506 = sasxport.get("Raw NHANES Data/BMX_D.XPT") %>%
  filter(seqn %in% unique(demo.0506$seqn)) %>% 
  filter(bmxbmi < 100) %>% ##remove the participants with extreme BMI
  select(seqn,bmxbmi)


####### Combined #######
df = na.omit(demo.0506 %>% inner_join(db.0506) %>% inner_join(bw.0506))


####### Final #######
nhanes_2005 = df %>%
  left_join(pa.0506) %>%
  mutate(seqn = as.integer(seqn), ridageyr = as.integer(ridageyr), riagendr = as.integer(riagendr), ridreth1 = as.integer(ridreth1), 
         wtint2yr = as.numeric(wtint2yr), bmxbmi = as.numeric(bmxbmi)) %>%
  mutate(age_group = ifelse(ridageyr<=44, "44 and younger", ifelse(ridageyr>=45 & ridageyr<=64, "45 to 64", "65 and older")),
         riagendr = factor(riagendr, levels = c(1,2), labels = c("Male","Female")),
         ridreth1 = factor(ridreth1, levels = c(1:5), labels = c("Mexican American","Other Hispanic","Non-Hispanic White",
                                                                 "Non-Hispanic Black","Other Race - Including Multi-Racial")),
         ridreth1_regroup = factor(case_when(ridreth1 %in% c("Mexican American","Other Hispanic") ~ 1,
                                             ridreth1 == "Non-Hispanic White" ~ 2,
                                             ridreth1 == "Non-Hispanic Black" ~ 3,
                                             ridreth1 == "Other Race - Including Multi-Racial" ~ 4),
                                   levels = c(1:4),
                                   labels = c("Hispanic","Non-Hispanic White",
                                              "Non-Hispanic Black","Other Race - Including Multi-Racial")),
         dmdeduc2 = factor(ifelse(dmdeduc2 %in% c(7,9), NA, dmdeduc2), levels = c(1:5), labels = c("Less Than 9th Grade","9-11th Grade","High School Grad/GED or Equivalent",
                                                                                                   "Some College or AA degree","College Graduate or above")),
         dmdeduc2_regroup = factor(case_when(dmdeduc2 %in% c("Less Than 9th Grade","9-11th Grade","High School Grad/GED or Equivalent") ~ 1,
                                             dmdeduc2 == "Some College or AA degree" ~ 2,
                                             dmdeduc2 == "College Graduate or above" ~ 3), 
                                   levels = c(1:3), 
                                   labels = c("High School Grad/GED or Lower","Some College or AA degree","College Graduate or above")),
         db = factor(db, levels = c(0,1), labels = c("No","Yes")),
         bmi_group = factor(case_when(bmxbmi<18.5 ~ 1,
                                      bmxbmi>=18.5 & bmxbmi<=24.9 ~ 2,
                                      bmxbmi>=25 & bmxbmi<=29.9 ~ 3,
                                      bmxbmi>=30 ~ 4),
                            levels = c(1:4),
                            labels = c("Underweight","Normal weight","Overweight","Obese"))) %>%
  mutate(riagendr.Female = ifelse(riagendr == "Female", 1, 0),
         ridreth1_regroup.Hispanic = ifelse(ridreth1 %in% c("Mexican American","Other Hispanic"), 1, 0),
         ridreth1_regroup.White = ifelse(ridreth1 == "Non-Hispanic White", 1, 0),
         ridreth1_regroup.Black = ifelse(ridreth1 == "Non-Hispanic Black", 1, 0),
         ridreth1_regroup.Other = ifelse(ridreth1 == "Other Race - Including Multi-Racial", 1, 0),
         db.1 = ifelse(db == "Yes", 1, 0)) %>%
  arrange(seqn,paxhour)
# length(unique(nhanes_2005$seqn))


save(nhanes_2005, file = "NHANES_2005.Rda")



