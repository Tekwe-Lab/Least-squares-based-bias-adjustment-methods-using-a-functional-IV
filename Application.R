
library(dplyr); library(tidyr)
library(ggplot2); library(ggpubr)
library(tableone); library(kableExtra)
library(foreach); library(doParallel)
library(survey)

# Load cleaned dataset
load("NHANES_2005.Rda")




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Exploratory Data Analysis
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Demographic
demo = unique(nhanes_2005[,c("seqn","ridageyr","age_group","riagendr","ridreth1","ridreth1_regroup","dmdeduc2","dmdeduc2_regroup","db","bmxbmi","bmi_group","wtint2yr")])
demo.svy = svydesign(ids=~seqn, weights=~wtint2yr, data=demo)
tab1.svy = svyCreateTableOne(data=demo.svy, vars=c("ridageyr","age_group","riagendr","ridreth1","ridreth1_regroup","dmdeduc2","dmdeduc2_regroup","db","bmxbmi"), strata = "bmi_group", addOverall = TRUE)
tab1.svy %>% kableone(showAllLevels = TRUE) %>% kable_styling("striped", full_width = F)


# BMI
## Figure
fig.bmi = ggplot(demo, aes(x = bmxbmi)) + 
  geom_density(color = "black", fill = "grey", alpha = 0.2) + 
  labs(x = "Body Mass Index (BMI)", y = "") +
  scale_y_continuous(limits = c(0,0.08), breaks = c(0,0.02,0.04,0.06,0.08)) +
  scale_x_continuous(breaks = c(18.5,25,30)) +
  geom_vline(xintercept = 18.5, linetype="dashed", color="red") +
  geom_vline(xintercept = 25, linetype="dashed", color="red") +
  geom_vline(xintercept = 30, linetype="dashed", color="red") +
  annotate(geom="text", x=14, y=0.035, label = "Underweight", color="burlywood4", angle=90) +
  annotate(geom="text", x=21.5, y=0.035, label = "Normal", color="burlywood4", angle=90) +
  annotate(geom="text", x=27.5, y=0.035, label = "Overweight", color="burlywood4", angle=90) +
  annotate(geom="text", x=50, y=0.035, label = "Obese", color="burlywood4", angle=90) +
  theme_classic()
fig.bmi

pdf(file = "Application_FigBMI.pdf")
fig.bmi
dev.off()


# Physical Activity
## Figure
pa.mean = nhanes_2005 %>%
  group_by(paxhour) %>%
  summarise(Weekday.mean = mean(Weekday, na.rm = TRUE), Weekend.mean = mean(Weekend, na.rm = TRUE))

fig.pa.weekday = ggplot() +
  geom_line(data=nhanes_2005, aes(x=paxhour, y=Weekday, group=seqn), color="grey") +
  geom_line(data=pa.mean, aes(x=paxhour, y=Weekday.mean), color="deeppink3", size=1.25) +
  labs(x="Time (hour)",y="Weekday Physical Intensity") +
  ggtitle("(a) Weekday") +
  theme(axis.text=element_text(size=10),axis.title=element_text(size=12,face="bold"),
        plot.caption = element_text(size=14,face="bold",hjust = 0.5)) +
  theme_classic() +
  theme(panel.spacing = unit(2, "lines")) +
  theme(strip.background = element_rect(fill="cadetblue2", size=1.5, linetype="solid"))
fig.pa.weekend = ggplot() +
  geom_line(data=nhanes_2005, aes(x=paxhour, y=Weekend, group=seqn), color="grey") +
  geom_line(data=pa.mean, aes(x=paxhour, y=Weekend.mean), color="deeppink3", size=1.25) +
  labs(x="Time (hour)",y="Weekend Physical Intensity") +
  ggtitle("(b) Weekend") +
  theme(axis.text=element_text(size=10),axis.title=element_text(size=12,face="bold"),
        plot.caption = element_text(size=14,face="bold",hjust = 0.5)) +
  theme_classic() +
  theme(panel.spacing = unit(2, "lines")) +
  theme(strip.background = element_rect(fill="cadetblue2", size=1.5, linetype="solid"))

fig.pa = ggarrange(plotlist = list(fig.pa.weekday, fig.pa.weekend), nrow = 1, ncol = 2)
fig.pa

pdf(file = "Application_FigPA.pdf", width = 6, height = 3)
fig.pa
dev.off()

## Check zero
### Weekday
nhanes_2005 %>%
  group_by(paxhour) %>%
  summarise(count_zeros = sum(Weekday == 0)) %>%
  mutate(percent_zeros = 100*count_zeros/length(unique(nhanes_2005$seqn))) %>%
  kable("html") %>%
  kable_styling(full_width = F)

### Weekend
nhanes_2005 %>%
  group_by(paxhour) %>%
  summarise(count_zeros = sum(Weekend == 0)) %>%
  mutate(percent_zeros = 100*count_zeros/length(unique(nhanes_2005$seqn))) %>%
  kable("html") %>%
  kable_styling(full_width = F)




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Functional Linear Regression
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Setup parallel backend to use many processors
# cores = detectCores()
cl = makeCluster(30) #not to overload your computer cores[1]-1
registerDoParallel(cl)


# Main model
source("Application_LinearFuncIV.R")
result = linear_analysis(df=nhanes_2005 %>% filter(paxhour >= 7 & paxhour <= 22), 
                         ID="seqn", Y="bmxbmi", Time="paxhour", FV="Weekday", IV="Weekend", 
                         EF=c("ridageyr","riagendr.Female","ridreth1_regroup.Black","ridreth1_regroup.Hispanic","ridreth1_regroup.Other","db.1"),
                         wt="wtint2yr")


# Bootstrap
##Re-sampling with replacement
nboot = 500
id_boot_list = 
  foreach (icount(nboot)) %dopar% {
    id_list = unique(nhanes_2005$seqn)
    id_boot = sample(id_list, length(id_list), replace = TRUE)
    id_boot = as.data.frame(id_boot); colnames(id_boot) = "seqn"
    
    id_boot
  }

result_boot = 
  foreach (id_boot = id_boot_list, .errorhandling = 'remove') %dopar% {
    source("Application_LinearFuncIV.R")
    res_boot = linear_analysis_boot(df=nhanes_2005 %>% filter(paxhour >= 7 & paxhour <= 22), 
                                    ID="seqn", Y="bmxbmi", Time="paxhour", FV="Weekday", IV="Weekend", 
                                    EF=c("ridageyr","riagendr.Female","ridreth1_regroup.Black","ridreth1_regroup.Hispanic","ridreth1_regroup.Other","db.1"),
                                    wt="wtint2yr",
                                    id_boot=id_boot)
    
    res_boot
  }


# Stop cluster
stopCluster(cl)


# Function for main analysis output rearrangement
main_output_func = function(result, result_boot){
  # Functional Variable #
  # Beta
  result_W = result$beta_W
  
  # Beta percent change vs. Naive method - By time
  beta_W_percent_change_bytime = result_W %>%
    mutate(MLM = 100 * (MLM-Naive) / Naive,
           LM = 100 * (LM-Naive) / Naive,
           SIMEX = 100 * (SIMEX-Naive) / Naive) %>%
    dplyr::select(Time, MLM, LM, SIMEX)
  
  # Beta percent change vs. Naive method - Averaged
  beta_W_percent_change_avg = data.frame(MLM = mean(abs(beta_W_percent_change_bytime$MLM)),
                                         LM = mean(abs(beta_W_percent_change_bytime$LM)),
                                         SIMEX = mean(abs(beta_W_percent_change_bytime$SIMEX)))
  
  # Bootstrap
  result_boot_W = NULL
  for (i in 1:length(result_boot)){
    beta_W = result_boot[[i]]$beta_W %>% mutate(nboot = i)
    result_boot_W = rbind(result_boot_W, beta_W)
  }
  
  boot.lcl = result_boot_W %>%
    group_by(Time) %>%
    summarise_at(vars("Naive","SIMEX","LM","MLM"), .funs = list(~quantile(., probs = 0.025))) %>%
    mutate_at(c("Naive","SIMEX","LM","MLM"), as.numeric)
  boot.ucl = result_boot_W %>%
    group_by(Time) %>%
    summarise_at(vars("Naive","SIMEX","LM","MLM"), .funs = list(~quantile(., probs = 0.975))) %>%
    mutate_at(c("Naive","SIMEX","LM","MLM"), as.numeric)
  
  beta_W = result_W %>%
    full_join(boot.lcl %>% rename_at(vars(-Time), function(x) paste0(x,"_lcl"))) %>% 
    full_join(boot.ucl %>% rename_at(vars(-Time), function(x) paste0(x,"_ucl"))) %>%
    dplyr::select(Time, 
                  MLM, MLM_lcl, MLM_ucl, LM, LM_lcl, LM_ucl,
                  SIMEX, SIMEX_lcl, SIMEX_ucl, Naive, Naive_lcl, Naive_ucl)
  
  # Figure b
  fig_b.mlm = ggplot(data=beta_W) +
    geom_ribbon(aes(x=Time, ymin=MLM_lcl, ymax=MLM_ucl), fill="lightgray") +
    geom_line(aes(x=Time, y=MLM_lcl, colour="95% Bootstrap C.I."), size=1, linetype="dotdash") +
    geom_line(aes(x=Time, y=MLM_ucl, colour="95% Bootstrap C.I."), size=1, linetype="dotdash") +
    geom_line(aes(x=Time, y=MLM, color="MULTI-2SLS"), lwd=1) +
    geom_line(aes(x=Time, y=Naive, color="Naive"), lwd=1) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    scale_x_continuous(limits=c(7, 22), breaks=c(7, 12, 17, 22)) +
    scale_y_continuous(limits=c(-4e-04, 4.5e-04), breaks=c(-4e-04, -2e-04, 0, 2e-04, 4e-04)) +
    labs(x="Time", y=expression(widehat(beta)[1] ~ "(t)")) +
    scale_colour_manual("Legend",
                        breaks = c("95% Bootstrap C.I.", "MULTI-2SLS", "PW-2SLS", "SIMEX", "Naive"),
                        values = c("95% Bootstrap C.I." ="cyan", "MULTI-2SLS"="purple", "PW-2SLS"="blue1", "SIMEX"="green4", "Naive"="red")) +
    ggtitle("(a) MULTI-2SLS vs. Naive") +
    theme_classic() + 
    theme(legend.position="none")
  fig_b.lm = ggplot(data=beta_W) +
    geom_ribbon(aes(x=Time, ymin=LM_lcl, ymax=LM_ucl), fill="lightgray") +
    geom_line(aes(x=Time, y=LM_lcl, colour="95% Bootstrap C.I."), size=1, linetype="dotdash") +
    geom_line(aes(x=Time, y=LM_ucl, colour="95% Bootstrap C.I."), size=1, linetype="dotdash") +
    geom_line(aes(x=Time, y=LM, color="PW-2SLS"), lwd=1) +
    geom_line(aes(x=Time, y=Naive, color="Naive"), lwd=1) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    scale_x_continuous(limits=c(7, 22), breaks=c(7, 12, 17, 22)) +
    scale_y_continuous(limits=c(-4e-04, 4.5e-04), breaks=c(-4e-04, -2e-04, 0, 2e-04, 4e-04)) +
    labs(x="Time", y=expression(widehat(beta)[1] ~ "(t)")) +
    scale_colour_manual("Legend",
                        breaks = c("95% Bootstrap C.I.", "MULTI-2SLS", "PW-2SLS", "SIMEX", "Naive"),
                        values = c("95% Bootstrap C.I." ="cyan", "MULTI-2SLS"="purple", "PW-2SLS"="blue1", "SIMEX"="green4", "Naive"="red")) +
    ggtitle("(b) PW-2SLS vs. Naive") +
    theme_classic() + 
    theme(legend.position="none")
  fig_b.simex = ggplot(data=beta_W) +
    geom_ribbon(aes(x=Time, ymin=SIMEX_lcl, ymax=SIMEX_ucl), fill="lightgray") +
    geom_line(aes(x=Time, y=SIMEX_lcl, colour="95% Bootstrap C.I."), size=1, linetype="dotdash") +
    geom_line(aes(x=Time, y=SIMEX_ucl, colour="95% Bootstrap C.I."), size=1, linetype="dotdash") +
    geom_line(aes(x=Time, y=SIMEX, color="SIMEX"), lwd=1) +
    geom_line(aes(x=Time, y=Naive, color="Naive"), lwd=1) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    scale_x_continuous(limits=c(7, 22), breaks=c(7, 12, 17, 22)) +
    scale_y_continuous(limits=c(-4e-04, 4.5e-04), breaks=c(-4e-04, -2e-04, 0, 2e-04, 4e-04)) +
    labs(x="Time", y=expression(widehat(beta)[1] ~ "(t)")) +
    scale_colour_manual("Legend",
                        breaks = c("95% Bootstrap C.I.", "MULTI-2SLS", "PW-2SLS", "SIMEX", "Naive"),
                        values = c("95% Bootstrap C.I." ="cyan", "MULTI-2SLS"="purple", "PW-2SLS"="blue1", "SIMEX"="green4", "Naive"="red")) +
    ggtitle("(c) SIMEX vs. Naive") +
    theme_classic() + 
    theme(legend.position="none")
  
  fig_b.mlm.chng = ggplot(data=beta_W_percent_change_bytime) +
    geom_segment(aes(x=Time, xend=Time, y=0, yend=MLM), linetype = "dashed", color = "purple") +
    geom_point(aes(x=Time, y=MLM)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    scale_x_continuous(limits=c(7, 22), breaks=c(7, 12, 17, 22)) +
    scale_y_continuous(limits=c(-2500, 6000), breaks=c(-2000, 0, 2000, 4000, 6000)) +
    labs(x="Time", y="Percent Change Compared to Naive Method") +
    ggtitle("(d) MULTI-2SLS vs. Naive") +
    theme_classic()
  fig_b.lm.chng = ggplot(data=beta_W_percent_change_bytime) +
    geom_segment(aes(x=Time, xend=Time, y=0, yend=LM), linetype = "dashed", color = "blue1") +
    geom_point(aes(x=Time, y=LM)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    scale_x_continuous(limits=c(7, 22), breaks=c(7, 12, 17, 22)) +
    scale_y_continuous(limits=c(-2500, 6000), breaks=c(-2000, 0, 2000, 4000, 6000)) +
    labs(x="Time", y="Percent Change Compared to Naive Method") +
    ggtitle("(e) PW-2SLS vs. Naive") +
    theme_classic()
  fig_b.simex.chng = ggplot(data=beta_W_percent_change_bytime) +
    geom_segment(aes(x=Time, xend=Time, y=0, yend=SIMEX), linetype = "dashed", color = "green4") +
    geom_point(aes(x=Time, y=SIMEX)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    scale_x_continuous(limits=c(7, 22), breaks=c(7, 12, 17, 22)) +
    scale_y_continuous(limits=c(-2500, 6000), breaks=c(-2000, 0, 2000, 4000, 6000)) +
    labs(x="Time", y="Percent Change Compared to Naive Method") +
    ggtitle("(f) SIMEX vs. Naive") +
    theme_classic()
  
  fig_b = ggarrange(plotlist = list(fig_b.mlm, fig_b.lm, fig_b.simex,
                                    fig_b.mlm.chng, fig_b.lm.chng, fig_b.simex.chng), 
                    nrow = 2, ncol = 3)
  
  # Figure b legend
  fig_b_legend = ggplot(data=beta_W) +
    geom_ribbon(aes(x=Time, ymin=LM_lcl, ymax=LM_ucl), fill="lightgray") +
    geom_line(aes(x=Time, y=LM_lcl, colour="95% Bootstrap C.I."), size=1, linetype="dotdash") +
    geom_line(aes(x=Time, y=LM_ucl, colour="95% Bootstrap C.I."), size=1, linetype="dotdash") +
    geom_line(aes(x=Time, y=MLM, color="MULTI-2SLS"), lwd=1) +
    geom_line(aes(x=Time, y=LM, color="PW-2SLS"), lwd=1) +
    geom_line(aes(x=Time, y=SIMEX, color="SIMEX"), lwd=1) +
    geom_line(aes(x=Time, y=Naive, color="Naive"), lwd=1) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    labs(x="Time", y=expression(widehat(beta)[1] ~ "(t)")) +
    scale_colour_manual("Legend",
                        breaks = c("95% Bootstrap C.I.", "MULTI-2SLS", "PW-2SLS", "SIMEX", "Naive"),
                        values = c("95% Bootstrap C.I." ="cyan", "MULTI-2SLS"="purple", "PW-2SLS"="blue1", "SIMEX"="green4", "Naive"="red"))
  fig_b.legend = get_legend(fig_b_legend, position = "bottom")
  
  
  # Error-free Covariates #
  # Beta
  result_EF = result$beta_EF %>% mutate(variable = row.names(.))
  
  # Bootstrap
  result_boot_EF = NULL
  for (i in 1:length(result_boot)){
    beta_EF = result_boot[[i]]$beta_EF %>% mutate(variable = row.names(.), nboot = i)
    result_boot_EF = rbind(result_boot_EF, beta_EF)
  }
  
  boot.lcl = result_boot_EF %>%
    group_by(variable) %>%
    summarise_at(vars("Naive","SIMEX","LM","MLM"), .funs = list(~quantile(., probs = 0.025))) %>%
    mutate_at(c("Naive","SIMEX","LM","MLM"), as.numeric)
  boot.ucl = result_boot_EF %>%
    group_by(variable) %>%
    summarise_at(vars("Naive","SIMEX","LM","MLM"), .funs = list(~quantile(., probs = 0.975))) %>%
    mutate_at(c("Naive","SIMEX","LM","MLM"), as.numeric)
  
  beta_EF = result_EF %>%
    full_join(boot.lcl %>% rename_at(vars(-variable), function(x) paste0(x,"_lcl"))) %>% 
    full_join(boot.ucl %>% rename_at(vars(-variable), function(x) paste0(x,"_ucl"))) %>%
    dplyr::select(variable, 
                  MLM, MLM_lcl, MLM_ucl, LM, LM_lcl, LM_ucl, 
                  SIMEX, SIMEX_lcl, SIMEX_ucl, Naive, Naive_lcl, Naive_ucl)
  
  
  return(list(beta_W = beta_W,
              beta_W_percent_change_bytime = beta_W_percent_change_bytime,
              beta_W_percent_change_avg = beta_W_percent_change_avg,
              fig_b = fig_b,
              fig_b.legend = fig_b.legend,
              
              beta_EF = beta_EF))
}

output = main_output_func(result = result, result_boot = result_boot)


# Result for functional variable
## Percent change vs. Naive method
kable(output$beta_W_percent_change_avg, "html") %>% kable_styling(full_width = F)

## Figure
fig.beta1.b = ggarrange(plotlist = list(output$fig_b), legend = "bottom", legend.grob = output$fig_b.legend)
fig.beta1.b

pdf(file = "Application_FigBeta1.b.pdf", width = 12, height = 10)
fig.beta1.b
dev.off()


# Result for error-free covariates
## Beta coefficient
kable(output$beta_EF, "html") %>% kable_styling(full_width = F)



