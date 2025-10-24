#Baseline Microfilarial Prevalence = BMP
time <- (seq(1, length(list2[[2]]), by=1)/(366)) + 1840 #run the list2 below before running this line
library(ggpubr)
library(patchwork) #to be able to: a2 + a6

#1. Savanes non-SIZ (Tone and Tandjoare perfectures)
#Now add the misplaced surveys of hyperendemic villages
#Mesoendemic
list2 <- lapply(means.examplesSNSO, `[[`, 2)
list22 <- lapply(means.examplesSNSR, `[[`, 2)
list222 <- lapply(means.examplesSNSP, `[[`, 2)
#I pushed the values a year behind of the surveyed mesoendemic village, as vector control started in 1976 (instead of 1977) and makes a considerable impact
a2 <- ggplot()+geom_line(aes(time,list2[[2]]*100), colour="lightblue", size=1)+geom_line(aes(time,list22[[2]]*100), colour="darkblue", size=1)+geom_line(aes(time,list222[[2]]*100), colour="lightblue", size=1)+ coord_cartesian(xlim = c(1972.75, 2027.28), ylim = c(4.52,95.5)) + theme_classic() +labs(title="d. Savanes non-SIZ - 50% BMP",  x ="Year", y = "Microfilarial prevalence (%)") + geom_errorbar(data = SavanesNo2, aes(ymin=CILower, ymax=CIUpper, x = Survey_year-1), width=.1, size=0.5)+ geom_point(data = SavanesNo2,aes(x = Survey_year-1, y = PREV_cr,colour=SavanesNo2$ENDE), size = 4)+scale_colour_manual(values=c("orange","dark red"))+   labs(colour=NULL)+geom_vline(xintercept=1977,colour="green") + geom_vline(xintercept=1994,colour="dark green")+geom_vline(xintercept=1991,colour="red")+geom_text(mapping=aes(x=1991, y=80, label="MDA start"), size=4, angle=90, vjust=-0.4, hjust=0)+geom_text(mapping=aes(x=1994, y=80, label="VC end"), size=4, angle=90, vjust=-0.4, hjust=0)+geom_text(mapping=aes(x=1977, y=80, label="VC start"), size=4, angle=90, vjust=-0.4, hjust=0)+theme(
  axis.title.x = element_text(size = 14),  # Increase x-axis title size
  axis.title.y = element_text(size = 14),  # Increase y-axis title size
  axis.text.x = element_text(size = 12),   # Increase x-axis text size
  axis.text.y = element_text(size = 12)    # Increase y-axis text size
)+scale_y_continuous(breaks = seq(0, 100, by = 10)) +  # Y-axis ticks every 10%
  scale_x_continuous(breaks = seq(1970, 2030, by = 5))  # X-axis ticks every 5 years

#Hyperendemic + 100% Vector control efficacy
library(dplyr)
#Hyper 100%VC
list9 <- lapply(means.examplesSNSO100, `[[`, 3)
list99 <- lapply(means.examplesSNSR100, `[[`, 3)
list999 <- lapply(means.examplesSNSP100, `[[`, 3)
a6 <- ggplot()+geom_line(aes(time,list9[[2]]*100), colour="lightblue", size=1)+geom_line(aes(time,list99[[2]]*100), colour="darkblue", size=1)+geom_line(aes(time,list999[[2]]*100), colour="lightblue", size=1)+ coord_cartesian(xlim = c(1972.75, 2027.28), ylim = c(4.52,95.5)) + theme_classic() +labs(title="e. Savanes non-SIZ - 70% BMP, 100% VC efficacy",  x ="Year", y = "Microfilarial prevalence (%)")+ geom_errorbar(data = SavanesHyperNo1, aes(ymin=CILower, ymax=CIUpper, x = Survey_year+3), width=.1, size=0.5) + geom_point(data = SavanesHyperNo1,aes(x = Survey_year+3, y = PREV_cr,colour=SavanesHyperNo1$ENDE), size = 4) +scale_colour_manual(values=c("dark red"))+   labs(colour=NULL)+geom_vline(xintercept=1977,colour="green") + geom_vline(xintercept=1994,colour="dark green")+geom_vline(xintercept=1991,colour="red")+geom_text(mapping=aes(x=1991, y=80, label="MDA start"), size=4, angle=90, vjust=-0.4, hjust=0)+geom_text(mapping=aes(x=1994, y=80, label="VC end"), size=4, angle=90, vjust=-0.4, hjust=0)+geom_text(mapping=aes(x=1977, y=80, label="VC start"), size=4, angle=90, vjust=-0.4, hjust=0)+theme(
  axis.title.x = element_text(size = 14),  # Increase x-axis title size
  axis.title.y = element_text(size = 14),  # Increase y-axis title size
  axis.text.x = element_text(size = 12),   # Increase x-axis text size
  axis.text.y = element_text(size = 12)    # Increase y-axis text size
)+scale_y_continuous(breaks = seq(0, 100, by = 10)) +  # Y-axis ticks every 10%
  scale_x_continuous(breaks = seq(1970, 2030, by = 5))  # X-axis ticks every 5 years

time <- (seq(1, length(list2[[2]]), by=1)/(366)) + 1840 #for plotting the model

#2. Savanes SIZ 
#100% VC efficacy + biannual
#Hypoendemicity
list3 <- lapply(means.examplesSSO100, `[[`, 1)
list33 <- lapply(means.examplesSSR100, `[[`, 1)
list333 <- lapply(means.examplesSSP100, `[[`, 1)
a3 <- ggplot()+geom_line(aes(time,list3[[2]]*100), colour="lightblue", size=1)+geom_line(aes(time,list33[[2]]*100), colour="darkblue", size=1)+geom_line(aes(time,list333[[2]]*100), colour="lightblue", size=1)+ coord_cartesian(xlim = c(1972.75, 2027.28), ylim = c(4.52,95.5)) + theme_classic() +labs(title="a. Savanes SIZ - 30% BMP, 100% VC efficacy, biannual MDA",  x ="Year", y = "Microfilarial prevalence (%)") + geom_errorbar(data = SavanesYes1B, aes(ymin=CILower, ymax=CIUpper, x = Survey_year), width=.1, size=0.5)+ geom_point(data = SavanesYes1B,aes(x = Survey_year, y = PREV_cr,colour=SavanesYes1B$ENDE), size = 4)+scale_colour_manual(values=c("yellow"))+   labs(colour=NULL)+geom_vline(xintercept=1977,colour="green") + geom_vline(xintercept=1993,colour="dark green")+geom_vline(xintercept=1991,colour="red")+geom_vline(xintercept=2003,colour="ORANGE")+geom_text(mapping=aes(x=1991, y=80, label="MDA start"), size=4, angle=90, vjust=-0.4, hjust=0)+geom_text(mapping=aes(x=1993, y=80, label="VC end"), size=4, angle=90, vjust=-0.4, hjust=0)+geom_text(mapping=aes(x=1977, y=80, label="VC start"), size=4, angle=90, vjust=-0.4, hjust=0)+geom_text(mapping=aes(x=2003, y=65, label="Biannual MDA start"), size=4, angle=90, vjust=-0.4, hjust=0)+ theme(
  axis.title.x = element_text(size = 14),  # Increase x-axis title size
  axis.title.y = element_text(size = 14),  # Increase y-axis title size
  axis.text.x = element_text(size = 12),   # Increase x-axis text size
  axis.text.y = element_text(size = 12)    # Increase y-axis text size
)+scale_y_continuous(breaks = seq(0, 100, by = 10)) +  # Y-axis ticks every 10%
  scale_x_continuous(breaks = seq(1970, 2030, by = 5))  # X-axis ticks every 5 years
#Mesoendemicity biannual
list5 <- lapply(means.examplesSSOB, `[[`, 2)
list55 <- lapply(means.examplesSSRB, `[[`, 2)
list555 <- lapply(means.examplesSSPB, `[[`, 2)
a5 <- ggplot()+geom_line(aes(time,list5[[2]]*100), colour="lightblue", size=1)+geom_line(aes(time,list55[[2]]*100), colour="darkblue", size=1)+geom_line(aes(time,list555[[2]]*100), colour="lightblue", size=1)+ coord_cartesian(xlim = c(1972.75, 2027.28), ylim = c(4.52,95.5)) + theme_classic() +labs(title="b. Savanes SIZ - 50% BMP, biannual MDA",  x ="Year", y = "Microfilarial prevalence (%)") + geom_errorbar(data = PANGA_dataset, aes(ymin=CILower, ymax=CIUpper, x = Survey_year), width=.1, size=0.5)+ geom_point(data = PANGA_dataset,aes(x = Survey_year, y = PREV_cr,colour=PANGA_dataset$ENDE), size = 4)+scale_colour_manual(values=c("orange"))+   labs(colour=NULL)+geom_vline(xintercept=1977,colour="green") + geom_vline(xintercept=1993,colour="dark green")+geom_vline(xintercept=1991,colour="red")+geom_text(mapping=aes(x=1991, y=80, label="MDA start"), size=4, angle=90, vjust=-0.4, hjust=0)+geom_text(mapping=aes(x=1993, y=80, label="VC end"), size=4, angle=90, vjust=-0.4, hjust=0)+geom_text(mapping=aes(x=1977, y=80, label="VC start"), size=4, angle=90, vjust=-0.4, hjust=0)+geom_text(mapping=aes(x=2003, y=65, label="Biannual MDA start"), size=4, angle=90, vjust=-0.4, hjust=0)+geom_vline(xintercept=2003,colour="ORANGE") +theme(
  axis.title.x = element_text(size = 14),  # Increase x-axis title size
  axis.title.y = element_text(size = 14),  # Increase y-axis title size
  axis.text.x = element_text(size = 12),   # Increase x-axis text size
  axis.text.y = element_text(size = 12)    # Increase y-axis text size
)+scale_y_continuous(breaks = seq(0, 100, by = 10)) +  # Y-axis ticks every 10%
  scale_x_continuous(breaks = seq(1970, 2030, by = 5))  # X-axis ticks every 5 years


#Mesoendemicity biannual + 100%VC
list4 <- lapply(means.examplesSSO100, `[[`, 2)
list44 <- lapply(means.examplesSSR100, `[[`, 2)
list444 <- lapply(means.examplesSSP100, `[[`, 2)
a4 <- ggplot()+geom_line(aes(time,list4[[2]]*100), colour="lightblue", size=1)+geom_line(aes(time,list44[[2]]*100), colour="darkblue", size=1)+geom_line(aes(time,list444[[2]]*100), colour="lightblue", size=1)+ coord_cartesian(xlim = c(1972.75, 2027.28), ylim = c(4.52,95.5)) + theme_classic() +labs(title="c. Savanes SIZ - 50% BMP, 100% VC efficacy, biannual MDA",  x ="Year", y = "Microfilarial prevalence (%)") + geom_errorbar(data = SavanesYes2100, aes(ymin=CILower, ymax=CIUpper, x = Survey_year), width=.1, size=0.5)+ geom_point(data = SavanesYes2100,aes(x = Survey_year, y = PREV_cr,colour=ENDE), size = 4)+scale_colour_manual(values=c("orange"))+   labs(colour=NULL)+geom_vline(xintercept=1977,colour="green") + geom_vline(xintercept=1993,colour="dark green")+geom_vline(xintercept=1991,colour="red")+geom_vline(xintercept=2003,colour="ORANGE")+geom_text(mapping=aes(x=1991, y=80, label="MDA start"), size=4, angle=90, vjust=-0.4, hjust=0)+geom_text(mapping=aes(x=1993, y=80, label="VC end"), size=4, angle=90, vjust=-0.4, hjust=0)+geom_text(mapping=aes(x=1977, y=80, label="VC start"), size=4, angle=90, vjust=-0.4, hjust=0)+geom_text(mapping=aes(x=2003, y=65, label="Biannual MDA start"), size=4, angle=90, vjust=-0.4, hjust=0) +theme(
  axis.title.x = element_text(size = 14),  # Increase x-axis title size
  axis.title.y = element_text(size = 14),  # Increase y-axis title size
  axis.text.x = element_text(size = 12),   # Increase x-axis text size
  axis.text.y = element_text(size = 12)    # Increase y-axis text size
)+scale_y_continuous(breaks = seq(0, 100, by = 10)) +  # Y-axis ticks every 10%
  scale_x_continuous(breaks = seq(1970, 2030, by = 5))  # X-axis ticks every 5 years

combined_plot <-  (a3+a5)/  (a4 +a2) /( a6+ plot_spacer())
theme_set(theme_minimal(base_family = "Arial"))
ggsave(filename = "Figure 2.pdf", plot = combined_plot, width = 15, height = 21, units = "in", dpi = 300)



#3. Kara SIZ 
#Hypoendemicity under biannual MDA
list1 <- lapply(means.examplesKSOB, `[[`, 1) #A list with 4 numeric columns (plus the fifth). We just need to select the second one, that corresponds to the prevalence
list11 <- lapply(means.examplesKSRB, `[[`, 1) #A list with 4 numeric columns (plus the fifth). We just need to select the second one, that corresponds to the prevalence
list111 <- lapply(means.examplesKSPB, `[[`, 1) #A list with 4 numeric columns (plus the fifth). We just need to select the second one, that corresponds to the prevalence
a1 <- ggplot()+geom_line(aes(time,list1[[2]]*100), colour="lightblue", size=1)+geom_line(aes(time,list11[[2]]*100), colour="darkblue", size=1)+geom_line(aes(time,list111[[2]]*100), colour="lightblue", size=1)+ coord_cartesian(xlim = c(1972.75, 2027.28), ylim = c(4.52,95.5)) + theme_classic() +labs(title="a. Kara SIZ - 30% BMP, biannual MDA",  x ="Year", y = "Microfilarial prevalence (%)") + geom_errorbar(data = KaraYes1, aes(ymin=CILower, ymax=CIUpper, x = Survey_year), width=.1, size=0.5)+ geom_point(data = KaraYes1,aes(x = Survey_year, y = PREV_cr, colour=KaraYes1$ENDE), size = 4) +scale_colour_manual(values=("yellow"))+   labs(colour=NULL)+geom_vline(xintercept=1977,colour="green") + geom_vline(xintercept=2007,colour="dark green")+geom_vline(xintercept=1991,colour="red")+geom_text(mapping=aes(x=1991, y=80, label="MDA start"), size=4, angle=90, vjust=-0.4, hjust=0)+geom_text(mapping=aes(x=2007, y=80, label="VC end"), size=4, angle=90, vjust=-0.4, hjust=0)+geom_text(mapping=aes(x=1977, y=80, label="VC start"), size=4, angle=90, vjust=-0.4, hjust=0)+geom_text(mapping=aes(x=2003, y=65, label="Biannual MDA start"), size=4, angle=90, vjust=-0.4, hjust=0)+geom_vline(xintercept=2003,colour="ORANGE") +theme(
  axis.title.x = element_text(size = 14),  # Increase x-axis title size
  axis.title.y = element_text(size = 14),  # Increase y-axis title size
  axis.text.x = element_text(size = 12),   # Increase x-axis text size
  axis.text.y = element_text(size = 12)    # Increase y-axis text size
)+scale_y_continuous(breaks = seq(0, 100, by = 10)) +  # Y-axis ticks every 10%
  scale_x_continuous(breaks = seq(1970, 2030, by = 5))  # X-axis ticks every 5 years

#Mesoendemicity under biannual MDA
list2 <- lapply(means.examplesKSOB, `[[`, 2)
list22 <- lapply(means.examplesKSRB, `[[`, 2)
list222 <- lapply(means.examplesKSPB, `[[`, 2)
a2 <- ggplot()+geom_line(aes(time,list2[[2]]*100), colour="lightblue", size=1)+geom_line(aes(time,list22[[2]]*100), colour="darkblue", size=1)+geom_line(aes(time,list222[[2]]*100), colour="lightblue", size=1)+ coord_cartesian(xlim = c(1972.75, 2027.28), ylim = c(4.52,95.5)) + theme_classic() +labs(title="b. Kara SIZ - 50% BMP, biannual MDA",  x ="Year", y = "Microfilarial prevalence (%)") + geom_errorbar(data = KaraYes2U, aes(ymin=CILower, ymax=CIUpper, x = Survey_year), width=.1, size=0.5)+ geom_point(data = KaraYes2U,aes(x = Survey_year, y = PREV_cr, colour=KaraYes2U$ENDE), size = 4) +scale_colour_manual(values="orange")+   labs(colour=NULL)+geom_vline(xintercept=1977,colour="green") + geom_vline(xintercept=2007,colour="dark green")+geom_vline(xintercept=1991,colour="red")+geom_text(mapping=aes(x=1991, y=80, label="MDA start"), size=4, angle=90, vjust=-0.4, hjust=0)+geom_text(mapping=aes(x=2007, y=80, label="VC end"), size=4, angle=90, vjust=-0.4, hjust=0)+geom_text(mapping=aes(x=1977, y=80, label="VC start"), size=4, angle=90, vjust=-0.4, hjust=0)+geom_text(mapping=aes(x=2003, y=65, label="Biannual MDA start"), size=4, angle=90, vjust=-0.4, hjust=0)+geom_vline(xintercept=2003,colour="ORANGE") + theme(
  axis.title.x = element_text(size = 14),  # Increase x-axis title size
  axis.title.y = element_text(size = 14),  # Increase y-axis title size
  axis.text.x = element_text(size = 12),   # Increase x-axis text size
  axis.text.y = element_text(size = 12)    # Increase y-axis text size
)+scale_y_continuous(breaks = seq(0, 100, by = 10)) +  # Y-axis ticks every 10%
  scale_x_continuous(breaks = seq(1970, 2030, by = 5))  # X-axis ticks every 5 years

#Hyperendemicity under biannual MDA
list3 <- lapply(means.examplesKSOB, `[[`, 3)
list33 <- lapply(means.examplesKSRB, `[[`, 3)
list333 <- lapply(means.examplesKSPB, `[[`, 3)
a3 <- ggplot()+geom_line(aes(time,list3[[2]]*100), colour="lightblue", size=1)+geom_line(aes(time,list33[[2]]*100), colour="darkblue", size=1)+geom_line(aes(time,list333[[2]]*100), colour="lightblue", size=1)+ coord_cartesian(xlim = c(1972.75, 2027.28), ylim = c(4.52,95.5)) + theme_classic() +labs(title="c. Kara SIZ - 70% BMP, biannual MDA",  x ="Year", y = "Microfilarial prevalence (%)")+ geom_errorbar(data = KaraYes3U, aes(ymin=CILower, ymax=CIUpper, x = Survey_year), width=.1, size=0.5) + geom_point(data = KaraYes3U,aes(x = Survey_year, y = PREV_cr, colour=KaraYes3U$ENDE), size = 4) +scale_colour_manual(values="dark red")+   labs(colour=NULL)+geom_vline(xintercept=1977,colour="green") + geom_vline(xintercept=2007,colour="dark green")+geom_vline(xintercept=1991,colour="red")+geom_text(mapping=aes(x=1991, y=80, label="MDA start"), size=4, angle=90, vjust=-0.4, hjust=0)+geom_text(mapping=aes(x=2007, y=80, label="VC end"), size=4, angle=90, vjust=-0.4, hjust=0)+geom_text(mapping=aes(x=1977, y=80, label="VC start"), size=4, angle=90, vjust=1.1, hjust=0)+geom_text(mapping=aes(x=2003, y=65, label="Biannual MDA start"), size=4, angle=90, vjust=-0.4, hjust=0)+geom_vline(xintercept=2003,colour="ORANGE") + theme(
  axis.title.x = element_text(size = 14),  # Increase x-axis title size
  axis.title.y = element_text(size = 14),  # Increase y-axis title size
  axis.text.x = element_text(size = 12),   # Increase x-axis text size
  axis.text.y = element_text(size = 12)    # Increase y-axis text size
)+scale_y_continuous(breaks = seq(0, 100, by = 10)) +  # Y-axis ticks every 10%
  scale_x_continuous(breaks = seq(1970, 2030, by = 5))  # X-axis ticks every 5 years

#Holoendemicity under biannual MDA
list4 <- lapply(means.examplesKSOB, `[[`, 4)
list44 <- lapply(means.examplesKSRB, `[[`, 4)
list444 <- lapply(means.examplesKSPB, `[[`, 4)
a4 <- ggplot()+geom_line(aes(time,list4[[2]]*100), colour="lightblue", size=1)+geom_line(aes(time,list44[[2]]*100), colour="darkblue", size=1)+geom_line(aes(time,list444[[2]]*100), colour="lightblue", size=1)+ coord_cartesian(xlim = c(1972.75, 2027.28), ylim = c(4.52,95.5)) + theme_classic() +labs(title="d. Kara SIZ - 90% BMP, biannual MDA",  x ="Year", y = "Microfilarial prevalence (%)") + geom_errorbar(data = KaraYes4, aes(ymin=CILower, ymax=CIUpper, x = Survey_year), width=.1, size=0.5)+ geom_point(data = KaraYes4,aes(x = Survey_year, y = PREV_cr, colour=KaraYes4$ENDE), size = 4)+scale_colour_manual(values="dark red")+   labs(colour=NULL)+geom_vline(xintercept=1977,colour="green") + geom_vline(xintercept=2007,colour="dark green")+geom_vline(xintercept=1991,colour="red")+geom_text(mapping=aes(x=1991, y=50, label="MDA start"), size=4, angle=90, vjust=-0.4, hjust=0)+geom_text(mapping=aes(x=2007, y=80, label="VC end"), size=4, angle=90, vjust=-0.4, hjust=0)+geom_text(mapping=aes(x=1977, y=50, label="VC start"), size=4, angle=90, vjust=1.1, hjust=0)+geom_text(mapping=aes(x=2003, y=65, label="Biannual MDA start"), size=4, angle=90, vjust=-0.4, hjust=0)+geom_vline(xintercept=2003,colour="ORANGE") + theme(
  axis.title.x = element_text(size = 14),  # Increase x-axis title size
  axis.title.y = element_text(size = 14),  # Increase y-axis title size
  axis.text.x = element_text(size = 12),   # Increase x-axis text size
  axis.text.y = element_text(size = 12)    # Increase y-axis text size
)+scale_y_continuous(breaks = seq(0, 100, by = 10)) +  # Y-axis ticks every 10%
  scale_x_continuous(breaks = seq(1970, 2030, by = 5))  # X-axis ticks every 5 years

combined_plot <- (a1+a2)/(a3+a4)
ggsave(filename = "Figure 3.pdf", plot = combined_plot, width = 15, height = 21*2/3, units = "in", dpi = 300)


#5. Centrale non-SIZ 
#hypoendemic under annual MDA
list1 <- lapply(means.examplesKNSO, `[[`, 1) #A list with 4 numeric columns (plus the fifth). We just need to select the second one, that corresponds to the prevalence
list11 <- lapply(means.examplesKNSR, `[[`, 1) #A list with 4 numeric columns (plus the fifth). We just need to select the second one, that corresponds to the prevalence
list111 <- lapply(means.examplesKNSP, `[[`, 1) #A list with 4 numeric columns (plus the fifth). We just need to select the second one, that corresponds to the prevalence
a1 <- ggplot()+geom_line(aes(time,list1[[2]]*100), colour="lightblue", size=1)+geom_line(aes(time,list11[[2]]*100), colour="darkblue", size=1)+geom_line(aes(time,list111[[2]]*100), colour="lightblue", size=1)+ coord_cartesian(xlim = c(1972.75, 2027.28), ylim = c(4.52,95.5)) + theme_classic() +labs(title="b. Centrale non-SIZ - 30% BMP",  x ="Year", y = "Microfilarial prevalence (%)") + geom_errorbar(data = CentraleNo1A, aes(ymin=CILower, ymax=CIUpper, x = Survey_year), width=.1, size=0.5)+ geom_point(data = CentraleNo1A,aes(x = Survey_year, y = PREV_cr, colour=CentraleNo1A$ENDE), size = 4)+scale_colour_manual(values="yellow")+   labs(colour=NULL) + geom_vline(xintercept=1989,colour="green") + geom_vline(xintercept=2002,colour="dark green")+geom_vline(xintercept=1991,colour="red")+geom_text(mapping=aes(x=1991, y=80, label="MDA start"), size=4, angle=90, vjust=-0.4, hjust=0)+geom_text(mapping=aes(x=2002, y=80, label="VC end"), size=4, angle=90, vjust=-0.4, hjust=0)+geom_text(mapping=aes(x=1989, y=80, label="VC start"), size=4, angle=90, vjust=-0.4, hjust=0)+ theme(
  axis.title.x = element_text(size = 14),  # Increase x-axis title size
  axis.title.y = element_text(size = 14),  # Increase y-axis title size
  axis.text.x = element_text(size = 12),   # Increase x-axis text size
  axis.text.y = element_text(size = 12)    # Increase y-axis text size
)+scale_y_continuous(breaks = seq(0, 100, by = 10)) +  # Y-axis ticks every 10%
  scale_x_continuous(breaks = seq(1970, 2030, by = 5))  # X-axis ticks every 5 years

#Mesoendemic under annual MDA
list2 <- lapply(means.examplesKNSO, `[[`, 2)
list22 <- lapply(means.examplesKNSR, `[[`, 2)
list222 <- lapply(means.examplesKNSP, `[[`, 2)
a2 <- ggplot()+geom_line(aes(time,list2[[2]]*100), colour="lightblue", size=1)+geom_line(aes(time,list22[[2]]*100), colour="darkblue", size=1)+geom_line(aes(time,list222[[2]]*100), colour="lightblue", size=1)+ coord_cartesian(xlim = c(1972.75, 2027.28), ylim = c(4.52,95.5)) + theme_classic() +labs(title="c. Centrale non-SIZ - 50% BMP",  x ="Year", y = "Microfilarial prevalence (%)") + geom_errorbar(data = CentraleNo2A, aes(ymin=CILower, ymax=CIUpper, x = Survey_year), width=.1, size=0.5)+ geom_point(data = CentraleNo2A,aes(x = Survey_year, y = PREV_cr, colour=CentraleNo2A$ENDE), size = 4) +scale_colour_manual(values="orange")+   labs(colour=NULL) + geom_vline(xintercept=1989,colour="green") + geom_vline(xintercept=2002,colour="dark green")+geom_vline(xintercept=1991,colour="red")+geom_text(mapping=aes(x=1991, y=80, label="MDA start"), size=4, angle=90, vjust=-0.4, hjust=0)+geom_text(mapping=aes(x=2002, y=80, label="VC end"), size=4, angle=90, vjust=-0.4, hjust=0)+geom_text(mapping=aes(x=1989, y=80, label="VC start"), size=4, angle=90, vjust=-0.4, hjust=0)+ theme(
  axis.title.x = element_text(size = 14),  # Increase x-axis title size
  axis.title.y = element_text(size = 14),  # Increase y-axis title size
  axis.text.x = element_text(size = 12),   # Increase x-axis text size
  axis.text.y = element_text(size = 12)    # Increase y-axis text size
)+scale_y_continuous(breaks = seq(0, 100, by = 10)) +  # Y-axis ticks every 10%
  scale_x_continuous(breaks = seq(1970, 2030, by = 5))  # X-axis ticks every 5 years

#Hyperendemicity under annual MDA
list3 <- lapply(means.examplesKNSO, `[[`, 3)
list33 <- lapply(means.examplesKNSR, `[[`, 3)
list333 <- lapply(means.examplesKNSP, `[[`, 3)
a3 <- ggplot()+geom_line(aes(time,list3[[2]]*100), colour="lightblue", size=1)+geom_line(aes(time,list33[[2]]*100), colour="darkblue", size=1)+geom_line(aes(time,list333[[2]]*100), colour="lightblue", size=1)+ coord_cartesian(xlim = c(1972.75, 2027.28), ylim = c(4.52,95.5)) + theme_classic() +labs(title="d. Centrale non-SIZ - 70% BMP",  x ="Year", y = "Microfilarial prevalence (%)") + geom_errorbar(data = CentraleNo3A, aes(ymin=CILower, ymax=CIUpper, x = Survey_year), width=.1, size=0.5)+ geom_point(data = CentraleNo3A,aes(x = Survey_year, y = PREV_cr, colour=CentraleNo3A$ENDE), size = 4) +scale_colour_manual(values="dark red")+   labs(colour=NULL) + geom_vline(xintercept=1989,colour="green") + geom_vline(xintercept=2002,colour="dark green")+geom_vline(xintercept=1991,colour="red")+geom_text(mapping=aes(x=1991, y=80, label="MDA start"), size=4, angle=90, vjust=1.1, hjust=0)+geom_text(mapping=aes(x=2002, y=80, label="VC end"), size=4, angle=90, vjust=-0.4, hjust=0)+geom_text(mapping=aes(x=1989, y=80, label="VC start"), size=4, angle=90, vjust=-0.4, hjust=0)+ theme(
  axis.title.x = element_text(size = 14),  # Increase x-axis title size
  axis.title.y = element_text(size = 14),  # Increase y-axis title size
  axis.text.x = element_text(size = 12),   # Increase x-axis text size
  axis.text.y = element_text(size = 12)    # Increase y-axis text size
)+scale_y_continuous(breaks = seq(0, 100, by = 10)) +  # Y-axis ticks every 10%
  scale_x_continuous(breaks = seq(1970, 2030, by = 5))  # X-axis ticks every 5 years

#6. Centrale SIZ 
#Hyperendemic under biannual MDA
list4 <- lapply(means.examplesKSOB, `[[`, 3)
list44 <- lapply(means.examplesKSRB, `[[`, 3)
list444 <- lapply(means.examplesKSPB, `[[`, 3)
a4 <- ggplot()+geom_line(aes(time,list4[[2]]*100), colour="lightblue", size=1)+geom_line(aes(time,list44[[2]]*100), colour="darkblue", size=1)+geom_line(aes(time,list444[[2]]*100), colour="lightblue", size=1)+ coord_cartesian(xlim = c(1972.75, 2027.28), ylim = c(4.52,95.5)) + theme_classic() +labs(title="a. Centrale SIZ - 70% BMP, biannual MDA",  x ="Year", y = "Microfilarial prevalence (%)") + geom_errorbar(data = CentraleYes3B, aes(ymin=CILower, ymax=CIUpper, x = Survey_year), width=.1, size=0.5)+ geom_point(data = CentraleYes3B,aes(x = Survey_year, y = PREV_cr, colour=CentraleYes3B$ENDE), size = 4) +scale_colour_manual(values=c("dark red"))+   labs(colour=NULL)+geom_vline(xintercept=1977,colour="green") + geom_vline(xintercept=2007,colour="dark green")+geom_vline(xintercept=1991,colour="red")+geom_vline(xintercept=2003,colour="ORANGE")+geom_text(mapping=aes(x=1991, y=80, label="MDA start"), size=4, angle=90, vjust=-0.4, hjust=0)+geom_text(mapping=aes(x=2007, y=80, label="VC end"), size=4, angle=90, vjust=-0.4, hjust=0)+geom_text(mapping=aes(x=1977, y=80, label="VC start"), size=4, angle=90, vjust=-0.4, hjust=0)+geom_text(mapping=aes(x=2003, y=65, label="Biannual MDA start"), size=4, angle=90, vjust=-0.4, hjust=0)+ theme(
  axis.title.x = element_text(size = 14),  # Increase x-axis title size
  axis.title.y = element_text(size = 14),  # Increase y-axis title size
  axis.text.x = element_text(size = 12),   # Increase x-axis text size
  axis.text.y = element_text(size = 12)    # Increase y-axis text size
)+scale_y_continuous(breaks = seq(0, 100, by = 10)) +  # Y-axis ticks every 10%
  scale_x_continuous(breaks = seq(1970, 2030, by = 5))  # X-axis ticks every 5 years

combined_plot <- (a4 + a1) /  (a2+a3)
ggsave(filename = "Figure 4.pdf", plot = combined_plot, width = 15, height = 21*2/3, units = "in", dpi = 300)




#7. Plateaux 
#Hypoendemic under annual MDA
list1 <- lapply(means.examplesPO, `[[`, 1) #A list with 4 numeric columns (plus the fifth). We just need to select the second one, that corresponds to the prevalence
list11 <- lapply(means.examplesPR, `[[`, 1) #A list with 4 numeric columns (plus the fifth). We just need to select the second one, that corresponds to the prevalence
list111 <- lapply(means.examplesPP, `[[`, 1) #A list with 4 numeric columns (plus the fifth). We just need to select the second one, that corresponds to the prevalence
a1 <- ggplot()+geom_line(aes(time,list1[[2]]*100), colour="lightblue", size=1)+geom_line(aes(time,list11[[2]]*100), colour="darkblue", size=1)+geom_line(aes(time,list111[[2]]*100), colour="lightblue", size=1)+ coord_cartesian(xlim = c(1972.75, 2027.28), ylim = c(4.52,95.5)) + theme_classic() +labs(title="a. Plateaux - 30% BMP",  x ="Year", y = "Microfilarial prevalence (%)") + geom_errorbar(data = PlateauxNo1A, aes(ymin=CILower, ymax=CIUpper, x = Survey_year), width=.1, size=0.5)+ geom_point(data = PlateauxNo1A,aes(x = Survey_year, y = PREV_cr, colour=PlateauxNo1A$ENDE), size=4) +geom_vline(xintercept=1989,colour="green") + geom_vline(xintercept=2002,colour="dark green")+geom_vline(xintercept=1991,colour="red")+scale_colour_manual(values=c("yellow"))+labs(colour=NULL)+geom_text(mapping=aes(x=1991, y=80, label="MDA start"), size=4, angle=90, vjust=--1.2, hjust=0)+geom_text(mapping=aes(x=2002, y=80, label="VC end"), size=4, angle=90, vjust=-0.4, hjust=0)+geom_text(mapping=aes(x=1989, y=80, label="VC start"), size=4, angle=90, vjust=-0.4, hjust=0)+ theme(
  axis.title.x = element_text(size = 14),  # Increase x-axis title size
  axis.title.y = element_text(size = 14),  # Increase y-axis title size
  axis.text.x = element_text(size = 12),   # Increase x-axis text size
  axis.text.y = element_text(size = 12)    # Increase y-axis text size
)+scale_y_continuous(breaks = seq(0, 100, by = 10)) +  # Y-axis ticks every 10%
  scale_x_continuous(breaks = seq(1970, 2030, by = 5))  # X-axis ticks every 5 years

#hypoendemic under biannual MDA
list5 <- lapply(means.examplesPOB, `[[`, 1) #A list with 4 numeric columns (plus the fifth). We just need to select the second one, that corresponds to the prevalence
list55 <- lapply(means.examplesPRB, `[[`, 1) #A list with 4 numeric columns (plus the fifth). We just need to select the second one, that corresponds to the prevalence
list555 <- lapply(means.examplesPPB, `[[`, 1) #A list with 4 numeric columns (plus the fifth). We just need to select the second one, that corresponds to the prevalence
a5 <- ggplot()+geom_line(aes(time,list5[[2]]*100), colour="lightblue", size=1)+geom_line(aes(time,list55[[2]]*100), colour="darkblue", size=1)+geom_line(aes(time,list555[[2]]*100), colour="lightblue", size=1)+ coord_cartesian(xlim = c(1972.75, 2027.28), ylim = c(4.52,95.5)) + theme_classic() +labs(title="b. Plateaux - 30% BMP, biannual MDA",  x ="Year", y = "Microfilarial prevalence (%)") + geom_errorbar(data = PlateauxNo1B, aes(ymin=CILower, ymax=CIUpper, x = Survey_year), width=.1, size=0.5)+ geom_point(data = PlateauxNo1B,aes(x = Survey_year, y = PREV_cr, colour=PlateauxNo1B$ENDE), size=4) +geom_vline(xintercept=1989,colour="green") + geom_vline(xintercept=2002,colour="dark green")+geom_vline(xintercept=1991,colour="red")+scale_colour_manual(values=c("yellow"))+labs(colour=NULL)+geom_text(mapping=aes(x=1991, y=80, label="MDA start"), size=4, angle=90, vjust=--1.2, hjust=0)+geom_text(mapping=aes(x=2002, y=80, label="VC end"), size=4, angle=90, vjust=-0.4, hjust=0)+geom_text(mapping=aes(x=1989, y=80, label="VC start"), size=4, angle=90, vjust=-0.4, hjust=0)+ geom_text(mapping=aes(x=2014, y=65, label="Biannual MDA start"), size=4, angle=90, vjust=--1.2, hjust=0)+geom_vline(xintercept=2014,colour="ORANGE") +theme(
  axis.title.x = element_text(size = 14),  # Increase x-axis title size
  axis.title.y = element_text(size = 14),  # Increase y-axis title size
  axis.text.x = element_text(size = 12),   # Increase x-axis text size
  axis.text.y = element_text(size = 12)    # Increase y-axis text size
)+scale_y_continuous(breaks = seq(0, 100, by = 10)) +  # Y-axis ticks every 10%
  scale_x_continuous(breaks = seq(1970, 2030, by = 5))  # X-axis ticks every 5 years

#mesoendemic under annual MDA
list2 <- lapply(means.examplesPO, `[[`, 2)
list22 <- lapply(means.examplesPR, `[[`, 2)
list222 <- lapply(means.examplesPP, `[[`, 2)
a2 <- ggplot()+geom_line(aes(time,list2[[2]]*100), colour="lightblue", size=1)+geom_line(aes(time,list22[[2]]*100), colour="darkblue", size=1)+geom_line(aes(time,list222[[2]]*100), colour="lightblue", size=1)+ coord_cartesian(xlim = c(1972.75, 2027.28), ylim = c(4.52,95.5)) + theme_classic() +labs(title="c. Plateaux - 50% BMP",  x ="Year", y = "Microfilarial prevalence (%)") + geom_errorbar(data = PlateauxNo2A, aes(ymin=CILower, ymax=CIUpper, x = Survey_year), width=.1, size=0.5)+ geom_point(data = PlateauxNo2A,aes(x = Survey_year, y = PREV_cr, colour=PlateauxNo2A$ENDE), size=4) +geom_vline(xintercept=1989,colour="green") + geom_vline(xintercept=2002,colour="dark green")+geom_vline(xintercept=1991,colour="red")+scale_colour_manual(values=c("orange","dark red"))+   labs(colour=NULL)+geom_text(mapping=aes(x=1991, y=80, label="MDA start"), size=4, angle=90, vjust=+1.2, hjust=0)+geom_text(mapping=aes(x=2002, y=80, label="VC end"), size=4, angle=90, vjust=-0.4, hjust=0)+geom_text(mapping=aes(x=1989, y=80, label="VC start"), size=4, angle=90, vjust=-0.4, hjust=0)+ theme(
  axis.title.x = element_text(size = 14),  # Increase x-axis title size
  axis.title.y = element_text(size = 14),  # Increase y-axis title size
  axis.text.x = element_text(size = 12),   # Increase x-axis text size
  axis.text.y = element_text(size = 12)    # Increase y-axis text size
)+scale_y_continuous(breaks = seq(0, 100, by = 10)) +  # Y-axis ticks every 10%
  scale_x_continuous(breaks = seq(1970, 2030, by = 5))  # X-axis ticks every 5 years

#meso under biannual MDA
PlateauxNo2B <- subset(PlateauxNo2B, !(VILNOM == "ADJABOULOUKOUKOPE" & CECI == 9999  &  POS_BASE == 21)) #Remove an unnusual value, same as the survey 1997.715
PlateauxNo2B <- subset(PlateauxNo2B, !(VILNOM == "KPODJI" )) #KPODJI is hyperendemic
list6 <- lapply(means.examplesPOB, `[[`, 2)
list66 <- lapply(means.examplesPRB, `[[`, 2)
list666 <- lapply(means.examplesPPB, `[[`, 2)
a6 <- ggplot()+geom_line(aes(time,list6[[2]]*100), colour="lightblue", size=1)+geom_line(aes(time,list66[[2]]*100), colour="darkblue", size=1)+geom_line(aes(time,list666[[2]]*100), colour="lightblue", size=1)+ coord_cartesian(xlim = c(1972.75, 2027.28), ylim = c(4.52,95.5)) + theme_classic() +labs(title="d. Plateaux - 50% BMP, biannual MDA",  x ="Year", y = "Microfilarial prevalence (%)") + geom_errorbar(data = PlateauxNo2B, aes(ymin=CILower, ymax=CIUpper, x = Survey_year), width=.1, size=0.5)+ geom_point(data = PlateauxNo2B,aes(x = Survey_year, y = PREV_cr, colour=PlateauxNo2B$ENDE), size=4) +geom_vline(xintercept=1989,colour="green") + geom_vline(xintercept=2002,colour="dark green")+geom_vline(xintercept=1991,colour="red")+scale_colour_manual(values=c("orange","dark red"))+   labs(colour=NULL)+geom_text(mapping=aes(x=1991, y=80, label="MDA start"), size=4, angle=90, vjust=+1.2, hjust=0)+geom_text(mapping=aes(x=2002, y=80, label="VC end"), size=4, angle=90, vjust=-0.4, hjust=0)+geom_text(mapping=aes(x=1989, y=80, label="VC start"), size=4, angle=90, vjust=-0.4, hjust=0)+ geom_text(mapping=aes(x=2014, y=65, label="Biannual MDA start"), size=4, angle=90, vjust=--1.2, hjust=0)+geom_vline(xintercept=2014,colour="ORANGE") + theme(
  axis.title.x = element_text(size = 14),  # Increase x-axis title size
  axis.title.y = element_text(size = 14),  # Increase y-axis title size
  axis.text.x = element_text(size = 12),   # Increase x-axis text size
  axis.text.y = element_text(size = 12)    # Increase y-axis text size
)+scale_y_continuous(breaks = seq(0, 100, by = 10)) +  # Y-axis ticks every 10%
  scale_x_continuous(breaks = seq(1970, 2030, by = 5))  # X-axis ticks every 5 years


#hyperendemic under annual MDA
list3 <- lapply(means.examplesPO, `[[`, 3)
list33 <- lapply(means.examplesPR, `[[`, 3)
list333 <- lapply(means.examplesPP, `[[`, 3)
a3 <- ggplot()+geom_line(aes(time,list3[[2]]*100), colour="lightblue", size=1)+geom_line(aes(time,list33[[2]]*100), colour="darkblue", size=1)+geom_line(aes(time,list333[[2]]*100), colour="lightblue", size=1)+ coord_cartesian(xlim = c(1972.75, 2027.28), ylim = c(4.52,95.5)) + theme_classic() +labs(title="e. Plateaux - 70% BMP",  x ="Year", y = "Microfilarial prevalence (%)") + geom_errorbar(data = PlateauxNo3A, aes(ymin=CILower, ymax=CIUpper, x = Survey_year), width=.1, size=0.5)+ geom_point(data = PlateauxNo3A,aes(x = Survey_year, y = PREV_cr, colour=PlateauxNo3A$ENDE), size=4)+ geom_vline(xintercept=1989,colour="green") + geom_vline(xintercept=2002,colour="dark green")+geom_vline(xintercept=1991,colour="red")+scale_colour_manual(values="dark red")+   labs(colour=NULL)+geom_text(mapping=aes(x=1991, y=80, label="MDA start"), size=4, angle=90, vjust=-0.4, hjust=0)+geom_text(mapping=aes(x=2002, y=80, label="VC end"), size=4, angle=90, vjust=-0.4, hjust=0)+geom_text(mapping=aes(x=1989, y=80, label="VC start"), size=4, angle=90, vjust=-0.4, hjust=0)+ theme(
  axis.title.x = element_text(size = 14),  # Increase x-axis title size
  axis.title.y = element_text(size = 14),  # Increase y-axis title size
  axis.text.x = element_text(size = 12),   # Increase x-axis text size
  axis.text.y = element_text(size = 12)    # Increase y-axis text size
)+scale_y_continuous(breaks = seq(0, 100, by = 10)) +  # Y-axis ticks every 10%
  scale_x_continuous(breaks = seq(1970, 2030, by = 5))  # X-axis ticks every 5 years

#hyper under biannual MDA
list7 <- lapply(means.examplesPOB, `[[`, 3)
list77 <- lapply(means.examplesPRB, `[[`, 3)
list777 <- lapply(means.examplesPPB, `[[`, 3)
a7 <- ggplot()+geom_line(aes(time,list7[[2]]*100), colour="lightblue", size=1)+geom_line(aes(time,list77[[2]]*100), colour="darkblue", size=1)+geom_line(aes(time,list777[[2]]*100), colour="lightblue", size=1)+ coord_cartesian(xlim = c(1972.75, 2027.28), ylim = c(4.52,95.5)) + theme_classic() +labs(title="f. Plateaux - 70% BMP, biannual MDA",  x ="Year", y = "Microfilarial prevalence (%)") + geom_errorbar(data = PlateauxNo3B, aes(ymin=CILower, ymax=CIUpper, x = Survey_year), width=.1, size=0.5)+ geom_point(data = PlateauxNo3B,aes(x = Survey_year, y = PREV_cr, colour=PlateauxNo3B$ENDE), size=4)+geom_vline(xintercept=1989,colour="green") + geom_vline(xintercept=2002,colour="dark green")+geom_vline(xintercept=1991,colour="red")+scale_colour_manual(values="dark red")+   labs(colour=NULL)+geom_text(mapping=aes(x=1991, y=80, label="MDA start"), size=4, angle=90, vjust=1.2, hjust=0)+geom_text(mapping=aes(x=2002, y=80, label="VC end"), size=4, angle=90, vjust=-0.4, hjust=0)+geom_text(mapping=aes(x=1989, y=80, label="VC start"), size=4, angle=90, vjust=-0.4, hjust=0)+geom_text(mapping=aes(x=2014, y=65, label="Biannual MDA start"), size=4, angle=90, vjust=-0.4, hjust=0)+geom_vline(xintercept=2014,colour="ORANGE")+ theme(
  axis.title.x = element_text(size = 14),  # Increase x-axis title size
  axis.title.y = element_text(size = 14),  # Increase y-axis title size
  axis.text.x = element_text(size = 12),   # Increase x-axis text size
  axis.text.y = element_text(size = 12)    # Increase y-axis text size
)+scale_y_continuous(breaks = seq(0, 100, by = 10)) +  # Y-axis ticks every 10%
  scale_x_continuous(breaks = seq(1970, 2030, by = 5))  # X-axis ticks every 5 years

combined_plot <- (a1 + a5) /  (a2+a6) /(a3 + a7)
ggsave(filename = "Figure 5.pdf", plot = combined_plot, width = 15, height = 21, units = "in", dpi = 300)



#8. Maritime - until 2030
#hypoendemic under annual MDA
list1 <- lapply(means.examplesMO, `[[`, 1) #A list with 4 numeric columns (plus the fifth). We just need to select the second one, that corresponds to the prevalence
list11 <- lapply(means.examplesMR, `[[`, 1) #A list with 4 numeric columns (plus the fifth). We just need to select the second one, that corresponds to the prevalence
list111 <- lapply(means.examplesMP, `[[`, 1) #A list with 4 numeric columns (plus the fifth). We just need to select the second one, that corresponds to the prevalence
a1 <- ggplot()+geom_line(aes(time,list1[[2]]*100), colour="light blue", size=1)+geom_line(aes(time,list11[[2]]*100), colour="darkblue", size=1)+geom_line(aes(time,list111[[2]]*100), colour="lightblue", size=1)+ coord_cartesian(xlim = c(1972.75, 2027.28), ylim = c(4.52,95.5)) + theme_classic() +labs(title="a. Maritime - 30% BMP",  x ="Year", y = "Microfilarial prevalence (%)") + geom_errorbar(data = MaritimeNo1, aes(ymin=CILower, ymax=CIUpper, x = Survey_year), width=.1, size=0.5)+ geom_point(data = MaritimeNo1,aes(x = Survey_year, y = PREV_cr, colour=MaritimeNo1$ENDE), size=4)+geom_vline(xintercept=1988,colour="green") + geom_vline(xintercept=2002,colour="dark green")+geom_vline(xintercept=1991,colour="red")+scale_colour_manual(values="yellow")+   labs(colour=NULL)+geom_text(mapping=aes(x=1991, y=80, label="MDA start"), size=4, angle=90, vjust=-0.4, hjust=0)+geom_text(mapping=aes(x=2002, y=80, label="VC end"), size=4, angle=90, vjust=-0.4, hjust=0)+geom_text(mapping=aes(x=1988, y=80, label="VC start"), size=4, angle=90, vjust=-0.4, hjust=0)+ theme(
  axis.title.x = element_text(size = 14),  # Increase x-axis title size
  axis.title.y = element_text(size = 14),  # Increase y-axis title size
  axis.text.x = element_text(size = 12),   # Increase x-axis text size
  axis.text.y = element_text(size = 12)    # Increase y-axis text size
)+scale_y_continuous(breaks = seq(0, 100, by = 10)) +  # Y-axis ticks every 10%
  scale_x_continuous(breaks = seq(1970, 2030, by = 5))  # X-axis ticks every 5 years

#hyper under annual MDA
list3 <- lapply(means.examplesMO, `[[`, 3)
list33 <- lapply(means.examplesMR, `[[`, 3)
list333 <- lapply(means.examplesMP, `[[`, 3)
a3 <- ggplot()+geom_line(aes(time,list3[[2]]*100), colour="lightblue", size=1)+geom_line(aes(time,list33[[2]]*100), colour="darkblue", size=1)+geom_line(aes(time,list333[[2]]*100), colour="lightblue", size=1)+ coord_cartesian(xlim = c(1972.75, 2027.28), ylim = c(4.52,95.5)) + theme_classic() +labs(title="b. Maritime - 70% BMP",  x ="Year", y = "Microfilarial prevalence (%)")+ geom_errorbar(data = MaritimeNo3, aes(ymin=CILower, ymax=CIUpper, x = Survey_year), width=.1, size=0.5) + geom_point(data = MaritimeNo3,aes(x = Survey_year, y = PREV_cr, colour=MaritimeNo3$ENDE), size=4) +geom_vline(xintercept=1988,colour="green") + geom_vline(xintercept=2002,colour="dark green")+geom_vline(xintercept=1991,colour="red")+scale_colour_manual(values="dark red")+   labs(colour=NULL)+geom_text(mapping=aes(x=1991, y=80, label="MDA start"), size=4, angle=90, vjust=-0.4, hjust=0)+geom_text(mapping=aes(x=2002, y=80, label="VC end"), size=4, angle=90, vjust=-0.4, hjust=0)+geom_text(mapping=aes(x=1988, y=80, label="VC start"), size=4, angle=90, vjust=-0.4, hjust=0)+ theme(
  axis.title.x = element_text(size = 14),  # Increase x-axis title size
  axis.title.y = element_text(size = 14),  # Increase y-axis title size
  axis.text.x = element_text(size = 12),   # Increase x-axis text size
  axis.text.y = element_text(size = 12)    # Increase y-axis text size
)+scale_y_continuous(breaks = seq(0, 100, by = 10)) +  # Y-axis ticks every 10%
  scale_x_continuous(breaks = seq(1970, 2030, by = 5))  # X-axis ticks every 5 years


#a1+a3
combined_plot <- (a1 + a3)
ggsave(filename = "Figure 6.pdf", plot = combined_plot, width = 15, height = 21/3, units = "in", dpi = 300)


