#Figures for the villages with no recorded baseline of Supplementary Material
#Baseline Prevalence = BMP

#1839 means the vector control started in 1976 (true for the surveyed vilalges with pre-control info of Savanes non-SIZ)
library(ggpubr)
library(patchwork) #to be able to: a2 + a6
#1. Savanes non-SIZ 
#hypoendemic and mesoendemic under annual MDA (Tone and Kpendjal prefecture)
list1 <- lapply(means.examplesSNSOB, `[[`, 1)
list11 <- lapply(means.examplesSNSRB, `[[`, 1)
list111 <- lapply(means.examplesSNSPB, `[[`, 1)
a1 <- ggplot()+geom_line(aes(time,list1[[2]]*100), colour="lightblue", size=1)+geom_line(aes(time,list11[[2]]*100), colour="darkblue", size=1)+geom_line(aes(time,list111[[2]]*100), colour="lightblue", size=1)+ coord_cartesian(xlim = c(1972.75, 2027.28), ylim = c(4.52,95.5)) + theme_classic() +labs(title="A: Savanes non-SIZ - 30% BMP",  x ="Year", y = "Microfilarial prevalence (%)") + geom_errorbar(data = SavanesNoEND1, aes(ymin=CILower, ymax=CIUpper, x = Survey_year-1), width=.1, size=0.5)+ geom_point(data = SavanesNoEND1,aes(x = Survey_year-1, y = PREV_cr,colour=SavanesNoEND1$ENDE), size = 4)+scale_colour_manual(values=c("yellow","dark red"))+   labs(colour=NULL)+geom_vline(xintercept=1977,colour="green") + geom_vline(xintercept=1994,colour="dark green")+geom_vline(xintercept=1991,colour="red")+geom_text(mapping=aes(x=1991, y=80, label="MDA start"), size=4, angle=90, vjust=-0.4, hjust=0)+geom_text(mapping=aes(x=1994, y=80, label="VC end"), size=4, angle=90, vjust=-0.4, hjust=0)+geom_text(mapping=aes(x=1977, y=80, label="VC start"), size=4, angle=90, vjust=-0.4, hjust=0)+theme(
  axis.title.x = element_text(size = 14),  # Increase x-axis title size
  axis.title.y = element_text(size = 14),  # Increase y-axis title size
  axis.text.x = element_text(size = 12),   # Increase x-axis text size
  axis.text.y = element_text(size = 12)    # Increase y-axis text size
)+scale_y_continuous(breaks = seq(0, 100, by = 10)) +  # Y-axis ticks every 10%
  scale_x_continuous(breaks = seq(1970, 2030, by = 5))  # X-axis ticks every 5 years
time <- (seq(1, length(list1[[2]]), by=1)/(366)) + 1840 


#2. Savanes SIZ 
#under biannual MDA with and without 100% vector control (VC) efficacy for hypo and meso that cannot be differentiated
list3 <- lapply(means.examplesSSO100, `[[`, 1)
list33 <- lapply(means.examplesSSR100, `[[`, 1)
list333 <- lapply(means.examplesSSP100, `[[`, 1)
list3333 <- lapply(means.examplesSSO100, `[[`, 2)
list33333 <- lapply(means.examplesSSR100, `[[`, 2)
list333333 <- lapply(means.examplesSSP100, `[[`, 2)
list9 <- lapply(means.examplesSSOB, `[[`, 1)
list99 <- lapply(means.examplesSSRB, `[[`, 1)
list999 <- lapply(means.examplesSSPB, `[[`, 1)
list9999 <- lapply(means.examplesSSOB, `[[`, 2)
list99999 <- lapply(means.examplesSSRB, `[[`, 2)
list999999 <- lapply(means.examplesSSPB, `[[`, 2)
a3 <- ggplot()+geom_line(aes(time,list3[[2]]*100), colour="pink", size=1)+geom_line(aes(time,list33[[2]]*100), colour="red", size=1)+geom_line(aes(time,list333[[2]]*100), colour="pink", size=1)+geom_line(aes(time,list9[[2]]*100), colour="lightblue", size=1)+geom_line(aes(time,list99[[2]]*100), colour="darkblue", size=1)+geom_line(aes(time,list999[[2]]*100), colour="lightblue", size=1)+ coord_cartesian(xlim = c(1972.75, 2027.28), ylim = c(4.52,95.5)) + theme_classic() +labs(title="A: Savanes SIZ - 30% BMP with biannual MDA, including 100% VC efficacy" , x ="Year", y = "Microfilarial prevalence (%)") + geom_errorbar(data = SavanesYesEND1, aes(ymin=CILower, ymax=CIUpper, x = Survey_year), width=.1, size=0.5)+ geom_point(data = SavanesYesEND1,aes(x = Survey_year, y = PREV_cr,colour=SavanesYesEND1$ENDE), size = 4)+scale_colour_manual(values=c("yellow"))+   labs(colour=NULL)+geom_vline(xintercept=1977,colour="green") + geom_vline(xintercept=1993,colour="dark green")+geom_vline(xintercept=1991,colour="red")+geom_vline(xintercept=2003,colour="ORANGE")+geom_text(mapping=aes(x=1991, y=80, label="MDA start"), size=4, angle=90, vjust=-0.4, hjust=0)+geom_text(mapping=aes(x=1993, y=80, label="VC end"), size=4, angle=90, vjust=-0.4, hjust=0)+geom_text(mapping=aes(x=1977, y=80, label="VC start"), size=4, angle=90, vjust=-0.4, hjust=0)+geom_text(mapping=aes(x=2003, y=65, label="Biannual MDA start"), size=4, angle=90, vjust=-0.4, hjust=0)+ theme(
  axis.title.x = element_text(size = 14),  # Increase x-axis title size
  axis.title.y = element_text(size = 14),  # Increase y-axis title size
  axis.text.x = element_text(size = 12),   # Increase x-axis text size
  axis.text.y = element_text(size = 12)    # Increase y-axis text size
)+scale_y_continuous(breaks = seq(0, 100, by = 10)) +  # Y-axis ticks every 10%
  scale_x_continuous(breaks = seq(1970, 2030, by = 5))  # X-axis ticks every 5 years
a9 <- ggplot()+geom_line(aes(time,list3333[[2]]*100), colour="pink", size=1)+geom_line(aes(time,list33333[[2]]*100), colour="red", size=1)+geom_line(aes(time,list333333[[2]]*100), colour="pink", size=1)+geom_line(aes(time,list9999[[2]]*100), colour="lightblue", size=1)+geom_line(aes(time,list99999[[2]]*100), colour="darkblue", size=1)+geom_line(aes(time,list999999[[2]]*100), colour="lightblue", size=1)+ coord_cartesian(xlim = c(1972.75, 2027.28), ylim = c(4.52,95.5)) + theme_classic() +labs(title="B: Savanes SIZ - 50% BMP with biannual MDA, including 100% VC efficacy" , x ="Year", y = "Microfilarial prevalence (%)") + geom_errorbar(data = SavanesYesEND2, aes(ymin=CILower, ymax=CIUpper, x = Survey_year), width=.1, size=0.5)+ geom_point(data = SavanesYesEND2,aes(x = Survey_year, y = PREV_cr,colour=SavanesYesEND2$ENDE), size = 4)+scale_colour_manual(values=c("orange"))+   labs(colour=NULL)+geom_vline(xintercept=1977,colour="green") + geom_vline(xintercept=1993,colour="dark green")+geom_vline(xintercept=1991,colour="red")+geom_vline(xintercept=2003,colour="ORANGE")+geom_text(mapping=aes(x=1991, y=80, label="MDA start"), size=4, angle=90, vjust=-0.4, hjust=0)+geom_text(mapping=aes(x=1993, y=80, label="VC end"), size=4, angle=90, vjust=-0.4, hjust=0)+geom_text(mapping=aes(x=1977, y=80, label="VC start"), size=4, angle=90, vjust=-0.4, hjust=0)+geom_text(mapping=aes(x=2003, y=65, label="Biannual MDA start"), size=4, angle=90, vjust=-0.4, hjust=0)+ theme(
  axis.title.x = element_text(size = 14),  # Increase x-axis title size
  axis.title.y = element_text(size = 14),  # Increase y-axis title size
  axis.text.x = element_text(size = 12),   # Increase x-axis text size
  axis.text.y = element_text(size = 12)    # Increase y-axis text size
)+scale_y_continuous(breaks = seq(0, 100, by = 10)) +  # Y-axis ticks every 10%
  scale_x_continuous(breaks = seq(1970, 2030, by = 5))  # X-axis ticks every 5 years

#hyperendemic under biannual MDA with and without 100% vector control (VC) efficacy
list4 <- lapply(means.examplesSSOB, `[[`, 3)
list44 <- lapply(means.examplesSSRB, `[[`, 3)
list444 <- lapply(means.examplesSSPB, `[[`, 3)
list4444 <- lapply(means.examplesSSO100, `[[`, 3)
list44444 <- lapply(means.examplesSSO100, `[[`, 3)
list444444 <- lapply(means.examplesSSO100, `[[`, 3)
a4 <- ggplot()+geom_line(aes(time,list4[[2]]*100), colour="lightblue", size=1)+geom_line(aes(time,list44[[2]]*100), colour="darkblue", size=1)+geom_line(aes(time,list444[[2]]*100), colour="lightblue", size=1)+geom_line(aes(time,list4444[[2]]*100), colour="pink", size=1)+geom_line(aes(time,list44444[[2]]*100), colour="red", size=1)+geom_line(aes(time,list444444[[2]]*100), colour="pink", size=1)+ coord_cartesian(xlim = c(1972.75, 2027.28), ylim = c(4.52,95.5)) + theme_classic() +labs(title="C: Savanes SIZ - 70% BMP with biannual MDA, including 100% VC efficacy",  x ="Year", y = "Microfilarial prevalence (%)") + geom_errorbar(data = SavanesYesEND3N, aes(ymin=CILower, ymax=CIUpper, x = Survey_year), width=.1, size=0.5)+ geom_point(data = SavanesYesEND3N,aes(x = Survey_year, y = PREV_cr,colour=ENDE), size = 4)+scale_colour_manual(values=c("brown"))+   labs(colour=NULL)+geom_vline(xintercept=1977,colour="green") + geom_vline(xintercept=1993,colour="dark green")+geom_vline(xintercept=1991,colour="red")+geom_vline(xintercept=2003,colour="ORANGE")+geom_text(mapping=aes(x=1991, y=80, label="MDA start"), size=4, angle=90, vjust=-0.4, hjust=0)+geom_text(mapping=aes(x=1993, y=80, label="VC end"), size=4, angle=90, vjust=-0.4, hjust=0)+geom_text(mapping=aes(x=1977, y=80, label="VC start"), size=4, angle=90, vjust=-0.4, hjust=0)+geom_text(mapping=aes(x=2003, y=65, label="Biannual MDA start"), size=4, angle=90, vjust=-0.4, hjust=0) +theme(
  axis.title.x = element_text(size = 14),  # Increase x-axis title size
  axis.title.y = element_text(size = 14),  # Increase y-axis title size
  axis.text.x = element_text(size = 12),   # Increase x-axis text size
  axis.text.y = element_text(size = 12)    # Increase y-axis text size
)+scale_y_continuous(breaks = seq(0, 100, by = 10)) +  # Y-axis ticks every 10%
  scale_x_continuous(breaks = seq(1970, 2030, by = 5))  # X-axis ticks every 5 years

combined_plotSIZ <- (a3 + a9) /  (a4+ plot_spacer()) 
#ggsave(filename = "USavanes RegionSIZNoBMP.png", plot = combined_plotSIZ, width = 15, height = 21*2/3, units = "in", dpi = 300)
combined_plotnoSIZ <- (a1+ plot_spacer()) 
#ggsave(filename = "USavanes RegionNoSIZNoBMP.png", plot = combined_plotnoSIZ, width = 15, height = 21*1/3, units = "in", dpi = 300)



#3. Kara SIZ (without baseline)
#Hypoendemicity under biannual MDA
list1 <- lapply(means.examplesKSOB, `[[`, 1) #A list with 4 numeric columns (plus the fifth). We just need to select the second one, that corresponds to the prevalence
list11 <- lapply(means.examplesKSRB, `[[`, 1) #A list with 4 numeric columns (plus the fifth). We just need to select the second one, that corresponds to the prevalence
list111 <- lapply(means.examplesKSPB, `[[`, 1) #A list with 4 numeric columns (plus the fifth). We just need to select the second one, that corresponds to the prevalence
a1 <- ggplot()+geom_line(aes(time,list1[[2]]*100), colour="lightblue", size=1)+geom_line(aes(time,list11[[2]]*100), colour="darkblue", size=1)+geom_line(aes(time,list111[[2]]*100), colour="lightblue", size=1)+ coord_cartesian(xlim = c(1972.75, 2027.28), ylim = c(4.52,95.5)) + theme_classic() +labs(title="A: Kara SIZ - 30% BMP, biannual MDA",  x ="Year", y = "Microfilarial prevalence (%)") + geom_errorbar(data = KaraYesEND1, aes(ymin=CILower, ymax=CIUpper, x = Survey_year), width=.1, size=0.5)+ geom_point(data = KaraYesEND1,aes(x = Survey_year, y = PREV_cr, colour=KaraYesEND1$ENDE), size = 4) +scale_colour_manual(values=color_group)+   labs(colour=NULL)+geom_vline(xintercept=1977,colour="green") + geom_vline(xintercept=2007,colour="dark green")+geom_vline(xintercept=1991,colour="red")+geom_text(mapping=aes(x=1991, y=80, label="MDA start"), size=4, angle=90, vjust=-0.4, hjust=0)+geom_text(mapping=aes(x=2007, y=80, label="VC end"), size=4, angle=90, vjust=-0.4, hjust=0)+geom_text(mapping=aes(x=1977, y=80, label="VC start"), size=4, angle=90, vjust=-0.4, hjust=0)+geom_text(mapping=aes(x=2003, y=65, label="Biannual MDA start"), size=4, angle=90, vjust=-0.4, hjust=0)+geom_vline(xintercept=2003,colour="ORANGE") +theme(
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
a2 <- ggplot()+geom_line(aes(time,list2[[2]]*100), colour="lightblue", size=1)+geom_line(aes(time,list22[[2]]*100), colour="darkblue", size=1)+geom_line(aes(time,list222[[2]]*100), colour="lightblue", size=1)+ coord_cartesian(xlim = c(1972.75, 2027.28), ylim = c(4.52,95.5)) + theme_classic() +labs(title="B: Kara SIZ - 50% BMP, biannual MDA",  x ="Year", y = "Microfilarial prevalence (%)") + geom_errorbar(data = KaraYesEND2, aes(ymin=CILower, ymax=CIUpper, x = Survey_year), width=.1, size=0.5)+ geom_point(data = KaraYesEND2,aes(x = Survey_year, y = PREV_cr, colour=KaraYesEND2$ENDE), size = 4) +scale_colour_manual(values="orange")+   labs(colour=NULL)+geom_vline(xintercept=1977,colour="green") + geom_vline(xintercept=2007,colour="dark green")+geom_vline(xintercept=1991,colour="red")+geom_text(mapping=aes(x=1991, y=80, label="MDA start"), size=4, angle=90, vjust=-0.4, hjust=0)+geom_text(mapping=aes(x=2007, y=80, label="VC end"), size=4, angle=90, vjust=-0.4, hjust=0)+geom_text(mapping=aes(x=1977, y=80, label="VC start"), size=4, angle=90, vjust=-0.4, hjust=0)+geom_text(mapping=aes(x=2003, y=65, label="Biannual MDA start"), size=4, angle=90, vjust=-0.4, hjust=0)+geom_vline(xintercept=2003,colour="ORANGE") + theme(
  axis.title.x = element_text(size = 14),  # Increase x-axis title size
  axis.title.y = element_text(size = 14),  # Increase y-axis title size
  axis.text.x = element_text(size = 12),   # Increase x-axis text size
  axis.text.y = element_text(size = 12)    # Increase y-axis text size
)+scale_y_continuous(breaks = seq(0, 100, by = 10)) +  # Y-axis ticks every 10%
  scale_x_continuous(breaks = seq(1970, 2030, by = 5))  # X-axis ticks every 5 years

#Hyperendemicity under biannual MDA
list3 <- lapply(means.examplesKNSOB, `[[`, 3)
list33 <- lapply(means.examplesKNSRB, `[[`, 3)
list333 <- lapply(means.examplesKNSPB, `[[`, 3)
a3 <- ggplot()+geom_line(aes(time,list3[[2]]*100), colour="lightblue", size=1)+geom_line(aes(time,list33[[2]]*100), colour="darkblue", size=1)+geom_line(aes(time,list333[[2]]*100), colour="lightblue", size=1)+ coord_cartesian(xlim = c(1972.75, 2027.28), ylim = c(4.52,95.5)) + theme_classic() +labs(title="C: Kara SIZ - 70% BMP, biannual MDA",  x ="Year", y = "Microfilarial prevalence (%)")+ geom_errorbar(data = KaraYesEND3, aes(ymin=CILower, ymax=CIUpper, x = Survey_year), width=.1, size=0.5) + geom_point(data = KaraYesEND3,aes(x = Survey_year, y = PREV_cr, colour=KaraYesEND3$ENDE), size = 4) +scale_colour_manual(values="dark red")+   labs(colour=NULL)+geom_vline(xintercept=1977,colour="green") + geom_vline(xintercept=2007,colour="dark green")+geom_vline(xintercept=1991,colour="red")+geom_text(mapping=aes(x=1991, y=80, label="MDA start"), size=4, angle=90, vjust=-0.4, hjust=0)+geom_text(mapping=aes(x=2007, y=80, label="VC end"), size=4, angle=90, vjust=-0.4, hjust=0)+geom_text(mapping=aes(x=1977, y=80, label="VC start"), size=4, angle=90, vjust=1.1, hjust=0)+geom_text(mapping=aes(x=2003, y=65, label="Biannual MDA start"), size=4, angle=90, vjust=-0.4, hjust=0)+geom_vline(xintercept=2003,colour="ORANGE") + theme(
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
a4 <- ggplot()+geom_line(aes(time,list4[[2]]*100), colour="lightblue", size=1)+geom_line(aes(time,list44[[2]]*100), colour="darkblue", size=1)+geom_line(aes(time,list444[[2]]*100), colour="lightblue", size=1)+ coord_cartesian(xlim = c(1972.75, 2027.28), ylim = c(4.52,95.5)) + theme_classic() +labs(title="D: Kara SIZ - 90% BMP, biannual MDA",  x ="Year", y = "Microfilarial prevalence (%)") + geom_errorbar(data = KaraYesEND4, aes(ymin=CILower, ymax=CIUpper, x = Survey_year), width=.1, size=0.5)+ geom_point(data = KaraYesEND4,aes(x = Survey_year, y = PREV_cr, colour=KaraYesEND4$ENDE), size = 4)+scale_colour_manual(values="dark red")+   labs(colour=NULL)+geom_vline(xintercept=1977,colour="green") + geom_vline(xintercept=2007,colour="dark green")+geom_vline(xintercept=1991,colour="red")+geom_text(mapping=aes(x=1991, y=50, label="MDA start"), size=4, angle=90, vjust=-0.4, hjust=0)+geom_text(mapping=aes(x=2007, y=80, label="VC end"), size=4, angle=90, vjust=-0.4, hjust=0)+geom_text(mapping=aes(x=1977, y=50, label="VC start"), size=4, angle=90, vjust=1.1, hjust=0)+geom_text(mapping=aes(x=2003, y=65, label="Biannual MDA start"), size=4, angle=90, vjust=-0.4, hjust=0)+geom_vline(xintercept=2003,colour="ORANGE") + theme(
  axis.title.x = element_text(size = 14),  # Increase x-axis title size
  axis.title.y = element_text(size = 14),  # Increase y-axis title size
  axis.text.x = element_text(size = 12),   # Increase x-axis text size
  axis.text.y = element_text(size = 12)    # Increase y-axis text size
)+scale_y_continuous(breaks = seq(0, 100, by = 10)) +  # Y-axis ticks every 10%
  scale_x_continuous(breaks = seq(1970, 2030, by = 5))  # X-axis ticks every 5 years

combined_plot <- (a1 + a2) /  (a3+a4)
#ggsave(filename = "UKara RegionNoBaseline.png", plot = combined_plot, width = 15, height = 21*2/3, units = "in", dpi = 300)



#5. Centrale non-SIZ without recorded baseline
#Hypoendemicity under annual MDA
list1 <- lapply(means.examplesKNSO, `[[`, 1) #A list with 4 numeric columns (plus the fifth). We just need to select the second one, that corresponds to the prevalence
list11 <- lapply(means.examplesKNSR, `[[`, 1) #A list with 4 numeric columns (plus the fifth). We just need to select the second one, that corresponds to the prevalence
list111 <- lapply(means.examplesKNSP, `[[`, 1) #A list with 4 numeric columns (plus the fifth). We just need to select the second one, that corresponds to the prevalence
a1 <- ggplot()+geom_line(aes(time,list1[[2]]*100), colour="lightblue", size=1)+geom_line(aes(time,list11[[2]]*100), colour="darkblue", size=1)+geom_line(aes(time,list111[[2]]*100), colour="lightblue", size=1)+ coord_cartesian(xlim = c(1972.75, 2027.28), ylim = c(4.52,95.5)) + theme_classic() +labs(title="A: Centrale non-SIZ - 30% BMP",  x ="Year", y = "Microfilarial prevalence (%)") + geom_errorbar(data = CentraleNoEND1, aes(ymin=CILower, ymax=CIUpper, x = Survey_year), width=.1, size=0.5)+ geom_point(data = CentraleNoEND1,aes(x = Survey_year, y = PREV_cr, colour=CentraleNoEND1$ENDE), size = 4)+scale_colour_manual(values=color_group)+   labs(colour=NULL) + geom_vline(xintercept=1989,colour="green") + geom_vline(xintercept=2002,colour="dark green")+geom_vline(xintercept=1991,colour="red")+geom_text(mapping=aes(x=1991, y=80, label="MDA start"), size=4, angle=90, vjust=-0.4, hjust=0)+geom_text(mapping=aes(x=2002, y=80, label="VC end"), size=4, angle=90, vjust=-0.4, hjust=0)+geom_text(mapping=aes(x=1989, y=80, label="VC start"), size=4, angle=90, vjust=-0.4, hjust=0)+ theme(
  axis.title.x = element_text(size = 14),  # Increase x-axis title size
  axis.title.y = element_text(size = 14),  # Increase y-axis title size
  axis.text.x = element_text(size = 12),   # Increase x-axis text size
  axis.text.y = element_text(size = 12)    # Increase y-axis text size
)+scale_y_continuous(breaks = seq(0, 100, by = 10)) +  # Y-axis ticks every 10%
  scale_x_continuous(breaks = seq(1970, 2030, by = 5))  # X-axis ticks every 5 years

#Mesoendemicity under annual MDA
list2 <- lapply(means.examplesKNSO, `[[`, 2)
list22 <- lapply(means.examplesKNSR, `[[`, 2)
list222 <- lapply(means.examplesKNSP, `[[`, 2)
a2 <- ggplot()+geom_line(aes(time,list2[[2]]*100), colour="lightblue", size=1)+geom_line(aes(time,list22[[2]]*100), colour="darkblue", size=1)+geom_line(aes(time,list222[[2]]*100), colour="lightblue", size=1)+ coord_cartesian(xlim = c(1972.75, 2027.28), ylim = c(4.52,95.5)) + theme_classic() +labs(title="B: Centrale non-SIZ - 50% BMP",  x ="Year", y = "Microfilarial prevalence (%)") + geom_errorbar(data = CentraleNoEND2, aes(ymin=CILower, ymax=CIUpper, x = Survey_year), width=.1, size=0.5)+ geom_point(data = CentraleNoEND2,aes(x = Survey_year, y = PREV_cr, colour=CentraleNoEND2$ENDE), size = 4) +scale_colour_manual(values="orange")+   labs(colour=NULL) + geom_vline(xintercept=1989,colour="green") + geom_vline(xintercept=2002,colour="dark green")+geom_vline(xintercept=1991,colour="red")+geom_text(mapping=aes(x=1991, y=80, label="MDA start"), size=4, angle=90, vjust=-0.4, hjust=0)+geom_text(mapping=aes(x=2002, y=80, label="VC end"), size=4, angle=90, vjust=-0.4, hjust=0)+geom_text(mapping=aes(x=1989, y=80, label="VC start"), size=4, angle=90, vjust=-0.4, hjust=0)+ theme(
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
a3 <- ggplot()+geom_line(aes(time,list3[[2]]*100), colour="lightblue", size=1)+geom_line(aes(time,list33[[2]]*100), colour="darkblue", size=1)+geom_line(aes(time,list333[[2]]*100), colour="lightblue", size=1)+ coord_cartesian(xlim = c(1972.75, 2027.28), ylim = c(4.52,95.5)) + theme_classic() +labs(title="C: Centrale non-SIZ - 70% BMP",  x ="Year", y = "Microfilarial prevalence (%)") + geom_errorbar(data = CentraleNoEND3, aes(ymin=CILower, ymax=CIUpper, x = Survey_year), width=.1, size=0.5)+ geom_point(data = CentraleNoEND3,aes(x = Survey_year, y = PREV_cr, colour=CentraleNoEND3$ENDE), size = 4) +scale_colour_manual(values="dark red")+   labs(colour=NULL) + geom_vline(xintercept=1989,colour="green") + geom_vline(xintercept=2002,colour="dark green")+geom_vline(xintercept=1991,colour="red")+geom_text(mapping=aes(x=1991, y=80, label="MDA start"), size=4, angle=90, vjust=1.1, hjust=0)+geom_text(mapping=aes(x=2002, y=80, label="VC end"), size=4, angle=90, vjust=-0.4, hjust=0)+geom_text(mapping=aes(x=1989, y=80, label="VC start"), size=4, angle=90, vjust=-0.4, hjust=0)+ theme(
  axis.title.x = element_text(size = 14),  # Increase x-axis title size
  axis.title.y = element_text(size = 14),  # Increase y-axis title size
  axis.text.x = element_text(size = 12),   # Increase x-axis text size
  axis.text.y = element_text(size = 12)    # Increase y-axis text size
)+scale_y_continuous(breaks = seq(0, 100, by = 10)) +  # Y-axis ticks every 10%
  scale_x_continuous(breaks = seq(1970, 2030, by = 5))  # X-axis ticks every 5 years

combined_plot <- (a1 + a2) /  (a3+ plot_spacer())
#ggsave(filename = "UUCentrale RegionNon-SIZNoBas.png", plot = combined_plot, width = 15, height = 21*2/3, units = "in", dpi = 300)


#6. Centrale SIZ without recorded baseline
#Hypoendemicity under biannual MDA
list1 <- lapply(means.examplesKSOB, `[[`, 1)
list11 <- lapply(means.examplesKSRB, `[[`, 1)
list111 <- lapply(means.examplesKSPB, `[[`, 1)
a1 <- ggplot()+geom_line(aes(time,list1[[2]]*100), colour="lightblue", size=1)+geom_line(aes(time,list11[[2]]*100), colour="darkblue", size=1)+geom_line(aes(time,list111[[2]]*100), colour="lightblue", size=1)+ coord_cartesian(xlim = c(1972.75, 2027.28), ylim = c(4.52,95.5)) + theme_classic() +labs(title="A: Centrale SIZ - 30% BMP, biannual MDA",  x ="Year", y = "Microfilarial prevalence (%)") + geom_errorbar(data = CentraleYesEND1, aes(ymin=CILower, ymax=CIUpper, x = Survey_year), width=.1, size=0.5)+ geom_point(data = CentraleYesEND1,aes(x = Survey_year, y = PREV_cr, colour=CentraleYesEND1$ENDE), size = 4) +scale_colour_manual(values=c(color_group))+   labs(colour=NULL)+geom_vline(xintercept=1977,colour="green") + geom_vline(xintercept=2007,colour="dark green")+geom_vline(xintercept=1991,colour="red")+geom_vline(xintercept=2003,colour="ORANGE")+geom_text(mapping=aes(x=1991, y=80, label="MDA start"), size=4, angle=90, vjust=-0.4, hjust=0)+geom_text(mapping=aes(x=2007, y=80, label="VC end"), size=4, angle=90, vjust=-0.4, hjust=0)+geom_text(mapping=aes(x=1977, y=80, label="VC start"), size=4, angle=90, vjust=-0.4, hjust=0)+geom_text(mapping=aes(x=2003, y=65, label="Biannual MDA start"), size=4, angle=90, vjust=-0.4, hjust=0)+ theme(
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
a3 <- ggplot()+geom_line(aes(time,list3[[2]]*100), colour="lightblue", size=1)+geom_line(aes(time,list33[[2]]*100), colour="darkblue", size=1)+geom_line(aes(time,list333[[2]]*100), colour="lightblue", size=1)+ coord_cartesian(xlim = c(1972.75, 2027.28), ylim = c(4.52,95.5)) + theme_classic() +labs(title="B: Centrale SIZ - 70% BMP, biannual MDA",  x ="Year", y = "Microfilarial prevalence (%)") + geom_errorbar(data = CentraleYesEND3, aes(ymin=CILower, ymax=CIUpper, x = Survey_year), width=.1, size=0.5)+ geom_point(data = CentraleYesEND3,aes(x = Survey_year, y = PREV_cr, colour=CentraleYesEND3$ENDE), size = 4) +scale_colour_manual(values=c("dark red"))+   labs(colour=NULL)+geom_vline(xintercept=1977,colour="green") + geom_vline(xintercept=2007,colour="dark green")+geom_vline(xintercept=1991,colour="red")+geom_vline(xintercept=2003,colour="ORANGE")+geom_text(mapping=aes(x=1991, y=80, label="MDA start"), size=4, angle=90, vjust=-0.4, hjust=0)+geom_text(mapping=aes(x=2007, y=80, label="VC end"), size=4, angle=90, vjust=-0.4, hjust=0)+geom_text(mapping=aes(x=1977, y=80, label="VC start"), size=4, angle=90, vjust=-0.4, hjust=0)+geom_text(mapping=aes(x=2003, y=65, label="Biannual MDA start"), size=4, angle=90, vjust=-0.4, hjust=0)+ theme(
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
a4 <- ggplot()+geom_line(aes(time,list4[[2]]*100), colour="lightblue", size=1)+geom_line(aes(time,list44[[2]]*100), colour="darkblue", size=1)+geom_line(aes(time,list444[[2]]*100), colour="lightblue", size=1)+ coord_cartesian(xlim = c(1972.75, 2027.28), ylim = c(4.52,95.5)) + theme_classic() +labs(title="C: Centrale SIZ - 90% BMP, biannual MDA",  x ="Year", y = "Microfilarial prevalence (%)") + geom_errorbar(data = CentraleYesEND4, aes(ymin=CILower, ymax=CIUpper, x = Survey_year), width=.1, size=0.5)+ geom_point(data = CentraleYesEND4,aes(x = Survey_year, y = PREV_cr, colour=CentraleYesEND4$ENDE), size = 4) +scale_colour_manual(values=c("dark red"))+   labs(colour=NULL)+geom_vline(xintercept=1977,colour="green") + geom_vline(xintercept=2007,colour="dark green")+geom_vline(xintercept=1991,colour="red")+geom_vline(xintercept=2003,colour="ORANGE")+geom_text(mapping=aes(x=1991, y=80, label="MDA start"), size=4, angle=90, vjust=-0.4, hjust=0)+geom_text(mapping=aes(x=2007, y=80, label="VC end"), size=4, angle=90, vjust=-0.4, hjust=0)+geom_text(mapping=aes(x=1977, y=80, label="VC start"), size=4, angle=90, vjust=-0.4, hjust=0)+geom_text(mapping=aes(x=2003, y=65, label="Biannual MDA start"), size=4, angle=90, vjust=-0.4, hjust=0)+ theme(
  axis.title.x = element_text(size = 14),  # Increase x-axis title size
  axis.title.y = element_text(size = 14),  # Increase y-axis title size
  axis.text.x = element_text(size = 12),   # Increase x-axis text size
  axis.text.y = element_text(size = 12)    # Increase y-axis text size
)+scale_y_continuous(breaks = seq(0, 100, by = 10)) +  # Y-axis ticks every 10%
  scale_x_continuous(breaks = seq(1970, 2030, by = 5))  # X-axis ticks every 5 years

combined_plot <- (a1 + a3) /  (a4+ plot_spacer())
#ggsave(filename = "UUCentrale RegionSIZNoBas.png", plot = combined_plot, width = 15, height = 21*2/3, units = "in", dpi = 300)


#7. Plateaux without pre-control recorded baseline
#hypoendemicity under annual MDA 
list1 <- lapply(means.examplesPO, `[[`, 1) #A list with 4 numeric columns (plus the fifth). We just need to select the second one, that corresponds to the prevalence
list11 <- lapply(means.examplesPR, `[[`, 1) #A list with 4 numeric columns (plus the fifth). We just need to select the second one, that corresponds to the prevalence
list111 <- lapply(means.examplesPP, `[[`, 1) #A list with 4 numeric columns (plus the fifth). We just need to select the second one, that corresponds to the prevalence
a1 <- ggplot()+geom_line(aes(time,list1[[2]]*100), colour="lightblue", size=1)+geom_line(aes(time,list11[[2]]*100), colour="darkblue", size=1)+geom_line(aes(time,list111[[2]]*100), colour="lightblue", size=1)+ coord_cartesian(xlim = c(1972.75, 2027.28), ylim = c(4.52,95.5)) + theme_classic() +labs(title="A: Plateaux - 30% BMP",  x ="Year", y = "Microfilarial prevalence (%)") + geom_errorbar(data = PlateauxNoEND1A, aes(ymin=CILower, ymax=CIUpper, x = Survey_year), width=.1, size=0.5)+ geom_point(data = PlateauxNoEND1A,aes(x = Survey_year, y = PREV_cr, colour=PlateauxNoEND1A$ENDE), size=4) +geom_vline(xintercept=1989,colour="green") + geom_vline(xintercept=2002,colour="dark green")+geom_vline(xintercept=1991,colour="red")+scale_colour_manual(values=color_group)+labs(colour=NULL)+geom_text(mapping=aes(x=1991, y=80, label="MDA start"), size=4, angle=90, vjust=--1.2, hjust=0)+geom_text(mapping=aes(x=2002, y=80, label="VC end"), size=4, angle=90, vjust=-0.4, hjust=0)+geom_text(mapping=aes(x=1989, y=80, label="VC start"), size=4, angle=90, vjust=-0.4, hjust=0)+ theme(
  axis.title.x = element_text(size = 14),  # Increase x-axis title size
  axis.title.y = element_text(size = 14),  # Increase y-axis title size
  axis.text.x = element_text(size = 12),   # Increase x-axis text size
  axis.text.y = element_text(size = 12)    # Increase y-axis text size
)+scale_y_continuous(breaks = seq(0, 100, by = 10)) +  # Y-axis ticks every 10%
  scale_x_continuous(breaks = seq(1970, 2030, by = 5))  # X-axis ticks every 5 years

#hypoendemicity under biannual MDA
list5 <- lapply(means.examplesPOB, `[[`, 1) #A list with 4 numeric columns (plus the fifth). We just need to select the second one, that corresponds to the prevalence
list55 <- lapply(means.examplesPRB, `[[`, 1) #A list with 4 numeric columns (plus the fifth). We just need to select the second one, that corresponds to the prevalence
list555 <- lapply(means.examplesPPB, `[[`, 1) #A list with 4 numeric columns (plus the fifth). We just need to select the second one, that corresponds to the prevalence
a5 <- ggplot()+geom_line(aes(time,list5[[2]]*100), colour="lightblue", size=1)+geom_line(aes(time,list55[[2]]*100), colour="darkblue", size=1)+geom_line(aes(time,list555[[2]]*100), colour="lightblue", size=1)+ coord_cartesian(xlim = c(1972.75, 2027.28), ylim = c(4.52,95.5)) + theme_classic() +labs(title="B: Plateaux - 30% BMP, biannual MDA",  x ="Year", y = "Microfilarial prevalence (%)") + geom_errorbar(data = PlateauxNoEND1B, aes(ymin=CILower, ymax=CIUpper, x = Survey_year), width=.1, size=0.5)+ geom_point(data = PlateauxNoEND1B,aes(x = Survey_year, y = PREV_cr, colour=PlateauxNoEND1B$ENDE), size=4) +geom_vline(xintercept=1989,colour="green") + geom_vline(xintercept=2002,colour="dark green")+geom_vline(xintercept=1991,colour="red")+scale_colour_manual(values=color_group)+labs(colour=NULL)+geom_text(mapping=aes(x=1991, y=80, label="MDA start"), size=4, angle=90, vjust=--1.2, hjust=0)+geom_text(mapping=aes(x=2002, y=80, label="VC end"), size=4, angle=90, vjust=-0.4, hjust=0)+geom_text(mapping=aes(x=1989, y=80, label="VC start"), size=4, angle=90, vjust=-0.4, hjust=0)+ geom_text(mapping=aes(x=2014, y=65, label="Biannual MDA start"), size=4, angle=90, vjust=--1.2, hjust=0)+geom_vline(xintercept=2014,colour="ORANGE") +theme(
  axis.title.x = element_text(size = 14),  # Increase x-axis title size
  axis.title.y = element_text(size = 14),  # Increase y-axis title size
  axis.text.x = element_text(size = 12),   # Increase x-axis text size
  axis.text.y = element_text(size = 12)    # Increase y-axis text size
)+scale_y_continuous(breaks = seq(0, 100, by = 10)) +  # Y-axis ticks every 10%
  scale_x_continuous(breaks = seq(1970, 2030, by = 5))  # X-axis ticks every 5 years

#mesoendemicity under annual MDA
list2 <- lapply(means.examplesPO, `[[`, 2)
list22 <- lapply(means.examplesPR, `[[`, 2)
list222 <- lapply(means.examplesPP, `[[`, 2)
a2 <- ggplot()+geom_line(aes(time,list2[[2]]*100), colour="lightblue", size=1)+geom_line(aes(time,list22[[2]]*100), colour="darkblue", size=1)+geom_line(aes(time,list222[[2]]*100), colour="lightblue", size=1)+ coord_cartesian(xlim = c(1972.75, 2027.28), ylim = c(4.52,95.5)) + theme_classic() +labs(title="C: Plateaux - 50% BMP",  x ="Year", y = "Microfilarial prevalence (%)") + geom_errorbar(data = PlateauxNoEND2A, aes(ymin=CILower, ymax=CIUpper, x = Survey_year), width=.1, size=0.5)+ geom_point(data = PlateauxNoEND2A,aes(x = Survey_year, y = PREV_cr, colour=PlateauxNoEND2A$ENDE), size=4) +geom_vline(xintercept=1989,colour="green") + geom_vline(xintercept=2002,colour="dark green")+geom_vline(xintercept=1991,colour="red")+scale_colour_manual(values=c("orange","dark red"))+   labs(colour=NULL)+geom_text(mapping=aes(x=1991, y=80, label="MDA start"), size=4, angle=90, vjust=+1.2, hjust=0)+geom_text(mapping=aes(x=2002, y=80, label="VC end"), size=4, angle=90, vjust=-0.4, hjust=0)+geom_text(mapping=aes(x=1989, y=80, label="VC start"), size=4, angle=90, vjust=-0.4, hjust=0)+ theme(
  axis.title.x = element_text(size = 14),  # Increase x-axis title size
  axis.title.y = element_text(size = 14),  # Increase y-axis title size
  axis.text.x = element_text(size = 12),   # Increase x-axis text size
  axis.text.y = element_text(size = 12)    # Increase y-axis text size
)+scale_y_continuous(breaks = seq(0, 100, by = 10)) +  # Y-axis ticks every 10%
  scale_x_continuous(breaks = seq(1970, 2030, by = 5))  # X-axis ticks every 5 years

#mesoendemicity under biannual MDA
PlateauxNo2B <- subset(PlateauxNo2B, !(VILNOM == "ADJABOULOUKOUKOPE" & CECI == 9999  &  POS_BASE == 21)) #Remove an unnusual value, same as the survey 1997.715
PlateauxNo2B <- subset(PlateauxNo2B, !(VILNOM == "KPODJI" )) #KPODJI is hyperendemic
list6 <- lapply(means.examplesPOB, `[[`, 2)
list66 <- lapply(means.examplesPRB, `[[`, 2)
list666 <- lapply(means.examplesPPB, `[[`, 2)
a6 <- ggplot()+geom_line(aes(time,list6[[2]]*100), colour="lightblue", size=1)+geom_line(aes(time,list66[[2]]*100), colour="darkblue", size=1)+geom_line(aes(time,list666[[2]]*100), colour="lightblue", size=1)+ coord_cartesian(xlim = c(1972.75, 2027.28), ylim = c(4.52,95.5)) + theme_classic() +labs(title="D: Plateaux - 50% BMP, biannual MDA",  x ="Year", y = "Microfilarial prevalence (%)") + geom_errorbar(data = PlateauxNoEND2B, aes(ymin=CILower, ymax=CIUpper, x = Survey_year), width=.1, size=0.5)+ geom_point(data = PlateauxNoEND2B,aes(x = Survey_year, y = PREV_cr, colour=PlateauxNoEND2B$ENDE), size=4) +geom_vline(xintercept=1989,colour="green") + geom_vline(xintercept=2002,colour="dark green")+geom_vline(xintercept=1991,colour="red")+scale_colour_manual(values=c("orange","dark red"))+   labs(colour=NULL)+geom_text(mapping=aes(x=1991, y=80, label="MDA start"), size=4, angle=90, vjust=+1.2, hjust=0)+geom_text(mapping=aes(x=2002, y=80, label="VC end"), size=4, angle=90, vjust=-0.4, hjust=0)+geom_text(mapping=aes(x=1989, y=80, label="VC start"), size=4, angle=90, vjust=-0.4, hjust=0)+ geom_text(mapping=aes(x=2014, y=65, label="Biannual MDA start"), size=4, angle=90, vjust=--1.2, hjust=0)+geom_vline(xintercept=2014,colour="ORANGE") + theme(
  axis.title.x = element_text(size = 14),  # Increase x-axis title size
  axis.title.y = element_text(size = 14),  # Increase y-axis title size
  axis.text.x = element_text(size = 12),   # Increase x-axis text size
  axis.text.y = element_text(size = 12)    # Increase y-axis text size
)+scale_y_continuous(breaks = seq(0, 100, by = 10)) +  # Y-axis ticks every 10%
  scale_x_continuous(breaks = seq(1970, 2030, by = 5))  # X-axis ticks every 5 years

#hyperendemicity under annual MDA
list3 <- lapply(means.examplesPO, `[[`, 3)
list33 <- lapply(means.examplesPR, `[[`, 3)
list333 <- lapply(means.examplesPP, `[[`, 3)
a3 <- ggplot()+geom_line(aes(time,list3[[2]]*100), colour="lightblue", size=1)+geom_line(aes(time,list33[[2]]*100), colour="darkblue", size=1)+geom_line(aes(time,list333[[2]]*100), colour="lightblue", size=1)+ coord_cartesian(xlim = c(1972.75, 2027.28), ylim = c(4.52,95.5)) + theme_classic() +labs(title="E: Plateaux - 70% BMP",  x ="Year", y = "Microfilarial prevalence (%)") + geom_errorbar(data = PlateauxNoEND3A, aes(ymin=CILower, ymax=CIUpper, x = Survey_year), width=.1, size=0.5)+ geom_point(data = PlateauxNoEND3A,aes(x = Survey_year, y = PREV_cr, colour=PlateauxNoEND3A$ENDE), size=4)+geom_vline(xintercept=1989,colour="green") + geom_vline(xintercept=2002,colour="dark green")+geom_vline(xintercept=1991,colour="red")+scale_colour_manual(values="dark red")+   labs(colour=NULL)+geom_text(mapping=aes(x=1991, y=80, label="MDA start"), size=4, angle=90, vjust=-0.4, hjust=0)+geom_text(mapping=aes(x=2002, y=80, label="VC end"), size=4, angle=90, vjust=-0.4, hjust=0)+geom_text(mapping=aes(x=1989, y=80, label="VC start"), size=4, angle=90, vjust=-0.4, hjust=0)+ theme(
  axis.title.x = element_text(size = 14),  # Increase x-axis title size
  axis.title.y = element_text(size = 14),  # Increase y-axis title size
  axis.text.x = element_text(size = 12),   # Increase x-axis text size
  axis.text.y = element_text(size = 12)    # Increase y-axis text size
)+scale_y_continuous(breaks = seq(0, 100, by = 10)) +  # Y-axis ticks every 10%
  scale_x_continuous(breaks = seq(1970, 2030, by = 5))  # X-axis ticks every 5 years

#hyperendemicity under biannual MDA
list7 <- lapply(means.examplesPOB, `[[`, 3)
list77 <- lapply(means.examplesPRB, `[[`, 3)
list777 <- lapply(means.examplesPPB, `[[`, 3)
a7 <- ggplot()+geom_line(aes(time,list7[[2]]*100), colour="lightblue", size=1)+geom_line(aes(time,list77[[2]]*100), colour="darkblue", size=1)+geom_line(aes(time,list777[[2]]*100), colour="lightblue", size=1)+ coord_cartesian(xlim = c(1972.75, 2027.28), ylim = c(4.52,95.5)) + theme_classic() +labs(title="F: Plateaux - 70% BMP, biannual MDA",  x ="Year", y = "Microfilarial prevalence (%)") + geom_errorbar(data = PlateauxNoEND3B, aes(ymin=CILower, ymax=CIUpper, x = Survey_year), width=.1, size=0.5)+ geom_point(data = PlateauxNoEND3B,aes(x = Survey_year, y = PREV_cr, colour=PlateauxNoEND3B$ENDE), size=4)+geom_vline(xintercept=1989,colour="green") + geom_vline(xintercept=2002,colour="dark green")+geom_vline(xintercept=1991,colour="red")+scale_colour_manual(values="dark red")+   labs(colour=NULL)+geom_text(mapping=aes(x=1991, y=80, label="MDA start"), size=4, angle=90, vjust=1.2, hjust=0)+geom_text(mapping=aes(x=2002, y=80, label="VC end"), size=4, angle=90, vjust=-0.4, hjust=0)+geom_text(mapping=aes(x=1989, y=80, label="VC start"), size=4, angle=90, vjust=-0.4, hjust=0)+geom_text(mapping=aes(x=2014, y=65, label="Biannual MDA start"), size=4, angle=90, vjust=-0.4, hjust=0)+geom_vline(xintercept=2014,colour="ORANGE")+ theme(
  axis.title.x = element_text(size = 14),  # Increase x-axis title size
  axis.title.y = element_text(size = 14),  # Increase y-axis title size
  axis.text.x = element_text(size = 12),   # Increase x-axis text size
  axis.text.y = element_text(size = 12)    # Increase y-axis text size
)+scale_y_continuous(breaks = seq(0, 100, by = 10)) +  # Y-axis ticks every 10%
  scale_x_continuous(breaks = seq(1970, 2030, by = 5))  # X-axis ticks every 5 years

combined_plot <- (a1 + a5) /  (a2+a6) /(a3 + a7)
#ggsave(filename = "UUPlateaux RegionNoBasAnieCorrected.png", plot = combined_plot, width = 15, height = 21, units = "in", dpi = 300)


#8. Maritime - until 2030 for villages with unknown baseline
#hypoendemicity under annual MDA
list1 <- lapply(means.examplesMO, `[[`, 1) #A list with 4 numeric columns (plus the fifth). We just need to select the second one, that corresponds to the prevalence
list11 <- lapply(means.examplesMR, `[[`, 1) #A list with 4 numeric columns (plus the fifth). We just need to select the second one, that corresponds to the prevalence
list111 <- lapply(means.examplesMP, `[[`, 1) #A list with 4 numeric columns (plus the fifth). We just need to select the second one, that corresponds to the prevalence
a1 <- ggplot()+geom_line(aes(time,list1[[2]]*100), colour="light blue", size=1)+geom_line(aes(time,list11[[2]]*100), colour="darkblue", size=1)+geom_line(aes(time,list111[[2]]*100), colour="lightblue", size=1)+ coord_cartesian(xlim = c(1972.75, 2027.28), ylim = c(4.52,95.5)) + theme_classic() +labs(title="A: Maritime - 30% BMP",  x ="Year", y = "Microfilarial prevalence (%)") + geom_errorbar(data = MaritimeNoEND1, aes(ymin=CILower, ymax=CIUpper, x = Survey_year), width=.1, size=0.5)+ geom_point(data = MaritimeNoEND1,aes(x = Survey_year, y = PREV_cr, colour=MaritimeNoEND1$ENDE), size=4)+geom_vline(xintercept=1988,colour="green") + geom_vline(xintercept=2002,colour="dark green")+geom_vline(xintercept=1991,colour="red")+scale_colour_manual(values=color_group)+   labs(colour=NULL)+geom_text(mapping=aes(x=1991, y=80, label="MDA start"), size=4, angle=90, vjust=-0.4, hjust=0)+geom_text(mapping=aes(x=2002, y=80, label="VC end"), size=4, angle=90, vjust=-0.4, hjust=0)+geom_text(mapping=aes(x=1988, y=80, label="VC start"), size=4, angle=90, vjust=-0.4, hjust=0)+ theme(
  axis.title.x = element_text(size = 14),  # Increase x-axis title size
  axis.title.y = element_text(size = 14),  # Increase y-axis title size
  axis.text.x = element_text(size = 12),   # Increase x-axis text size
  axis.text.y = element_text(size = 12)    # Increase y-axis text size
)+scale_y_continuous(breaks = seq(0, 100, by = 10)) +  # Y-axis ticks every 10%
  scale_x_continuous(breaks = seq(1970, 2030, by = 5))  # X-axis ticks every 5 years

#mesoendemicity under annual MDA
list2 <- lapply(means.examplesMO, `[[`, 2)
list22 <- lapply(means.examplesMR, `[[`, 2)
list222 <- lapply(means.examplesMP, `[[`, 2)
a2 <- ggplot()+geom_line(aes(time,list2[[2]]*100), colour="lightblue", size=1)+geom_line(aes(time,list22[[2]]*100), colour="darkblue", size=1)+geom_line(aes(time,list222[[2]]*100), colour="lightblue", size=1)+ coord_cartesian(xlim = c(1972.75, 2027.28), ylim = c(4.52,95.5)) + theme_classic() +labs(title="B: Maritime - 50% BMP",  x ="Year", y = "Microfilarial prevalence (%)")+ geom_errorbar(data = MaritimeNoEND2, aes(ymin=CILower, ymax=CIUpper, x = Survey_year), width=.1, size=0.5) + geom_point(data = MaritimeNoEND2,aes(x = Survey_year, y = PREV_cr, colour=MaritimeNoEND2$ENDE), size=4) +geom_vline(xintercept=1988,colour="green") + geom_vline(xintercept=2002,colour="dark green")+geom_vline(xintercept=1991,colour="red")+scale_colour_manual(values="orange")+   labs(colour=NULL)+geom_text(mapping=aes(x=1991, y=80, label="MDA start"), size=4, angle=90, vjust=-0.4, hjust=0)+geom_text(mapping=aes(x=2002, y=80, label="VC end"), size=4, angle=90, vjust=-0.4, hjust=0)+geom_text(mapping=aes(x=1988, y=80, label="VC start"), size=4, angle=90, vjust=-0.4, hjust=0)+ theme(
  axis.title.x = element_text(size = 14),  # Increase x-axis title size
  axis.title.y = element_text(size = 14),  # Increase y-axis title size
  axis.text.x = element_text(size = 12),   # Increase x-axis text size
  axis.text.y = element_text(size = 12)    # Increase y-axis text size
)+scale_y_continuous(breaks = seq(0, 100, by = 10)) +  # Y-axis ticks every 10%
  scale_x_continuous(breaks = seq(1970, 2030, by = 5))  # X-axis ticks every 5 years

#hyperendemicity under annual MDA
list3 <- lapply(means.examplesMO, `[[`, 3)
list33 <- lapply(means.examplesMR, `[[`, 3)
list333 <- lapply(means.examplesMP, `[[`, 3)
a3 <- ggplot()+geom_line(aes(time,list3[[2]]*100), colour="lightblue", size=1)+geom_line(aes(time,list33[[2]]*100), colour="darkblue", size=1)+geom_line(aes(time,list333[[2]]*100), colour="lightblue", size=1)+ coord_cartesian(xlim = c(1972.75, 2027.28), ylim = c(4.52,95.5)) + theme_classic() +labs(title="C: Maritime - 70% BMP",  x ="Year", y = "Microfilarial prevalence (%)")+ geom_errorbar(data = MaritimeNoEND3, aes(ymin=CILower, ymax=CIUpper, x = Survey_year), width=.1, size=0.5) + geom_point(data = MaritimeNoEND3,aes(x = Survey_year, y = PREV_cr, colour=MaritimeNoEND3$ENDE), size=4) +geom_vline(xintercept=1988,colour="green") + geom_vline(xintercept=2002,colour="dark green")+geom_vline(xintercept=1991,colour="red")+scale_colour_manual(values="dark red")+   labs(colour=NULL)+geom_text(mapping=aes(x=1991, y=80, label="MDA start"), size=4, angle=90, vjust=-0.4, hjust=0)+geom_text(mapping=aes(x=2002, y=80, label="VC end"), size=4, angle=90, vjust=-0.4, hjust=0)+geom_text(mapping=aes(x=1988, y=80, label="VC start"), size=4, angle=90, vjust=-0.4, hjust=0)+ theme(
  axis.title.x = element_text(size = 14),  # Increase x-axis title size
  axis.title.y = element_text(size = 14),  # Increase y-axis title size
  axis.text.x = element_text(size = 12),   # Increase x-axis text size
  axis.text.y = element_text(size = 12)    # Increase y-axis text size
)+scale_y_continuous(breaks = seq(0, 100, by = 10)) +  # Y-axis ticks every 10%
  scale_x_continuous(breaks = seq(1970, 2030, by = 5))  # X-axis ticks every 5 years

combined_plot <- (a1 + a2)/  (a3+ plot_spacer())
#ggsave(filename = "UMaritime RegionNoBas.png", plot = combined_plot, width = 15, height = 21*2/3, units = "in", dpi = 300)