setwd("~/Desktop/NYU_Projects/selectiontest_tips/analysis/sameQTL/focus4eQTL/")

library("plyr")
#require(data.table)
library('gridExtra')
library(cowplot)
library(ggplot2)
library(reshape2)

allinfo = read.table("tips.pos.allinfo.txt", header = T, sep="\t")

#osJAZ1 (JAP -- hapA)
i = 8
tip = allinfo$TIP[i]
fname=paste("selStat/win25/Selstats.",tip,".25SNPwin.JAP.txt", sep="")
file2 = read.table(fname, header=TRUE, sep="\t")
file2$Pos <- (file2$start+(file2$end-file2$start)/2)/1000000
fname=paste("H12/outputH12/",tip,".JAP.hapA_H12out_50_25.txt", sep="")
file2_ab = read.table(fname, header = F, sep="\t")
names(file2_ab) = c("Pos", "start", "stop", "k", "haploFreq", "indID", "H1", "H2", "H12", "H2/H1", "H123")

insert = allinfo[which(allinfo$TIP == tip),]$Pos/1000000
xlabel = paste(allinfo[which(allinfo$TIP == tip),]$Chr," (",tip,")",sep="")
st = allinfo$StartPos_Sel[i]/1000000
en = allinfo$EndPos_Sel[i]/1000000

fname=paste("permutationtest/",tip,".JAP.perm.txt",sep="")
perm = read.table(fname, header = TRUE, sep = "\t")
pi = mean(perm$pipr)
h12 = mean(perm$h12pr)

newpi2 <- ggplot(file2) + 
  geom_point(aes(x=Pos,y=pi_pr), color="darkolivegreen3", size=2) +
  geom_hline(yintercept = pi, color="black", size=1, linetype="dashed") + xlim(st, en) + 
  geom_vline(xintercept=insert, color="black", linetype = "dotted") + ylab("pi") + theme_bw() +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_rect(color = 'gray', size = 0.25), 
        axis.text.x = element_blank(), axis.title.x = element_blank(),
        axis.text.y = element_text(face="bold", size=10),axis.title.y = element_text(size=14),
        axis.ticks.length=unit(0.3, "cm"))

H122 <- ggplot() + 
  geom_point(data = file2_ab, aes(x=Pos/1000000,y=H12), color="darkolivegreen3", size=2) +
  geom_hline(yintercept = h12, color="black", size=1, linetype="dashed") + xlim(st, en) + 
  geom_vline(xintercept=insert, color="black", linetype = "dotted") + ylab("H12") + theme_bw() +
  xlab(xlabel) + 
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_rect(color = 'gray', size = 0.25), 
        axis.text.x = element_text(face="bold", size=10), axis.title.x = element_text(size=14),
        axis.text.y = element_text(face="bold", size=10), axis.title.y = element_text(size=14),
        axis.ticks.length=unit(0.3, "cm"))

#aligned1 <- align_plots(newpi2, H122, align = "v")
plot_grid(newpi2, H122, ncol = 1, align = "v")

#osJAZ1 (JAP -- hapB)
i = 8
tip = allinfo$TIP[i]
fname=paste("selStat/win25/Selstats.",tip,".25SNPwin.JAP.txt", sep="")
file2 = read.table(fname, header=TRUE, sep="\t")
file2$Pos <- (file2$start+(file2$end-file2$start)/2)/1000000
fname=paste("H12/outputH12/",tip,".JAP.hapB_H12out_50_25.txt", sep="")
file2_ab = read.table(fname, header = F, sep="\t")
names(file2_ab) = c("Pos", "start", "stop", "k", "haploFreq", "indID", "H1", "H2", "H12", "H2/H1", "H123")

insert = allinfo[which(allinfo$TIP == tip),]$Pos/1000000
xlabel = paste(allinfo[which(allinfo$TIP == tip),]$Chr," (",tip,")",sep="")
st = allinfo$StartPos_Sel[i]/1000000
en = allinfo$EndPos_Sel[i]/1000000

fname=paste("permutationtest/",tip,".JAP.perm.txt",sep="")
perm = read.table(fname, header = TRUE, sep = "\t")
pi = mean(perm$piab)
h12 = mean(perm$h12ab)

newpi2 <- ggplot(file2) + 
  geom_point(aes(x=Pos,y=pi_ab), color="coral3", size=2) +
  geom_hline(yintercept = pi, color="black", size=1, linetype="dashed") + xlim(st, en) + 
  geom_vline(xintercept=insert, color="black", linetype = "dotted") + ylab("pi") + theme_bw() +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_rect(color = 'gray', size = 0.25), 
        axis.text.x = element_blank(), axis.title.x = element_blank(),
        axis.text.y = element_text(face="bold", size=10),axis.title.y = element_text(size=14),
        axis.ticks.length=unit(0.3, "cm"))

H122 <- ggplot() + 
  geom_point(data = file2_ab, aes(x=Pos/1000000,y=H12), color="coral3", size=2) +
  geom_hline(yintercept = h12, color="black", size=1, linetype="dashed") + xlim(st, en) + 
  geom_vline(xintercept=insert, color="black", linetype = "dotted") + ylab("H12") + theme_bw() +
  xlab(xlabel) + 
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_rect(color = 'gray', size = 0.25), 
        axis.text.x = element_text(face="bold", size=10), axis.title.x = element_text(size=14),
        axis.text.y = element_text(face="bold", size=10), axis.title.y = element_text(size=14),
        axis.ticks.length=unit(0.3, "cm"))

#aligned1 <- align_plots(newpi2, H122, align = "v")
plot_grid(newpi2, H122, ncol = 1, align = "v")


#osJAZ1 (IND -- hapA)
i = 8
tip = allinfo$TIP[i]
fname=paste("selStat/win25/Selstats.",tip,".25SNPwin.IND.txt", sep="")
file2 = read.table(fname, header=TRUE, sep="\t")
file2$Pos <- (file2$start+(file2$end-file2$start)/2)/1000000
fname=paste("H12/outputH12/",tip,".IND.hapA_H12out_50_25.txt", sep="")
file2_ab = read.table(fname, header = F, sep="\t")
names(file2_ab) = c("Pos", "start", "stop", "k", "haploFreq", "indID", "H1", "H2", "H12", "H2/H1", "H123")

insert = allinfo[which(allinfo$TIP == tip),]$Pos/1000000
xlabel = paste(allinfo[which(allinfo$TIP == tip),]$Chr," (",tip,")",sep="")
st = allinfo$StartPos_Sel[i]/1000000
en = allinfo$EndPos_Sel[i]/1000000

fname=paste("permutationtest/",tip,".IND.perm.txt",sep="")
perm = read.table(fname, header = TRUE, sep = "\t")
pi = mean(perm$pipr)
h12 = mean(perm$h12pr)

newpi2 <- ggplot(file2) + 
  geom_point(aes(x=Pos,y=pi_pr), color="darkcyan", size=2) +
  geom_hline(yintercept = pi, color="black", size=1, linetype="dashed") + xlim(st, en) + 
  geom_vline(xintercept=insert, color="black", linetype = "dotted") + ylab("pi") + theme_bw() +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_rect(color = 'gray', size = 0.25), 
        axis.text.x = element_blank(), axis.title.x = element_blank(),
        axis.text.y = element_text(face="bold", size=10),axis.title.y = element_text(size=14),
        axis.ticks.length=unit(0.3, "cm"))

H122 <- ggplot() + 
  geom_point(data = file2_ab, aes(x=Pos/1000000,y=H12), color="darkcyan", size=2) +
  geom_hline(yintercept = h12, color="black", size=1, linetype="dashed") + xlim(st, en) + 
  geom_vline(xintercept=insert, color="black", linetype = "dotted") + ylab("H12") + theme_bw() +
  xlab(xlabel) + 
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_rect(color = 'gray', size = 0.25), 
        axis.text.x = element_text(face="bold", size=10), axis.title.x = element_text(size=14),
        axis.text.y = element_text(face="bold", size=10), axis.title.y = element_text(size=14),
        axis.ticks.length=unit(0.3, "cm"))

#aligned1 <- align_plots(newpi2, H122, align = "v")
plot_grid(newpi2, H122, ncol = 1, align = "v")


#osJAZ1 (IND -- hapB)
i = 8
tip = allinfo$TIP[i]
fname=paste("selStat/win25/Selstats.",tip,".25SNPwin.IND.txt", sep="")
file2 = read.table(fname, header=TRUE, sep="\t")
file2$Pos <- (file2$start+(file2$end-file2$start)/2)/1000000
fname=paste("H12/outputH12/",tip,".IND.hapB_H12out_50_25.txt", sep="")
file2_ab = read.table(fname, header = F, sep="\t")
names(file2_ab) = c("Pos", "start", "stop", "k", "haploFreq", "indID", "H1", "H2", "H12", "H2/H1", "H123")

insert = allinfo[which(allinfo$TIP == tip),]$Pos/1000000
xlabel = paste(allinfo[which(allinfo$TIP == tip),]$Chr," (",tip,")",sep="")
st = allinfo$StartPos_Sel[i]/1000000
en = allinfo$EndPos_Sel[i]/1000000

fname=paste("permutationtest/",tip,".IND.perm.txt",sep="")
perm = read.table(fname, header = TRUE, sep = "\t")
pi = mean(perm$piab)
h12 = mean(perm$h12ab)

newpi2 <- ggplot(file2) + 
  geom_point(aes(x=Pos,y=pi_ab), color="darkgoldenrod1", size=2) +
  geom_hline(yintercept = pi, color="black", size=1, linetype="dashed") + xlim(st, en) + 
  geom_vline(xintercept=insert, color="black", linetype = "dotted") + ylab("pi") + theme_bw() +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_rect(color = 'gray', size = 0.25), 
        axis.text.x = element_blank(), axis.title.x = element_blank(),
        axis.text.y = element_text(face="bold", size=10),axis.title.y = element_text(size=14),
        axis.ticks.length=unit(0.3, "cm"))

H122 <- ggplot() + 
  geom_point(data = file2_ab, aes(x=Pos/1000000,y=H12), color="darkgoldenrod1", size=2) +
  geom_hline(yintercept = h12, color="black", size=1, linetype="dashed") + xlim(st, en) + 
  geom_vline(xintercept=insert, color="black", linetype = "dotted") + ylab("H12") + theme_bw() +
  xlab(xlabel) + 
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_rect(color = 'gray', size = 0.25), 
        axis.text.x = element_text(face="bold", size=10), axis.title.x = element_text(size=14),
        axis.text.y = element_text(face="bold", size=10), axis.title.y = element_text(size=14),
        axis.ticks.length=unit(0.3, "cm"))

#aligned1 <- align_plots(newpi2, H122, align = "v")
plot_grid(newpi2, H122, ncol = 1, align = "v")


#osGAP (JAP -- hapA)
i = 6
tip = allinfo$TIP[i]
fname=paste("selStat/win25/Selstats.",tip,".25SNPwin.JAP.txt", sep="")
file2 = read.table(fname, header=TRUE, sep="\t")
file2$Pos <- (file2$start+(file2$end-file2$start)/2)/1000000
fname=paste("H12/outputH12/",tip,".JAP.present_H12out_50_25.txt", sep="")
file2_ab = read.table(fname, header = F, sep="\t")
names(file2_ab) = c("Pos", "start", "stop", "k", "haploFreq", "indID", "H1", "H2", "H12", "H2/H1", "H123")

insert = allinfo[which(allinfo$TIP == tip),]$Pos/1000000
xlabel = paste(allinfo[which(allinfo$TIP == tip),]$Chr," (","OsGAP",")",sep="")
st = allinfo$StartPos_Sel[i]/1000000
en = allinfo$EndPos_Sel[i]/1000000

fname=paste("permutationtest/",tip,".JAP.perm.txt",sep="")
perm = read.table(fname, header = TRUE, sep = "\t")
pi = mean(perm$pipr)
h12 = mean(perm$h12pr)

newpi2 <- ggplot(file2) + 
  geom_point(aes(x=Pos,y=pi_pr), color="darkolivegreen3", size=2) +
  geom_hline(yintercept = pi, color="black", size=1, linetype="dashed") + xlim(st, en) + 
  geom_vline(xintercept=insert, color="black", linetype = "dotted") + ylab("pi") + theme_bw() +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_rect(color = 'gray', size = 0.25), 
        axis.text.x = element_blank(), axis.title.x = element_blank(),
        axis.text.y = element_text(face="bold", size=10),axis.title.y = element_text(size=14),
        axis.ticks.length=unit(0.3, "cm"))

H122 <- ggplot() + 
  geom_point(data = file2_ab, aes(x=Pos/1000000,y=H12), color="darkolivegreen3", size=2) +
  geom_hline(yintercept = h12, color="black", size=1, linetype="dashed") + xlim(st, en) + 
  geom_vline(xintercept=insert, color="black", linetype = "dotted") + ylab("H12") + theme_bw() +
  xlab(xlabel) + 
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_rect(color = 'gray', size = 0.25), 
        axis.text.x = element_text(face="bold", size=10), axis.title.x = element_text(size=14),
        axis.text.y = element_text(face="bold", size=10), axis.title.y = element_text(size=14),
        axis.ticks.length=unit(0.3, "cm"))

plot_grid(newpi2, H122, ncol = 1, align = "v")


#osGAP (JAP -- hapB)
#Shows no sign of selection
i = 7
tip = allinfo$TIP[i]
fname=paste("selStat/win25/Selstats.",tip,".25SNPwin.JAP.txt", sep="")
file2 = read.table(fname, header=TRUE, sep="\t")
file2$Pos <- (file2$start+(file2$end-file2$start)/2)/1000000
fname=paste("H12/outputH12/",tip,".JAP.present_H12out_50_25.txt", sep="")
file2_ab = read.table(fname, header = F, sep="\t")
names(file2_ab) = c("Pos", "start", "stop", "k", "haploFreq", "indID", "H1", "H2", "H12", "H2/H1", "H123")

insert = allinfo[which(allinfo$TIP == tip),]$Pos/1000000
xlabel = paste(allinfo[which(allinfo$TIP == tip),]$Chr," (","OsGAP",")",sep="")
st = allinfo$StartPos_Sel[i]/1000000
en = allinfo$EndPos_Sel[i]/1000000

fname=paste("permutationtest/",tip,".JAP.perm.txt",sep="")
perm = read.table(fname, header = TRUE, sep = "\t")
pi = mean(perm$pipr)
h12 = mean(perm$h12pr)

newpi2 <- ggplot(file2) + 
  geom_point(aes(x=Pos,y=pi_pr), color="coral3", size=2) +
  geom_hline(yintercept = pi, color="black", size=1, linetype="dashed") + xlim(st, en) + 
  geom_vline(xintercept=insert, color="black", linetype = "dotted") + ylab("pi") + theme_bw() +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_rect(color = 'gray', size = 0.25), 
        axis.text.x = element_blank(), axis.title.x = element_blank(),
        axis.text.y = element_text(face="bold", size=10),axis.title.y = element_text(size=14),
        axis.ticks.length=unit(0.3, "cm"))

H122 <- ggplot() + 
  geom_point(data = file2_ab, aes(x=Pos/1000000,y=H12), color="coral3", size=2) +
  geom_hline(yintercept = h12, color="black", size=1, linetype="dashed") + xlim(st, en) + 
  geom_vline(xintercept=insert, color="black", linetype = "dotted") + ylab("H12") + theme_bw() +
  xlab(xlabel) + 
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_rect(color = 'gray', size = 0.25), 
        axis.text.x = element_text(face="bold", size=10), axis.title.x = element_text(size=14),
        axis.text.y = element_text(face="bold", size=10), axis.title.y = element_text(size=14),
        axis.ticks.length=unit(0.3, "cm"))

plot_grid(newpi2, H122, ncol = 1, align = "v")



#osGAP (IND -- hapA)
#shows sel but cannot be trusted since N=4 (too few ind)
i = 6
tip = allinfo$TIP[i]
fname=paste("selStat/win25/Selstats.",tip,".25SNPwin.IND.txt", sep="")
file2 = read.table(fname, header=TRUE, sep="\t")
file2$Pos <- (file2$start+(file2$end-file2$start)/2)/1000000
fname=paste("H12/outputH12/",tip,".IND.present_H12out_50_25.txt", sep="")
file2_ab = read.table(fname, header = F, sep="\t")
names(file2_ab) = c("Pos", "start", "stop", "k", "haploFreq", "indID", "H1", "H2", "H12", "H2/H1", "H123")

insert = allinfo[which(allinfo$TIP == tip),]$Pos/1000000
xlabel = paste(allinfo[which(allinfo$TIP == tip),]$Chr," (","OsGAP",")",sep="")
st = allinfo$StartPos_Sel[i]/1000000
en = allinfo$EndPos_Sel[i]/1000000

fname=paste("permutationtest/",tip,".IND.perm.txt",sep="")
perm = read.table(fname, header = TRUE, sep = "\t")
pi = mean(perm$pipr)
h12 = mean(perm$h12pr)

newpi2 <- ggplot(file2) + 
  geom_point(aes(x=Pos,y=pi_pr), color="darkcyan", size=2) +
  geom_hline(yintercept = pi, color="black", size=1, linetype="dashed") + xlim(st, en) + 
  geom_vline(xintercept=insert, color="black", linetype = "dotted") + ylab("pi") + theme_bw() +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_rect(color = 'gray', size = 0.25), 
        axis.text.x = element_blank(), axis.title.x = element_blank(),
        axis.text.y = element_text(face="bold", size=10),axis.title.y = element_text(size=14),
        axis.ticks.length=unit(0.3, "cm"))

H122 <- ggplot() + 
  geom_point(data = file2_ab, aes(x=Pos/1000000,y=H12), color="darkcyan", size=2) +
  geom_hline(yintercept = h12, color="black", size=1, linetype="dashed") + xlim(st, en) + 
  geom_vline(xintercept=insert, color="black", linetype = "dotted") + ylab("H12") + theme_bw() +
  xlab(xlabel) + 
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_rect(color = 'gray', size = 0.25), 
        axis.text.x = element_text(face="bold", size=10), axis.title.x = element_text(size=14),
        axis.text.y = element_text(face="bold", size=10), axis.title.y = element_text(size=14),
        axis.ticks.length=unit(0.3, "cm"))

plot_grid(newpi2, H122, ncol = 1, align = "v")


#osGAP (IND -- hapB)
#Shows no sign of selection
i = 7
tip = allinfo$TIP[i]
fname=paste("selStat/win25/Selstats.",tip,".25SNPwin.IND.txt", sep="")
file2 = read.table(fname, header=TRUE, sep="\t")
file2$Pos <- (file2$start+(file2$end-file2$start)/2)/1000000
fname=paste("H12/outputH12/",tip,".IND.present_H12out_50_25.txt", sep="")
file2_ab = read.table(fname, header = F, sep="\t")
names(file2_ab) = c("Pos", "start", "stop", "k", "haploFreq", "indID", "H1", "H2", "H12", "H2/H1", "H123")

insert = allinfo[which(allinfo$TIP == tip),]$Pos/1000000
xlabel = paste(allinfo[which(allinfo$TIP == tip),]$Chr," (","OsGAP",")",sep="")
st = allinfo$StartPos_Sel[i]/1000000
en = allinfo$EndPos_Sel[i]/1000000

fname=paste("permutationtest/",tip,".IND.perm.txt",sep="")
perm = read.table(fname, header = TRUE, sep = "\t")
pi = mean(perm$pipr)
h12 = mean(perm$h12pr)

newpi2 <- ggplot(file2) + 
  geom_point(aes(x=Pos,y=pi_pr), color="darkgoldenrod1", size=2) +
  geom_hline(yintercept = pi, color="black", size=1, linetype="dashed") + xlim(st, en) + 
  geom_vline(xintercept=insert, color="black", linetype = "dotted") + ylab("pi") + theme_bw() +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_rect(color = 'gray', size = 0.25), 
        axis.text.x = element_blank(), axis.title.x = element_blank(),
        axis.text.y = element_text(face="bold", size=10),axis.title.y = element_text(size=14),
        axis.ticks.length=unit(0.3, "cm"))

H122 <- ggplot() + 
  geom_point(data = file2_ab, aes(x=Pos/1000000,y=H12), color="darkgoldenrod1", size=2) +
  geom_hline(yintercept = h12, color="black", size=1, linetype="dashed") + xlim(st, en) + 
  geom_vline(xintercept=insert, color="black", linetype = "dotted") + ylab("H12") + theme_bw() +
  xlab(xlabel) + 
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_rect(color = 'gray', size = 0.25), 
        axis.text.x = element_text(face="bold", size=10), axis.title.x = element_text(size=14),
        axis.text.y = element_text(face="bold", size=10), axis.title.y = element_text(size=14),
        axis.ticks.length=unit(0.3, "cm"))

plot_grid(newpi2, H122, ncol = 1, align = "v")



#osMPH1 (JAP -- hapA)
i = 4
tip = allinfo$TIP[i]
fname=paste("selStat/win25/Selstats.",tip,".25SNPwin.JAP.txt", sep="")
file2 = read.table(fname, header=TRUE, sep="\t")
file2$Pos <- (file2$start+(file2$end-file2$start)/2)/1000000
fname=paste("H12/outputH12/",tip,".JAP.present_H12out_50_25.txt", sep="")
file2_ab = read.table(fname, header = F, sep="\t")
names(file2_ab) = c("Pos", "start", "stop", "k", "haploFreq", "indID", "H1", "H2", "H12", "H2/H1", "H123")

insert = allinfo[which(allinfo$TIP == tip),]$Pos/1000000
xlabel = paste(allinfo[which(allinfo$TIP == tip),]$Chr," (","OsMPH1",")",sep="")
st = allinfo$StartPos_Sel[i]/1000000
en = allinfo$EndPos_Sel[i]/1000000

fname=paste("permutationtest/",tip,".JAP.perm.txt",sep="")
perm = read.table(fname, header = TRUE, sep = "\t")
pi = mean(perm$pipr)
h12 = mean(perm$h12pr)

newpi2 <- ggplot(file2) + 
  geom_point(aes(x=Pos,y=pi_pr), color="darkolivegreen3", size=2) +
  geom_hline(yintercept = pi, color="black", size=1, linetype="dashed") + xlim(st, en) + 
  geom_vline(xintercept=insert, color="black", linetype = "dotted") + ylab("pi") + theme_bw() +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_rect(color = 'gray', size = 0.25), 
        axis.text.x = element_blank(), axis.title.x = element_blank(),
        axis.text.y = element_text(face="bold", size=10),axis.title.y = element_text(size=14),
        axis.ticks.length=unit(0.3, "cm"))

H122 <- ggplot() + 
  geom_point(data = file2_ab, aes(x=Pos/1000000,y=H12), color="darkolivegreen3", size=2) +
  geom_hline(yintercept = h12, color="black", size=1, linetype="dashed") + xlim(st, en) + 
  geom_vline(xintercept=insert, color="black", linetype = "dotted") + ylab("H12") + theme_bw() +
  xlab(xlabel) + 
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_rect(color = 'gray', size = 0.25), 
        axis.text.x = element_text(face="bold", size=10), axis.title.x = element_text(size=14),
        axis.text.y = element_text(face="bold", size=10), axis.title.y = element_text(size=14),
        axis.ticks.length=unit(0.3, "cm"))

plot_grid(newpi2, H122, ncol = 1, align = "v")


#osMPH1 (IND -- hapA)
i = 4
tip = allinfo$TIP[i]
fname=paste("selStat/win25/Selstats.",tip,".25SNPwin.IND.txt", sep="")
file2 = read.table(fname, header=TRUE, sep="\t")
file2$Pos <- (file2$start+(file2$end-file2$start)/2)/1000000
fname=paste("H12/outputH12/",tip,".IND.present_H12out_50_25.txt", sep="")
file2_ab = read.table(fname, header = F, sep="\t")
names(file2_ab) = c("Pos", "start", "stop", "k", "haploFreq", "indID", "H1", "H2", "H12", "H2/H1", "H123")

insert = allinfo[which(allinfo$TIP == tip),]$Pos/1000000
xlabel = paste(allinfo[which(allinfo$TIP == tip),]$Chr," (","OsMPH1",")",sep="")
st = allinfo$StartPos_Sel[i]/1000000
en = allinfo$EndPos_Sel[i]/1000000

fname=paste("permutationtest/",tip,".IND.perm.txt",sep="")
perm = read.table(fname, header = TRUE, sep = "\t")
pi = mean(perm$pipr)
h12 = mean(perm$h12pr)

newpi2 <- ggplot(file2) + 
  geom_point(aes(x=Pos,y=pi_pr), color="darkcyan", size=2) +
  geom_hline(yintercept = pi, color="black", size=1, linetype="dashed") + xlim(st, en) + 
  geom_vline(xintercept=insert, color="black", linetype = "dotted") + ylab("pi") + theme_bw() +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_rect(color = 'gray', size = 0.25), 
        axis.text.x = element_blank(), axis.title.x = element_blank(),
        axis.text.y = element_text(face="bold", size=10),axis.title.y = element_text(size=14),
        axis.ticks.length=unit(0.3, "cm"))

H122 <- ggplot() + 
  geom_point(data = file2_ab, aes(x=Pos/1000000,y=H12), color="darkcyan", size=2) +
  geom_hline(yintercept = h12, color="black", size=1, linetype="dashed") + xlim(st, en) + 
  geom_vline(xintercept=insert, color="black", linetype = "dotted") + ylab("H12") + theme_bw() +
  xlab(xlabel) + 
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_rect(color = 'gray', size = 0.25), 
        axis.text.x = element_text(face="bold", size=10), axis.title.x = element_text(size=14),
        axis.text.y = element_text(face="bold", size=10), axis.title.y = element_text(size=14),
        axis.ticks.length=unit(0.3, "cm"))

plot_grid(newpi2, H122, ncol = 1, align = "v")

