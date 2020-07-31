##
## 20 May 2020
## Calculation of percentages of (output) cells from historyOUT.txt (event-based output file)
##
## Simulation: Scenario 2. Affinity-based CD40 signal.  (BLIMP1 decision rule)
## 
## Queries from Antoine van Kampen and Elena Merino Tejero
##
## Version 2.0
##    Changes 
##       -min/max in plot from df.Dead.Ouy dataframe to ensure that color scaling will
##        be the same in all plots.
##       -jitter only in vertical direction
##       -Subsets that are plotted don't have overlap anymore (eg Ag+ do not contain PC)
##       -axis labels improved. Title of plot removed.
##       -ProcessEvents.Rmd is no more required. Reading of data and preprocessing of data 
##       -   occurs now in this script
##       -Replaced ggsave with tiff()/dev.off()
##
## Version 3.0 (12 June 2020)
##    Changes
##       - Elena has re-run simulations. THerefore queries were updated.
##       - Plots were updated (small changes)
##       - Single-cell lineage plots were changed and improved
##
## Version 3.1  (13 June 2020)
##    -Some cleaning up of the code and documentation
##
##
## Version 3.11 (13 June 2020)
##    - This file is a copy of ProcessEvents_Scenario2_Fixed_v3.1.R
##    - Adapted for Scenario 2, Affinity-based signal
##
## Version 3.11-Unix
##    -This script is a copy of Version 3.11 and contains some modification to run
##    -it on Unix (Virtual Machine, Ubuntu). This was necessary because the 
##    -historyOUT.txt file is about 47Gbyte and could not be processed on a laptop.
##    -In principle it should be possible to run on Ukkepuk if no other processes consume much 
##    -RAM

#
#clear workspace
#

if(TRUE) { #set to TRUE if data needs to be read again
  
rm(list=ls())

#library(RColorBrewer)
#library(furrr)
#library(tictoc)
#library(scales)
library(dplyr)
library(readr)
library(ggnewscale)
library(viridis)
library(tidyverse)
library(forcats)    #for reordering in the lineage plots
#


#Parameters for saving tiff files
#
#This will make a file that is 2400 x 2400 pixels, with an 300
#resolution gives you 2400/300 = 8 inches.
#
pwidth=2400
pheight=2400
punits="px"
psize=12
pres=300 #dpi

#Prefix to be used in the file name of all figures
prefix="S2-Affinity-"

#
## End plotting parameters


########
#Scenario 2. BLIMP1 as decision rule. Fixed affinity signal
########
#setwd("F:/Dropbox/_toprint/Discussion/Elena-Kopie/4- 9June2020 (Elena) - Copy/output_BLIMP1/bcinflow-affinity-cd40-50-PC-BLIMP1-8-MC-IamAgHigh_all events/Timecourse_scData")

## Lines copied from populationDynamics.Rmd (script of Elena)
## This is done to make this script independend of populationDynamics.RmD
#
cat("test1\n")
df.historyOut.raw <- read_delim("historyOut.txt","\t", escape_double = FALSE, trim_ws = TRUE)
cat("test2\n")

# DATA PREPROCESSING: Calculating number of time points per ID and change in BLIMP1 
# df.Dead.Out includes the dead cells and the output cells
df.Dead.Out <-df.historyOut.raw %>% 
  distinct(.keep_all = FALSE)%>% # Remove of born and divided events
  group_by(ID) %>%  
  arrange(time)%>%
  mutate(
    n_t_points = n() # number of time points per cell ID and time point
  )%>%
  ungroup()

df.Dead.Out$cellType <- gsub("1", "CBs",df.Dead.Out$cellType)
df.Dead.Out$cellType <- gsub("2", "CCs",df.Dead.Out$cellType)
df.Dead.Out$cellType <- gsub("6", "PCs",df.Dead.Out$cellType)
df.Dead.Out$cellType <- gsub("7", "MBCs",df.Dead.Out$cellType)
df.Dead.Out$IamAgHigh <- gsub("1", "OCs",df.Dead.Out$IamAgHigh)
df.Dead.Out$IamAgHigh <- gsub("0", "GC cells",df.Dead.Out$IamAgHigh)
#
##
## END COPY

} #END IF()

###
#    CALCULATION OF PERCENTAGES
###

###
# The next code is to determine the number of cells in each class.
# From these numbers the percentages are calculated in an Excel sheet
###
tmp <- df.Dead.Out %>% select(time, MID, ID,cellType,cellState, IamAgHigh,BLIMP1,Event) #just make copy of data.

###
# COUNT EVENTS
##
tmp1 = tmp  %>% distinct(ID, Event,.keep_all=TRUE) 
tmp2 = table(tmp1$Event)
print(tmp2)
###
#    PERCENTAGES FOR BLIMP1+
###

#Ag+BLIMP1+  --> PC; 13,110
p=tmp  %>% filter(Event=="become_plasma" & IamAgHigh=="OCs" & BLIMP1>8) %>% summarize(n()) 
#tmp1=tmp  %>% filter(Event=="become_plasma" & IamAgHigh=="OCs" & BLIMP1>8) 
cat("Ag+BLIMP1+ / PC ","\n")
print(p)

#Ag+BLIMP1+  --> NOT PC (=plasmablast); 50767
p=tmp  %>% filter(Event!="become_plasma" & IamAgHigh=="OCs" & BLIMP1>8) %>% select(ID) %>% distinct() %>% summarize(n()) 
#tmp2=tmp  %>% filter(Event!="become_plasma" & IamAgHigh=="OCs" & BLIMP1>8) %>% distinct(ID,.keep_all=TRUE) 
cat("Ag+BLIMP1+ / Not PC ","\n")
print(p)

#Ag-BLIMP1+ --> PC; 25,574
p=tmp  %>% filter(Event=="become_plasma" & IamAgHigh=="GC cells" & BLIMP1>8) %>% summarize(n())
#tmp3=tmp  %>% filter(Event=="become_plasma" & IamAgHigh=="GC cells" & BLIMP1>8) 
cat("Ag-BLIMP1+ / PC ","\n")
print(p)

#Ag-BLIMP1+ --> NOT PC (=plasmablast); 25,859
p=tmp  %>% filter(Event!="become_plasma" & IamAgHigh=="GC cells" & BLIMP1>8) %>% select(ID) %>% distinct() %>% summarize(n()) 
#tmp4=tmp  %>% filter(Event!="become_plasma" & IamAgHigh=="GC cells" & BLIMP1>8)  %>% distinct(ID,.keep_all=TRUE)  
cat("Ag-BLIMP1+ /Not PC","\n")
print(p)

#Note: in the above 4 queries some cells (IDs) occur in multiple sets. For example ID=4365 first has a born event, followed
#by a become_plasma event. Thus this ID is counted as Ag+BLIMP1+ in the OUTPUT and non-OUTPUT cell set.
#Therefore, if you count the number 

#The size of the joint set (tmp6) = 115310:
#tmp6=rbind(tmp1,tmp2,tmp3,tmp4)

#but the distinct IDs in this set = 93,500:
#tmp6 %>% distinct(ID)

#which is equal to the total number of BLIMP1+ cells.

#Total BLIMP1+; 93,500
p=tmp %>% filter(BLIMP1>8) %>% select(ID) %>% distinct() %>% summarize(n()) 
#tmp5=tmp %>% filter(BLIMP1>8) %>% distinct(ID,.keep_all=TRUE)  
cat("BLIMP1+ ","\n")
print(p)
                 
###
#    PERCENTAGES FOR BLIMP1-
###

#Ag+BLIMP1-  --> MBC; 0
p=tmp  %>% filter(Event=="become_memory" & IamAgHigh=="OCs" & BLIMP1<=8) %>% summarize(n())
#tmp1=tmp  %>% filter(Event=="become_memory" & IamAgHigh=="OCs" & BLIMP1<=8) 
cat("Ag+BLIMP1- / MBC ","\n")
print(p)

#Ag+BLIMP1-  --> NOT MBC (=CB);  135
p=tmp  %>% filter(Event!="become_memory" & IamAgHigh=="OCs" & BLIMP1<=8) %>% summarize(n())
#tmp2=tmp  %>% filter(Event!="become_memory" & IamAgHigh=="OCs" & BLIMP1<=8) 
cat("Ag+BLIMP1- / Not MBC ","\n")
print(p)

#Ag-BLIMP1- --> CB; this should be zero. Just to check  
p=tmp  %>% filter((Event=="become_plasma" | Event=="become_memory") & IamAgHigh=="GC cells" & BLIMP1<=8) %>% summarize(n())
#tmp3=tmp  %>% filter((Event=="become_plasma" | Event=="become_memory") & IamAgHigh=="GC cells" & BLIMP1<=8) 
cat("Ag-BLIMP1- / should be zero ","\n")
print(p)

#Ag-BLIMP1- --> CB; 161484 (not output)
p=tmp  %>% filter(Event!="become_plasma" & Event!="become_memory" & IamAgHigh=="GC cells" & BLIMP1<=8) %>% select(ID) %>% distinct() %>% summarize(n())
#tmp4=tmp  %>% filter(Event!="become_plasma" & Event!="become_memory" & IamAgHigh=="GC cells" & BLIMP1<=8) %>% distinct(ID,.keep_all = TRUE)
cat("Ag-BLIMP1- ","\n")
print(p)

#BLIMP1-; 161,484
p=tmp %>% filter(BLIMP1<=8) %>% select(ID) %>% distinct() %>% summarize(n()) 
#tmp5=tmp %>% filter(BLIMP1<=8) %>% distinct(ID,.keep_all=TRUE)  
cat("BLIMP1- ","\n")
print(p)

#The size of the joint set (tmp6) = 161,619
#tmp6=rbind(tmp1,tmp2,tmp3,tmp4)

#but the distinct IDs in this set = 161,484
#tmp6 %>% distinct(ID)

#Ag+; 50767
p=tmp %>% filter(IamAgHigh=="OCs") %>% distinct(ID) %>% summarize(n())
cat("Ag+ ","\n")
print(p)

#Ag-; 173862
p=tmp %>% filter(IamAgHigh=="GC cells") %>% distinct(ID) %>% summarize(n())
cat("Ag- ","\n")
print(p)

#All cells; 200,947
p=tmp %>%  distinct(ID) %>% summarize(n())
cat("All cells\n")
print(p)

###
#    FIGURES
###

##
# The code below plots subsets of the data. Basically, the selected subsets are
# in agreement with the subsets for determining the counts in each subset.
#
##

##
# First make all the subsets that we want to plot
##

#select all cells 
all <- df.Dead.Out %>% distinct(ID,.keep_all = TRUE)
  
#Select all Ag+ cells; 
AgH <- df.Dead.Out %>% filter(IamAgHigh=="OCs") %>% distinct(ID,.keep_all = TRUE)

#Select all Ag+ cells but not PC or MBC; Again to avoid bias in density of points
AgH_noOut <-df.Dead.Out %>% filter(Event!="become_plasma" & Event!="become_memory" & IamAgHigh=="OCs") %>% distinct(ID,.keep_all = TRUE)

# Select PCs 
pc <-df.Dead.Out  %>% filter(Event=="become_plasma") %>% distinct(ID,.keep_all = TRUE)

#Select Ag+BLIMP1- cells  --> MBC
mbc <- df.Dead.Out %>% filter(Event=="become_memory") %>% distinct(ID,.keep_all = TRUE)
  
#select all cells except PC, MBC, and Ag+
#Remaining is plot with ggplot. I do not plot 'All' since this subset also contains the
#other subsets and, therefore, would give a biased view in terms of the density of the
#points in the plot.
remaining <- df.Dead.Out %>% 
  filter(IamAgHigh!="OCs" &
            Event!= 'become_plasma' &
            Event != 'become_mbc') %>% distinct(ID,.keep_all = TRUE)


##
# Next do the actuall plotting
# All figures are written to a file
##

if(FALSE) {
  
#PLOT BLIMP vs AFFINITY (TIME gradient)
max=round(max(df.Dead.Out$time),0) # Taking the max/min from the tmp data frame ensures that all scales have the same range in all plots
min=round(min(df.Dead.Out$time),0)

tiff(paste(prefix,"DZ-cells-TimeGradient.tiff",sep=''), res=pres,width = pwidth, height = pheight, units = punits, pointsize = psize)
cat("tiff\n")
plot=ggplot(data=remaining,aes(x=Affinity, y=BLIMP1)) +
  geom_jitter(size=0.5,height=0,width=0.01,alpha=0.05) +
  geom_jitter(data=AgH_noOut,aes(x=Affinity,y=BLIMP1,color=time),height=0,width=0.05,alpha=0.05,size=1) +
  scale_color_gradient(name="Time (Ag+)",low="#05D9F6", high="#5011D1",limits=c(min,max),breaks=round(c(min,max/2,max),0))+ 
  new_scale_color() +
  geom_jitter(data=pc,aes(x=Affinity,y=BLIMP1,color=time),height=0,width=0.05,alpha=0.05, size=1) +
  scale_color_gradient(name="Time (PC)", low="white",high="red",limits=c(min,max),breaks=round(c(min,max/2,max),0)) +
  new_scale_color() +
  geom_jitter(data=mbc,aes(x=Affinity,y=BLIMP1,color=time),height=0,width=0.05,alpha=1, size=1.0) +
  scale_color_gradient(name="Time (MBC)", low="#66FF00",high="#006600",limits=c(min,max),breaks=round(c(min,max/2,max),0)) +
  geom_hline(yintercept=8,linetype="dashed",color="black")+
  ylim(-0.1,10)+ #if set to zero it will generate warnings about missing values
  theme(axis.title=element_text(size=14),axis.text=element_text(size=14))+
  labs(x="Affinity",y=expression("BLIMP1  ["~10^-8 ~"M]"), title="")+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
print(plot)
dev.off()

#PLOT BLIMP vs AFFINITY (Recycling gradient)
max=round(max(df.Dead.Out$Recyclings),2)
min=round(min(df.Dead.Out$Recyclings),2)

tiff(paste(prefix,"DZ-cells-RecyclingGradient.tiff",sep=''), res=pres, width = pwidth, height = pheight, units = punits, pointsize = psize)
plot=ggplot(data=remaining,aes(x=Affinity, y=BLIMP1)) +
  geom_jitter(size=0.5,width=0.01,alpha=0.05) +
  geom_jitter(data=AgH,aes(x=Affinity,y=BLIMP1,color=Recyclings),width=0.05,alpha=0.05,size=1) +
  scale_color_gradient(name="Recyclings (Ag+)",low="#05D9F6", high="#5011D1",limits=c(min,max),breaks=c(min,max/2,max))+ 
  new_scale_color() +
  geom_jitter(data=pc,aes(x=Affinity,y=BLIMP1,color=Recyclings),width=0.05,alpha=0.05, size=1) +
  scale_color_gradient(name="Recyclings (PC)", low="white",high="red",limits=c(min,max),breaks=c(min,max/2,max)) +
  new_scale_color() +
  geom_jitter(data=mbc,aes(x=Affinity,y=BLIMP1,color=Recyclings),width=0.05,alpha=1, size=1.0) +
  scale_color_gradient(name="Recyclings (MBC)", low="#66FF00",high="#006600",limits=c(min,max),breaks=c(min,max/2,max)) +
  geom_hline(yintercept=8,linetype="dashed",color="black")+
  ylim(-0.1,10)+ #if set to zero it will generate warnings about missing values
  theme(axis.title=element_text(size=14),axis.text=element_text(size=14))+
  labs(x="Affinity",y=expression("BLIMP1  ["~10^-8 ~"M]"), title="")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
print(plot)
dev.off()


#PLOT BLIMP vs AFFINITY (BCL6 gradient)
max=round(max(df.Dead.Out$BCL6),2)
min=round(min(df.Dead.Out$BCL6),2)

tiff(paste(prefix,"DZ-cells-BCL6Gradient.tiff",sep=''), res=pres, width = pwidth, height = pheight, units = punits, pointsize = psize)
plot=ggplot(data=remaining,aes(x=Affinity, y=BLIMP1)) +
  #geom_point(size=0.5, alpha=0.05)+ 
  geom_jitter(size=0.5,width=0.01,alpha=0.05) +
  geom_jitter(data=AgH,aes(x=Affinity,y=BLIMP1,color=BCL6),width=0.05,alpha=0.05,size=1) +
  scale_color_gradient(name="BCL6 (Ag+)",low="#05D9F6", high="#5011D1",limits=c(min,max),breaks=c(min,max/2,max))+ 
  new_scale_color() +
  geom_jitter(data=pc,aes(x=Affinity,y=BLIMP1,color=BCL6),width=0.05,alpha=0.05, size=1) +
  scale_color_gradient(name="BCL6 (PC)", low="white",high="red",limits=c(min,max),breaks=c(min,max/2,max)) +
  new_scale_color() +
  geom_jitter(data=mbc,aes(x=Affinity,y=BLIMP1,color=BCL6),width=0.05,alpha=1, size=1.0) +
  scale_color_gradient(name="BCL6 (MBC)", low="#66FF00",high="#006600",limits=c(min,max),breaks=c(min,max/2,max)) +
  geom_hline(yintercept=8,linetype="dashed",color="black")+
  ylim(-0.1,10)+ #if set to zero it will generate warnings about missing values
  theme(axis.title=element_text(size=14),axis.text=element_text(size=14))+
  labs(x="Affinity",y=expression("BLIMP1  ["~10^-8 ~"M]"), title="")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
print(plot)
dev.off()


#PLOT BLIMP vs AFFINITY (IRF4 gradient)
max=round(max(df.Dead.Out$IRF4),2)
min=round(min(df.Dead.Out$IRF4),2)

tiff(paste(prefix,"DZ-cells-IRF4Gradient.tiff",sep=''), res=pres, width = pwidth, height = pheight, units = punits, pointsize = psize)
plot=ggplot(data=remaining,aes(x=Affinity, y=BLIMP1)) +
  geom_jitter(size=0.5,width=0.01,alpha=0.05) +
  geom_jitter(data=AgH,aes(x=Affinity,y=BLIMP1,color=IRF4),width=0.05,alpha=0.05,size=1) +
  scale_color_gradient(name="IRF4 (Ag+)",low='#05D9F6', high='#5011D1',limits=c(min,max),breaks=c(min,max/2,max))+ 
  new_scale_color() +
  geom_jitter(data=pc,aes(x=Affinity,y=BLIMP1,color=IRF4),width=0.05,alpha=0.05, size=1) +
  scale_color_gradient(name="IRF4 (PC)", low="white",high="red",limits=c(min,max),breaks=c(min,max/2,max)) +
  new_scale_color() +
  geom_jitter(data=mbc,aes(x=Affinity,y=BLIMP1,color=IRF4),width=0.05,alpha=1, size=1.0) +
  scale_color_gradient(name="IRF4 (MBC)", low="#66FF00",high="#006600",limits=c(min,max),breaks=c(min,max/2,max)) +
  geom_hline(yintercept=8,linetype="dashed",color="black")+
  ylim(-0.1,10)+ #if set to zero it will generate warnings about missing values
  theme(axis.title=element_text(size=14),axis.text=element_text(size=14))+
  labs(x="Affinity",y=expression("BLIMP1  ["~10^-8 ~"M]"), title="")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
print(plot)
dev.off()


#####
# Select subsets of CCs
#####

cc <- df.Dead.Out %>%    #not plotted.
  filter(cellType =='CCs')
          
unselected <- df.Dead.Out %>% 
  filter(Event =='unselected')  %>% 
  select(time, ID, cellType, IamAgHigh, Event, BCL6, IRF4, BLIMP1, Affinity, Recyclings)

die <- df.Dead.Out %>% 
  filter(Event =='die')  %>% 
  select(time, ID, cellType, IamAgHigh, Event, BCL6, IRF4, BLIMP1, Affinity, Recyclings)

FDCsel <- df.Dead.Out %>% 
  filter(Event =='FDC_selected')  %>% 
  select(time, ID, cellType, IamAgHigh, Event, BCL6, IRF4, BLIMP1, Affinity, Recyclings)

StopS <- df.Dead.Out %>% 
  filter(Event =='stop_signaling_TC')  %>% 
  select(time, ID, cellType, IamAgHigh, Event, BCL6, IRF4, BLIMP1, Affinity, Recyclings)

RemainingCC <- df.Dead.Out %>% 
  filter((  Event =='die' |
            Event =='catch_FDC' |
            Event == 'start_contact_TC' |
            Event == 'start_signaling_TC' |
            Event == 'stop_contact_TC'))  %>% 
  select(time, ID, cellType, IamAgHigh, Event, BCL6, IRF4, BLIMP1, Affinity, Recyclings)

RemainingCC2 <- df.Dead.Out %>% #This also leaves out the dead cells
  filter((    Event =='catch_FDC' |
              Event == 'start_contact_TC' |
              Event == 'start_signaling_TC' |
              Event == 'stop_contact_TC'))  %>% 
  select(time, ID, cellType, IamAgHigh, Event, BCL6, IRF4, BLIMP1, Affinity, Recyclings)

RemainingCC3 <- df.Dead.Out %>% #no dead cells
  filter((    Event =='catch_FDC' |
              Event == 'start_contact_TC' |
              Event == 'start_signaling_TC' |
              Event == 'stop_contact_TC'))  %>% 
  select(time, ID, cellType, IamAgHigh, Event, BCL6, IRF4, BLIMP1, Affinity, Recyclings)

##
# Do the actual plotting
##

#PLOT BLIMP vs AFFINITY (TIME gradient) for different events of CC (exclude dead cells)
max=round(max(df.Dead.Out$time),0)
min=round(min(df.Dead.Out$time),0)

tiff(paste(prefix,"LZ-cells-TimeGradient.tiff",sep=''), res=pres, width = pwidth, height = pheight, units = punits, pointsize = psize)
plot=ggplot(data=RemainingCC,aes(x=Affinity, y=BLIMP1)) +
  geom_jitter(size=0.5,height=0,width=0.01,alpha=0.05) +
  geom_jitter(data=unselected,aes(x=Affinity,y=BLIMP1,color=time),height=0,width=0.05,alpha=0.05,size=1) +
  scale_color_gradient(name="Time (Unsel.)",low="#05D9F6", high="#5011D1",limits=c(min,max),breaks=round(c(min,max/2,max),0))+ 
  new_scale_color() +
  geom_jitter(data=FDCsel,aes(x=Affinity,y=BLIMP1,color=time),height=0,width=0.05,alpha=0.05, size=1) +
  scale_color_gradient(name="Time (FDC sel.)", low="white",high="red",limits=c(min,max),breaks=round(c(min,max/2,max),0)) +
  new_scale_color() +
  geom_jitter(data=StopS,aes(x=Affinity,y=BLIMP1,color=time),height=0,width=0.05,alpha=1, size=1.0) +
  scale_color_gradient(name="Time (Stop T)", low="#66FF00",high="#006600",limits=c(min,max),breaks=round(c(min,max/2,max),0)) +
  theme(axis.title=element_text(size=14),axis.text=element_text(size=14))+
  labs(x="Affinity",y=expression("BLIMP1  ["~10^-8 ~"M]"), title="")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
print(plot)
dev.off()

#PLOT BLIMP vs AFFINITY (TIME gradient) for dead cells only
max=round(max(df.Dead.Out$time),0)
min=round(min(df.Dead.Out$time),0)

tiff(paste(prefix,"LZ-cells-Dead-TimeGradient.tiff",sep=''), res=pres, width = pwidth, height = pheight, units = punits, pointsize = psize)
plot=ggplot(data=RemainingCC3,aes(x=Affinity, y=BLIMP1)) +
  geom_jitter(size=0.5,height=0,width=0.01,alpha=0.05) +
  geom_jitter(data=die,aes(x=Affinity,y=BLIMP1,color=time),height=0,width=0.05,alpha=0.05,size=1) +
  scale_color_gradient(name="Time (Dead)",low="yellow", high="gold",limits=c(min,max),breaks=round(c(min,max/2,max),0))+ 
  theme(axis.title=element_text(size=14),axis.text=element_text(size=14))+
  labs(x="Affinity",y=expression("BLIMP1  ["~10^-8 ~"M]"), title="")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
print(plot)
dev.off()

#PLOT BLIMP vs AFFINITY (TIME gradient) for all events
#Note: ensure that the plot window is large enough if ggsave is used. Otherwise part of the scales will
#not be visable. 
max=round(max(df.Dead.Out$time),0)
min=round(min(df.Dead.Out$time),0)

tiff(paste(prefix,"LZ-cells-All-TimeGradient.tiff",sep=''), res=pres,width = pwidth, height = pheight, units = punits, pointsize = psize)
plot=ggplot(data=RemainingCC2,aes(x=Affinity, y=BLIMP1)) +
  geom_jitter(size=0.5,height=0,width=0.01,alpha=0.05) +
  geom_jitter(data=die,aes(x=Affinity,y=BLIMP1,color=time),height=0,width=0.05,alpha=0.05,size=1) +
  scale_color_gradient(name="Time (Dead)",low="yellow", high="gold",limits=c(min,max),breaks=round(c(min,max/2,max),0))+ 
  new_scale_color() +
  geom_jitter(data=unselected,aes(x=Affinity,y=BLIMP1,color=time),height=0,width=0.05,alpha=0.05,size=1) +
  scale_color_gradient(name="Time (Unsel.)",low="#05D9F6", high="#5011D1",limits=c(min,max),breaks=round(c(min,max/2,max),0))+ 
  new_scale_color() +
  geom_jitter(data=FDCsel,aes(x=Affinity,y=BLIMP1,color=time),height=0,width=0.05,alpha=0.05, size=1) +
  scale_color_gradient(name="Time (FDC sel.)", low="white",high="red",limits=c(min,max),breaks=round(c(min,max/2,max),0)) +
  new_scale_color() +
  geom_jitter(data=StopS,aes(x=Affinity,y=BLIMP1,color=time),height=0,width=0.05,alpha=1, size=1.0) +
  scale_color_gradient(name="Time (Stop T)", low="#66FF00",high="#006600",limits=c(min,max),breaks=round(c(min,max/2,max),0)) +
  theme(axis.title=element_text(size=14),axis.text=element_text(size=14))+
  labs(x="Affinity",y=expression("BLIMP1  ["~10^-8 ~"M]"), title="")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
print(plot)
dev.off()


} #end if

#
#Other plots
#

if(TRUE) {  #make óther plots

  #
  # Output cell related plots
  # 
  
  # All cells vs Output cells
  df.do=df.Dead.Out
  
  tmp1=df.do  %>% filter(Event=="become_plasma" | Event == 'become_memory') 
  tmp1=cbind(tmp1,count=rep(1,dim(tmp1)[1]))
  #tmp1=tmp1  %>% filter(Event=="become_plasma" | Event=='become_memory')  #select output cells
  tmp1$cumulative = cumsum(tmp1$count/sum(tmp1$count))     #cummulative cell counts (scaled from 0 to 1)
  tmp1=tmp1 %>% mutate(cPC=cumsum(ifelse(Event=="become_plasma",1,0)))   #cumulative count of PCs
  tmp1=tmp1 %>% mutate(cMBC=cumsum(ifelse(Event=="become_memory",1,0))) #cumulative count of MBCs
  tmp1$cPC=tmp1$cPC/sum(tmp1$count) #scale at same range of cumulative cell counts
  tmp1$cMBC=tmp1$cMBC/sum(tmp1$count) 
  
  d=dim(tmp1)[1]
  
  All=df.do %>%
    mutate(day=cut_width(time,boundary=0,width=24)) %>% group_by(day) %>% count(Event=(Event=="born")) %>% filter(Event==TRUE)
  All=All[c(-1,-2),]  #remove first 2 days, because no output cells on first two days
  
  Out=df.do %>%
    mutate(day=cut_width(time,boundary=0,width=24)) %>% group_by(day) %>% 
    count(Event=(Event=="become_plasma")) %>% filter(Event==TRUE)
  All=cbind(All,Out)
  All=cbind(Day=c(3:(2+dim(All)[1])),All)
  colnames(All)=c('Day','day','Event','n1','day2','Event2', 'n2')
  
  r=round(cor(All$n1,All$n2),3)
  tiff(paste(prefix,"PC-OUT-Correlation.tiff",sep=''), res=pres,width = pwidth, height = pheight, units = punits, pointsize = psize)
  plot(All$n1,All$n2,type='b',col='red',xlab="Born cells",ylab="PC cells",main=paste("r=",r))
  dev.off()
  
 tiff(paste(prefix,"PC-OUT-timeprofile.tiff",sep=''), res=pres,width = pwidth, height = pheight, units = punits, pointsize = psize)
 plot=ggplot(All,aes(Day)) +
    geom_line(aes(y=n1),color='blue') +
    geom_line(aes(y=n2),color='darkgreen') +
    labs(x="Time (days)",y="#PCs",title=paste("r=",r))+
    theme(axis.title=element_text(size=14),axis.text=element_text(size=14))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
print(plot)
dev.off()


All=df.do %>%
  mutate(day=cut_width(time,boundary=0,width=24)) %>% group_by(day) %>% 
  count(Event=(Event=="born")) %>% filter(Event==TRUE)
All=All[c(-1,-2,-12,-13,-14,-15),]  #remove days on which no MBC is produced

Out=df.do %>%
  mutate(day=cut_width(time,boundary=0,width=24)) %>% group_by(day) %>% 
  count(Event=(Event=="become_memory")) %>% filter(Event==TRUE)
All=cbind(All,Out)
All=cbind(Day=c(2:(1+dim(All)[1])),All)
colnames(All)=c('Day','day','Event','n1','day2','Event2', 'n2')

r=round(cor(All$n1,All$n2),3)
tiff(paste(prefix,"MBC-OUT-Correlation.tiff",sep=''), res=pres,width = pwidth, height = pheight, units = punits, pointsize = psize)
plot(All$n1,All$n2,type='b',col='red',xlab="Born cells",ylab="MBC",main=paste("r=",r))
dev.off()

tiff(paste(prefix,"MBC-OUT-timeprofile.tiff",sep=''), res=pres,width = pwidth, height = pheight, units = punits, pointsize = psize)
plot=ggplot(All,aes(Day)) +
  geom_line(aes(y=n1),color='blue') +
  geom_line(aes(y=n2),color='darkgreen') +
  labs(x="Time (days)",y="#MBCs",title=paste("r=",r))+
  theme(axis.title=element_text(size=14),axis.text=element_text(size=14))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
print(plot)
dev.off()



tiff(paste(prefix,"PC-cumulative.tiff",sep=''), res=pres,width = pwidth, height = pheight, units = punits, pointsize = psize)
ymax=1500  #take max of histogram
plot=ggplot(data=tmp1,aes(x=time/24)) + 
    geom_histogram(binwidth = 1, fill='lightgray', col="white") +
    geom_line(aes(x=time/24,y=cumulative*ymax, colour="All"), lwd=1) +
    geom_line(aes(x=time/24,y=cPC*ymax, colour="PC"), lwd=1) +
    geom_line(aes(x=time/24,y=cMBC*ymax, colour="MBC"), lwd=1) +
    scale_color_manual(name='Cumulative',values=c('All'= "Black",'PC' = "Red",'MBC' = "Green"))+
    scale_y_continuous(name = 'Number of output cells', #limits=c(0,2500),
                       sec.axis = sec_axis(~./ymax,name="Cummulative output (%)"))+
    geom_hline(yintercept=c(.25,.50,.75,1)*ymax,linetype="dashed", color = "black")+
    labs(x="Time (days)",title="")+
    theme(axis.title=element_text(size=14),axis.text=element_text(size=14))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
print(plot)
dev.off()


tmp2=tmp1  %>% filter(Event=="become_plasma")  #select output cells
ac=round(sort(unique(log10(df.do$Affinity))),0) #determine affinity classes
t=table(log10(tmp2$Affinity))
min_affinity=as.numeric(unlist(dimnames(t))[1])
tiff(paste(prefix,"Hist-PC-affinity.tiff",sep=''), res=pres,width = pwidth, height = pheight, units = punits, pointsize = psize)
hist(log10(tmp2$Affinity),col='lightgreen',xlim=c(min(ac),max(ac)),xlab='log10 Affinity',main="Plasma cells",cex.lab=1.3,cex.axis=1.3)
abline(v=min_affinity,lty=2,col='red')
dev.off()

tmp2=tmp1  %>% filter(Event=="become_memory")  #select output cells
ac=round(sort(unique(log10(df.do$Affinity))),0) #determine affinity classes
t=table(log10(tmp2$Affinity))
min_affinity=as.numeric(unlist(dimnames(t))[1])
tiff(paste(prefix,"Hist-MBC-affinity.tiff",sep=''), res=pres,width = pwidth, height = pheight, units = punits, pointsize = psize)
hist(log10(tmp2$Affinity),col='lightgreen',xlim=c(min(ac),max(ac)),xlab='log10 Affinity',main="MBCs",cex.lab=1.3,cex.axis=1.3)
abline(v=min_affinity,lty=2,col='red')
dev.off()

#
# HISTOGRAMS
#

#Number of recyclings of MBCs
tiff(paste(prefix,"Hist-MBC-Recyclings.tiff",sep=''), res=pres,width = pwidth, height = pheight, units = punits, pointsize = psize)
hist(mbc$Recyclings,col='lightgreen',xlab='Recycles',main="Recycling of MBC",cex.lab=1.3,cex.axis=1.3)
dev.off()

#Number of recyclings of PCs
tiff(paste(prefix,"Hist-PC-Recyclings.tiff",sep=''), res=pres,width = pwidth, height = pheight, units = punits, pointsize = psize)
hist(pc$Recyclings,col='lightgreen',xlab='Recycles',main="Recycling of PC",cex.lab=1.3,cex.axis=1.3)
dev.off()

#Number of recyclings of Ag+ but not output cell
tiff(paste(prefix,"Hist-AgH-NoOut-Recyclings.tiff",sep=''), res=pres,width = pwidth, height = pheight, units = punits, pointsize = psize)
hist(AgH_noOut$Recyclings,col='lightgreen',xlab='Recycles',main="Recycling of Ag+ (no output)",cex.lab=1.3,cex.axis=1.3)
dev.off()


} #END IF


##
# Various plots for Single cell lineages
##
if(FALSE) {
## Single cell lineage
q=as.numeric(237813
             )  #PC #provide starting cell (=end of lineage)

lineage=c(q)
info=data.frame(df.Dead.Out %>% filter(ID==q))  #Compile all data of the lineage in data.frame.

#Trace back to founder cell by finding corresponding MID (parent cell)
while(q != -1) {
  mid = df.Dead.Out %>% filter(ID==q) %>% select(MID)
  mid=mid[1,]
  all = df.Dead.Out %>% filter(ID==q) 
  info=rbind(all,info)
  lineage=c(as.numeric(mid),lineage)
  q=as.numeric(mid)
}

#Extend the data frame 'info' with abbreviations of the various effects. To be used as labels in the plots
info=cbind(info,lab=info$Event)
for(i in 1:dim(info)[1]) {
  info[i,'lab']=''
  if (info[i,'Event']=='born') {info[i,'lab']='b'}
  if (info[i,'Event']=='divide') {info[i,'lab']='d'}
  if (info[i,'Event']=='recycling') {info[i,'lab']='r'}
  if (info[i,'Event']=='catch_FDC' | info[i,'Event']=='FDC_selected' | info[i,'Event']=='unselected' ) {info[i,'lab']='F'}
  if (info[i,'Event']=='start_contact_TC' | info[i,'Event']=='start_signaling_TC' | 
      info[i,'Event']=='stop_contact_TC'  | info[i,'Event']=='stop_signaling_TC' ) {info[i,'lab']='T'}
  if (info[i,'Event']=='become_plasma') {info[i,'lab']='P'}
  if (info[i,'Event']=='become_memory') {info[i,'lab']='M'}
  
}

lineage

#Use other graphic parameters for lineage plots
#This will make a file that is 2400 x 2400 pixels, with an 300
#resolution gives you 2400/300 = 8 inches.
#
pwidth=3200
pheight=pwidth
punits="px"
psize=12
pres=300 #dpi


#Plot lineage in Affinity vs BLIMP1 space
start=0
end=504
info2 <- info  %>% filter(time>=start & time<=end & lab!='F' & lab!='T') #exclude certain events if necessary
info2 <- info  %>% filter(time>=start & time<=end) #include everything between the given times


#Event vs Time, IDs, BLIMP1
#
#BLIMP level colored 
max=round(max(df.Dead.Out$BLIMP1),2) #Ensure same scaling for BLIMP1 in all plots
min=round(min(df.Dead.Out$BLIMP1),2) #Ensure same scaling for BLIMP1 in all plots

tiff(paste(prefix,"scLineage-EventTime.tiff",sep=''), res=pres,width = pwidth, height = pheight, units = punits, pointsize = psize)
plot=info2 %>% mutate(Event=fct_relevel(Event,"divide","born", "unselected", "catch_FDC", "FDC_selected", "start_contact_TC", "start_signaling_TC","stop_signaling_TC",
                                   "stop_contact_TC","recycling")) %>%
  ggplot(mapping=aes(x=time, y=Event)) +
  #ggplot(data=info2,aes(x=Event, y=time, color=BLIMP1)) +
  geom_point(aes(color=BLIMP1,size=Affinity)) +
  scale_size_continuous(range=c(1,6),limits=c(0,1)) +
  scale_color_gradient(name="BLIMP1",low='blue', high='red',limits=c(min,max),breaks=c(min,max/2,max))+ 
  geom_path(aes(group=ID),color="gray") +  #only draw lines between cells with same ID
  #geom_text(aes(label=lab),vjust=-0.5, hjust=1,size=5)+  #plot label of point
  #geom_text(aes(label=ID),vjust=-0.5, hjust=0,size=3)+   #plot ID of point
  #geom_text(aes(label=nFDCcontacts),vjust=-0.5, hjust=0,size=3)+   #plot ID of point
  geom_text(data=filter(info2, Event=="FDC_selected"),mapping=aes(label=nFDCcontacts),
            vjust=3.0, hjust=-0.1,color="black",size=3)+   #plot number of FDC contacts
  theme(axis.title=element_text(size=14),axis.text=element_text(size=14))+
  scale_y_discrete(labels=c("born"="Born","unselected"="Unselected","FDC_selected"="FDC sel." ,"catch_FDC" = "Catch Ag",
                            "start_contact_TC"="Tfh contact","start_signaling_TC"="Tfh signal",
                            "stop_signaling_TC"="Stop Tfh sign.","stop_contact_TC"="Stop Tfh cont.","recycling"="Recycle",
                            "divide"="Divide")) +
  labs(x="Time (hours)",y="Event", title="")+
  scale_x_continuous(sec.axis=dup_axis(~./24,name="Time (days)"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
print(plot)
dev.off()

#BLIMP1 vs Time, labels and IDs
#
#Events plotted with different colors, labels, and labels
tiff(paste(prefix,"scLineage-TimeBlimp1-labels.tiff",sep=''), res=pres,width = pwidth, height = pheight, units = punits, pointsize = psize)
plot=ggplot(data=info2,aes(x=time, y=BLIMP1, color=Event)) +
  geom_point(aes(color=Event),size=2) +
  geom_text(aes(label=lab),vjust=-0.5, hjust=1,size=5)+  #plot label of point
  geom_text(aes(label=ID),vjust=-0.5, hjust=0,size=3)+   #plot ID of point
  geom_line(color='black') +   
  geom_hline(yintercept=8,linetype="dashed",color="black")+
  #scale_y_continuous(trans = 'log10') + 
  ylim(-0.1,10)+ #if set to zero it will generate warnings about missing values
  theme(axis.title=element_text(size=14),axis.text=element_text(size=14))+
  labs(x="Time (hours)",y=expression("BLIMP1  ["~10^-8 ~"M]"), title="")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
print(plot)
dev.off()

#BLIMP1 vs Time, and affintiy
#
#Events plotted with different colors
#Point size indicates affinity
tiff(paste(prefix,"scLineage-TimeBlimp1-affinity.tiff",sep=''), res=pres,width = pwidth, height = pheight, units = punits, pointsize = psize)
plot=ggplot(data=info,aes(x=time, y=BLIMP1, color=Event)) +
  geom_point(aes(color=Event,size=Affinity)) +
  geom_line(color='black') +   
  geom_hline(yintercept=8,linetype="dashed",color="black")+
  ylim(-0.1,10)+ #if set to zero it will generate warnings about missing values
  theme(axis.title=element_text(size=14),axis.text=element_text(size=14))+
  labs(x="Time (hours)",y=expression("BLIMP1  ["~10^-8 ~"M]"), title="")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
print(plot)
dev.off()

#Affinity vs Time, events
#
#Events in different colors
tiff(paste(prefix,"scLineage-TimeAffinity-Events.tiff",sep=''), res=pres,width = pwidth, height = pheight, units = punits, pointsize = psize)
plot=ggplot(data=info,aes(x=time, y=Affinity, color=Event)) +
  geom_point(aes(color=Event),size=4) +
  geom_line(color='red') +   
  ylim(-0.1,1.0)+ #if set to zero it will generate warnings about missing values
  theme(axis.title=element_text(size=14),axis.text=element_text(size=14))+
  labs(x="Time (hours)",y=expression("Affinity (a.u.)"), title="")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
print(plot)
dev.off()


} # END IF()

