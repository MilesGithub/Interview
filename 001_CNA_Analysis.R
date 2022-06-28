
library(ggplot2)

data_dir<-'/.mounts/labs/reimandlab/private/projects/PPCG_CNA_pipeline/data/'
figure_dir<-'/.mounts/labs/reimandlab/private/projects/PPCG_CNA_pipeline/figures/'

#QC
SCNA_QC_Summary<-read.delim(paste(data_dir,'SCNA_QC_Summary_01_June_2020.csv',sep=""), sep=",")
SCNA_QC_Summary_table<-data.frame(table(SCNA_QC_Summary$SCNA_QC_Assessment))

#Donors and Samples
Sample_Donor<-read.table(paste(data_dir,'Sample_Donor_Tissue_Origin_01_June_2020.csv',sep=""), sep=",", header=T, fill = TRUE)
Sample_Donor_table<-data.frame(table(Sample_Donor$Tissue_Origin))

#Remove BPH, Normal, Recurrence, Unsatisfactory
Sample_Donor_filtered<-Sample_Donor[!(Sample_Donor$Tissue_Origin %in% c('BPH', 'Normal', 'Recurrence')),]
Sample_Donor_filtered<-Sample_Donor_filtered[!(Sample_Donor_filtered$PPCG_Sample_ID %in% SCNA_QC_Summary[SCNA_QC_Summary$SCNA_QC_Assessment == 'Unsatisfactory_SCNA_call',]$PPCG_Sample_ID), ]

#Get Donor and Sample lists for Metastasis and Primary  
PPCG_Donor_ID_Metastasis<-unique(Sample_Donor[Sample_Donor$Tissue_Origin == 'Metastasis',]$PPCG_Donor_ID)
PPCG_Sample_ID_Metastasis<-Sample_Donor[Sample_Donor$PPCG_Donor_ID %in% PPCG_Donor_ID_Metastasis,]$PPCG_Sample_ID

PPCG_Donor_ID_Primary<-unique(Sample_Donor$PPCG_Donor_ID)
PPCG_Donor_ID_Primary<-PPCG_Donor_ID_Primary[!(PPCG_Donor_ID_Primary %in% PPCG_Donor_ID_Metastasis)]
PPCG_Sample_ID_Primary<-Sample_Donor[Sample_Donor$PPCG_Donor_ID %in% PPCG_Donor_ID_Primary,]$PPCG_Sample_ID

#Get Number of Samples per Donor
PPCG_donors<-unique(Sample_Donor$PPCG_Donor_ID)
PPCG_donor_sample_num<-data.frame(unique(Sample_Donor$PPCG_Donor_ID))
PPCG_donor_sample_num$Primary_sample_num<-0
PPCG_donor_sample_num$Metastasis_sample_num<-0
colnames(PPCG_donor_sample_num)<-c('PPCG.donor.ID', 'Primary_num', 'Metastasis_num')

for (i in 1:length(PPCG_donors)){
  print(i)
  donor<-PPCG_donors[i]
  Primary_sub<-Sample_Donor[Sample_Donor$PPCG_Donor_ID == donor & Sample_Donor$Tissue_Origin == 'Primary',]
  Metastasis_sub<-Sample_Donor[Sample_Donor$PPCG_Donor_ID == donor & Sample_Donor$Tissue_Origin == 'Metastasis',]  
  Primary_num<-nrow(Primary_sub)
  Metastasis_num<-nrow(Metastasis_sub)
  PPCG_donor_sample_num[PPCG_donor_sample_num$PPCG.donor.ID == donor, 2]<-Primary_num
  PPCG_donor_sample_num[PPCG_donor_sample_num$PPCG.donor.ID == donor, 3]<-Metastasis_num
}

PPCG_donor_sample_num$all_samples<-PPCG_donor_sample_num$Primary_num + PPCG_donor_sample_num$Metastasis_num
PPCG_donor_sample_num<-PPCG_donor_sample_num[order(-PPCG_donor_sample_num$all_samples),]
PPCG_donor_sample_num_Primary<-PPCG_donor_sample_num[,c(1,2)]
PPCG_donor_sample_num_Primary$type<-'Primary'
colnames(PPCG_donor_sample_num_Primary)<-c('donor_id', 'Freq', 'type')
PPCG_donor_sample_num_Metastasis<-PPCG_donor_sample_num[,c(1,3)]
PPCG_donor_sample_num_Metastasis$type<-'Metastasis'
colnames(PPCG_donor_sample_num_Metastasis)<-c('donor_id', 'Freq', 'type')
PPCG_donor_sample_num_figure<-rbind(PPCG_donor_sample_num_Primary, PPCG_donor_sample_num_Metastasis)

PPCG_donor_sample_num_figure$donor_id <- factor(PPCG_donor_sample_num_figure$donor_id, levels=PPCG_donor_sample_num$PPCG.donor.ID)
PPCG_donor_sample_num_figure_Primary<-PPCG_donor_sample_num_figure[as.character(PPCG_donor_sample_num_figure$donor_id) %in% PPCG_Donor_ID_Primary,]
PPCG_donor_sample_num_figure_Metastasis<-PPCG_donor_sample_num_figure[as.character(PPCG_donor_sample_num_figure$donor_id) %in% PPCG_Donor_ID_Metastasis,]

#Samples per Donor figures
#Primary
p<-ggplot(PPCG_donor_sample_num_figure_Primary, aes(x=donor_id, y=Freq, fill=type)) +
  geom_bar(stat="identity", fill="#008080") +
  ylab('Number of Samples')+
  xlab('Donor')+
  scale_y_continuous(limits = c(0, 8), breaks = (seq(0,8,by = 1)))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1),
        panel.background = element_blank(), axis.line = element_line(color = "#2d3f44"),
        legend.title = element_blank(),
        legend.position = "none")

pdf(paste(c(figure_dir, '001_Samples_per_Donor_Primary.pdf'), collapse=''), width=100, height= 4)
print(p)
dev.off()

#Metastasis
p<-ggplot(PPCG_donor_sample_num_figure_Metastasis, aes(x=donor_id, y=Freq, fill=type)) +
  geom_bar(stat="identity") +
  ylab('Number of Samples')+
  xlab('Donor')+
  scale_y_continuous(limits = c(0, 10), breaks = (seq(0,10,by = 1)))+
  scale_fill_manual(values=c("#c83771", "#008080"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1),
        panel.background = element_blank(), axis.line = element_line(color = "#2d3f44"),
        legend.title = element_blank(),
        legend.position = "none")

pdf(paste(c(figure_dir, '002_Samples_per_Donor_Metastasis.pdf'), collapse=''), width=20, height= 4)
print(p)
dev.off()


#Get CNA data
file_dir<-'D:/projects/PPCG_project/data/Subclonal_SCNA_01_June_2020/'
files<-list.files(file_dir)
PPCG_CNAs<-data.frame(matrix(nrow=0,ncol=69))
for(i in 1:length(files)){
  print(i)
  file<-files[i]
  tab<-read.delim(paste(file_dir,file,sep=""))
  tab$file<-file
  sample_id<-substr(file, start = 1, stop = 13)
  Sample_Donor_sub<-Sample_Donor[Sample_Donor$PPCG_Sample_ID == sample_id,]
  donor_id<-Sample_Donor_sub$PPCG_Donor_ID
  tab$donor_id<-donor_id
  tab$sample_id<-sample_id
  PPCG_CNAs<-rbind(PPCG_CNAs,tab)
}

PPCG_CNAs$Seg.CN<-PPCG_CNAs$nMaj1_A + PPCG_CNAs$nMin1_A
PPCG_CNAs$chromosome<-paste('chr',PPCG_CNAs$chr, sep='')
PPCG_CNAs$Num.Markers<-NA
PPCG_CNAs_sub<-PPCG_CNAs[,c(69,70,72,2,3,73,71)]
colnames(PPCG_CNAs_sub)<-c("Donor", "Sample", "Chromosome", "Start.bp", "End.bp", "Num.Markers", "Seg.CN")

PPCG_CNA_Metastasis<-PPCG_CNAs_sub[PPCG_CNAs_sub$Donor %in% PPCG_Donor_ID_Metastasis,]
PPCG_CNA_Primary<-PPCG_CNAs_sub[PPCG_CNAs_sub$Donor %in% PPCG_Donor_ID_Primary,]

#Get avg CNA for each sample
Primary_samples_CNA_avg<-data.frame(PPCG_Sample_ID_Primary,0,0)
Metastasis_samples_CNA_avg<-data.frame(PPCG_Sample_ID_Metastasis,0,0)
colnames(Primary_samples_CNA_avg)<-c('Sample_id', 'calculated_avg_CN', 'avg_CN')
colnames(Metastasis_samples_CNA_avg)<-c('Sample_id', 'calculated_avg_CN', 'avg_CN')

for(i in 1:length(PPCG_Sample_ID_Primary)){
  print(i)
  Primary_sample<-Primary_samples_CNA_avg[i,1]
  sample_CNAs<-PPCG_CNA_Primary[PPCG_CNA_Primary$Sample == Primary_sample,]
  seg_cn<-as.numeric(as.character(sample_CNAs$Seg.CN))
  seg_cn<-seg_cn[!is.na(seg_cn)]
  sample_CNA_mean<-median(seg_cn)
  Primary_samples_CNA_avg[i,2]<-sample_CNA_mean
  if(sample_CNA_mean > 2){
    Primary_samples_CNA_avg[i,3]<-4
  }else{
    Primary_samples_CNA_avg[i,3]<-2
  }
}

for(i in 1:length(PPCG_Sample_ID_Metastasis)){
  print(i)
  Metastasis_sample<-Metastasis_samples_CNA_avg[i,1]
  sample_CNAs<-PPCG_CNA_Metastasis[PPCG_CNA_Metastasis$Sample == Metastasis_sample,]
  seg_cn<-as.numeric(as.character(sample_CNAs$Seg.CN))
  seg_cn<-seg_cn[!is.na(seg_cn)]
  sample_CNA_mean<-median(seg_cn)
  Metastasis_samples_CNA_avg[i,2]<-sample_CNA_mean
  if(sample_CNA_mean > 2){
    Metastasis_samples_CNA_avg[i,3]<-4
  }else{
    Metastasis_samples_CNA_avg[i,3]<-2
  }
}

#Get relative CNA
PPCG_CNA_Primary$relative_CN<-0
for(i in 1:nrow(PPCG_CNA_Primary)){
  print(i)
  Primary_sample<-PPCG_CNA_Primary[i,]
  Sample<-Primary_sample$Sample
  seg_cn<-Primary_sample$Seg.CN
  sample_sub<-Primary_samples_CNA_avg[Primary_samples_CNA_avg$PPCG_Sample_ID_Primary == Sample,]
  avg_CN<-sample_sub$avg_CN
  PPCG_CNA_Primary[i,]$relative_CN<-seg_cn -avg_CN
}

PPCG_CNA_Metastasis$relative_CN<-0
for(i in 1:nrow(PPCG_CNA_Metastasis)){
  print(i)
  Metastasis_sample<-PPCG_CNA_Metastasis[i,]
  Sample<-Metastasis_sample$Sample
  seg_cn<-Metastasis_sample$Seg.CN
  sample_sub<-Metastasis_samples_CNA_avg[Metastasis_samples_CNA_avg$PPCG_Sample_ID_Metastasis == Sample,]
  avg_CN<-sample_sub$avg_CN
  PPCG_CNA_Metastasis[i,]$relative_CN<-seg_cn -avg_CN
}

PPCG_CNA_Metastasis_CNs<-PPCG_CNA_Metastasis[abs(as.numeric(as.character(PPCG_CNA_Metastasis$relative_CN))) > 0,]
PPCG_CNA_Primary_CNs<-PPCG_CNA_Primary[abs(as.numeric(as.character(PPCG_CNA_Primary$relative_CN))) > 0,]

#Get number of CNAs a Sample
Sample_Donor_Primary<-Sample_Donor[Sample_Donor$Tissue_Origin == 'Primary',]
Sample_Donor_Primary$num_AMP_CNA<-NA
Sample_Donor_Primary$num_DEL_CNA<-NA
Sample_Donor_Primary$num_all_CNA<-NA

Sample_Donor_Metastasis<-Sample_Donor[Sample_Donor$Tissue_Origin == 'Metastasis',]
Sample_Donor_Metastasis$num_AMP_CNA<-NA
Sample_Donor_Metastasis$num_DEL_CNA<-NA
Sample_Donor_Metastasis$num_all_CNA<-NA

for(i in 1:nrow(Sample_Donor_Primary)){
  print(i)
  sample<-Sample_Donor_Primary[i,]
  sample_id<-sample$PPCG_Sample_ID
  sample_CNAs<-PPCG_CNA_Primary_CNs[PPCG_CNA_Primary_CNs$Sample == sample_id,]
  sample_AMP_CNAs<-sample_CNAs[sample_CNAs$relative_CN > 0,]
  sample_DEL_CNAs<-sample_CNAs[sample_CNAs$relative_CN < 0,]
  sample_AMP_CNAs<-sample_AMP_CNAs[!is.na(sample_AMP_CNAs$Donor),]
  sample_DEL_CNAs<-sample_DEL_CNAs[!is.na(sample_DEL_CNAs$Donor),]
  num_AMP_CNA<-nrow(sample_AMP_CNAs)
  num_DEL_CNA<-nrow(sample_DEL_CNAs)
  Sample_Donor_Primary[i,]$num_AMP_CNA<-num_AMP_CNA
  Sample_Donor_Primary[i,]$num_DEL_CNA<-num_DEL_CNA
  Sample_Donor_Primary[i,]$num_all_CNA<-num_DEL_CNA + num_AMP_CNA
  
}

for(i in 1:nrow(Sample_Donor_Metastasis)){
  print(i)
  sample<-Sample_Donor_Metastasis[i,]
  sample_id<-sample$PPCG_Sample_ID
  sample_CNAs<-PPCG_CNA_Metastasis_CNs[PPCG_CNA_Metastasis_CNs$Sample == sample_id,]
  sample_AMP_CNAs<-sample_CNAs[sample_CNAs$relative_CN > 0,]
  sample_DEL_CNAs<-sample_CNAs[sample_CNAs$relative_CN < 0,]
  sample_AMP_CNAs<-sample_AMP_CNAs[!is.na(sample_AMP_CNAs$Donor),]
  sample_DEL_CNAs<-sample_DEL_CNAs[!is.na(sample_DEL_CNAs$Donor),]
  num_AMP_CNA<-nrow(sample_AMP_CNAs)
  num_DEL_CNA<-nrow(sample_DEL_CNAs)
  Sample_Donor_Metastasis[i,]$num_AMP_CNA<-num_AMP_CNA
  Sample_Donor_Metastasis[i,]$num_DEL_CNA<-num_DEL_CNA
  Sample_Donor_Metastasis[i,]$num_all_CNA<-num_DEL_CNA + num_AMP_CNA
  
}

#Choose One sample per donor
Sample_Donor_Primary_filtered<-Sample_Donor_Primary[Sample_Donor_Primary$PPCG_Sample_ID %in% Sample_Donor_filtered$PPCG_Sample_ID,]
Sample_Donor_Metastasis_filtered<-Sample_Donor_Metastasis[Sample_Donor_Metastasis$PPCG_Sample_ID %in% Sample_Donor_filtered$PPCG_Sample_ID,]

PPCG_Sample_Primary_sub<-c()
for(i in 1:length(PPCG_Donor_ID_Primary)){
  print(i)
  Donor_ID<-PPCG_Donor_ID_Primary[i]
  samples<-Sample_Donor_Primary_filtered[Sample_Donor_Primary_filtered$PPCG_Donor_ID == Donor_ID,]
  if(nrow(samples) == 1){
    PPCG_Sample_Primary_sub<-c(PPCG_Sample_Primary_sub, samples$PPCG_Sample_ID)
  }else{
    Sample_Donor_Primary_sub<-Sample_Donor_Primary_filtered[Sample_Donor_Primary_filtered$PPCG_Sample_ID %in% samples$PPCG_Sample_ID,]
    Sample_Donor_Primary_sub<-Sample_Donor_Primary_sub[order(-Sample_Donor_Primary_sub$num_all_CNA),]
    PPCG_Sample_Primary_sub<-c(PPCG_Sample_Primary_sub, Sample_Donor_Primary_sub[1,]$PPCG_Sample_ID)
  }
}

PPCG_Sample_Metastasis_sub<-c()
for(i in 1:length(PPCG_Donor_ID_Metastasis)){
  print(i)
  Donor_ID<-PPCG_Donor_ID_Metastasis[i]
  samples<-Sample_Donor_Metastasis_filtered[Sample_Donor_Metastasis_filtered$PPCG_Donor_ID == Donor_ID,]
  if(nrow(samples) == 1){
    PPCG_Sample_Metastasis_sub<-c(PPCG_Sample_Metastasis_sub, samples$PPCG_Sample_ID)
  }else{
    Sample_Donor_Metastasis_sub<-Sample_Donor_Metastasis_filtered[Sample_Donor_Metastasis_filtered$PPCG_Sample_ID %in% samples$PPCG_Sample_ID,]
    Sample_Donor_Metastasis_sub<-Sample_Donor_Metastasis_sub[order(-Sample_Donor_Metastasis_sub$num_all_CNA),]
    PPCG_Sample_Metastasis_sub<-c(PPCG_Sample_Metastasis_sub, Sample_Donor_Metastasis_sub[1,]$PPCG_Sample_ID)
  }
}

# CNA Size Distribution
PPCG_CNA_Primary_CNs_sub<-PPCG_CNA_Primary_CNs[PPCG_CNA_Primary_CNs$Sample %in% PPCG_Sample_Primary_sub,]
PPCG_CNA_Primary_CNs_AMP<-PPCG_CNA_Primary_CNs_sub[PPCG_CNA_Primary_CNs_sub$relative_CN > 0,]
PPCG_CNA_Primary_CNs_AMP$type<-'AMP'
PPCG_CNA_Primary_CNs_DEL<-PPCG_CNA_Primary_CNs_sub[PPCG_CNA_Primary_CNs_sub$relative_CN < 0,]
PPCG_CNA_Primary_CNs_DEL$type<-'DEL'

PPCG_CNA_Metastasis_CNs_sub<-PPCG_CNA_Metastasis_CNs[PPCG_CNA_Metastasis_CNs$Sample %in% PPCG_Sample_Metastasis_sub,]
PPCG_CNA_Metastasis_CNs_AMP<-PPCG_CNA_Metastasis_CNs_sub[PPCG_CNA_Metastasis_CNs_sub$relative_CN > 0,]
PPCG_CNA_Metastasis_CNs_AMP$type<-'AMP'
PPCG_CNA_Metastasis_CNs_DEL<-PPCG_CNA_Metastasis_CNs_sub[PPCG_CNA_Metastasis_CNs_sub$relative_CN < 0,]
PPCG_CNA_Metastasis_CNs_DEL$type<-'DEL'

PPCG_CNA_Primary_CNs_AMP_table<-data.frame(table(PPCG_CNA_Primary_CNs_AMP$Donor))
PPCG_CNA_Primary_CNs_AMP_table$type<-'AMP'
PPCG_CNA_Primary_CNs_DEL_table<-data.frame(table(PPCG_CNA_Primary_CNs_DEL$Donor))
PPCG_CNA_Primary_CNs_DEL_table$type<-'DEL'

PPCG_CNA_Metastasis_CNs_AMP_table<-data.frame(table(PPCG_CNA_Metastasis_CNs_AMP$Donor))
PPCG_CNA_Metastasis_CNs_AMP_table$type<-'AMP'
PPCG_CNA_Metastasis_CNs_DEL_table<-data.frame(table(PPCG_CNA_Metastasis_CNs_DEL$Donor))
PPCG_CNA_Metastasis_CNs_DEL_table$type<-'DEL'

#Primary AMP
PPCG_CNA_Primary_CNs_AMP$width<-PPCG_CNA_Primary_CNs_AMP$End.bp - PPCG_CNA_Primary_CNs_AMP$Start.bp
PPCG_CNA_Primary_CNs_AMP$log_width<-log(PPCG_CNA_Primary_CNs_AMP$width)
p<-ggplot(PPCG_CNA_Primary_CNs_AMP, aes(x=as.numeric(as.character(log_width)))) + 
  geom_density(fill="#ff7f2a")+
  ylab('Density')+
  xlab('log CNA width')+
  geom_vline(aes(xintercept=11.51293),color="#cccccc", linetype="dashed", size=0.5)+
  geom_vline(aes(xintercept=13.81551),color="#cccccc", linetype="dashed", size=0.5)+
  geom_vline(aes(xintercept=16.1181),color="#cccccc", linetype="dashed", size=0.5)+
  scale_x_continuous(breaks=c(11.51293,13.81551,16.1181),labels=c("100kb", "1Mb", "10Mb"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.background = element_blank(), axis.line = element_line(color = "#2d3f44"))

pdf(paste(c(figure_dir, '003_Primary_AMP_Size_Distribution.pdf'), collapse=''), width=3, height= 2)
print(p)
dev.off()

#Primary DEL
PPCG_CNA_Primary_CNs_DEL$width<-PPCG_CNA_Primary_CNs_DEL$End.bp - PPCG_CNA_Primary_CNs_DEL$Start.bp
PPCG_CNA_Primary_CNs_DEL$log_width<-log(PPCG_CNA_Primary_CNs_DEL$width)
p<-ggplot(PPCG_CNA_Primary_CNs_DEL, aes(x=as.numeric(as.character(log_width)))) + 
  geom_density(fill="#3771c8")+
  ylab('Density')+
  xlab('log CNA width')+
  geom_vline(aes(xintercept=11.51293),color="#cccccc", linetype="dashed", size=0.5)+
  geom_vline(aes(xintercept=13.81551),color="#cccccc", linetype="dashed", size=0.5)+
  geom_vline(aes(xintercept=16.1181),color="#cccccc", linetype="dashed", size=0.5)+
  scale_x_continuous(breaks=c(11.51293,13.81551,16.1181),labels=c("100kb", "1Mb", "10Mb"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.background = element_blank(), axis.line = element_line(color = "#2d3f44"))
pdf(paste(c(figure_dir, '004_Primary_DEL_Size_Distribution.pdf'), collapse=''), width=3, height= 2)
print(p)
dev.off()

#Metastasis AMP
PPCG_CNA_Metastasis_CNs_AMP$width<-PPCG_CNA_Metastasis_CNs_AMP$End.bp - PPCG_CNA_Metastasis_CNs_AMP$Start.bp
PPCG_CNA_Metastasis_CNs_AMP$log_width<-log(PPCG_CNA_Metastasis_CNs_AMP$width)
p<-ggplot(PPCG_CNA_Metastasis_CNs_AMP, aes(x=as.numeric(as.character(log_width)))) + 
  geom_density(fill="#ff7f2a")+
  ylab('Density')+
  xlab('log CNA width')+
  geom_vline(aes(xintercept=11.51293),color="#cccccc", linetype="dashed", size=0.5)+
  geom_vline(aes(xintercept=13.81551),color="#cccccc", linetype="dashed", size=0.5)+
  geom_vline(aes(xintercept=16.1181),color="#cccccc", linetype="dashed", size=0.5)+
  scale_x_continuous(breaks=c(11.51293,13.81551,16.1181),labels=c("100kb", "1Mb", "10Mb"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.background = element_blank(), axis.line = element_line(color = "#2d3f44"))

pdf(paste(c(figure_dir, '005_Metastasis_AMP_Size_Distribution.pdf'), collapse=''), width=3, height= 2)
print(p)
dev.off()

#Metastasis DEL
PPCG_CNA_Metastasis_CNs_DEL$width<-PPCG_CNA_Metastasis_CNs_DEL$End.bp - PPCG_CNA_Metastasis_CNs_DEL$Start.bp
PPCG_CNA_Metastasis_CNs_DEL$log_width<-log(PPCG_CNA_Metastasis_CNs_DEL$width)
p<-ggplot(PPCG_CNA_Metastasis_CNs_DEL, aes(x=as.numeric(as.character(log_width)))) + 
  geom_density(fill="#3771c8")+
  ylab('Density')+
  xlab('log CNA width')+
  geom_vline(aes(xintercept=11.51293),color="#cccccc", linetype="dashed", size=0.5)+
  geom_vline(aes(xintercept=13.81551),color="#cccccc", linetype="dashed", size=0.5)+
  geom_vline(aes(xintercept=16.1181),color="#cccccc", linetype="dashed", size=0.5)+
  scale_x_continuous(breaks=c(11.51293,13.81551,16.1181),labels=c("100kb", "1Mb", "10Mb"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.background = element_blank(), axis.line = element_line(color = "#2d3f44"))
pdf(paste(c(figure_dir, '006_Metastasis_DEL_Size_Distribution.pdf'), collapse=''), width=3, height= 2)
print(p)
dev.off()

#Frequcy of CNAs
#Primary
PPCG_CNA_Primary_CNs_all<-rbind(PPCG_CNA_Primary_CNs_AMP,PPCG_CNA_Primary_CNs_DEL)
PPCG_CNA_Primary_CNs_all_table<-data.frame(table(PPCG_CNA_Primary_CNs_all$Donor))
PPCG_CNA_Primary_CNs_all_table<-PPCG_CNA_Primary_CNs_all_table[order(-PPCG_CNA_Primary_CNs_all_table$Freq),]
PPCG_CNA_Primary_CNs_all_table$Var1<-factor(PPCG_CNA_Primary_CNs_all_table$Var1)
PPCG_CNA_Primary_CNs_table<-rbind(PPCG_CNA_Primary_CNs_AMP_table,PPCG_CNA_Primary_CNs_DEL_table)
PPCG_CNA_Primary_CNs_table$Var1<-factor(PPCG_CNA_Primary_CNs_table$Var1, levels=PPCG_CNA_Primary_CNs_all_table$Var1)

p<-ggplot(PPCG_CNA_Primary_CNs_table, aes(x=Var1, y=Freq, fill=type)) +
  geom_bar(stat="identity") +
  ylab('Number of CNAs')+
  xlab('Donor')+
  scale_y_continuous(limits = c(0, 640), breaks = (seq(0,650,by = 50)))+
  scale_fill_manual(values=c("#ff7f2a", "#3771c8"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1),
        panel.background = element_blank(), axis.line = element_line(color = "#2d3f44"))

pdf(paste(c(figure_dir, '007_Primary_CNA_Frequency.pdf'), collapse=''), width=100, height= 10)
print(p)
dev.off()

#Metastasis
PPCG_CNA_Metastasis_CNs_all<-rbind(PPCG_CNA_Metastasis_CNs_AMP,PPCG_CNA_Metastasis_CNs_DEL)
PPCG_CNA_Metastasis_CNs_all_table<-data.frame(table(PPCG_CNA_Metastasis_CNs_all$Donor))
PPCG_CNA_Metastasis_CNs_all_table<-PPCG_CNA_Metastasis_CNs_all_table[order(-PPCG_CNA_Metastasis_CNs_all_table$Freq),]
PPCG_CNA_Metastasis_CNs_all_table$Var1<-factor(PPCG_CNA_Metastasis_CNs_all_table$Var1)
PPCG_CNA_Metastasis_CNs_table<-rbind(PPCG_CNA_Metastasis_CNs_AMP_table,PPCG_CNA_Metastasis_CNs_DEL_table)
PPCG_CNA_Metastasis_CNs_table$Var1<-factor(PPCG_CNA_Metastasis_CNs_table$Var1, levels=PPCG_CNA_Metastasis_CNs_all_table$Var1)

p<-ggplot(PPCG_CNA_Metastasis_CNs_table, aes(x=Var1, y=Freq, fill=type)) +
  geom_bar(stat="identity") +
  ylab('Number of CNAs')+
  xlab('Donor')+
  scale_y_continuous(limits = c(0, 640), breaks = (seq(0,650,by = 50)))+
  scale_fill_manual(values=c("#ff7f2a", "#3771c8"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1),
        panel.background = element_blank(), axis.line = element_line(color = "#2d3f44"))

pdf(paste(c(figure_dir, '008_Metastasis_CNA_Frequency.pdf'), collapse=''), width=15, height= 10)
print(p)
dev.off()

#GISTIC Input
PPCG_CNA_Metastasis_GISTIC<-PPCG_CNA_Metastasis[,c(2:6,8)]
colnames(PPCG_CNA_Metastasis_GISTIC)<-c("Sample", "Chromosome", "Start.bp", "End.bp", "Num.Markers", "Seg.CN")
PPCG_CNA_Primary_GISTIC<-PPCG_CNA_Primary[,c(2:6,8)]
colnames(PPCG_CNA_Primary_GISTIC)<-c("Sample", "Chromosome", "Start.bp", "End.bp", "Num.Markers", "Seg.CN")

write.table(PPCG_CNA_Metastasis_GISTIC, paste(c(data_dir, 'PPCG_CNA_Metastasis_GISTIC_2.tsv'), collapse=''), row.names=F, quote=F, sep='\t')
write.table(PPCG_CNA_Primary_GISTIC, paste(c(data_dir, 'PPCG_CNA_Primary_GISTIC_2.tsv'), collapse=''), row.names=F, quote=F, sep='\t')




