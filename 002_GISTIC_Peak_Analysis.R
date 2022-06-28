
library(Gviz)

gistic_dir<-'/.mounts/labs/reimandlab/private/projects/PPCG_CNA_pipeline/data/'
data_dir<-'/.mounts/labs/reimandlab/private/projects/PPCG_CNA_pipeline/data/'
figure_dir<-'/.mounts/labs/reimandlab/private/projects/PPCG_CNA_pipeline/figures/'

#Run GISTIC using BASH code in GISTIC directory /.mounts/labs/reimandlab/private/projects/PPCG_CNA_GISTIC/
#basedir=`pwd`/PPCG_results_Metastasis
#segfile=`pwd`/PPCG_GISTIC_files/PPCG_CNA_Metastasis_GISTIC.tsv
#refgenefile=`pwd`/refgenefiles/hg19.mat

#./gistic2 -b $basedir -seg $segfile -refgene $refgenefile -genegistic 1 -smallmem 1 -broad 1 -brlen 0.5 -conf 0.90 -armpeel 1 -savegene 1 -gcm extreme

#basedir=`pwd`/PPCG_results_Primary
#segfile=`pwd`/PPCG_GISTIC_files/PPCG_CNA_Primary_GISTIC.tsv
#refgenefile=`pwd`/refgenefiles/hg19.mat

#./gistic2 -b $basedir -seg $segfile -refgene $refgenefile -genegistic 1 -smallmem 1 -broad 1 -brlen 0.5 -conf 0.90 -armpeel 1 -savegene 1 -gcm extreme


table<-data.frame(table(Sample_Donor_filtered$Tissue_Origin))

#Get GISTIC Results
Metastasis_GISTIC<-read.delim(paste(c(gistic_dir, 'PPCG_results_Metastasis/all_lesions.conf_90.txt'), collapse=''))
Primary_GISTIC<-read.delim(paste(c(gistic_dir, 'PPCG_results_Primary/all_lesions.conf_90.txt'), collapse=''))

#Create GISTIC Peak GRanges
#Seperate GISTIC Output into AMP and DEL
Primary_GISTIC_AMP<-Primary_GISTIC[1:25,]
Primary_GISTIC_DEL<-Primary_GISTIC[26:57,]

Metastasis_GISTIC_AMP<-Metastasis_GISTIC[1:34,]
Metastasis_GISTIC_DEL<-Metastasis_GISTIC[35:70,]

Primary_GISTIC_AMP<-Primary_GISTIC_AMP[,1:3]
Primary_GISTIC_DEL<-Primary_GISTIC_DEL[,1:3]

Metastasis_GISTIC_AMP<-Metastasis_GISTIC_AMP[,1:3]
Metastasis_GISTIC_DEL<-Metastasis_GISTIC_DEL[,1:3]

Primary_GISTIC_AMP$chr<-NA
Primary_GISTIC_AMP$start<-NA
Primary_GISTIC_AMP$end<-NA
Primary_GISTIC_AMP$width<-NA

for(i in 1:nrow(Primary_GISTIC_AMP)){
  print(i)
  peak<-Primary_GISTIC_AMP[i,]
  list<-unlist(strsplit(peak$Wide.Peak.Limits, "[(:-]"))
  chr<-list[1]
  start<-as.numeric(as.character(list[2]))
  end<-as.numeric(as.character(list[3]))
  Primary_GISTIC_AMP[i,]$chr<-chr
  Primary_GISTIC_AMP[i,]$start<-start
  Primary_GISTIC_AMP[i,]$end<-end
  Primary_GISTIC_AMP[i,]$width<-end-start
}

Primary_GISTIC_DEL$chr<-NA
Primary_GISTIC_DEL$start<-NA
Primary_GISTIC_DEL$end<-NA
Primary_GISTIC_DEL$width<-NA

for(i in 1:nrow(Primary_GISTIC_DEL)){
  print(i)
  peak<-Primary_GISTIC_DEL[i,]
  list<-unlist(strsplit(peak$Wide.Peak.Limits, "[(:-]"))
  chr<-list[1]
  start<-as.numeric(as.character(list[2]))
  end<-as.numeric(as.character(list[3]))
  Primary_GISTIC_DEL[i,]$chr<-chr
  Primary_GISTIC_DEL[i,]$start<-start
  Primary_GISTIC_DEL[i,]$end<-end
  Primary_GISTIC_DEL[i,]$width<-end-start
}

Metastasis_GISTIC_AMP$chr<-NA
Metastasis_GISTIC_AMP$start<-NA
Metastasis_GISTIC_AMP$end<-NA
Metastasis_GISTIC_AMP$width<-NA

for(i in 1:nrow(Metastasis_GISTIC_AMP)){
  print(i)
  peak<-Metastasis_GISTIC_AMP[i,]
  list<-unlist(strsplit(peak$Wide.Peak.Limits, "[(:-]"))
  chr<-list[1]
  start<-as.numeric(as.character(list[2]))
  end<-as.numeric(as.character(list[3]))
  Metastasis_GISTIC_AMP[i,]$chr<-chr
  Metastasis_GISTIC_AMP[i,]$start<-start
  Metastasis_GISTIC_AMP[i,]$end<-end
  Metastasis_GISTIC_AMP[i,]$width<-end-start
}

Metastasis_GISTIC_DEL$chr<-NA
Metastasis_GISTIC_DEL$start<-NA
Metastasis_GISTIC_DEL$end<-NA
Metastasis_GISTIC_DEL$width<-NA

for(i in 1:nrow(Metastasis_GISTIC_DEL)){
  print(i)
  peak<-Metastasis_GISTIC_DEL[i,]
  list<-unlist(strsplit(peak$Wide.Peak.Limits, "[(:-]"))
  chr<-list[1]
  start<-as.numeric(as.character(list[2]))
  end<-as.numeric(as.character(list[3]))
  Metastasis_GISTIC_DEL[i,]$chr<-chr
  Metastasis_GISTIC_DEL[i,]$start<-start
  Metastasis_GISTIC_DEL[i,]$end<-end
  Metastasis_GISTIC_DEL[i,]$width<-end-start
}

#Create GISTIC Peak GRanges
PPCG_CNA_Primary_CNs_AMP<-PPCG_CNA_Primary_CNs_AMP[!is.na(PPCG_CNA_Primary_CNs_AMP$Donor),]
PPCG_CNA_Primary_CNs_AMP_grange<-GRanges(Rle(PPCG_CNA_Primary_CNs_AMP$Chromosome), ranges = IRanges(start = PPCG_CNA_Primary_CNs_AMP$Start.bp, end = PPCG_CNA_Primary_CNs_AMP$End.bp))
PPCG_CNA_Primary_CNs_AMP_grange$sample_id<-PPCG_CNA_Primary_CNs_AMP$Sample
PPCG_CNA_Primary_CNs_DEL<-PPCG_CNA_Primary_CNs_DEL[!is.na(PPCG_CNA_Primary_CNs_DEL$Donor),]
PPCG_CNA_Primary_CNs_DEL_grange<-GRanges(Rle(PPCG_CNA_Primary_CNs_DEL$Chromosome), ranges = IRanges(start = PPCG_CNA_Primary_CNs_DEL$Start.bp, end = PPCG_CNA_Primary_CNs_DEL$End.bp))
PPCG_CNA_Primary_CNs_DEL_grange$sample_id<-PPCG_CNA_Primary_CNs_DEL$Sample

PPCG_CNA_Metastasis_CNs_AMP<-PPCG_CNA_Metastasis_CNs_AMP[!is.na(PPCG_CNA_Metastasis_CNs_AMP$Donor),]
PPCG_CNA_Metastasis_CNs_AMP_grange<-GRanges(Rle(PPCG_CNA_Metastasis_CNs_AMP$Chromosome), ranges = IRanges(start = PPCG_CNA_Metastasis_CNs_AMP$Start.bp, end = PPCG_CNA_Metastasis_CNs_AMP$End.bp))
PPCG_CNA_Metastasis_CNs_AMP_grange$sample_id<-PPCG_CNA_Metastasis_CNs_AMP$Sample
PPCG_CNA_Metastasis_CNs_DEL<-PPCG_CNA_Metastasis_CNs_DEL[!is.na(PPCG_CNA_Metastasis_CNs_DEL$Donor),]
PPCG_CNA_Metastasis_CNs_DEL_grange<-GRanges(Rle(PPCG_CNA_Metastasis_CNs_DEL$Chromosome), ranges = IRanges(start = PPCG_CNA_Metastasis_CNs_DEL$Start.bp, end = PPCG_CNA_Metastasis_CNs_DEL$End.bp))
PPCG_CNA_Metastasis_CNs_DEL_grange$sample_id<-PPCG_CNA_Metastasis_CNs_DEL$Sample

#Create GISTIC Peak Figures
#load data for figures
all_gene_ids<-readRDS(paste(c(data_dir, 'all_gene_ids.rds'), collapse=''))
all_gene_ids_cgc<-all_gene_ids[all_gene_ids$cgc == TRUE,]

PRAD_peakCalls<-read.delim(paste(c(data_dir, 'PRAD_peakCalls.txt'), collapse=''))
PRAD_ATAC_seq_grange<-GRanges(Rle(PRAD_peakCalls$seqnames), ranges = IRanges(start = PRAD_peakCalls$start, end = PRAD_peakCalls$end))
PRAD_ATAC_seq_grange$score<-PRAD_peakCalls$score

#Primary AMP GISTIC Peaks
for(i in 1:nrow(Primary_GISTIC_AMP)){
  print(i)
  gistic_peak<-Primary_GISTIC_AMP[i,]
  gistic_chr<-gistic_peak$chr
  gistic_start<-as.numeric(as.character(gistic_peak$start))
  gistic_end<-as.numeric(as.character(gistic_peak$end))
  band<-gistic_peak$Descriptor
  
  gistic_grange<-GRanges(Rle(gistic_chr), ranges = IRanges(start = gistic_start, end = gistic_end))
  
  #Figure Width
  start<-as.numeric(as.character(gistic_start)) - 3000000
  end<-as.numeric(as.character(gistic_end)) + 2000000
  chr<-as.character(gistic_chr)
  grange<-GRanges(Rle(chr), ranges = IRanges(start = start, end = end))
  
  #Chromosome Figure and Gene Tracks 
  gtrack <- GenomeAxisTrack()
  gistic_track<-AnnotationTrack(gistic_grange, genome="hg19", name="gistic", fill="#ffcc00")
  biomartTrack2 <- BiomartGeneRegionTrack(rotation.title=0, genome="hg19", chromosome=chr, start=start, end=end, stacking="squish", collapseTranscripts = "meta", name="ENSEMBL", fill="#910000", filters=list(biotype="protein_coding", ensembl_gene_id = as.character(all_gene_ids_cgc$ensembl_gene_id)))
  biomartTrack1 <- BiomartGeneRegionTrack(rotation.title=0,genome="hg19", chromosome=chr, start=start, end=end, stacking="squish", collapseTranscripts = "meta", name="ENSEMBL", fill="#ebb52e", filters=list(biotype="protein_coding"))
  biomartTrack1 <- biomartTrack1[!(gene(biomartTrack1) %in% all_gene_ids_cgc$ensembl_gene_id  )]
  itrack <- IdeogramTrack(genome="hg19", chromosome=chr)
  
  #TCGA PRAD ATAC-Seq
  x<-as.data.frame(findOverlaps(grange, PRAD_ATAC_seq_grange))
  ATAC_seq_grange_sub<-PRAD_ATAC_seq_grange[x$subjectHits]
  ATAC_start<-as.numeric(as.character(start(ATAC_seq_grange_sub)))
  ATAC_end<-as.numeric(as.character(end(ATAC_seq_grange_sub)))
  ATAC_dat<-as.numeric(as.character(ATAC_seq_grange_sub$score))
  ATAC_track <- DataTrack(data = ATAC_dat, start = ATAC_start,end = ATAC_end, chromosome = chr, genome ="hg19",name = "ATAC_peak", type="histogram")
  tracks<-c(itrack, gtrack, biomartTrack1, biomartTrack2, ATAC_track, gistic_track)
  
  #CNAs Overlaping GISTIC Peak
  x<-as.data.frame(findOverlaps(grange, PPCG_CNA_Primary_CNs_AMP_grange))
  PPCG_CNA_sub_df<-data.frame(PPCG_CNA_Primary_CNs_AMP_grange[x$subjectHits])
  PPCG_CNA_sub_df<-PPCG_CNA_sub_df[order(-PPCG_CNA_sub_df$width),]
  
  all_CNA_patients<-unique(PPCG_CNA_sub_df$sample_id)
  all_CNA_patients<-rev(all_CNA_patients)
  
  x<-as.data.frame(findOverlaps(gistic_grange, PPCG_CNA_Primary_CNs_AMP_grange))
  gistic_CNA_patients<-unique(all_CNA_patients[all_CNA_patients %in% unique(PPCG_CNA_Primary_CNs_AMP_grange[x$subjectHits]$sample_id)])
  
  if(length(gistic_CNA_patients) > 0){
    for(j in 1:length(gistic_CNA_patients)){
      sample_id<-gistic_CNA_patients[j]
      patient_CNA<-PPCG_CNA_Primary_CNs_AMP_grange[PPCG_CNA_Primary_CNs_AMP_grange$sample_id == sample_id]
      x<-as.data.frame(findOverlaps(grange, patient_CNA))
      CNA<-patient_CNA[x$subjectHits]
      CNA<-reduce(CNA)
      track <- AnnotationTrack(CNA, name=sample_id, rotation.title=0, genome="hg19", background.title ="grey", fontcolor.group="black", fill="#ff6600",showTitle=TRUE)
      tracks<-c(tracks, track)
    }    
  }
  
  #CNAs Not Overlaping GISTIC Peak
  other_CNA_patients<-unique(all_CNA_patients[!(all_CNA_patients %in% gistic_CNA_patients)])
  
  if(length(other_CNA_patients) > 0){
    for(j in 1:length(other_CNA_patients)){
      sample_id<-other_CNA_patients[j]
      patient_CNA<-PPCG_CNA_Primary_CNs_AMP_grange[PPCG_CNA_Primary_CNs_AMP_grange$sample_id == sample_id]
      x<-as.data.frame(findOverlaps(grange, patient_CNA))
      CNA<-patient_CNA[x$subjectHits]
      CNA<-reduce(CNA)
      track <- AnnotationTrack(CNA, name=sample_id, rotation.title=0, genome="hg19", background.title ="grey", fontcolor.group="black", fill="#ffb380",showTitle=TRUE)
      tracks<-c(tracks, track)
    }    
  }
  
  pdf(paste(c(figure_dir, 'GISTIC_Peaks/Primary_AMP/', band, "_", chr, "_", start, "_", end, ".pdf"), collapse=''), width=10, height= 30)
  plotTracks(tracks, from=start, to=end, title.width=4, cex.title=0.3, transcriptAnnotation = "symbol", main=paste(c(chr, ':', start,'-', end), collapse=' '), cex.main=1)
  dev.off()
  
}

#Primary DEL GISTIC Peaks
for(i in 1:nrow(Primary_GISTIC_DEL)){
  print(i)
  gistic_peak<-Primary_GISTIC_DEL[i,]
  gistic_chr<-gistic_peak$chr
  gistic_start<-as.numeric(as.character(gistic_peak$start))
  gistic_end<-as.numeric(as.character(gistic_peak$end))
  band<-gistic_peak$Descriptor
  
  gistic_grange<-GRanges(Rle(gistic_chr), ranges = IRanges(start = gistic_start, end = gistic_end))
  
  start<-as.numeric(as.character(gistic_start)) - 5000000
  end<-as.numeric(as.character(gistic_end)) + 5000000
  chr<-as.character(gistic_chr)
  grange<-GRanges(Rle(chr), ranges = IRanges(start = start, end = end))
  
  gtrack <- GenomeAxisTrack()
  gistic_track<-AnnotationTrack(gistic_grange, genome="hg19", name="gistic", fill="#ffcc00")
  biomartTrack2 <- BiomartGeneRegionTrack(rotation.title=0, genome="hg19", chromosome=chr, start=start, end=end, stacking="squish", collapseTranscripts = "meta", name="ENSEMBL", fill="#910000", filters=list(biotype="protein_coding", ensembl_gene_id = as.character( all_gene_ids_cgc$ensembl_gene_id)))
  biomartTrack1 <- BiomartGeneRegionTrack(rotation.title=0,genome="hg19", chromosome=chr, start=start, end=end, stacking="squish", collapseTranscripts = "meta", name="ENSEMBL", fill="#ebb52e", filters=list(biotype="protein_coding"))
  biomartTrack1 <- biomartTrack1[!(gene(biomartTrack1) %in% all_gene_ids_cgc$ensembl_gene_id  )]
  itrack <- IdeogramTrack(genome="hg19", chromosome=chr)

  x<-as.data.frame(findOverlaps(grange, PRAD_ATAC_seq_grange))
  ATAC_seq_grange_sub<-PRAD_ATAC_seq_grange[x$subjectHits]
  ATAC_start<-as.numeric(as.character(start(ATAC_seq_grange_sub)))
  ATAC_end<-as.numeric(as.character(end(ATAC_seq_grange_sub)))
  ATAC_dat<-as.numeric(as.character(ATAC_seq_grange_sub$score))
  ATAC_track <- DataTrack(data = ATAC_dat, start = ATAC_start,end = ATAC_end, chromosome = chr, genome ="hg19",name = "ATAC_peak", type="histogram")
  tracks<-c(itrack, gtrack, biomartTrack1, biomartTrack2, ATAC_track, gistic_track)
  
  x<-as.data.frame(findOverlaps(grange, PPCG_CNA_Primary_CNs_DEL_grange))
  PPCG_CNA_sub_df<-data.frame(PPCG_CNA_Primary_CNs_DEL_grange[x$subjectHits])
  PPCG_CNA_sub_df<-PPCG_CNA_sub_df[order(-PPCG_CNA_sub_df$width),]
  
  all_CNA_patients<-unique(PPCG_CNA_sub_df$sample_id)
  all_CNA_patients<-rev(all_CNA_patients)
  
  x<-as.data.frame(findOverlaps(gistic_grange, PPCG_CNA_Primary_CNs_DEL_grange))
  gistic_CNA_patients<-unique(all_CNA_patients[all_CNA_patients %in% unique(PPCG_CNA_Primary_CNs_DEL_grange[x$subjectHits]$sample_id)])
  
  if(length(gistic_CNA_patients) > 0){
    for(j in 1:length(gistic_CNA_patients)){
      sample_id<-gistic_CNA_patients[j]
      patient_CNA<-PPCG_CNA_Primary_CNs_DEL_grange[PPCG_CNA_Primary_CNs_DEL_grange$sample_id == sample_id]
      x<-as.data.frame(findOverlaps(grange, patient_CNA))
      CNA<-patient_CNA[x$subjectHits]
      CNA<-reduce(CNA)
      track <- AnnotationTrack(CNA, name=sample_id, rotation.title=0, genome="hg19", background.title ="grey", fontcolor.group="black", fill="#003380",showTitle=TRUE)
      tracks<-c(tracks, track)
    }    
  }
  
  other_CNA_patients<-unique(all_CNA_patients[!(all_CNA_patients %in% gistic_CNA_patients)])
  
  if(length(other_CNA_patients) > 0){
    for(j in 1:length(other_CNA_patients)){
      sample_id<-other_CNA_patients[j]
      patient_CNA<-PPCG_CNA_Primary_CNs_DEL_grange[PPCG_CNA_Primary_CNs_DEL_grange$sample_id == sample_id]
      x<-as.data.frame(findOverlaps(grange, patient_CNA))
      CNA<-patient_CNA[x$subjectHits]
      CNA<-reduce(CNA)
      track <- AnnotationTrack(CNA, name=sample_id, rotation.title=0, genome="hg19", background.title ="grey", fontcolor.group="black", fill="#5599ff",showTitle=TRUE)
      tracks<-c(tracks, track)
    }    
  }
  
  pdf(paste(c(figure_dir, 'GISTIC_Peaks/Primary_DEL/', band, "_", chr, "_", start, "_", end, ".pdf"), collapse=''), width=10, height= 30)
  plotTracks(tracks, from=start, to=end, title.width=4, cex.title=0.3, transcriptAnnotation = "symbol", main=paste(c(chr, ':', start,'-', end), collapse=' '), cex.main=1)
  dev.off()
  
}

#Metastasis AMP GISTIC Peaks
for(i in 1:nrow(Metastasis_GISTIC_AMP)){
  print(i)
  gistic_peak<-Metastasis_GISTIC_AMP[i,]
  gistic_chr<-gistic_peak$chr
  gistic_start<-as.numeric(as.character(gistic_peak$start))
  gistic_end<-as.numeric(as.character(gistic_peak$end))
  band<-gistic_peak$Descriptor
  
  gistic_grange<-GRanges(Rle(gistic_chr), ranges = IRanges(start = gistic_start, end = gistic_end))
  
  start<-as.numeric(as.character(gistic_start)) - 3000000
  end<-as.numeric(as.character(gistic_end)) + 2000000
  chr<-as.character(gistic_chr)
  grange<-GRanges(Rle(chr), ranges = IRanges(start = start, end = end))
  
  gtrack <- GenomeAxisTrack()
  gistic_track<-AnnotationTrack(gistic_grange, genome="hg19", name="gistic", fill="#ffcc00")
  biomartTrack2 <- BiomartGeneRegionTrack(rotation.title=0, genome="hg19", chromosome=chr, start=start, end=end, stacking="squish", collapseTranscripts = "meta", name="ENSEMBL", fill="#910000", filters=list(biotype="protein_coding", ensembl_gene_id = as.character(all_gene_ids_cgc$ensembl_gene_id)))
  biomartTrack1 <- BiomartGeneRegionTrack(rotation.title=0,genome="hg19", chromosome=chr, start=start, end=end, stacking="squish", collapseTranscripts = "meta", name="ENSEMBL", fill="#ebb52e", filters=list(biotype="protein_coding"))
  biomartTrack1 <- biomartTrack1[!(gene(biomartTrack1) %in% all_gene_ids_cgc$ensembl_gene_id  )]
  itrack <- IdeogramTrack(genome="hg19", chromosome=chr)
  
  x<-as.data.frame(findOverlaps(grange, PRAD_ATAC_seq_grange))
  ATAC_seq_grange_sub<-PRAD_ATAC_seq_grange[x$subjectHits]
  ATAC_start<-as.numeric(as.character(start(ATAC_seq_grange_sub)))
  ATAC_end<-as.numeric(as.character(end(ATAC_seq_grange_sub)))
  ATAC_dat<-as.numeric(as.character(ATAC_seq_grange_sub$score))
  ATAC_track <- DataTrack(data = ATAC_dat, start = ATAC_start,end = ATAC_end, chromosome = chr, genome ="hg19",name = "ATAC_peak", type="histogram")
  tracks<-c(itrack, gtrack, biomartTrack1, biomartTrack2, ATAC_track, gistic_track)

  x<-as.data.frame(findOverlaps(grange, PPCG_CNA_Metastasis_CNs_AMP_grange))
  PPCG_CNA_sub_df<-data.frame(PPCG_CNA_Metastasis_CNs_AMP_grange[x$subjectHits])
  PPCG_CNA_sub_df<-PPCG_CNA_sub_df[order(-PPCG_CNA_sub_df$width),]
  
  all_CNA_patients<-unique(PPCG_CNA_sub_df$sample_id)
  all_CNA_patients<-rev(all_CNA_patients)
  
  x<-as.data.frame(findOverlaps(gistic_grange, PPCG_CNA_Metastasis_CNs_AMP_grange))
  
  gistic_CNA_patients<-unique(all_CNA_patients[all_CNA_patients %in% unique(PPCG_CNA_Metastasis_CNs_AMP_grange[x$subjectHits]$sample_id)])
  
  if(length(gistic_CNA_patients) > 0){
    for(j in 1:length(gistic_CNA_patients)){
      sample_id<-gistic_CNA_patients[j]
      patient_CNA<-PPCG_CNA_Metastasis_CNs_AMP_grange[PPCG_CNA_Metastasis_CNs_AMP_grange$sample_id == sample_id]
      x<-as.data.frame(findOverlaps(grange, patient_CNA))
      CNA<-patient_CNA[x$subjectHits]
      CNA<-reduce(CNA)
      track <- AnnotationTrack(CNA, name=sample_id, rotation.title=0, genome="hg19", background.title ="grey", fontcolor.group="black", fill="#ff6600",showTitle=TRUE)
      tracks<-c(tracks, track)
    }    
  }
  
  other_CNA_patients<-unique(all_CNA_patients[!(all_CNA_patients %in% gistic_CNA_patients)])
  
  if(length(other_CNA_patients) > 0){
    for(j in 1:length(other_CNA_patients)){
      sample_id<-other_CNA_patients[j]
      patient_CNA<-PPCG_CNA_Metastasis_CNs_AMP_grange[PPCG_CNA_Metastasis_CNs_AMP_grange$sample_id == sample_id]
      x<-as.data.frame(findOverlaps(grange, patient_CNA))
      CNA<-patient_CNA[x$subjectHits]
      CNA<-reduce(CNA)
      track <- AnnotationTrack(CNA, name=sample_id, rotation.title=0, genome="hg19", background.title ="grey", fontcolor.group="black", fill="#ffb380",showTitle=TRUE)
      tracks<-c(tracks, track)
    }    
  }
  
  pdf(paste(c(figure_dir, 'GISTIC_Peaks/Metastasis_AMP/', band, "_", chr, "_", start, "_", end, ".pdf"), collapse=''), width=10, height= 15)
  plotTracks(tracks, from=start, to=end, title.width=4, cex.title=0.3, transcriptAnnotation = "symbol", main=paste(c(chr, ':', start,'-', end), collapse=' '), cex.main=1)
  dev.off()
  
}

#Metastasis DEL GISTIC Peaks
for(i in 1:nrow(Metastasis_GISTIC_DEL)){
  print(i)
  gistic_peak<-Metastasis_GISTIC_DEL[i,]
  gistic_chr<-gistic_peak$chr
  gistic_start<-as.numeric(as.character(gistic_peak$start))
  gistic_end<-as.numeric(as.character(gistic_peak$end))
  band<-gistic_peak$Descriptor
  
  gistic_grange<-GRanges(Rle(gistic_chr), ranges = IRanges(start = gistic_start, end = gistic_end))
  
  start<-as.numeric(as.character(gistic_start)) - 5000000
  end<-as.numeric(as.character(gistic_end)) + 5000000
  chr<-as.character(gistic_chr)
  grange<-GRanges(Rle(chr), ranges = IRanges(start = start, end = end))
  
  gtrack <- GenomeAxisTrack()
  gistic_track<-AnnotationTrack(gistic_grange, genome="hg19", name="gistic", fill="#ffcc00")
  biomartTrack2 <- BiomartGeneRegionTrack(rotation.title=0, genome="hg19", chromosome=chr, start=start, end=end, stacking="squish", collapseTranscripts = "meta", name="ENSEMBL", fill="#910000", filters=list(biotype="protein_coding", ensembl_gene_id = as.character( all_gene_ids_cgc$ensembl_gene_id)))
  biomartTrack1 <- BiomartGeneRegionTrack(rotation.title=0,genome="hg19", chromosome=chr, start=start, end=end, stacking="squish", collapseTranscripts = "meta", name="ENSEMBL", fill="#ebb52e", filters=list(biotype="protein_coding"))
  biomartTrack1 <- biomartTrack1[!(gene(biomartTrack1) %in% all_gene_ids_cgc$ensembl_gene_id  )]
  itrack <- IdeogramTrack(genome="hg19", chromosome=chr)

  x<-as.data.frame(findOverlaps(grange, PRAD_ATAC_seq_grange))
  ATAC_seq_grange_sub<-PRAD_ATAC_seq_grange[x$subjectHits]
  ATAC_start<-as.numeric(as.character(start(ATAC_seq_grange_sub)))
  ATAC_end<-as.numeric(as.character(end(ATAC_seq_grange_sub)))
  ATAC_dat<-as.numeric(as.character(ATAC_seq_grange_sub$score))
  ATAC_track <- DataTrack(data = ATAC_dat, start = ATAC_start,end = ATAC_end, chromosome = chr, genome ="hg19",name = "ATAC_peak", type="histogram")
  tracks<-c(itrack, gtrack, biomartTrack1, biomartTrack2, ATAC_track, gistic_track)
  
  x<-as.data.frame(findOverlaps(grange, PPCG_CNA_Metastasis_CNs_DEL_grange))
  PPCG_CNA_sub_df<-data.frame(PPCG_CNA_Metastasis_CNs_DEL_grange[x$subjectHits])
  PPCG_CNA_sub_df<-PPCG_CNA_sub_df[order(-PPCG_CNA_sub_df$width),]
  
  all_CNA_patients<-unique(PPCG_CNA_sub_df$sample_id)
  all_CNA_patients<-rev(all_CNA_patients)
  
  x<-as.data.frame(findOverlaps(gistic_grange, PPCG_CNA_Metastasis_CNs_DEL_grange))
  gistic_CNA_patients<-unique(all_CNA_patients[all_CNA_patients %in% unique(PPCG_CNA_Metastasis_CNs_DEL_grange[x$subjectHits]$sample_id)])
  
  if(length(gistic_CNA_patients) > 0){
    for(j in 1:length(gistic_CNA_patients)){
      sample_id<-gistic_CNA_patients[j]
      patient_CNA<-PPCG_CNA_Metastasis_CNs_DEL_grange[PPCG_CNA_Metastasis_CNs_DEL_grange$sample_id == sample_id]
      x<-as.data.frame(findOverlaps(grange, patient_CNA))
      CNA<-patient_CNA[x$subjectHits]
      CNA<-reduce(CNA)
      track <- AnnotationTrack(CNA, name=sample_id, rotation.title=0, genome="hg19", background.title ="grey", fontcolor.group="black", fill="#003380",showTitle=TRUE)
      tracks<-c(tracks, track)
    }    
  }
  
  other_CNA_patients<-unique(all_CNA_patients[!(all_CNA_patients %in% gistic_CNA_patients)])
  
  if(length(other_CNA_patients) > 0){
    for(j in 1:length(other_CNA_patients)){
      sample_id<-other_CNA_patients[j]
      patient_CNA<-PPCG_CNA_Metastasis_CNs_DEL_grange[PPCG_CNA_Metastasis_CNs_DEL_grange$sample_id == sample_id]
      x<-as.data.frame(findOverlaps(grange, patient_CNA))
      CNA<-patient_CNA[x$subjectHits]
      CNA<-reduce(CNA)
      track <- AnnotationTrack(CNA, name=sample_id, rotation.title=0, genome="hg19", background.title ="grey", fontcolor.group="black", fill="#5599ff",showTitle=TRUE)
      tracks<-c(tracks, track)
    }    
  }
  
  pdf(paste(c(figure_dir, 'GISTIC_Peaks/Metastasis_DEL/', band, "_", chr, "_", start, "_", end, ".pdf"), collapse=''), width=10, height= 15)
  plotTracks(tracks, from=start, to=end, title.width=4, cex.title=0.3, transcriptAnnotation = "symbol", main=paste(c(chr, ':', start,'-', end), collapse=' '), cex.main=1)
  dev.off()
  
}




