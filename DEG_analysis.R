
#########TCGA
##step1 download data
# https://xenabrowser.net/datapages/?dataset=TCGA-BLCA.htseq_counts.tsv&host=https%3A%2F%2Fgdc.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443
# https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-BLCA.GDC_phenotype.tsv.gz
# https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-BLCA.htseq_counts.tsv.gz
# https://gdc-hub.s3.us-east-1.amazonaws.com/download/gencode.v22.annotation.gene.probeMap

# 01cancer，11paracancerous

if (!file.exists("./Transcriptome/TCGA/TCGA-BLCA.Rdata")) {
  gzfile<- "./Transcriptome/TCGA/TCGA-BLCA.htseq_counts.tsv.gz"
  download.file("https://tcga.xenahubs.net/download/TCGA.BLCA.sampleMap/HiSeqV2.gz",
                destfile = gzfile) ##downalod 
  library(R.utils)
  gunzip(gzfile, remove=F)
  library(data.table)
  raw_data<- fread("./Transcriptome/TCGA/TCGA-BLCA.htseq_counts.tsv", 
                   header = T)
  raw_data<- as.data.frame(raw_data)
  raw_data[1:5,1:3]
  rownames(raw_data)<- raw_data[,1]
  raw_data<- raw_data[, -1]
  raw_data<- 2^raw_data - 1
  raw_data<- ceiling(raw_data) ##60488
  pick_raw<- apply(raw_data, 1, function(x) {
    sum(x==0)<10
  }) 

  raw_data<- raw_data[pick_raw,] ##17920
 
  a=read.table('./Transcriptome/TCGA/gencode.v22.annotation.gene.probeMap',header = T)
  head(a)
  ids=a[match(rownames(raw_data),a$id),1:2] #17920
  head(ids)
  tmp=table(ids$gene)
  head(ids[ids$gene %in% names(tmp[tmp==2]),])
  ids[ids$gene %in% head(names(tmp[tmp==2])),]
  colnames(ids)=c('probe_id','symbol')  
  ids=ids[ids$symbol != '',]
  dat=raw_data
  ids=ids[ids$probe_id %in%  rownames(dat),]

  dat=dat[ids$probe_id,] ##17917
  ids$median=apply(dat,1,median) 

  ids=ids[order(ids$symbol,ids$median,decreasing = T),]
  ids=ids[!duplicated(ids$symbol),]
  dat=dat[ids$probe_id,] 
  colnames(dat) <- colnames(raw_data)
  dat <- dat[,!(colnames(dat) %in% c('TCGA-BL-A0C8-01B', 'TCGA-BL-A13I-01B', 'TCGA-BL-A13J-01B'))]
  colname <- gsub("(.+\\-[[:digit:]]+)([[:alpha:]]+)",'\\1',colnames(dat))
  tmp <- table(colname)
  colnames(dat) <- colname
  save(dat, file='./Transcriptome/TCGA/TCGA-BLCA.Rdata')
} else {
  load('./Transcriptome/TCGA/TCGA-BLCA.Rdata')
}


##step2 grouping by special clinical information
if (!file.exists('./Transcriptome/TCGA/TCGA-BLCA_phenotype.Rdata')) {
  gzfile<- './Transcriptome/TCGA/TCGA-BLCA_phenotype.tsv.gz'
  download.file('https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-BLCA.GDC_phenotype.tsv.gz',
                destfile=gzfile)
  phenotype<- read.table(gzfile,
                         header =T,
                         sep = '\t',
                         quote = '\"',

  group <- gsub("(.+)-([[:digit:]])",'\\2', phenotype$sampleID)

   phenotype$group <- ifelse(group %in% "11", "normal", "tumor") 
 
  phenotype[249:260,c(1,140)]
  save(phenotype, file="./Transcriptome/TCGA/TCGA-BLCA_phenotype.Rdata")
} else{
  load("./Transcriptome/TCGA/TCGA-BLCA_phenotype.Rdata")
}


  

#####
  group_info <- phenotype[,c(1,140)] #436
  tmp <- table(group_info$sampleID)  #436

  group_info <-group_info[match(colname,group_info$sampleID),] #427
  group_info[1:3,]
  dat[1:3,group_info$sampleID[425:427]]
  colnames(group_info) <- c('sampleID', 'type')
  table(group_info$type)
  group_info$type <-factor(group_info$type,levels=c("normal","tumor"))


library(limma)
library(edgeR)
DGE_limma<- DGEList(counts = dat, group =group_info$type)

keep_gene<- rowSums(cpm(DGE_limma) >1) >= 2  
DGE_limma<- DGE_limma[keep_gene, , keep.lib.sizes=F]
DGE_limma<- calcNormFactors(DGE_limma)

design<- model.matrix(~0 + group_info$type) 
rownames(design)<- colnames(DGE_limma)
colnames(design)<- levels(group_info$type)


##transform RNA-seq data ready for linear modeling
v<- voom(DGE_limma, design , plot = T, normalize="quantile")
fit<- lmFit(v, design)
contrasts <- paste(rev(levels(group_info$type)), collapse = '-')
cont.martrix<- makeContrasts(contrasts = c("tumor-normal"), levels = design)
fit2<- contrasts.fit(fit, cont.martrix)
fit2<- eBayes(fit2) 
DGE_limma<- topTable(fit2, coef = contrasts, n=Inf)
DGE_limma<- na.omit(DGE_limma)
head(DGE_limma)


logFC_cutoff <- with(DGE_limma, mean(abs(logFC))+2*sd(abs(logFC)))

DGE_limma$change_LFCM2SD<- as.factor(ifelse(DGE_limma$P.Value < 0.05 & abs(DGE_limma$logFC
) > logFC_cutoff,  ifelse(DGE_limma$logFC > logFC_cutoff, "UP", "DOWN"), "NOT"))
save(raw_degs1, DGE_edgeR,DGE_limma,file="./Transcriptome/TCGA/TCGA_DEG_TumorVSnormal.RData")

# 1.4 
tj <- data.frame(DESeq2=as.integer(table(raw_degs1$change_LFCM2SD)),
                 edgeR=as.integer(table(DGE_edgeR$change_LFCM2SD)),
                 limma=as.integer(table(DGE_limma$change_LFCM2SD)),
                 row.names = c( 'DOWN',   'NOT',    'UP' ))

up <- intersect(intersect(rownames(raw_degs1)[raw_degs1$change_LFCM2SD=='UP'],rownames(DGE_edgeR)[DGE_edgeR$change_LFCM2SD=='UP']),
                rownames(DGE_limma)[DGE_limma$change_LFCM2SD=='UP']) #171
down <- intersect(intersect(rownames(raw_degs1)[raw_degs1$change_LFCM2SD=='DOWN'],rownames(DGE_edgeR)[DGE_edgeR$change_LFCM2SD=='DOWN']),
                rownames(DGE_limma)[DGE_limma$change_LFCM2SD=='DOWN']) #453


dat <- log2(cpm(dat)+1) 

library(ggplot2)
pca <- prcomp(t(dat), scale=F)  
pca <- pca[["x"]]
pca <- as.data.frame(pca)

NT <- sapply(strsplit(rownames(pca),"-"),"[",4)

NT <-ifelse(NT %in% "11", "N", "T") 
group <- ArchR::ArchRPalettes[2][[1]][c(14,16)]
names(group) <- c("N","T")

ggplot(pca, aes(PC1,PC2, color = NT)) + geom_point(size=2,alpha=1) + 
  scale_color_manual(values = group)+
  stat_ellipse(data=pca,aes(x=PC1,y=PC2,fill=NT,color=NT),
               geom = "polygon",alpha=0.2,level=0.96,type="t",linetype = 0,show.legend = F)+
  scale_fill_manual(values = group)+theme_bw() +
  theme(
    legend.title = element_blank(),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.text = element_text(color = 'black',size = 15, face = 'plain'),
    axis.title = element_text(color = 'black',size = 15, face = 'plain'),
    axis.ticks = element_line(color = 'black')) + 
  geom_vline(xintercept = 0,linetype="dashed")+geom_hline(yintercept = 0,linetype="dashed")

ggsave("./Transcriptome/TCGA/Fig1TCGA_PCA.pdf")


ggplot(pca, aes(PC1,PC2, color = NT)) + 
  stat_ellipse(data=pca,aes(x=PC1,y=PC2,fill=NT,color=NT),
               geom = "polygon",alpha=0,level=0.96,type="t",linetype = 1,show.legend = F)+
  scale_fill_manual(values = group)+theme_bw() 



 nrDEG_down <- rownames(raw_degs1)[raw_degs1$change_LFCM2SD=='DOWN'] 
 nrDEG_up<- rownames(raw_degs1)[raw_degs1$change_LFCM2SD=='UP']
 choose_gene<- c(nrDEG_down,nrDEG_up)
 choose_matrix<- dat[choose_gene,]
 choose_matrix<- t(scale(t(choose_matrix))) 
 quantile(choose_matrix, 
          c(0.01, 0.99)) 

 choose_matrix[choose_matrix > 2.5]= 2.5
 choose_matrix[choose_matrix < -1.5] = -1.5
 
 group_list <-group_info$type
 names(group_list) <- group_info$sampleID
 group_list <- sort(group_list)
 annotation_col<- as.data.frame(factor(group_list))

 choose_matrix <- choose_matrix[,rownames(annotation_col)]
 filename<- paste("./Transcriptome/TCGA/TCGA_limma_heatmap_allDEG_logFC4.pdf",
                  sep = "", collapse = NULL)
 pheatmap(fontsize = 6, choose_matrix, annotation_col = annotation_col, 
          show_colnames = F, cluster_cols = F, annotation_legend = F,
          filename = filename)
 
 pdf("./Transcriptome/TCGA/TCGA_allDEG_Heatmap_DEGseq.pdf",width = 10,height = 15)
 col <- RColorBrewer::brewer.pal(11,"RdBu")[11:1]
 Heatmap(choose_matrix, 
         col = circlize::colorRamp2(seq(-2,2,length.out = 11),col),
         cluster_rows = T,
         cluster_columns = T,
         show_column_names = F,
         show_row_names = F,
         use_raster=F,
         show_row_dend = T,
         na_col = "grey85",
         clustering_method_rows = "ward.D",
         column_split = annotation_col$`factor(group_list)`,
         column_gap = unit(1, "mm"),
         row_gap = unit(0, "mm"),
         column_title = NULL) 
 dev.off()
 save(list = ls(), file = './Transcriptome/TCGA/TCGA_allDEG.RData' )

 nrDEG_down <- rownames(raw_degs1)[raw_degs1$change_LFCM2SD=='DOWN'] 
 nrDEG_up<- rownames(raw_degs1)[raw_degs1$change_LFCM2SD=='UP']
 choose_gene<- c(nrDEG_down,nrDEG_up)
 choose_matrix<- dat[choose_gene,]
 tcga_signatures=data.frame(
   type=c(rep("normal",length(nrDEG_down)),rep("tumor",length(nrDEG_up))),
   gene=c(nrDEG_down,nrDEG_up)
 )
 write.table(tcga_signatures, "./Transcriptome/TCGA/tcga_signatures.txt",sep = '\t',quote = F,) 
 