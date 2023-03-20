#!$HOME/bin/Rscript

library(ggplot2)
library(scales)
library(reshape2)
library(eulerr)
library(gridExtra)
library(dplyr)

HOME <- getwd()

libs <- c("DZ.20210304_RRS","DZ.20210304_RRS.cgh1")
track_types <- c("22G")
path <- paste0(HOME,"/results/")

transcriptome <- read.table(paste0(HOME,"/reference/WS230.reporter/ce_WS230.coding_transcript.exon.merge.gene.bed"),sep="\t",header=F)

draw_browser_3panel <- function(WBGene,tracks_object,fix_y,boost_22G,transcriptome,lib_names,gene_name){
  
  gene.id <- read.table(paste0(HOME,"/reference/master.v0.22G.mrna.anti.txt"),sep="\t",header=T)
  # translate <- as.character(gene.id$Gene.ID[which(gene.id$Gene.ID==WBGene | gene.id$Gene.name==WBGene | gene.id$Locus.name==WBGene)])
  
  row <- transcriptome[which(transcriptome[,4]%in%WBGene),]
  
  row_edited <- row
  row_edited[1,2] <- 134
  row_edited[1,10] <- 4
  row_edited[1,11] <- c("166,154,164,412")
  row_edited[1,12] <- c("134,351,556,771")

  row <- row_edited

  chr <- as.character(row[1,1])
  range <- row[1,2]:row[1,3]
  range=range+1
  strand <- as.character(row[1,6])
  
  tracks_ls <- list()
  
  for (t in 1:length(tracks_object)){
    df.tmp <- tracks_object[[t]]
    df.sub <- df.tmp[which(tracks_object[[t]][,1]==chr),]
    tracks_ls[[t]] <- bedgraph_to_depth(df.sub,range)
  }
  
  range=range-134
  
  rect_x <- as.numeric(as.vector(unlist(strsplit(x=as.character(row[1,12]),split=","))))+range[1]
  rect_y <- rep(1,times=length(rect_x))
  rect_w <- as.numeric(as.vector(unlist(strsplit(x=as.character(row[1,11]),split=","))))
  
  df.rect <- data.frame(rect_x,
                        rect_y,
                        rect_w)
  
  

  df <- data.frame(range=range)
  df_ls <- list()
  for (i in 1:length(lib_names)){
    df_tmp <- cbind(df,tracks_ls[[i]])
    colnames(df_tmp) <- c("range","value")
    df_tmp$variable <- lib_names[i]
    df_ls[[i]] <- df_tmp
  }
  
  df <- rbind(df_ls[[1]],df_ls[[2]])
  if (length(df_ls) > 2){
    for (x in 3:length(df_ls)){
      df <- rbind(df,df_ls[[x]])
    }
  }
  
  df
  
  df$variable <- factor(df$variable,
                        levels=lib_names,
                        ordered=T)
  
  if (!is.na(fix_y)){
    browser_max <- fix_y
    adjusted_min=-1*browser_max/15
    adjusted_max=-1*browser_max/5
  } else {
    browser_max <- max(df$value)+0.05*max(df$value)
    adjusted_min=-1*max(df$value)/15
    adjusted_max=-1*max(df$value)/5
  }
  
  if (!is.na(boost_22G)){
    df$value[which(df$type=="22G")] <- df$value[which(df$type=="22G")]*boost_22G
  }
  
  browser <- ggplot(data=df,aes(x=range,y=value))+
    facet_wrap(~variable,ncol=1)+
    theme_bw(base_size=15)+theme(aspect.ratio=0.2,
                                 legend.position = 0,
                                 panel.grid = element_blank(),
                                 legend.title = element_blank(),
                                 strip.background = element_blank(),
                                 text = element_text(family = "sans",color="black"))+
    geom_col(fill="blue")+
    # scale_fill_manual(values=c("red","blue"))+
    coord_cartesian(ylim=c((adjusted_min+adjusted_max),browser_max))+
    labs(y="Normalized read depth",x=gene_name)
  
  pirna_coor <- data.frame("21ur-gfp",
                           266+134,
                           286+134)
  
  gene <- ggplot(data=df.rect,aes(xmin=rect_x,xmax=rect_x+rect_w,ymin=rect_y,ymax=rect_y+1))+
    theme_void()+theme(legend.position = 0,
                       legend.title = element_blank())+
    scale_x_discrete(breaks=seq(from=df.rect$rect_x[1],to=df.rect$rect_x[nrow(df.rect)],by=100))+
    coord_cartesian(xlim=c(df.rect$rect_x[1],(df.rect$rect_x[nrow(df.rect)]+df.rect$rect_w[nrow(df.rect)])))+
    geom_rect(fill="black",color="white")+
    geom_rect(data=pirna_coor,aes(xmin=pirna_coor[,2],xmax=pirna_coor[,3],ymin=-0.25,ymax=0.5),fill="darkgreen",color="darkgreen")+
    geom_abline(slope=0,intercept=1.5)
  gene.grob <- ggplotGrob(gene)
  
  browser+
    annotation_custom(
      grob=gene.grob,
      xmin=range[1],
      xmax=range[length(range)],
      ymin=adjusted_min,
      ymax=adjusted_max
    )
}

tracks <- list()
for (y in 1:length(libs)){
  for (x in 1:length(track_types)){
    if (track_types[x]=="22G"){
      tracks[[(y+length(track_types)-2)+x]] <- read.table(paste(HOME,"/results/",libs[y],"/",libs[y],".inserts.v0.22G.norm.anti.bedgraph",sep="")
                                                          ,sep="\t",header=F)
    } else if (track_types[x]=="26G"){
      tracks[[(y+length(track_types)-2)+x]] <- read.table(paste(HOME,"/results/",libs[y],"/",libs[y],".inserts.v0.26G.norm.anti.bedgraph",sep=""),
                                                          sep="\t",header=F)
    }
  }
}


bedgraph_to_depth <- function(df.sub,range){
  select <- df.sub[which(between(df.sub$V2,min(range),max(range))),]
  
  depth <- rep(x=0,times=length(range))
  
  for (i in 1:nrow(select)){
    depth[which(range==select[i,2]):which(range==select[i,3])] <- select[i,4]
  }
  
  depth
}



cgh1_reporter <- draw_browser_3panel(WBGene="cdk-1::gfp",
                    tracks_object=tracks,
                    fix_y=25000,
                    boost_22G=NA,
                    transcriptome=transcriptome,
                    lib_names <- c("Wild type","cgh-1"),
                    gene_name="gfp")

ggsave(cgh1_reporter,path=path,filename="cgh1_reporter.png",
       device="png",dpi=300,height=4,width=5,scale=2)
ggsave(cgh1_reporter,path=path,filename="cgh1_reporter.pdf",
      device="pdf",dpi=300,height=4,width=5,scale=2)




transcriptome_RNAi <- read.table(paste0(HOME,"/reference/WS230.mex5.gfp/ce_WS230.coding_transcript.exon.merge.gene.bed"),sep="\t",header=T)

draw_browser_RNAi <- function(WBGene,tracks_object,fix_y,boost_22G,transcriptome,lib_names,gene_name){
  
  gene.id <- read.table(paste0(HOME,"/reference/master.v0.22G.mrna.anti.txt"),sep="\t",header=T)
  # translate <- as.character(gene.id$Gene.ID[which(gene.id$Gene.ID==WBGene | gene.id$Gene.name==WBGene | gene.id$Locus.name==WBGene)])
  
  row <- transcriptome[which(transcriptome[,4]%in%WBGene),]
  
  chr <- as.character(row[1,1])
  range <- row[1,2]:row[1,3]
  range=range+1
  strand <- as.character(row[1,6])
  
  tracks_ls <- list()
  
  for (t in 1:length(tracks_object)){
    df.tmp <- tracks_object[[t]]
    df.sub <- df.tmp[which(tracks_object[[t]][,1]==chr),]
    tracks_ls[[t]] <- bedgraph_to_depth(df.sub,range)
  }
    
  rect_x <- as.numeric(as.vector(unlist(strsplit(x=as.character(row[1,12]),split=","))))+range[1]
  rect_y <- rep(1,times=length(rect_x))
  rect_w <- as.numeric(as.vector(unlist(strsplit(x=as.character(row[1,11]),split=","))))
  
  df.rect <- data.frame(rect_x,
                        rect_y,
                        rect_w)
  
  
  
  df <- data.frame(range=range)
  df_ls <- list()
  for (i in 1:length(lib_names)){
    df_tmp <- cbind(df,tracks_ls[[i]])
    colnames(df_tmp) <- c("range","value")
    df_tmp$variable <- lib_names[i]
    df_ls[[i]] <- df_tmp
  }
  
  df <- rbind(df_ls[[1]],df_ls[[2]])
  if (length(df_ls) > 2){
    for (x in 3:length(df_ls)){
      df <- rbind(df,df_ls[[x]])
    }
  }
  
  df
  
  df$variable <- factor(df$variable,
                        levels=lib_names,
                        ordered=T)
  
  if (!is.na(fix_y)){
    browser_max <- fix_y
    adjusted_min=-1*browser_max/15
    adjusted_max=-1*browser_max/5
  } else {
    browser_max <- max(df$value)+0.05*max(df$value)
    adjusted_min=-1*max(df$value)/15
    adjusted_max=-1*max(df$value)/5
  }
  
  if (!is.na(boost_22G)){
    df$value[which(df$type=="22G")] <- df$value[which(df$type=="22G")]*boost_22G
  }
  
  browser <- ggplot(data=df,aes(x=range,y=value))+
    facet_wrap(~variable,ncol=1)+
    theme_bw(base_size=15)+theme(aspect.ratio=0.2,
                                 legend.position = 0,
                                 panel.grid = element_blank(),
                                 legend.title = element_blank(),
                                 strip.background = element_blank(),
                                 text = element_text(family = "sans",color="black"))+
    geom_col(fill="blue")+
    coord_cartesian(ylim=c((adjusted_min+adjusted_max),browser_max))+
    labs(y="Normalized read depth",x=gene_name)

  
  gene <- ggplot(data=df.rect,aes(xmin=rect_x,xmax=rect_x+rect_w,ymin=rect_y,ymax=rect_y+1))+
    theme_void()+theme(legend.position = 0,
                       legend.title = element_blank())+
    scale_x_discrete(breaks=seq(from=df.rect$rect_x[1],to=df.rect$rect_x[nrow(df.rect)],by=100))+
    coord_cartesian(xlim=c(df.rect$rect_x[1],(df.rect$rect_x[nrow(df.rect)]+df.rect$rect_w[nrow(df.rect)])))+
    geom_rect(fill="black",color="white")+
    geom_abline(slope=0,intercept=1.5)
  gene.grob <- ggplotGrob(gene)
  
  browser+
    annotation_custom(
      grob=gene.grob,
      xmin=range[1],
      xmax=range[length(range)],
      ymin=adjusted_min,
      ymax=adjusted_max
    )
}


bedgraph_to_depth <- function(df.sub,range){
  select <- df.sub[which(between(df.sub$V2,min(range),max(range))),]
  
  depth <- rep(x=0,times=length(range))
  
  for (i in 1:nrow(select)){
    
    if (max(which(range==select[i,3]),0)==0){
      select[i,3] <- max(range)
    }
    
    depth[which(range==select[i,2]):which(range==select[i,3])] <- select[i,4]
  }
  
  depth
}

libs <- c("DZ.20220708_wt.mex5.gfp.GFP.P0","DZ.20220708_wt.mex5.gfp.GFP.F1","DZ.20220708_wt.mex5.gfp.GFP.F2",
          "DZ.20220708_cgh1.mex5.gfp.GFP.P0","DZ.20220708_cgh1.mex5.gfp.GFP.F1","DZ.20220708_cgh1.mex5.gfp.GFP.F2")
track_types <- c("22G","22G","22G","22G","22G","22G")

tracks <- list()
for (y in 1:length(libs)){
    tracks[[y]] <- read.table(paste(HOME,"/results/bedgraph/",libs[y],".inserts.v0.22G.norm.anti.bedgraph",sep=""),sep="\t",header=F)  
}

cgh1_RNAi <- draw_browser_RNAi(WBGene="gfp",
                                     tracks_object=tracks,
                                     fix_y=100,
                                     boost_22G=NA,
                                     transcriptome=transcriptome_RNAi,
                                     lib_names <- c("Wild type P0","Wild type F1","Wild type F2",
                                                    "cgh-1 P0","cgh-1 F1","cgh-1 F2"),
                                     gene_name="gfp")


ggsave(cgh1_RNAi,path=path,filename="cgh1_RNAi.png",
       device="png",dpi=300,height=4,width=5,scale=2)
ggsave(cgh1_RNAi,path=path,filename="cgh1_RNAi.pdf",
       device="pdf",dpi=300,height=4,width=5,scale=2)



cgh1_RNAi_full <- draw_browser_RNAi(WBGene="gfp",
                               tracks_object=tracks,
                               fix_y=NA,
                               boost_22G=NA,
                               transcriptome=transcriptome_RNAi,
                               lib_names <- c("Wild type P0","Wild type F1","Wild type F2",
                                              "cgh-1 P0","cgh-1 F1","cgh-1 F2"),
                               gene_name="gfp")


ggsave(cgh1_RNAi_full,path=path,filename="cgh1_RNAi_full.png",
       device="png",dpi=300,height=4,width=5,scale=2)
ggsave(cgh1_RNAi_full,path=path,filename="cgh1_RNAi_full.pdf",
       device="pdf",dpi=300,height=4,width=5,scale=2)


plot.meta <- function(data,data_sd,genes,class,type,fig,labels){
  df <- data.frame(NULL)
  for (x in 1:length(data)){
    gene.rows <- data[[x]][data[[x]][,1]%in%genes,]
    sum.row <- colSums(gene.rows[,2:ncol(gene.rows)])
    df <- rbind(df,sum.row)
  }
  
  sum.df <- data.frame(gene.length=1:(ncol(gene.rows)-1),
                       t(df))
  
  colnames(sum.df) <- c("gene.length",names(data))
  sum.df.melt <- melt(sum.df,id.vars="gene.length")
  colnames(sum.df.melt) <- c("gene.length","library","sum")
  
  df_sd <- data.frame(NULL)
  for (x in 1:length(data)){
    gene.rows <- data_sd[[x]][data_sd[[x]][,1]%in%genes,]
    row_sd <- vector(length=0L)
    for (z in 2:ncol(gene.rows)){
      col_sd <- sqrt(sum(gene.rows[,z]^2))
      row_sd[z-1] <- col_sd
    }
    df_sd <- rbind(df_sd,row_sd)
  }
  
  sum.df_sd <- data.frame(gene.length=1:(ncol(gene.rows)-1),
                          t(df_sd))
  
  colnames(sum.df_sd) <- c("gene.length",names(data))
  sum.df_sd.melt <- melt(sum.df_sd,id.vars="gene.length")
  colnames(sum.df_sd.melt) <- c("gene.length","library","sum")
  
  ggplot(data=sum.df.melt,aes(x=gene.length,y=sum,color=library))+
    theme_bw()+theme(aspect.ratio=0.33,legend.title = element_blank(),panel.grid = element_blank(),legend.position = "top")+
    scale_color_manual(values=c("blue","red","green","yellow","purple","orange"),labels=c(labels))+
    geom_line()+
    geom_ribbon(data=sum.df_sd.melt,aes(x=gene.length,ymin=sum.df.melt$sum-sum,ymax=sum.df.melt$sum+sum,fill=library),alpha=0.25,color=NA)+
    scale_fill_manual(values=c("blue","red","green","yellow","purple","orange"),labels=c(labels))+
    labs(x="Gene Length (%)",y="Total RPM",title=paste(class,sep=""))
  
}

avg_rpm <- function(cov1,cov2){
  colnames(cov1) <- c("Gene.ID",
                      paste("r1_",c(1:100),sep=""))
  colnames(cov2) <- c("Gene.ID",
                      paste("r2_",c(1:100),sep=""))
  
  cov_merge <- merge(cov1,cov2,by="Gene.ID")
  
  cov_avg <- data.frame(V1=cov_merge$Gene.ID)
  for (x in 1:100){
    new_col <- (cov_merge[paste("r1_",x,sep="")]+cov_merge[paste("r2_",x,sep="")])/2
    cov_avg[,paste("V",(x+1),sep="")] <- new_col
  }
  
  cov_avg
}

sd_rpm <- function(cov1,cov2){
  colnames(cov1) <- c("Gene.ID",
                      paste("r1_",c(1:100),sep=""))
  colnames(cov2) <- c("Gene.ID",
                      paste("r2_",c(1:100),sep=""))
  
  cov_merge <- merge(cov1,cov2,by="Gene.ID")
  
  cov_sd <- data.frame(V1=cov_merge$Gene.ID)
  for (x in 1:100){
    new_col <- apply(data.frame(r1=cov_merge[paste("r1_",x,sep="")],
                                r2=cov_merge[paste("r2_",x,sep="")]),
                     1,
                     sd)
    cov_sd[,paste("V",(x+1),sep="")] <- new_col
  }
  
  cov_sd
}


avg_rpm_3 <- function(cov1,cov2,cov3){
  colnames(cov1) <- c("Gene.ID",
                      paste("r1_",c(1:100),sep=""))
  colnames(cov2) <- c("Gene.ID",
                      paste("r2_",c(1:100),sep=""))
  colnames(cov3) <- c("Gene.ID",
                      paste("r3_",c(1:100),sep=""))
  
  cov_merge <- merge(cov1,cov2,by="Gene.ID")
  cov_merge <- merge(cov_merge,cov3,by="Gene.ID")
  
  cov_avg <- data.frame(V1=cov_merge$Gene.ID)
  for (x in 1:100){
    new_col <- (cov_merge[paste("r1_",x,sep="")]+cov_merge[paste("r2_",x,sep="")]+cov_merge[paste("r3_",x,sep="")])/3
    cov_avg[,paste("V",(x+1),sep="")] <- new_col
  }
  
  cov_avg
}

sd_rpm_3 <- function(cov1,cov2,cov3){
  colnames(cov1) <- c("Gene.ID",
                      paste("r1_",c(1:100),sep=""))
  colnames(cov2) <- c("Gene.ID",
                      paste("r2_",c(1:100),sep=""))
  colnames(cov3) <- c("Gene.ID",
                      paste("r3_",c(1:100),sep=""))
  
  cov_merge <- merge(cov1,cov2,by="Gene.ID")
  cov_merge <- merge(cov_merge,cov3,by="Gene.ID")
  
  cov_sd <- data.frame(V1=cov_merge$Gene.ID)
  for (x in 1:100){
    new_col <- apply(data.frame(r1=cov_merge[paste("r1_",x,sep="")],
                                r2=cov_merge[paste("r2_",x,sep="")],
                                r3=cov_merge[paste("r3_",x,sep="")]),
                     1,
                     sd)
    cov_sd[,paste("V",(x+1),sep="")] <- new_col
  }
  
  cov_sd
}

data_22G <- read.table(paste0(HOME,"/results/master/master.v0.22G.mrna.anti.txt"),header=T,sep="\t")
all.WBGene <- data_22G$Gene.ID

wt_wago4_WBGene <- vector(length=0L)

rpm_cutoff=5
fold_cutoff=2
arb=sort(unique(data_22G$DZ.20211015_wt.csr1.ip.r1),decreasing=F)[2]

for (x in 1:nrow(data_22G)){
  
  if (data_22G$DZ.20211015_wt.wago4.ip.r1[x]>=rpm_cutoff &
      ((data_22G$DZ.20211015_wt.wago4.ip.r1[x]+arb)/(data_22G$DZ.20211015_wt.wago4.input.r1[x]+arb))>=fold_cutoff){
    wt_wago4_WBGene <- c(wt_wago4_WBGene,data_22G$Gene.ID[x])
  }
  
}




meta_rpm_wt_wago4_r1 <- read.table(paste0(HOME,"/cov/DZ.20211015_wt.wago4.ip.r1/DZ.20211015_wt.wago4.ip.r1_22G_rpm.txt"),header=F,sep="\t")
meta_rpm_cgh1_wago4_r1 <- read.table(paste0(HOME,"/cov/DZ.20211015_cgh1.wago4.ip.r1/DZ.20211015_cgh1.wago4.ip.r1_22G_rpm.txt"),header=F,sep="\t")

gene.info <- read.table(paste0(HOME,"/reference/master.v0.22G.mrna.anti.txt"),sep="\t",header=T)
csr.WBGene <- gene.info$Gene.ID[which(gene.info$CSR.target==T)]

meta_rpm_wt_wago4_sd <- meta_rpm_wt_wago4_r1
meta_rpm_wt_wago4_sd[,2:101] <- 0
meta_rpm_cgh1_wago4_sd <- meta_rpm_cgh1_wago4_r1
meta_rpm_cgh1_wago4_sd[,2:101] <- 0


data.wago <- list( meta_rpm_wt_wago4_r1,
             meta_rpm_cgh1_wago4_r1)
names(data.wago) <- c("Wild type WAGO-4 IP","cgh-1 WAGO-4 IP")

data.wago_sd <- list(meta_rpm_wt_wago4_sd,
                meta_rpm_cgh1_wago4_sd)
names(data.wago_sd) <- c("Wild type WAGO-4 IP","cgh-1 WAGO-4 IP")

gg_meta_wago4.IP_wago4 <- plot.meta(data=data.wago,
                                   data_sd=data.wago_sd,
                                   genes=wt_wago4_WBGene,
                                   class=paste("WAGO-4 targets n=",length(wt_wago4_WBGene),sep=""),
                                   type="22G",
                                   labels=c("Wild type WAGO-4 IP","cgh-1 WAGO-4 IP"))

ggsave(gg_meta_wago4.IP_wago4,path=path,filename="gg_meta_wago4.IP_wago4.png",
       device="png",dpi=300,width=10,height=3)

ggsave(gg_meta_wago4.IP_wago4,path=path,filename="gg_meta_wago4.IP_wago4.pdf",
       device="pdf",dpi=300,width=10,height=3)

