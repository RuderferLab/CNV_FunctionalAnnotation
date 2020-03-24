args = commandArgs(trailingOnly=TRUE)

print(paste("CNV File: ",args[1],sep=" "))
cnv_file <- args[1]
data = read.table(cnv_file,header=T)
data$SVLen = data$end - data$start 
all = data


# coding proportion
data.exon = read.table(paste(cnv_file,"exon",sep="-"))
data.gene = read.table("Ensemblv75/Ensemblv75_merged_exons",header=T)
data.exon.sum = aggregate(data.exon$V3,list(data.exon$V1,data.exon$V2),sum)
data.gene = data.gene[,c(1,6)]
data.exon.sum = merge(data.exon.sum,data.gene,by.x="Group.2",by.y="Geneid")
data.exon.sum$exon = data.exon.sum$x/data.exon.sum$Length
data.exon.sum = data.exon.sum[,c(2,1,5)]

# enhancer proportion
data.enh = read.table(paste(cnv_file,"PEC-enhancer-predictions",sep="-"))
data.enh.sum = aggregate(data.enh$V3,list(data.enh$V1,data.enh$V2),sum)

# in tad
data.tad.gene = read.table(paste(cnv_file,"gene-tad",sep="-"))
data.tad.sv = read.table(paste(cnv_file,"tad",sep="-"))
data.intad = merge(data.tad.sv,data.tad.gene,by="V2")
data.intad = unique(data.intad[,c(2,3)])
data.intad$intad = 1

# transcript proportion
data.trans = read.table(paste(cnv_file,"gene",sep="-"))

# promoter
data.pro = read.table(paste(cnv_file,"promoters",sep="-"))

data.annot = merge(data.trans,data.exon.sum,by.x=c("V1","V2"),by.y=c("Group.1","Group.2"),all.x=T,all.y=T)
data.annot = merge(data.annot,data.pro,by=c("V1","V2"),all.x=T,all.y=T)
data.annot = merge(data.annot,data.enh.sum,by.x=c("V1","V2"),by.y=c("Group.1","Group.2"),all.x=T,all.y=T)
data.annot = merge(data.annot,data.intad,by.x=c("V1","V2"),by.y=c("V1.x","V1.y"),all.x=T)
names(data.annot) = c("id","gene","trans","ExonProp","pro.sum","enh.sum","inTad")
data.annot[is.na(data.annot$ExonProp),]$ExonProp = 0
data.annot[is.na(data.annot$pro.sum),]$pro.sum = 0
data.annot[is.na(data.annot$enh.sum),]$enh.sum = 0
data.annot[is.na(data.annot$trans),]$trans = 0
data.annot[is.na(data.annot$inTad),]$inTad = 0

data = merge(data,data.annot,by="id")

print("Reading in models...")
load("CNV-models.Rdata")


print("Applying models...")
data$pred_exp = NA
data[data$SVType == "DEL",]$pred_exp = predict(pfc.model.del,data[data$SVType == "DEL",])
data[data$SVType == "DUP",]$pred_exp = predict(pfc.model.dup,data[data$SVType == "DUP",])

data$sum = data$ExonProp + data$enh.sum + data$pro.sum
data$cre = data$enh.sum + data$pro.sum
data$gene.cnt = 0
data[data$ExonProp > 0,]$gene.cnt = 1

ginfo = read.table("Ensemblv75/Ensemblv75-hgnc.convert",header=T,sep="\t")
names(ginfo) = c("gene","symbol")

print("Adding constraint information...")
# gnomAD constraint (https://storage.cloud.google.com/gnomad-public/release/2.1/ht/constraint/constraint.txt.bgz)
p = read.table("release_2.1_ht_constraint_constraint.txt",header=T)
p = p[p$canonical == "true",]
p = p[,c(1,6,22)]
names(p)[1] = "symbol"

p$oe_z = (p$oe_lof - mean(p$oe_lof,na.rm=T))/sd(p$oe_lof,na.rm=T)
p$oe_z = p$oe_z * -1

# windsorize at 5%
p[p$oe_z < -1.548 & !is.na(p$oe_z),]$oe_z = -1.548
p$oe.norm = (p$oe_z + abs(min(p$oe_z,na.rm=T)))/(max(p$oe_z,na.rm=T)+abs(min(p$oe_z,na.rm=T)))

p = merge(p,ginfo,by="symbol",all.x=T)

p = p[,-1]
p = unique(p)

pf = as.data.frame(unique(p[!is.na(p$gene),]$gene))
names(pf) = "gene"
oe = aggregate(p$oe.norm,list(p$gene),mean)
names(oe) = c("gene","oe")

pf = merge(pf,oe,all.x=T,by="gene")

data = merge(data,pf,by="gene",all.x=T)
data$wt = data$pred_exp * data$oe


print("Creating CNV level regulatory disruption scores...")
cnv.ngene = aggregate(data$gene.cnt,list(data$id),sum,na.rm=T)
cnv.exon = aggregate(data$ExonProp,list(data$id),sum,na.rm=T)
cnv.enh = aggregate(data$enh.sum,list(data$id),sum,na.rm=T)
cnv.pro = aggregate(data$pro.sum,list(data$id),sum,na.rm=T)
cnv.exp = aggregate(data$pred_exp,list(data$id),sum,na.rm=T)
cnv.wt = aggregate(data[!is.na(data$wt),]$wt,list(data[!is.na(data$wt),]$id),sum,na.rm=T)

names(cnv.ngene) = c("id","ngenes")
names(cnv.exon) = c("id","exon")
names(cnv.enh) = c("id","enh")
names(cnv.pro) = c("id","pro")
names(cnv.exp) = c("id","exp")
names(cnv.wt) = c("id","wt")

cnv = merge(all,cnv.ngene,by="id",all.x=T)
cnv = merge(cnv,cnv.exon,by="id",all.x=T)
cnv = merge(cnv,cnv.enh,by="id",all.x=T)
cnv = merge(cnv,cnv.pro,by="id",all.x=T)
cnv = merge(cnv,cnv.exp,by="id",all.x=T)
cnv = merge(cnv,cnv.wt,by="id",all.x=T)

cnv[is.na(cnv)] = 0

write.table(cnv,paste(cnv_file,"reg_disruption",sep="."),quote=F,row.names=F,col.names=T)