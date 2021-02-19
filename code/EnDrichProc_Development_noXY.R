source("/group/card2/Evangelyn_Sim/Transcriptome_chromatin_human/Sequencing_ATAC_RNA/20180726_hATACseq_MF/R/1.1kbTSS/7.integration/EnDrich.R")
reactome<-GMT2DF("/group/card2/Evangelyn_Sim/Transcriptome_chromatin_human/Sequencing_ATAC_RNA/20180726_hATACseq_MF/R/1.1kbTSS/7.integration/c2.cp.reactome.v5.0.symbols.gmt")
msigdb<-GMT2DF("/group/card2/Evangelyn_Sim/Transcriptome_chromatin_human/Sequencing_ATAC_RNA/20180726_hATACseq_MF/R/1.1kbTSS/7.integration/msigdb.v5.0.symbols.gmt")
encodedb<-read.table("/group/card2/Evangelyn_Sim/Transcriptome_chromatin_human/Sequencing_ATAC_RNA/20180726_hATACseq_MF/R/1.1kbTSS/7.integration/ENCODE_TFBS.gmt.txt" ,sep="\t")

#compare atac to rna differences
xdat<-read.table("Development_noXY.jn.rnk",header=T,row.names=1)
xdat<-as.matrix(xdat)
zdat<-xdat

reactomeres<-EnDrichMANOVA(zdat, reactome)
reactomeres$p.adjustMANOVA<-p.adjust(reactomeres$pMANOVA,"fdr")
reactomeres<- reactomeres[order(reactomeres$pMANOVA),]
reactomeres$minAbsS<-apply(reactomeres[,4:6],1, function(zz){min(abs(zz))})
write.table(reactomeres, "Development_noXY.EnDrichReactomeResult.txt",sep="\t",col.names=T,row.names=F, quote=F)
ss<-apply(zdat,2,rank)
pdf("Development_noXY.EnDrichReactomeResults.pdf") ; plot2DSets(ss,reactome,reactomeres,resrows=1:500) ; dev.off()

msigres<-EnDrichMANOVA(zdat, msigdb)
msigres$p.adjustMANOVA<-p.adjust(msigres$pMANOVA,"fdr")
msigres<- msigres[order(msigres$pMANOVA),]
msigres$minAbsS<-apply(msigres[,4:6],1, function(zz){min(abs(zz))})
write.table(msigres, "Development_noXY.EndrichMSigDBResult.txt",sep="\t",col.names=T,row.names=F, quote=F)
ss<-apply(zdat,2,rank)
pdf("Development_noXY.EndrichMSigDBResults.pdf") ; plot2DSets(ss,msigdb,msigres,resrows=1:500) ; dev.off()

encoderes<-EnDrichMANOVA(zdat, encodedb)
encoderes$p.adjustMANOVA<-p.adjust(encoderes$pMANOVA,"fdr")
encoderes<- encoderes[order(encoderes$pMANOVA),]
encoderes$minAbsS<-apply(encoderes[,4:6],1, function(zz){min(abs(zz))})
write.table(encoderes,"Development_noXY.EnDrichEncodeResult.txt",sep="\t",col.names=T,row.names=F, quote=F)
ss<-apply(zdat,2,rank)
pdf("Development_noXY.EnDrichEncodeResults.pdf") ; plot2DSets(ss,encodedb,encoderes,resrows=1:500) ; dev.off()

save.image("Development_noXY.EnDrichAnalysis.Rdata")

