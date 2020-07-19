


#genome-derived features generation
#we build a package for the genome-derived features 
#devtools::install_github("ZhenWei10/m6ALogisticModel")
print("*********")
library(m6ALogisticModel)
library(SummarizedExperiment)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(BSgenome.Hsapiens.UCSC.hg19)
library(fitCons.UCSC.hg19)
library(phastCons100way.UCSC.hg19)

matureSE <- SummarizedExperiment()
Trainingtable <-read.table("/home/sdb2/jing/mature/A549/generealtest/test1.txt",stringsAsFactors = F)
Trainingtable <- Trainingtable[-1,]
Training.predictor.data <-GRanges(seqnames = Trainingtable[,2],
                                  ranges = IRanges(start = as.numeric(Trainingtable[,4]),width = 1)
                                  ,strand = Trainingtable[,6])

print(Training.predictor.data)

rowRanges(matureSE) <- Training.predictor.data

Additional_features_hg19 = list(
  HNRNPC_eCLIP = eCLIP_HNRNPC_gr,
  YTHDC1_TREW = YTHDC1_TREW_gr,
  YTHDF1_TREW = YTHDF1_TREW_gr,
  YTHDF2_TREW = YTHDF2_TREW_gr,
  miR_targeted_genes = miR_targeted_genes_grl,
  TargetScan = TargetScan_hg19_gr,
  Verified_miRtargets = verified_targets_gr,
  METTL3_TREW = METTL3_TREW,
  METTL14_TREW = METTL14_TREW,
  WTAP_TREW = WTAP_TREW,
  METTL16_CLIP = METTL16_CLIP,
  ALKBH5_PARCLIP = ALKBH5_PARCLIP,
  FTO_CLIP = FTO_CLIP,
  FTO_eCLIP = FTO_eCLIP
  
  
)
matureFE <- predictors_annot(se = matureSE,
                             txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,
                             bsgnm = Hsapiens,
                             fc = fitCons.UCSC.hg19,
                             pc = phastCons100way.UCSC.hg19,
                             struct_hybridize = Struc_hg19,
                             feature_lst = Additional_features_hg19,
                             hk_genes_list = HK_hg19_eids,
                             genes_ambiguity_method = "average",
                             standardization=F)

Training.predictor.Gendata <- mcols(matureFE)
write.table(Training.predictor.Gendata,"/home/sdb2/jing/mature/A549/generealtest/Test1Fea.txt")


