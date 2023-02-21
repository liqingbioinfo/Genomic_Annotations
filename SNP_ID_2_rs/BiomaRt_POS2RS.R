########
##Run out of time due to Internet connection
#######

BiocManager::install("biomaRt")

library(biomaRt)
#Check database version 
#ensembl <- useEnsembl(biomart = "snps")
#searchDatasets(mart = ensembl, pattern = "hsapiens")

#            dataset                                                                            description
#           hsapiens_snp         Human Short Variants (SNPs and indels excluding flagged variants) (GRCh38.p13)
#       hsapiens_snp_som Human Somatic Short Variants (SNPs and indels excluding flagged variants) (GRCh38.p13)
#     hsapiens_structvar                                                 Human Structural Variants (GRCh38.p13)
# hsapiens_structvar_som                                         Human Somatic Structural Variants (GRCh38.p13)

## Use the default ENSEMBL Variation Mart & Human dataset
snpMart = useEnsembl(biomart = "snps", 
                     dataset = "hsapiens_snp",
                     mirror = "uswest")

eqtls=read.csv("E:/PHD/Projects/Proteomics/eQTLs_colocalization/Prostate.v8.EUR.signif_pairs.biomaRt.csv")

coords<-eqtls$coords[1:5000]
#> [1] "1:10020:10020" "1:10039:10039"

## Submit the query
getBM(attributes = c('refsnp_id', 'chr_name', 'chrom_start', 'chrom_end', 'allele'),
      filters = c('chromosomal_region'), 
      values = coords, 
      mart = snpMart)  