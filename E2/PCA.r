setwd("C:/Users/Brise/OneDrive - Universidade de Vigo/Sequentia/E2/Variants_file.vcf")
library("SNPRelate")
vcf.fn<-"C:/Users/Brise/OneDrive - Universidade de Vigo/Sequentia/E2/Variants_file.vcf/Variants_file.vcf"
snpgdsVCF2GDS(vcf.fn, "VCF.gds",  method="biallelic.only")
genofile <- openfn.gds("VCF.gds")
VCF_pca<-snpgdsPCA(genofile, autosome.only=FALSE)

plotdat = data.frame(sample=VCF_pca$sample.id,
                     PC1=VCF_pca$eigenvect[,1],
                     PC2=VCF_pca$eigenvect[,2],
                     PC3=VCF_pca$eigenvect[,3])

ggplot(plotdat, aes(x=PC1, y=PC2)) +
  geom_point(show.legend = FALSE) +
  theme_bw()
