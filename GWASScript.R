
setwd("/Users/Emix/Desktop/snpanalysis")
library(detectRUNS)
library(tidyverse)
library(janitor)
library(qqman)


#Step 1: Quality Control of SNPs and Samples
#missigness per dog < 0.05 
#minor allele frequency > 0.01 
#missingness per marker < 0.05 
#removal of genetically identical donors (kinship values exceeding 0.354)
#no MT, no sex, no unlocalized SNPs
#for each breed: created ped and bed files after QC 

#for GS
system("./plink2 --bfile GSD_fs --dog --autosome --mind 0.05 --geno 0.05 --maf 0.01 --king-cutoff 0.354 --export ped --out GS_afterQC ")
system("./plink --file GS_afterQC --dog --autosome --make-bed --out GS_afterQC ")


#for LR 
system("./plink2 --bfile LR_fs --dog --autosome --mind 0.05 --geno 0.05 --maf 0.01 --king-cutoff 0.354 --export ped --out LR_afterQC ")
system("./plink --file LR_afterQC --dog --autosome --make-bed --out LR_afterQC ")

#need to change tab-delimited ped to space-delimited with plink 1.9 ped 
system("./plink --file GS_afterQC --dog --recode --out GSafterQC_notab")
system("./plink --file LR_afterQC --dog --recode --out LRafterQC_notab")


#Step 2: Calculate genomic inbreeding coefficients by calculating runs of homozygosity with both 
#sliding window and consecutive SNP methods. 

#for GS
gsd_slideruns <- slidingRUNS.run(genotypeFile = "GSafterQC_notab.ped" , 
                               mapFile = "GSafterQC_notab.map", 
                               minSNP = 2, 
                               minLengthBps = 100, 
                               maxGap = 500000, 
                               maxOppRun = 0, 
                               maxMissRun = 0)

gsd_consecutiveRuns <- consecutiveRUNS.run( genotypeFile ="GSafterQC_notab.ped",
                                            mapFile = "GSafterQC_notab.map",
                                            minSNP = 2,
                                            maxGap = 500000,
                                            minLengthBps = 100,
                                            maxOppRun = 0,
                                            maxMissRun = 0)

gs_slide_summaryList <- summaryRuns(
  runs = gsd_slideruns, mapFile = "GSafterQC_notab.map", genotypeFile = "GSafterQC_notab.ped")

gs_consecutive_summaryList <- summaryRuns(
  runs = gsd_consecutiveRuns, mapFile = "GSafterQC_notab.map", genotypeFile = "GSafterQC_notab.ped")

#for LR 
lr_slideruns <- slidingRUNS.run(genotypeFile = "LRafterQC_notab.ped" , 
                                 mapFile = "LRafterQC_notab.map", 
                                 minSNP = 2, 
                                 minLengthBps = 100, 
                                 maxGap = 500000, 
                                 maxOppRun = 0, 
                                 maxMissRun = 0)

lr_consecutiveRuns <- consecutiveRUNS.run( genotypeFile ="LRafterQC_notab.ped",
                                            mapFile = "LRafterQC_notab.map",
                                            minSNP = 2,
                                            maxGap = 500000,
                                            minLengthBps = 100,
                                            maxOppRun = 0,
                                            maxMissRun = 0)

lr_slide_summaryList <- summaryRuns(
  runs = lr_slideruns, mapFile = "LRafterQC_notab.map", genotypeFile = "LRafterQC_notab.ped")

lr_consecutive_summaryList <- summaryRuns(
  runs = lr_consecutiveRuns, mapFile = "LRafterQC_notab.map", genotypeFile = "LRafterQC_notab.ped")

#calculate genomic inbreeding coefficient (FROH) genome wide for each 

 total_Froh_gs <- mean(mean(gs_slide_summaryList$result_Froh_genome_wide[,4]) + 
   mean(gs_consecutive_summaryList$result_Froh_genome_wide[,4]))
 
 total_Froh_lr <- mean(mean(lr_slide_summaryList$result_Froh_genome_wide[,4]) + 
                         mean(lr_consecutive_summaryList$result_Froh_genome_wide[,4]))
 Froh_results <- capture.output('Froh for GS:', print(total_Froh_gs), 'Froh for LR:', print(total_Froh_lr))
 writeLines(Froh_results, con = file("froh_results.txt"))
 #Step 3: Compare obtained genomic inbreeding coefficients using Welch's two sample t-test. 
 #use consecutive method bc significantly more runs detected  
 breeds <- c('GS', 'LR')
 png("FROHcomparison.png")
 boxplot(gs_consecutive_summaryList$result_Froh_genome_wide[,4], 
         lr_consecutive_summaryList$result_Froh_genome_wide[,4], 
         names = breeds , ylab = 'FROH')
 dev.off()
 
 t_testresults <- capture.output(t.test(gs_consecutive_summaryList$result_Froh_genome_wide[,4], lr_consecutive_summaryList$result_Froh_genome_wide[,4]))
 writeLines(t_testresults, con = file("t_testresults.txt"))
 #therefore no significant difference between breeds 
 
 
 #Step 4: Perform PCA and visualize results. 
 #for GS
 system("./plink --file GSafterQC_notab --dog --pca --out GS_PCA")
 gs_eigenvalues <- read_delim("GS_PCA.eigenval", delim = " ", col_names = F)
 gs_eigenvectors <- read_delim("GS_PCA.eigenvec", delim = " ", col_names = F)
 gs_eigen_percent <- round(((gs_eigenvalues)/sum(gs_eigenvalues))*100,2)  

 #for LR 
 system("./plink --file LRafterQC_notab --dog --pca --out LR_PCA")
 lr_eigenvalues <- read_delim("LR_PCA.eigenval", delim = " ", col_names = F)
 lr_eigenvectors <- read_delim("LR_PCA.eigenvec", delim = " ", col_names = F)
 lr_eigen_percent <- round(((lr_eigenvalues)/sum(lr_eigenvalues))*100,2)  
 
 #bring in phenotype information for PCA visualization 
 #GS
 gsd_phenotype <- read.table("GSD_fs_pheno.txt", sep = "", header = TRUE, stringsAsFactors = FALSE)
 gsd_phenotype_t <- data.frame(t(gsd_phenotype)) %>%  
   row_to_names(row_number = 1) %>% 
   mutate(FID = 'GermanShepherd') %>% 
   rename(IID = DogID) %>% 
   relocate(FID)
 #save phenotype information 
 write_delim(gsd_phenotype_t, "gsd_phenotype.txt", delim = " ")
 #merge with eigeninfo
 gs_eigen_phenotype <- merge(gs_eigenvectors, gsd_phenotype_t, 
                                 by.x = 'X2', by.y = 'IID')
 
 
#LR
 lr_phenotype <- read.table("LR_fs_pheno.txt", sep = "", header = TRUE, stringsAsFactors = FALSE)
 lr_phenotype_t <- data.frame(t(lr_phenotype)) %>%  
   row_to_names(row_number = 1) %>% 
   mutate(FID = 'LabradorRetriever') %>% 
   rename(IID = DogID) %>% 
   relocate(FID)
 #save phenotype information 
 write_delim(lr_phenotype_t, "lr_phenotype.txt", delim = " ")
 #merge with eigeninfo
 lr_eigen_phenotype <- merge(lr_eigenvectors, lr_phenotype_t, 
                             by.x = 'X2', by.y = 'IID')
 
 
 #plot PCAs with color denoting qualified or not 
 png("gs_pca.png")
 ggplot(data = gs_eigen_phenotype) + 
   geom_point(mapping = aes(x = X3, y=X4, color = qualification), size = 3, show.legend = TRUE)+
   geom_hline(yintercept = 0, linetype = 'dotted') + 
   geom_vline(xintercept = 0, linetype = 'dotted') +
   labs(title = 'PCA of german shepherd population', 
        x = paste0("Principal component 1(", gs_eigen_percent[1,1], " %)"), 
        y = paste0("Principal component 2(", gs_eigen_percent[2,1], " %)")) + 
   theme_minimal()
dev.off() 

png("lr_pca.png")
ggplot(data = lr_eigen_phenotype) + 
   geom_point(mapping = aes(x = X3, y=X4, color = qualification), size = 3, show.legend = TRUE)+
   geom_hline(yintercept = 0, linetype = 'dotted') + 
   geom_vline(xintercept = 0, linetype = 'dotted') +
   labs(title = 'PCA of labrador retriever population', 
        x = paste0("Principal component 1(", lr_eigen_percent[1,1], " %)"), 
        y = paste0("Principal component 2(", lr_eigen_percent[2,1], " %)")) + 
   theme_minimal()
 dev.off()
 
 #Population structure clearly detected in GS, slight population structure in LR
 
 #Step 5: Run a GWAS with gemma for each phenotype: 7 behavioral traits and the qualification outcome 
 #and save the significant SNPS detected. 
 
 #for Qualification (binary): 
 #GS
 #prepare bed input file 
 system(str_c("./plink2 --bfile GS_afterQC --dog --pheno gsd_phenotype.txt --pheno-name ", 'Qualification', 
              ' --make-bed --out gsd_gemma_qualif'))
 #estimate a relatedness matrix with bed files 
 system("gemma -bfile gsd_gemma_qualif -gk 1 -o gsd_RelMat_qualif")  
 system("gemma -bfile gsd_gemma_qualif -k ./output/gsd_RelMat_qualif.cXX.txt -lmm 2 -o gsd_qualif_GWASresults.lmm")
 
 #read results in 
 gsd_qualif_results <- read_table("./output/gsd_qualif_GWASresults.lmm.assoc.txt")
 #compute the Bonferroni threshold -> will be same for all association tests for GSD
 bonferroni_gsd<- -log10(0.05/ nrow(gsd_qualif_results)) 
 #manhattan plot 
 png("GWAS_gsd_qualif.png")
 manhattan(gsd_qualif_results,chr="chr",bp="ps",p="p_lrt",snp="rs",genomewideline=bonferroni_gsd)
 dev.off()
 #qqplot
 png("GWAS_gsd_qualif_qq.png")
 qq(gsd_qualif_results$p_lrt)
 dev.off()
 
 #significant SNPS above threshold on chromosome 10 -> saving these SNPs to table
 gsd_sig_qualif <- gsd_qualif_results %>% 
   mutate(negLogP = -log10(p_lrt)) %>% 
   select(chr, rs, p_lrt, negLogP) %>% 
   filter(negLogP > 5)
 
 write_csv(gsd_sig_qualif, "gsd_sig_qualif.csv")
 
 #LR 
 system(str_c("./plink2 --bfile LR_afterQC --dog --pheno lr_phenotype.txt --pheno-name ", 'Qualification', 
              ' --make-bed --out lr_gemma_qualif'))
 #estimate a relatedness matrix with bed files 
 system("gemma -bfile lr_gemma_qualif -gk 1 -o lr_RelMat_qualif")  
 system("gemma -bfile lr_gemma_qualif -k ./output/lr_RelMat_qualif.cXX.txt -lmm 2 -o lr_qualif_GWASresults.lmm")
 #read results in 
 lr_qualif_results <- read_table("./output/lr_qualif_GWASresults.lmm.assoc.txt")
 #compute the Bonferroni threshold -> will be same for all association tests for LR
 bonferroni_lr<- -log10(0.05/ nrow(lr_qualif_results))
 #manhattan plot
 png("GWAS_lr_qualif.png")
 manhattan(lr_qualif_results,chr="chr",bp="ps",p="p_lrt",snp="rs",genomewideline=bonferroni_lr)
 dev.off()
 #qqplot
 png("GWAS_lr_qualif_qq.png")
 qq(lr_qualif_results$p_lrt)
 dev.off()
 
 #significant SNPS above threshold on chromosome 23 -> saving these SNPs to table
 lr_sig_qualif <- lr_qualif_results %>% 
   mutate(negLogP = -log10(p_lrt)) %>% 
   select(chr, rs, p_lrt, negLogP) %>% 
   filter(negLogP > 5)
 
 write_csv(lr_sig_qualif, "lr_sig_qualif.csv")
 
 
 
 #for Activity: 
 #GS
 system(str_c("./plink2 --bfile GS_afterQC --dog --pheno gsd_phenotype.txt --pheno-name ", 'Activity', 
              ' --make-bed --out gsd_gemma_active'))
 #estimate a relatedness matrix with bed files 
 system("gemma -bfile gsd_gemma_active -gk 1 -o gsd_RelMat_active")  
 system("gemma -bfile gsd_gemma_active -k ./output/gsd_RelMat_active.cXX.txt -lmm 2 -o gsd_active_GWASresults.lmm")
 
 #read results in 
 gsd_active_results <- read_table("./output/gsd_active_GWASresults.lmm.assoc.txt")
 #use the Bonferroni threshold for gsd 
 #manhattan plot 
 png("GWAS_gsd_active.png")
 manhattan(gsd_active_results,chr="chr",bp="ps",p="p_lrt",snp="rs",genomewideline=bonferroni_gsd)
 dev.off()
 #qqplot
 png("GWAS_gsd_active_qq.png")
 qq(gsd_active_results$p_lrt)
 dev.off()
 
 #no significant SNPs detected -> none above threshold 
 
 #LR 
 system(str_c("./plink2 --bfile LR_afterQC --dog --pheno lr_phenotype.txt --pheno-name ", 'Activity', 
              ' --make-bed --out lr_gemma_active'))
 #estimate a relatedness matrix with bed files 
 system("gemma -bfile lr_gemma_active -gk 1 -o lr_RelMat_active")  
 system("gemma -bfile lr_gemma_active -k ./output/lr_RelMat_active.cXX.txt -lmm 2 -o lr_active_GWASresults.lmm")
 #read results in 
 lr_active_results <- read_table("./output/lr_active_GWASresults.lmm.assoc.txt")
 #use the Bonferroni threshold for lr
 #manhattan plot 
 png("GWAS_lr_active.png")
 manhattan(lr_active_results,chr="chr",bp="ps",p="p_lrt",snp="rs",genomewideline=bonferroni_lr)
 dev.off()
 #qqplot
 png("GWAS_lr_active_qq.png")
 qq(lr_active_results$p_lrt)
 dev.off()
 
 #significant SNPS above threshold on chromosome 3 -> saving these SNPs to table
 lr_sig_active <- lr_active_results %>% 
   mutate(negLogP = -log10(p_lrt)) %>% 
   select(chr, rs, p_lrt, negLogP) %>% 
   filter(negLogP > 5)
 
 write_csv(lr_sig_active, "lr_sig_active.csv")
 
 
 
 #for Independency: 
 #GS
 system(str_c("./plink2 --bfile GS_afterQC --dog --pheno gsd_phenotype.txt --pheno-name ", 'Independency', 
              ' --make-bed --out gsd_gemma_independent'))
 #estimate a relatedness matrix with bed files 
 system("gemma -bfile gsd_gemma_independent -gk 1 -o gsd_RelMat_independent")  
 system("gemma -bfile gsd_gemma_independent -k ./output/gsd_RelMat_independent.cXX.txt -lmm 2 -o gsd_independent_GWASresults.lmm")
 
 #read results in 
 gsd_independent_results <- read_table("./output/gsd_independent_GWASresults.lmm.assoc.txt")
 #use the Bonferroni threshold for gsd 
 #manhattan plot 
 png("GWAS_gsd_indepen.png")
 manhattan(gsd_independent_results,chr="chr",bp="ps",p="p_lrt",snp="rs",genomewideline=bonferroni_gsd)
 dev.off()
 #qqplot
 png("GWAS_gsd_indepen_qq.png")
 qq(gsd_independent_results$p_lrt)
 dev.off()
 
 #significant SNPS above threshold on chromosome 15 -> saving these SNPs to table
 gsd_sig_independent <- gsd_independent_results %>% 
   mutate(negLogP = -log10(p_lrt)) %>% 
   select(chr, rs, p_lrt, negLogP) %>% 
   filter(negLogP > 5)
 
 write_csv(gsd_sig_independent, "gsd_sig_independent.csv")
 
 #LR 
 system(str_c("./plink2 --bfile LR_afterQC --dog --pheno lr_phenotype.txt --pheno-name ", 'Independency', 
              ' --make-bed --out lr_gemma_independent'))
 #estimate a relatedness matrix with bed files 
 system("gemma -bfile lr_gemma_independent -gk 1 -o lr_RelMat_independent")  
 system("gemma -bfile lr_gemma_independent -k ./output/lr_RelMat_independent.cXX.txt -lmm 2 -o lr_independent_GWASresults.lmm")
 #read results in 
 lr_independent_results <- read_table("./output/lr_independent_GWASresults.lmm.assoc.txt")
 #use the Bonferroni threshold for lr 
 #manhattan plot 
 png("GWAS_lr_indepen.png")
 manhattan(lr_independent_results,chr="chr",bp="ps",p="p_lrt",snp="rs",genomewideline=bonferroni_lr)
 dev.off()
 #qqplot
 png("GWAS_lr_indepen_qq.png")
 qq(lr_independent_results$p_lrt)
 dev.off()
 
 #significant SNPS above threshold on chromosome 8 and 15 -> saving these SNPs to table
 lr_sig_independent <- lr_independent_results %>% 
   mutate(negLogP = -log10(p_lrt)) %>% 
   select(chr, rs, p_lrt, negLogP) %>% 
   filter(negLogP > 5)
 
 write_csv(lr_sig_independent, "lr_sig_independent.csv")
 
 
 
 #for Concentration: 
 #GS
 system(str_c("./plink2 --bfile GS_afterQC --dog --pheno gsd_phenotype.txt --pheno-name ", 'Concentration', 
              ' --make-bed --out gsd_gemma_concen'))
 #estimate a relatedness matrix with bed files 
 system("gemma -bfile gsd_gemma_concen -gk 1 -o gsd_RelMat_concen")  
 system("gemma -bfile gsd_gemma_concen -k ./output/gsd_RelMat_concen.cXX.txt -lmm 2 -o gsd_concen_GWASresults.lmm")
 
 #read results in 
 gsd_concen_results <- read_table("./output/gsd_concen_GWASresults.lmm.assoc.txt")
 #use the Bonferroni threshold for gsd 
 #manhattan plot 
 png("GWAS_gsd_concen.png")
 manhattan(gsd_concen_results,chr="chr",bp="ps",p="p_lrt",snp="rs",genomewideline=bonferroni_gsd)
 dev.off()
 #qqplot
 png("GWAS_gsd_concen_qq.png")
 qq(gsd_concen_results$p_lrt)
 dev.off()
 
 #no significant SNPs detected -> none above threshold 
 
 #LR 
 system(str_c("./plink2 --bfile LR_afterQC --dog --pheno lr_phenotype.txt --pheno-name ", 'Concentration', 
              ' --make-bed --out lr_gemma_concen'))
 #estimate a relatedness matrix with bed files 
 system("gemma -bfile lr_gemma_concen -gk 1 -o lr_RelMat_concen")  
 system("gemma -bfile lr_gemma_concen -k ./output/lr_RelMat_concen.cXX.txt -lmm 2 -o lr_concen_GWASresults.lmm")
 #read results in 
 lr_concen_results <- read_table("./output/lr_concen_GWASresults.lmm.assoc.txt")
 #use the Bonferroni threshold for lr 
 #manhattan plot 
 png("GWAS_lr_concen.png")
 manhattan(lr_concen_results,chr="chr",bp="ps",p="p_lrt",snp="rs",genomewideline=bonferroni_lr)
 dev.off()
 #qqplot
 png("GWAS_lr_concen_qq.png")
 qq(lr_concen_results$p_lrt)
 dev.off()
 
 #significant SNPS above threshold on chromosome 7 -> saving these SNPs to table
 lr_sig_concen <- lr_concen_results %>% 
   mutate(negLogP = -log10(p_lrt)) %>% 
   select(chr, rs, p_lrt, negLogP) %>% 
   filter(negLogP > 5)
 
 write_csv(lr_sig_concen, "lr_sig_concen.csv")
 
 
 
 #for Friendliness to humans: 
 #GS
 system(str_c("./plink2 --bfile GS_afterQC --dog --pheno gsd_phenotype.txt --pheno-name ", 'Friendliness_to_humans', 
              ' --make-bed --out gsd_gemma_friend'))
 #estimate a relatedness matrix with bed files 
 system("gemma -bfile gsd_gemma_friend -gk 1 -o gsd_RelMat_friend")  
 system("gemma -bfile gsd_gemma_friend -k ./output/gsd_RelMat_friend.cXX.txt -lmm 2 -o gsd_friend_GWASresults.lmm")
 
 #read results in 
 gsd_friend_results <- read_table("./output/gsd_friend_GWASresults.lmm.assoc.txt")
 #use the Bonferroni threshold for gsd 
 #manhattan plot 
 png("GWAS_gsd_friend.png")
 manhattan(gsd_friend_results,chr="chr",bp="ps",p="p_lrt",snp="rs",genomewideline=bonferroni_gsd)
 dev.off()
 #qqplot
 png("GWAS_gsd_friend_qq.png")
 qq(gsd_friend_results$p_lrt)
 dev.off()
 
 #no significant SNPs detected -> none above threshold 
 
 #LR 
 system(str_c("./plink2 --bfile LR_afterQC --dog --pheno lr_phenotype.txt --pheno-name ", 'Friendliness_to_humans', 
              ' --make-bed --out lr_gemma_friend'))
 #estimate a relatedness matrix with bed files 
 system("gemma -bfile lr_gemma_friend -gk 1 -o lr_RelMat_friend")  
 system("gemma -bfile lr_gemma_friend -k ./output/lr_RelMat_friend.cXX.txt -lmm 2 -o lr_friend_GWASresults.lmm")
 #read results in 
 lr_friend_results <- read_table("./output/lr_friend_GWASresults.lmm.assoc.txt")
 #use the Bonferroni threshold for lr 
 #manhattan plot 
 png("GWAS_lr_friend.png")
 manhattan(lr_friend_results,chr="chr",bp="ps",p="p_lrt",snp="rs",genomewideline=bonferroni_lr)
 dev.off()
 #qqplot
 png("GWAS_lr_friend_qq.png")
 qq(lr_friend_results$p_lrt)
 dev.off()
 
 #significant SNPS above threshold on chromosome 15 -> saving these SNPs to table
 lr_sig_friend <- lr_friend_results %>% 
   mutate(negLogP = -log10(p_lrt)) %>% 
   select(chr, rs, p_lrt, negLogP) %>% 
   filter(negLogP > 5)
 
 write_csv(lr_sig_friend, "lr_sig_friend.csv")
 
 
 
 #for Tolerance to dogs: 
 #GS
 system(str_c("./plink2 --bfile GS_afterQC --dog --pheno gsd_phenotype.txt --pheno-name ", 'Tolerance_to_dogs', 
              ' --make-bed --out gsd_gemma_toler'))
 #estimate a relatedness matrix with bed files 
 system("gemma -bfile gsd_gemma_toler -gk 1 -o gsd_RelMat_toler")  
 system("gemma -bfile gsd_gemma_toler -k ./output/gsd_RelMat_toler.cXX.txt -lmm 2 -o gsd_toler_GWASresults.lmm")
 
 #read results in 
 gsd_toler_results <- read_table("./output/gsd_toler_GWASresults.lmm.assoc.txt")
 #use the Bonferroni threshold for gsd 
 #manhattan plot 
 png("GWAS_gsd_toler.png")
 manhattan(gsd_toler_results,chr="chr",bp="ps",p="p_lrt",snp="rs",genomewideline=bonferroni_gsd)
 dev.off()
 #qqplot
 png("GWAS_gsd_toler_qq.png")
 qq(gsd_toler_results$p_lrt)
 dev.off()
 
 #significant SNPS above threshold on chromosome 9, 17, 29 -> saving these SNPs to table
 #chromosome 9 SNPS are VERY suggestive 
 gsd_sig_toler <- gsd_toler_results %>% 
   mutate(negLogP = -log10(p_lrt)) %>% 
   select(chr, rs, p_lrt, negLogP) %>% 
   filter(negLogP > 5)
 
 write_csv(gsd_sig_toler, "gsd_sig_toler.csv")
 
 #LR 
 system(str_c("./plink2 --bfile LR_afterQC --dog --pheno lr_phenotype.txt --pheno-name ", 'Tolerance_to_dogs', 
              ' --make-bed --out lr_gemma_toler'))
 #estimate a relatedness matrix with bed files 
 system("gemma -bfile lr_gemma_toler -gk 1 -o lr_RelMat_toler")  
 system("gemma -bfile lr_gemma_toler -k ./output/lr_RelMat_toler.cXX.txt -lmm 2 -o lr_toler_GWASresults.lmm")
 #read results in 
 lr_toler_results <- read_table("./output/lr_toler_GWASresults.lmm.assoc.txt")
 #use the Bonferroni threshold for lr 
 #manhattan plot 
 png("GWAS_lr_toler.png")
 manhattan(lr_toler_results,chr="chr",bp="ps",p="p_lrt",snp="rs",genomewideline=bonferroni_lr)
 dev.off()
 #qqplot
 png("GWAS_lr_toler_qq.png")
 qq(lr_toler_results$p_lrt)
 dev.off()
 
 #significant SNPS above threshold on chromosome 18 -> saving these SNPs to table
 lr_sig_toler <- lr_toler_results %>% 
   mutate(negLogP = -log10(p_lrt)) %>% 
   select(chr, rs, p_lrt, negLogP) %>% 
   filter(negLogP > 5)
 
 write_csv(lr_sig_toler, "lr_sig_toler.csv")
 
 
 
 #for Boldness: 
 #GS
 system(str_c("./plink2 --bfile GS_afterQC --dog --pheno gsd_phenotype.txt --pheno-name ", 'Boldness', 
              ' --make-bed --out gsd_gemma_bold'))
 #estimate a relatedness matrix with bed files 
 system("gemma -bfile gsd_gemma_bold -gk 1 -o gsd_RelMat_bold")  
 system("gemma -bfile gsd_gemma_bold -k ./output/gsd_RelMat_bold.cXX.txt -lmm 2 -o gsd_bold_GWASresults.lmm")
 
 #read results in 
 gsd_bold_results <- read_table("./output/gsd_bold_GWASresults.lmm.assoc.txt")
 #use the Bonferroni threshold for gsd 
 #manhattan plot 
 png("GWAS_gsd_bold.png")
 manhattan(gsd_bold_results,chr="chr",bp="ps",p="p_lrt",snp="rs",genomewideline=bonferroni_gsd)
 dev.off()
 #qqplot
 png("GWAS_gsd_bold_qq.png")
 qq(gsd_bold_results$p_lrt)
 dev.off()
 
 #no significant SNPs detected -> none above threshold 
 
 #LR 
 system(str_c("./plink2 --bfile LR_afterQC --dog --pheno lr_phenotype.txt --pheno-name ", 'Boldness', 
              ' --make-bed --out lr_gemma_bold'))
 #estimate a relatedness matrix with bed files 
 system("gemma -bfile lr_gemma_bold -gk 1 -o lr_RelMat_bold")  
 system("gemma -bfile lr_gemma_bold -k ./output/lr_RelMat_bold.cXX.txt -lmm 2 -o lr_bold_GWASresults.lmm")
 #read results in 
 lr_bold_results <- read_table("./output/lr_bold_GWASresults.lmm.assoc.txt")
 #use the Bonferroni threshold for lr 
 #manhattan plot 
 png("GWAS_lr_bold.png")
 manhattan(lr_bold_results,chr="chr",bp="ps",p="p_lrt",snp="rs",genomewideline=bonferroni_lr)
 dev.off()
 #qqplot
 png("GWAS_lr_bold_qq.png")
 qq(lr_bold_results$p_lrt)
 dev.off()
 
 #significant SNPS above threshold on chromosome 9 and 21 -> saving these SNPs to table
 lr_sig_bold <- lr_bold_results %>% 
   mutate(negLogP = -log10(p_lrt)) %>% 
   select(chr, rs, p_lrt, negLogP) %>% 
   filter(negLogP > 5)
 
 write_csv(lr_sig_bold, "lr_sig_bold.csv")
 
 
 
 #for Interest in dummy: 
 #GS
 system(str_c("./plink2 --bfile GS_afterQC --dog --pheno gsd_phenotype.txt --pheno-name ", 'Interest_in_the_dummy', 
              ' --make-bed --out gsd_gemma_interest'))
 #estimate a relatedness matrix with bed files 
 system("gemma -bfile gsd_gemma_interest -gk 1 -o gsd_RelMat_interest")  
 system("gemma -bfile gsd_gemma_interest -k ./output/gsd_RelMat_interest.cXX.txt -lmm 2 -o gsd_interest_GWASresults.lmm")
 
 #read results in 
 gsd_interest_results <- read_table("./output/gsd_interest_GWASresults.lmm.assoc.txt")
 #use the Bonferroni threshold for gsd 
 #manhattan plot 
 png("GWAS_gsd_interest.png")
 manhattan(gsd_interest_results,chr="chr",bp="ps",p="p_lrt",snp="rs",genomewideline=bonferroni_gsd)
 dev.off()
 #qqplot
 png("GWAS_gsd_interest_qq.png")
 qq(gsd_interest_results$p_lrt)
 dev.off()
 
 #significant SNPS above threshold on chromosome 24 -> saving these SNPs to table
 gsd_sig_interest <- gsd_interest_results %>% 
   mutate(negLogP = -log10(p_lrt)) %>% 
   select(chr, rs, p_lrt, negLogP) %>% 
   filter(negLogP > 5)
 
 write_csv(gsd_sig_interest, "gsd_sig_interest.csv")
 
 #LR 
 system(str_c("./plink2 --bfile LR_afterQC --dog --pheno lr_phenotype.txt --pheno-name ", 'Interest_in_the_dummy', 
              ' --make-bed --out lr_gemma_interest'))
 #estimate a relatedness matrix with bed files 
 system("gemma -bfile lr_gemma_interest -gk 1 -o lr_RelMat_interest")  
 system("gemma -bfile lr_gemma_interest -k ./output/lr_RelMat_interest.cXX.txt -lmm 2 -o lr_interest_GWASresults.lmm")
 #read results in 
 lr_interest_results <- read_table("./output/lr_interest_GWASresults.lmm.assoc.txt")
 #use the Bonferroni threshold for lr 
 #manhattan plot 
 png("GWAS_lr_interest.png")
 manhattan(lr_interest_results,chr="chr",bp="ps",p="p_lrt",snp="rs",genomewideline=bonferroni_lr)
 dev.off()
 #qqplot
 png("GWAS_lr_interest_qq.png")
 qq(lr_interest_results$p_lrt)
 dev.off()
 
 #significant SNPS above threshold on chromosome 23 -> saving these SNPs to table
 lr_sig_interest <- lr_interest_results %>% 
   mutate(negLogP = -log10(p_lrt)) %>% 
   select(chr, rs, p_lrt, negLogP) %>% 
   filter(negLogP > 5)
 
 write_csv(lr_sig_interest, "lr_sig_interest.csv")
 
 