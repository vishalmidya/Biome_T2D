


# Negative mode

bwqs_neg_pfas_met <-   read.csv("C:/Users/midyav01/OneDrive - The Mount Sinai Hospital/MSSM Projects/BIOME-PFAS pilot/T2D/Paper T2D/neg/negmode_meta_vs_pfas_bwqs.csv")
d_lm_status <- read.csv("C:/Users/midyav01/OneDrive - The Mount Sinai Hospital/MSSM Projects/BIOME-PFAS pilot/T2D/Paper T2D/neg/mixedmodel_exwas_t2d_meta_all.csv")

bb <- bwqs_neg_pfas_met[bwqs_neg_pfas_met$p.value < 0.05,c("refmet_name")]
tt <- d_lm_status[d_lm_status$p.value < 0.05, c("refmet_name")]
cc <- intersect(bb,tt) 

ft <-bwqs_neg_pfas_met[bwqs_neg_pfas_met$refmet_name %in% cc,c("refmet_name","mean","p.value","q.value")]
ft <- ft[order(ft$refmet_name),]

dt <- d_lm_status[d_lm_status$refmet_name %in% cc,c("refmet_name","Value","p.value","q.value")]
dt <- dt[order(dt$refmet_name),]

common_tab <- data.frame(refmet_name = dt[,"refmet_name"], mean.pfas.met = ft$mean, pval.pfas.met = ft$p.value, 
                         qval.pfas.met = ft$q.value, mean.met.t2d = dt$Value, pval.met.t2d = dt$p.value,
                         qval.met.t2d = dt$q.value)

common_tab

#               refmet_name mean.pfas.met pval.pfas.met qval.pfas.met mean.met.t2d pval.met.t2d qval.met.t2d
# 1    5-Hydroxy-tryptophan     0.3394450  7.232089e-04  1.062669e-02   -0.3158309  0.031736747    0.7303621
# 2          Glucoheptulose     0.2048479  1.256501e-02  8.180441e-02   -0.3400293  0.025262938    0.7303621
# 3 Sulfolithocholylglycine     0.5170315  4.024093e-09  1.142842e-06   -0.4648498  0.001309185    0.3718085



# Negative

em_comp_pfas_meta <- read.delim("C:/Users/midyav01/OneDrive - The Mount Sinai Hospital/MSSM Projects/BIOME-PFAS pilot/T2D/Paper T2D/neg/result_pfas_meta/tables/ListOfEmpiricalCompounds.tsv")
colnames(em_comp_pfas_meta)[4:5] <- c("compounds_mummichog","compound_names_mummichog")


em_comp_meta_t2d <- read.delim("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\BIOME-PFAS pilot\\T2D\\Paper T2D\\neg\\result_meta_t2d\\tables\\ListOfEmpiricalCompounds.tsv")
colnames(em_comp_meta_t2d)[4:5] <- c("compounds_mummichog","compound_names_mummichog")


em_pathways_pfas_meta  <- read.delim("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\BIOME-PFAS pilot\\T2D\\Paper T2D\\neg\\result_pfas_meta\\tables\\mcg_pathwayanalysis_.tsv")
colnames(em_pathways_pfas_meta)[5:7] <- c("overlap_EID","overlap_features_ID","overlap_features_names")

em_pathways_meta_t2d  <- read.delim("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\BIOME-PFAS pilot\\T2D\\Paper T2D\\neg\\result_meta_t2d\\tables\\mcg_pathwayanalysis_.tsv")
colnames(em_pathways_meta_t2d)[5:7] <- c("overlap_EID","overlap_features_ID","overlap_features_names")


cutoff <- 0.05

common_paths_neg <- unique(c(em_pathways_pfas_meta$pathway[em_pathways_pfas_meta$p.value < cutoff], 
                             em_pathways_meta_t2d$pathway[em_pathways_meta_t2d$p.value < cutoff]))

common_paths_neg

intersect(common_paths_pos,common_paths_neg)



