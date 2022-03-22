
# upload datasets:
moffitt_df      <- read.table(...)
collisson_df    <- read.table(...)
bailey_df       <- read.table(...)
icgc_array_df   <- read.table(...)
badea_df        <- read.table(...)
tcga_df         <- read.table(...)
yang_df         <- read.table(...)
sandhu_df       <- read.table(...)

# add column with batch 
moffitt_df$batch      <- rep(1, nrow(moffitt_df     ))
collisson_df$batch    <- rep(2, nrow(collisson_df   ))
bailey_df$batch       <- rep(3, nrow(bailey_df      ))
icgc_array_df$batch   <- rep(4, nrow(icgc_array_df  ))
badea_df$batch        <- rep(5, nrow(badea_df       ))
tcga_df$batch         <- rep(6, nrow(tcga_df        ))
yang_df$batch         <- rep(7, nrow(yang_df        ))
sandhu_df$batch       <- rep(8, nrow(sandhu_df     ))

# traspose to have genes on the rows
moffitt_df      <- as.data.frame(t(moffitt_df))
collisson_df    <- as.data.frame(t(collisson_df))
bailey_df       <- as.data.frame(t(bailey_df))
icgc_array_df   <- as.data.frame(t(icgc_array_df))
badea_df        <- as.data.frame(t(badea_df))
tcga_df         <- as.data.frame(t(tcga_df))
yang_df         <- as.data.frame(t(yang_df))
sandhu_df       <- as.data.frame(t(sandhu_df))

moffitt_df$genes      <- rownames(moffitt_df)      
collisson_df$genes    <- rownames(collisson_df)    
bailey_df$genes       <- rownames(bailey_df)      
icgc_array_df$genes   <- rownames(icgc_array_df)  
badea_df$genes        <- rownames(badea_df)       
tcga_df$genes         <- rownames(tcga_df)        
yang_df$genes         <- rownames(yang_df)     
sandhu_df$genes       <- rownames(sandhu_df)      

# combine all the datasets
datasets = list(moffitt_df, collisson_df, bailey_df, icgc_array_df, badea_df, tcga_df, yang_df, sandhu_df)
require(plyr)
df_all <- join_all(datasets, by = 'genes', type="inner")
rownames(df_all) <- df_all$genes
df_all$genes <- NULL
df_all <- as.data.frame(t(df_all))

# save data NOT corrected
write.csv(df_all, "/nfs/home/users/mlautizi/PDAC/data/all_datasets_notCorrected.csv")


# remove batch effect
batch <- df_all$batch
df_all$batch <- NULL
df2 <- as.data.frame(sapply(df_all, function(x) as.numeric(as.character(x))))
rownames(df2) <- rownames(df_all)
df2 <- as.data.frame(t(df2))
df_all_corr <- limma::removeBatchEffect(df2, batch)
df_all_corr <- as.data.frame(t(df_all_corr))
df_all_corr$batch <- batch
# save dataset corrected
write.csv(df_all_corr, "/nfs/home/users/mlautizi/PDAC/data/all_datasets_batchCorrected.csv")
