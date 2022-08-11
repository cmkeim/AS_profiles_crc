########################################################################
# Investigate DSG (Different Spliced Genes acc to Liu et al., 2018)
# Chris Keim 1.7.2022
########################################################################


library(dplyr)
library(tidyr)

# Load the crc psi dataset (merged from COAD and READ)
crc_psi <- read.csv("C:/Users/cmk_8/OneDrive/Dokumente/Master ACLS/Master Thesis/Semester 6/data/PSI_COAD_and_READ.tsv", header = T,  sep="\t")

# ------------------------------------------------------------------------------------------------
# This attempt is calculating the p-value (t-test) for all patients (cancer, normal)
# for a splice event in a given gene (neglects the fact, that several STRs are found)
# this can be done via the "as_id" which gives each AS in a gene a specific identifier
# To make the code work, lines only containing NA's, must be removed first
str_in_as_merge_crc_overlap <- read.csv('C:/Users/cmk_8/OneDrive/Dokumente/Master ACLS/Master Thesis/Semester 6/data/str_in_as_merge_crc_overlap.tsv', header = T,  sep="\t")
str_in_as_merge_crc_overlap <- str_in_as_merge_crc_overlap  %>% 
  filter(!db_match == "1175191:CTTTCT,CTCTGT,CTCTGT,CTCTGT,CTCTCT,CTCTGT,CTTTCT:0.88:3") %>%
  filter(!db_match == "1176667:GT,GT,GT,GC,GT,GT,GT,CT,TT,GT,GT,GT,GT,GT,TT,GT,GT,GT,GT,GT,GT,GT,GT,GT,GT,GT,GT,GT,CT,GT:0.92:13") %>%
  filter(!db_match == ".:TCTCTG,TCTCTG,TCTCTG,TCTCTG,TCTCTG,TCTCTG:1.0:6")
# Remove all entries where the STR length (nucleoditdes) are bigger then 200
# because some STRs are wrongly liftovered & longest "real" detected STR is 114nt
# str_in_as_merge_crc_overlap <- str_in_as_merge_crc_overlap %>% mutate(str_len = end-start+1, as_length = end_as-start_as+1) %>% 
#   filter(!str_len > 200)

# Correct the bo_overlap by adding one
str_in_as_merge_crc_overlap["bp_overlap"] = str_in_as_merge_crc_overlap["bp_overlap"]+1

# Functions to Run t-test and handle the error while occuring when NA is in the data
run_test <- function(x, y) {
  res <- try(t.test(x, y)$p.value)
  
  if(grepl(pattern = "Error", x = res)){
    return(NA)
  } else {
    res
  }
}

# Same as above, just to ad DF --> no more necessary
run_test2 <- function(x, y) {
  res <- try(t.test(x, y)$parameter)
  
  if(grepl(pattern = "Error", x = res)){
    return(NA)
  } else {
    res
  }
}

# COmpare All cancer tisues (#648 cases) against all Normal tissues (#51)
# and calculate its p-values (t-test)
df_t = str_in_as_merge_crc_overlap %>% select(contains(c("as_id","cancer_type", "symbol", "TCGA"))) %>% distinct() %>%
  rowwise() %>%
  mutate(values_n = list(c_across(contains("TCGA") & contains("Norm"))),
         values_c = list(c_across(contains("TCGA") & !contains("Norm")))) %>%
  select(as_id, cancer_type, symbol, values_n, values_c) %>%
  mutate(p_val = run_test(unlist(values_n), unlist(values_c)), 
         DF = run_test2(unlist(values_n), unlist(values_c))) %>%
  filter(!is.na(p_val))

df_t$p_adj <- p.adjust(df_t$p_val, "BH") 
df_t_fdr <- df_t %>% filter(p_adj <= 0.05)


# Trial to calculate the Log2FC from the beginning on
# df_t2 = str_in_as_merge_crc_overlap %>% select(contains(c("as_id","cancer_type", "symbol", "TCGA"))) %>% distinct() %>%
#   mutate(psi_mean_c = rowMeans(select(., contains("TCGA") & !contains("Norm")), na.rm=T), 
#          psi_mean_n = rowMeans(select(., contains("TCGA") & contains("Norm")), na.rm=T),
#          log2FC = log2(psi_mean_c/psi_mean_n)) %>% 
#   rowwise() %>%
#   mutate(values_n = list(c_across(contains("TCGA") & contains("Norm"))),
#          values_c = list(c_across(contains("TCGA") & !contains("Norm")))) %>%
#   select(as_id, cancer_type, symbol, values_n, values_c, log2FC) %>%
#   mutate(p_val = run_test(unlist(values_n), unlist(values_c)), 
#          DF = run_test2(unlist(values_n), unlist(values_c))) %>%
#   filter(!is.na(p_val))
#
# test <- df_t_fdr %>% mutate(mean_c = mean(unlist(values_c),na.rm =TRUE),
#                             mean_n = mean(unlist(values_n), na.rm = TRUE),
#                             log2FC = log2(mean_c/mean_n))
# test$p_adj <- p.adjust(test$p_val, "BH")
# test <- df_t2 %>% filter(p_adj <= 0.05) %>% filter(log2FC > 1 | log2FC < -1, log2FC != Inf)

# -------------------------------------------------------------------------------------
# Filter original df "str_in_as_merge_crc_overlap" to the found signficant
# p-values (adjusted) --> mapping via "as_id"

df_t_sig <- str_in_as_merge_crc_overlap %>% right_join(df_t_fdr, by= c("as_id", "cancer_type","symbol")) %>% distinct() 

# df_t_sig %>% distinct() %>%  select(splice_type) %>% group_by(splice_type) %>% summarize(count=n())
# df_t_sig %>% select(symbol) %>% unique()

# -------------------------------------------------------------------------------------
# filter out the most top 50 AS events according to p-values adjusted
# -------------------------------------------------------------------------------------
get_top_50 <- function(x) {
  res_top = x[order(x$p_adj),][1:50,]
  return(res_top)
}
df_t_top = get_top_50(df_t_sig)

# Here, for the attempt where log2FC was directly calculated in lines 106-123
# get_top_50log <- function(x) {
#   res_top <- x %>% arrange(desc(log2FC)) %>% head(25)
#   res_bot <- x %>% arrange(desc(log2FC)) %>% tail(25)
#   res = rbind(res_top, res_bot)
#   return(res)
# }
# test_top =  get_top_50log(test_sig)


# ----------------------------------------------------------------------------------------
# Create a Overview table with top50 genes with STRs, AS, and additional info
#------------------------------------------------------------------------------------------
crc_psi_top50 <- filter(crc_psi, as_id %in% unique(df_t_top$as_id)) 

crc_psi_temp <- df_t_top %>% select(-c(24:702), -X,cancer_type) %>% 
  right_join(crc_psi_top50, by= c("as_id",  "symbol", "cancer_type")) %>%
  mutate(psi_mean_c = rowMeans(select(., contains("TCGA") & !contains("Norm")), na.rm=T), 
         psi_mean_n = rowMeans(select(., contains("TCGA") & contains("Norm")), na.rm=T),
         log2FC = log2(psi_mean_c/psi_mean_n))

# READ types will be removed, because if look at top50 most significant splice events, non is found in READ
# crc_psi_table <- crc_psi_temp %>%  select(symbol, splice_type.x, exons, str, str_len, as_length, p_adj, cancer_type, psi_mean_c, psi_mean_n, log2FC) %>%
#   filter(cancer_type != "READ")
# write.table(crc_psi_table, file="crc_psi_table.csv", sep=",", quote=FALSE, row.names=FALSE)
# crc_psi_all <- filter(crc_psi, as_id %in% unique(df_t_sig$as_id))

#------------------------------------------------------------------------------------------------------
# Identify no the significant STR-containing differently expressed AS genes 
# that are highly enrichted at the borders of ASE
#------------------------------------------------------------------------------------------------------

repeat_overlap_table <- read.csv("C:/Users/cmk_8/OneDrive/Dokumente/Master ACLS/Master Thesis/Semester 6/data/repeat_overlap_table.tsv", header = T,  sep="\t")
repeat_table_sig <- filter(repeat_overlap_table, symbol %in% unique(df_t_sig$symbol))
# Write the repeat_table_sig to a file with unique gene symbol names, as an import file for g:Profiler over representation
#write.table(unique(repeat_table_sig$symbol), file="sig_repeat_table.csv", sep=",", quote=FALSE, row.names=FALSE)
# repeat_table_sig_export <- repeat_table_sig %>% select(1:8, 12:15)
# write.table(repeat_table_sig_export, file="repeat_table_sig.csv", sep=",", quote=FALSE, row.names=FALSE)


str_dsg <- filter(str_in_as_merge_crc_overlap, as_id %in% unique(df_t_sig$as_id))
str_dsg_top50 <- filter(str_in_as_merge_crc_overlap, as_id %in% unique(df_t_top$as_id))
# str_dsg_top500 <- filter(str_in_as_merge_crc_overlap, as_id %in% unique(df_t_top500$as_id))
# str_dsg_log2FC <- filter(str_in_as_merge_crc_overlap, as_id %in% unique(df_t2_top$as_id))
# str_dsg_log2FC_all <- filter(str_in_as_merge_crc_overlap, as_id %in% unique(df_t2_sig$as_id))


# str_dsg_unique <- unique(str_dsg$as_id)
# str_dsg_top = filter(str_in_as_merge_crc_overlap, as_id %in% unique(df_t_top$as_id))
# str_dsg_symbol <- filter(str_in_as_merge_crc_overlap, symbol %in% unique(df_t_top$symbol))
# sort(unique(str_dsg_top50$symbol)) #sort(unique(str_dsg_paired_top$symbol)) 
# sort(unique(str_dsg$symbol))  #sort(unique(str_dsg_paired$symbol)) 
# sort(unique(str_dsg_top500$symbol))  #sort(unique(str_dsg_paired$symbol)) 
# sort(unique(str_dsg_log2FC_all$symbol)) 
# # sort(unique(str_dsg_symbol$symbol)); sort(unique(str_dsg_paired_symbol$symbol)) 
# write.table(unique(str_dsg_log2FC$symbol), file="str_dsg_log2FC.tsv", sep="\t", quote=FALSE, row.names=FALSE)
# write.table(unique(str_dsg_log2FC_all$symbol), file="str_dsg_log2FC_all.tsv", sep="\t", quote=FALSE, row.names=FALSE)

# ---------------------------------------------------------------------------------------
# check top sign. regulated genes  For all cancer patients
# While adding "cancer_type to merge column, duplicate entries of READ will be removed
# ---------------------------------------------------------------------------------------

df_top4 <- str_dsg %>% select(contains(c("splice_type","symbol", "TCGA"))) %>%
  mutate(psi_mean_c = rowMeans(select(., contains("TCGA") & !contains("Norm")), na.rm=T), 
         psi_mean_n = rowMeans(select(., contains("TCGA") & contains("Norm")), na.rm=T)) %>%
  group_by(splice_type, symbol) %>%
  summarise(cancer = mean(psi_mean_c), normal = mean(psi_mean_n), 
            sd_psi_c = sd(psi_mean_c), sd_psi_n = sd(psi_mean_n), .groups="keep") %>% 
  pivot_longer(., cols = c(cancer, normal), names_to = "tissue", values_to = "mean_psi") %>%
  group_by(splice_type, tissue) %>%
  summarise(Mean = mean(mean_psi),
            Nb = n()) %>%
  rowwise() %>%
  mutate(Labels = as.character(splice_type)) %>%
  mutate(Labels = paste0(Labels,"\nn = ",Nb))

# Same data but with Standard error of the mean
# df_top5 <- crc_psi_temp %>% filter(cancer_type != "READ") %>%
#   select(contains(c("splice_type.x","symbol", "psi_mean"))) %>% 
#   group_by(splice_type.x, symbol) %>% 
#   summarise(cancer = mean(psi_mean_c), normal = mean(psi_mean_n), 
#             sd_psi_c = sd(psi_mean_c), sd_psi_n = sd(psi_mean_n), .groups="keep") %>%
#   pivot_longer(., cols = c(cancer, normal), names_to = "tissue", values_to = "mean_psi") %>% 
#   group_by(splice_type.x, tissue) %>%
#   summarise(Mean = mean(mean_psi),
#             SEM = sd(mean_psi) / sqrt(n()),
#             Nb = n()) 


library(ggplot2)
ggplot(df_top4, aes(fill=tissue, y=Mean, x=splice_type)) + 
  geom_bar(position=position_dodge(width=0.55), width = 0.5, stat='identity') +
  # geom_errorbar(aes(ymin = Mean-SEM, ymax = Mean+SEM), width = 0.2, position=position_dodge(0.55)) +
  theme_bw() + 
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(.90, .85),
        legend.title = element_text(size=24),
        legend.text = element_text(size=24),
        # plot.title = element_text(hjust=0.5, size=14, face='bold'),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title = element_text(size=24, face = "bold"),
        axis.text = element_text(face = "bold"),
        strip.text = element_text(size = 24, face = "bold"),
        text = element_text(size=24)) +
  ylim(0,1)+
  scale_fill_manual('Tissue', values=c('red', 'blue')) +
  labs(x='Splice Types', y='Average Alternative Splicing', fill="Tissue")

## --------------------------------
## Overall STR count in CRC per Splice type and repeat length
# df_temp_overll <- str_dsg_top %>% select(splice_type, repeat_length) %>%
#   group_by(splice_type, repeat_length) %>%
#   summarise(str_count = n(), .groups="keep")
# df_temp_paired <- str_dsg_paired_top %>% select(splice_type, repeat_length) %>%
#   group_by(splice_type, repeat_length) %>%
#   summarise(str_count = n(), .groups="keep")
# df_temp_paired2 <- str_dsg_paired %>% select(splice_type, repeat_length) %>%
#   group_by(splice_type, repeat_length) %>%
#   summarise(str_count = n(), .groups="keep")
top50 <- str_dsg_top50 %>% select(splice_type, repeat_length) %>%
  group_by(splice_type, repeat_length) %>%
  summarise(str_count = n(), .groups="keep")


ggplot(top50, aes(y=str_count, x=repeat_length, fill=splice_type)) + 
  geom_bar(position=position_dodge(width=0.55), width = 0.5, stat='identity') +
  # facet_wrap(~cancer_type, scales="free")+
  theme_bw() +
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(.92, .75),
        legend.title=element_blank(),
        legend.text = element_text(size=10, face="bold"),
        plot.title = element_text(hjust=0.5, size=10, face='bold'),
        axis.title = element_text(size=10, face = "bold"),
        axis.text = element_text(face = "bold")) +
  # scale_y_continuous(breaks=seq(0,30,by=5))+
  scale_x_continuous(breaks=seq(1,7,by=1))+
  labs(x='Repeat Length', y='STR count', title='Overall STR count in CRC per Splice Type and STR repeat length')+
  scale_fill_manual(values=c("#66CCFF", "#33CC99",
                             "#3333FF", "#FFCC33",
                             "#CC3333", "#FF6666","#FF6666"))




































