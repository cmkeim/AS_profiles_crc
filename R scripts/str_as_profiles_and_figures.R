######################################################################
# This file provides basic Profile overview of the 
# AS and STR landscape in CRC
# 10.7.2022
#####################################################################
#Load the necessary packages
library(dplyr)
library(tidyr) 
library(ggplot2)

# Load the Mapped date file ()
str_in_as_merge_crc_overlap <- read.csv('C:/Users/cmk_8/OneDrive/Dokumente/Master ACLS/Master Thesis/Semester 6/data/str_in_as_merge_crc_overlap.tsv', header = T,  sep="\t")
str_in_as_merge_crc_overlap <- str_in_as_merge_crc_overlap  %>% 
  filter(!db_match == "1175191:CTTTCT,CTCTGT,CTCTGT,CTCTGT,CTCTCT,CTCTGT,CTTTCT:0.88:3") %>%
  filter(!db_match == "1176667:GT,GT,GT,GC,GT,GT,GT,CT,TT,GT,GT,GT,GT,GT,TT,GT,GT,GT,GT,GT,GT,GT,GT,GT,GT,GT,GT,GT,CT,GT:0.92:13") %>%
  filter(!db_match == ".:TCTCTG,TCTCTG,TCTCTG,TCTCTG,TCTCTG,TCTCTG:1.0:6")

# correct overlap cooridnates by +1
str_in_as_merge_crc_overlap["bp_overlap"] = str_in_as_merge_crc_overlap["bp_overlap"]+1

# Load the COC_PSI data (merged from COAD and READ)
crc_psi <- read.csv("C:/Users/cmk_8/OneDrive/Dokumente/Master ACLS/Master Thesis/Semester 6/data/PSI_COAD_and_READ.tsv", header = T,  sep="\t")
# read_psi <- read.csv('C:/Users/cmk_8/OneDrive/data/PSI_download_READ.txt', sep="\t")
# coad_psi <- read.csv("C:/Users/cmk_8/OneDrive/data/PSI_download_COAD.txt", sep="\t")
# =============================================================================================================
# ------------------  Overview Summary table for STR ASE  ---------------------------------------------
# =============================================================================================================
str_in_ase <- str_in_as_merge_crc_overlap %>% 
  mutate(ase_length = (end_as-start_as+1), str_length = end-start+1) %>%
  select(splice_type, symbol, as_id, str, repeat_length, ase_length, str_length) %>% 
  distinct() %>%  
  # filter(ase_length < 100) %>% 
  # filter(ase_length >= 100 & ase_length < 500) %>% 
  # filter(ase_length >= 500 & ase_length < 1000) %>% 
  # filter(ase_length >= 1000) %>% 
  group_by(splice_type) %>%
  summarise(count = n(), 
            av_ase_len = round(mean(ase_length)), ase_sd = round(sd(ase_length)), 
            av_str_len = round(mean(str_length),2), str_sd = round(sd(str_length),2)) %>%View()
  
# Count the STRs overall
global_most_freq_STR <- str_in_as_merge_crc_overlap %>% 
    filter(!splice_type == "ME") %>% # ME exon contains only 1 STR
    select(splice_type, end_as, start_as, end, start, str) %>% distinct() %>% 
    mutate(as_length = (end_as-start_as+1), 
           str_length = end-start+1) %>%
    group_by(splice_type,str) %>% 
    summarise(count = n(), 
              av_ase_len = round(mean(as_length)), ase_sd = round(sd(as_length)), 
              av_str_len = round(mean(str_length),2), str_sd = round(sd(str_length),2))  
write.table(global_most_freq_STR[order(global_most_freq_STR$count),], file="global_most_frequent_strs_ordered.csv", sep=",", quote=FALSE, row.names=FALSE)

# Same for the Fully contained  
fully_contained_strs <- df1 %>% 
  group_by(splice_type) %>%
  summarise(count = n(), 
            av_ase_len = round(mean(as_length)), ase_sd = round(sd(as_length)), 
            av_str_len = round(mean(str_length),2), str_sd = round(sd(str_length),2)) %>%View()
# count the strs in each event
most_frequent_strs <- df1 %>% 
  group_by(splice_type, str) %>%
  summarise(count = n(), 
            av_ase_len = round(mean(as_length)), ase_sd = round(sd(as_length)), 
            av_str_len = round(mean(str_length),2), str_sd = round(sd(str_length),2)) %>%
# write.table(most_frequent_strs, file="most_frequent_strs.csv", sep=",", quote=FALSE, row.names=FALSE)
# write.table(most_frequent_strs[order(most_frequent_strs$count),], file="most_frequent_strs_ordered.csv", sep=",", quote=FALSE, row.names=FALSE)

partially_overlapping_strs <- df4 %>% 
  group_by(splice_type) %>%
  summarise(count = n(), 
            av_ase_len = round(mean(as_length)), ase_sd = round(sd(as_length)), 
            av_str_len = round(mean(str_length),2), str_sd = round(sd(str_length),2)) %>%View()

most_frequent_strs_border <- df4 %>% 
  group_by(splice_type, str) %>%
  summarise(count = n(), 
            av_ase_len = round(mean(as_length)), ase_sd = round(sd(as_length)), 
            av_str_len = round(mean(str_length),2), str_sd = round(sd(str_length),2))
# write.table(most_frequent_strs_border[order(most_frequent_strs_border$count),], file="most_frequent_strs_border_ordered.csv", sep=",", quote=FALSE, row.names=FALSE)
# most_frequent_strs_border2 <- df4 %>% 
#   group_by(splice_type, str) %>%
#   summarise(str_count = n(), 
#             av_ase_len = round(mean(as_length)), ase_sd = round(sd(as_length)), 
#             av_str_len = round(mean(str_length),2), str_sd = round(sd(str_length),2)) %>%
#   group_by(splice_type) %>%
#   summarise(count = n()) %>% View()


# -------------------------------------------------------------------------------------------------
#◘ ASE overview across CRC
# -------------------------------------------------------------------------------------------------
total_ase = 35957 #crc_psi %>% select(as_id) %>% distinct()
str_ase = 7934    #str_in_as_merge_crc_overlap %>% select(as_id) %>% distinct() %>% View()
fully_contained = 7583 #df1 %>% select(as_id) %>% distinct() %>% View()
partially_overlapping = 554 # df4 %>% select(as_id) %>% distinct() %>% View()
fully_covering = 19 #str_in_as_merge_crc_overlap %>%  filter(!is.na(str_overlap_start) & !is.na(str_overlap_end)) %>% select(splice_type, symbol, as_id) %>% distinct()
filter(df4, as_id %in% df1$as_id) %>% distinct()  # 218


names = rbind("Total ASE in CRC", "Total STR containing ASE", "ASE with fcSTRs",
              "ASE with bcSTRs", "ASE fully covered by STR")
data = rbind(35957, 7934, 7583, 554, 19)
overview=data.frame(names, data)
ggplot(overview, aes(y=reorder(names, data), x=data))+
  geom_bar(stat="identity") +
  geom_text(aes(label=data), hjust=-0.15,  color="black", size=6.5)+
  theme_bw()+
  theme(
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.y=element_blank(),
    axis.ticks.x=element_blank(),
    axis.text.y=element_text(size = 30),
    legend.position="none",
    axis.title = element_blank(),
    axis.title.x = element_blank(),
    axis.text = element_blank(),
    strip.text = element_blank(),
    text = element_text(size = 30))+
  xlim(0,38000)+
  scale_fill_grey(start=0.8,end =0.8)
  
# =============================================================================================================
# ------------------ AS in CRC globally (w/o STR containing data) ---------------------------------------------
# =============================================================================================================

# How many READ and COAD patients are at hand
# read_p <- read_psi %>% select(., contains("TCGA"))
# read_pn <- read_psi %>% select(., contains("TCGA") & contains("Norm"))     
# coad_p <- coad_psi %>% select(., contains("TCGA"))
# coad_pn <- coad_psi %>% select(., contains("TCGA") & contains("Norm"))     

# How many patients and how many AS in genes overall in CRC
as_in_crc <- crc_psi %>% select(as_id) %>% distinct() # unique AS per gene
as_in_crc_pie <- as_in_crc  %>% group_by(splice_type) %>% summarise(count=n()) %>% mutate(perc=count/sum(count)*100)
as_in_crc_pie$data_set <- "all"

as_in_crc_pie %>%
  ggplot(aes(x=data_set, y=count, fill=splice_type)) +
  geom_bar(stat="identity") +
  coord_polar("y", start=0) +
  geom_col(color = "black") +
  theme_void() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        text = element_text(size = 24)) +
  labs(fill="Splice Type", title = "ASE in CRC by splice type") 

# as_in_crc_tot <- crc_psi %>% select(symbol, as_id) %>% unique()                    # number of genes with AS events
# write.table(as_in_crc$symbol, file="as_genes_coad_read.tsv", sep="\t", quote=FALSE, row.names=FALSE)


crc_plot <- crc_psi %>% select(splice_type, symbol, as_id) %>% distinct() %>% 
  group_by(symbol, splice_type) %>%
  summarise(as_count_gene = n(), .groups = "keep") %>%
  group_by(splice_type) %>%
  summarise(AS_Exons = sum(as_count_gene), AS_Genes = n()) %>%
  pivot_longer(-splice_type, names_to = "subgroup", values_to = "subgroup_count")


ggplot(crc_plot, aes(y=subgroup_count, x=splice_type, fill = subgroup)) + 
  geom_bar(position=position_dodge(width=0.55), width = 0.5, stat='identity') +
  theme_bw() + 
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(.92, .9),
        legend.title=element_blank(),
        legend.text = element_text(size=12),
        plot.title = element_text(hjust=0.5, size=12, face='bold'),
        axis.title = element_text(size=12, face = "bold"),
        axis.text = element_text(face = "bold")) +
  scale_y_continuous(breaks=seq(0,26000,by=2000)) +
  labs(x='Splice Types', y='Alternative splicing events', title='AS event counts per splice type and per gene') +
  scale_fill_brewer(palette="Paired")

# Figure X - Overview of AS types in CRC (without STRs)
# pie_crc <- crc_psi %>% select(-cancer_type) %>% distinct() %>%
#   group_by(splice_type) %>% summarise(count=n()) %>% mutate(perc=count/sum(count)*100)
# pie_crc$data_set <- "all"
# 
# pie_crc %>% ggplot(aes(x=data_set, y=count, fill=splice_type)) +
#   geom_bar(stat="identity") +
#   coord_polar("y", start=0) +
#   geom_col(color = "black") +
#   theme_void() +
#   theme(axis.text = element_blank(),
#         axis.ticks = element_blank(),
#         panel.grid = element_blank(),
#         axis.title = element_blank()) +
#   labs(fill="Splice Type", title = "ASE in CRC by splice type") 

# =============================================================================================================
# ----------------------------- STR containing AS events in CRC   ---------------------------------------------
# =============================================================================================================
# How many genes and AS are identified in CRC samples
# str_in_as_merge_crc_overlap %>% select(symbol) %>% distinct() 
# str_in_as_merge_crc_overlap %>% select(as_id, exon) %>% distinct() 
# amount of STRs detected: observations from df1 (fcSTRs) + df4 (bcSTRs)

# Count how many STR containing AS events occuring per spliceType and per Gene-symbol
str_in_as_plot <- str_in_as_merge_crc_overlap %>% 
  select(splice_type, symbol,as_id) %>% distinct() %>% 
  group_by(symbol, splice_type) %>%
  summarise(as_count_gene = n(), .groups = "keep") %>%
  group_by(splice_type) %>%
  summarise(AS_Exons = sum(as_count_gene), AS_Genes = n()) %>%
  pivot_longer(-splice_type, names_to = "subgroup_str", values_to = "subgroup_count_str") %>%
  select(-splice_type)

r = cbind(crc_plot, str_in_as_plot) %>% 
  select(splice_type, subgroup, subgroup_count, subgroup_str, subgroup_count_str) %>%
  mutate(overall = subgroup_count - subgroup_count_str, str_fraction = subgroup_count_str) %>% 
  select(-c("subgroup_count", "subgroup_count_str", "subgroup_str")) %>% 
  pivot_longer(-c("splice_type", "subgroup"), names_to = "subgroups", values_to = "subgroup_count")

# Figure 3 - STR containing AS events per alternatively spliced Exons (ASE) and alternatively spliced genes (ASG)
ggplot(r) + 
  geom_bar(aes(x=splice_type, y=subgroup_count, fill = subgroups), position  = 'stack',  stat='identity', width = 0.5)+  
  facet_wrap(~subgroup, nrow = 2)+
  theme_bw() + 
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(.93, .92),
        legend.title= element_text(size=14, face = "bold"),
        legend.text = element_text(size=14),
        # plot.title = element_text(hjust=0.5, size=12, face='bold'),
        axis.title = element_text(size=14, face = "bold"),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.text = element_text(size=14, color="black"),
        strip.text = element_text(size = 14, face = "bold")) +
  scale_y_continuous(breaks=seq(0,14000,by=2000))+
  labs(x='Splice Types', y='Alternative splicing events', fill="AS in CRC") +
  scale_fill_brewer(palette="Paired")



# =============================================================================================================
# ----------------------------- Density coverage of STRs in AS (how many of AS are STRs)  -------------------------------
# =============================================================================================================
df1 <- str_in_as_merge_crc_overlap %>%
  filter(is.na(str_overlap_start) & is.na(str_overlap_end), !splice_type == "ME") %>% # ME exon contains only 1 STR
  mutate(as_length = (end_as-start_as+1), 
         str_length = end-start+1, 
         max_str_start_loc = as_length-str_length+1, # max. starting location
         str_inAS_start_pos = (start-start_as+1),
         str_center_loc = (str_inAS_start_pos/max_str_start_loc)) %>%
  select(chromosome, str, repeat_length, str_length, exon, as_id, splice_type, bp_overlap, as_length, max_str_start_loc, str_inAS_start_pos, 
         start, start_as, str_center_loc, str_overlap_start, str_overlap_end,symbol) %>% distinct()

# Figure 5a - Location of STRs in AES by different splice types
ggplot(df1, aes(x = str_center_loc, fill = splice_type)) +
  geom_density(alpha = 0.7)+
  theme_bw() +
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank(),
        legend.position="none",
        axis.title = element_text(size=40, face = "bold"),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.text = element_text(size=24, color="black"),
        strip.text = element_text(size = 40, face = "bold"),
        text = element_text(size = 28)) +
  xlim(0,1)+
  labs(x='Normalized STR location', y='density')+
  facet_wrap(~splice_type)
  # facet_wrap(~repeat_length)


# ---------------------------------------------------------------------------------------
# Sum up length of all exons and STRs for all AS events
# Figure 4a - 'Total STR coverage in alternatively spliced Exons'
df2 <- df1 %>% group_by(exon, symbol, splice_type, repeat_length, as_id) %>% distinct() %>% 
  summarize(str_len_sum = sum(str_length), 
            as_len_sum = as_length, 
            coverage = str_len_sum/as_len_sum,.groups = "keep") 

ggplot(df2, aes(x = coverage, fill=splice_type)) +
  geom_density(alpha = 0.7)+
  theme_bw() +
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank(),
        legend.position = "none",
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title = element_text(size=40, face = "bold"),
        axis.text = element_text(size=24, color="black"),
        strip.text = element_text(size = 40, face = "bold"),
        text = element_text(size = 28)) +
  xlim(0,0.15)+
  labs(x='Total STR coverage over all ASE', y='density')+
  facet_wrap(~splice_type)

# Figure 4b - Pie chart of Total STR coverage in alternatively spliced Exons
pie_df2 <- df1 %>% group_by(splice_type) %>% summarise(count=n()) %>% mutate(perc=count/sum(count)*100)
pie_df2$data_set <- "all"

pie_df2 %>% ggplot(aes(x=data_set, y=count, fill=splice_type)) +
  geom_bar(stat="identity") +
  coord_polar("y", start=0) +
  geom_col(color = "black") +
  theme_void() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        text = element_text(size=24)) +
  labs(fill="Splice Type", title = "Total STR coverage in AES by splice type") 



ggplot(df2, aes(x=as_len_sum, y=str_len_sum, color=splice_type)) +
  geom_point()+
  ylim(0,7500)+
  geom_smooth(method="lm", se=T, level=0.95, aes(color="green"))
  facet_wrap(~splice_type)

plot(str_len_sum ~ as_len_sum, df2)
data.lm = lm(str_len_sum ~ as_len_sum, df2%>% filter(splice_type=="ES"))
 # --> H0 rejected that all groups are equal (p=0.015 < 0.05)
library(car)
Anova(data.lm, type = 2) # Unbalanced datasets (not completely randomized)
summary(data.lm)


# Paritally Overlap (bcSTRs)
df4 <- str_in_as_merge_crc_overlap %>%
  filter(!is.na(str_overlap_start) | !is.na(str_overlap_end), !splice_type == "ME", 
         is.na(str_overlap_start) | is.na(str_overlap_end)) %>% # ME exon contains only 1 STR
  mutate(as_length = (end_as-start_as+1), 
         str_length = end-start+1, 
         coverage = str_length/as_length,
         max_str_start_loc = case_when(!is.na(str_overlap_start) & is.na(str_overlap_end) ~ bp_overlap/as_length,
                                       is.na(str_overlap_start) & !is.na(str_overlap_end) ~ 1-bp_overlap/as_length))  %>%  # max. starting location
  select(str, coverage,repeat_length, str_length, exon, as_id, splice_type, bp_overlap, as_length, max_str_start_loc,
         start, start_as, str_overlap_start, str_overlap_end,symbol) %>% distinct()

# 
# df4<- df1 %>% filter(!is.na(str_overlap_start) | !is.na(str_overlap_end)) %>% 
#   group_by(as_id) %>%
#   summarize(str_len_sum = sum(bp_overlap), 
#             as_len_sum = sum(as_length), 
#             coverage = str_len_sum/as_len_sum, 
#             splice_type = paste(unique(splice_type)), .groups = "keep") 

# Figure 6a - STR coverage @ exon borders in ASE
df4 %>% filter(coverage < 0.99) %>% ggplot(aes(x = max_str_start_loc, fill=splice_type)) +
  geom_density(alpha = 0.7)+
  theme_bw() +
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank(),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title = element_text(size=40, face = "bold"),
        axis.text = element_text(size=24, color="black"),
        strip.text = element_text(size = 40, face = "bold"),
        text = element_text(size = 28)) +
  xlim(0,1)+
  labs(x='Normlized STR location', y='density')+
  facet_wrap(~splice_type)

# df4 %>% filter(coverage > 0.25) %>% ggplot(aes(x = max_str_start_loc, fill=splice_type)) +
#   geom_density(alpha = 0.7)+
#   theme_bw() +
#   theme(plot.background = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         legend.position = "none",
#         axis.ticks.y=element_blank(),
#         axis.text.y=element_blank(),
#         axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
#         axis.title = element_text(size=14, face = "bold"),
#         axis.text = element_text(size=14, color="black"),
#         strip.text = element_text(size = 14, face = "bold")) +
#   xlim(0,1)+
#   labs(x='Normlized STR location', y='density')+
#   facet_wrap(~splice_type)

# Figure 6b - Pie chart of STR border propensity in ASE
pie_df4 <- df4 %>% group_by(splice_type) %>% summarise(count=n()) %>% mutate(perc=count/sum(count)*100)
pie_df4$data_set <- "all"

pie_df4 %>% ggplot(aes(x=data_set, y=count, fill=splice_type)) +
  geom_bar(stat="identity") +
  coord_polar("y", start=0) +
  geom_col(color = "black") +
  theme_void() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        text = element_text(size=24)) +
  labs(fill="Splice Type", title = "STR coverage in AES @ Exon borders by splice type") 


# =============================================================================================================
# -----------------------------  Density coverage of STRs in AS by STR-Repeat-Length  -------------------------------
# =============================================================================================================

# Figure 7 - Overall STR coverage in ASE per Repeat_length
# --> shows similar pattern like for ~splice_type facet wrap
ggplot(df2, aes(x = coverage, fill=splice_type)) +
  geom_density(alpha = 0.7)+
  theme_bw() +
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank(),
        axis.title = element_text(size=24, face = "bold"),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.text = element_text(size=24, color="black"),
        strip.text = element_text(size = 24, face = "bold"),
        text = element_text(size = 24)) +
  xlim(0,0.2)+
  labs(x='Total STR coverage', y='density', fill="Splice Type")+
  facet_wrap(~repeat_length)

# Figure 8 - STR coverage of fully contained STRs in ASE per Repeat-length
df1 %>% group_by(splice_type, repeat_length) %>%
  mutate(count = n()) %>% ungroup() %>%
  filter(count >= 3) %>%
  ggplot(aes(x = str_center_loc, fill = splice_type)) +
  geom_density(alpha = 0.7)+
  theme_bw() +
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank(),
        axis.title = element_text(size=40, face = "bold"),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.text = element_text(size=24, color="black"),
        strip.text = element_text(size = 40, face = "bold"),
        text = element_text(size = 28)) +
  xlim(0,1)+
  labs(x='Normalized STR location', y='density', fill="Splice Type")+
  facet_wrap(~repeat_length)

# Select those STRs which are showing clear peaks in density plots
df1_table <- df1 %>% 
  group_by(splice_type, repeat_length) %>%
  mutate(count = n()) %>% ungroup() %>%
  filter(count >= 3) %>%
  filter(repeat_length == 6, splice_type == "AA") %>%
  select(str, repeat_length, str_length, exon, as_id, splice_type, bp_overlap, as_length, max_str_start_loc,
         start, start_as, str_overlap_start, str_overlap_end, symbol, count)


# Figure 9 - STR coverage @ exon borders in ASE per Repeat-length
# filter the df4 where less then 3 ASE are observed
df4 %>% group_by(splice_type, repeat_length) %>%
  mutate(count = n()) %>% ungroup() %>%
  filter(count >= 3) %>% 
  ggplot(aes(x = max_str_start_loc, fill=splice_type)) +
  geom_density(alpha = 0.7)+
  theme_bw() +
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank(),
        axis.title = element_text(size=40, face = "bold"),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.text = element_text(size=24, color="black"),
        strip.text = element_text(size = 40, face = "bold"),
        text = element_text(size = 40)) +
  xlim(0,1)+
  labs(x='Normalized STR location', y='density', fill="Splice Type")+
  facet_wrap(~repeat_length)

# Select those STRs which are showing clear peaks in density plots
df4_table <- df4 %>% 
  group_by(splice_type, repeat_length) %>%
  mutate(count = n()) %>% ungroup() %>%
  filter(count >= 3) %>%
  filter(repeat_length == 3 & splice_type == "AT" |
         repeat_length == 4 & splice_type == "RI" | 
         repeat_length == 6 & splice_type == "AP") %>%
  select(-coverage, str, repeat_length, str_length, exon, as_id, splice_type, bp_overlap, as_length, max_str_start_loc,
                            start, start_as, str_overlap_start, str_overlap_end, symbol, count)

# combine those selected STRs with a clear peak
repeat_overlap_table <- rbind(df1_table, df4_table)
#write.table(repeat_overlap_table, file="repeat_overlap_table.tsv", sep="\t", quote=FALSE, row.names=FALSE)

# =============================================================================================================
# -----------------------------  g:Profiler Process  -------------------------------
# =============================================================================================================

# Reformat the g:Profiler output --> change .csv to a .tsv file - just adapt the path and file-name for reformatting
# data <- read.csv('C:/Users/cmk_8/OneDrive/Dokumente/Master ACLS/Master Thesis/Semester 6/data/gProfiler_hsapiens_repeat_length_significant.csv', header = T,  sep=",")
# write.table(data, "gProfiler_hsapiens_repeat_length_significant.tsv", sep="\t", row.names=FALSE, quote=FALSE)

# Investigate Over-representation of top50 significantly STR-containing DEAS genes
data_df <- read.csv('C:/Users/cmk_8/OneDrive/Dokumente/Master ACLS/Master Thesis/Semester 6/data/gProfiler_top50_reformat_out.tsv', header = T,  sep="\t", quote="")
plot_df <- data_df %>% arrange(log2fc)
fold_levels <- factor(plot_df$term_name, ordered=TRUE)
plot_df$term_name <- factor(plot_df$term_name, levels=fold_levels)

plot_df %>% group_by(source) %>% slice_max(n=25, order_by=negative_log10_of_adjusted_p_value) %>% ungroup() %>% 
  ggplot(aes(y=reorder(term_name, log2fc), x=log2fc)) +
  geom_bar(stat="identity") +
  theme_classic() +
  xlab("log2(Fold change)") +
  theme(axis.title.y = element_blank(),
        axis.title = element_text(size=24, face = "bold"),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.text = element_text(size=24, color="black"),
        strip.text = element_text(size = 24, face = "bold"),
        text = element_text(size = 24)) +
  facet_wrap(~source, scales="free_y", nrow = 2)


# Investigate Over-representation of all significantly STR-containing DEAS genes
data_df1 <- read.csv('C:/Users/cmk_8/OneDrive/Dokumente/Master ACLS/Master Thesis/Semester 6/data/gProfiler_allDSG_woAmbiguous_reformat_out.tsv', header = T,  sep="\t", quote="")
plot_df1 <- data_df1 %>% arrange(log2fc)
fold_levels <- factor(plot_df1$term_name, ordered=TRUE)
plot_df1$term_name <- factor(plot_df1$term_name, levels=fold_levels)

plot_df1 %>% group_by(source) %>% slice_max(n=15, order_by=negative_log10_of_adjusted_p_value) %>% ungroup() %>% 
  ggplot(aes(y=reorder(term_name, log2fc), x=log2fc)) +
  geom_bar(stat="identity") +
  theme_classic() +
  xlab("log2(Fold change)") +
  theme(axis.title.y = element_blank(),
        axis.title = element_text(size=24, face = "bold"),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.text = element_text(size=24, color="black"),
        strip.text = element_text(size = 24, face = "bold"),
        text = element_text(size = 24)) +
  facet_wrap(~source, scales="free_y", nrow = 2)

data_df1 %>% select(-c(5:10)) %>%
  write.table("enrichment_table_all.csv", sep=",", quote=FALSE, row.names=FALSE)


# Investigate the g:Profiler over representation for STRs with high densities at the border
data_df2 <- read.csv('C:/Users/cmk_8/OneDrive/Dokumente/Master ACLS/Master Thesis/Semester 6/data/gProfiler_hsapiens_repeat_length_significant_out.tsv', header = T,  sep="\t", quote="")
plot_df2 <- data_df2 %>% arrange(log2fc)
fold_levels <- factor(plot_df2$term_name, ordered=TRUE)
plot_df2$term_name <- factor(plot_df2$term_name, levels=fold_levels)

plot_df2 %>% group_by(source) %>% slice_max(n=15, order_by=negative_log10_of_adjusted_p_value) %>% ungroup() %>% 
  ggplot(aes(y=reorder(term_name, log2fc), x=log2fc)) +
  geom_bar(stat="identity") +
  theme_classic() +
  xlab("log2(Fold change)") +
  theme(axis.title.y = element_blank(),
        axis.title = element_text(size=24, face = "bold"),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.text = element_text(size=24, color="black"),
        strip.text = element_text(size = 24, face = "bold"),
        text = element_text(size = 24)) +
  facet_wrap(~source, scales="free_y", nrow = 2)

# =============================================================================================================
# ----------------------------- Venn Diagram  -------------------------------
# =============================================================================================================

library("VennDiagram")
grid.newpage()                                        # Move to new plotting page
draw.triple.venn(area1 = 19,                          # Add name to each set
                 area2 = 554,
                 area3 = 7583,
                 n12 = 1,
                 n13 = 4,
                 n23 = 218,
                 n123 = 1,
                 fill = c("pink", "green","orange"),
                 col="black",
                 # category = c("fullyCovered", "bcSTRs", "fcSTRs"),
                 alpha=1,
                 cex=c(3))

# =============================================================================================================
# ----------------------------- STR loci Variation data   -------------------------------
# =============================================================================================================

# str_variation = read.csv('C:/Users/cmk_8/OneDrive/Dokumente/Master ACLS/Master Thesis/Semester 6/20220527_locus_variation_no_groups.csv', header=T, sep=",")





