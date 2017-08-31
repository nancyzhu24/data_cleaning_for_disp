hcopen_pool <- dbPool(drv = "PostgreSQL",
                      host = "shiny.hc.local",
                      dbname = "hcopen",
                      user = "hcreader",
                      password = "canada1")
hcopen <- src_pool(hcopen_pool)


#date on which the master table is created:
data_date   <- "20160630" 

# removing these columns first results in a faster join
cv_bcpnn <- hcopen %>% tbl("PT_IC_161027") %>%
  select(drug_code, event_effect, count, expected_count, median_IC, LB95_IC, UB95_IC)
cv_gps <- hcopen %>% tbl("PT_GPS_161121") %>%
  select(drug_code, event_effect, postH0, postE_lambda, Q0.05_lambda)
cv_rfet <- hcopen %>% tbl("PT_RFET_161101") %>%
  select(drug_code, event_effect, midRFET, RFET)%>%as.data.frame()
cv_prr <- hcopen %>% tbl("PT_PRR_160927") %>%
  select(drug_code, event_effect, PRR, LB95_PRR, UB95_PRR)
cv_ror <- hcopen %>% tbl("PT_ROR_160927") %>%
  select(drug_code, event_effect, ROR, LB95_ROR, UB95_ROR)
cv_rrr <- hcopen %>% tbl("PT_RRR_160927") %>%
  select(drug_code, event_effect, RRR, LB95_RRR, UB95_RRR)
cv_counts <- hcopen %>% tbl("PT_counts_161103") %>%
  select(drug_code, event_effect, drug_margin, event_margin)


# Grab PT, HLT and SOC Terms from MedDra
soc_all <- hcopen %>% tbl("meddra") %>%
  select(PT_Term,HLT_Term,SOC_Term)
cv_soc <- soc_all %>%
  select(PT_Term,SOC_Term) %>%
  rename(event_effect = PT_Term) %>%
  as.data.frame()
cv_soc <- cv_soc[!duplicated(cv_soc),]

# Concatenate all SOC Terms
cat_pt <- cv_soc %>%
  group_by(event_effect) %>%
  summarize(SOC_Term = paste(unique(SOC_Term), collapse = ", "))


# Grab PT, HLT, quarter info
reportid_all <- hcopen %>% tbl("cv_drug_rxn_meddra") %>%
  select(PT_NAME_ENG, HLT_NAME_ENG,quarter)
cv_reportid <- reportid_all %>%
  select(PT_NAME_ENG, quarter) %>%
  rename(event_effect = PT_NAME_ENG) %>%
  as.data.frame()
cv_reportid <- cv_reportid[!duplicated(cv_reportid),]

# OPTION1: save var report id(a lot bigger)
# a <- cv_reportid %>%
#  group_by(event_effect,quarter) %>%
#  summarize(REPORT_ID = paste(unique(REPORT_ID), collapse = ", "))
# OPTION2: save var quarter
cat_id_pt <- cv_reportid %>%
  group_by(event_effect) %>%
  summarize(quarter = paste(unique(quarter), collapse = ", "))

## question: what is cat_id_pt for? Not seen anywhere in disp app, also summarize pt count from 1972 instead of 2006??????
master_table_pt<- cv_bcpnn %>%
  inner_join(cv_counts, by = c("drug_code", "event_effect")) %>%
  inner_join(cv_gps, by = c("drug_code", "event_effect")) %>%
  inner_join(cv_rfet, by = c("drug_code", "event_effect")) %>%
  inner_join(cv_prr, by = c("drug_code", "event_effect")) %>%
  inner_join(cv_ror, by = c("drug_code", "event_effect")) %>%
  inner_join(cv_rrr, by = c("drug_code", "event_effect")) %>%
  left_join(cat_pt,"event_effect" = "event_effect",copy = TRUE) %>%
  left_join(cat_id_pt, "event_effect" = "event_effect", copy = TRUE) %>%
  select(1:4, 24:25, 8:9, 5:23)%>%as.data.frame()
# already pull entire table to display at onset of app, might as well do everything locally (faster)

master_table_pt <- master_table_pt[!duplicated(master_table_pt),]


hlt_bcpnn <- hcopen %>% tbl("HLT_IC_161027") %>%
  select(drug_code, event_effect, count, expected_count, median_IC, LB95_IC, UB95_IC)
hlt_gps <- hcopen %>% tbl("HLT_GPS_161121") %>%
  select(drug_code, event_effect, postH0, postE_lambda, Q0.05_lambda)
hlt_rfet <- hcopen %>% tbl("HLT_RFET_161101") %>%
  select(drug_code, event_effect, midRFET, RFET)
hlt_prr <- hcopen %>% tbl("HLT_PRR_160927") %>%
  select(drug_code, event_effect, PRR, LB95_PRR, UB95_PRR)
hlt_ror <- hcopen %>% tbl("HLT_ROR_160927") %>%
  select(drug_code, event_effect, ROR, LB95_ROR, UB95_ROR)
hlt_rrr <- hcopen %>% tbl("HLT_RRR_160927") %>%
  select(drug_code, event_effect, RRR, LB95_RRR, UB95_RRR)
hlt_counts <- hcopen %>% tbl("HLT_counts_161103") %>%
  select(drug_code, event_effect, drug_margin, event_margin)
hlt_soc <- soc_all %>%
  select(HLT_Term,SOC_Term) %>%
  rename(event_effect = HLT_Term) %>%
  as.data.frame()
hlt_soc <- hlt_soc[!duplicated(hlt_soc),]
cat_hlt <- hlt_soc %>%
  group_by(event_effect) %>%
  summarize(SOC_Term = paste(unique(SOC_Term),collapse = ", "))

hlt_reportid <- reportid_all %>%
  select(HLT_NAME_ENG, quarter) %>%
  rename(event_effect = HLT_NAME_ENG) %>%
  as.data.frame()
hlt_reportid <- hlt_reportid[!duplicated(hlt_reportid),]

# OPTION1: save var report id(a lot bigger)
# new_id_hlt <- hlt_reportid %>%
#  group_by(event_effect,quarter) %>%
#  summarize(REPORT_ID = paste(unique(REPORT_ID), collapse = ", "))
# OPTION2: save var quarter
cat_id_hlt <- hlt_reportid %>%
  group_by(event_effect) %>%
  summarize(quarter = paste(unique(quarter), collapse = ", "))

master_table_hlt <- hlt_bcpnn %>%
  inner_join(hlt_counts, by = c("drug_code", "event_effect")) %>%
  inner_join(hlt_gps, by = c("drug_code", "event_effect")) %>%
  inner_join(hlt_rfet, by = c("drug_code", "event_effect")) %>%
  inner_join(hlt_prr, by = c("drug_code", "event_effect")) %>%
  inner_join(hlt_ror, by = c("drug_code", "event_effect")) %>%
  inner_join(hlt_rrr, by = c("drug_code", "event_effect")) %>%
  left_join(cat_hlt, "event_effect" = "event_effect", copy = TRUE) %>%
  left_join(cat_id_hlt, "event_effect" = "event_effect", copy = TRUE) %>%
  select(1:4, 24:25, 8:9, 5:23) %>%
  as.data.frame()

master_table_hlt <- master_table_hlt[!duplicated(master_table_hlt),]

#assign name with date for uploading to PostgreSQL
master_table_pt_name<-paste0("master_table_pt_",data_date)
master_table_hlt_name<-paste0("master_table_hlt_",data_date)

testconnect <- dbConnect(drv = "PostgreSQL", host = "shiny.hc.local", user = "hcwriter", dbname = "hcopen", password = "canada2")
dbWriteTable(testconnect,master_table_pt_name,master_table_pt, row.names = FALSE)
dbWriteTable(testconnect,master_table_hlt_name,master_table_hlt, row.names = FALSE)