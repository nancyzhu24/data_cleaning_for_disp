#MAP drug table to cv_drug_report to recreate cv_drug_rxn table:
library(lubridate)

hcopen<- dbPool(drv=RPostgreSQL::PostgreSQL(),
                host = "shiny.hc.local",
                dbname = "hcopen",
                user = "hcreader",
                password = "canada1")


cv_report<-tbl(hcopen,'cv_reports_20160630')%>%select(REPORT_ID,DATINTRECEIVED_CLEAN,ends_with('ENG'),AGE_GROUP_CLEAN)%>%as.data.table()
cv_reaction<-tbl(hcopen,'cv_reactions_20160630')%>%select(REPORT_ID,PT_NAME_ENG,SOC_NAME_ENG,MEDDRA_VERSION)%>%as.data.table()
cv_drug_report<-tbl(hcopen,'cv_report_drug_20160630')%>%
                filter(DRUGINVOLV_ENG=="Suspect")%>% #use only suspected drug in report
                select(REPORT_ID,DRUGNAME)%>%
                as.data.table()

cv_report%<>%select(-AGE_UNIT_ENG,-WEIGHT_UNIT_ENG,-HEIGHT_UNIT_ENG)

cv_drug_rxn_new<-left_join(cv_report,cv_reaction)%>%
  left_join(cv_drug_report%>%filter())%>%
  left_join(final_list)%>%
  filter(!is.na(ing))%>% #data loss due to unmatched drugnames Eg: DRUGNAME:PRIVIGEN
  mutate(month=floor_date(DATINTRECEIVED_CLEAN,unit="month"))%>%
  mutate(quarter=quarter(DATINTRECEIVED_CLEAN,with_year=TRUE))%>%
  select(REPORT_ID,DRUGNAME,ing,DATINTRECEIVED_CLEAN,REPORT_TYPE_ENG:MEDDRA_VERSION,quarter,month,category)

rm(cv_reaction,cv_report,cv_drug_report)
# testconnect <- dbConnect(drv = "PostgreSQL", host = "shiny.hc.local", user = "hcwriter", dbname = "hcopen", password = "canada2")
# dbWriteTable(testconnect, "cv_drug_rxn_new",cv_drug_rxn_new, row.names = FALSE)
# dbDisconnect(testconnect)

cv_drug_rxn_new<-tbl(hcopen,'cv_drug_rxn_new')%>%as.data.table()
meddra<-tbl(hcopen,'meddra')%>%select(PT_Term,HLT_Term)%>%distinct()%>%as.data.table()
cv_drug_rxn_meddra_new<-left_join(cv_drug_rxn_new,meddra,by=c('PT_NAME_ENG'='PT_Term'))

colnames(cv_drug_rxn_meddra_new)[18]<-"HLT_NAME_ENG"


#Calculate for mastertable:
source("~/CV_app_Git/cvapps/data_processing/stats_functions.R")

#refer to main.calcs document for Baysiean analysis:

cv_drug_rxn_new_2006 <- cv_drug_rxn_meddra_new %>%
  filter(!is.na(PT_NAME_ENG)) %>%
  filter(quarter >= 2006.1)
count_df <- cv_drug_rxn_new_2006 %>%
  group_by(ing, PT_NAME_ENG) %>%    # coerce to data frame to drop grouped and tbl attributes
  dplyr::summarise(count = n_distinct(REPORT_ID)) %>% as.data.frame()

count_df_hlt <- cv_drug_rxn_new_2006 %>%
  group_by(ing, HLT_NAME_ENG) %>%    # coerce to data frame to drop grouped and tbl attributes
  dplyr::summarise(count = n_distinct(REPORT_ID)) %>% as.data.frame() %>%
  filter(!is.na(HLT_NAME_ENG))

input_df <- as.PhViD_HCSC(count_df)
input_df_hlt <- as.PhViD_HCSC(count_df_hlt)



PT_count <- cbind(input_df$L, input_df$data) %>%
  rename(drug_code = ing, event_effect = PT_NAME_ENG,
         count = n11, drug_margin = n1., event_margin = n.1)

HLT_count <- cbind(input_df_hlt$L, input_df_hlt$data) %>%
    rename(drug_code = ing, event_effect = HLT_NAME_ENG,
           count = n11, drug_margin = n1., event_margin = n.1)

DATA <- input_df_hlt$data
N <- input_df_hlt$N
L <- input_df_hlt$L
n11 <- DATA[,1]
n1. <- DATA[,2] # marginal drug counts
n.1 <- DATA[,3] # marginal AE counts
n10 <- n1. - n11
n01 <- n.1 - n11
n00 <- N - (n11+n10+n01)
expected_count <- n1. * n.1 / N

PRR <- (n11 / (n11 + n10)) / (n01 / (n01 + n00))
logPRR <- log(PRR)
var_logPRR <- 1/n11 - 1/(n11 + n10) + 1/n01 - 1/(n01 + n00)
LB95_logPRR <- qnorm(0.025,logPRR,sqrt(var_logPRR))
UB95_logPRR <- qnorm(0.975,logPRR,sqrt(var_logPRR))
LB95_PRR <- exp(LB95_logPRR)
UB95_PRR <- exp(UB95_logPRR)
PRR_result <- data.frame(drug_code = L[[1]], event_effect = L[[2]],
                         count = n11, expected_count,
                         PRR, LB95_PRR, UB95_PRR,
                         logPRR, LB95_logPRR, UB95_logPRR,
                         var_logPRR,
                         stringsAsFactors = FALSE)
rm("PRR", "logPRR", "var_logPRR", "LB95_logPRR", "UB95_logPRR",
   "LB95_PRR", "UB95_PRR")

ROR <- n11 * n00 /(n10 * n01)
logROR <- log(ROR)
var_logROR <- 1/n11 + 1/n10 + 1/n01 + 1/n00
LB95_logROR <- qnorm(0.025,logROR,sqrt(var_logROR))
UB95_logROR <- qnorm(0.975,logROR,sqrt(var_logROR))
LB95_ROR <- exp(LB95_logROR)
UB95_ROR <- exp(UB95_logROR)
ROR_result <- data.frame(drug_code = L[[1]], event_effect = L[[2]],
                         count = n11, expected_count,
                         ROR, LB95_ROR, UB95_ROR,
                         logROR, LB95_logROR, UB95_logROR,
                         var_logROR,
                         stringsAsFactors = FALSE)
rm("ROR", "logROR", "var_logROR", "LB95_logROR", "UB95_logROR",
   "LB95_ROR", "UB95_ROR")

RRR <- (n11*N) / (n1.*n.1)
logRRR <- log(RRR)
var_logRRR <- 1/n11 - 1/n1. + 1/n.1 - 1/N
LB95_logRRR <- qnorm(0.025,logRRR,sqrt(var_logRRR))
UB95_logRRR <- qnorm(0.975,logRRR,sqrt(var_logRRR))
LB95_RRR <- exp(LB95_logRRR)
UB95_RRR <- exp(UB95_logRRR)
RRR_result <- data.frame(drug_code = L[[1]], event_effect = L[[2]],
                         count = n11, expected_count,
                         RRR, LB95_RRR, UB95_RRR,
                         logRRR, LB95_logRRR, UB95_logRRR,
                         var_logRRR,
                         stringsAsFactors = FALSE)
rm("RRR", "logRRR", "var_logRRR", "LB95_logRRR", "UB95_logRRR",
   "LB95_RRR", "UB95_RRR")
rm(n.1,n00,n01,n1.,n10,n11,L,DATA)


RFET_result <- RFET_HCSC(input_df_hlt)

BCPNN_result <- BCPNN_HCSC(input_df_hlt, MC = TRUE, NB.MC = 10000)
BCPNN_result %<>% dplyr::select(drug_code = `drug code`, event_effect = `event effect`,
                                count, expected_count = `expected count`, median_IC,
                                LB95_IC = `Q_0.025(log(IC))`, UB95_IC = `Q_0.975(log(IC))`,
                                postH0, `n11/E (RRR)` = `n11/E`)

GPS_result1 <- GPS(input_df_hlt, RANKSTAT = 2)
GPS_quantile <- as.data.frame(GPS_result1$ALLSIGNALS) %>% arrange(drug, event)
GPS_result2 <- GPS(input_df_hlt, RANKSTAT = 3)
GPS_postE <- as.data.frame(GPS_result2$ALLSIGNALS) %>% arrange(drug, event)
GPS_final <- data.frame(
  drug_code = GPS_quantile$drug,
  event_effect = GPS_quantile$event,
  count = GPS_quantile$count,
  expected_count = GPS_quantile$`expected count`,
  `n11/E` = GPS_quantile$`n11/E`,
  postH0 = GPS_quantile$postH0,
  postE_lambda = GPS_postE$`post E(Lambda)`,
  Q0.05_lambda = GPS_quantile$`Q_0.05(lambda)`
)

rm(GPS_result1,GPS_quantile,GPS_result2,GPS_postE)

soc_all <- tbl(hcopen,'meddra')%>%
  dplyr::select(PT_Term,HLT_Term,SOC_Term)
cv_soc <- soc_all %>%
  dplyr::select(HLT_Term,SOC_Term) %>%
  rename(event_effect = HLT_Term) %>%
  as.data.frame()
cv_soc <- cv_soc[!duplicated(cv_soc),]

cat_pt <- cv_soc %>%
  group_by(event_effect) %>%
  summarize(SOC_Term = paste(unique(SOC_Term), collapse = ", "))


cv_reportid <- cv_drug_rxn_new_2006 %>%
  dplyr::select(HLT_NAME_ENG, quarter) %>%
  rename(event_effect = HLT_NAME_ENG)%>%
  filter(event_effect!='')%>%
  distinct()

cat_id_pt <- cv_reportid %>%
  group_by(event_effect) %>%
  summarize(quarter = paste(unique(quarter), collapse = ", "))


mastertable_hlt_new<-BCPNN_result[,1:7] %>%
  inner_join(HLT_count[,c(1,2,4,5)], by = c("drug_code", "event_effect")) %>%
  inner_join(GPS_final[,c(1,2,6,7,8)], by = c("drug_code", "event_effect")) %>%
  inner_join(RFET_result[,c(1,2,5,6)], by = c("drug_code", "event_effect")) %>%
  inner_join(PRR_result[,c(1,2,5,6,7)], by = c("drug_code", "event_effect")) %>%
  inner_join(ROR_result[,c(1,2,5,6,7)], by = c("drug_code", "event_effect")) %>%
  inner_join(RRR_result[,c(1,2,5,6,7)], by = c("drug_code", "event_effect")) %>%
  inner_join(cat_pt,by='event_effect')%>%
  inner_join(cat_id_pt,by='event_effect')%>%
  dplyr::select(1:4, 24:25, 8:9, 5:23)

# testconnect <- dbConnect(drv = "PostgreSQL", host = "shiny.hc.local", user = "hcwriter", dbname = "hcopen", password = "canada2")
# dbWriteTable(testconnect, "cv_drug_rxn_meddra_new",cv_drug_rxn_meddra_new, row.names = FALSE)
# dbWriteTable(testconnect,'master_table_pt_170831',mastertable_pt_new,row.names=FALSE)
# dbWriteTable(testconnect,'master_table_hlt_170831',mastertable_hlt_new,row.names=FALSE)
# dbDisconnect(testconnect)