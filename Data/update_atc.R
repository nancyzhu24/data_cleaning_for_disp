#updating atc_codes on postgrel:
source('whoatc.R')

#for single update:
extractATC<-function(atc_prefix,code){
output_v<-vector()
for (i in seq_along(code)){
  output_v[i]<-paste0(atc_prefix,code[i])
}
  
output_list<-list()
  for (i in seq_along(output_v)){
    output_list[[i]]<-whoatc(output_v[i])
  }

output_df<-do.call(rbind.data.frame,output_list)%>%distinct()
}

# #extracting ATC index:
# V09XX<-c("AA","AB","AX","BA","CA","CX","DA","DB","DX","EA","EB","EX","FX","GA","GB","GX",
#          "HA","HB","HX","IA","IB","IX","XA","XX")
# V09_df<-extractATC("V09",V09XX)
# 
# V10XX<-c("AA","AX","BX","XA","XX")
# V10_df<-extractATC("V10",V10XX)
# 
# L01XE<-c("35","36","37","38","39","40")
# L01XE_df<-extractATC("L01XE",L01XE)
# 
# L01X<-c("X51","X52","X53","Y","C23","C24","C26")
# L01X<-extractATC("L01X",L01X)
# 
# L04A<-c("C13","C14","A36","A37")
# L04A<-extractATC("L04A",L04A)
# 
# L03AA16<-whoatc("L03AA16")

###############################################################
#For bulk update:
#generate pseudo ATC codes to search for ATC hierachy:
# level1<-c('A','B','C','D','G','H','J','L','M','N','P','R','S','V')
# level2<-paste0("0",c(1:9))
# level2<-c(level2,c(10:20))
# 
# level3<-c('A','B','C','D','E','F','G','H','I','J','X')
# 
# code<-list()
# for(i in seq_along(level1)){
#   code[[i]]<-paste0(level1[i],level2)
# }
# 
# 
# code<-unlist(code)
# 
# code_f<-list()
# for(i in seq_along(code)){
#   code_f[[i]]<-paste0(code[i],level3)
# }
# 
# code_f<-unlist(code_f)%>%sort()%>%as.data.frame()
# 
# code_f<-read.csv('code.csv',header=F,stringsAsFactors = F)
# 
# atccode<-list()
# for(i in seq_along(code_f[[1]])){
#   atccode[[i]]<-whoatc(code_f$V1[i])$code
# }
# 
# atccode2<-unlist(atccode)%>%as.data.frame()%>%distinct()
# atccode2%<>%filter(grepl(".{5}",atccode2$.))
# 
# updated_atc<-read.csv('updated_atc.csv',stringsAsFactors = F)
# #find missing atc_code4 in updated_atc
# 
# dif<-setdiff(atccode2$.,updated_atc$atc_code4)
# 
# remove<-c("A02AG","A02AH","A02AX","A03DC","A03EA","A03ED", "A07BB","A07CA","V07AA","V07AB","V07AC","V07AD",
#           "V07AN","V07AR","V07AS","V07AT","V07AV","V07AX","V07AY","V07AZ"	)
# dif<-dif[!dif%in%remove]
# 
# atc2017<-list()
# for(i in seq_along(dif)){
#   atc2017[[i]]<-whoatc(dif[i])
# }



#transform into table:
table<-rbind(V09_df,V10_df,L01XE_df,L01X,L04A,L03AA16)%>%distinct()
atc_code<-table[nchar(table$code)==7,]
atc_index<-anti_join(table,atc_code)
colnames(atc_code)<-c("atc_code5","atc_desc5")
atc_table<-atc_code%>%mutate(atc_code1=substr(atc_code$atc_code5,1,1))%>%
                     mutate(atc_code2=substr(atc_code$atc_code5,1,3))%>%
                     mutate(atc_code3=substr(atc_code$atc_code5,1,4))%>%
                      mutate(atc_code4=substr(atc_code$atc_code5,1,5))%>%
                      left_join(atc_index,by=c("atc_code1"="code"))%>%
                      rename(atc_desc1=desc)%>%
                      left_join(atc_index,by=c("atc_code2"="code"))%>%
                      rename(atc_desc2=desc)%>% 
                      left_join(atc_index,by=c("atc_code3"="code"))%>%
                       rename(atc_desc3=desc)%>%
                      left_join(atc_index,by=c("atc_code4"="code"))%>%
                      rename(atc_desc4=desc)

#pull table from hcref to update:                     
# hcref_pool <- DBI::dbConnect(drv = "PostgreSQL",
#                              host = "shiny.hc.local",
#                              dbname = "hcref",
#                              user = "hcreader",
#                              password = "canada1")
# 
# 
# atc<-hcref_pool%>%
#   tbl("atc_codes")%>%
#   as.data.frame()

atc<-bind_rows(atc,atc_table)%>%distinct()
write.csv(updated_atc,'updated_atc.csv',row.names = F)
testconnect <- dbConnect(drv = "PostgreSQL", host = "shiny.hc.local", user = "hcwriter", dbname = "hcref", password = "canada2")
dbWriteTable(testconnect,"atc_codes_20170718",atc, row.names = FALSE)
                      