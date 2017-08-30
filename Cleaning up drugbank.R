library(tidyr)


#clean up drugbank for mapping

drugbank<-fread('drugbank.csv',header=T,sep=',',stringsAsFactors = F)
drugbank<-drugbank%>%select(atc_code,name,synonyms)%>%mutate(atc=str_extract_all(atc_code,"\\([[:alnum:]]{7}\\)"))
drugbank[drugbank==""]<-NA


drugbank$atc<-sapply(drugbank$atc,function(x)paste(x,collapse='|'))


drugbank<-drugbank%>%
  mutate(atc=gsub("\\(|\\)","",atc))%>%
  mutate(synonyms=strsplit(synonyms,"\n"))%>%
  unnest(synonyms)%>%
  mutate(name=tolower(name),synonyms=tolower(synonyms))%>%
  mutate(name=gsub("[[:punct:]]"," ",name),synonyms=gsub("[[:punct:]]"," ",synonyms))%>%
  mutate(name=gsub("\\s+"," ",name),synonyms=gsub("\\s+"," ",synonyms))%>%
  mutate(name=trimws(name),synonyms=trimws(synonyms))%>%
  rename(active_ingredient_name=name)

drugbank$atc[drugbank$atc=="NA"]<-NA  

#Load data from updated WHO ATC code list and clean the table:                               
atc<-fread('updated_atc.csv',header=T,sep=',',stringsAsFactors = F)


#concatenate ing with multiple ATC into one row:
atc_final<-atc%>%
  mutate(atc_desc5=gsub("combinations.*|in combination.*|in combination with.*","",atc_desc5))%>%
  mutate(atc_desc5=gsub(" and .*","",atc_desc5))%>%#for drug combinations, concatenate ATC code
  mutate(atc_desc5=gsub("[[:punct:]]"," ",atc_desc5))%>%  #removing punctuations, trimming extra space
  mutate(atc_desc5=gsub("\\s+"," ",atc_desc5))%>%
  mutate(atc_desc5=tolower(str_trim(atc_desc5)))%>%
  group_by(atc_desc5)%>%
  mutate(atc_code5=paste(unique(atc_code5),collapse="|"))%>%
  `[`(c(1,2))%>%
  filter(atc_desc5!="")%>%
  distinct()
 



#combine atc_final with drugbank for the matching database:
atc_code<-rep(NA,nrow(atc_final))
synonyms<-rep(NA,nrow(atc_final))
atc_final<-data.frame(atc_code,atc_final,synonyms)%>%
  select(atc_code,active_ingredient_name=atc_desc5,atc=atc_code5,synonyms)
  

#Add synonyms for "aluminium":
atc_final[grep("aluminium",atc_final$active_ingredient_name),4]<-atc_final[grep("aluminium",atc_final$active_ingredient_name),2]
atc_final$synonyms<-gsub("aluminium","aluminum",atc_final$synonyms)


#Identify drugs from WHO ATC updated list that does not exist in DrugBank:
drugbank_unique<-distinct(drugbank,active_ingredient_name)
drugbank_unique_s<-distinct(drugbank,synonyms)
atc_join<-anti_join(atc_final,drugbank_unique)%>%anti_join(drugbank_unique_s,by=c("active_ingredient_name"="synonyms"))

#Supplement the drugbank with extra input for Mapping:
drugbank<-bind_rows(drugbank,atc_join)%>%distinct(active_ingredient_name,atc,synonyms,.keep_all=T)%>%filter(!is.na(atc))

#input ATC code for fosaprepitant:
drugbank[grep("fosaprepitant",drugbank$active_ingredient_name),3]<-"A04AD12"

rm(atc_code,synonyms,drugbank_unique,drugbank_unique_s,atc_join,atc_final,atc)