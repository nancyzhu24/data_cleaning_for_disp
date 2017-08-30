library(dplyr)
library(stringr)
library(magrittr)
library(tidyr)
library(DBI)
library(pool)
library(data.table)
library(stringdist)

#this package was developed to do ngram clustering in R, the same way done in OpenRefine:
devtools::install_github("ChrisMuir/refinr")
library(refinr)
  


#download NHP data:
# source('~/nhp/NHP import.R')

#clean up drugbank and ATC code datatable:
source('Cleaning up drugbank.R')

# link to Adminer:
cv <- dbPool(drv=RPostgreSQL::PostgreSQL(),
             host     = "rest.hc.local",
             dbname   = "cv",
             user     = "dbuijs",
             password = "carlito27",
             options = "-c search_path=cv_20170331")


table_name_cv_drug_product_ingredients <- "cv_drug_product_ingredients"


#Extract drug_product_id and active_ingredient_name from the table:
cv_drug_product_ingredients <- tbl(cv, table_name_cv_drug_product_ingredients)%>%
                                as.data.frame()%>%`[`(c(2,5))
                
#matching NHP medicinal & non-medicinal ingredients:

#clean up nhp_med names: remove all punctuation, extra spaces:                                                   
# nhp_med<-nhptables[[2]]
nhp_med<-fread('nhp_med.csv',header=TRUE,sep=",",stringsAsFactors = F)
nhp_med[nhp_med==""]<-NA
nhp_med<-nhp_med[!is.na(nhp_med$proper_name),]  ## remove NAs in proper name

nhp_med%<>% mutate(proper_name=tolower(str_trim(proper_name)),common_name=tolower(str_trim(common_name)))%<>%
  mutate(proper_name=gsub("[[:punct:]]"," ",proper_name),common_name=gsub("[[:punct:]]"," ",common_name))%<>%
  mutate(proper_name=gsub("\\s+"," ",proper_name),common_name=gsub("\\s+"," ",common_name))%>%
  distinct(proper_name,common_name)

#create an index number for unique common name
nhp_med$index<-group_indices(nhp_med,proper_name)


#Clean up cv_drug_product_ingredient table:
#remove ws, brackets at the end in active_ingredient_name
cv_drug_product_ingredients<-cv_drug_product_ingredients%>%
  mutate(ing=sub("\\(.*\\)$","",active_ingredient_name)) %>%
                                mutate(ing=gsub("[[:punct:]]"," ",ing))%>%
                                mutate(ing=gsub("\\s+"," ",ing))%>%
                                 mutate(ing=trimws(ing))

#create a template with both drug_product_id and ing name for later table joining
cv_drug_product_ingredients_template<-cv_drug_product_ingredients%>%
  filter(active_ingredient_name!="",active_ingredient_name!=0)%>%
  distinct()

#extract list of vaccines into a separate list:
cv_vaccine<-cv_drug_product_ingredients[grepl("vaccine|vaccines",cv_drug_product_ingredients$active_ingredient_name),]%>%
            distinct()

#create a template with both drug_product_id and ing name for later table joining
cv_drug_product_ingredients<-anti_join(cv_drug_product_ingredients,cv_vaccine)%>%
                             filter(active_ingredient_name!="",active_ingredient_name!=0)%>%
                             distinct(ing)%>%
                             rename(active_ingredient_name=ing)


#using ngram find active ingredient names that are different only due to word transposition, white space difference: 
#Eg: 17b estradiol vs. estradiol17b will be in the same cluster
#cluster the result

refin<-cv_drug_product_ingredients$active_ingredient_name%>%
           key_collision_merge()%>%
           n_gram_merge(numgram=2,
                        edit_threshold=1,
                        edit_dist_weights=c(d = 0.33, i = 0.33, s = 1, t = 0.5),
                        bus_suffix=FALSE)

ngram_result<-data_frame(active_ingredient_name=cv_drug_product_ingredients$active_ingredient_name,cluster=refin)
              

#first "Easy" match: with both name and synonyms from drugbank and names from nhp
cv_ing_atc<-left_join(ngram_result,drugbank[,c(2,3)])%>%
            left_join(drugbank[,c(3,4)],by=c("active_ingredient_name"="synonyms"))%>%
            left_join(nhp_med[,c(1,3)],by=c("active_ingredient_name"="proper_name"))%>%
            left_join(nhp_med[,c(2,3)],by=c("active_ingredient_name"="common_name"))%>%
            mutate(atc=coalesce(atc.x,atc.y))%>%
            mutate(nhp=coalesce(index.x,index.y))%>%
            select(-atc.x,-atc.y,-index.x,-index.y)%>%
            distinct()

#Assign atc/nhp code for cluster: (assign the atc/nhp code to first non-NA value in a cluster group)

assign_cluster<-function(x){
  x%>%group_by(cluster)%>%
    mutate(cluster_n=sum(is.na(atc)))%>%
    mutate(atc=ifelse(cluster_n>0,atc[!is.na(atc)][1],atc))%>%
    mutate(nhp=ifelse(cluster_n>0,nhp[!is.na(nhp)][1],nhp))%>%
    ungroup()%>%
    select(-cluster_n)
}

cv_ing_atc<-assign_cluster(cv_ing_atc)
                         



cv_ing_atc_NA<-cv_ing_atc%>%filter(is.na(atc))

#create a new column with salt names,words like: herb, powder, dry, leaf, extract removed:
words<-c("dried","dry","powder","herb","leaf","extract","oil" ,"of","essential","seed","bark","flower","juice","root","stem",
         "bisulfate","sulfate","sulphate","diacetate","acetates","dihydrochloride","dihydrogencitrate","adipate",
         "hydrochloride","nos","dihcl","hcl","concentrate","disodium","dimeglumine","meglumine",
         "sodium$","hemitartrate","bitartrate","tartrate","mesylate","dimesilate","mesilate","besilate","borate",
         "salts","salt","not specified","succinate","tetrahydrate","syrup","chlorhydrate","camsilate","trihydrate",
         "tetrahydrate","bromohydrate","hexahydrate","heptahydrate","pentahydrate",
         "monohydrate","hemihydrate","dihydrate","dehydrated","otc","decahydrate","hydrated",
         "hydrate","orthophosphate","pyridoxalphosphate","pyrophosphate","biphosphate","phosphate",
         "hydrobromide","butylbromide","bromide","aminoxide","recombinant","fruit","prepared",
         "citrated","citrate","resinate","hemifumarate","hemifumarat","fumarate","fumarat","dipotassium","aspartate",
         "hemisuccinate","succinate","saccharate","stearate","acetate","nitrate","gluconate","maleate",
         "potassium","chloride","methylsulfate","tromethamine","propionate","furoate","protamine",
         "embonate","combo","amide","implant","hemodialysis","hemodialys","hydrogenated","hydrogen","fraction",
         "blinded","injection","refined","magnesium","calcium","water","\\d{1} water","spp","pamoate",
         "fluid","coated","liposome","sesquihydrate","fermentated","racemic","terephthalate","with.*","monohydrochloride",
         "anhydrous")

words<-sapply(words,function(x)paste0('\\b',x,'\\b'))  #for exact word match

cv_ing_atc_NA$ing<-str_replace_all(cv_ing_atc_NA$active_ingredient_name,paste(words,collapse='|'),"")
  
 
cv_ing_atc_NA%<>% mutate(ing=trimws(ing))%<>%
  select(-atc,-nhp)

#change all "vit" to "vitamin" in cv_ing_atc:
cv_ing_atc_NA%<>%mutate(ing=str_replace_all(ing,"vit\\s","vitamin "))

#match again:
cv_ing_atc_NA<-cv_ing_atc_NA%>%
  left_join(drugbank[,c(2,3)],by=c("ing"="active_ingredient_name"))%>%
  left_join(drugbank[,c(3,4)],by=c("ing"="synonyms"))%>%
  mutate(atc=coalesce(atc.x,atc.y))%>%
  left_join(nhp_med[,c(1,3)],by=c("ing"="proper_name"))%>%
  left_join(nhp_med[,c(2,3)],by=c("ing"="common_name"))%>%
  mutate(nhp=coalesce(index.x,index.y))%>%
  select(-atc.x,-atc.y,-index.x,-index.y)%>%
  distinct()


#update cv_ing_atc with new values
cv_ing_atc$ing<-"NA"  #create ing column in cv_ing_atc for matching purpose, rearrange columns
cv_ing_atc%<>%select(colnames(cv_ing_atc_NA))
cv_ing_atc[match(cv_ing_atc_NA$active_ingredient_name,cv_ing_atc$active_ingredient_name),]<-cv_ing_atc_NA

#Again, assign for the cluster:
cv_ing_atc<-assign_cluster(cv_ing_atc)

cv_ing_atc_NA<-cv_ing_atc%>%filter(is.na(atc)&is.na(nhp))

#matching unmatched result to drugbank using jaro-wrinkle distance: 
#to improve matching accuracy, need to expand the valid_list to foreign drugs and chemicals(PubChem).
#Currently, only drug list from WHO ATC and drugBank is used, not sufficient to guarantee accuracy

test_list<-cv_ing_atc_NA$ing
valid_list1<-c(unique(drugbank$active_ingredient_name),unique(nhp_med$proper_name),
               unique(drugbank$synonyms[!is.na(drugbank$synonyms)]),
               unique(nhp_med$common_name[!is.na(nhp_med$proper_name)]))

name_compare<-data.frame(test=test_list,active_ingredient_name=valid_list1[amatch(test_list,valid_list1,method="jw",p=0.1)],
                         stringsAsFactors = F)%>%
              filter(test!="")

#assuming 100% accurary from jaro wrinkler distance string matching: update cv_ing_atc table
name_compare<-left_join(name_compare,distinct(drugbank[,c(2,3)]))%>%
              left_join(distinct(drugbank[!is.na(drugbank$synonyms),c(4,3)]),by=c("active_ingredient_name"="synonyms"))%>%
              left_join(distinct(nhp_med[,c(1,3)]),by=c("active_ingredient_name"="proper_name"))%>%
              left_join(distinct(nhp_med[!is.na(nhp_med$common_name),c(2,3)]),by=c("active_ingredient_name"="common_name"))%>%
              mutate(atc=coalesce(atc.x,atc.y),nhp=coalesce(index.x,index.y))%>%
              select(-atc.x,-atc.y,-index.x,-index.y)%>%
              filter(!is.na(atc)|!is.na(nhp))%>%
              distinct()

cv_ing_atc[match(name_compare$test,cv_ing_atc$ing),'atc']<-name_compare$atc
cv_ing_atc[match(name_compare$test,cv_ing_atc$ing),'nhp']<-name_compare$nhp

#assign value within cluster:
cv_ing_atc<-cv_ing_atc%>%select(-ing)%>%assign_cluster()%>%distinct()
            


#final_unmatch<-cv_ing_atc%>%filter(is.na(atc),is.na(nhp))%>%filter(active_ingredient_name!="")
#match back to drug_product_id
cv_drug_product_ingredients_atc<-left_join(cv_drug_product_ingredients_template,cv_ing_atc,by=c('ing'='active_ingredient_name'))%>%
                                  select(-cluster,-ing)

#sensitivity analysis for disproportionality analysis in Canada Vigilance database:

#add column in cv_drug_product_ingredient: category
cv_drug_product_ingredients_atc$category<-"NA"
cv_drug_product_ingredients_atc[cv_drug_product_ingredients_atc$drug_product_id %in% cv_vaccine$drug_product_id,'category']<-'vaccine'
#ATC code starts with J07 are also vaccines
cv_drug_product_ingredients_atc[grepl("^J07",cv_drug_product_ingredients_atc$atc),'category']<-'vaccine'


cv_drug_product_ingredients_atc[!is.na(cv_drug_product_ingredients_atc$atc)&cv_drug_product_ingredients_atc$category=="NA",'category']<-'drug'
cv_drug_product_ingredients_atc[!is.na(cv_drug_product_ingredients_atc$nhp)&cv_drug_product_ingredients_atc$category=="NA",'category']<-'NHP'


#Mapping to drugname in cv_drug_product_ingredients and ing name from DrugBank
cv_drug<-tbl(cv,"cv_drug_product_ingredients")%>%select(drug_product_id,drugname)%>%as.data.table()%>%distinct()
cv_drug_product_ingredients_atc%<>%left_join(cv_drug)%<>%filter(category!="NA")


#separate lists into Drug,NHP and Vaccine and transform the table into final format:
drug<-cv_drug_product_ingredients_atc%>%
      filter(category=="drug")%>%
      left_join(drugbank[,c(2,3)],by='atc')%>%
      select(DRUGNAME=drugname,ing=active_ingredient_name.y,category)%>%
      distinct()

NHP<-cv_drug_product_ingredients_atc%>%
  filter(category=="NHP")%>%
  left_join(nhp_med[,c(1,3)],by=c('nhp'='index'))%>%
  select(DRUGNAME=drugname,ing=proper_name,category)%>%
  distinct()

vaccine<-cv_drug_product_ingredients_atc%>%
  filter(category=="vaccine")%>%
  select(DRUGNAME=drugname,ing=active_ingredient_name,category)%>%
  distinct()

final_list<-bind_rows(drug,NHP,vaccine)

#remove intermediate tables:
rm(name_compare,test_list,valid_list1,cv_drug_product_ingredients_template,cv_ing_atc_NA,refin,ngram_result,words,
   drug,NHP,vaccine,cv_drug_product_ingredients,cv_drug,cv_ing_atc,drugbank,nhp_med,table_name_cv_drug_product_ingredients,cv_vaccine)






#use partial string match (match strings in active_ingredient_name that contains strings in proper_name:
#try 3 methods:
#-grep for full partial match
#problem: duplicates ie: insulin zinc match with both zinc and insulin

#-jaro-winkler string distance:
#problem: accuracy not guranteed

#-levenshtein Distance, similar to jaro-winkler method:

#grep:
# ing_to_match<-sapply(drugbank$active_ingredient_name,function(x)paste0('\\b',x,'\\b')) #match exact word, not a substring
# idx2<-sapply(ing_to_match,grep,final_unmatch$active_ingredient_name)
# idx1<-sapply(seq_along(idx2),function(i)rep(i,length(idx2[[i]]))) #create matching index for binding the table
# 
# match2<-data.frame(test=final_unmatch$active_ingredient_name[unlist(idx2),drop=F],
#                    valid=drugbank$active_ingredient_name[unlist(idx1),drop=F],
#                    stringsAsFactors = F)%>% distinct()
# 
# #from match2 list, remove rows with common chemical names: Eg: calcium, potassium... for improved accuracy
# #list generated from valid names with the most count number: top20
# #match2%>%group_by(valid)%>%summarise(count=n())%>%arrange(desc(count))%>%head(n=20)
# names_to_remove<-c("calcium","potassium","phenol","ethanol","magnesium","salicylic acid","protamine","acetic acid","interferon",
#                     "ammonium chloride","aluminium","alum","ethanolamine","cina")
# 
# match2%<>%filter(!valid%in%names_to_remove)
# 
# #for active_ingredients with multiple matches, keep the one with more characters(avoid substring match)
# match2_1<-match2%>%group_by(test)%<>%filter(n()>1)
# match2_1%<>%group_by(test)%<>%filter(nchar(valid)==max(nchar(valid)))                           
# match2<-left_join(match2,match2_1,by="test")%>%mutate(valid=coalesce(valid.y,valid.x))%>%select(-valid.y,-valid.x)%>%distinct()
# 
# #Map ATC code to match2:
# match2<-left_join(match2,drugbank[,c(2,3)],by=c("valid"="active_ingredient_name"))%>%distinct()
# 
# cv_ing_atc[match(match2$test,cv_ing_atc$active_ingredient_name),'atc']<-match2$atc
# 
# #Mapping NHP
# nhp_to_match<-sapply(nhp_med$proper_name,function(x)paste0('\\b',x,'\\b'))
# idx4<-sapply(nhp_to_match,grep,final_unmatch$active_ingredient_name)
# idx3<-sapply(seq_along(idx4),function(i)rep(i,length(idx4[[i]]))) #create matching index for binding the table
# 
# match3<-data.frame(test=final_unmatch$active_ingredient_name[unlist(idx4),drop=F],
#                    valid=nhp_med$proper_name[unlist(idx3),drop=F],
#                    stringsAsFactors = F)%>% distinct()
# 
# names_to_remove_2<-c("calcium","potassium","sodium","pollen","magnesium","phosphate","citrus",
#                      "bovine","aluminum","aluminium","chloride","urea","zinc","pollen extract","citrate")
# 
# match3%<>%filter(!valid%in%names_to_remove_2)
# 
# #for active_ingredients with multiple matches, keep the one with more characters(avoid substring match)
# match3_1<-match3%>%group_by(test)%<>%filter(n()>1)
# match3_1%<>%group_by(test)%<>%filter(nchar(valid)==max(nchar(valid)))                           
# match3<-left_join(match3,match3_1,by="test")%>%mutate(valid=coalesce(valid.y,valid.x))%>%select(-valid.y,-valid.x)%>%distinct()
# 
# #Map ATC code to match3:
# match3<-left_join(match3,nhp_med[,c(1,3)],by=c("valid"="proper_name"))%>%distinct()
# 
# cv_ing_atc[match(match3$test,cv_ing_atc$active_ingredient_name),'nhp']<-match3$index
# 
# 
# 

# 
# #levenshtein distance:
# #final_unmatch$partials= as.character(sapply(final_unmatch$active_ingredient_name, agrep,drugbank$active_ingredient_name,max.distance =0.1,value=T))
# 
# #results comparable to jaro-winkler, but significantly slower.