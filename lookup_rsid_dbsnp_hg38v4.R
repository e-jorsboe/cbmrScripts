##############
## SCRIPT FOR ASSIGNING RSID to variants, using either DANMAC5 annotation (rsID) or the REGENIE pipeline annotations (avsnp150)
## or from DBSNP (using this file /data/preprocessed/genetics/RSID_map/RSID_map.HG38_19_18.tsv)
##################

library(data.table)

args <- commandArgs(trailingOnly=TRUE)
inputname<-args[1]
idcolumn<-args[2]
anyRSID<-as.logical(as.integer(args[3]))
map<-fread2("/data/preprocessed/genetics/RSID_map/RSID_map.HG38_19_18.tsv")
input<-fread2(inputname)

## All the rows with duplicated IDs have duplicated RSIDs
map<-map[!duplicated(map$ID_Original),]
rownames(map)<-map$ID_Original

map<-map[,c("ID_Original","RSID_dbsnp","HG19_CHR","HG19_POS")]
map<-map[ map$ID_Original%in%input[,idcolumn],]

## this lookup is supposed to be faster according to:
## https://stackoverflow.com/questions/29399259/what-is-the-fastest-way-to-select-a-row-of-a-data-frame-with-a-specific-id-value
setkey(setDT(map), "ID_Original")
fun2<-function(x){ if(x%%10000==0){print(x)}; return(unlist(map[.(input[,idcolumn][x]),"RSID_dbsnp"])) }

fun19chr<-function(x){ if(x%%10000==0){print(x)}; return(unlist(map[.(input[,idcolumn][x]),"HG19_CHR"])) }
fun19pos<-function(x){ if(x%%10000==0){print(x)}; return(unlist(map[.(input[,idcolumn][x]),"HG19_POS"])) }

res<-sapply(1:nrow(input),fun2)
input$RSID_dbsnp<-unlist(res)
input$RSID_dbsnp[ input$RSID_dbsnp==""]<-NA

res19chr<-sapply(1:nrow(input),fun19chr)
input$HG19_CHR<-unlist(res19chr)
input$HG19_CHR[ input$HG19_CHR==""]<-NA

res19pos<-sapply(1:nrow(input),fun19pos)
input$HG19_POS<-unlist(res19pos)
input$HG19_POS[ input$HG19_POS==""]<-NA


suffix<-tail(unlist(strsplit(inputname,"\\.")),1)

if(anyRSID){
        input$any_RSID<-apply(input[ ,c("rsID","avsnp150","RSID_dbsnp")],1,function(x) na.omit(x)[1])
        inputname2<-sub(suffix,paste0("withHg19.noDups.RSID_dbsnp.anyRSID.",suffix),inputname)
} else{
        inputname2<-sub(suffix,paste0("withHg19.noDups.RSID_dbsnp.",suffix),inputname)
}

print(suffix)
print(inputname2)

print("This many rows with a new dbsnp RSID:")
print(sum(!is.na(input$RSID_dbsnp)))
if(anyRSID){
        print("This many rows with a new RSID from: rsID, avsnp150, dbsnp (any of them)")
        print(sum(!is.na(input$any_RSID)))
}
write.table(input,inputname2,col.names=T,row.names=F,qu=F)
