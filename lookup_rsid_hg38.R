library(data.talbe)

args <- commandsArgs(trailingOnly=TRUE)
inputname<-args[1]
idcolumn<-args[2]
mapname<-args[3]
map<-fread(mapname,data.table=F)
input<-fread(inputname,data.table=F)

rownames(map)<-map$ID_Original

map<-map[,c("ID_Original","RSID_dbsnp")]
map<-map[ map$ID_Original%in%input$SNP,]

# this look us supposed to be faster according to 
# https://stackoverflow.com/questions/29399259/what-is-the-fastest-way-to-select-a-row-of-a-data-frame-with-a-specific-id-value
setkey(setDT(map), "ID_Original")
fun2<-function(x){ if(x%%10000==0){print(x)}; return(unlist(map[.(input$SNP[x]),"RSID_dbsnp"]))}

res<-sapply(1:nrow(input),fun2)

input$RSID_dbsnp<-unlist(res)
input$RSID_dbsnp[ input$RSID_dbsnp==""]<-NA

input$any_RSID<-apply(input[ ,c("rsID","avsnp150","RSID_dbsnp")],1,function(x) na.omit(x)[1])
suffix<-tail(unlist(strsplit(inputname,"\\.")),1)

inputname2<-sub(suffix,paste0("noDups.RSID_dbsnp.anyRSID.",suffix),inputname)

print(suffix)
print(inputname2)

print("This many rows with a new dbsnp RSID:")
print(sum(!is.na(input$RSID_dbsnp)))
print("This many rows with a new RSID from: rsID, avsnp150, dbsnp (any of them)")
print(sum(!is.na(input$any_RSID)))

write.table(input,inputname2,col.names=T,row.names=F,qu=F)
  
                      

                      
