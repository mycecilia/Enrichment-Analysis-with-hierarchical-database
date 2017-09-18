##' Parsed obo file with customized set of GO terms
##' Integrate the obo file into R analysis pipeline

GOList<-unique(Pv_GO[which(as.character(toupper(Pv_GO$Gene.Name))%in%as.character(Sig_input$name)),])
write.table(unique(GOList$GO.ID),file="GOList_IVDMD_Modules.txt",quote=FALSE,row.names = FALSE,col.names = FALSE)
# This Perl script parse the GO database OBO format
system("perl ./ONTO-PELforEnrichment.pl -f GO.obo -t GOList_IVDMD_Modules.txt",intern = TRUE)

# Root terms need to be excluded from the results
hfile<-read.table(file=">>hgo.txt",header=TRUE)
hfile<-hfile[-which(as.character(hfile$progeny)==as.character(hfile$query)),]
