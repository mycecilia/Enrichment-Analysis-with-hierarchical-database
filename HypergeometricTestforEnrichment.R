###'
###'Use in combination with EnrichmentAnalysis.R
###'
###The enrichment function starts here
#Sig_input is in the Cytoscape PINBPA format
#Pv_Annot is the Pv_KO in the phytozome biomart format
#hfile is the Path_KO in the parsed annotation htxt format
HypergeometricTestforEnrichment<-function(Sig_input, Pv_Annot,hfile,Trait_GAR){
  Enrichment.hgeo<-c(NA,NA,NA,NA,NA,NA,NA,NA,NA)
  #Select one set of significant module (SM) as input genes
  for (j in 1: max(Sig_input$ModuleID)){
    Sig_module<-Sig_input[which(Sig_input$ModuleID == j),]
    #Read the phytozome gene annotation
    Pv_caterms<-Pv_Annot[,c("Gene.Name","Term.ID")]
    Pv_caterms<-Pv_caterms[which(rowSums(is.na(Pv_caterms)) != 1),] 
    Pv_caterms<-unique(Pv_caterms) #Empty rows and other annotations were eliminated
    #Total number of reference genes is the total tested GWAS genes filtered by Pv_caterms (Ko_number & Pv_gene database)
    Gene_total<-Pv_caterms[which(Pv_caterms$Gene.Name%in%Trait_GAR$Gene),]
    #  Gene_total<-Pv_caterms[which(Pv_caterms$Gene.Name%in%unique(annot$gene)),]
    Ntotal<-length(unique(Gene_total$Gene.Name))
    #
    Sig_total<-Sig_module[which(as.character(Sig_module$Gene.Name)%in%toupper(as.character(Pv_caterms$Gene.Name))),]
    for (i in 1:nrow(Sig_total)) {
      Sig_total$Term.ID[i]<-Pv_caterms[which(toupper(as.character(Pv_caterms$Gene.Name))==as.character(Sig_total$Gene.Name[i])),]$Term.ID
    }
    #Initial value of the result table
    
    #Class 1; DE; n1+
    n1_<-length(Sig_total$Gene.Name)
    #Class 2; Non-DE; n2+
    n2_<-Ntotal-n1_
    ##'Annotation database filtered by 
    ##'
    Cate_total<-hfile[which(as.character(hfile[,1]) %in% unique(Gene_total$Term.ID)),] #Originated from Annotation htxt database
    Cate_Sig<-hfile[which(as.character(hfile[,1])%in%unique(Sig_total$Term.ID)),]
    Cate_name <- NA
    for (i in 2:(ncol(Cate_Sig)/2)) {
      Cate_name<-unique(c(Cate_name,as.character(Cate_Sig[,2*i])))
    }
    Cate_name<-na.omit(Cate_name)
    
    #Get the name of the testing category
    for (m in 1:length(Cate_name))
      {
      cate_name<-Cate_name[m]
      #Get the level of heirarchy file for the Cate_name
      Cond=0
      for (i in 1:(ncol(Cate_Sig)/2)) {
        Cond <- nrow(Cate_Sig[which(as.character(Cate_Sig[,2*i])%in%as.character(cate_name)),])
        print(i)
        print(Cond)
        if (Cond >= 1) {FindHierarchy=i; break}
      }
      if (FindHierarchy <= (ncol(Cate_total)/2 -2)) 
        { 
        #Within a category
        n_1 <- nrow(Cate_total[which(Cate_total[,2*FindHierarchy]%in%cate_name),])
        #Not category
        n_2 <- Ntotal-n_1
        #Both significant and within a category
        n11 <- nrow(Sig_total[which(as.character(Sig_total$Term.ID)%in%as.character(Cate_Sig[which(Cate_Sig[,2*FindHierarchy]%in%cate_name),][,1])),])
        n11_name<-paste(as.vector(Sig_total[which(Sig_total$Term.ID%in%Cate_Sig[which(Cate_Sig[,2*FindHierarchy]%in%cate_name),][,1]),]$Gene.Name), collapse = "/")
        #Calculate the hypergeometric one-sided enrichment p-value
        Pvalue<-phyper(n11,n_1,n_2,n1_,lower.tail = FALSE) + 0.5*dhyper(n11,n_1,n_2,n1_)
        Enrichment.hgeo<-rbind(Enrichment.hgeo, 
                               c(j,
                                 FindHierarchy,
                                 cate_name,
                                 n11_name,
                                 n11,
                                 Ntotal,
                                 n1_,
                                 n_1,
                                 Pvalue)) #The values go into the result data table
      }
    }
  }
  Enrichment.hgeo<-as.data.frame(na.omit(Enrichment.hgeo))
  colnames(Enrichment.hgeo)<-c("ModuleID","Node","Node.Name","Gene.Name","#SignificantInNode","#TotalGene","#ofSignificance","#inNode","P")
  Enrichment.hgeo$P_BH<-p.adjust(as.numeric(as.character(Enrichment.hgeo$P)), method="BH")
  return(Enrichment.hgeo)
}
