##' Integrate KEGG database into analysis pipeline in R
##'
# InputFiles is a list of .keg file names for parsing
InputFiles<-read.table(file="Parse_Input_Files.txt",header=FALSE,sep = ".")
for (i in 12:nrow(InputFiles)) {
  InputArg<-paste("python KEGG_parser.py", paste(InputFiles[i,1],InputFiles[i,2],sep = "."), "csv")
  system(InputArg)
  OutputArg<-paste(InputFiles[i,1],"xls",sep=".")
  nam<-paste("CompileHfile_",i,sep="")
  
  # Read the parsed files into R in batches
  assign(nam, read.xls(OutputArg,sheet=1,header=TRUE))
}
