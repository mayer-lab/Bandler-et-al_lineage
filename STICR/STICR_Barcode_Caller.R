library(data.table)
library(parallel)
library(Matrix)
library(ggplot2)
library(rlist)

######Functions######

Barcode_Collapse_Function<-function(x){
  VBCs<-as.character(working.cell.split[[x]]$barcode)
  return(VBCs)
}

Barcode_Calling_Function<-function(x){
  working.set<-working.cell.split[[x]]
  barcode.set<-as.character(working.set$barcode)
  hits<-list()
  for (i in barcode.set){
    hits<-list.append(hits,grep(i,Full_VBC.Single.ls, fixed = TRUE))
  }
  result<-ifelse(length(unique(unlist(hits)))>1,0, unique(unlist(hits)) )
  return(result)
} #Multiplet result will return a value of -1, Single Cell Clones will return a value of -2

Multi_VBC_count_function<-function(x){
  working.set<-working.cell.split[[x]]
  barcode.set<-as.character(working.set$barcode)
  VBC_count<-length(unlist(barcode.set))
  return(VBC_count)
} #checking to see what % of cells belonging to multi-VBC clones contain all VBCs

UMI_Count_Matching_function<-function(x){
  working.subset<-working.lib[working.lib$CBC == x,]
  UMI_Result<-ifelse(Barcodes.aligned.df$VBC[Barcodes.aligned.df$CBC==x]==0,0,working.subset$UMI_Count[working.subset$CBC==x])
  return(UMI_Result)
}

######Import and Format Data for Analysis######

#Import VBC calls
Full_BC.df<-readRDS("/STICR_Barcode_Extractor_outputFolder/Final_Barcodes.tsv") 

#Import "Final_Barcodes.tsv" file generated from the STICR Barcode Extractor Script as a data.frame named "Full_BC.df"
#Add an additional column called "Library" that describes the 10X transcriptomic library this STICR Barcode data was generated from. Value should be a character string.
#Dataframes from multiple libraries can be concatenated together into data.frame named "Full_BC.df".
#Final data.frame should contain 4 columns with names: "CBC","barcode","UMI_Count","Library"


#Create data.frame of "valid" CBCs- i.e. those in scRNA-seq transcriptomic set passing QC cutoffs
Clone_Seq_BC.df<-#valid CBC set
#Clone_Seq_BC.df should have 2 columns titled: "Library" and "CBC"
#Library values should be character strings and match values in Full_BC.df$Library
#CBC values should be character strings and match formating used in Full_BC.df$CBC. Of note, depending on how datasets were processed/merged in Seurat, CBC values might have a "-N" appended to the end that you will want to match.

#Subset VBC dataframe by CBCs detected in Transcriptomic set 
Full_BC.df<-Full_BC.df[Full_BC.df$CBC %in% Clone_Seq_BC.df$CBC,]

#Retain CBC/STICR combinations with at least 5 UMI count
Full_BC.df<-Full_BC.df[Full_BC.df$UMI_Count>4,] 


######Set variables######
numCores<-20 #set number of cores for multithreading,
y<-5 #Ratio of Most Prevelent STICR UMI counts vs. Second-most Prevelant STICR UMI counts used for calling "Dominant Barcode" when multiple STICR barcodes are recovered from a single CBC. Default is 5. Potential causes of multiple barcodes per CBC include transcript "leak", chimeric PCR, and superinfection.

outpath<-"/out/path" #change to output path
UMI.cutoff<-9 #Minimum UMI count threshold for a barcode to be assigned to a CBC. "Dominant" barcodes that do not exceed this threshold are not retained for further analysis.
library.set<-unique(Full_BC.df$Library) #create character vector of Library names





for (i in library.set){
  working.lib<-data.frame(Full_BC.df[Full_BC.df$Library == i,]) #Subset VBC Count Matrix by Library
  VBC_Barcode_Index<-data.frame(VBC=sort(unique(working.lib[,2])),Index=1:length(unique(working.lib[,2]))) #Create a numerical VBC index
  VBC_Single_Cell.vec<-names(table(working.lib$barcode)[table(working.lib$barcode)==1]) #create character vector of VBCs that appear in a single cell
  VBC_Barcode_Index$Single_Cell_Clone<-ifelse(VBC_Barcode_Index$VBC %in% VBC_Single_Cell.vec,1,0) #Add a metadata column to numerical VBC index indicating whether a VBC appears in just one cell (value=1) or not (value=0)
  VBC_Barcode_Index.comb<-data.frame(matrix(unlist(combn(VBC_Barcode_Index$Index[VBC_Barcode_Index$Single_Cell_Clone==0],2,simplify=FALSE)),ncol=2, byrow=TRUE)) #make all unique pairwise between STICR barcodes found in multiple cells
  MOI_1_VBC<-as.character(VBC_Barcode_Index$VBC[VBC_Barcode_Index$Single_Cell_Clone==0]) #Make vector of STICR barcodes found in more than one cell
  working.cell.split<-split(working.lib,working.lib$CBC) #split data.frame by CBC
  CBC_VBC_Calls<-lapply(names(working.cell.split),Barcode_Collapse_Function)
  names(CBC_VBC_Calls)<-names(working.cell.split)
  Full_VBC.ls<-MOI_1_VBC #Copy STICR barcode list to vector of different name for addition of other elements (i.e. STICR barcodes appearing in only one cell) later in the script. 
  Full_VBC.Single.ls<-c(Full_VBC.ls,VBC_Single_Cell.vec)
  Barcodes.aligned<-mclapply(1:length(working.cell.split), Barcode_Calling_Function, mc.cores = numCores) #For each CBC, returns all associated STICR barcodes. If only one STICR barcode is present, returns the Index value of that STICR barcode as defined in VBC_Barcode_Index. If multiple STICR barcodes found, returns a value of 0.
  Barcodes.aligned.df<-data.frame(CBC=names(working.cell.split),VBC=unlist(Barcodes.aligned)) #Turn matched STICR barcodes/CBCs into a dataframe
  Barcodes.aligned.df$VBC_Count<-as.numeric(lapply(Barcodes.aligned.df$CBC,Multi_VBC_count_function)) # Returns count of unique STICR barcodes associated with each CBC.
  Barcodes.aligned.df$Tier<-ifelse(Barcodes.aligned.df$VBC==0,0,1) #Creates new column called "Tier" by assigned CBCs that have only a single STICR barcode a Tier of 1, while those with multiple STICR barcodes get assigned a Tier of 0 and will be analyzed further below. Tier 0 STICR barcode/CBC pairings might be updated further depending on whether one of the associated STICR barcodes is determined to be "Dominant".
  #make final STICR barcode key
  Final.VBC.Index.df<-data.frame(VBC_Final=1:length(Full_VBC.Single.ls),Type=c(rep("MOI_1",length(MOI_1_VBC)),rep("SingleCellClone",length(VBC_Single_Cell.vec))))
  Final.VBC.Index.df$Type<-as.character(Final.VBC.Index.df$Type)
  #Attempt to call high confidence "dominant" barcodes from cells with multiple VBCs
  Multiplet.df<-Barcodes.aligned.df[Barcodes.aligned.df$VBC == 0,] #Subset Barcodes.aligned.df to create a data.frame called Multiplet.df containing CBCs with multiple STICR barcodes.  
  if (dim(Multiplet.df)[1]>0){
    Dom.ls<-list() #Creates an empty list to store STICR Barcode Indices/UMI counts for each CBC
    for (j in unique(Multiplet.df$CBC)){
      working.set<-working.cell.split[[j]]
      barcode.set<-as.character(working.set$barcode)
      hits<-c()
      for (x in barcode.set){
        hits[x]<-grep(x,Full_VBC.Single.ls, fixed = TRUE)
      }
      combined.matches<-data.frame(VBC=hits)
      combined.matches$UMI_Count<-working.set$UMI_Count[working.set$barcode %in% rownames(combined.matches)]
      rownames(combined.matches)<-c()
      VBC.unique<-data.frame(VBC=unique(combined.matches$VBC))
      VBC.unique$UMI_Count<- lapply(VBC.unique$VBC, function(x) max(combined.matches$UMI_Count[combined.matches$VBC == x]))
      Dom.ls[[j]]<-VBC.unique
    } #Fills Dom.ls with STICR Barcode Indices/UMI counts for each CBC
    Tier2_exploration<-data.frame(CBC=names(Dom.ls)) #make data.frame to hold results of STICR barcode "Dominance" analysis.
    Tier2_exploration$CBC<-as.character(Tier2_exploration$CBC)
    Tier2.results.ls<-list()
    for (k in Tier2_exploration$CBC){
      current.set<-Dom.ls[[k]]
      current.set$UMI_Count<-as.numeric(current.set$UMI_Count)
      Max.UMI<-max(current.set$UMI_Count)
      Sum.UMI<-sum(current.set$UMI_Count)
      Second.UMI<-as.numeric(current.set$UMI_Count)[order(as.numeric(current.set$UMI_Count))][length(current.set$UMI_Count)-1]
      Dom.bc<-current.set$VBC[current.set$UMI_Count==Max.UMI]
      results.merged<-paste(Max.UMI,Second.UMI,Sum.UMI,Dom.bc,sep="-")
      Tier2.results.ls[k]<-ifelse(length(results.merged)==1,results.merged, "0-0-0-0")
    } # For each cell, returns the UMI count from the most prevelant STICR barcode (Max.UMI), the UMI count fromt he second most prevelant STICR barcode (Second.UMI), Total UMI count of all STICR barcodes associated with cell (Sum.UMI), and the STICR barcode Index of the most prevelant barcode (Dom.bc), all separated by "-".  Will output "0-0-0-0" in the event 2 barcodes occur with the exact same UMI
    Tier2.results.df <- data.frame(do.call(rbind, str_split(Tier2.results.ls,"-"))) #removes "-" delimiter
    Tier2.results.df<-cbind(Tier2_exploration,Tier2.results.df) #adds matching CBC for Tier2.results.df in a column called "CBC".
    names(Tier2.results.df)[2:5]<-c("Max","Second","Total","DOM_BC") #add column names
    Tier2.results.df$Max<-as.numeric(as.character(Tier2.results.df$Max)) #formating
    Tier2.results.df$Second<-as.numeric(as.character(Tier2.results.df$Second)) #formating
    Tier2.results.df$Total<-as.numeric(as.character(Tier2.results.df$Total)) #formating
    Tier2.results.df$Prop_Dom<-Tier2.results.df$Max/Tier2.results.df$Total #calculates the proportion of most frequent barcode counts vs. all barcode counts from cell.
    Tier2.results.df$Prop_Second<-Tier2.results.df$Max/Tier2.results.df$Second #calculates the proportion of most frequent barcode count vs. second-most frequent barcode count.
    UMI.val<-ifelse(dim(Tier2.results.df[Tier2.results.df$Max==0,])[1]>0,"yes","no") #spot check for equal UMI CBC/STICR barcode combinations
    Equal.UMI.set.df<-Tier2.results.df[Tier2.results.df$Max==0,] #create subset Tier2 results with equal UMI
    Tier2.results.df<-Tier2.results.df[!Tier2.results.df$CBC %in% Equal.UMI.set.df$CBC,] #remove subset Tier2 results with equal UMI
    Tier2.results.df$DOM_BC<-as.numeric(as.character(Tier2.results.df$DOM_BC)) #Convert DOM_BC from factor into interger for VBC type matching
    Tier2.Viral.Subtype.vec<-c() #create empty vector in which to add VBC type in following line
    for (l in 1:dim(Tier2.results.df)[1]){
      VBC.key<-Tier2.results.df[l,5]
      Output.type<-Final.VBC.Index.df[Final.VBC.Index.df$VBC_Final==VBC.key,2]
      Tier2.Viral.Subtype.vec[l]<-Output.type
    } #add VBC type (ie. MOI_1, Superinfection, etc.) for dominant VBC within in each cell
    Tier2.results.df$VBC.type<-Tier2.Viral.Subtype.vec #add VBC type to Tier2 data.frame
    #make ggplot of Prop_Second vs Prop_Dom
    VBC_Jaccard.dom<-ggplot(Tier2.results.df, aes(x=Prop_Second, y=Prop_Dom )) + geom_point() + theme_bw() + geom_vline(xintercept = 5, linetype="dotted", color = "red", size=1.5) +xlim(0,50) +ylim(0,1) #make a dot plot depicting both the Dominant barcode vs. all STICR barcode ratio and the Dominant vs. Second STICR barcode ratio for each cell.
    title<-paste0("Multiple_VBC_plot_",i,".pdf") #Create a title for the histogram complete with file extension (.pdf)
    pdf(title, width = 4, height = 4) 
    print(VBC_Jaccard.dom + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                               panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + ggtitle(title) + xlab("Jaccard Similarity") + xlab("Log2(Count)"))
    dev.off()
    #Call Tier 2 cells (cells with multiple STICR barcodes, but one of which is "Dominant" meaning that it is equal to or greater than the previously defined ratio of most common STICR barcode count vs. second-most common STICR barcode count)
    Tier2.results.df$Tier2<-ifelse(Tier2.results.df$Prop_Second >=y,2,0) #Standard cutoff for Tier2 inclusion (First/Second VBC UMI ratio >5)
    Tier2.results.df$DOM_BC<-as.numeric(Tier2.results.df$DOM_BC)
    #Add back in CBCs that have multiple STICR barcodes with equal UMI counts
    if (UMI.val=="yes"){
      Equal.UMI.set.df$VBC.type<-"Same_UMI_Count"
      Equal.UMI.set.df$Tier2<-0 #Modify the set of potential Tier2 CBCs with equal counts of first/second VBCs in order to merge with rest of Tier2 DF
      Tier2.results.df<-rbind(Equal.UMI.set.df,Tier2.results.df)
    }
    #Add Multiplet entry to Final.VBC.Index.df
    Final.VBC.Index.df<-rbind(c(0,"Multiplet"),Final.VBC.Index.df)
    Final.VBC.Index.df$VBC_Final<-as.numeric(Final.VBC.Index.df$VBC_Final)
    #Update Barcode DF with CBCs resolved as being Tier2
    for (m in Barcodes.aligned.df$CBC[Barcodes.aligned.df$Tier==0]){
      position<-as.numeric(rownames(Barcodes.aligned.df[Barcodes.aligned.df$CBC== m, ]))
      Barcodes.aligned.df[position,2]<-ifelse(Tier2.results.df$Tier2[Tier2.results.df$CBC == m] ==2, Tier2.results.df$DOM_BC[Tier2.results.df$CBC == m], 0)
      Barcodes.aligned.df[position,4]<-ifelse(Tier2.results.df$Tier2[Tier2.results.df$CBC == m]==2, 2, 0)
    }
  }
  #Add VBC type to Final.VBC.Index.df
  Viral.Subtype.vec<-c()
  for (n in 1:dim(Barcodes.aligned.df)[1]){
    VBC.key<-Barcodes.aligned.df[n,2]
    Output.type<-Final.VBC.Index.df[Final.VBC.Index.df$VBC_Final==VBC.key,2]
    Viral.Subtype.vec[n]<-Output.type
  }
  Barcodes.aligned.df$VBC.type<-Viral.Subtype.vec
  Barcodes.aligned.df$CBC<-as.character(Barcodes.aligned.df$CBC)
  #add a multiple entry in index.df if not present
  if (!"Multiplet" %in% Final.VBC.Index.df$Type){
    Final.VBC.Index.df<-rbind(Final.VBC.Index.df,c(0,"Multiplet"))
  }
  #Add VBC sequence to Index
  Final.VBC.Index.df$Sequence<-c("none",Full_VBC.Single.ls)
  #Add VBC sequence to aligned.df
  Barcodes.aligned.df$Sequence<-Final.VBC.Index.df$Sequence[match(Barcodes.aligned.df$VBC,Final.VBC.Index.df$VBC_Final)] 
  #Run UMI Count Function
  Barcodes.aligned.df$UMI<-unlist(lapply(Barcodes.aligned.df$CBC,UMI_Count_Matching_function))
  #Apply UMI Cutoff
  Barcodes.aligned.df$VBC.type<-ifelse(Barcodes.aligned.df$UMI<UMI.cutoff,"LOW_UMI",Barcodes.aligned.df$VBC.type)
  Barcodes.aligned.df$VBC<-ifelse(Barcodes.aligned.df$UMI<UMI.cutoff,0,Barcodes.aligned.df$VBC)
  Barcodes.aligned.df$Tier<-ifelse(Barcodes.aligned.df$UMI<UMI.cutoff,0,Barcodes.aligned.df$Tier)
  #Save Files
  save.path<-paste(outpath,i,sep="/")
  Barcodes.aligned.df.path<-paste0(save.path,"_Barcodes.aligned.rds")
  Full_VBC.Single.ls.path<-paste0(save.path,"_Full_VBC.Single.ls.rds")
  Final.VBC.Index.df.path<-paste0(save.path,"_Final.VBC.Index.df.rds")
  saveRDS(Barcodes.aligned.df,Barcodes.aligned.df.path)
  saveRDS(Full_VBC.Single.ls,Full_VBC.Single.ls.path)
  saveRDS(Final.VBC.Index.df,Final.VBC.Index.df.path)
  rm(Barcodes.aligned.df)
  rm(Full_VBC.Single.ls)
  rm(Final.VBC.Index.df)
  rm(Tier2.results.df)
  rm(VBC_Jaccard.dom)
  rm(Equal.UMI.set.df)
  rm(Multiplet.df)
  rm(UMI.val)
}


#Make Final CBC/VBC Set
Final.File.set<-list.files("/out/path/")
Final.File.set.aligned<-grep("aligned",Final.File.set, value = T)
Final.File.set.aligned<-paste0("/out/path/",Final.File.set.aligned)
Barcode_Call_Combination_function<-function(x){
  working.set<-readRDS(x) #import RDS file
  working.table<-data.frame(table(working.set$VBC)) #make a data.frame tabulating the number of cells assigned each barcode (indexed)
  working.table$Var1<-as.numeric(as.character(working.table$Var1)) #Convert the barcode count from factor to numeric
  Clone.set<-working.table$Var1[working.table$Freq>1]
  Clone.set<-as.numeric(as.character(Clone.set))
  Clone.set<-Clone.set[!Clone.set==0] #remove VBC 0 which refers to cells in which a dominant VBC cannot be determined
  working.table$Call<-ifelse(working.table$Var1==0,"Multiple",ifelse(working.table$Var1 %in% Clone.set,"Clone","Single.Call"))
  working.set$Final.VBC.call<-working.table$Call[match(working.set$VBC,working.table$Var1)]
  return(working.set)
}
Full_VBC_CBC.ls<-lapply(Final.File.set.aligned,Barcode_Call_Combination_function)
Full_VBC_CBC.df<-do.call(rbind.data.frame, Full_VBC_CBC.ls)
Full_VBC_CBC.df$Library<-sub("\\_.*", "", Full_VBC_CBC.df$CBC)
Full_VBC_CBC.df<-Full_VBC_CBC.df[,c(1,3,4,6,7,8,9)]
saveRDS(Full_VBC_CBC.df,"Results.rds")





