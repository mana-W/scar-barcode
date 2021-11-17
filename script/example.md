#example
#data is from UMI_CB/CB_UMI
#fa is ref file
#cutsite is a file define each sgRNA start and end positon

```
library(Slin)
data<-read.table("UMI_CB/CB_UMI",sep="\t",header=T)
ref<-ReadFasta("V3.fasta")
cutsite<-read.table("V3.cutSites",col.names = c("indx","start","end"))
scarref<-ReadCutsite(cutsite)
scarref_all<-ReadCutsite(cutsite,reftype="ALL")
```

#array identify
```
scarinfo<-FindScar(data=data,scarfull=ref,scar=cutsite,indel.coverage="ALL",type="test",cln=8)
scarform<-INDELChangeForm(scarinfo,scarref = scarref_all,cln=4)
scarinfo$Scar<-scarform
```
#define scar pattern for each cell
```
cellsinfo<-
```

#pattern visualization 
