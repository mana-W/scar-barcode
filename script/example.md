example<br />
data is from UMI_CB/CB_UMI<br />
fa is ref file<br />
cutsite is a file define each sgRNA start and end positon<br />
celltype.tsv is a file include cell barcode and its' annotations, header: Cell.BC Cell.type

```
library(Slin)
data<-read.table("UMI_CB/CB_UMI",sep="\t",header=T)
ref<-ReadFasta("V3.fasta")
cutsite<-read.table("V3.cutSites",col.names = c("indx","start","end"))
scarref<-ReadCutsite(cutsite)
scarref_all<-ReadCutsite(cutsite,reftype="ALL")
celltype<-read.table("celltype.tsv",header=T,stringsAsFactors=F)
```

#array identify<br />
```
scarinfo<-FindScar(data=data,scarfull=ref,scar=cutsite,indel.coverage="ALL",type="test",cln=8)
scarform<-INDELChangeForm(scarinfo,scarref = scarref_all,cln=4)
scarinfo$Scar<-scarform
```
#define scar pattern for each cell<br />
```
cellsinfo<-INDELIdents(scarinfo,scarref=scarref_all,scarfull=ref,scar=cutsite,method.use="umi.num",indel.coverage="ALL",cln=4)
```
tree reconstruct
```
tag<-TagDataProcess(cellsinfo$info,Cells=celltype)
treeinfo<-BuildTagTree(tag,Cells=celltype)
```

#pattern visualization <br />
