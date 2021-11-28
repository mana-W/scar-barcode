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
scarref_all<-ReadCutsite(cutsite,reftype="All")
celltype<-read.table("celltype.tsv",header=T,stringsAsFactors=F)
```

#array identify<br />
```
scarinfo<-FindScar(data=data,scarfull=ref,scar=cutsite,indel.coverage="All",type="test",cln=8)
scarform<-INDELChangeForm(scarinfo,scarref = scarref_all,cln=4)
scarinfo$Scar<-scarform
```
#define scar pattern for each cell<br />
```
cellsinfo<-INDELIdents(scarinfo,scarref=scarref_all,scarfull=ref,scar=cutsite,method.use="umi.num",indel.coverage="All",cln=4)
```

#pattern visualization <br />
```
IndelPlot(cellsinfo = cellsinfo,scar=cutsite,indel.coverage="All")
```
<p align="center">
<img src="https://github.com/mana-W/scar-barcode/blob/main/image/Indel_pattern.png" width = "620" height = "450" align=center />
</p >
Then <br />

tree reconstruct and plot
```
tag<-TagDataProcess(cellsinfo$info,Cells=celltype)
treeinfo<-BuildTagTree(tag,Cells=celltype)
```
```
plotinfo<-PlotTagTree(treeinfo = treeinfo,data.extract = "T",annotation = "T")
plotinfo$p
```
<p align="center">
<img src="https://github.com/mana-W/scar-barcode/blob/main/image/tree_tag_pattern.png" width = "400" height = "500" align=center />
</p >

Or <br />

```
plotinfo<-PlotTagTree(treeinfo = treeinfo,data.extract = "T",annotation = "F")
plotinfo$p
```
<p align="center">
<img src="https://github.com/mana-W/scar-barcode/blob/main/image/tree_pattern.png" width = "400" height = "500" align=center />
</p >
