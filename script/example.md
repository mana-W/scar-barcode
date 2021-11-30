# Example<br />
data is from UMI_CB/CB_UMI<br />
fa is ref file<br />
cutsite is a file define each sgRNA start and end positon<br />
celltype.tsv is a file include cell barcode and its' annotations, header: Cell.BC Cell.type

### Data importing
```
library(LinTInd)
data<-read.table("UMI_CB/CB_UMI",sep="\t",header=T)
ref<-ReadFasta("V3.fasta")
cutsite<-read.table("V3.cutSites",col.names = c("indx","start","end"))
scarref<-ReadCutsite(cutsite)
scarref_all<-ReadCutsite(cutsite,reftype="All")
celltype<-read.table("celltype.tsv",header=T,stringsAsFactors=F)
```

### Array identify<br />
Alignment

```
scarinfo<-FindIndel(data=data,scarfull=ref,scar=cutsite,indel.coverage="All",type="test",cln=8)
scarform<-IndelForm(scarinfo,scarref = scarref_all,cln=4)
scarinfo$Scar<-scarform
```
Define scar pattern for each cell<br />
```
cellsinfo<-IndelIdents(scarinfo,scarref=scarref_all,scarfull=ref,scar=cutsite,method.use="umi.num",indel.coverage="All",cln=4)
```

Pattern visualization <br />
```
IndelPlot(cellsinfo = cellsinfo,scar=cutsite,indel.coverage="All")
```
<p align="center">
<img src="https://github.com/mana-W/scar-barcode/blob/main/image/Indel_pattern.png" width = "620" height = "450" align=center />
</p >
<br />


### Indel extracted
```
tag<-TagProcess(cellsinfo$info,Cells=celltype)
```

### Tree reconstruct 
```
treeinfo<-BuildTree(tag,Cells=celltype)
```

### Visualization

**Similarity of each pair of clusters**
```
tag_dist=TagDist(tag,method = "spearman")
```
<p align="center">
<img src="https://github.com/mana-W/scar-barcode/blob/main/image/Indel.png" width = "500" height = "300" align=center />
</p >

<p align="center">
<img src="https://github.com/mana-W/scar-barcode/blob/main/image/cluster_similarity.png" width = "490" height = "450" align=center />
</p >

***Visualization for tree***
```
plotinfo<-PlotTree(treeinfo = treeinfo,data.extract = "T",annotation = "T")
plotinfo$p
```
<p align="center">
<img src="https://github.com/mana-W/scar-barcode/blob/main/image/tree_tag_pattern.png" width = "400" height = "500" align=center />
</p >

Or <br />

```
plotinfo<-PlotTree(treeinfo = treeinfo,data.extract = "T",annotation = "F")
plotinfo$p
```
<p align="center">
<img src="https://github.com/mana-W/scar-barcode/blob/main/image/tree_pattern.png" width = "400" height = "500" align=center />
</p >
