# Barcode diversity estimation
To run the following scripts. You would need to install both R and Python. 

**1. Extract and cluster barcodes from raw reads**
extract the barcode from raw read2 sequencing files using [Bartender](https://github.com/LaoZZZZZ/bartender-1.1)

```sh
$ bartender_extractor_com -f CloneSeqBarcodePlasmid200512_.conc.R2.fastq -o CloneSeq200512_extracted -q ? -p CTGA[2]CTG[2]ACT[2]GAC[2]TGA[2]CTG[2]ACT[2]GAC[2]GACT
```
Cluster the barcodes
```sh
bartender_single_com -f CloneSeq200512_extracted_barcode.txt -o CloneSeq200512_barcode -d 3 z -1
```

**2. Analyze hamming distance of barcodes**
QC of library and making sure bacteria library has high hamming distance
- they should all be around 11 and up
- when plotted, a histogram of the hamming distance should reveal a population of around 5 or higher
- run the following scripts in R

**Output hamming distance stats**
```
raw <- read.table(file = "CloneSeq200512_extracted_barcode.txt", sep = ",")
barcodes <- raw %>% dplyr::pull(V1)
subset<- barcodes[1:1000]
subset <- as.character(subset)

stats = analyse.barcodes(subset, metric = c("hamming", "seqlev", "levenshtein"), cores = detectCores()/2, cost_sub = 1, cost_indel = 1) # check stats 
stats
```

Convert hamming distance to vector for plotting

```
dm <- stringdistmatrix(subset, subset, method="hamming")
colnames(dm) <- subset
rownames(dm) <- subset
dm_v <- c(dm)
remove <- c(0)
dm_v <- dm_v[!dm_v %in% remove]
dm_v <- as.data.frame(dm_v)

ggplot(data = dm_v, aes(x = dm_v)) +
    geom_histogram(aes(y=..count../1000), binwidth = 1, fill="seagreen4",alpha=0.9) +
    labs(x = "Barcode Hamming Distance") + 
    labs(y= "Frequency(x1000)") +
    theme_classic()

```





