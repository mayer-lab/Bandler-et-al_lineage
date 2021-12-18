#  Created by May Ho on 18.09.19.
#  Notes for script

## The steps to preprocess FASTQ files of barcodes are as follows:
# 1. Trim and select reads
# 2. Remove short reeads and re-pair
# 3. Identify cell barcodes from R1 sequencing file/Seurat
# 4. Add the cell barcodes from step 3 to R2 reads 
# I use different versions of the script depending on whether I am trimming a virus library vs the PiggyBac library, the 

# Load program paths (change to your path)
export PATH=$PATH:/opt/bbmap # to load onto server
export PATH=$PATH:/Users/mayho/Downloads/samtools-1.9


### 1.Trimming and selecting

# Replace in1 and in2 with your R1 and R2 of your retrieved barcode file
bbduk.sh in1=CloneSeqBarcode200528__S3_R1_001.fastq.gz in2=CloneSeqBarcode200528__S3_R2_001.fastq.gz k=17 literal=GACTCTGGCTCACAAAT ktrim=r out=stdout.fq int=f skipr1 maq=15 |\
bbduk.sh in=stdin.fq literal=CTGA k=4 restrictleft=4 ktrim=l out1=trim_R1.fastq out2=trim_R2.fastq outm1=dis_R1.fastq outm2=dis_R2.fastq int=t skipr1

### 2. Remove short reads and re-pair

# Remove short reads and re-pair R1 and R2 (for libraries using piggybac)

bbduk.sh in=trim_R2.fastq out=stdout.fq minlength=37 maxlength=37 | \
repair.sh in1=stdin.fq in2=trim_R1.fastq out1=pair_R1.fastq out2=pair_R2.fastq repair

### 3. identify correct cell barcodes

# You may have to adjust the settings for umi_tools whitelist depending on the file for more info see: 
umi_tools whitelist --stdin $datadir/CloneSeqBarcode200528__S3_R1_001.fastq.gz \ ## change to your R1 file
--bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNN \
--knee-method=density \
--expect-cells=3000 \ # change number to the expected cell recovery
--plot-prefix=expect_whitelist \
--log2stderr > whitelist.txt

### 4. Add the cell barcodes from step 3 to R2 reads 

umi_tools extract --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNN \
--stdin pair_R1.fastq \
--stdout pair_R1_extracted.fastq.gz \
--read2-in pair_R2.fastq \
--read2-out=barcode_extracted.fastq \ ## final extracted barcodes
--filter-cell-barcode \
--error-correct-cell \
--whitelist=whitelist.txt;








