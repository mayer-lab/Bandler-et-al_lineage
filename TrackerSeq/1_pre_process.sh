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

# write path to input (where your files are) and output (where you want to output files) directories
datadir=/datastore_share/bioinformatics/datasets/sec/P350_210225_ED_IV_barcode/DualIndex/conc.fastq
output=/data/mayerlab/mayho/TrackerSeq_demo

### 1.Trimming and selecting

# Replace in1 and in2 with your R1 and R2 of your retrieved barcode file
bbduk.sh in1=$datadir/210303_B.conc.R1.fastq.gz in2=$datadir/210303_B.conc.R2.fastq.gz k=17 literal=GACTCTGGCTCACAAAT ktrim=r out=stdout.fq int=f skipr1 maq=15 |\
bbduk.sh in=stdin.fq literal=CTGA k=4 restrictleft=4 ktrim=l out1=$output/trim_R1.fastq out2=$output/trim_R2.fastq outm1=$output/dis_R1.fastq outm2=$output/dis_R2.fastq int=t skipr1

### 2. Remove short reads and re-pair

# Remove short reads and re-pair R1 and R2 (for libraries using piggybac)

bbduk.sh in=$output/trim_R2.fastq out=stdout.fq minlength=37 maxlength=37 | \
repair.sh in1=stdin.fq in2=$output/trim_R1.fastq out1=$output/pair_R1.fastq out2=$output/pair_R2.fastq repair

### 3. identify correct cell barcodes

# You may have to adjust the settings for umi_tools whitelist depending on the file for more info see: 
umi_tools whitelist --stdin $datadir/210303_B.conc.R1.fastq.gz \
--bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNN \
--knee-method=density \
--expect-cells=3000 \
--plot-prefix=expect_whitelist \
--log2stderr > $output/MUC28072_whitelist.txt

### 4. Add the cell barcodes from step 3 to R2 reads 

umi_tools extract --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNN \
--stdin $output/pair_R1.fastq \
--stdout $output/pair_R1_extracted.fastq.gz \
--read2-in $output/pair_R2.fastq \
--read2-out=$output/MUC28072_barcode_extracted.fastq \
--filter-cell-barcode \
--error-correct-cell \
--whitelist=$output/MUC28072_whitelist.txt








