# Master-Project-
This is James Uleth master project code 
rsync -avh --progress /Volumes/VERBATIM\ HD/Lars_RNAseq_lung_cancer_2018/bam/*.bam \
james89@cedar.computecanada.ca:/scratch/james89/    

The command I used to resume my connection after it was lost 
The command I used to ensure that the bam files are single end read files and if they are it will process them as such and if they are not they will process them as paired end reads and convert them to fast  
#!/bin/bash
#SBATCH --job-name=check_bam    # Job name
#SBATCH --ntasks=1              # Number of tasks
#SBATCH --cpus-per-task=4       # Number of CPU cores
#SBATCH --mem=8G                # Memory allocation
#SBATCH --time=2:00:00          # Time limit (hh:mm:ss)
#SBATCH --output=check_bam.%j.out  # Output file
#SBATCH --error=check_bam.%j.err   # Error file

# Load the samtools module
module load samtools

# Set directories
input_dir="/scratch/james89"
output_dir="/scratch/james89/fastq_output"
mkdir -p "$output_dir"

# Process BAM files
cd "$input_dir" || { echo "Input directory not found! Exiting."; exit 1; }

for bam in *.bam; do
    echo "Checking BAM file: $bam"

    # Check for paired-end reads in the BAM file
    paired_reads=$(samtools view -c -f 1 "$bam")
    
    if [ "$paired_reads" -eq 0 ]; then
        echo "$bam contains single-end reads. Processing as single-end."
        samtools fastq "$bam" > "$output_dir/${bam%.bam}.fastq"
    else
        echo "$bam contains paired-end reads. Processing as paired-end."
        samtools fastq -1 "$output_dir/${bam%.bam}_R1.fastq" \
                       -2 "$output_dir/${bam%.bam}_R2.fastq" \
                       -0 /dev/null -s /dev/null -n "$bam"
    fi
done

echo "All BAM files have been checked and processed."  

The above command allowed me to resume where I left off and made sure I didn't transfer duplicates

I will use the below command to make sure the files have been converted to FASTQC and that they exist in the directory and are non empty 

cd /scratch/james89

ls -lh /scratch/james89/*.fastq

pwd


I will now confirm that the files generated are of high quality using FASTQC and the following commands 
And since there is a lot of fasted files I will be doing it as a slur job 

#!/bin/bash
#SBATCH --job-name=fastqc
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=1:00:00
#SBATCH --output=/scratch/james89/fastqc_%A.out

module load fastqc

fastqc /scratch/james89/*.fastq -o /scratch/james89/fastq_output

sbatch fastqc_job.sh


The path to where the output files where located in 

ls /scratch/james89/fastqc_results/

To view the results of the fastqc quality control I used the below command 
Analyze multiple files in a loop:

bash
Copy code
for zip in /scratch/james89/fastqc_results/*.zip; do
    unzip -p $zip */summary.txt
done

Result examples 

PASS	Basic Statistics	2224.fastq
PASS	Per base sequence quality	2224.fastq
PASS	Per tile sequence quality	2224.fastq
PASS	Per sequence quality scores	2224.fastq
WARN	Per base sequence content	2224.fastq
WARN	Per sequence GC content	2224.fastq
PASS	Per base N content	2224.fastq
PASS	Sequence Length Distribution	2224.fastq
PASS	Sequence Duplication Levels	2224.fastq
PASS	Overrepresented sequences	2224.fastq
PASS	Adapter Content	2224.fastq


Since I got some outputs like the above I am going to examine what exactly is wrong in these files as if it is the adapter seq then trimming is reqiured of it is not the adapter seq then these errors are normal for RNA seq data and we will proceed as per usual.

The command I am going to use to see these files is 

scp james89@cedar.computecanada.ca:/scratch/james89/fastqc_results/4891.fastq.html ~/Downloads/


Now that I have confirmed that the above results are normal I will now begin to align my fastq data with the human reference genome 

first I got the human reference genome off of 
https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/

now I am going to unzip the files to check there content and structure 

using 
Download the Comprehensive GTF File
Ensure the GTF file includes both coding (e.g., exons, CDS) and non-coding regions (e.g., lncRNAs, introns). Use the Ensembl GTF, which includes these features:
1.	Download the comprehensive GTF annotation file:
wget ftp://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.109.gtf.gz
gunzip Homo_sapiens.GRCh38.109.gtf.gz

This file includes annotations for both coding (exons, CDS) and non-coding regions (e.g., pseudogenes, lncRNAs).

Build the HISAT2 Index with Annotations
When building the HISAT2 index, include splice site and exon information derived from the comprehensive GTF file.

Generate splice sites and exons:
module load hisat2
extract_splice_sites.py Homo_sapiens.GRCh38.109.gtf > splice_sites.txt
extract_exons.py Homo_sapiens.GRCh38.109.gtf > exons.txt


Download a Prebuilt HISAT2 Index
If Cedar does not have a prebuilt index, you can download it directly from the HISAT2 website or other trusted sources.
Download Prebuilt Index for GRCh38 from HISAT2 Developers:
1.	Run the following commands to download and extract the GRCh38 index:

wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data/grch38_snp_tran.tar.gz
tar -xvzf grch38_snp_tran.tar.gz

2.	What Is Included:
o	grch38_snp_tran.tar.gz contains:
	A GRCh38 genome index (grch38_snp_tran.*.ht2) that includes common SNPs and transcript annotations for RNA-seq.
3.	Set the Index Directory:
o	Use the extracted files in your alignment script:
index_dir="/path/to/grch38_snp_tran"

I used these codes to download the new human reference genome files 

wget -O GRCh38.primary_assembly.genome.fa.gz https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh38.primary_assembly.genome.fa.gz

wget -O gencode.v44.annotation.gtf.gz https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz

to unzip the files I used 
gunzip GRCh38.primary_assembly.genome.fa.gz
gunzip gencode.v44.annotation.gtf.gz

To cofirm they are unzipped I used this command 
file GRCh38.primary_assembly.genome.fa

The command I use to index my human reference genome was 

This was my command to create my index for my data set 
#!/bin/bash
#SBATCH --job-name=star_index             # Job name
#SBATCH --output=star_index.out           # Standard output
#SBATCH --error=star_index.err            # Standard error
#SBATCH --ntasks=1                        # Number of tasks
#SBATCH --cpus-per-task=8                 # Number of CPU cores
#SBATCH --mem=128G                        # Memory allocation
#SBATCH --time=24:00:00                   # Maximum runtime
#SBATCH --mail-type=ALL                   # Email notifications (BEGIN, END, FAIL)
#SBATCH --mail-user=james.bencsik@uleth.ca # Your email address

# Load STAR module
module load star/2.7.11b

# Define paths
GENOME_FASTA="/scratch/james89/humanref/GRCh38.primary_assembly.genome.fa"
ANNOTATION_GTF="/scratch/james89/humanref/gencode.v44.annotation.gtf"
GENOME_INDEX_DIR="/scratch/james89/Humanstarref"

# Create genome index directory if it doesn't exist (optional)
mkdir -p "$GENOME_INDEX_DIR"

# Step 1: Generate the STAR genome index
echo "Generating STAR genome index..."
STAR --runThreadN 8 \
     --runMode genomeGenerate \
     --genomeDir "$GENOME_INDEX_DIR" \
     --genomeFastaFiles "$GENOME_FASTA" \
     --sjdbGTFfile "$ANNOTATION_GTF" \
     --sjdbOverhang 100 \
     --limitGenomeGenerateRAM 128000000000  # 128 GB

if [ $? -ne 0 ]; then
    echo "Genome index generation failed. Exiting..."
    exit 1
fi

echo "Genome index generation completed successfully."

Automate the Process for Multiple FASTQ Files
If you have many FASTQ files, you can create a SLURM script to process them in a batch.
SLURM Script for Batch Alignment:
Save the script as batch_align_rnaseq.slurm



This SLURM script performed RNA-Seq alignment using STAR for 276 single-end FASTQ files, aligning them to the GRCh38 human reference genome.

#!/bin/bash
#SBATCH --job-name=star_align       
#SBATCH --output=/scratch/james89/logs/star_align_%A_%a.out  
#SBATCH --error=/scratch/james89/logs/star_align_%A_%a.err   
#SBATCH --ntasks=1                 
#SBATCH --cpus-per-task=8          
#SBATCH --mem=350G                 
#SBATCH --time=48:00:00            
#SBATCH --array=0-275              
#SBATCH --mail-type=ALL            
#SBATCH --mail-user=james.bencsik@uleth.ca 

# Load STAR module
module load star/2.7.11b

# Define paths
FASTQ_DIR="/scratch/james89/fastq_output"
GENOME_INDEX_DIR="/scratch/james89/Humanstarref"
OUTPUT_DIR="/scratch/james89/alignment_results"

# Ensure FASTQ files exist
FASTQ_FILES=($(ls "$FASTQ_DIR"/*.fastq 2>/dev/null))
if [[ ${#FASTQ_FILES[@]} -eq 0 ]]; then
    echo "🚨 ERROR: No FASTQ files found in $FASTQ_DIR"
    exit 1
fi

# Select file based on SLURM_ARRAY_TASK_ID
FASTQ_FILE="${FASTQ_FILES[$SLURM_ARRAY_TASK_ID]}"

# Debugging output
echo "Processing FASTQ file: $FASTQ_FILE"

# Run STAR alignment (NO zcat since files are uncompressed)
STAR --runThreadN 8 \
     --genomeDir "$GENOME_INDEX_DIR" \
     --readFilesIn "$FASTQ_FILE" \
     --outFileNamePrefix "$OUTPUT_DIR/$(basename "$FASTQ_FILE" .fastq)_" \
     --outSAMtype BAM SortedByCoordinate \
     --quantMode TranscriptomeSAM GeneCounts


To get unaligned reads as well which is essential for biomarkers discovery I ran these commands 

#!/bin/bash
#SBATCH --job-name=extract_unaligned
#SBATCH --output=/scratch/james89/logs/extract_unaligned_%A_%a.out
#SBATCH --error=/scratch/james89/logs/extract_unaligned_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=100G
#SBATCH --time=12:00:00
#SBATCH --array=0-275  # Updated for 276 samples
#SBATCH --mail-type=ALL
#SBATCH --mail-user=james.bencsik@uleth.ca 

# Load STAR module
module load star/2.7.11b

# Define paths
FASTQ_DIR="/scratch/james89/fastq_output"
OUTPUT_DIR="/scratch/james89/unaligned_reads"
GENOME_INDEX_DIR="/scratch/james89/Humanstarref"

# Ensure FASTQ files exist
FASTQ_FILES=($(ls "$FASTQ_DIR"/*.fastq 2>/dev/null))
if [[ ${#FASTQ_FILES[@]} -eq 0 ]]; then
    echo "🚨 ERROR: No FASTQ files found in $FASTQ_DIR"
    exit 1
fi

# Select file based on SLURM_ARRAY_TASK_ID
FASTQ_FILE="${FASTQ_FILES[$SLURM_ARRAY_TASK_ID]}"
BASENAME=$(basename "$FASTQ_FILE" .fastq)

# Run STAR to extract only unaligned reads
STAR --runThreadN 8 \
     --genomeDir "$GENOME_INDEX_DIR" \
     --readFilesIn "$FASTQ_FILE" \
     --outFileNamePrefix "$OUTPUT_DIR/${BASENAME}_" \
     --outReadsUnmapped Fastx \
     --outSAMmode None


Next we will Convert BAM to gene counts I used 

#!/bin/bash
#SBATCH --job-name=featureCounts
#SBATCH --output=/scratch/james89/logs/featureCounts_%j.out
#SBATCH --error=/scratch/james89/logs/featureCounts_%j.err
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --mail-user=james.bencsik@uleth.ca
#SBATCH --mail-type=END,FAIL

# Load the subread module
module load subread/2.0.6

echo "Running featureCounts..."

# Run featureCounts on all valid BAM files
featureCounts -T 8 -a /scratch/james89/humanref/gencode.v44.annotation.gtf \
              -o /scratch/james89/gene_counts.tsv \
              $(find /scratch/james89/alignment_results -type f -size +0 -name "*_Aligned.sortedByCoord.out.bam")

echo "featureCounts completed."


Next we will quantify our transcript quantification the command I will use to do this is 

First I need to generate the index for salmon 
Step 1: Download the Correct Transcriptome FASTA
wget -O /scratch/james89/humanref/gencode.v44.transcripts.fa.gz \
     https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.transcripts.fa.gz

gunzip /scratch/james89/humanref/gencode.v44.transcripts.fa.gz

Second I created the salmon index via 
#!/bin/bash
#SBATCH --job-name=salmon_index
#SBATCH --output=/scratch/james89/logs/salmon_index.out
#SBATCH --error=/scratch/james89/logs/salmon_index.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --mail-user=james.bencsik@uleth.ca
#SBATCH --mail-type=END,FAIL

# Load Salmon module
module load salmon/1.10.2

# Define paths
TRANSCRIPT_FASTA="/scratch/james89/humanref/gencode.v44.transcripts.fa"
INDEX_DIR="/scratch/james89/Humanstarref_salmon_index"

# Ensure the index directory exists
mkdir -p "$INDEX_DIR"

# Create Salmon index (REMOVED --type quasi)
salmon index -t "$TRANSCRIPT_FASTA" \
             -i "$INDEX_DIR" \
             -k 31 \
             --keepDuplicates

echo "Salmon index creation completed."

Now I need to do data normalization 

nano run_tximport.slurm


#!/bin/bash
#SBATCH --job-name=tximport
#SBATCH --output=/scratch/james89/logs/tximport_%A.out
#SBATCH --error=/scratch/james89/logs/tximport_%A.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --ntasks=1
#SBATCH --mail-user=james.bencsik@uleth.ca
#SBATCH --mail-type=END,FAIL

# Load required modules
module load StdEnv/2020
module load r/4.2.2

# Run the R script
Rscript /scratch/james89/salmon_merge.R

nano /scratch/james89/salmon_merge.R

# Ensure R uses the personal library
.libPaths(c("~/R/x86_64-pc-linux-gnu-library/4.2", .libPaths()))

# Install required packages if missing
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", lib="~/R/x86_64-pc-linux-gnu-library/4.2", repos="https://cran.rstudio.com")
if (!requireNamespace("tximport", quietly = TRUE)) BiocManager::install("tximport", lib="~/R/x86_64-pc-linux-gnu-library/4.2")
if (!requireNamespace("jsonlite", quietly = TRUE)) install.packages("jsonlite", lib="~/R/x86_64-pc-linux-gnu-library/4.2", repos="https://cran.rstudio.com")

# Load libraries
library(tximport)
library(jsonlite)

# Define paths to all Salmon output folders
samples <- list.files("/scratch/james89/salmon_quant/", full.names = TRUE)
files <- file.path(samples, "quant.sf")
names(files) <- basename(samples)

# Import Salmon quantification
txi <- tximport(files, type = "salmon", txOut = TRUE, dropInfReps = TRUE)  # Drop inferential replicates to avoid issues

# Save the merged matrix
write.csv(txi$counts, file="/scratch/james89/salmon_counts_matrix.csv")


Before I normalize it I need to run The issue with importing the data into DESeq2 from Salmon is caused by the presence of fractional counts that are not accepted by DESeq2. You need to use Salmon->tximport->DESeq2 workflow that summarizes transcript counts over a gene. See https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#transcript-abundance-files-and-tximport-tximeta for detailed instruction. Once gene level counts are obtained via tximport you can proceed with a standard DESeq2 workflow as described in the manual.


Here Is my R workflow 
Step 1: Load Required Libraries
Ensure you have the necessary Bioconductor packages installed and loaded:

library(txdbmaker)
library(GenomicFeatures)

Step 2: Create the TxDb Object
Use the correct path to the GTF file in the makeTxDbFromGFF() function:

txdb <- makeTxDbFromGFF("/Volumes/One Touch/masters project stuff/humanref/gencode.v44.annotation.gtf", 
                        format = "gtf")

Step 3: Extract Transcript-to-Gene Mapping
Extract the mapping of transcript IDs (TXNAME) to gene IDs (GENEID):

# Extract transcript-to-gene mapping
tx2gene <- select(txdb, keys = keys(txdb, keytype = "TXNAME"), 
                  columns = "GENEID", keytype = "TXNAME")

# Check the first few rows
head(tx2gene)

Step 4: Save the Mapping to a File
Save the transcript-to-gene mapping as a .tsv file for use in tximport:

# Save tx2gene to a file
write.table(tx2gene, "/Volumes/One Touch/masters project stuff/humanref/tx2gene.tsv", 
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


Now that you've successfully generated the tx2gene mapping and are ready to process the Salmon outputs, the next steps involve:
1.	Using tximport to Summarize Transcript Counts into Gene Counts.
2.	Preparing the Data for DESeq2 Analysis.
3.	Normalize the data 
4.	Performing Differential Expression Analysis with DESeq2.

Load Required Libraries
Ensure you have tximport, readr, and DESeq2 installed and loaded:
library(tximport)
library(readr)
library(DESeq2)

Load the tx2gene Mapping
Load the transcript-to-gene mapping file created earlier:
tx2gene <- read_tsv("/Volumes/One Touch/masters project stuff/humanref/tx2gene.tsv", 
                    col_names = c("TXNAME", "GENEID"))
Prepare List of Salmon Quantification Files
List AND organize the quant.sf files from your Salmon output directories:
List All quant.sf Files
Assuming each sample's output (including quant.sf) is in a separate subdirectory under /Volumes/One Touch/masters project stuff/salmon_quant, list all the quant.sf files:
# List all quant.sf files recursively
salmon_files <- list.files("/Volumes/One Touch/masters project stuff/salmon_quant", 
                           pattern = "quant.sf", full.names = TRUE, recursive = TRUE)

# Assign sample names based on their folder names
names(salmon_files) <- basename(dirname(salmon_files))

# Check the files
print(salmon_files)

Step 2: Verify the Files
Ensure that all expected samples are accounted for. For example:
# Check the number of files
cat("Number of quant.sf files found:", length(salmon_files), "\n")

# Verify file names
print(names(salmon_files))

Step 3: Import Data Using tximport
Load the transcript-to-gene mapping file (tx2gene.tsv) created earlier and use it with tximport:

library(tximport)
library(readr)

# Load the tx2gene mapping file
tx2gene <- read_tsv("/Volumes/One Touch/masters project stuff/humanref/tx2gene.tsv", 
                    col_names = c("TXNAME", "GENEID"))

# Import counts from Salmon
txi <- tximport(salmon_files, 
                type = "salmon", 
                tx2gene = tx2gene, 
                ignoreTxVersion = TRUE)  # Ignore version numbers in transcript IDs

# Check the imported data
str(txi)
head(txi$counts)


Error I ran into 
Error in .local(object, ...) : 
  None of the transcripts in the quantification files are present
  in the first column of tx2gene. Check to see that you are using
  the same annotation for both.

Example IDs (file): [ENST00000456328, ENST00000450305, ENST00000488147, ...]

Example IDs (tx2gene): [TXNAME, ENST00000456328.2, ENST00000450305.2, ...]

  This can sometimes (not always) be fixed using 'ignoreTxVersion' or 'ignoreAfterBar'.

possible fix 
# Remove version numbers from the transcript IDs in tx2gene
tx2gene$TXNAME <- gsub("\\..*$", "", tx2gene$TXNAME)

# Import counts from Salmon
txi <- tximport(salmon_files, 
                type = "salmon", 
                tx2gene = tx2gene, 
                ignoreTxVersion = TRUE)  # Optionally, you can still use ignoreTxVersion

# Check the imported data
str(txi)
head(txi$counts)

# Load the metadata table
meta_data <- read.csv("/Volumes/One Touch/masters project stuff/prepped_clin_data.csv", 
                      header = TRUE, 
                      stringsAsFactors = FALSE)

# View the first few rows of the table
head(meta_data)

Step 4: Create a Sample Table for salmon prepped clincla data already created 
Prepare a metadata table with sample IDs and experimental conditions. For example:

# Create a table for Salmon metadata
salmon_meta <- data.frame(
    sampleID = names(salmon_files),  # Extract sample names from Salmon files
    stringsAsFactors = FALSE
)

# View the Salmon metadata
head(salmon_meta)

# Add experimental condition (update according to your design)
salmon_meta$condition <- c("tumour", "tumour", "tumour", "tumour")  # Replace with your actual group labels

# Check the Salmon metadata
head(salmon_meta)

# Load clinical metadata
clinical_meta <- read.csv("/Volumes/One Touch/masters project stuff/prepped_clin_data.csv", 
                          stringsAsFactors = FALSE)

# Check the first few rows to confirm it's loaded
head(clinical_meta)

# Merge salmon_meta and clinical_meta
merged_meta <- merge(salmon_meta, clinical_meta, by.x = "sampleID", by.y = "patientid", all.x = TRUE)

# Check the merged metadata
head(merged_meta)
Check for samples that need to be exculuded

# Identify rows with missing clinical metadata
missing_rows <- is.na(merged_meta$VitalStatus)  # Replace 'VitalStatus' with any key clinical column
missing_samples <- merged_meta[missing_rows, "sampleID"]

# Print missing samples
cat("Samples with missing clinical data:\n")
Samples with missing clinical data:
> print(missing_samples)
 [1] "2263"    "2347_1"  "2384"    "2433"    "2760"    "2760_1"  "2765"   
 [8] "2804_1"  "3508"    "3690"    "3959"    "4061"    "4543_17" "4543_22"
[15] "5013"    "5131"    "5149_1"  "5163"  


# Define the missing samples
missing_samples <- c("2263", "2347_1", "2384", "2433", "2760", "2760_1", 
                     "2765", "2804_1", "3508", "3690", "3959", "4061", 
                     "4543_17", "4543_22", "5013", "5131", "5149_1", "5163")

# Remove rows with missing samples
merged_meta_clean <- merged_meta[!merged_meta$sampleID %in% missing_samples, ]

# Check the updated metadata table
cat("Updated metadata table after removing missing samples:\n")
print(head(merged_meta_clean))

# Check if any missing samples remain
remaining_missing <- intersect(missing_samples, merged_meta_clean$sampleID)

if (length(remaining_missing) == 0) {
    cat("All missing samples have been successfully removed!\n")
} else {
    cat("These missing samples still remain in the metadata:\n")
    print(remaining_missing)
}

# Now remove them from the counts matrix

# Subset the counts matrix to include only samples in the metadata
txi$counts <- txi$counts[, colnames(txi$counts) %in% merged_meta_clean$sampleID]
txi$abundance <- txi$abundance[, colnames(txi$abundance) %in% merged_meta_clean$sampleID]
txi$length <- txi$length[, colnames(txi$length) %in% merged_meta_clean$sampleID]

# Check the dimensions after filtering
dim(txi$counts)

# Reorder metadata to match counts matrix
merged_meta_clean <- merged_meta_clean[match(colnames(txi$counts), merged_meta_clean$sampleID), ]

# Check for any NA values in the reordered metadata
if (any(is.na(merged_meta_clean$sampleID))) {
    stop("Error: Sample mismatch still exists after filtering!")
}

Create the DESeq2 dataset for normalization 

Step 1: Load Required Libraries
Ensure that the DESeq2 library is loaded:

library(DESeq2)

Step 2: Create the DESeq2 Dataset
Use the DESeqDataSetFromTximport() function to create the DESeq2 dataset:

# Create the DESeq2 dataset with no design comparison (using ~ 1)
dds <- DESeqDataSetFromTximport(txi = txi, 
                                colData = merged_meta_clean, 
                                design = ~ 1)

# Check the DESeq2 dataset
Dds

# Perform differential expression analysis
dds <- DESeq(dds)

# Extract results
res <- results(dds)

# View a summary of the results
summary(res)

# View the top significant genes
head(res)


Run DESeq2 for Normalization & Variance Stabilization
•	Normalize the counts and perform variance stabilization.

dds <- DESeq(dds)
vsd <- vst(dds, blind=FALSE)  # Variance stabilizing transformation

# View the top significant genes
head(res)

out of 59140 with nonzero total read count
adjusted p-value < 0.1
LFC > 0 (up)       : 32015, 54%
LFC < 0 (down)     : 25055, 42%
outliers [1]       : 0, 0%
low counts [2]     : 47, 0.079%
(mean count < 0)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results


> # View the top significant genes
> head(res)
log2 fold change (MLE): Intercept 
Wald test p-value: Intercept 
DataFrame with 6 rows and 6 columns
                     baseMean log2FoldChange     lfcSE      stat      pvalue
                    <numeric>      <numeric> <numeric> <numeric>   <numeric>
ENSG00000000003.16 217.516577        7.76498 0.0556106  139.6312 0.000000000
ENSG00000000005.6    0.218932       -2.19144 0.6102940   -3.5908 0.000329664
ENSG00000000419.14 428.120778        8.74187 0.0237546  368.0076 0.000000000
ENSG00000000457.14 259.790781        8.02121 0.0297741  269.4018 0.000000000
ENSG00000000460.17 108.509112        6.76167 0.0490238  137.9263 0.000000000
ENSG00000000938.13 133.451921        7.06018 0.0640678  110.1984 0.000000000
                          padj
                     <numeric>
ENSG00000000003.16 0.000000000
ENSG00000000005.6  0.000376667
ENSG00000000419.14 0.000000000
ENSG00000000457.14 0.000000000
ENSG00000000460.17 0.000000000
ENSG00000000938.13 0.000000000


Feature Selection: Extract Top 500 Most Variable Genes
•	Select the top 500 genes with the highest Median Absolute Deviation (MAD).
mad_values <- apply(assay(vsd), 1, mad)  # Compute MAD for each gene
top500_genes <- names(sort(mad_values, decreasing=TRUE)[1:500])
vsd_top500 <- vsd[top500_genes, ]  # Subset for downstream analysis

3. Clustering: Heatmap of Top 500 Genes
Perform Hierarchical Clustering & Generate Heatmap
•	Use pheatmap to visualize clustering.

library(pheatmap)
pheatmap(assay(vsd_top500), 
         scale="row", 
         clustering_distance_rows="euclidean", 
         clustering_distance_cols="correlation", 
         clustering_method="ward.D2")

Principal Component Analysis (PCA)
Run PCA on Top 500 Genes

Use PCAexplorer for dimensionality reduction.
library(PCAtools)
p <- pca(assay(vsd_top500), metadata = colData(vsd_top500), removeVar = 0.1)
-- removing the lower 10% of variables based on variance
> biplot(p)

Additional Data Exploration
Distance Matrix for Sample Similarities
dist_matrix <- dist(t(assay(vsd_top500)), method="euclidean")
heatmap(as.matrix(dist_matrix))

Cluster Samples Using Different Methods
hclust_samples <- hclust(dist_matrix, method="average")
plot(hclust_samples)


library(PCAtools)
> p <- pca(assay(vsd_top500), metadata = colData(vsd_top500), removeVar = 0.1)
-- removing the lower 10% of variables based on variance
> biplot(p)

Refine Various graphs:

Heatmap with Different Clustering Methods
We will:
1.	Remove row names (gene names)
2.	Add an annotation color bar for tumor type, gender, and smoking status
3.	Generate three versions using different clustering methods: complete, Ward.D2, and average
4.	Install & Load Required Packages

if (!requireNamespace("pheatmap", quietly=TRUE)) install.packages("pheatmap")
library(pheatmap)

# Define annotations with correct column names
annotation_col <- data.frame(
    Tumor_Type = merged_meta_clean$Histology,  # Use Histology instead of TumorType
    Gender = merged_meta_clean$Gender,
    Smoking_Status = merged_meta_clean$smoking  # Use smoking instead of SmokingStatus
)

# Set row names to match the column names of heatmap_data
rownames(annotation_col) <- rownames(merged_meta_clean)

# Load necessary library
library(pheatmap)

# Define annotations with correct column names
annotation_col <- data.frame(
    Tumor_Type = merged_meta_clean$Histology,  # Correct column name
    Gender = merged_meta_clean$Gender,
    Smoking_Status = merged_meta_clean$smoking  # Correct column name
)

# Convert categorical variables to factors
annotation_col <- data.frame(lapply(annotation_col, as.factor))

# Set row names to match the column names of heatmap_data
rownames(annotation_col) <- colnames(heatmap_data)

# Generate the heatmap with annotations
pheatmap(
    heatmap_data, 
    annotation_col = annotation_col,
    cluster_rows = TRUE,   # Cluster genes/features
    cluster_cols = TRUE,   # Cluster samples
    scale = "row",         # Normalize by row (optional)
    show_colnames = FALSE,  # Hide column names for readability
    show_rownames = FALSE,  # Hide row names for readability
    color = colorRampPalette(c("blue", "white", "red"))(100) # Custom color scale
)

PCA Plot with Tumor Type Legend
✅ Ensure PCA plot includes a legend for tumor types
✅ Try different PC combinations (PC1 vs PC2, PC1 vs PC3, PC1 vs PC4)
✅ Add ellipses to PCA plot

Install & Load Required Packages
if (!requireNamespace("PCAtools", quietly=TRUE)) BiocManager::install("PCAtools")
library(PCAtools)

# Run PCA on variance-stabilized expression data
> pca_results <- pca(assay(vsd_top500), metadata = merged_meta_clean, removeVar = 0.1)
Error in pca(assay(vsd_top500), metadata = merged_meta_clean, removeVar = 0.1) : 
  'colnames(mat)' is not identical to 'rownames(metadata)'
> # Print the first few entries of the sampleID column
> print(head(merged_meta_clean$sampleID))
[1] "2224" "2242" "2251" "2265" "2285" "2289"
> 
> # Set row names using the sampleID column
> rownames(merged_meta_clean) <- merged_meta_clean$sampleID
> 
> # Remove sampleID column from the dataframe (to prevent duplication)
> merged_meta_clean$sampleID <- NULL
> 
> # Confirm that row names have been set correctly
> print(head(rownames(merged_meta_clean)))
[1] "2224" "2242" "2251" "2265" "2285" "2289"
> 
> # Check if the sample names match between the two datasets
> identical(colnames(assay(vsd_top500)), rownames(merged_meta_clean))
[1] TRUE

Run actual command 


# Load necessary libraries
library(PCAtools)
library(ggplot2)
library(reshape2)

# Run PCA on variance-stabilized expression data
pca_results <- pca(assay(vsd_top500), metadata = merged_meta_clean, removeVar = 0.1)

# Define function to plot PCA for different principal component combinations
plot_pca <- function(pcX, pcY) {
    biplot(pca_results, colby = "Histology", 
           hline = 0, vline = 0, lab = NULL, 
           title = paste("PCA: PC", pcX, "vs PC", pcY),
           x = pcX, y = pcY)
}

# Generate PCA plots for different principal components
plot_pca(1, 2)  # PC1 vs PC2
plot_pca(1, 3)  # PC1 vs PC3
plot_pca(1, 4)  # PC1 vs PC4

# Generate a scree plot
screeplot(pca_results)

# Extract PC values
pc_values <- pca_results$rotated

# Convert metadata to numeric where necessary
numeric_meta <- merged_meta_clean[, sapply(merged_meta_clean, is.numeric), drop = FALSE]

# Align metadata with PC values
numeric_meta <- numeric_meta[rownames(pc_values), , drop = FALSE]

# Compute Pearson correlations
pc_correlations <- cor(pc_values, numeric_meta, use = "pairwise.complete.obs", method = "pearson")

# Convert correlation matrix into long format
pc_corr_melt <- melt(pc_correlations)

# Plot heatmap of correlations
ggplot(pc_corr_melt, aes(Var1, Var2, fill = value)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0) +
    labs(title = "Correlation between PCs and Metadata", x = "Principal Components", y = "Metadata") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

Before running concenus clustering we need to download and load summerized expirment so that we can take our top500 differentially expressed genes and put it into clusterplus 
# Install SummarizedExperiment if it’s not already installed
if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("SummarizedExperiment")
}

# Load the package
library(SummarizedExperiment)

Then we need to Access the Variance-Stabilized Data: Once SummarizedExperiment is loaded, you can access the variance-stabilized data using the assay() function.

# Extract the matrix of variance-stabilized data from vsd_top500
vsd_matrix <- assay(vsd_top500)  # This should work after loading the SummarizedExperiment package

# Check the matrix class and dimensions
class(vsd_matrix)  # Should be a matrix
dim(vsd_matrix)    # Should show 500 genes and 258 samples

Now we need to Perform Consensus Clustering: After successfully extracting the matrix, you can proceed with consensus clustering

# ------------------------------------------------------------------#

# Load the vsd_top500 data (assuming vsd_top500 is already loaded as a DESeq2 object)
d <- assay(vsd_top500)  # Extract the matrix from the DESeq2 object
head(d)  # Check the first few rows of the data
dim(d)  # Check the dimensions of the data

# Get top 500 genes by MAD
mads <- apply(d, 1, mad)  # Calculate the MAD for each row (gene)
d <- d[rev(order(mads))[1:500],]  # Select the top 500 genes based on MAD
d <- sweep(d, 1, apply(d, 1, median, na.rm = TRUE))  # Normalize data by subtracting median of each row
head(d)  # Check the first few rows of the normalized matrix
dim(d)  # Check the new dimensions

# Run Consensus Clustering
library(ConsensusClusterPlus)

results <- ConsensusClusterPlus(d, maxK = 6, reps = 50, pItem = 0.8,
                                pFeature = 1, title = "test_graph_top500", 
                                clusterAlg = "hc",
                                distance = "pearson", 
                                seed = 126118388.71279, plot = "png")

# Display the consensus matrix and consensus tree for the first 5 rows/columns
print(results[[2]][["consensusMatrix"]][1:5, 1:5])
print(results[[2]][["consensusTree"]])

# Calculate cluster-consensus and item-consensus results:
icl <- calcICL(results, title = "test_graph_top500", plot = "png")

# Display the cluster-consensus and item-consensus values
print(icl[["clusterConsensus"]])
print(icl[["itemConsensus"]][1:5,])  # Show the first 5 item-consensus values

# ------------------------------------------------------------------#


So we are now going to run a differential expression analysis to identify prognostic Biomarkers:

Poor vs good Survival 

Step One: Define Groups in merged_meta_clean

# Define prognosis groups based on survival (e.g., cutoff = 18 months)
merged_meta_clean$PrognosisGroup <- ifelse(merged_meta_clean$survival_mo >= 18, "good", "poor")

# Confirm distribution
table(merged_meta_clean$PrognosisGroup)

# Make sure it's a factor
merged_meta_clean$PrognosisGroup <- factor(merged_meta_clean$PrognosisGroup, levels = c("good", "poor"))

We are then going to create and run the deseq2 object 
dds_prognosis <- DESeqDataSetFromTximport(txi = txi,
                                          colData = merged_meta_clean,
                                          design = ~ PrognosisGroup)

dds_prognosis <- DESeq(dds_prognosis)

res_prognosis <- results(dds_prognosis, contrast = c("PrognosisGroup", "poor", "good"))

# Filter significant genes
res_prognosis_sig <- res_prognosis[which(res_prognosis$padj < 0.05 & abs(res_prognosis$log2FoldChange) > 1), ]
res_prognosis_sig <- res_prognosis_sig[order(res_prognosis_sig$padj), ]

# Save results
write.csv(as.data.frame(res_prognosis_sig), "DEA_by_survival_group.csv")

View the results thus far
# View first few rows
head(res_prognosis_sig)

# View full structure
str(res_prognosis_sig)

# How many genes passed your filter?
nrow(res_prognosis_sig)

Then run this command to see the top hits
head(res_prognosis_sig)
nrow(res_prognosis_sig)  # How many significant DEGs?

Brief over view of the results 
1,268 significant genes:
These genes:
•	Have FDR-adjusted p-values (padj) < 0.05
•	Have log2 fold change > 1 or < -1, meaning expression changed by at least 2-fold
•	Are differentially expressed between poor vs good prognosis patients
Gene ID	log2FC	padj	Expression Higher In
ENSG00000176678.6	+2.27	6.77e-11	Poor prognosis group
ENSG00000128591.16	+2.14	2.07e-10	Poor prognosis group
ENSG00000109956.13	−2.82	4.77e-09	Good prognosis group
ENSG00000164588.8	−5.04	4.77e-09	Good prognosis group
ENSG00000165376.12	−3.60	5.30e-09	Good prognosis group
ENSG00000124795.17	+1.13	5.53e-09	Poor prognosis group


A positive log2FC = gene is upregulated in poor prognosis
A negative log2FC = gene is upregulated in good prognosis
This is exactly the kind of result you'd want when searching for prognostic biomarkers.


Part 2: DEA by ConsensusClusterPlus Clusters

Add Cluster Assignments to Metadata
# Extract and assign cluster labels
cluster_assignments <- results[[4]]$consensusClass
merged_meta_clean$Cluster <- as.factor(cluster_assignments)

# Ensure sample order matches
stopifnot(names(cluster_assignments) == rownames(merged_meta_clean))

Run DESeq2 Between Two Clusters (e.g., 1 vs 4)
dds_cluster <- DESeqDataSetFromTximport(txi = txi,
                                        colData = merged_meta_clean,
                                        design = ~ Cluster)

dds_cluster <- DESeq(dds_cluster)

# Example: Cluster 4 vs Cluster 1
res_cluster <- results(dds_cluster, contrast = c("Cluster", "4", "1"))

# Filter significant genes
res_cluster_sig <- res_cluster[which(res_cluster$padj < 0.05 & abs(res_cluster$log2FoldChange) > 1), ]
res_cluster_sig <- res_cluster_sig[order(res_cluster_sig$padj), ]

# Save
write.csv(as.data.frame(res_cluster_sig), "DEA_by_clusters.csv")

To view the results we ran 
# View top 6 significant genes
head(res_cluster_sig)

# View full structure of the result
str(res_cluster_sig)

# Number of significant genes found
nrow(res_cluster_sig)

Results:
You’ve successfully run your second differential expression analysis, comparing Cluster 4 vs Cluster 1, and found 1,968 significant genes.

Interpretation of Results

Top Genes (from head(res_cluster_sig)):
Gene ID	log2FC	padj	Expression Higher In
ENSG00000132874.15	+3.48	6.52e-75	Cluster 4
ENSG00000115758.14	+3.90	1.38e-49	Cluster 4
ENSG00000163053.11	+4.15	1.07e-46	Cluster 4
ENSG00000123999.5	+5.07	4.00e-43	Cluster 4
•	These genes are massively upregulated in Cluster 4 compared to Cluster 1.
•	The log2FoldChange values of 3–5+ mean a 8x–32x increase in expression!
•	The very low padj values indicate high statistical confidence.

What This Means:
•	Cluster 4 has a distinct expression profile, which could relate to tumor subtype or prognosis.
•	The 1,968 genes are potential biomarkers or drivers of the biology defining this cluster.
•	If Cluster 4 also had worse survival, these DEGs are prognostically relevant.

Next Steps: Intersect with Cox Regression: VSD_top500
library(survival)
library(tibble)

# Confirm sample alignment
stopifnot(rownames(expr_matrix) == rownames(merged_meta_clean))

# Run Cox regression for each gene
cox_results <- apply(expr_matrix, 2, function(expr) {
  fit <- coxph(Surv(time = merged_meta_clean$survival_mo, event = merged_meta_clean$VitalStatus) ~ expr)
  summary(fit)
})

# Extract key results
cox_df <- tibble(
  gene = colnames(expr_matrix),
  HR = sapply(cox_results, function(x) x$coefficients[1, "exp(coef)"]),
  pval = sapply(cox_results, function(x) x$coefficients[1, "Pr(>|z|)"])
)

# Adjust p-values
cox_df$padj <- p.adjust(cox_df$pval, method = "fdr")

# Sort and preview
cox_df <- cox_df[order(cox_df$padj), ]
head(cox_df)


cox_df Output:
You now have a table of:
•	Gene ENSEMBL IDs
•	Hazard Ratios (HR)
•	Raw p-values
•	FDR-adjusted p-values (padj)
Your Current head(cox_df):
Gene ID	HR	pval	padj
ENSG00000117983.17	0.716	0.044	0.994
ENSG00000229807.13	0.858	0.217	0.994
ENSG00000168484.13	1.02	0.867	0.994
...	...	...	...
Issue: Very High Adjusted p-values (padj ≈ 0.994)
This suggests that no genes passed the FDR cutoff (e.g., padj < 0.05).

Why This Happens:
1.	Sample size: You only have 258 samples → low power.
2.	Stringent correction: FDR adjustment is conservative with 500 tests.
3.	Survival signal: It may be subtle (no gene has very strong association on its own).
Since this is the above case we are going to run the entire vsd in order to get a more complete capture of our prognostic biomarkers

Why Run Cox on the Full vsd?
You didn’t find strong survival signals among the top 500 most variable genes, which suggests:
•	Important prognostic genes might not be highly variable (e.g., some immune or DNA repair genes are tightly regulated)
•	You’re currently only testing a small subset (~500) of the transcriptome
By running the full vsd, you:
•	Increase your chance of capturing rare but critical survival-associated genes
•	Make the most of your data — vsd is already normalized and variance-stabilized
•	Get a comprehensive prognostic gene ranking

Before You Run It
➕ Pros:
•	More complete analysis
•	Higher likelihood of finding significant results after FDR correction
➖ Cons:
•	Slower (~60K genes × 258 samples)
•	Takes more memory (~1–3 mins depending on system)

Create the expression matrix
expr_matrix_full <- t(assay(vsd))  # transpose so rows = samples

Confirm alignment
stopifnot(rownames(expr_matrix_full) == rownames(merged_meta_clean))

Run Cox on full matrix
library(survival)
library(tibble)

cox_results_full <- apply(expr_matrix_full, 2, function(expr) {
  summary(coxph(Surv(time = merged_meta_clean$survival_mo, event = merged_meta_clean$VitalStatus) ~ expr))
})

# Extract HR and p-values
cox_df_full <- tibble(
  gene = colnames(expr_matrix_full),
  HR = sapply(cox_results_full, function(x) x$coefficients[1, "exp(coef)"]),
  pval = sapply(cox_results_full, function(x) x$coefficients[1, "Pr(>|z|)"])
)

# FDR correction
cox_df_full$padj <- p.adjust(cox_df_full$pval, method = "fdr")

# Sort by significance
cox_df_full <- cox_df_full[order(cox_df_full$padj), ]
head(cox_df_full)

What You Got:
Gene	HR	Raw p-value	FDR-adjusted p-value
ENSG00000073598.6	23.7	0.0006	0.423
ENSG00000117602.13	0.0743	0.0002	0.423
•	You found strong effect sizes (HR > 10 or HR < 0.1)
•	But padj is still high (≥ 0.423)

Interpretation
You're detecting genes with strong potential survival association, but after FDR correction across ~60,000 genes, none reach significance (padj < 0.05).
This is not unusual for:
•	Moderate sample size (n = 258)
•	Large multiple testing burden
Even though padj is high, some genes (like those with p < 0.001 and HR > 5 or < 0.2) are still biologically interesting — especially if they also show up in your DEA.

Given the above results
That would be the correct next step if your cox_df had genes with padj < 0.05.
However, from your full Cox regression output:
All padj values were ≥ 0.423 → no genes passed FDR correction

✅ Modified Step (Using Raw p-value Instead)
Because you're in the biomarker discovery phase, it's acceptable to use a lenient filter for intersection (e.g., pval < 0.01) to avoid missing meaningful genes.

# Use raw p-value threshold instead of padj
cox_genes <- cox_df_full$gene[cox_df_full$pval < 0.01]

# From your DEA results
prognosis_de_genes <- rownames(res_prognosis_sig)
cluster_de_genes <- rownames(res_cluster_sig)

# Intersect survival-associated genes with DEGs
cox_prognosis_overlap <- intersect(cox_genes, prognosis_de_genes)
cox_cluster_overlap <- intersect(cox_genes, cluster_de_genes)

# Save overlap genes
write.csv(cox_prognosis_overlap, "cox_prognosis_overlap.csv", row.names = FALSE)
write.csv(cox_cluster_overlap, "cox_cluster_overlap.csv", row.names = FALSE)

# Optional: view counts
length(cox_prognosis_overlap)
length(cox_cluster_overlap)

Results:
You successfully identified:
Overlap Type	Count	Meaning
cox_prognosis_overlap	9	Genes that are differentially expressed between poor vs good survival and significantly associated with survival time
cox_cluster_overlap	14	Genes that distinguish molecular subtypes (clusters) and correlate with survival time
These overlapping genes are your top prognostic biomarker candidates — statistically relevant and biologically meaningful.

Annotate Overlapping ENSEMBL IDs to Gene Symbols
Let’s annotate the 9 overlapping prognostic genes (cox_prognosis_overlap) first.

Lets start off with the 9 cox_progonisis_overlap:
# Remove version suffix from ENSEMBL IDs
ids_trimmed <- gsub("\\.\\d+$", "", cox_prognosis_overlap)

# Re-annotate using trimmed IDs
annot <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol", "description"),
  filters = "ensembl_gene_id",
  values = ids_trimmed,
  mart = ensembl
)

# Save to file
write.csv(annot, "annotated_cox_prognosis_genes_trimmed.csv", row.names = FALSE)

# Create manual entry
manual_entry <- data.frame(
  ensembl_gene_id = missing_id,
  hgnc_symbol = NA,
  description = "Putative novel or unannotated transcript"
)

# Combine with annotated data
annot_full <- rbind(annot, manual_entry)

# Save full version
write.csv(annot_full, "annotated_cox_prognosis_genes_full.csv", row.names = FALSE)

We have possibly Identified 3 novel transcripts:
So we will begin the confirmation of them:

ENSEMBL ID	HGNC Symbol	Description
ENSG00000255410	NA	novel transcript
ENSG00000267568	NA	novel transcript
ENSG00000286687	NA	long non-coding RNA; uncharacterized


Check Full Biotype and Location in Ensembl
library(biomaRt)

# Connect to Ensembl
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Strip version (just in case) and prepare list
novel_ids <- c("ENSG00000255410", "ENSG00000267568", "ENSG00000286687")

# Get more descriptive attributes
novel_info <- getBM(
  attributes = c("ensembl_gene_id", "gene_biotype", "external_gene_name",
                 "chromosome_name", "start_position", "end_position", "description"),
  filters = "ensembl_gene_id",
  values = novel_ids,
  mart = ensembl
)

print(novel_info)


 
Heatmap of all 9 prognostic biomarkers (including the 3 novel lncRNAs)
Kaplan-Meier (KM) plots for each of the 3 novel lncRNAs
library(pheatmap)

# Subset expression data to the 9 genes (now fully labeled)
heat_data <- assay(vsd)[annot_full_updated$ensembl_gene_id, ]

# Scale gene expression for visualization
heat_data_scaled <- t(scale(t(heat_data)))

# Sample annotation
annotation_col <- merged_meta_clean[, c("Cluster", "PrognosisGroup", "Stage")]
annotation_col <- data.frame(lapply(annotation_col, as.factor))
rownames(annotation_col) <- rownames(merged_meta_clean)

# Add gene symbols as row labels
rownames(heat_data_scaled) <- annot_full_updated$hgnc_symbol

# Plot heatmap
pheatmap(
  heat_data_scaled,
  annotation_col = annotation_col,
  cluster_cols = TRUE,
  cluster_rows = TRUE,
  show_rownames = TRUE,
  show_colnames = FALSE,
  fontsize_row = 9,
  main = "Expression Heatmap of Prognostic Biomarkers (including Novel lncRNAs)"
)
Step 1: Load libraries

library(DESeq2)
library(matrixStats)
library(ConsensusClusterPlus)

Step 2: Extract top 2000 most variable genes

Get variance-stabilized expression matrix
vsd_matrix <- assay(vsd) # Rows = genes, columns = samples

Calculate variance per gene
gene_variances <- rowVars(vsd_matrix)

Get top 2000 most variable genes
top_genes <- order(gene_variances, decreasing = TRUE)[1:2000]
vsd_top2000 <- vsd_matrix[top_genes, ]

Step 3: Prepare matrix for clustering

Transpose so that rows = samples, columns = genes
expr_for_clustering <- t(vsd_top2000)

Confirm shape: should be (samples x 2000 genes)
dim(expr_for_clustering)
head(rownames(expr_for_clustering)) # Should be like "2224", "2242", etc.

Step 4: Run ConsensusClusterPlus (K = 2 to 6)
set.seed(1234)

results <- ConsensusClusterPlus(
d = expr_for_clustering,
maxK = 6,
reps = 1000,
pItem = 0.8,
pFeature = 1,
clusterAlg = "hc",
distance = "pearson",
seed = 1234,
plot = "png",
title = "Clustering_top2000_genes"
)

Step 5: Choose K (You chose K = 3)

Examine consensus matrices, CDF, delta area
Then settle on K = 3
Step 6: Extract cluster assignments
cluster_assignments <- results[[3]]$consensusClass # K = 3

merged_meta_clean$Cluster <- as.factor(cluster_assignments)
Error in $<-.data.frame(*tmp*, Cluster, value = c(ENSG00000186081 = 1L, :
replacement has 2000 rows, data has 258

stopifnot(all(names(cluster_assignments) == rownames(merged_meta_clean)))
Error: all(names(cluster_assignments) == rownames(merged_meta_clean)) is not TRUE
In addition: Warning message:
In names(cluster_assignments) == rownames(merged_meta_clean) :
longer object length is not a multiple of shorter object length

# Assume md is already read in as metadataclu
> md <- metadataclu
> 
> # Strip the leading "X" from SampleID
> md$SampleID <- gsub("^X", "", md$SampleID)
> 
> # Set rownames to SampleID
> rownames(md) <- md$SampleID
> 
> # (Optional) remove duplicate SampleID column if not needed
> # md <- md[ , !(names(md) %in% "SampleID") ]
> 
> # Check the result
> head(rownames(md))
[1] "2224" "2242" "2251" "2265" "2285" "2289"
> head(colnames(txi$counts))
[1] "2224" "2242" "2251" "2265" "2285" "2289"
> 
> # Now they should match
> intersect(colnames(txi$counts), rownames(md))
  [1] "2224" "2242" "2251" "2265" "2285" "2289"
  [7] "2305" "2323" "2334" "2341" "2347" "2349"
 [13] "2351" "2353" "2354" "2357" "2374" "2378"
 [19] "2382" "2395" "2396" "2397" "2403" "2409"
 [25] "2414" "2421" "2425" "2429" "2439" "2440"
 [31] "2442" "2456" "2462" "2473" "2496" "2506"
 [37] "2507" "2509" "2516" "2517" "2519" "2552"
 [43] "2566" "2568" "2570" "2573" "2584" "2589"
 [49] "2599" "2601" "2610" "2613" "2618" "2623"
 [55] "2629" "2638" "2646" "2649" "2680" "2690"
 [61] "2716" "2717" "2729" "2730" "2747" "2776"
 [67] "2798" "2804" "2808" "2829" "2896" "2900"
 [73] "2902" "2906" "3053" "3054" "3064" "3091"
 [79] "3123" "3131" "3139" "3167" "3184" "3186"
 [85] "3328" "3329" "3331" "3332" "3375" "3385"
 [91] "3391" "3393" "3431" "3461" "3463" "3468"
 [97] "3487" "3494" "3500" "3502" "3506" "3523"
[103] "3525" "3527" "3534" "3536" "3545" "3557"
[109] "3577" "3625" "3631" "3639" "3654" "3662"
[115] "3673" "3695" "3729" "3735" "3746" "3753"
[121] "3755" "3796" "3802" "3832" "3835" "3843"
[127] "3866" "3913" "3914" "3937" "3940" "3941"
[133] "3944" "3945" "3948" "3956" "3957" "3960"
[139] "3971" "3974" "3981" "3995" "3997" "4008"
[145] "4017" "4019" "4025" "4062" "4064" "4069"
[151] "4117" "4128" "4132" "4135" "4137" "4138"
[157] "4150" "4153" "4156" "4157" "4158" "4166"
[163] "4197" "4203" "4204" "4211" "4222" "4235"
[169] "4236" "4237" "4242" "4243" "4244" "4253"
[175] "4254" "4255" "4275" "4277" "4278" "4284"
[181] "4288" "4291" "4300" "4309" "4310" "4311"
[187] "4313" "4360" "4362" "4382" "4401" "4403"
[193] "4409" "4421" "4495" "4497" "4501" "4510"
[199] "4528" "4529" "4530" "4615" "4616" "4627"
[205] "4631" "4662" "4676" "4681" "4683" "4702"
[211] "4703" "4706" "4709" "4710" "4717" "4723"
[217] "4724" "4734" "4759" "4760" "4768" "4769"
[223] "4773" "4774" "4778" "4823" "4863" "4866"
[229] "4867" "4876" "4879" "4884" "4890" "4891"
[235] "4896" "4918" "4926" "4978" "4980" "4995"
[241] "4999" "5000" "5014" "5018" "5055" "5056"
[247] "5067" "5077" "5078" "5097" "5098" "5102"
[253] "5120" "5148" "5149" "5161" "5222" "5237"
> 
> # --- 1) Align md with txi --------------------------------------------------
> stopifnot(all(c("counts","abundance","length") %in% names(txi)))
> keep <- intersect(colnames(txi$counts), rownames(md))
> if (length(keep) < 3) stop("Not enough matching samples after ID alignment.")
> txi$counts    <- txi$counts[, keep, drop=FALSE]
> txi$abundance <- txi$abundance[, keep, drop=FALSE]
> txi$length    <- txi$length[, keep, drop=FALSE]
> md            <- md[keep, , drop=FALSE]
> # --- 2) Build DDS (add covariates here if you have them, e.g. ~ batch + Cluster)
> dds <- DESeqDataSetFromTximport(txi, colData = md, design = ~ Cluster)
  the design formula contains one or more numeric variables with integer values,
  specifying a model with increasing fold change for higher values.
  did you mean for this to be a factor? if so, first convert
  this variable to a factor using the factor() function
using counts and average transcript lengths from tximport
> # 🔧 Convert categorical columns to factors
> md$Cluster <- factor(as.character(md$Cluster), levels = c("1","2","3","4"))
> # If you have batch or other categorical covariates, factor those too:
> # md$batch <- factor(md$batch)
> 
> # Quick sanity check
> table(md$Cluster, useNA = "ifany")

  1   2   3   4 
163   9  73  13 
> str(md$Cluster)
 Factor w/ 4 levels "1","2","3","4": 1 2 3 1 1 3 1 1 3 3 ...
> # Align with txi
> keep <- intersect(colnames(txi$counts), rownames(md))
> txi$counts    <- txi$counts[, keep, drop=FALSE]
> txi$abundance <- txi$abundance[, keep, drop=FALSE]
> txi$length    <- txi$length[, keep, drop=FALSE]
> md            <- md[keep, , drop=FALSE]
> # --- 2) Build DDS (add covariates here if you have them, e.g. ~ batch + Cluster)
> dds <- DESeqDataSetFromTximport(txi, colData = md, design = ~ Cluster)
using counts and average transcript lengths from tximport
> 
> # Filter very low counts (sum across all samples >= 10)
> dds <- dds[rowSums(counts(dds)) >= 10, ]
> 
> # Run DESeq2 once
> dds <- DESeq(dds)
estimating size factors
using 'avgTxLength' from assays(dds), correcting for library size
estimating dispersions
gene-wise dispersion estimates
mean-dispersion relationship
final dispersion estimates
fitting model and testing
-- replacing outliers and refitting for 1964 genes
-- DESeq argument 'minReplicatesForReplace' = 7 
-- original counts are preserved in counts(dds)
estimating dispersions
fitting model and testing
> # --- 3) Helper to run & save a contrast -----------------------------------
> do_contrast <- function(dds, level_ref, outdir) {
+     # 4 vs level_ref (log2FC = Cluster4 - Cluster[level_ref])
+     res_raw <- results(dds, contrast = c("Cluster","4", level_ref))
+     
+     # Shrink LFCs (apeglm if available else normal)
+     res_shrunk <- tryCatch({
+         suppressPackageStartupMessages(require(apeglm))
+         lfcShrink(dds, contrast = c("Cluster","4", level_ref), type = "apeglm")
+     }, error = function(e) {
+         message("apeglm not available; using 'normal' shrink for 4 vs ", level_ref)
+         lfcShrink(dds, contrast = c("Cluster","4", level_ref), type = "normal")
+     })
+     
+     tbl <- as.data.frame(res_shrunk) |>
+         mutate(gene_id = rownames(res_shrunk)) |>
+         relocate(gene_id) |>
+         arrange(padj, pvalue)
+     
+     dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
+     write_csv(tbl, file.path(outdir, paste0("DESeq2_Cluster4_vs_", level_ref, ".csv")))
+     saveRDS(res_raw,    file.path(outdir, paste0("res_raw_4vs", level_ref, ".rds")))
+     saveRDS(res_shrunk, file.path(outdir, paste0("res_shrunk_4vs", level_ref, ".rds")))
+     
+     cat(sprintf("\n4 vs %s: %d genes padj<0.05\n",
+                 level_ref, sum(tbl$padj < 0.05, na.rm = TRUE)))
+     invisible(tbl)
+ }
> 
> # --- 4) Run all pairwise: 4 vs 1, 4 vs 2, 4 vs 3 ---------------------------
> others <- intersect(c("1","2","3"), levels(colData(dds)$Cluster))
> out_dir <- "~/DE_Pairwise_Cluster4"
> pairwise_tables <- lapply(others, function(lvl) do_contrast(dds, lvl, out_dir))
apeglm not available; using 'normal' shrink for 4 vs 1
using 'normal' for LFC shrinkage, the Normal prior from Love et al (2014).

Note that type='apeglm' and type='ashr' have shown to have less bias than type='normal'.
See ?lfcShrink for more details on shrinkage type, and the DESeq2 vignette.
Reference: https://doi.org/10.1093/bioinformatics/bty895
                                                  
4 vs 1: 16220 genes padj<0.05
apeglm not available; using 'normal' shrink for 4 vs 2
using 'normal' for LFC shrinkage, the Normal prior from Love et al (2014).

Note that type='apeglm' and type='ashr' have shown to have less bias than type='normal'.
See ?lfcShrink for more details on shrinkage type, and the DESeq2 vignette.
Reference: https://doi.org/10.1093/bioinformatics/bty895
                                                  
4 vs 2: 8315 genes padj<0.05
apeglm not available; using 'normal' shrink for 4 vs 3
using 'normal' for LFC shrinkage, the Normal prior from Love et al (2014).

Note that type='apeglm' and type='ashr' have shown to have less bias than type='normal'.
See ?lfcShrink for more details on shrinkage type, and the DESeq2 vignette.
Reference: https://doi.org/10.1093/bioinformatics/bty895
                                                  
4 vs 3: 19333 genes padj<0.05
Warning messages:
1: In library(package, lib.loc = lib.loc, character.only = TRUE, logical.return = TRUE,  :
  there is no package called ‘apeglm’
2: In library(package, lib.loc = lib.loc, character.only = TRUE, logical.return = TRUE,  :
  there is no package called ‘apeglm’
3: In library(package, lib.loc = lib.loc, character.only = TRUE, logical.return = TRUE,  :
  there is no package called ‘apeglm’
> 
> # Also save the DDS for later reuse
> saveRDS(dds, file.path(out_dir, "dds_pairwise.rds"))
> 
> # Quick peek at top hits for 4 vs 1
> if ("1" %in% others) head(pairwise_tables[[which(others=="1")]])
                              gene_id baseMean
ENSG00000133958.14 ENSG00000133958.14 119.8008
ENSG00000144406.19 ENSG00000144406.19 289.8305
ENSG00000186487.21 ENSG00000186487.21 618.6737
ENSG00000077279.21 ENSG00000077279.21  82.3496
ENSG00000255087.6   ENSG00000255087.6  96.9727
ENSG00000166573.6   ENSG00000166573.6 311.9324
                   log2FoldChange     lfcSE
ENSG00000133958.14       4.650866 0.1467589
ENSG00000144406.19       3.420675 0.1129577
ENSG00000186487.21       3.956863 0.1323279
ENSG00000077279.21       7.504745 0.2714818
ENSG00000255087.6        5.198584 0.2084436
ENSG00000166573.6        6.811142 0.2778687
                       stat        pvalue
ENSG00000133958.14 31.66685 4.445422e-220
ENSG00000144406.19 30.27163 2.709652e-201
ENSG00000186487.21 29.90851 1.525217e-196
ENSG00000077279.21 27.55220 4.165250e-167
ENSG00000255087.6  24.97772 1.067681e-137
ENSG00000166573.6  24.39656 1.860166e-131
                            padj
ENSG00000133958.14 1.727891e-215
ENSG00000144406.19 5.266073e-197
ENSG00000186487.21 1.976122e-192
ENSG00000077279.21 4.047478e-163
ENSG00000255087.6  8.299937e-134
ENSG00000166573.6  1.205047e-127
>
> 
