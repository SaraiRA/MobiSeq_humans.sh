#! /bin/bash

## Load the modules required to run this pipeline on the high performance cluster
module load python/v2.7.12
module load perl/v5.24.0
module load R/v3.4.1
module load java/v1.8.0_131
module load mapDamage/v2.0.6
module load AdapterRemoval/v2.2.2
module load bwa/v0.7.15
module load htslib/v1.6
module load samtools/v1.6
module load paleomix/v1.2.12
module load preseq/v2.0
module load bedtools/v2.26.0
module load kentTools/v22032018
module load picard/v2.13.2
module load agplus/v1.0
module load cutadapt/v1.11
module load angsd/v0.921
module load ngsTools
module load fastme/v2.1.5
module load RAxML-NG/v0.5.1b
module load mapDamage/v2.0.9

## Define the variables required for the pipeline. Most of these are directory
PROJECT=/groups/hologenomics/sarai/data/human_mobiSeq
HUMANGENOME=/groups/hologenomics/data/genomes/human/hg38.UCSC.fasta
FASTQ=$PROJECT/fastqs
MAP=$PROJECT/mapping
ADAPTER=$PROJECT/adapRem
MATCH=$PROJECT/matches
MAPDAMAGE=$PROJECT/mapDamage


##################################### MATCHES #####################################
# Have a look of all the matches of the primer in the reference genome 
# The matches in the genome itself is not used downstream, since there are plots
# of other places in the genome where these primers are present, but these loci
# are not present in the reference genome. We use another technique to get all
# the possible matches later.

echo "Generating matches."
mkdir -p $MATCH
cd $MATCH

## For each primer, make afile with the primer sequrence.
if [ ! -e AluV1.txt ]; then echo "AluV1 CGATTACAGGCGTGAGCCACCGCGCC " > AluV1.txt; fi
if [ ! -e AluV2.txt ]; then echo "AluV2 CATGAGCCACCGCGCCCGGC" > AluV2.txt; fi

## For each primer, find the matches in the appropriate reference genome.
if [ ! -e human_AluV1.bed ]; then
  python $PROJECT/code/estNumberPrimers.py -f $HUMANGENOME -s AluV1.txt -o human_AluV1.txt -b human_AluV1.bed
fi
if [ ! -e human_AluV2.bed ]; then
  python $PROJECT/code/estNumberPrimers.py -f $HUMANGENOME -s AluV2.txt -o human_AluV2.txt -b human_AluV2.bed
fi
echo "Done generating matches."


################################## INFORMATION OF READS ##################################
# Obtain the number of reads for every fastq file
## NO run with the rest of the pipeline

##Count number total of reads in fastq.gz files 
#bc command is used for command line calculator
echo $(zcat PRI_HVMQ_P56_Aluv1__MobiSeq_lib_S1_L004_R1_001.fastq.gz|wc -l)/4 | bc

#loop, ls list the samples 
# If you have problems running the loop, do it in one line:  
# for file in $(ls *R1_001.fastq.gz); do echo $file; echo $(zcat $file |wc -l)/4 | bc; done 
for file in $(ls *R1_001.fastq.gz)
do 
	echo $file  
	echo $(zcat $file |wc -l)/4 | bc
done



##################################### ADAPTERS #####################################

# This script is for unpair, single data

####### Adapter Removal

## Remove the illumina adapters P5 and P7 (AdapterRemoval). 
# Yeah, here we are removing the adapters

#INPUT: fastq files

echo "Filter for TE target and remove adapters"
# -p: no error if existing 
mkdir -p $ADAPTER && cd $ADAPTER

## ALUV1
ALUV1ADAP=$ADAPTER/ALUV1
mkdir -p $ALUV1ADAP && cd $ALUV1ADAP
# Check if the adapter removal is done, if yes then this file must exist.
# If it does, skip this step, if not then perform cutadapt steps.
if [ ! -e .adap.done ]; then
	# For each fastq pair, do cutadapt
	for f in $FASTQ/Aluv1/*_L004_R1_001.fastq.gz
		do
		# Figure out the name of the sample.
		bn=$(basename $f _L004_R1_001.fastq.gz)
    		
		## Run Adapter removal upon the successful completion of cutadapt
		# collapse: Combined into a single consensus sequence pair aligments
		# Output: output_paired.collapsed containing merged reads,
		# output_paired.collapsed.truncated containing merged reads that have been trimmed due to the --trimns or --trimqualities options.
		# The sequencing primers are not specified, since AdapterRemoval knows the Illumina sequencing primers.
   		echo "AdapterRemoval --qualitybase 33 --qualitybase-output 33 --qualitymax 45 --gzip --mm 3 --minlength 25 --trimns --trimqualities --minquality 10 --basename ${bn}_noAdap --file1 $f"
  	done | xsbatch -c 1 --mem-per-cpu=2G -R -J ALUV1adap --
 	touch .adap.done
fi

## ALUV2
ALUV2ADAP=$ADAPTER/ALUV2
mkdir -p $ALUV2ADAP && cd $ALUV2ADAP
if [ ! -e .adap.done ]; then
	for f in $FASTQ/Aluv2/*_L004_R1_001.fastq.gz
  		do
    		bn=$(basename $f _L004_R1_001.fastq.gz)
    		echo "AdapterRemoval --qualitybase 33 --qualitybase-output 33 --qualitymax 45 --gzip --mm 3 --minlength 25 --trimns --trimqualities --minquality 10 --basename ${bn}_noAdap --file1 $f"
  		done | xsbatch -c 1 --mem-per-cpu=2G -R -J ALUV2adap --
  		touch .adap.done
fi


####### Keep the sequences with the TE target #################

## Verify that all the sequences have the TE target (cutadapt).
# Filter all the sequences that does not end in the 5' of R2 with the TE target, allowing 5% mismatches.
# IMPORTANT: We are not removing the primer sequence, just the sequences that does not have the primer ;)


## ALUV1
cd $ALUV1ADAP
# Check if the adapter removal is done, if yes then this file must exist.
# If it does, skip this step, if not then perform cutadapt steps.
if [ ! -e .cut.done ]; then
	# For each fastq pair, do cutadapt
	for f in $ALUV1ADAP/*_noAdap.truncated.gz;
		do
		# Figure out the name of the sample.
		bn=$(basename $f _noAdap.truncated.gz)
		## Run cutadapt to verify that all the sequences have the TE target. 
    		# -a : check for sequence on the 3' end of read1. It is N, so any base will do.
		# -G : check for sequence on the 5' end of read2. It is the primer for ALUV1, ^ means it should be in the star. 
		# --discard-untrimmed : discard sequences that dont have the adapters
  		# --no-trim: do not remove the primer sequences.
    		# --no-indels: do not account for short insertions and deletions.
    		# -e : mismatch rate, here it is 5%, 1 mismatch in the primer sequence.
    		# -o : output name for read 1
  		# -p : output name for read 2
    		# $f and ${f/R1/R2}: read 1 and read 2 input files respectively.
   		echo "mv $f ${bn}_noAdap.fastq.gz; cutadapt -a GGCGCGGTGGCTCACGCCTGTAATCG$ --discard-untrimmed --no-trim --no-indels -e 0.05 -o ${bn}_noAdap_primer.gz ${bn}_noAdap.fastq.gz "
  	done | xsbatch -c 1 --mem-per-cpu=2G -R -J ALUV1cut --
 	touch .cut.done
fi


## ALUV2
cd $ALUV2ADAP
# Check if the adapter removal is done, if yes then this file must exist.
# If it does, skip this step, if not then perform cutadapt steps.
if [ ! -e .cut.done ]; then
	# For each fastq pair, do cutadapt
	for f in $ALUV2ADAP/*_noAdap.truncated.gz;
		do
		# Figure out the name of the sample.
		bn=$(basename $f _noAdap.truncated.gz)
		## Run cutadapt to verify that all the sequences have the TE target. 
    		# -a : check for sequence on the 3' end of read1. It is N, so any base will do.
		# -G : check for sequence on the 5' end of read2. It is the primer for ALUV1, ^ means it should be in the star. 
		# --discard-untrimmed : discard sequences that dont have the adapters
  		# --no-trim: do not remove the primer sequences.
    		# --no-indels: do not account for short insertions and deletions.
    		# -e : mismatch rate, here it is 5%, 1 mismatch in the primer sequence.
    		# -o : output name for read 1
  		# -p : output name for read 2
    		# $f and ${f/R1/R2}: read 1 and read 2 input files respectively.
   		echo "mv $f ${bn}_noAdap.fastq.gz; cutadapt -a GCCGGGCGCGGTGGCTCATG$ --discard-untrimmed --no-trim --no-indels -e 0.05 -o ${bn}_noAdap_primer.gz ${bn}_noAdap.fastq.gz "
  	done | xsbatch -c 1 --mem-per-cpu=2G -R -J ALUV2cut --
 	touch .cut.done
fi



################################## MAPPING AND MARK DUPLICATES ##########################################
# CAREFUL: Names were changed
# Map to the appropriate reference genomes using BWA mem. 
# Use mem so that we can soft clip the primer sequences in the beginning of the read.

# Sort the mapping by coordinates and mark duplicates using samtools 
echo "Map to reference genome"
# -p: no error if existing 
mkdir -p $MAP && cd $MAP

#INPUT: *_primer_noAdap.collapsed.gz and  *_primer_noAdap.collapsed.truncated.gz

## ALUV1
ALUV1BAM=$MAP/ALUV1
mkdir -p $ALUV1BAM
cd $ALUV1BAM
if [ ! -e .map.done ]; then
	## Run BWA mem to map
	# Discard any reads that are orphans of a pair, so keep only valid pairs
  	for f in $ADAPTER/ALUV1/*_noAdap_primer.fq;
  		do
    		bn=$(basename $f _noAdap_primer.fq)
    		# Run bwa mem and then sort then sort the bam by coordinates.
		# M: mark shorter split hits as secondary
    		echo "(bwa aln $HUMANGENOME $ADAPTER/ALUV1/${bn}_noAdap_primer.fq > ${bn}_noAdap_primer_hg38.UCSC.sai; bwa samse -f ${bn}_noAdap_primer_hg38.UCSC.sam $HUMANGENOME ${bn}_noAdap_primer_hg38.UCSC.sai $ADAPTER/ALUV1/${bn}_noAdap_primer.fq; samtools view -u ${bn}_noAdap_primer_hg38.UCSC.sam | samtools fixmate -r -p -m - - | samtools sort - | samtools markdup -r - ${bn}_noAdap_primer_hg38.UCSC.markdup.bam  )"
  	done | xsbatch -c 1 --mem-per-cpu=10G -J aluv1 -R --max-array-jobs=10 --
  touch .map.done
fi



## ALUV2
ALUV2BAM=$MAP/ALUV2
mkdir -p $ALUV2BAM
cd $ALUV2BAM
if [ ! -e .map.done ]; then
	## Run BWA mem to map
	# Discard any reads that are orphans of a pair, so keep only valid pairs
  	for f in $ADAPTER/ALUV2/*_noAdap_primer.fq;
  		do
    		bn=$(basename $f _noAdap_primer.fq)
    		# Run bwa mem and then sort then sort the bam by coordinates.
		# M: mark shorter split hits as secondary
    		echo "(bwa aln $HUMANGENOME $ADAPTER/ALUV2/${bn}_noAdap_primer.fq > ${bn}_noAdap_primer_hg38.UCSC.sai; bwa samse -f ${bn}_noAdap_primer_hg38.UCSC.sam $HUMANGENOME ${bn}_noAdap_primer_hg38.UCSC.sai $ADAPTER/ALUV2/${bn}_noAdap_primer.fq; samtools view -u ${bn}_noAdap_primer_hg38.UCSC.sam | samtools fixmate -r -p -m - - | samtools sort - | samtools markdup -r - ${bn}_noAdap_primer_hg38.UCSC.markdup.bam )"
  	done | xsbatch -c 1 --mem-per-cpu=10G -J aluv2 -R --max-array-jobs=10 --
  touch .map.done
fi



################################## MAP DAMAGE ##########################################

mapDamage -i $bam -r $reference  -rescale

echo "Running map Damage"
# -p: no error if existing 
mkdir -p $MAPDAMAGE && cd $MAPDAMAGE

## ALUV1
ALUV1DAM=$MAPDAMAGE/ALUV1
mkdir -p $ALUV1DAM && cd $ALUV1DAM
# Check if the adapter removal is done, if yes then this file must exist.
# If it does, skip this step, if not then perform cutadapt steps.
if [ ! -e .dam.done ]; then
	# For each fastq pair, do cutadapt
	for f in $MAP/ALUV1/*.markdup.bam; 
		do
		# Figure out the name of the sample.
		bn=$(basename $f .markdup.bam )
    		
		## Run Map Damage
		# parameter can be optionally used to rescale quality scores of likely damaged positions in the reads. A new BAM file is constructed by downscaling quality values for misincorporations likely due to ancient DNA damage according to their initial qualities, position in reads and damage patterns.
   		echo "mapDamage -i $f -r $HUMANGENOME --rescale "
  	done | xsbatch -c 1 --mem-per-cpu=2G -R -J ALUV1dam --
 	touch .dam.done
fi


## ALUV2
ALUV2DAM=$MAPDAMAGE/ALUV2
mkdir -p $ALUV2DAM && cd $ALUV2DAM
# Check if the adapter removal is done, if yes then this file must exist.
# If it does, skip this step, if not then perform cutadapt steps.
if [ ! -e .dam.done ]; then
	# For each fastq pair, do cutadapt
	for f in $MAP/ALUV2/*.markdup.bam 
		do
		# Figure out the name of the sample.
		bn=$(basename $f .markdup.bam )
    		
		## Run Map Damage
		# parameter can be optionally used to rescale quality scores of likely damaged positions in the reads. A new BAM file is constructed by downscaling quality values for misincorporations likely due to ancient DNA damage according to their initial qualities, position in reads and damage patterns.
   		echo "mapDamage -i $f -r $HUMANGENOME --rescale "
  	done | xsbatch -c 1 --mem-per-cpu=2G -R -J ALUV2dam --
 	touch .dam.done
fi


xsbatch -c 1 --mem-per-cpu=2G -J paleomix --paleomix bam_pipeline run paloemix_P56.yaml


##### WE ARE HERE
