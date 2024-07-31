*varAmpliCNV*

varAmpliCNV is both an online/off-line tool for CNV calling from amplicon based targeted resequencing (TR) data enriched via Haloplex technologies. The highlighting features of the tool are:

  - Assigns reads directly to the amplicon coordinates to prune out aspecific alignments
  - Utilizes non-parametric model based on PCA/MDS to control the variance present in the data
  - Utilizes position specific dependencies between amplicons to filter out false positive CNV segments
  - Provides visualization plots to select CNV segments



### Requirements
 - R (>3.2.1 and <3.5.0 ), python 2.7
 - R packages: DNAcopy; ggplot2, Matrix, gtools
 - bedtools v2.28.0
 - All the BAM files
 - Amplicon design file (tab separated)
 - SNP-Amplicon-Gene File (tab separated)
 - Region of Interest (ROI) file annotated with exon number, gene name (tab separated)
 - Gender information per sample for every batches


### Stage 1: Preprocessing step 
1. Processing input Amplicon design file to remove SNP coordinates and duplicate coordinates. Also adding Region number to input region of interest (ROI) design file. 
```sh
$ python processAmpliconRegion.py -a  <RAW-amplicon-bed-file> -s <snp-amplicon-file> -r <ROI-design-file> -b <path-to-bedtools> -o <path-to-output-directory> 
 *Example:
$ python processAmpliconRegion.py -a /home/shared_data_medgen_aorta/VCNV/TAAD/BED/37328-1448381652_Amplicons.bed -s /home/shared_data_medgen_aorta/VCNV/TAAD/BED/snpAmpGene.bed -o /home/shared_data_medgen_aorta/VCNV/TAAD/BED/ -b /opt/NGS/binaries/BedTools/2.28.0/bin/bedtools -r /home/shared_data_medgen_aorta/VCNV/TAAD/BED/37328-1448381652_Regions.bed
```
2. Getting GC content per amplicon
	* Should be pre-computed in the variant calling step itself inorder to reduce the computational time.
```sh
$ python getGCAmplicon.py -a <path-to-amplicon-bed-file> -o <path-to-GCContent-output-file> -b <two-bit-binary> -f <two-bit-fasta-file>
 * Example:
$ python getGCAmplicon.py -a /home/shared_data_medgen_aorta/VCNV/TAAD/BED/AmpRmSNPRmDup.bed -o /home/shared_data_medgen_aorta/VCNV/TAAD/BED/GCContent.csv -b /opt/NGS/binaries/twoBitToFa/default/twoBitToFa -f /opt/NGS/References/hg19/2bit/hg19.2bit 
```
### Stage 2: Assigning reads uniquely to Amplicons
1.  The BAM files for a given batch are parsed to retrieve reads and assign uniquely to the amplicons based on exact matching to their start and end coordinates
* For every batch inside DATA folder for every amplicon and corresponding read count a file is created with sample name
```sh
$ python processBAM.py -i <input-BAM-directory> -a <amplicon-bed-file> -o <out-directory-path> -b <analysis-batch-name>
```
Example command for batch BATCH\_1:
```sh
$/opt/software/python2.7.5/bin/python processBAM.py -i /home/shared_data_medgen_aorta/TAAD/53 -a /home/shared_data_medgen_aorta/VCNV/TAAD/BED/AmpRmSNPRmDup.bed -o /home/shared_data_medgen_aorta/VCNV/TAAD -b BATCH_1
```  

### Stage 3: Predicting the CNV per analysis batch
Requirements
- Batch specific amplicon-by-sample Read Count matrix
- Unique, sorted amplicon design file
- GC content file obtained for given amplicon design file
- Gender details per sample
```sh
$ Rscript varAmpliCNV.R -i <sample-amplicon-count-Matrix> -o <output-directory> -b <amplicon-bed-file> -c <gc-content-file> -r <ROI-file> -s <gender-file> -a <batch Name/Number> -p <proprtion-of-variance for Auto/Sex: default:0.80> -n < DS=>0 or AOF=>1 Flag; default:0>
```
Example for processing BATCH\_1 of TAAD samples with direct segmentation or DS flag = 0:
```sh
$ /opt/software/R/R-3.2.1/bin/Rscript varAmpliCNV.R -i /home/shared_data_medgen_aorta/VCNV/TAAD/BATCH_1/MAT_RC/AmpCountMat.RData -o /home/shared_data_medgen_aorta/VCNV/TAAD/ -b /home/shared_data_medgen_aorta/VCNV/TAAD/BED/AmpRmSNPRmDup.bed -r /home/shared_data_medgen_aorta/VCNV/TAAD/BED/sortMergeROIGene.bed -c /home/shared_data_medgen_aorta/VCNV/TAAD/BED/GCContent.csv -s /home/shared_data_medgen_aorta/VCNV/TAAD/ALL_RUN_Sample_Gender.txt -a BATCH_1 -p 0.80 -n 0
```
Example for processing BATCH\_1 of TAAD samples using amplicon overlap filtering (DS-AOF) approach or DS flag = 1:
```sh
$ /opt/software/R/R-3.2.1/bin/Rscript varAmpliCNV.R -i /home/shared_data_medgen_aorta/VCNV/TAAD/BATCH_1/MAT_RC/AmpCountMat.RData -o /home/shared_data_medgen_aorta/VCNV/TAAD/ -b /home/shared_data_medgen_aorta/VCNV/TAAD/BED/AmpRmSNPRmDup.bed -r /home/shared_data_medgen_aorta/VCNV/TAAD/BED/sortMergeROIGene.bed -c /home/shared_data_medgen_aorta/VCNV/TAAD/BED/GCContent.csv -s /home/shared_data_medgen_aorta/VCNV/TAAD/ALL_RUN_Sample_Gender.txt -a BATCH_1 -p 0.80 -n 1
```

### Things To Do
- Add little more documentation to the code/functions
- Documentation explaining input file format types related to Amplicon design files.


