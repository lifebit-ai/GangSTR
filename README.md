# GangSTR

NB: Many of the inputs/flags in the original GangSTR documentation will not work! Instead use the paramters used in the example below. If any paramters from the original GangSTR documentation (which is included below the example) are required they can be implemented by modifying the main.nf script. For any queries about implementing more paramters (eg --genomewide) you can email me at: phil@lifebit.ai.

## Example command that can be run over Deploit
Example/test command that can be run on [Deploit](https://deploit.lifebit.ai/) using publically available test data. The data can be imported from the AWS S3 bucket [s3://lifebit-featured-datasets/](https://s3.console.aws.amazon.com/s3/buckets/lifebit-featured-datasets/pipelines/RepeatExpansion/?region=eu-west-1&tab=overview) 
```        
nextflow run lifebit-ai/GangSTR --genome RepeatExpansion/RepeatExpansion/Reference/hs37d5.fa
                                --bam RepeatExpansion/RepeatExpansion/Bams/HG00472.mapped.ILLUMINA.bwa.CHS.exome.20121211.bam
                                --regions RepeatExpansion/RepeatExpansion/GangSTR/hs37_ver8.bed 
                                --nonuniform 
                                --ref RepeatExpansion/RepeatExpansion/ExpansionHunter/repeat-specs/grch37
```

* Please ensure the index file (.bam.bai) is located within the same folder. For example, for the test data the index file HG00472.mapped.ILLUMINA.bwa.CHS.exome.20121211.bam.bai is already located within RepeatExpansion/RepeatExpansion/Bams folder and so nothing more needs to be done.
* **--nonuniform** was used due to non-uniform coverage in alignment file in the test data. This option will likely not be needed for larger bam files espeiallly full-genome sequences
*  **Warning** The execution time was 1hr 30mins on a m2.2xlarge (medium cost saving spot instance)
* The output should be found in the format of output.vcf


<br />
<br />
<br />
<br />
<br />



## Original GangSTR Documentation

GangSTR is a tool for genome-wide profiling tandem repeats from short reads. A key advantage of GangSTR over existing tools (e.g. [lobSTR](https://github.com/mgymrek/lobstr-code) or [hipSTR](https://github.com/tfwillems/HipSTR)) is that it can handle repeats that are longer than the read length.

GangSTR takes aligned reads (BAM) and a set of repeats in the reference genome as input and outputs a VCF file containing genotypes for each locus.

**bioRxiv preprint** manuscript: http://bit.ly/gangstr-preprint

For questions on installation or usage, please open an issue, submit a pull request, or contact Nima Mousavi (mousavi@ucsd.edu).

[Download](#download) | [Install](#install) | [Usage](#usage) | [File formats](#formats) | [References](#references)

<a name="download"></a>
## Download

The latest GangSTR release is available on the [releases page](https://github.com/gymreklab/GangSTR/releases).

For a list of TR references available, see [references](#references) below. 

<a name="install"></a>

## Basic Install

GangSTR requires third party packages [nlopt](https://nlopt.readthedocs.io/en/latest/), [gsl](https://www.gnu.org/software/gsl/doc/html/index.html), and [htslib](http://www.htslib.org//). The built-in script `install-gangstr.sh` installs these for you before compiling and installing GangSTR. Both UNIX and Mac OSX are supported.

If you are running as root:
```
tar -xzvf GangSTR-X.X.tar.gz
cd GangSTR-X.X
sudo ./install-gangstr.sh
```

If you are installing locally (e.g. on a cluster where you don't have root access):
```
tar -xzvf GangSTR-X.X.tar.gz
cd GangSTR-X.X
./install-gangstr.sh PREFIX
```


where `PREFIX` is a place you have write permissions. In most cases this will be your home directory, e.g. `$HOME`. If you install locally, make sure `$PREFIX/bin` is on your `PATH`.


Typing `GangSTR --help` should show a help message if GangSTR was successfully installed.

<a name="usage"></a>
## Usage
To run GangSTR using default parameters use the following command:
```
GangSTR --bam file.bam 
        --ref ref.fa 
        --regions regions.bed 
        --out outprefix 
```
Required parameters:
* **--bam** Alignment file (.bam)
* **--ref** Refererence genome (.fa)
* **--regions** Target TR loci (.bed)
* **--out** Output prefix

Additional general options:
* **--genomewide** Run GangSTR in genome-wide mode. This mode has more stringent filtering steps to prevent false positive in genome-wide profiling.

Options for different sequencing settings
* **--readlength \<int\>** Preset read length (default: extract from alignments if not provided)
* **--coverage \<float\>** Preset average coverage, should be set for targeted data (default: calculate if not provided)
* **--nonuniform** Indicates non-uniform coverage in alignment file (i.e., used for exome sequencing). Using this flag removes the likelihood term corresponding to FRR count.

Advanced parameters for likelihood model:
* **--frrweight \<float\>** Reset weight for FRR class in likelihood model (default 0.5)
* **--spanweight \<float\>** Reset weight for Spanning class in likelihood model (default 1.0)
* **--enclweight \<float\>** Reset weight for Enclosing class in likelihood model (default 1.0)
* **--flankweight \<float\>** Reset weight for Flanking class in likelihood model (default 1.0)
* **--ploidy [0,1]** Haploid (1) or diploid (2) genotyping (default 2)
* **--useofftarget** Extract off-target FRRs based on the off-target regions provided in the regions file.
* **--insertmean \<float\>** Fragment length mean (default: calculate if not provided)
* **--insertsdev \<float\>** Fragment length standard deviation (default: calculate if not provided)
* **--insertmax \<float\>** Maximum allowed fragment length (default: no filtering based on fragment length)
* **--readprobmode** Only use read probabilities in likelihood model (ignore class probability)
* **--numbstrap \<int\>** Number of bootstrap samples for calculating confidence intervals (default 100)

Parameters for local realignment:
* **--minscore \<int\>** Minimun alignment score for accepting reads (default 75)
* **--minmatch \<int\>** Minimum matching basepairs required at the edge of the repeat region to accept flanking and enclosing reads (default 5)

Stutter model parameters:
* **--stutterup \<float\>** Stutter insertion probability (default 0.0364653)
* **--stutterdown \<float\>**	Stutter deletion probability (default: 0.0428387)
* **--stutterprob \<float\>**	Stutter step size parameter (default: 0.818913)

Parameters for more detailed info about each locus:
* **--output-readinfo** Output a file containing extracted read information
* **--output-bootstraps** Output a file containing bootstrap samples

Additional optional parameters:
* **-h,--help** display help screen
* **--seed** Random number generator initial seed
* **-v,--verbose** Print progress information (major steps)
* **--very** Print detailed progress information
* **--version** Print out the version of this software

<a name="formats"></a>
## File formats

GangSTR takes as input a BAM file of short read alignments, a reference set of TRs, and a reference genome, and outputs genotypes in a VCF file. Each of these formats is described below.

### BAM (`--bam`)
GangSTR requires a [BAM](https://samtools.github.io/hts-specs/SAMv1.pdf) file produced by an indel-sensitive aligner. The BAM file must be sorted and indexed e.g. by using `samtools sort` and `samtools index`. GangSTR currently only processes a single sample at a time.

### FASTA Reference genome (`--ref`)
You must input a reference genome in FASTA format. This must be the same reference build used to align the sequences in the BAM file.

### TR regions (`--regions`)
GangSTR requires a reference set of regions to genotype. This is a BED-like file with the following columns:

1. The name of the chromosome on which the STR is located
2. The start position of the STR on its chromosome
3. The end position of the STR on its chromosome
4. The motif length
5. The repeat motif

An optional 6th column may contain a comma-separated list of off-target regions for each TR. These are regions where misaligned reads for a given TR may be incorrectly mapped.

Below is an example file which contains 5 TR loci. Standard references for hg19 and GRCh38 can be obtained [below](#references).
**NOTE: The table header is for descriptive purposes. The BED file should not have a header**

| **CHROM** | **START** | **END** | **MOTIF_LEN** | **MOTIF** | **OFFTARGET (optional)** |
|-----------|-----------|---------|----------------|----------|------------|
|chr1	|10689	|10700|	5	|CGCGC|	|
| chr1  |  28589  | 28603  | 1 |      T    |   |
|chr4  |  11173|   11194  | 11   |   CGCCGGCGCGG |    |
|chr4   | 150889 | 150909 | 2    |   TG      ||
| chr19 | 45770205 | 45770264	| 3	| CAG	|chr2:163338502-163338506,chr3:197333949-197333955,chr6:16327632-16327646,chr6:170561926-170561931,chr7:122288209-122288215,chr8:133055822-133055827,chr11:28310883-28310888,chr17:4887671-4887677,chr18:55586148-55586165,chr19:13207866-13207871 |

### VCF (output)
For more information on VCF file format, see the [VCF spec](http://samtools.github.io/hts-specs/VCFv4.2.pdf). In addition to standard VCF fields, GangSTR adds custom fields described below.

#### INFO fields

INFO fields contain aggregated statistics about each TR. The following custom fields are added:

| **FIELD** | **DESCRIPTION** |
|-----------|------------------|
| END | End position of the TR |
| RU| Repeat motif | 
| REF| Reference copy number (number of repeat units| 

#### FORMAT fields
FORMAT fields contain information specific to each genotype call. The following custom fields are added:

| **FIELD** | **DESCRIPTION** |
|-----------|------------------|
| GB | Base pair length differences of genotype from reference for each allele |
| CI| 95% confidence intervals for each allele | 
| RC| Number of reads in each class (enclosing, spanning, FRR, flanking)| 
| Q| Minimum negative likelihood| 
| INS| Insert size mean and stddev at the locus| 


<a name="references"></a>
## GangSTR reference files

The following lists available references created using Tandem Repeats Finder. We update the reference periodically with additional loci or annotation changes. Please see the Changelog file for details. Note references must be unzipped before using with GangSTR. 

| **Reference build** | **Version** | **Link** |
| --------------------| ------------|----------|
| hg19 | ver8 | [hg19_ver8.bed.gz](https://s3.amazonaws.com/gangstr/hg19_ver8.bed.gz) | 
| hs37 | ver8 | [hs37_ver8.bed.gz](https://s3.amazonaws.com/gangstr/hs37_ver8.bed.gz) |
| hg38 | ver5 | [hg38_ver5.bed.gz](https://s3.amazonaws.com/gangstr/hg38_ver5.bed.gz) |

The references below contain pre-defined off-target loci for target pathogenic loci (hg38 coordinates):

| **Locus** | **Link** |
| ------| ------|
| SCA1 | [SCA1_hg38.bed](https://s3.amazonaws.com/gangstr/SCA1_hg38.bed) |
| SCA2 | [SCA2_hg38.bed](https://s3.amazonaws.com/gangstr/SCA2_hg38.bed) |
| SCA3 | [SCA3_hg38.bed](https://s3.amazonaws.com/gangstr/SCA3_hg38.bed) |
| SCA6 | [SCA6_hg38.bed](https://s3.amazonaws.com/gangstr/SCA6_hg38.bed) |
| SCA7 | [SCA7_hg38.bed](https://s3.amazonaws.com/gangstr/SCA7_hg38.bed) |
| SCA8 | [SCA8_hg38.bed](https://s3.amazonaws.com/gangstr/SCA8_hg38.bed) |
| SCA12 | [SCA12_hg38.bed](https://s3.amazonaws.com/gangstr/SCA12_hg38.bed) |
| SCA17 | [SCA17_hg38.bed](https://s3.amazonaws.com/gangstr/SCA17_hg38.bed) |
| HTT | [HTT_hg38.bed](https://s3.amazonaws.com/gangstr/HTT_hg38.bed) |
| DM1 | [DM1_hg38.bed](https://s3.amazonaws.com/gangstr/DM1_hg38.bed) |


