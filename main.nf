params.genome = "s3://repeat-expansion/Reference/hs37d5.fa"
params.bam = "s3://repeat-expansion/Bams/HG00457.mapped.ILLUMINA.bwa.CHS.exome.20121211.bam"
params.ref= "s3://repeat-expansion/ExpansionHunter/repeat-specs/grch37"
params.regions = "s3://repeat-expansion/GangSTR/hs37_ver8.bed"

genome_file = file(params.genome)
genome_index = file(params.genome+".fai")
bam_file = file(params.bam)
bai_file = file(params.bam+".bai")
ref_dir = file(params.ref)
regions_file = file(params.regions)

process gangstr {
	publishDir 'results'

	input:
	file('aln.bam') from bam_file
	file('aln.bam.bai') from bai_file
	file('genome.fa') from genome_file
	file('genome.fa.fai') from genome_index
	file('ref') from ref_dir
  	file('regions.bed') from regions_file

	output:
	file('output.*') into results
	
	script:
	"""	
  	GangSTR \
  	--bam aln.bam \
  	--ref-fasta genome.fa \
  	--regions regions.bed \
  	--out output
	"""
}

workflow.onComplete {
	println ( workflow.success ? "\nGangSTR is done!" : "Oops .. something went wrong" )
}
    
