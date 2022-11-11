#!/usr/bin/env nextflow

/*
 * SET UP CONFIGURATION VARIABLES
 */
bam = Channel
		.fromPath(params.bam)
		.ifEmpty { exit 1, "${params.bam} not found.\nPlease specify --bam option (--bam bamfile)"}

Channel
		.fromPath(params.genome)
		.ifEmpty { exit 1, "${params.genome} not found.\nPlease specify --genome option (--genome fastafile)"}
		.into { fastaToFai; fastaToGangSTR }

if(params.fai){
	Channel
			.fromPath(params.fai)
			.ifEmpty { exit 1, "${params.fai} not found.\nMake sure your the file exists (--fai faifile) or remove the fai option for it be automatically generated"}
}

bed = Channel
    .fromPath(params.bed)
    .ifEmpty { exit 1, "${params.bed} not found.\nPlease specify --bed option (--bed bedfile)"}

if(params.nonuniform) {
	extraflags = '--nonuniform'
} else {
	extraflags = ""
}

// Header log info
log.info """=======================================================
		GangSTR
======================================================="""
def summary = [:]
summary['Pipeline Name']    = 'GangSTR'
summary['Bam file']         = params.bam
summary['Bed file']         = params.bed
summary['Reference genome'] = params.genome
summary['Output dir']       = params.outdir
summary['Working dir']      = workflow.workDir
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "========================================="


process preprocess_bam{

  tag "${bam}"
	container 'broadinstitute/gatk:4.3.0.0'

  input:
  file bam from bam

  output:
  set file("ready/${bam}"), file("ready/${bam}.bai") into completeChannel

  script:
  """
	echo '$task.container'
  mkdir ready
  [[ `samtools view -H ${bam} | grep '@RG' | wc -l`   > 0 ]] && { mv $bam ready;}|| { gatk  AddOrReplaceReadGroups \
	-I ${bam[0]} -O ready/${bam[0]} --RGID ${params.rgid} --RGLB ${params.rglb} --RGPL ${params.rgpl} --RGPU ${params.rgpu} --RGSM ${params.rgsm};}
  cd ready ;samtools index ${bam};
  """
}

if(!params.fai) {
  process preprocess_genome {

			tag "${fasta}"
			container 'lifebitai/preprocessingvctools'

      input:
      file fasta from fastaToFai

      output:
      file("${fasta}.fai") into fai

      script:
      """
      samtools faidx $fasta
      """
  }
}

process gangstr {
	publishDir "${params.outdir}", mode: 'copy'

	container = 'lifebitai/gangstr'

	input:
	set file(bam), file(bai) from completeChannel
	file fasta from fastaToGangSTR
	file fai from fai
	file bed from bed

	output:
	file('output*') into results

	script:
	"""
  	GangSTR \
  	--bam ${bam} \
  	--ref ${fasta} \
  	--regions ${bed} \
  	--out output ${extraflags}
	"""
}

workflow.onComplete {
	println ( workflow.success ? "\nGangSTR is done!" : "Oops .. something went wrong" )
}
