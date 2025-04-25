nextflow.enable.dsl=2

params.samplesheet = "samplesheet.csv"
params.genomeDir   = "reference/star_index"
params.gtf         = "reference/gencode.v44.annotation.gtf"
params.outdir      = "results"

workflow {
  samples = Channel
    .fromPath(params.samplesheet)
    .splitCsv(header: true)
    .map { row -> tuple(row.sample_id, row.condition, file(row.fastq)) }

  samples | fastqc | star

  star.out.collect() | featurecounts
}

// ────── FASTQC ──────
process fastqc {
  tag "$sample_id"
  conda "envs/rnaseq.yml"

  input:
  tuple val(sample_id), val(condition), path(fastq)

  output:
  path "results/fastqc/${sample_id}_fastqc.zip" into fastqc_out

  script:
  '''
  mkdir -p results/fastqc
  fastqc $fastq -o results/fastqc
  '''
}

// ────── STAR ALIGNMENT ──────
process star {
  tag "$sample_id"
  conda "envs/rnaseq.yml"

  input:
  tuple val(sample_id), val(condition), path(fastq)

  output:
  tuple val(sample_id), path("results/star/${sample_id}_Aligned.sortedByCoord.out.bam") into star_out

  script:
  '''
  mkdir -p results/star
  STAR --runThreadN 4 \
       --genomeDir ${params.genomeDir} \
       --readFilesIn $fastq \
       --outFileNamePrefix results/star/${sample_id}_ \
       --outSAMtype BAM SortedByCoordinate
  '''
}

// ────── FEATURECOUNTS ──────
process featurecounts {
  tag "featureCounts"
  conda "envs/rnaseq.yml"

  input:
  set sample_id_bam_pairs from star_out.collect()

  output:
  path "results/counts/gene_counts.txt"

  script:
  def bams = sample_id_bam_pairs.collect{ it[1] }.join(' ')
  '''
  mkdir -p results/counts
  featureCounts -T 4 \
    -a ${params.gtf} \
    -o results/counts/gene_counts.txt \
    $bams
  '''
}

