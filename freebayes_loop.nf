nextflow.enable.dsl=2

// Freebayes-Loop v1.0
// Region-based parallel variant calling with freebayes for pooled HPV samples

params.input_file = 'start_range.bed'
params.bam_path = null
params.publish_dir = "./results"
params.min_del = 10
params.ref = null
params.cpu = 4
params.mem = 40
params.help = false

def helpMessage() {
    log.info """
    =========================================
     Freebayes-Loop v1.0
     Region-based parallel variant calling
    =========================================

    Usage:
      nextflow run freebayes_loop.nf --bam_path <bam_list> --ref <reference.fa> [options]

    Required:
      --bam_path        Path to text file listing BAM file paths (one per line)
      --ref             Path to reference FASTA file (must be indexed)

    Optional:
      --input_file      BED file with genomic regions (default: start_range.bed)
      --publish_dir     Output directory (default: ./results)
      --min_del         Minimum deletion length in bp to report (default: 10)
      --cpu             Number of CPUs per process (default: 4)
      --mem             Memory in GB per process (default: 40)
      --help            Show this help message

    Output:
      merge.vcf.gz      Merged VCF with variants from all regions
      merge.vcf.gz.bed  BED file of detected large deletions
      vcf.log           Log of merge operation
    """.stripIndent()
}

if (params.help) {
    helpMessage()
    exit 0
}

if (!params.bam_path) {
    log.error "ERROR: --bam_path is required. Provide a file listing BAM paths."
    helpMessage()
    exit 1
}

if (!params.ref) {
    log.error "ERROR: --ref is required. Provide a reference FASTA file."
    helpMessage()
    exit 1
}

process freebayes {
    cpus params.cpu
    memory "${params.mem} GB"
    publishDir "${params.publish_dir}", mode: 'copy', overwrite: true

    input:
        tuple val(chr), val(start), val(end)

    output:
        path("variant_${start}-${end}.vcf.gz") , emit: vcf
        path("variant_${start}-${end}.vcf.gz.csi") , emit: index

    """
    freebayes \
        --fasta-reference $params.ref \
        --region '$chr':$start-$end \
        --ploidy 4 --pvar 0 --min-alternate-fraction 0.03 --min-alternate-count 1 -m 55 -q 13 --read-indel-limit 2 --use-best-n-alleles 4 \
        --min-coverage 200 --pooled-discrete --pooled-continuous --haplotype-length 0  \
        --report-monomorphic  \
        --bam-list $params.bam_path | bgzip -c > variant_${start}-${end}.vcf.gz

    bcftools index -f variant_${start}-${end}.vcf.gz
    """
}

process merge_vcfs {
    cpus params.cpu
    memory "${params.mem} GB"
    publishDir "${params.publish_dir}", mode: 'copy', overwrite: true

    input:
        path vcf_gz
        path vcf_index

    output:
        path "merge.vcf.gz"
        path "merge.vcf.gz.bed" , emit: del_bed
        path "vcf.log"

    script:
    """
    count=\$(ls $vcf_gz | wc -l)
    echo "Number of values: \$count"

    if [[ \$count -gt 1 ]]; then
        echo "Multiple input files detected = \$count :"$vcf_gz >> vcf.log
        bcftools concat $vcf_gz --allow-overlaps --output merge.vcf.gz
        python ${projectDir}/scripts/vcf2bed_03.py $params.min_del merge.vcf.gz 30 > merge.vcf.gz.bed
    else
        echo "Single input files detected = \$count :"$vcf_gz >> vcf.log
        mv $vcf_gz merge.vcf.gz
        python ${projectDir}/scripts/vcf2bed_03.py $params.min_del merge.vcf.gz 30 > merge.vcf.gz.bed
    fi
    """
}

workflow variant_calling {
    take:
        bed_input

    main:
        bin = bed_input
            | splitText { it.trim() }
            | splitCsv(header: ['chr', 'start', 'end'], sep: ' ' )
            | map { row-> tuple(row.chr, row.start, row.end) }

        freebayes(bin)
        merge_vcfs(freebayes.out.vcf.collect(), freebayes.out.index.collect())

    emit:
        merge_vcfs.out.del_bed
}

workflow {
    variant_calling(Channel.fromPath(params.input_file))
}
