include { SAMTOOLS_COVERAGE } from '../../../modules/local/samtools/coverage/main'
include { SAMTOOLS_INDEX    } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_VIEW     } from '../../../modules/nf-core/samtools/view/main'
include { BCFTOOLS_MPILEUP  } from '../../../modules/nf-core/bcftools/mpileup/main'
include { BCFTOOLS_INDEX    } from '../../../modules/nf-core/bcftools/index/main'
include { BCFTOOLS_CALL     } from '../../../modules/nf-core/bcftools/call/main'

workflow BAM_SIMULATE {

    take:
    ch_bam // channel: [ val(meta), [ bam ] ]
    ch_region
    ch_ref

    main:

    ch_versions = Channel.empty()

    SAMTOOLS_VIEW ( ch_bam )
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions.first())

    SAMTOOLS_INDEX ( SAMTOOLS_SORT.out.bam )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    emit:
    // TODO nf-core: edit emitted channels
    bam      = SAMTOOLS_SORT.out.bam           // channel: [ val(meta), [ bam ] ]
    bai      = SAMTOOLS_INDEX.out.bai          // channel: [ val(meta), [ bai ] ]
    csi      = SAMTOOLS_INDEX.out.csi          // channel: [ val(meta), [ csi ] ]

    versions = ch_versions                     // channel: [ versions.yml ]
}

