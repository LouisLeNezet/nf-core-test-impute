#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SAMPLE_CHECK                           } from '../subworkflows/local/input_check.nf'
include { REGION_CHECK                           } from '../subworkflows/local/input_check.nf'
include { DEPTH_CHECK                            } from '../subworkflows/local/input_check.nf'
include { PANEL_CHECK                            } from '../subworkflows/local/input_check.nf'
include { BAM_SIMULATE                           } from '../subworkflows/local/bam_simulate.nf'
include { GET_PANEL                              } from '../subworkflows/local/get_panel.nf'
include { COMPUTE_GL as GL_TRUTH                 } from '../subworkflows/local/compute_gl.nf'
include { COMPUTE_GL as GL_EMUL                  } from '../subworkflows/local/compute_gl.nf'
include { VCF_IMPUTE_GLIMPSE                     } from '../subworkflows/nf-core/vcf_impute_glimpse/main.nf'
include { MULTIPLE_IMPUTE_GLIMPSE2               } from '../subworkflows/nf-core/multiple_impute_glimpse2/main'
include { GLIMPSE_CONCORDANCE                    } from '../modules/nf-core/glimpse/concordance/main.nf'
include { GUNZIP                                 } from '../modules/nf-core/gunzip/main'
include { ADD_COLUMNS                            } from '../modules/local/add_columns.nf'
include { CONCATENATE                            } from '../modules/local/concatenate.nf'

workflow TESTQUALITY {
    SAMPLE_CHECK(params.input)
    REGION_CHECK(params.region)
    DEPTH_CHECK(params.depth)
    PANEL_CHECK(params.panel)

    BAM_SIMULATE(
        SAMPLE_CHECK.out.bam,
        REGION_CHECK.out.region,
        DEPTH_CHECK.out.depth
    )

    GET_PANEL(
        PANEL_CHECK.out.panel,
        REGION_CHECK.out.region,
        params.chr_rename
    )

    GL_TRUTH(
        BAM_SIMULATE.out.bam_region
            .combine(BAM_SIMULATE.out.bam_region_index, by:0),
        REGION_CHECK.out.region,
        GET_PANEL.out.panel_sites,
        GET_PANEL.out.panel_tsv
    )

    GL_EMUL(
        BAM_SIMULATE.out.bam_emul
            .combine(BAM_SIMULATE.out.bam_emul_index, by:0),
        REGION_CHECK.out.region,
        GET_PANEL.out.panel_sites,
        GET_PANEL.out.panel_tsv
    )

    impute_input = GL_EMUL.out.vcf
                    .combine(GL_EMUL.out.tbi, by:0)
                    .combine(Channel.of([[]]))
                    .map{meta, vcf, index, sample -> [meta, vcf, index, sample, meta.region]}

    VCF_IMPUTE_GLIMPSE(impute_input,
                        GET_PANEL.out.panel_phased
                            .combine(GET_PANEL.out.panel_phased_index, by:0),
                        Channel.of([["map": "map1"], params.map]).collect())

    //MULTIPLE_IMPUTE_GLIMPSE2(impute_input,
    //                    GET_PANEL.out.panel_phased
    //                        .combine(GET_PANEL.out.panel_phased_index, by:0),
    //                    Channel.of([[]]).collect(),
    //                    Channel.of([[],[],[]]).collect(),
    //                    "sequential")

    glimpse1_vcf = VCF_IMPUTE_GLIMPSE.out.merged_variants
                                    .map{meta, vcf -> [meta + ["soft":"glimpse1"], vcf]}

    //glimpse2_vcf = MULTIPLE_IMPUTE_GLIMPSE2.out.merged_variants
    //                                .map{meta, vcf -> [meta + ["soft":"glimpse2"], vcf]}

    //all_vcf = glimpse1_vcf.concat(glimpse2_vcf)

    ch_concordance = glimpse1_vcf
                    .map{metaIRRPDS, vcf -> [metaIRRPDS.subMap(["id","ref","region","panel"]), metaIRRPDS, vcf]}
                    .combine(GL_TRUTH.out.vcf, by:0)
                    .map{metaIRRP, metaIRRPDS, emul, truth -> [metaIRRPDS.subMap(["ref","region","panel"]), metaIRRPDS, emul, truth]}
                    .combine(GET_PANEL.out.panel_sites.map{metaIpRRP, vcf ->
                                [metaIpRRP.subMap(["ref","region","panel"]),metaIpRRP,vcf]},
                            by:0)
                    .map{metaIRRP,metaIRRPDS, emul, truth, metaIpRRP, freq ->
                        [metaIRRPDS, metaIRRPDS.region, freq, truth, emul]}

    GLIMPSE_CONCORDANCE ( ch_concordance, [], [], [])
    GUNZIP(GLIMPSE_CONCORDANCE.out.errors_grp)
    ADD_COLUMNS(GUNZIP.out.gunzip)

    CONCATENATE(ADD_COLUMNS.out.txt
                    .map{meta, txt -> [["id":"TestQuality"], txt]}
                    .groupTuple()
    )
}
