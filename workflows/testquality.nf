#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { VCF_IMPUTE_GLIMPSE  } from '../subworkflows/nf-core/vcf_impute_glimpse/main.nf'
include { GLIMPSE_CONCORDANCE } from '../modules/nf-core/glimpse/concordance/main.nf'

workflow TESTQUALITY {
    input_vcf = Channel.of([
        [ id:'input', single_end:false ], // meta map
        file("/groups/dog/llenezet/test-datasets/data/ind/12559.chr38.1x.vcf.gz",
            checkIfExists: true),
        file("/groups/dog/llenezet/test-datasets/data/ind/12559.chr38.1x.vcf.gz.csi",
            checkIfExists: true),
        "38"
    ])
    ref_panel = Channel.of([
        [ id:'reference', single_end:false ], // meta map
        file("/groups/dog/llenezet/test-datasets/data/panel/DVDBC.chr38.bcf",
            checkIfExists: true),
        file("/groups/dog/llenezet/test-datasets/data/panel/DVDBC.chr38.bcf.csi",
            checkIfExists: true)
    ]).collect()
    map = Channel.of([[]]).collect()
    sample = Channel.of([[]]).collect()
    VCF_IMPUTE_GLIMPSE ( input_vcf, ref_panel, [], sample )

    allele_freq = [file("/groups/dog/llenezet/test-datasets/data/panel/DVDBC.chr38.sites.vcf.gz",checkIfExists:true),
                   file("/groups/dog/llenezet/test-datasets/data/panel/DVDBC.chr38.sites.vcf.gz.csi",checkIfExists:true)]
    
    truth = [file("/groups/dog/llenezet/test-datasets/data/ind/12559.chr38.bcf",checkIfExists:true),
            file("/groups/dog/llenezet/test-datasets/data/ind/12559.chr38.bcf.csi",checkIfExists:true)]
    
    list_inputs = Channel.of(["38", allele_freq[0], truth[0]])
                    .combine(VCF_IMPUTE_GLIMPSE.out.merged_variants.map{it[1]}.collect().map{it[0]})
                    .collect()
    concordance_input=Channel.of([[ id:'input', single_end:false ]]).combine(list_inputs)

    GLIMPSE_CONCORDANCE ( concordance_input, [], [], []) // meta, Region, Frequencies, Truth, Estimate, minPROB, minDP, bins
}