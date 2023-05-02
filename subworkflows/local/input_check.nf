//
// Check input samplesheet and get bam channels
//

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'

workflow SAMPLE_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    SAMPLESHEET_CHECK ( samplesheet )
        .csv
        .splitCsv ( header:true, sep:',' )
        .map { create_bam_channel(it) }
        .set { bam }

    emit:
    bam                                       // channel: [ val(meta), [ bam, bai ] ]
    versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}

workflow REGION_CHECK {
    take:
    regionsheet // file: /path/to/samplesheet.csv

    main:
    Channel.fromPath ( regionsheet )
        .splitCsv ( header:true, sep:',' )
        .map { create_region_channel(it) }
        .set { region }

    emit:
    region                                    // channel: [meta, fasta, region ]
}

workflow DEPTH_CHECK {
    take:
    depthsheet // file: /path/to/depthsheet.csv

    main:
    Channel.fromPath ( depthsheet )
        .splitCsv ( header:true, sep:',' )
        .map { create_depth_channel(it) }
        .set { depth }

    emit:
    depth                                    // channel: [ depth ]
}
workflow PANEL_CHECK {
    take:
    panelsheet // file: /path/to/panelsheet.csv

    main:
    Channel.fromPath ( panelsheet )
        .splitCsv ( header:true, sep:',' )
        .map { create_panel_channel(it) }
        .set { panel }

    emit:
    panel                                    // channel: [meta,  panel, index ]
}

// Function to get list of [ meta, [ bam ] ]
def create_bam_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id         = row.sample_id
    meta.ref        = row.ref_id

    // add path(s) of the bam file(s) to the meta map
    def bam_meta = []
    bam_meta = [ meta, [file(row.bam)], [file(row.bam_index)] ]
    return bam_meta
}

// Function to get list of [ region(chr:start-end) ]
def create_region_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.ref        = row.ref
    meta.region     = "$row.chr:$row.start-$row.end"
    // colapse regions in chr:start-end format
    def region_meta = []
    region_meta = [meta, file(row.fasta), "$row.chr:$row.start-$row.end"]
    return region_meta
}

// Function to get list of panel
def create_panel_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.ref        = row.ref_id
    meta.panel      = row.panel_id
    meta.id         = row.panel_id

    // add path(s) of the vcf file(s) to the meta map
    def panel_meta = []
    panel_meta = [meta, [file(row.file)], [file(row.file_index)]]
    return panel_meta
}

// Function to get list of [ depth ]
def create_depth_channel(LinkedHashMap row) {
    def depth_meta = []
    depth_meta = "$row.depth"
    return depth_meta
}
