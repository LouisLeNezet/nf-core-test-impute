include { GLIMPSE_CHUNK as CHUNK_PANEL      } from '../../modules/nf-core/glimpse/chunk/main'
include { SHAPEIT5_PHASECOMMON              } from '../../modules/nf-core/shapeit5/phasecommon/main'  
include { SHAPEIT5_LIGATE                   } from '../../modules/nf-core/shapeit5/ligate/main'
include { BCFTOOLS_INDEX as VCF_INDEX1      } from '../../modules/nf-core/bcftools/index/main.nf'
include { BCFTOOLS_INDEX as VCF_INDEX2      } from '../../modules/nf-core/bcftools/index/main.nf'


workflow PHASE_PANEL {
    take:
    ch_vcf   // channel: [ [id, ref], vcf ]
    ch_ref      // channel (mandatory): [ meta, vcf, csi ]
    ch_scaffold  // channel  (optional): path to scaffold
    ch_map    // channel  (optional): path to map

    main:
    ch_versions = Channel.empty()

    input_chunk = ch_vcf
                    .map{meta, vcf, index ->
                        [meta, vcf, index, meta.region]}

    CHUNK_PANEL ( input_chunk )
    ch_versions = ch_versions.mix(CHUNK_PANEL.out.versions)

    chunk_output = CHUNK_PANEL.out.chunk_chr
                                .splitCsv(header: ['ID', 'Chr', 'RegionIn', 'RegionOut', 'Size1', 'Size2'], sep: "\t", skip: 0)
                                .map { meta, it -> [meta, it["RegionIn"]]}

    phase_input = ch_vcf
                    .combine(chunk_output, by:0)
                    .combine(Channel.of([[]]))
                    .map{meta, vcf, index, region, pedigree ->
                            [meta + ["region_chunk": region], vcf, index, region, pedigree]}

    SHAPEIT5_PHASECOMMON ( phase_input,
                            ch_ref,
                            ch_scaffold,
                            ch_map ) // [meta, vcf, index, regionin, regionout, regionindex, ref, ref_index, map, sample_infos]
    ch_versions = ch_versions.mix(SHAPEIT5_PHASECOMMON.out.versions.first())

    VCF_INDEX1(SHAPEIT5_PHASECOMMON.out.phased_variant)
    ch_versions = ch_versions.mix(VCF_INDEX1.out.versions.first())

    ligate_input = SHAPEIT5_PHASECOMMON.output.phased_variant
                                        .map{metaIRPRRc, vcf ->
                                         [metaIRPRRc.subMap(["ref","panel","id","region"]),vcf]}
                                        .groupTuple()
                                        .combine(VCF_INDEX1.out.csi
                                            .map{metaIRPRRc, csi ->
                                                [metaIRPRRc.subMap(["ref","panel","id","region"]), csi]}
                                            .groupTuple(),
                                            by:0)
    ligate_input_sorted = ligate_input
        .map{meta, vcf, csi ->
            [meta,
            vcf.sort { a, b ->
                def aStart = a.getName().split(':')[-1].split('-')[0].toInteger()
                def bStart = b.getName().split(':')[-1].split('-')[0].toInteger()
                aStart <=> bStart
            },
            csi]
        }.view()
    SHAPEIT5_LIGATE(ligate_input_sorted) 
    ch_versions = ch_versions.mix(SHAPEIT5_LIGATE.out.versions.first())

    VCF_INDEX2(SHAPEIT5_LIGATE.out.merged_variants)
    ch_versions = ch_versions.mix(VCF_INDEX2.out.versions.first())

    emit:
    panel_phased        = SHAPEIT5_LIGATE.out.merged_variants
    panel_phased_index  = VCF_INDEX2.out.csi
    versions            = ch_versions      // channel: [ versions.yml ]
}