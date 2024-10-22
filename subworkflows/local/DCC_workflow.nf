include { STAR_ALIGN as DCC_1ST_PASS       } from '../../modules/nf-core/star/align/main'
include { STAR_ALIGN as DCC_2ND_PASS       } from '../../modules/nf-core/star/align/main'
include { SJDB as DCC_SJDB                 } from '../../modules/local/star/sjdb/main'
include { STAR_ALIGN as DCC_MATE1_1ST_PASS } from '../../modules/nf-core/star/align/main'
include { STAR_ALIGN as DCC_MATE1_2ND_PASS } from '../../modules/nf-core/star/align/main'
include { SJDB as DCC_MATE1_SJDB           } from '../../modules/local/star/sjdb/main'
include { STAR_ALIGN as DCC_MATE2_1ST_PASS } from '../../modules/nf-core/star/align/main'
include { STAR_ALIGN as DCC_MATE2_2ND_PASS } from '../../modules/nf-core/star/align/main'
include { SJDB as DCC_MATE2_SJDB           } from '../../modules/local/star/sjdb/main'
include { DCC                              } from '../../modules/local/dcc/dcc/main'
include { DCC_FILTER                       } from '../../modules/local/dcc/filter/main'

workflow DCC_WORKFLOW {
    take:
    reads
    gtf
    star_index
    bsj_reads

    main:
    
    // DCC WORKFLOW

    DCC_1ST_PASS( reads, star_index, gtf, true, '', '' )
    DCC_SJDB( DCC_1ST_PASS.out.tab.map{ meta, tab -> return [ tab ] }.collect(), bsj_reads )
    DCC_2ND_PASS( reads, star_index, DCC_SJDB.out.sjtab, true, '', '' )

    mate1 = reads.map{ meta, reads -> return [meta, reads[0] ] }
    DCC_MATE1_1ST_PASS( mate1, star_index, gtf, true, '', '' )
    DCC_MATE1_SJDB( DCC_MATE1_1ST_PASS.out.tab.map{ meta, tab -> return [ tab ] }.collect(), bsj_reads )
    DCC_MATE1_2ND_PASS( mate1, star_index, DCC_MATE1_SJDB.out.sjtab, true, '', '' )

    mate2 = reads.map{ meta, reads -> return [ meta, reads[1] ] }
    DCC_MATE2_1ST_PASS( mate2, star_index, gtf, true, '', '' )
    DCC_MATE2_SJDB( DCC_MATE2_1ST_PASS.out.tab.map{ meta, tab -> return [ tab ] }.collect(), bsj_reads )
    DCC_MATE2_2ND_PASS( mate2, star_index, DCC_MATE2_SJDB.out.sjtab, true, '', '' )

    dcc_stage = DCC_2ND_PASS.out.junction.join( DCC_MATE1_2ND_PASS.out.junction, remainder: true ).join( DCC_MATE2_2ND_PASS.out.junction, remainder: true )
    dcc = dcc_stage.map{ it -> def meta = it[0]; if( meta.single_end ){ return [ it[0], it[1], [], [] ] } else { return it } }.view()
    DCC( dcc, fasta, gtf )
    DCC_FILTER( DCC.out.txt.map{ meta, txt -> meta.tool = "dcc"; return [ meta, txt ] }, bsj_reads )

    ch_versions = ch_versions.mix(DCC_MATE1_1ST_PASS.out.versions)
    ch_versions = ch_versions.mix(DCC.out.versions)

    emit:
    DCC_FILTER.out.matrix
    versions
}