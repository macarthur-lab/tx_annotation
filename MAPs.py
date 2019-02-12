from tx_annotation import *


def import_and_vep_exac_maps_info(f, out_f):
    """

    :param f: In gs://gnomad-berylc/tx-annotation/hail2/MAPs/data/
    ExAC.info.for.maps.adding.consequence.no.rows.txt or
    ExAC.info.for.maps.splice.2.txt

    :param out_f: Out file in gs://gnomad-berylc/tx-annotation/hail2/MAPs/data/
    e.g. ExAC.info.for.maps.splice.vepped.mt"
    :return:
    """
    rt = hl.import_table(f)
    rt = rt.annotate(variant=rt.chrom + ':' + rt.pos + ":" + rt.ref + ":" + rt.alt)
    rt = rt.annotate(** hl.parse_variant(rt.variant))
    rt = rt.key_by(rt.locus, rt.alleles)
    mt = hl.MatrixTable.from_rows_table(rt, partition_key='locus')
    mt = mt.repartition(2500)
    annotated_mt = hl.vep(mt, vep_config)
    annotated_mt.write(out_f)




#Annotation pipeline for stop gained
mt, gtex = read_tx_annotation_tables("gs://gnomad-berylc/tx-annotation/hail2/MAPs/data/ExAC.info.for.maps.adding.consequence.no.rows.vepped.mt",  gtex_v7_tx_summary_mt_path)

stop_gained_tx_table = tx_annotate_mt(mt, gtex, filter_to_csqs=all_coding_csqs)

v7_gene_maximums_kt = hl.read_table("gs://gnomad-berylc/tx-annotation/hail2/data/GTEx.v7.max_expression_per_gene_per_tissue.031318.kt")

stop_gained_tx_proportions = get_expression_proportion(stop_gained_tx_table, v7_tissues_to_drop, v7_gene_maximums_kt)



#Annotation pipeline for splice
mt, gtex = read_tx_annotation_tables("gs://gnomad-berylc/tx-annotation/hail2/MAPs/data/ExAC.info.for.maps.splice.vepped.mt",  gtex_v7_tx_summary_mt_path)

splice_tx_table = tx_annotate_mt(mt, gtex, filter_to_csqs=all_coding_csqs)
v7_gene_maximums_kt = hl.read_table("gs://gnomad-berylc/tx-annotation/hail2/data/GTEx.v7.max_expression_per_gene_per_tissue.031318.kt")
splice_tx_proportions = get_expression_proportion(splice_tx_table, v7_tissues_to_drop, v7_gene_maximums_kt)
splice_tx_proportions = splice_tx_proportions.drop(splice_tx_proportions.expression_proportion_dict)
splice_tx_proportions.write("gs://gnomad-berylc/tx-annotation/hail2/MAPs/data/ExAC.info.for.maps.splice.proportions.tsv.bgz")


