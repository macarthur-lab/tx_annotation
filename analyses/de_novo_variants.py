from tx_annotation import *

rt = hl.import_table("gs://gnomad-berylc/tx-annotation/hail2/DeNovoSignal/data/updated_includeASC_Hamdam/ASD_ddid_Updated.hamdam.ASC.011418.txt")
rt = rt.annotate(variant=rt.CHROM + ':' + rt.POSITION + ":" + rt.REF + ":" + rt.ALT)
rt = rt.annotate(** hl.parse_variant(rt.variant))
rt = rt.key_by(rt.locus, rt.alleles)
mt = hl.MatrixTable.from_rows_table(rt)
mt = mt.repartition(10)
annotated_mt = hl.vep(mt, vep_config)
annotated_mt.write("gs://gnomad-public/papers/2019-tx-annotation/results/de_novo_variants/asd_ddid_de_novos.vepped.021819.mt")