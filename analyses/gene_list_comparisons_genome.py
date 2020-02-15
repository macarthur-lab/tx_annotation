from tx_annotation import *

mt, gtex = read_tx_annotation_tables(gnomad_genomes_release_mt_path, gtex_v7_tx_summary_ht_path, "ht")
mt = mt.filter_rows(hl.len(mt.filters) == 0)

oe_genes = import_gene_list(constraint, gene_column="gene", oe_threshold=0.35)
hi_genes = import_gene_list(curated_haploinsufficient_genes, gene_column="ENSGID", ensg=True)
recessive_genes = import_gene_list(recessive_disease_genes, gene_column="gene")

out_dir = "gs://gnomad-public/papers/2019-tx-annotation/results/gene_list_comparisons/"

#gnomAD

# HI genes
mt_gnomad_hi = tx_annotate_mt(mt, gtex,"proportion",
                              filter_to_csqs=lof_csqs,
                              filter_to_genes=hi_genes, gene_column_in_mt="gene_id")
mt_gnomad_hi = mt_gnomad_hi.filter_rows(~hl.is_missing(mt_gnomad_hi.tx_annotation))
mt_gnomad_hi = pull_out_worst_from_tx_annotate(mt_gnomad_hi)
mt_gnomad_hi.rows().export("%sHI_genes.gnomad.genomes.r2.1.tx_annotated.021619.tsv.bgz" %out_dir)

# AR genes

mt_gnomad_recessive = tx_annotate_mt(mt, gtex, "proportion", filter_to_genes=recessive_genes,
                                     gene_column_in_mt="gene_symbol",
                                     filter_to_homs=True, filter_to_csqs = lof_csqs)

mt_gnomad_recessive = mt_gnomad_recessive.filter_rows(~hl.is_missing(mt_gnomad_recessive.tx_annotation))
mt_gnomad_recessive = pull_out_worst_from_tx_annotate(mt_gnomad_recessive)
mt_gnomad_recessive = mt_gnomad_recessive.drop("vep", "tx_annotation")
mt_gnomad_recessive.rows().export("%sAR_genes.gnomad.genomes.r2.1.tx_annotated.021519.tsv.bgz" %out_dir)

# LOEUF genes

# pLoFs
mt_gnomad_pli_lof = tx_annotate_mt(mt, gtex, "proportion",
                                   filter_to_genes=oe_genes, gene_column_in_mt="gene_symbol",
                                   filter_to_csqs=lof_csqs)
mt_gnomad_pli_lof = mt_gnomad_pli_lof.filter_rows(~hl.is_missing(mt_gnomad_pli_lof.tx_annotation))
mt_gnomad_pli_lof = mt_gnomad_pli_lof.select_rows('tx_annotation')
mt_gnomad_pli_lof = pull_out_worst_from_tx_annotate(mt_gnomad_pli_lof)
mt_gnomad_pli_lof = mt_gnomad_pli_lof.drop("tx_annotation")

mt_gnomad_pli_lof.rows().export("%sloeuf_genes.plof.gnomad.genomes.r2.1.tx_annotated.021519.tsv.bgz" %out_dir)

# synonymous variants
mt_gnomad_pli_syn = tx_annotate_mt(mt, gtex, "proportion",
                                   filter_to_genes=oe_genes, gene_column_in_mt="gene_symbol",
                                   filter_to_csqs=syn_csqs)
mt_gnomad_pli_syn = mt_gnomad_pli_syn.filter_rows(~hl.is_missing(mt_gnomad_pli_syn.tx_annotation))
mt_gnomad_pli_syn = mt_gnomad_pli_syn.annotate_rows(tx_annotation=mt_gnomad_pli_syn.tx_annotation.map(fix_loftee_beta_nonlofs))
mt_gnomad_pli_syn = mt_gnomad_pli_syn.select_rows('tx_annotation')
mt_gnomad_pli_syn = pull_out_worst_from_tx_annotate(mt_gnomad_pli_syn)
mt_gnomad_pli_syn = mt_gnomad_pli_syn.drop("tx_annotation")
mt_gnomad_pli_syn.rows().export("%sloeuf_genes.syn.gnomad.genomes.r2.1.tx_annotated.021519.tsv.bgz" %out_dir)
