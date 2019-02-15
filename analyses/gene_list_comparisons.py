from tx_annotation import *

clinvar_mt, gtex = read_tx_annotation_tables(clinvar_ht_path, gtex_v7_tx_summary_mt_path, "ht")
mt, gtex = read_tx_annotation_tables(gnomad_release_mt_path, gtex_v7_tx_summary_mt_path, "ht")
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
mt_gnomad_hi.rows().export("%sHI_genes.gnomad.exomes.r2.1.tx_annotated.021519.tsv.bgz" %out_dir)

# AR genes
recessive_genes = import_gene_list(blekhman_berg, gene_column="gene")

mt_gnomad_recessive = tx_annotate_mt(mt, gtex, "proportion", filter_to_genes=recessive_genes,
                                     gene_column_in_mt="gene_symbol",
                                     filter_to_homs=True, filter_to_csqs = lof_csqs)

mt_gnomad_recessive = mt_gnomad_recessive.filter_rows(~hl.is_missing(mt_gnomad_recessive.tx_annotation))
mt_gnomad_recessive = pull_out_worst_from_tx_annotate(mt_gnomad_recessive)
mt_gnomad_recessive = mt_gnomad_recessive.drop("vep", "tx_annotation")
mt_gnomad_recessive.rows().export("%sAR_genes.gnomad.exomes.r2.1.tx_annotated.021519.tsv.bgz" %out_dir)

# LOEUF genes

# pLoFs
mt_gnomad_pli_lof = tx_annotate_mt(mt, gtex, "proportion",
                                   filter_to_genes=oe_genes, gene_column_in_mt="gene_symbol",
                                   filter_to_csqs=lof_csqs)
mt_gnomad_pli_lof = mt_gnomad_pli_lof.filter_rows(~hl.is_missing(mt_gnomad_pli_lof.tx_annotation))
mt_gnomad_pli_lof = mt_gnomad_pli_lof.select_rows('tx_annotation')
mt_gnomad_pli_lof = pull_out_worst_from_tx_annotate(mt_gnomad_pli_lof)
mt_gnomad_pli_lof = mt_gnomad_pli_lof.drop("tx_annotation")

mt_gnomad_pli_lof.rows().export("%sloeuf_genes.plof.gnomad.exomes.r2.1.tx_annotated.021519.tsv.bgz" %out_dir)

# synonymous variants
mt_gnomad_pli_syn = tx_annotate_mt(mt, gtex, "proportion",
                                   filter_to_genes=oe_genes, gene_column_in_mt="gene_symbol",
                                   filter_to_csqs=syn_csqs)
mt_gnomad_pli_syn = mt_gnomad_pli_syn.filter_rows(~hl.is_missing(mt_gnomad_pli_syn.tx_annotation))
mt_gnomad_pli_syn = mt_gnomad_pli_syn.annotate_rows(tx_annotation=mt_gnomad_pli_syn.tx_annotation.map(fix_loftee_beta_nonlofs))
mt_gnomad_pli_syn = mt_gnomad_pli_syn.select_rows('tx_annotation')
mt_gnomad_pli_syn = pull_out_worst_from_tx_annotate(mt_gnomad_pli_syn)
mt_gnomad_pli_syn = mt_gnomad_pli_syn.drop("tx_annotation")
mt_gnomad_pli_syn.rows().export("%sloeuf_genes.syn.gnomad.exomes.r2.1.tx_annotated.021519.tsv.bgz" %out_dir)

# ClinVar

# HI genes
mt_clinvar_hi = tx_annotate_mt(clinvar_mt, gtex,"proportion",
                               filter_to_csqs=lof_csqs, filter_to_genes=hi_genes,
                               gene_column_in_mt="gene_id")
mt_clinvar_hi = mt_clinvar_hi.filter_rows(~hl.is_missing(mt_clinvar_hi.tx_annotation))
mt_clinvar_hi = pull_out_worst_from_tx_annotate(mt_clinvar_hi)
mt_clinvar_hi = mt_clinvar_hi.annotate_rows(**mt_clinvar_hi.info)
mt_clinvar_hi = mt_clinvar_hi.drop("vep", "tx_annotation","info")
mt_clinvar_hi.rows().export("%sHI_genes.clinvar.alleles.single.b37.tx_annotated.021519.tsv.bgz" %out_dir)

# AR genes
mt_clinvar_ar = tx_annotate_mt(clinvar_mt, gtex,"proportion",
                                            filter_to_csqs=lof_csqs,
                                            filter_to_genes=recessive_genes, gene_column_in_mt="gene_symbol")

mt_clinvar_ar = mt_clinvar_ar.filter_rows(~hl.is_missing(mt_clinvar_ar.tx_annotation))
mt_clinvar_ar = pull_out_worst_from_tx_annotate(mt_clinvar_ar)
mt_clinvar_ar = mt_clinvar_ar.annotate_rows(**mt_clinvar_ar.info)
mt_clinvar_ar = mt_clinvar_ar.drop("vep", "tx_annotation","info")
mt_clinvar_ar.rows().export("%sAR_genes.clinvar.alleles.single.b37.tx_annotated.021519.tsv.bgz" %out_dir)

