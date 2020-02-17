from tx_annotation import *

# Ran this in three parts

'''
# Part 1 : formatting isoform files
rsem_tpms = "gs://gnomad-public/papers/2019-tx-annotation/data/HBDR_fetal_RNAseq/HBDR.RSEM.tissue_names.tpm.tsv.gz"

get_gtex_summary(rsem_tpms, "gs://gnomad-public/papers/2019-tx-annotation/data/HBDR_fetal_RNAseq/HBDR.tissue_names.medians.021520.mt")

get_gene_expression("gs://gnomad-public/papers/2019-tx-annotation/data/HBDR_fetal_RNAseq/HBDR.tissue_names.medians.021520.mt",
                    "gs://gnomad-public/papers/2019-tx-annotation/data/HBDR_fetal_RNAseq/HBDR.gene_expression_per_gene_per_tissue.021520.ht")


# At this point updated hbdr_fetal_tissue_summary_ht_path, hbdr_fetal_tissue_gene_maximums_ht_path in tx_annotation_resources.py
# Part 2 : annotate context and identify max pext per gene

context_hbdr_annotated = "gs://gnomad-public/papers/2019-tx-annotation/data/HBDR_fetal_RNAseq/all.possible.snvs.tx_annotated.HBDR.021520.ht"
context_hbdr_max_per_gene = "gs://gnomad-public/papers/2019-tx-annotation/data/HBDR_fetal_RNAseq/all.genes.max.pext.HBDR.021520.tsv.bgz"

mt, hbdr = read_tx_annotation_tables(context_ht_path, hbdr_fetal_tissue_summary_ht_path, 'ht')

mt_annotated = tx_annotate_mt(mt, hbdr, "proportion",
                              tissues_to_filter=None,
                              gene_maximums_ht_path=hbdr_fetal_tissue_gene_maximums_ht_path,
                              filter_to_csqs=all_coding_csqs)

mt_annotated.write(context_hbdr_annotated, overwrite=True)

print("Wrote RSEM Brain Cortex HT, now going to identify maximums")

identify_maximum_pext_per_gene(context_hbdr_annotated, context_hbdr_max_per_gene)
'''


# Part 3 : ClinVar anntoation for supplemental figure

clinvar_mt, hbdr = read_tx_annotation_tables(clinvar_ht_path, hbdr_fetal_tissue_summary_ht_path, "ht")
mt, hbbr = read_tx_annotation_tables(gnomad_release_mt_path, hbdr_fetal_tissue_summary_ht_path, "ht")
mt = mt.filter_rows(hl.len(mt.filters) == 0)

genome_mt, hbdr = read_tx_annotation_tables(gnomad_genomes_release_mt_path, hbdr_fetal_tissue_summary_ht_path, "ht")
genome_mt = genome_mt.filter_rows(hl.len(genome_mt.filters) == 0)


hi_genes = import_gene_list(curated_haploinsufficient_genes, gene_column="ENSGID", ensg=True)
out_dir = "gs://gnomad-public/papers/2019-tx-annotation/results/fetal_RNAseq/"

# gnomAD exomes
mt_gnomad_hi = tx_annotate_mt(mt, hbdr, "proportion",
                              tissues_to_filter=None,
                              gene_maximums_ht_path=hbdr_fetal_tissue_gene_maximums_ht_path,
                              filter_to_csqs=lof_csqs,
                              filter_to_genes=hi_genes,
                              gene_column_in_mt="gene_id")

mt_gnomad_hi = mt_gnomad_hi.filter_rows(~hl.is_missing(mt_gnomad_hi.tx_annotation))
mt_gnomad_hi = pull_out_worst_from_tx_annotate(mt_gnomad_hi)
mt_gnomad_hi.rows().export("%sHBDR_HI_genes.gnomad.exomes.r2.1.tx_annotated.021520.tsv.bgz" %out_dir)

# ClinVar
mt_clinvar_hi = tx_annotate_mt(clinvar_mt, hbdr, "proportion",
                               tissues_to_filter=None,
                               gene_maximums_ht_path=hbdr_fetal_tissue_gene_maximums_ht_path,
                               filter_to_csqs=lof_csqs, filter_to_genes=hi_genes,
                               gene_column_in_mt="gene_id")

mt_clinvar_hi = mt_clinvar_hi.filter_rows(~hl.is_missing(mt_clinvar_hi.tx_annotation))
mt_clinvar_hi = pull_out_worst_from_tx_annotate(mt_clinvar_hi)
mt_clinvar_hi = mt_clinvar_hi.annotate_rows(**mt_clinvar_hi.info)
mt_clinvar_hi = mt_clinvar_hi.drop("vep", "tx_annotation","info")
mt_clinvar_hi.rows().export("%sHBDR_HI_genes.clinvar.alleles.single.b37.tx_annotated.021520.tsv.bgz" %out_dir)

# gnomAD genomes
mt_genome_gnomad_hi = tx_annotate_mt(genome_mt, hbdr, "proportion",
                                     tissues_to_filter=None,
                                     gene_maximums_ht_path=hbdr_fetal_tissue_gene_maximums_ht_path,
                                     filter_to_csqs=lof_csqs,
                                     filter_to_genes=hi_genes,
                                     gene_column_in_mt="gene_id")

mt_genome_gnomad_hi = mt_genome_gnomad_hi.filter_rows(~hl.is_missing(mt_genome_gnomad_hi.tx_annotation))
mt_genome_gnomad_hi = pull_out_worst_from_tx_annotate(mt_genome_gnomad_hi)
mt_genome_gnomad_hi.rows().export("%sHBDR_HI_genes.gnomad.genomes.r2.1.tx_annotated.021520.tsv.bgz" %out_dir)