'''To show that pext adds value on top of Polyphen, slightly editing the tx_annotate script to include Polyphen information
I've commented in the tx_annotate function where this happens'''




'''
Annotating all of gnomAD, but with polyphen in there

This is a terrible terrible way to do this, and I promise to the bioinformatics Gods that I will fix it

from tx_annotation import *
from datetime import datetime


def get_expression_proportion_polyphen(tx_table, tissues_to_filter, gene_maximum_ht):

    if tissues_to_filter:
        print("Filtering tissues:", tissues_to_filter)
        tx_table = tx_table.drop(*tissues_to_filter)

    remaining_tissue_columns = list(
        set(tx_table.row) - {'locus', 'alleles','csq', 'ensg', 'symbol','lof', 'lof_flag','polyphen'})

    tx_table = tx_table.annotate(
        tx_expression=
        {tissue_id: tx_table[tissue_id] for tissue_id in remaining_tissue_columns})

    tx_table = tx_table.key_by('ensg').join(gene_maximum_ht.key_by("ensg"))

    expression_proportion_table = tx_table.annotate(
        expression_proportion_dict=
        {tissue_id: tx_table.tx_expression[tissue_id] / tx_table.gene_maximum_dict[tissue_id] for tissue_id in remaining_tissue_columns})

    columns_to_drop = list(set(expression_proportion_table.row) - {'locus', 'alleles','csq', 'ensg',
                                                                   'symbol','lof', 'lof_flag','polyphen','expression_proportion_dict'})

    expression_proportion_table = expression_proportion_table.drop(*columns_to_drop)

    expression_proportion_table = expression_proportion_table.annotate(
        **{tissue_id: expression_proportion_table.expression_proportion_dict[tissue_id] for tissue_id in remaining_tissue_columns})

    expression_proportion_table = expression_proportion_table.annotate(
        mean_proportion=hl.mean( hl.filter(
            lambda e : ~hl.is_nan(e),  [expression_proportion_table[tissue_id] for tissue_id in remaining_tissue_columns]), filter_missing=True))

    expression_proportion_table = expression_proportion_table.drop(
        expression_proportion_table.expression_proportion_dict).key_by(
        'locus', 'alleles', 'ensg')

    return expression_proportion_table


def tx_annotate_mt_polyphen(mt, gtex, tx_annotation_type,
                   tissues_to_filter = v7_tissues_to_drop, gene_maximums_ht_path = gtex_v7_gene_maximums_ht_path,
                   filter_to_csqs=all_coding_csqs, filter_to_genes=None, gene_column_in_mt=None, filter_to_homs=False,
                   out_tx_annotation_tsv=None, out_tx_annotation_ht=None):


    """
    Annotate variants in the input MatrixTable with transcript-based expression values accross GTEx. Returns Table.
    :param MatrixTable mt: Input variant file
    :param MatrixTable gtex: Input GTEx summary MatrixTable, must have transcript_id column to key by
    :param str tx_annotation_type: One of ["expression", "proportion"]. Select proportion if you'd like the
    tx_annotation values to be normalized by max expression of the gene
    :param None or list filter_to_csqs: Default None. If you'd like to filter the mt before annotating
    (decreases time) feed in a list or set of consequence terms.
    :param str gene_column_in_mt: Must be set if filter_to_genes != None.
    Column in matrix table that contains gene information within vep.transcript_consequences.
    often ["gene_id", "gene_symbol"]
    :param None or list filter_to_csqs: Default None. If you'd like to filter the mt before annotating
    (decreases time) feed in a list or set of consequence terms.
    Example = ["stop_gained","splice_donor_variant", "splice_acceptor_variant","frameshift_variant"]
    :param None or str out_tx_annotation_tsv: Default None.
    If you'd like to write out the results table as a tsv, provide a tsv path
    :param None or str out_tx_annotation_ht: Default None.
    If you'd like to write out the results table as a Hail 0.2 table, provide a .ht path
    :param bool filter_to_homs: Default False
    If True, filter to variants with at least one homozygote in dataset
    :return: Table with columns: variant, worst_csq, ensg, LOFTEE LOF, LOFTEE LOF Flag, transcript-aware expression
    by GTEx Tissue
    :rtype: Table with variants annotated with transcript-aware tissue expression
    """

    #check_inputs(**locals())

    gtex_table = gtex.key_by("transcript_id")

    #mt = process_consequences(mt, penalize_flags=False)
    mt_exploded = mt.distinct_by_row()
    mt_exploded = mt_exploded.annotate_rows(vep=mt_exploded.vep.annotate(
        transcript_consequences=mt_exploded.vep.transcript_consequences.map(add_most_severe_consequence_to_consequence)))

    # Explode the mt for the transcript consequences to be able to key by transcript ID
    mt_exploded = mt_exploded.explode_rows(mt_exploded.vep.transcript_consequences)

    mt_kt = mt_exploded.rows()
    # Currently testing removal of protein coding transcripts
    mt_kt = mt_kt.filter(mt_kt.vep.transcript_consequences.biotype == "protein_coding")

    if filter_to_genes:
        print("Filtering to genes of interest")
        mt_kt = filter_table_to_gene_list(mt_kt, filter_to_genes, gene_column_in_mt)

    if filter_to_csqs:
        print("Filtering to csqs in %s" % (",".join(filter_to_csqs)))
        mt_kt = filter_table_to_csqs(mt_kt, filter_to_csqs)

    if filter_to_homs:
        print("Filtering to variants with at least 1 homozygote sample in dataset")
        #mt_kt = mt_kt.filter(mt_kt.info.Hom[mt_kt.a_index - 1] > 0)
        idx = mt_kt.globals.freq_index_dict['gnomad']
        mt_kt = mt_kt.filter(mt_kt.freq[idx].homozygote_count >= 1)

    # Annotate mt with the gtex values (ie. join them)
    mt_kt = mt_kt.annotate(tx_data=gtex_table[mt_kt.vep.transcript_consequences.transcript_id])

    # Group by gene, worst_csq and variant, and do a pairwise-sum
    grouped_table = (
        mt_kt.group_by(csq=mt_kt.vep.transcript_consequences.most_severe_consequence,
                       ensg=mt_kt.vep.transcript_consequences.gene_id,
                       symbol=mt_kt.vep.transcript_consequences.gene_symbol,
                       locus=mt_kt.locus,
                       alleles=mt_kt.alleles,
                       lof=mt_kt.vep.transcript_consequences.lof,
                       polyphen=mt_kt.vep.transcript_consequences.polyphen_prediction,
                       lof_flag=mt_kt.vep.transcript_consequences.lof_flags).aggregate(tx_annotation=hl.agg.array_sum(mt_kt.tx_data.agg_expression)))

    # Expand the columns from the arrays and add tissues as headers
    #tissue_ids = gtex.tissue.collect()
    # Since gtex no longer has .tissue just a new way to do this, i probably want to save it as a global at some point
    tissue_ids = sorted([y.tissue for y in gtex.values.take(1)[0]])
    d = {tiss: i for i, tiss in enumerate(tissue_ids)}

    tx_annotation_table = grouped_table.annotate(
        **{tissue_id.replace("-", "_").replace(" ", "_").replace("(", "_").replace(")", "_"):
               grouped_table.tx_annotation[d[tissue_id]] for tissue_id in tissue_ids})

    tx_annotation_table = tx_annotation_table.drop(tx_annotation_table.tx_annotation)

    # First of all do you want proportions or expression?
    if tx_annotation_type == "proportion":
        print("Returning expression proportion")
        gene_maximums_ht = hl.read_table(gene_maximums_ht_path)
        tx_annotation_table = get_expression_proportion_polyphen(tx_annotation_table, tissues_to_filter, gene_maximums_ht)

    #You can write the output that is exploded by variants-ensg-csq-symbol-LOFTEE-LOFTEEflag
    # and has a value for each tissue as column, either as a TSV or a KT

    if out_tx_annotation_tsv:
        print("Writing tsv file to %s" %out_tx_annotation_tsv)
        tx_annotation_table.export(out_tx_annotation_tsv)

    if out_tx_annotation_ht:
        print("Writing Table to %s" % out_tx_annotation_ht)
        tx_annotation_table.write(out_tx_annotation_ht)

    tx_annotation_table = tx_annotation_table.key_by(tx_annotation_table.locus, tx_annotation_table.alleles)
    tx_annotation_table = tx_annotation_table.collect_by_key('tx_annotation')
    mt = mt.annotate_rows(**tx_annotation_table[mt.locus, mt.alleles])

    return mt


start = datetime.now()

mt, gtex = read_tx_annotation_tables(gnomad_release_mt_path, gtex_v7_tx_summary_ht_path, 'ht')

mt_annotated = tx_annotate_mt_polyphen(mt, gtex,"proportion",
                                       gene_maximums_ht_path = gtex_v7_gene_maximums_ht_path,
                                       filter_to_csqs=all_coding_csqs)

finished_annotation = datetime.now()
print("Finished annotation ",finished_annotation)

mt_annotated.write("gs://gnomad-public/papers/2019-tx-annotation/data/gnomad_release_annotated/gnomad.exomes.r2.1.1.sites.tx_annotated.withpolyphen.021520.ht",overwrite = True)

stop = datetime.now()
print("Finished writing",stop)
print("total time", stop - start)

'''
'''
Running MAPs

# For code to run, need to git clone onto gnomad_lof on cluster
'''
import sys
sys.path.append('/home/hail/gnomad_lof')

from gnomad_hail import *
from gnomad_hail.resources.sample_qc import *
from gnomad_hail.utils.plotting import *
from constraint_utils import *
from tx_annotation import *

# For code to run, need to git clone onto gnomad_lof on cluster
import sys
sys.path.append('/home/hail/gnomad_lof')

from gnomad_hail import *
from gnomad_hail.resources.sample_qc import *
from gnomad_hail.utils.plotting import *
from constraint_utils import *
from tx_annotation import *

def load_tx_expression_data(tx_ht):
    tx_ht = tx_ht.rows()

    def process_expression_data(csq_expression):
        exprs_to_drop = ['ensg', 'csq', 'symbol', 'lof', 'lof_flag', 'mean_proportion', 'polyphen']
        expression_data = csq_expression.drop(*exprs_to_drop)
        all_tissues = list(expression_data.values())
        expression_data_list = list(zip(list(expression_data), all_tissues))
        brain_tissues = [x[1] for x in expression_data_list if 'Brain' in x[0]]
        return csq_expression.select('ensg', 'csq', 'symbol', 'lof', 'lof_flag','polyphen',
                                     mean_expression=hl.mean(hl.filter(lambda e: ~hl.is_nan(e), all_tissues), filter_missing=True),
                                     mean_brain_expression=hl.mean(hl.filter(lambda f: ~hl.is_nan(f), brain_tissues), filter_missing=True),
                                     Brain_Cortex=csq_expression.Brain_Cortex
                                     )

    return tx_ht.annotate(tx_annotation=tx_ht.tx_annotation.map(process_expression_data))

context_ht = hl.read_table(context_ht_path)

# Import and process gnomad 2.1.1 transcript annotation
ht = hl.read_matrix_table('gs://gnomad-public/papers/2019-tx-annotation/data/gnomad_release_annotated/gnomad.exomes.r2.1.1.sites.tx_annotated.withpolyphen.021520.ht')
ht = ht.filter_rows(~hl.is_missing(ht.tx_annotation))
ht = ht.annotate_rows(tx_annotation = ht.tx_annotation.map(fix_loftee_beta_nonlofs))
ht = load_tx_expression_data(ht)
ht = hl.MatrixTable.from_rows_table(ht)
ht = pull_out_worst_from_tx_annotate(ht)

# Only consider variants that pass RF
ht = ht.rows()
ht = ht.filter(hl.len(ht.filters) == 0)
context = context_ht[ht.key]
ht = ht.annotate(context=context.context, methylation=context.methylation)
ht = prepare_ht(ht, trimer=True, annotate_coverage=False)

# Prepare MAPS data
even_breaks = [0.999, 0.995, 0.99, 0.98] + list(map(lambda x: x/40, range(39, -1, -1)))

ht = ht.filter(ht.freq[0].AN > 125748 * 0.8 * 2)
mutation_ht = hl.read_table(mutation_rate_ht_path)


# Only consider LOFTEE HC pLoFs, missense and synonymous
ht = ht.annotate(keep = hl.case(missing_false=True)
                 .when(ht.csq == "missense_variant", "keep").default('filter'))

ht = ht.annotate(keep = hl.case(missing_false=True)
                 .when((ht.csq == "stop_gained") &(ht.lof == 'HC'), "keep")
                 .when((ht.csq == "splice_donor_variant") &(ht.lof == 'HC'), "keep")
                 .when((ht.csq == "splice_acceptor_variant" ) &(ht.lof == 'HC'), "keep")
                 .when((ht.csq == "missense_variant") & (ht.polyphen == 'possibly_damaging'), "keep")
                 .when((ht.csq == "missense_variant") & (ht.polyphen == 'probably_damaging'), "keep")
                 .when((ht.csq == "missense_variant") & (ht.polyphen == 'benign'), "keep")
                 .when(ht.csq == "synonymous_variant", "keep").default('filter'))

ht = ht.filter(ht.keep == "keep")



# Group pLoFs, remember can't calculate MAPs on frameshifts (no mutational model)
ht = ht.annotate(worst_csq = hl.case(missing_false=True)
                 .when(ht.csq == "stop_gained", "pLoF")
                 .when(ht.csq == "splice_donor_variant", "pLoF")
                 .when(ht.csq == "splice_acceptor_variant", "pLoF")
                 .when(ht.csq == "missense_variant", "missense_variant")
                 .when(ht.csq == "synonymous_variant", "synonymous_variant").default('irrev_var'),
                 lof = ht.lof)


print("finished processing")

constraint = hl.read_table(constraint_ht_path)
constraint = constraint.rename({"gene": "symbol"})
constraint = constraint.key_by("symbol")
ht = ht.key_by("symbol")

ht_constraint = ht.annotate(constraint_bin = constraint[ht.symbol].oe_lof_upper_bin,
                            constraint_value = constraint[ht.symbol].oe_lof_upper)

# Addded in filtering for max pext low genes
genes_to_filter = hl.import_table("gs://gnomad-public/papers/2019-tx-annotation/data/GRCH37_hg19/max_pext_low_genes.021520.tsv", force = True)
genes_to_filter = genes_to_filter.key_by('symbol')

ht_constraint = ht_constraint.filter(~hl.is_defined(genes_to_filter[ht_constraint.key]))


def run_maps_constraint_binexport(f, write, mut_ht = mutation_ht):
    m = maps(f, mut_ht, ['constraint_bin', 'polyphen'])
    m.export(write)

# run_maps_constraint_binexport(ht_constraint,
#                                "gs://gnomad-berylc/tx-annotation/hail2/reviewer_response/polyphen/try5.withpolyphen.all.maps.no.expression.eachplofcategoryexpression.070219.tsv.bgz")


oe_constraint_bin_below_01 = ht_constraint.filter(ht_constraint.mean_expression < 0.1)
run_maps_constraint_binexport(oe_constraint_bin_below_01,
                             "gs://gnomad-public/papers/2019-tx-annotation/results/maps/maps.withpolyphen.low.021520.tsv.bgz")
print('wrote low')

oe_constraint_bin_above_09 = ht_constraint.filter(ht_constraint.mean_expression > 0.9)
run_maps_constraint_binexport(oe_constraint_bin_above_09,
                              "gs://gnomad-public/papers/2019-tx-annotation/results/maps/maps.withpolyphen.high.021520.tsv.bgz")

print('wrote high')

oe_constraint_bin_between =  ht_constraint.filter((ht_constraint.mean_expression <= 0.9) & (ht_constraint.mean_expression >= 0.1))
run_maps_constraint_binexport(oe_constraint_bin_between,
                              "gs://gnomad-public/papers/2019-tx-annotation/results/maps/maps.withpolyphen.medium.021520.tsv.bgz")
