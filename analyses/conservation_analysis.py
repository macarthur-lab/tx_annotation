from tx_annotation import *


phylocsf = hl.import_table(phylocsf_file_path, impute = True)
all_baselevel_ht = hl.read_table(all_baselevel_ht_path)

# Remove tissues for consistency and recalculate the mean
all_baselevel_ht = all_baselevel_ht.drop(*v7_tissues_to_drop)
tissues = set(all_baselevel_ht.row) - {'ensg', 'symbol', 'locus', 'mean_proportion', 'mean_prop_correct'}
all_baselevel_ht = all_baselevel_ht.annotate(mean_prop_conservation=hl.mean(
    hl.filter(lambda e: ~hl.is_nan(e),[all_baselevel_ht[tissue_id] for tissue_id in tissues]), filter_missing=True))

# Define regions of high and low conservation, and filter remaning regions
phylocsf = phylocsf.annotate(conservation_type = hl.case(missing_false=True)
                             .when(phylocsf.max_score > 1000, "high")
                             .when(phylocsf.max_score < -100, "low")
                             .default('filter'))
phylocsf = phylocsf.filter(phylocsf.conservation_type != "filter")

# Make intervals to filter baselevel file
phylocsf = phylocsf.annotate(chrom = phylocsf.chromosome_name.replace("chr",""))

phylocsf = phylocsf.annotate(
    interval = hl.interval(hl.locus(phylocsf.chrom, phylocsf.start_coordinate), hl.locus(phylocsf.chrom, phylocsf.end_coordinate)),
    interval_name = phylocsf.chrom + ":" + hl.str(phylocsf.start_coordinate) + "-" + hl.str(phylocsf.end_coordinate) )
phylocsf = phylocsf.key_by(phylocsf.interval)

# Annotate baselevel expression file with whether bases are in the intervals
all_baselevel_ht = all_baselevel_ht.annotate(**phylocsf[all_baselevel_ht.locus])
all_baselevel_ht= all_baselevel_ht.filter(hl.is_defined(all_baselevel_ht.conservation_type), keep=True)

# Get mean pext in these intervals
mean_proportion_in_interval = (all_baselevel_ht.group_by(symbol = all_baselevel_ht.symbol,
                                                         ensg = all_baselevel_ht.ensg,
                                                         enst = all_baselevel_ht.transcript_id,
                                                         interval = all_baselevel_ht.interval_name,
                                                         conservation_type = all_baselevel_ht.conservation_type).
                               aggregate(mean_of_mean_pext =
                                        hl.agg.filter(~hl.is_nan(all_baselevel_ht.mean_prop_conservation),
                                                      hl.agg.mean(all_baselevel_ht.mean_prop_conservation))))

mean_proportion_in_interval.export("gs://gnomad-public/papers/2019-tx-annotation/results/conservation.phylocsf.vs.pext.021219.tsv.bgz")