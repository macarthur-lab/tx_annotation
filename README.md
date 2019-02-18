# Transcript expression-aware annotation 

Welcome to the repository for the "Transcript expression-aware annotation improves rare variant discovery and interpretation" manuscript. Here we'll outline how to get transcript expression values for your variant file and isoform expression expression matrix of interest, and outline the commands and code to recreate analyses in the pre-print. 

## Applying transcript expression aware annotation to your own dataset

You will need 
  1) A variant file that has columns for chrom, pos, ref and alt
  2) An isoform expression matrix 
  3) The ability to use Hail locally or on a cloud platform 
  
You can have additional columns in your variant file, which will be maintained, and only new columns of transcript-expression annotation will be added. Your isoform expression matrix must start with two columns : 1. transcript id and 2. gene_id. The remaining columns can be any samples or tissues. If you have biological replicates, they should be numbered with a '.; delimiter (e.g. MuscleSkeletal.1, MuscleSkeletal.2, MuscleSkeletal.3). 
 
Instructions to set up Hail are laid out in the [Hail docs](https://hail.is/docs/0.2/)

If you're unable to set up Hail in your local environment, we have also release the pext values for every possible SNV in the genome: gs://gnomad-public/papers/2019-tx-annotation/pre_computed/all.possible.snvs.tx_annotated.021819.tsv.bgz

More information about this files format is below.

Please be aware that while we don't expect any issues, the files may be iterated upon until publication of the manuscript! 

We will walk through an example of annotating *de novo* variants in autism and developmental delay / intellectual disability. 

###### 1) Prepare the variant file 
The variant file is available at : gs://gnomad-public/papers/2019-tx-annotation/data/asd_ddid_de_novos.txt


###### 2) Prepare the isoform expression file 

###### 3) Add pext values

## Analyses in manuscript 

###### Conservation and pext
d
