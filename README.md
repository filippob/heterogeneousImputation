# heterogeneousImputation
### Imputation experiments on multi-breed/accession/variety plant and animal populations

The objective is to assess the effect of **population structure** on the **accuracy of imputing missing SNP genotypes**.
Data from commercial SNP-arrays used to genotype multiple cattle and sheep breeds are available.

The within-breed imputation accuracy will be compared to the multi-breed (multiple breeds in the training and testing datasets) and across-breed (one breed to impute the other) accuracy.

Two main scenarios are envisaged: 

1. gap-filling: after genotyping, a small proportion of locus-sample cells are left uncalled. These missing genotypes usually need to be imputed before moving on to downstream analyses.
2. mixed geneotyping strategies: to optimize costs, typically a fraction of the samples is genotyped with a high-density SNP array, while the remaining samples are genotyped with a low density array. Here the objective is to impute from low to high density (more challenging scenario)

Go to this repository's [wiki](https://github.com/filippob/heterogeneousImputation/wiki)
