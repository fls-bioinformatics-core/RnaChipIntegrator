What it does
============

``RnaChipIntegrator`` performs integrated analyses of feature data (any set of
genomic features, for example canonical genes or CpG islands) with ‘peak’ data
(a set of regions, for example ChIP peaks). For each peak it reports the
'nearest' features to each peak, and the 'nearest' peaks for each feature.

For both comparisons only peaks and features on the same chromosome are
considered, and by default distances are calculated from the feature TSS
position to the peak edge which is closest. It is also possible to perform
the analyses using the nearest of either the TSS or TES positions.
