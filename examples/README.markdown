RnaChipIntegrator: examples
===========================

This directory contains three example data files:

 * `ExpressionData.txt` has sample expression data
 * `ChIP_summits.txt` has sample ChIP peak summit data
 * `ChIP_regions.txt` has sample ChIP peak region data

RnaChipIntegrator can be run on these data files as follows:

    RnaChipIntegrator --window=130000 ExpressionData.txt ChIP_regions.txt
    RnaChipIntegrator --window=130000 ExpressionData.txt ChIP_summits.txt

--