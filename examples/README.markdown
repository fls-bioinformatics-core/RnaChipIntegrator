RnaChipIntegrator: examples
===========================

This directory contains three example data files:

 * `ExpressionData.txt` has sample expression data
 * `ChIP_peaks.txt` has sample ChIP peak summit data
 * `ChIP_edges.txt` has sample ChIP peak region data

RnaChipIntegrator can be run on these data files as follows:

    RnaChipIntegrator.py --window=130000 ExpressionData.txt ChIP_edges.txt
    RnaChipIntegrator.py --window=130000 ExpressionData.txt ChIP_peaks.txt

--