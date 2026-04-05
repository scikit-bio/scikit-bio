# Desert Soil Biocrust Wetting Dataset

This dataset contains paired microbiome and metabolome data derived from biocrust samples after a laboratory wetting event. Adapted from:

- https://github.com/biocore/mmvec/tree/master/examples/soils

The dataset was used as a demonstration in the MMvec paper:

- Morton, J. T., Aksenov, A. A., Nothias, L. F., Foulds, J. R., Quinn, R. A., Badri, M. H., ... & Knight, R. (2019). Learning representations of microbe–metabolite interactions. Nature Methods, 16(12), 1306-1314.

The original dataset was published in:

- Swenson, T. L., Karaoz, U., Swenson, J. M., Bowen, B. P., & Northen, T. R. (2018). Linking soil biology and chemistry in biological soil crust using isolate exometabolomics. Nature Communications, 9(1), 19.

The dataset contains 19 samples:
- 4 time points: 3 min (immediately after wetting), 9 hours, 18 hours, 42 hours, and 49.5 hours.
- 4 biocrust successional stages: early, early-mid, late-mid, and late.

The dataset contains two omics layers:
- Microbes (n=466), generated using ribosomal protein L15 (rplO) genotyping.
- Metabolites (n=85), generated using liquid chromatography-mass spectrometry (LC/MS).
