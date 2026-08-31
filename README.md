# IntraTalker

This package enables condition- and cell-type-specific analysis of transcription factor activities and generates intracellular signaling Networks. The results can be combined with ligand-receptor predictions and analyzed differentially using CrossTalkeR. 

To perform receptor perturbation simulation, use the Python package [IntraTalkerpy](https://github.com/CostaLab/IntraTalkerpy)

A **full tutorial** including an example running in R can be found in the IntraTalkerpy documentation: <br>
[Il1rn KO Case Study](https://costalab.github.io/IntraTalkerpy/tutorials/il1rn_bone_marrow/index.html).

## Abstract

Single-cell sequencing has advanced the study of cell-cell communication, yet most methods focus on intercellular ligand-
receptor interactions while neglecting downstream intracellular signalling cascades and the possibility that downstream
target genes themselves encode ligands, thereby propagating communication across multiple cells. We present In-
traTalker+CrossTalkeR that combines intracellular (IntraTalker) and intercellular (CrossTalkeR) signalling from multimodal
single-cell data. IntraTalker infers cell-type-specific transcription factor activities and constructs receptomes that link receptors
to downstream target genes, which are then integrated with ligand-receptor predictions in CrossTalkeR. To prioritize signalling
receptors, the framework performs in silico receptor perturbation.

## Install

IntraTalker depends on two Bioconductor packages, which are not available from
CRAN and therefore have to be installed with `BiocManager`:

| Package | Used for |
| --- | --- |
| [ComplexHeatmap](https://bioconductor.org/packages/ComplexHeatmap) | Heatmaps of transcription factor activities |
| [scran](https://bioconductor.org/packages/scran) | Marker detection in the differential transcription factor analysis |

```r
# Install BiocManager if necessary
install.packages("BiocManager")

BiocManager::install(c("ComplexHeatmap", "scran"))
```

All remaining dependencies are on CRAN and are installed automatically.

Then install the package from GitHub with devtools:

```r
# Install devtools if necessary
install.packages("devtools")

# Install the package
devtools::install_github("CostaLab/IntraTalker")
```

or with the remotes package:

```r
# Install remotes if necessary
install.packages("remotes")

# Install the package
remotes::install_github("CostaLab/IntraTalker")
```

Both installers read the `biocViews` field in `DESCRIPTION` and add the
Bioconductor repositories automatically, so they can also resolve
`ComplexHeatmap` and `scran` on their own. Installing them up front with
`BiocManager` is the more reliable route, as it pins them to the Bioconductor
release matching your R version.
