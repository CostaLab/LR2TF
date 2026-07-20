# IntraTalker

This package enables condition- and cell-type-specific analysis of transcription factor activities and generates intracellular signaling Networks. The results can be combined with ligand-receptor predictions and analyzed differentially using CrossTalkeR. 

To perform receptor perturbation simulation, use the Python package [IntraTalkerpy](https://github.com/CostaLab/IntraTalkerpy)

## Abstract

Single-cell sequencing has advanced the study of cell-cell communication, yet most methods focus on intercellular ligand-
receptor interactions while neglecting downstream intracellular signalling cascades and the possibility that downstream
target genes themselves encode ligands, thereby propagating communication across multiple cells. We present In-
traTalker+CrossTalkeR that combines intracellular (IntraTalker) and intercellular (CrossTalkeR) signalling from multimodal
single-cell data. IntraTalker infers cell-type-specific transcription factor activities and constructs receptomes that link receptors
to downstream target genes, which are then integrated with ligand-receptor predictions in CrossTalkeR. To prioritize signalling
receptors, the framework performs in silico receptor perturbation.

## Install

Install the package from GitHub with devtools:

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
