# LipidMS 4.0

**Lipid annotation from LC-MS data using an xcms-based workflow**

LipidMS is an R package for the annotation of lipids from liquid chromatography–mass spectrometry (LC-MS/MS) data, supporting data-dependent acquisition (DDA) and data-independent acquisition (DIA/SWATH) workflows.

**Version 4.0 introduces a major redesign of the LipidMS workflow**, separating chromatographic data processing from lipid annotation. Rather than implementing its own peak detection and data processing algorithms, LipidMS now builds on the established [xcms](https://bioconductor.org/packages/xcms/) framework for LC-MS data processing, allowing the package to focus on lipid-specific annotation strategies.

This integration provides greater flexibility in data processing, improves interoperability with the Bioconductor ecosystem, and facilitates the development of reproducible lipidomics workflows.

## Installation

LipidMS 4.0 can be installed directly from GitHub. We recommend using an up-to-date version of R and installing the required Bioconductor dependencies before installing LipidMS.

```r
# Install BiocManager and remotes if needed
install.packages(c("BiocManager", "remotes"))

# Install Bioconductor dependencies
BiocManager::install(c("xcms", "MsExperiment", "Spectra", "BiocParallel", "CAMERA"))

# Install LipidMS from GitHub
remotes::install_github("maialba3/LipidMS")
```

Once installed, load the package:

```r
library(LipidMS)
```

To launch the interactive Shiny application (still in development):

```r
LipidMSapp()
```

## What's new in LipidMS 4.0?

### Native integration with xcms

The most significant change in LipidMS 4.0 is the adoption of `XcmsExperiment` as the central data structure.

Chromatographic peak detection, retention time alignment, and feature grouping are performed using xcms, while LipidMS provides the downstream lipid annotation workflow.

This architecture offers several advantages:

- **Flexible data processing:** users can select among the peak detection, alignment, and feature grouping algorithms available in xcms, adapting the workflow to different instruments, chromatographic methods, and experimental designs.
- **Reduced code redundancy:** by relying on an established LC-MS processing framework, LipidMS avoids duplicating chromatographic processing functionality.
- **Improved interoperability:** annotated data remain compatible with xcms and other Bioconductor packages, facilitating downstream statistical analysis and integration with other omics datasets.
- **Reproducibility:** chromatographic processing results, lipid annotations, and associated metadata can be retained within a common data structure.
- **Easier maintenance and extensibility:** separating data processing from annotation allows lipid-specific methods to evolve independently of chromatographic processing algorithms.

### Redesigned lipid annotation workflow

LipidMS 4.0 introduces a modular annotation strategy that combines accurate-mass matching, adduct relationships, isotope information, and MS/MS fragmentation evidence.

The main steps include:

1. **Isotope annotation:** identification of isotopic patterns to distinguish monoisotopic peaks from isotope-related signals.
2. **MS1 candidate generation:** matching experimental precursor masses against theoretical lipid databases, considering lipid classes, molecular compositions, and adducts.
3. **Adduct grouping:** identification of coeluting signals potentially corresponding to different adducts of the same lipid.
4. **MS/MS data extraction:** extraction of fragmentation information associated with lipid candidates, using acquisition-specific strategies for DDA and DIA/SWATH data.
5. **Rule-based lipid annotation:** evaluation of class-specific and fatty acyl chain-specific fragmentation evidence, including fragmentation rules and relative fragment intensities.
6. **Annotation reporting:** incorporation of lipid annotations and supporting information into the original `XcmsExperiment` object.

### Support for different acquisition strategies

LipidMS 4.0 supports lipid annotation from DDA and DIA/SWATH acquisitions, including datasets containing samples acquired using different acquisition strategies.

MS1 information is used to generate lipid candidates, while MS/MS evidence is extracted from the relevant acquisitions to support lipid identification.

### Integrated annotation results

LipidMS stores annotation results within the `XcmsExperiment` object, preserving the relationship between chromatographic signals and their lipid assignments.

Annotation information is organized at different levels:

- **Chromatographic peaks:** lipid assignments, adduct information, and annotation confidence levels are stored in `chromPeakData()`.
- **Grouped features:** when feature grouping has been performed, lipid annotation summaries are added to `featureDefinitions()`.
- **Annotation metadata:** intermediate results, MS/MS evidence, annotation parameters, and fragmentation rules are stored in `metadata(x)$LipidMS`.

This approach facilitates downstream analysis without requiring separate annotation result objects.

## LipidMS 4.0 workflow

The analysis is divided into two complementary stages.

**1. LC-MS data processing — xcms**

Raw LC-MS data are processed using xcms to perform chromatographic peak detection and, when appropriate, retention time alignment and feature grouping.

**2. Lipid annotation — LipidMS**

The resulting `XcmsExperiment` object is used as input for LipidMS, which performs isotope annotation, MS1 candidate generation, adduct grouping, MS/MS evidence extraction, and rule-based lipid annotation.

The annotated object can subsequently be used for lipid-level data exploration, feature quantification, statistical analysis, and downstream integration.

## Main annotation function

The lipid annotation workflow is implemented through:

`annotateLipids()`

This function takes an `XcmsExperiment` object and applies the LipidMS annotation workflow using configurable annotation parameters, lipid databases, fragmentation rules, and adduct definitions.

LipidMS returns the updated `XcmsExperiment` object, retaining the original chromatographic processing results together with the lipid annotations.

## Why xcms?

xcms is a widely used, open-source framework for mass spectrometry data processing within Bioconductor.

By building on xcms, LipidMS 4.0 combines an established and flexible chromatographic processing framework with lipid-specific annotation algorithms.

**The goal is not to replace xcms, but to complement it:** xcms provides the chromatographic features, and LipidMS provides their lipid-specific interpretation.

This separation makes LipidMS more adaptable to different LC-MS platforms and experimental workflows while concentrating development efforts on lipid annotation.

