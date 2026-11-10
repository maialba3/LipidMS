# LipidMS 4.0.0000

## Major changes

- **Redesigned workflow based on xcms.** LipidMS now uses `XcmsExperiment` objects as its primary data structure, separating chromatographic data processing from lipid annotation.
- **Chromatographic processing delegated to xcms.** Peak detection, retention time alignment, and chromatographic peak grouping are performed using `xcms`, eliminating the need for a separate chromatographic processing workflow within LipidMS.
- **Unified annotation workflow.** The new `annotateLipids()` function integrates isotope annotation, MS1 candidate generation, adduct grouping, MS/MS spectral extraction, and rule-based lipid annotation.

## Lipid annotation

- Added MS1-based lipid candidate generation using theoretical lipid databases.
- Added grouping of coeluting peaks corresponding to different adducts.
- Integrated isotope annotation into the annotation workflow.
- Added rule-based matching of experimental MS/MS spectra against expected lipid fragments.
- Added support for lipid annotation using DDA, DIA, and SWATH acquisition data.
- Added lipid annotation at different levels of structural detail, depending on the available MS/MS evidence.

## Data structures and results

- Lipid annotations are stored within the input `XcmsExperiment` object.
- Peak-level annotations are incorporated into `chromPeakData`.
- Feature-level annotations can be added to `featureDefinitions` when chromatographic peaks have been grouped.
- Intermediate results, annotation parameters, fragmentation rules, and other relevant information are retained in `metadata(x)$LipidMS`.

## Lipid databases and fragmentation rules

- Added functionality for generating theoretical lipid databases using `createLipidDbs()`.
- Added configurable fragmentation rules for MS/MS-based lipid annotation.
- Added support for user-defined lipid databases, fragmentation rules, and adduct tables.

## Additional improvements

- Added support for parallel processing through `BiocParallel`.
- Improved integration with the Bioconductor mass spectrometry ecosystem.

## Breaking changes

- LipidMS 4.0 introduces a substantially redesigned workflow that is not directly interchangeable with previous versions.
- Chromatographic processing must now be performed with `xcms` before lipid annotation.
- The main annotation workflow expects an `XcmsExperiment` object rather than the input structures used in previous versions.
- Existing scripts based on earlier LipidMS workflows may require adaptation to the new API.
