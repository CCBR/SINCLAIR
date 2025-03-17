## SINCLAIR development version

- Overhaul the CLI to use python rather than bash, which introduces breaking changes (#61, @kelly-sovacool).
  - Create a script (`bin/sinclair`) to provide an interface to the CLI that works out-of-the-box without the need to install the python package with pip. (#80, @kelly-sovacool)
- Bug fixes
  - Fix biowulf module syntax in `conf/modules.config`. (#81, @epehrsson)
  - Fix filtering thresholds and use filtered object for downstream steps in `SEURAT_PROCESS`. (#81, @epehrsson)
  - Fix seurat object and group assignment in `SEURAT_MERGE`. (#81, @epehrsson)
  - Correctly map h5 files from cellranger with their input fastq files during preprocessing. (#82, @kelly-sovacool)
  - Fix the $SLURM_JOB_ID variable name for biowulf. (#82, @kelly-sovacool)
  - Fix file paths for test dataset. (#82, @kelly-sovacool)
  - Add the tex module for preprocess and merge, which produce PDF files from R Markdown. (#82, @kelly-sovacool)
  - Use same number of PCs for merged object clustering as for integration. (#85, @epehrsson)
  - Add LIGER UMAP to integration report. (#85, @epehrsson)
  - Set all default parameters in `nextflow.config`. (#85, @epehrsson)
    - Previously, some parameters were set in `conf/process_params.config`, but we found this confusing, so we consolidated them to the main `nextflow.config` file.
  - Allow sample IDs to contain hyphens. (#94, @wong-nw)
  - Places SCVI batch correction on hold until a later release (#69, @wong-nw)
- New features
  - Allows users to determine what variables to regress out. (#55, @slsevilla)
- Documentation improvements
  - The docs website now has a drop-down menu to switch between different versions. (#103, @kelly-sovacool)

## SINCLAIR v0.2.0

### What's Changed

- feature: add precommmit
- fix: issue #49
- fix: issue #48
- doc: notes on wrapper
- feature: add stubruns
- fix: issue #43
- doc: add workflow image

## SINCLAIR v0.1.1

### What's Changed

- Feature/gex patch by @slsevilla in #35
- Feature/documentation by @slsevilla in #36
- Feature/documentation by @slsevilla in #37
- updated docs, gex_patch by @slsevilla in #38

## SINCLAIR v0.1.0

This is the first release of SINCLAIR 🎉

### What's Changed

- ci: sync source pipeline @slsevilla in #3
- feature: add R scripts by @abdallahamr in #8
- ci: auto add issues & PRs by @kelly-sovacool in #27
- feature: complete feature/gex pipeline by @slsevilla in #34

### New Contributors

- @abdallahamr made their first contribution in #8
- @kelly-sovacool made their first contribution in #27
