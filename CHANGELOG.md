## SINCLAIR development version

- Improve error messages in SAMPLESHEET_CHECK processes. (#143, @kelly-sovacool)

## SINCLAIR 0.3.1

- Minor documentation improvements. (#135, @kelly-sovacool)
- SINCLAIR is now archived in Zenodo and can be cited with the DOI [10.5281/zenodo.15283503](https://doi.org/10.5281/zenodo.15283503). (#136, @kelly-sovacool)
- Fix bug where nextflow schema was not included in the installation. (#138, @kelly-sovacool)

## SINCLAIR 0.3.0

- **Breaking change**: GEX is now the default workflow. The `-entry` argument is no longer used. (#129, @kelly-sovacool)

### New features

- Allows users to determine what variables to regress out. (#55, @slsevilla)
- Overhaul the CLI to use python rather than bash, which introduces **breaking changes** (#61, @kelly-sovacool).
  - Create a script (`bin/sinclair`) to provide an interface to the CLI that works out-of-the-box without the need to install the python package with pip. (#80, @kelly-sovacool)
- Use `nextflow run -resume` by default, or turn it off with `sinclair run --forceall`. (#110, @kelly-sovacool)
- Add `--output` argument for `sinclair init` and `sinclair run`. (#110, @kelly-sovacool)
  - If not provided, commands are run in the current working directory.
  - This is equivalent to the nextflow `$launchDir` constant.
- Set the `publish_dir_mode` nextflow option to `link` by default. (#110, @kelly-sovacool)
- Set the `process.cache` nextflow option to `deep` by default rather than lenient on biowulf. (#110, @kelly-sovacool)
- Before launching the pipeline run:
  - The nextflow preview is printed. (#117, @kelly-sovacool)
  - The nextflow parameters are validated. (#127, @kelly-sovacool)

### Bug fixes

- Fix biowulf module syntax in `conf/modules.config`. (#81, @epehrsson)
- Fix filtering thresholds and use filtered object for downstream steps in `SEURAT_PROCESS`. (#81, @epehrsson)
- Fix seurat object and group assignment in `SEURAT_MERGE`. (#81, @epehrsson)
- Correctly map h5 files from cellranger with their input fastq files during preprocessing. (#82, @kelly-sovacool)
- Fix the $SLURM_JOB_ID variable name for biowulf. (#82, @kelly-sovacool)
- Fix file paths for test dataset. (#82, @kelly-sovacool)
- Use same number of PCs for merged object clustering as for integration. (#85, @epehrsson)
- Add LIGER UMAP to integration report. (#85, @epehrsson)
- Set all default parameters in `nextflow.config`. (#85, @epehrsson)
  - Previously, some parameters were set in `conf/process_params.config`, but we found this confusing, so we consolidated them to the main `nextflow.config` file.
- Allow sample IDs to contain hyphens. (#94, @wong-nw)
- Disable SCVI batch correction. (#109, @wong-nw)
  - This feature is on hold until a later release.
- LIGER now runs with 50 PCs by default instead of 20 (#109, @wong-nw)
- Output all R Markdown documents as HTML rather than PDF. (#112, @kelly-sovacool)
- Make sure values in the contrasts sheet are treated as strings. (#133, @kelly-sovacool)

### Documentation improvements

- The docs website now has a drop-down menu to switch between different versions. (#103, @kelly-sovacool)
- Fix broken image link. (#108, @wong-nw)
- Now using the readthedocs theme. (#111, @kelly-sovacool)
- Update to include quickstart and clarify existing documentation (#122, @wong-nw)

## SINCLAIR 0.2.0

### What's Changed

- feature: add precommmit
- fix: issue #49
- fix: issue #48
- doc: notes on wrapper
- feature: add stubruns
- fix: issue #43
- doc: add workflow image

## SINCLAIR 0.1.1

### What's Changed

- Feature/gex patch by @slsevilla in #35
- Feature/documentation by @slsevilla in #36
- Feature/documentation by @slsevilla in #37
- updated docs, gex_patch by @slsevilla in #38

## SINCLAIR 0.1.0

This is the first release of SINCLAIR 🎉

### What's Changed

- ci: sync source pipeline @slsevilla in #3
- feature: add R scripts by @abdallahamr in #8
- ci: auto add issues & PRs by @kelly-sovacool in #27
- feature: complete feature/gex pipeline by @slsevilla in #34

### New Contributors

- @abdallahamr made their first contribution in #8
- @kelly-sovacool made their first contribution in #27
