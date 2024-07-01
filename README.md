# Evaluation of Clinical Decision Support Tools in IBD

This repository includes the following files, which are needed to
reproduce the output of the [analysis
notebook](https://dahhamalsoud.github.io/analysis-notebooks/report-cdst-evaluation-leuven.html)
for the paper
`Calibration, Clinical Utility and Specificity of Clinical Decision Support Tools in Inflammatory Bowel Disease`:

- renv.lock, .Rprofile, renv/settings.json and renv/activate.R: These
  files are necessary to reproduce the same versions of the packages
  used in the original analysis.

- helpers.R: This script should be sourced at the beginning of the
  analysis to produce some customized functions.

- variable-selection-instability.R: This script reproduces the bootstrap
  analysis exploring the stability of variable selection in clinical
  prediction models.
