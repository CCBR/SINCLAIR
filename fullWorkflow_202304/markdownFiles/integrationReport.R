args <- commandArgs(trailingOnly = TRUE)

mergedRDS <- as.character(args[1])
integratedRDS <- as.character(args[2])
harmonyGroupRDS <- as.character(args[3])
harmonySampleRDS <- as.character(args[4])
rpcaRDS <- as.character(args[5])
imageDir <- as.character(args[6])
outputDir <- as.character(args[7])
resolution <- as.character(args[8]) # comma-delimited string of resolutions
contrasts <- as.character(args[9]) #
groups <- as.character(args[10]) # boolean: are groups included
citeseq <- as.character(args[11]) # boolean: are citeseq results included

print(args)
# stop()


if (groups == "YES") {
  groups <- TRUE
} else {
  groups <- FALSE
}
if (citeseq == "YES") {
  citeseq <- TRUE
} else {
  citeseq <- FALSE
}

rmarkdown::render("workflow/scripts/integrationReports.Rmd",
  output_file = paste0("QC_Report_", contrasts, ".html"), output_dir = outputDir,
  params = list(
    mergedObj = mergedRDS, integratedObj = integratedRDS, rpcaObj = rpcaRDS, harmonyGroupObj = harmonyGroupRDS, harmonySampleObj = harmonySampleRDS,
    imageDir = imageDir, resolution = resolution, contrasts = contrasts,
    groups = groups, citeseq = citeseq
  )
)
