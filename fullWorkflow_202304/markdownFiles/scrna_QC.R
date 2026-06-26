args <- commandArgs(trailingOnly = TRUE)

inputDir <- as.character(args[1])
outputDir <- as.character(args[2])
sample <- as.character(args[3])
resolution <- as.character(args[4])
species <- as.character(args[5])
annotRef <- as.character(args[6])

# sample=strsplit(dir,".Rdat")[[1]][1]

rmarkdown::render("workflow/scripts/scRNA_QC.Rmd",
  output_file = paste0("QC_Report_", sample, ".html"), output_dir = outputDir,
  params = list(
    sample = sample, resolution = resolution,
    imageDir = paste0("QC/images/", sample), ref = species, outFile = paste0("filtered/", sample, ".rds"),
    annot = annotRef
  )
)
