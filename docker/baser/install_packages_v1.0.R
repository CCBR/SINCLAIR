scRNA_handle_packages <- function(pkg_df) {
  for (rowid in rownames(pkg_df)) {
    pkg <- pkg_df[rowid, "package"]
    source <- pkg_df[rowid, "source"]
    version <- pkg_df[rowid, "version"]
    gh_name <- pkg_df[rowid, "gh_name"]

    need_install <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(need_install) != 0) {
      print(paste0("Installing: ", pkg))
      if (source == "bc") BiocManager::install(pkg)
      if (source == "cr") {
        install.packages(pkg,
          version = version, repos = "http://cran.us.r-project.org",
          local = FALSE
        )
      }
      if (source == "gh") remotes::install_github(gh_name, version = version, local = FALSE)
    }
  }
}

# read in package info
pkg_df <- read.csv("/opt2/Rpack_v1.0.config")

# for each package check installation, if present then load library
scRNA_handle_packages(pkg_df)
