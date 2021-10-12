render_report = function(contrast, dataPath) {
  rmarkdown::render(
    "dea_table_template.Rmd", params = list(
      contrast = contrast,
      data = dataPath
    ),
    output_file = paste0(contrast, ".html")
  )
}

fit <- readRDS("fit_dea.Rds")

for (contrast_group in colnames(fit)[-1]) {
  render_report(contrast=contrast_group, data='fit_dea.Rds')
}

