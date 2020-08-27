data("sample_matrix", package = "xts")
sample_xts <- xts::xts(sample_matrix[,], order.by = as.Date(rownames(sample_matrix)))
usethis::use_data(sample_xts, overwrite = T)

