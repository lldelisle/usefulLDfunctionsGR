context("Rtracklayer Basic functions")
library(usefulLDfunctionsGR)
test_that("getTSSinUCSCFormatFromEnsemblGTF gives what is expected", {
  tests_dir <- system.file("tests", package = "usefulLDfunctionsGR")
  test_gtf <- file.path(tests_dir,
                        "Homo_sapiens.GRCh38.95_491firstLines.gtf.gz")
  load(file.path(tests_dir, "tss.Rdata"))
  expect_equal(getTSSinUCSCFormatFromEnsemblGTF(test_gtf), expectedResult)
})
