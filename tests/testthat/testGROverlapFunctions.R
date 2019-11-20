context("GR overlap functions")
library(usefulLDfunctionsGR)
test_that("overlap2GR gives what is expected", {
  gr1 <- GenomicRanges::GRanges(seqnames = "chr1",
                                ranges = IRanges::IRanges(start = c(1, 11, 199),
                                                          end = c(10, 12, 200)),
                                score = c(20, 30, 100))
  gr2 <- GenomicRanges::GRanges(seqnames = "chr1",
                                ranges = IRanges::IRanges(start = c(5, 101),
                                                          end = c(15, 102)),
                                score = c(1, 2))
  expect_equal(overlap2GR(list(first = gr1, second = gr2)),
               data.frame(first = c(1, 2),
                          second = rep(1, 2)))
})

test_that("findAllOverlap gives what is expected", {
  gr1 <- GenomicRanges::GRanges(seqnames = "chr1",
                                ranges = IRanges::IRanges(start = c(1, 11, 199),
                                                          end = c(10, 12, 200)),
                                score = c(20, 30, 100))
  gr2 <- GenomicRanges::GRanges(seqnames = "chr1",
                                ranges = IRanges::IRanges(start = c(5, 101),
                                                          end = c(15, 102)),
                                score = c(1, 2))
  gr3 <- GenomicRanges::GRanges(seqnames = "chr1",
                                ranges = IRanges::IRanges(start = c(1, 100),
                                                          end = c(4, 120)),
                                score = c(10, 20))
  grL <- list(first = gr1, second = gr2, third = gr3)
  my.dataframe <- data.frame(first = c(1, 2, NA, 1, 1),
                             second = c(1, 1, 2, NA, 1),
                             third = c(1, NA, 2, 1, NA))
  expect_equal(findAllOverlap(grL), my.dataframe)
  all.dataframes <- list(first____second = data.frame(first = 1:2,
                                                      second = c(1, 1)),
                         first____third = data.frame(first = 1, third = 1),
                         second____third = data.frame(second = 2, third = 2),
                         first____second____third = my.dataframe)
  expect_equal(findAllOverlap(grL, returnAllMergedCalculated = T),
               all.dataframes)
})

test_that("filterByScore gives what is expected", {
  gr1 <- GenomicRanges::GRanges(seqnames = "chr1",
                                ranges = IRanges::IRanges(start = c(1, 11,
                                                                    199),
                                                          end = c(10, 12,
                                                                  200)),
                                score = c(20, 30, 100))
  gr2 <- GenomicRanges::GRanges(seqnames = "chr1",
                                ranges = IRanges::IRanges(start = c(5, 101),
                                                          end = c(15, 102)),
                                score = c(60, 90))
  gr3 <- GenomicRanges::GRanges(seqnames = "chr1",
                                ranges = IRanges::IRanges(start = c(1, 100),
                                                          end = c(4, 120)),
                                score = c(1, 2))
  grL <- list(first = gr1, second = gr2, third = gr3)
  my.dataframe <- data.frame(first = c(1, NA, 2, 3),
                             second = c(1:2, NA, NA),
                             third = c(1:2, NA, NA))
  expect_equal(filterByScore(grL, useNormScore = T), my.dataframe)
  all.merges.my.dataframe <- data.frame(first = c(1, 2, NA, 1, 1),
                                        second = c(1, 1, 2, NA, 1),
                                        third = c(1, NA, 2, 1, NA))
  expect_equal(filterByScore(grL, mAll = all.merges.my.dataframe,
                             useNormScore = T),
               my.dataframe)

  my.dataframe2 <- data.frame(first = c(1, NA, 2:3),
                              second = c(NA, 2:1, NA),
                              third = c(1:2, NA, NA))

  expect_equal(filterByScore(grL, useNormScore = F),
               my.dataframe2)

})

test_that("findOverlapAsClusters gives what is expected", {
  gr1 <- GenomicRanges::GRanges(seqnames = "chr1",
                                ranges = IRanges::IRanges(start = c(1, 11, 199),
                                                          end = c(10, 12, 200)),
                                score = c(20, 30, 100))
  gr2 <- GenomicRanges::GRanges(seqnames = "chr1",
                                ranges = IRanges::IRanges(start = c(5, 101),
                                                          end = c(15, 102)),
                                score = c(1, 2))
  gr3 <- GenomicRanges::GRanges(seqnames = "chr1",
                                ranges = IRanges::IRanges(start = c(1, 100),
                                                          end = c(4, 120)),
                                score = c(10, 20))
  grL <- list(first = gr1, second = gr2, third = gr3)
  my.dataframe <- data.frame(id = 1:3)
  my.dataframe$first <- list(c(1, 2), 3, numeric(0))
  my.dataframe$second <- list(1, numeric(0), 2)
  my.dataframe$third <- list(1, numeric(0), 2)
  my.dataframe$id <- NULL
  expect_equal(findOverlapAsClusters(grL), my.dataframe)
})
