% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/myGRoverlapFunctions.R
\name{overlap2GR}
\alias{overlap2GR}
\title{Get a dataframe from findOverlaps from a GRangesList with 2 items where colnames match the names of the items}
\usage{
overlap2GR(my2GRsToOverlap)
}
\arguments{
\item{my2GRsToOverlap}{a GRangesList with 2 items to overlap}
}
\value{
a dataframe with 2 columns with colnames identicals to names of the input with the indices of the overlaps
}
\description{
Get a dataframe from findOverlaps from a GRangesList with 2 items where colnames match the names of the items
}
\examples{
gr1 <- GenomicRanges::GRanges(
seqnames = "chr1",
ranges = IRanges::IRanges(start = c(1, 11, 199),
                         end = c(10, 12, 200)),
score = c(20, 30, 100))
gr2 <- GenomicRanges::GRanges(
seqnames = "chr1",
ranges = IRanges::IRanges(start = c(5, 101),
                         end = c(15, 102)),
score = c(1, 2))
overlap2GR(list(first = gr1, second = gr2))
}
