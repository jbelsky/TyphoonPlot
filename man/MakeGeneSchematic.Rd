% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_typhoon_plot_functions.R
\name{MakeGeneSchematic}
\alias{MakeGeneSchematic}
\title{Creates a schematic displaying gene locations in a given chromosomal region.}
\usage{
MakeGeneSchematic(feature_chr, feature_start, feature_end, cex_title = 1,
  bg_type = "white", fwd_gene_col = "gray", rev_gene_col = "gray",
  proteinCoding = T, geneName = T, omit_genes = NA)
}
\arguments{
\item{feature_chr}{Character, chromosome (e.g. "1")}

\item{feature_start}{Numeric, genomic start position}

\item{feature_end}{Numeric, genomic end position}

\item{cex_title}{Numeric, specifies the cex expansion factor for the title}

\item{bg_type}{Character, whether to make the background "white" or "transparent"
(Default: "white").}

\item{fwd_gene_col, rev_gene_col}{Character, specifies the color for both
the Watson (+) and Crick (-) genes (Default: "gray" for both).}

\item{proteinCoding}{Boolean, should only protein-coding genes be included as
opposed to including all open-reading-frames (Default: T)}

\item{geneName}{Boolean, should the gene name be displayed (Default: T)}

\item{omit_genes}{Character vector, should any gene names be excluded from the
plot.  This is useful if there are overlapping genes. (Default: NA)}
}
\value{
A schematic showing gene locations in a given chromosomal location
}
\description{
Each gene is designated as a gray box, with genes on the Watson (+) strand
displayed in the top row and genes on the Crick (-) strand displayed on
the bottom row.  Orientation of the gene name also denotes transcription
direction.  Currently this function will only make a gene schematic for
yeast genes in the sacCer2/SGD R61 genome version.
}
\examples{
 MakeGeneSchematic("1", 40000, 41000)
}
