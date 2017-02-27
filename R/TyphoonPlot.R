#' TyphoonPlot: A package for creating a TyphoonPlot
#'
#' The TyphoonPlot package provides helper processing and graphics functions to
#' display paired-end next-generation sequencing data in a helpful visualization
#' format
#'
#' @section Required Packages:
#' This package requires the following packages:
#' GenomeInfoDb (>= 1.2.4),
#' GenomicRanges (>= 1.18.4),
#' gplots (>= 2.16.0),
#' IRanges (>= 2.0.1),
#' Rsamtools (>= 1.18.3),
#' S4Vectors (>= 0.4.0),
#' splus2R (>= 1.2.0)
#'
#' @section TyphoonPlot Processing Functions:
#'  GetTyphoonPlotMat
#'
#' @section TyphoonPlot Graphing Functions:
#'  DensDotPlot
#'  MakeGeneSchematic
#'  MakeIndividualTyphoonPlot
#'  MakeTyphoonPlotAxes
#'  MergeCoveragePlot
#'  SetChromatinSchematic
#'  PlotNucleosome
#'  PlotNucleosomeModel
#'
#' @section TyphoonPlot Density Smoothing Functions:
#'  ConvertPairedReadBAMToGR
#'  GetDensityPeaks
#'  GetMNaseFeatureDensity
#'
#' @docType package
#' @name TyphoonPlot
NULL
