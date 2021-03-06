# Ensure that the options(stringsAsFactors = FALSE)
old_options <- options(stringsAsFactors = FALSE)
on.exit(options(old_options), add = TRUE)

#' Obtain the TyphoonPlot matrix needed for plotting.
#'
#' Given an indexed bam file and chromosomal genomic coordinates, this function
#' will return a matrix displaying the paired-end, sequenced reads in
#' two dimensions.  Note that this function will plot reads at half their
#' span to reduce the overlap between adjacent reads.  For example, a 100 bp
#' paired-end read spanning from position 201 to 300 will fill a 50 bp
#' x-coordinate region from position 226 to 275.
#'
#' @param bam_file_name Character, a bam file and its associated .bai index file
#' @param chr Character, chromosome from the bam file
#' @param start_pos Numeric, genomic start position
#' @param end_pos Numeric, genomic end position
#' @return A frag_length x position matrix (250 rows x position columns) that
#'    can be visualized with DensDotPlot
#' @examples
#' \dontrun{GetTyphoonPlotMat("/home/jab112/dm242.bam", "1", 10000, 15000)}
GetTyphoonPlotMat = function(bam_file_name, chr, start_pos, end_pos){

	# Load the bam_file
	chr.gr = ConvertPairedReadBAMToGR(bam_file_name, chr, start_pos, end_pos)

	# Subset on reads less than 250 bp
	chr.gr = chr.gr[BiocGenerics::width(chr.gr) <= 250]

	# Set the query GR
	query.gr = GenomicRanges::GRanges(seqnames = chr,
			                              ranges = IRanges::IRanges(start = start_pos, end = end_pos)
			                             )

	# Subset on the indices that overlap with the query.gr
	idx = S4Vectors::subjectHits(IRanges::findOverlaps(query.gr, chr.gr))
	chr.gr = chr.gr[idx]

	# Get a new typhoon_plot_gr that will contain the positions of the "half width" reads
	chr_tp.gr = chr.gr

  # Update the start and end coordinates of the "half width" reads
	IRanges::ranges(chr_tp.gr) = IRanges::IRanges(start = BiocGenerics::start(chr.gr) + round(BiocGenerics::width(chr.gr) / 4),
				                                        end = BiocGenerics::end(chr.gr) - round(BiocGenerics::width(chr.gr) / 4)
				                                       )

	# Set up the matrix
	mat.m = matrix(0, nrow = 250, ncol = end_pos - start_pos + 1)
	colnames(mat.m) = start_pos:end_pos

	# Create the sub feature
	mat.gr = GenomicRanges::GRanges(seqnames = chr,
                                  ranges = IRanges::IRanges(start = start_pos:end_pos, width = 1)
                                 )

	# Iterate through each fragment width
	for(i in 20:250){

		# Get the reads that have a particular fragment width
		idx = which(BiocGenerics::width(chr.gr) == i)

		if(any(idx)){

			# Count the overlaps with the mat.gr
			mat.m[i,] = IRanges::countOverlaps(mat.gr, chr_tp.gr[idx])

		}

	}

	# Return the mat.m
	return(mat.m)

}


GetVPlotMat = function(bam_file_name, chr, pos, win){

	# Set the start and end pos
	start_pos = pos - win
	end_pos = pos + win

	# Load the bam_file
	chr.gr = ConvertPairedReadBAMToGR(bam_file_name, chr, start_pos, end_pos)

	# Subset on reads less than 250 bp
	chr.gr = chr.gr[BiocGenerics::width(chr.gr) <= 250]

	# Set the query GR
	query.gr = GenomicRanges::GRanges(seqnames = chr,
			                  ranges = IRanges::IRanges(start = start_pos, end = end_pos)
			                 )

	# Subset on the indices that overlap with the query.gr
	idx = S4Vectors::subjectHits(IRanges::findOverlaps(query.gr, chr.gr))
	chr.gr = chr.gr[idx]

	# Get a new typhoon_plot_gr that will contain the positions of the "half width" reads
	chr_tp.gr = chr.gr

	# Update the start and end coordinates of the "half width" reads
	IRanges::ranges(chr_tp.gr) = IRanges::IRanges(start = BiocGenerics::start(chr.gr) + round(BiocGenerics::width(chr.gr) / 2),
				                      width = 1
				                     )

	# Set up the matrix
	mat.m = matrix(0, nrow = 250, ncol = end_pos - start_pos + 1)
	colnames(mat.m) = start_pos:end_pos

	# Create the sub feature
	mat.gr = GenomicRanges::GRanges(seqnames = chr,
					ranges = IRanges::IRanges(start = start_pos:end_pos, width = 1)
				       )

	# Iterate through each fragment width
	for(i in 20:250){

		# Get the reads that have a particular fragment width
		idx = which(BiocGenerics::width(chr.gr) == i)

		if(any(idx)){

			# Count the overlaps with the mat.gr
			mat.m[i,] = IRanges::countOverlaps(mat.gr, chr_tp.gr[idx])

		}

	}

	# Return the mat.m
	return(mat.m)

}

#' Obtain the aggregate VPlot matrix over a set of chromosomal coordinates (e.g. TF-binding sites)
#'
#' Given a data frame listing chromosomal coordinates, an indexed bam file, 
#' and a bp window, this function will return a matrix displaying the paired-end, 
#' sequenced reads in a standard V-plot (e.g. at the midpoint-fragment length
#' location in a matrix.
#'
#' @param data.df Data Frame, a feature data frame containing the "chr", "mid", and "strand" components for each location in the dataset
#' @param bam_file_name Character, a bam file and its associated .bai index file
#' @param win Numeric, number of bp to extend the V-plot from the midpoint location
#' @return A frag_length x position matrix (250 rows x (2*win + 1) columns) that
#'    can be visualized with DensDotPlot
#' @examples
#' \dontrun{GetAggregateVPlot(abf1_yeast_sacCer2.df, "/home/jab112/dm242.bam", 500)}
GetAggregateVPlot = function(data.df, bam_file_name, win){

    # Initialize the aggregate V-plot matrix
    mat.m = matrix(0, nrow = 250, ncol = (2*win) + 1)
    colnames(mat.m) = -win:win

    # Iterate through each feature
    for(i in 1:nrow(data.df)){

	# Get the chromosome coordinates and strand
	chr = as.character(data.df$chr[i])
	pos = data.df$mid[i]
	strand = as.character(data.df$strand[i])

	# Get the V-Plot matrix for the position
	cur_mat.m = GetVPlotMat(bam_file_name, chr, pos, win)

	# Flip the matrix if on negative strand
	if(strand == "-"){
	    cur_mat.m = cur_mat.m[,ncol(cur_mat.m):1]
	}

	# Append to mat.m
	mat.m = mat.m + cur_mat.m

    }

    # Return the mat.m
    return(mat.m)


}



#' Makes a heatmap of data within a matrix
#'
#' This function takes as input a TyphoonPlot matrix created from
#' \code{\link{GetTyphoonPlotMat}} and uses the \code{\link[graphics]{image}}
#' function to create the heatmap.
#'
#' @param dot.m Matrix, the 2-dimensional matrix to plot
#' @param z_min Numeric, the minimum intensity of the plot (default: 0)
#' @param z_max Numeric, the maximum intensity of the plot.  If not specified,
#'  by default the maximum will be set to the value at the 95th percentile (default: NA)
#' @param low_col Character, the minimum color (default: "white")
#' @param med_col Character, the middle color (optional, default: "")
#' @param high_col Character, the maximum color (default: "blue")
#' @param num_colors Numeric, the number of colors between low_col and high_col
#' @param plot_title,x_label,y_label Character (optional, default: "" for all)
#' @param x_axt,y_axt Character, "s" (default) or "n"
#' @param plot_title_line Numeric, (optional, default: NA)
#' @param plot_box Boolean, (default: TRUE)
#' @return A heatmap plot utilizing the image function
#' @examples
#' \dontrun{DensDotPlot(dot.m = typhoon.m, z_max = 50)}
DensDotPlot = function(dot.m, z_min = 0, z_max = NA,
			 low_col = "white", med_col = "", high_col = "blue", num_colors = 100,
		   plot_title = "", plot_title_line = NA,
       x_label = "", y_label = "",
			 x_axt = "s", y_axt = "s",
			 plot_box = TRUE){

  # Check that dot.m is a matrix
  if(!is.matrix(dot.m)){
    stop("dot.m is not a matrix, DensDotPlot requires a matrix!")
  }

  if(is.na(z_max)){
    # Make the z_max equivalent to the 95th percentile if z_max is not specified
    z_max = quantile(as.vector(dot.m), probs = 0.95)

  }

	# For points that are either above or below z_max or z_min respectively, set them to
	# the z_min and z_max (otherwise, plot shows arbitrary colors
  dot.m[which(dot.m >= z_max)] = z_max
  dot.m[which(dot.m <= z_min)] = z_min

	# Get the current column names
	if(!is.null(colnames(dot.m))){

    position.v = as.numeric(colnames(dot.m))

	}else{

    position.v = 1:ncol(dot.m)

	}

  # Get the row values
  row.v = 1:nrow(dot.m)

  # Make the colorpanel
	if(med_col == ""){
	  make_colorpanel = gplots::colorpanel(n = num_colors, low = low_col, high = high_col)
  }else{
    make_colorpanel = gplots::colorpanel(n = num_colors, low = low_col, med = med_col, high = high_col)
	}

	# Make the heatmap utilizing the parameters specified above
  image(position.v, row.v, t(dot.m), col = make_colorpanel, zlim = c(z_min, z_max),
        xlab = x_label, ylab = y_label, xaxt = x_axt, yaxt = y_axt, bty = "n"
       )

	# Set the title
	title(main = plot_title, line = plot_title_line)

	# Add a box around the plot
	if(plot_box){
        	box(which = "plot", lty = "solid")
	}

}

#' Creates a schematic displaying gene locations in a given chromosomal region.
#'
#' Each gene is designated as a gray box, with genes on the Watson (+) strand
#' displayed in the top row and genes on the Crick (-) strand displayed on
#' the bottom row.  Orientation of the gene name also denotes transcription
#' direction.  Currently this function will only make a gene schematic for
#' yeast genes in the sacCer2/SGD R61 genome version.
#'
#' @param feature_chr Character, chromosome (e.g. "1")
#' @param feature_start Numeric, genomic start position
#' @param feature_end Numeric, genomic end position
#' @param cex_title Numeric, specifies the cex expansion factor for the title
#' @param bg_type Character, whether to make the background "white" or "transparent"
#'  (Default: "white").
#' @param fwd_gene_col,rev_gene_col Character, specifies the color for both
#'  the Watson (+) and Crick (-) genes (Default: "gray" for both).
#' @param proteinCoding Boolean, should only protein-coding genes be included as
#'  opposed to including all open-reading-frames (Default: T)
#' @param geneName Boolean, should the gene name be displayed (Default: T)
#' @param omit_genes Character vector, should any gene names be excluded from the
#'  plot.  This is useful if there are overlapping genes. (Default: NA)
#'
#' @return A schematic showing gene locations in a given chromosomal location
#'
#' @examples
#'  MakeGeneSchematic("1", 40000, 41000)
MakeGeneSchematic = function(feature_chr, feature_start, feature_end,
			       cex_title = 1, bg_type = "white",
             fwd_gene_col = "gray", rev_gene_col = "gray",
			       proteinCoding = T, geneName = T, omit_genes = NA
			       ){

	# Set up the plot
	plot(0, 0, type = "n", bty = "n", bg = bg_type,
	     xlim = c(feature_start, feature_end), xaxs = "i", xaxt = "n",
	     ylim = c(0, 1), yaxs = "i", yaxt = "n",
	     ann = F
	    )

	# Subset only on protein coding if selected
	if(proteinCoding){

		idx = which(as.character(yeast_gene.df$name) != as.character(yeast_gene.df$sgd_name))

		yeast_gene.df = yeast_gene.df[idx,]

	}

	# Omit any gene if necessary
	if(any(!is.na(omit_genes))){

		yeast_gene.df = yeast_gene.df[-which(yeast_gene.df$sgd_name %in% omit_genes),]

	}

	# Convert to a GenomicRanges object
	gene.gr = GenomicRanges::GRanges(seqnames = yeast_gene.df$chr,
			                             ranges = IRanges::IRanges(start = yeast_gene.df$start, end = yeast_gene.df$end),
			                             strand = yeast_gene.df$strand
			                            )
	names(gene.gr) = yeast_gene.df$name

	# Create the feature gr
	feature.gr = GenomicRanges::GRanges(seqnames = feature_chr,
			                                ranges = IRanges::IRanges(start = feature_start, end = feature_end)
			                               )

  # Find the overlaps
  overlaps.hits = GenomicRanges::findOverlaps(feature.gr, gene.gr)

	if(any(S4Vectors::subjectHits(overlaps.hits))){

    # Get the subjectHits
    subject_hits.v = S4Vectors::subjectHits(overlaps.hits)

		# Enter in the genes
		for(i in 1:length(subject_hits.v)){
			PlotGene(yeast_gene.df[subject_hits.v[i],], y_low = 0, y_high = 1,
				  feature_start, feature_end, cex_title, geneName, x_pos_title = 50, fwd_gene_col, rev_gene_col)
		}

	}

}

#' Helper function to plot individual genes as rectangles.
#'
#' Helper function for \code{\link{MakeGeneSchematic}} to plot individual
#' genes.  Each gene is designated as a gray box, with genes on the Watson (+)
#' strand displayed in the top row and genes on the Crick (-) strand displayed
#' on the bottom row.  Orientation of the gene name also denotes transcription
#' direction.
#'
#' @param gene.v Vector, selected row from \code{yeast_gene.df} feature file
#'  containing information about the gene.
#' @param y_low Numeric, the bottom \emph{y} rectangular coordinate.
#' @param y_high Numeric, the top \emph{y} rectangular coordinate.
#' @param x_start Numeric, the left \emph{x} rectangular coordinate of the plot.
#' @param x_end Numeric, the right \emph{x} rectangular coordinate of the plot.
#' @param x_pos_title Numeric, the number of inset bp specifying the gene name.
#' @inheritParams MakeGeneSchematic
#'
#' @return Rectangle schematic depicting gene coordinates
#'
#' @examples
#'  PlotGene(gene.v = yeast_gene.df[1,], y_low = 0.5, y_high = 1,
#'           x_start = 40000, x_end = 41000, cex_title = 1,
#'           geneName = T, x_pos_title = 50,
#'           fwd_gene_color = "gray", rev_gene_color = "gray"
#'          )
PlotGene = function(gene.v, y_low, y_high, x_start, x_end,
                    cex_title, geneName, x_pos_title = 50,
                    fwd_gene_color, rev_gene_color
                   ){

	# Get y_mid
	y_mid = (y_high + y_low) / 2

	# Add in the text
	if(gene.v$strand == "+"){

		# Make the rectangle
		rect(gene.v$start, y_mid + 0.1, gene.v$end, y_high - 0.1, col = fwd_gene_color)

		if(geneName){
			if(gene.v$start >= x_start){
				text(x = gene.v$start + x_pos_title, y = y_high - 0.15, adj = c(0, 1),
				     labels = gene.v$sgd_name, font = 3, cex = cex_title)
			}else{
				text(x = gene.v$end - x_pos_title, y = y_high - 0.15, adj = c(1, 1),
				     labels = gene.v$sgd_name, font = 3, cex = cex_title)
			}
		}
	}else{

		# Make the rectangle
		rect(gene.v$start, y_low + 0.1, gene.v$end, y_mid - 0.1, col = rev_gene_color)

		if(geneName){
			if(gene.v$end <= xend){
				text(x = gene.v$end - x_pos_title, y = y_low + 0.15, adj = c(0, 1),
				     labels = gene.v$sgd_name, srt = 180, font = 3, cex = cex_title)
			}else{
				text(x = gene.v$start + x_pos_title, y = y_low + 0.15, adj = c(1, 1),
				     labels = gene.v$sgd_name, srt = 180, font = 3, cex = cex_title)
			}
		}
	}

}

# Set up the schematic section
SetChromatinSchematic = function(x_start = 0, x_end = 1, y_start = 0, y_end = 1){

	plot(0, 0, type = "n", bty = "n",
	     xlim = c(x_start, x_end), xaxs = "i", xaxt = "n",
	     ylim = c(y_start, y_end), yaxs = "i", yaxt = "n",
	     ann = F
	    )

}

# Enter in the axes
MakeTyphoonPlotAxes = function(x_start, x_end, x_step = 500, x_divide = 1000, x_format = "f", x_digits = 1,
								               addXAxis = T, addXLabel = T, addXTick = T, x_axis_line = NA, x_axis_pos = NA,
								               addYAxis = T, addYLabel = T
								              ){

	# Get the axes to plot
	x_axes = seq((x_start - (x_start %% x_step)), x_end + x_step, x_step)

	# Make the WT typhoon plot
	if(addXAxis){

		# Check if the label should be added
		x_label = F
		if(addXLabel){
			x_label = formatC(x_axes / x_divide, digits = x_digits, format = x_format, big.mark = ",")
		}

		# Modify the start and end axes by 0.5
		x_axes[1] = x_axes[1] - 0.5
		x_axes[length(x_axes)] = x_axes[length(x_axes)] + 0.5

		axis(1, at = x_axes, labels = x_label,
		     line = x_axis_line, pos = x_axis_pos, tick = addXTick
		    )
		axis(1, at = x_axes - x_step/2, labels = F,
		     line = x_axis_line, pos = x_axis_pos,
			 tick = addXTick, tcl = 0.625 * par()$tcl
		    )
	}

	if(addYAxis){

		# Check if the label should be added
		y_label = F
		if(addYLabel){
			y_label = c(0, 100, 200)
		}

		axis(2, at = c(0.5, 100, 200), labels = y_label)
		axis(2, at = c(50, 150, 250.5), labels = F, tcl = 0.625 * par()$tcl)

	}

}




MakeIndividualTyphoonPlot = function(plot_file_name, bam_file_name, chr, start_pos, end_pos){

  # Set up the plot
  png(file = plot_file_name, width = 10, height = 7, units = "in", res = 300)

  # Close the screens
  close.screen(all.screens = T)

  # Set up the plot screening matrix
  plot_scr.m = matrix(c(0, 1, 0.9, 1,
                        0, 1, 0.8, 0.9,
                        0, 1, 0, 0.8
                       ), ncol = 4, byrow = T
                     )

  # Split the screen
  plot_scr.s = split.screen(plot_scr.m)

  # Make the GeneSchematic
  screen(plot_scr.s[1])
  par(mar = c(0, 4.5, 0, 2))
  MakeGeneSchematic(chr, start_pos, end_pos)

  # Make the Nucleosome Schematic
  screen(plot_scr.s[2])
  par(mar = c(0, 4.5, 0, 2))
  PlotNucleosomeModel(bam_file_name, chr, start_pos, end_pos)

  # Create the DensDotPlot
  screen(plot_scr.s[3])
  par(mar = c(4.5, 4.5, 0, 2))

  # Make the TyphoonPlot
  DensDotPlot(GetTyphoonPlotMat(bam_file_name, chr, start_pos, end_pos),
              x_axt = "n", y_axt = "n",
              y_lab = "Fragment length (bp)"
             )
  title(xlab = paste("Chr ", as.roman(chr), " position (kb)", sep = ""), family = "serif", line = 2.5)
  MakeTyphoonPlotAxes(x_start = start_pos, x_end = end_pos)

  # Close the device
  close.screen(all.screens = T)
  dev.off()

}


































MergeCoveragePlot = function(g1.m, g2.m, z_max = 10, plot_bty = "o"){

  # Get the flank
  flank = (ncol(g1.m) - 1)/2

  # Get the frag_high
  frag_high = nrow(g1.m)

  g1.v = as.vector(t(g1.m))
  g2.v = as.vector(t(g2.m))

  g1.v[which(g1.v > z_max)] = z_max
  g2.v[which(g2.v > z_max)] = z_max
  g1.v[which(g1.v < 0)] = 0
  g2.v[which(g2.v < 0)] = 0

  # Create the coordinates for the graph
  x_left = rep((-flank:flank) - 0.5, frag_high)
  y_low = rep((1:frag_high) - 0.5, each = (2 * flank + 1))

  # Create the plot
  plot(0, 0, type = "n", bty = plot_bty,
       xlim = c(-flank, flank), xaxs = "i", xaxt = "n",
       ylim = c(0.5, frag_high + 0.5), yaxs = "i", yaxt = "n",
       ann = F
  )

  # Add in the colors
  rect(x_left, y_low, x_left + 1, y_low + 1,
       col = rgb(g1.v, g2.v, 0, z_max, maxColorValue = z_max), border = NA
      )

}
