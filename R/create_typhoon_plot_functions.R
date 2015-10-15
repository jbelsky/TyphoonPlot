# Ensure that the options(stringsAsFactors = FALSE)
old_options <- options(stringsAsFactors = FALSE)
on.exit(options(old_options), add = TRUE)

#' Get the TyphoonPlot matrix
#'
#' @param bam_file_name character, a bam file and its associated .bai index file
#' @param chr character, chromosome from the bam file
#' @param start_pos numeric, genomic start position
#' @param end_pos numeric, genomic end position
#' @return A frag_length x position matrix (250 rows x position columns) that
#'  can be visualized with DensDotPlot
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

#' Makes a heatmap of data within a matrix
#'
#' @param dot.m matrix, the 2-dimensional matrix to plot
#' @param z_min numeric, the minimum intensity of the plot (default: 0)
#' @param z_max numeric, the maximum intensity of the plot (default: 50)
#' @param low_col character, the minimum color (default: "white")
#' @param med_col character, the middle color (optional, default: "")
#' @param high_col character, the maximum color (default: "blue")
#' @param num_colors numeric, the number of colors between low_col and high_col
#' @param plot_title,x_label,y_label character (optional, default: "" for all)
#' @param x_axt,y_axt character, "s" (default) or "n"
#' @param plot_title_line numeric, (optional, default: NA)
#' @param plot_box boolean, (default: TRUE)
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

# Set up the gene
MakeGeneSchematic = function(feature_chr, feature_start, feature_end,
			       y_low = 0, y_high = 1, cex_title = 1, bg_type = "white",
             fwd_gene_col = "gray", rev_gene_col = "gray",
			       proteinCoding = T, geneName = T, omit_genes = NA, x_pos_title = 50,
			       gene_file_name = "/data/illumina_pipeline/scripts/feature_files/yeast/gene_tables/sacCer2_ucsc_sgdGeneTable.csv"
			      ){

	# Set up the plot
	plot(0, 0, type = "n", bty = "n", bg = bg_type,
	     xlim = c(feature_start, feature_end), xaxs = "i", xaxt = "n",
	     ylim = c(0, 1), yaxs = "i", yaxt = "n",
	     ann = F
	    )

	# Load the gene dataframe
	gene.df = read.csv(gene_file_name)

	# Subset only on protein coding if selected
	if(proteinCoding){

		idx = which(as.character(gene.df$name) != as.character(gene.df$sgd_name))

		gene.df = gene.df[idx,]

	}

	# Omit any gene if necessary
	if(any(!is.na(omit_genes))){

		gene.df = gene.df[-which(gene.df$sgd_name %in% omit_genes),]

	}

	# Convert to a GenomicRanges object
	gene.gr = GenomicRanges::GRanges(seqnames = gene.df$chr,
			  ranges = IRanges(start = gene.df$start, end = gene.df$end),
			  strand = gene.df$strand
			 )
	names(gene.gr) = gene.df$name

	# Create the feature gr
	feature.gr = GenomicRanges::GRanges(seqnames = feature_chr,
			     ranges = IRanges(start = feature_start, end = feature_end)
			    )

  # return(list(a = feature.gr, b = gene.gr))

	# Find the overlaps
  overlaps.hits = GenomicRanges::findOverlaps(feature.gr, gene.gr)

	if(any(subjectHits(overlaps.hits))){

    # Get the subjectHits
    subject_hits.v = subjectHits(overlaps.hits)

		# Enter in the genes
		for(i in 1:length(subject_hits.v)){
			plot_gene(gene.df[subject_hits.v[i],], y_low, y_high,
				  feature_start, feature_end, cex_title, geneName, x_pos_title, fwd_gene_col, rev_gene_col)
		}

	}

}

# Plot gene
plot_gene = function(gene.v, ylow, yhigh, xstart, xend, cex_title, geneName, x_pos_title = 50,
                     fwd_gene_color, rev_gene_color
                     ){

	# Get y_mid
	ymid = (yhigh + ylow) / 2

	# Add in the text
	if(gene.v$strand == "+"){

		# Make the rectangle
		rect(gene.v$start, ymid + 0.1, gene.v$end, yhigh - 0.1, col = fwd_gene_color)

		if(geneName){
			if(gene.v$start >= xstart){
				text(x = gene.v$start + x_pos_title, y = yhigh - 0.15, adj = c(0, 1),
				     labels = gene.v$sgd_name, font = 3, cex = cex_title)
			}else{
				text(x = gene.v$end - x_pos_title, y = yhigh - 0.15, adj = c(1, 1),
				     labels = gene.v$sgd_name, font = 3, cex = cex_title)
			}
		}
	}else{

		# Make the rectangle
		rect(gene.v$start, ylow + 0.1, gene.v$end, ymid - 0.1, col = rev_gene_color)

		if(geneName){
			if(gene.v$end <= xend){
				text(x = gene.v$end - x_pos_title, y = ylow + 0.15, adj = c(0, 1),
				     labels = gene.v$sgd_name, srt = 180, font = 3, cex = cex_title)
			}else{
				text(x = gene.v$start + x_pos_title, y = ylow + 0.15, adj = c(1, 1),
				     labels = gene.v$sgd_name, srt = 180, font = 3, cex = cex_title)
			}
		}
	}

}

# Set up the schematic section
set_chromatin_schematic = function(x_start = 0, x_end = 1, y_start = 0, y_end = 1){

	plot(0, 0, type = "n", bty = "n",
	     xlim = c(x_start, x_end), xaxs = "i", xaxt = "n",
	     ylim = c(y_start, y_end), yaxs = "i", yaxt = "n",
	     ann = F
	    )

}

# Enter in the axes
make_typhoon_plot_axes = function(x_start, x_end, x_step = 500, x_divide = 1000, x_format = "f", x_digits = 1,
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
