# Name of functions contained in this functions file
# FUNCTION NAME			FUNCTION DEFINITION
# get_typhoon_plot_mat  	function(chr, start_pos, end_pos, bam_file_name)
# read_in_paired_bam 		function(bam_file_name, chr, type = "GR")
# dens_dot_plot 		function(dot.m, z_min = 0, z_max = 100, ...)
# make_gene_schematic 		function(feature_chr, feature_start, feature_end, ...)
# plot_gene 			function(gene.v, ylow, yhigh, xstart, xend, cex_title, geneName, x_pos_title = 50){
# set_chromatin_schematic	function(x_start = 0, x_end = 1, y_start = 0, y_end = 1){

# Load the libraries
# library(GenomicRanges)
# library(Rsamtools)
# library(gplots)

# Ensure that the options(stringsAsFactors = FALSE)
old_options <- options(stringsAsFactors = FALSE)
on.exit(options(old_options), add = TRUE)

# Get the process_bam_functions
# source("/data/illumina_pipeline/scripts/process_bam_scripts/process_bam_functions.R")

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
	chr.gr = convert_paired_read_bam_to_gr(bam_file_name, chr, start_pos, end_pos)

	# Subset on reads less than 250 bp
	chr.gr = chr.gr[width(chr.gr) <= 250]

	# Set the query GR
	query.gr = GRanges(seqnames = chr,
			   ranges = IRanges(start = start_pos, end = end_pos)
			  )

	# Subset on the indices that overlap with the query.gr
	idx = subjectHits(findOverlaps(query.gr, chr.gr))
	chr.gr = chr.gr[idx]

	# Get a new typhoon_plot_gr that will contain the positions of the "half width" reads
	chr_tp.gr = chr.gr

	# Update the start and end coordinates of the "half width" reads
	ranges(chr_tp.gr) = IRanges(start = start(chr.gr) + round(width(chr.gr) / 4),
				    end = end(chr.gr) - round(width(chr.gr) / 4)
				   )

	# Set up the matrix
	mat.m = matrix(0, nrow = 250, ncol = end_pos - start_pos + 1)
	colnames(mat.m) = start_pos:end_pos

	# Create the sub feature
	mat.gr = GRanges(seqnames = chr, ranges = IRanges(start = start_pos:end_pos, width = 1))

	# Iterate through each fragment width
	for(i in 20:250){

		# Get the reads that have a particular fragment width
		idx = which(width(chr.gr) == i)

		if(any(idx)){

			# Count the overlaps with the mat.gr
			mat.m[i,] = countOverlaps(mat.gr, chr_tp.gr[idx])

		}

	}

	# Scale the mat to 20,000,000 reads
	mat.m = mat.m / (get_total_read_number(bam_file_name, as.character(1:16)) / 20E6)

	# Return the mat.m
	return(mat.m)

}

get_typhoon_plot_mat_java = function(bam_file_name, chr, start, end){

	# Set the jar_file
	jar_file = "/data/illumina_pipeline/scripts/java_scripts/jar_files/get-typhoon-plot-R.jar"

	# Create the temporary file
	temp_file = paste(".", paste(chr, start, end, sep = "_"), sep = "")

	cat("The temp file is\t", temp_file, "\n")

	# Create the output string
	output_str = paste("java -jar", jar_file, bam_file_name, chr, start, end, temp_file)

	# Run the command
	system(command = output_str)

	# Read in the temp_file
	cat("Reading in\t", temp_file, "...\n")
	mat.df = read.csv(temp_file, header = T, row.names = 1, check.names = F)

	# Remove the temporary file
	file.remove(temp_file)

	# Convert into matrix
	mat.m = as.matrix(mat.df)

	# Return the mnase_dens.v
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
DensDotPlot = function(dot.m, z_min = 0, z_max = 100,
			 low_col = "white", med_col = "", high_col = "blue", num_colors = 100,
		   plot_title = "", plot_title_line = NA,
       x_label = "", y_label = "",
			 x_axt = "s", y_axt = "s",
			 plot_box = TRUE){

  # Check that dot.m is a matrix
  if(!is.matrix(dot.m)){
    stop("dot.m is not a matrix, DensDotPlot requires a matrix!")
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
make_gene_schematic = function(feature_chr, feature_start, feature_end,
			       y_low = 0, y_high = 1, cex_title = 1, bg_type = "white",
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
	gene.gr = GRanges(seqnames = gene.df$chr,
			  ranges = IRanges(start = gene.df$start, end = gene.df$end),
			  strand = gene.df$strand
			 )
	names(gene.gr) = gene.df$name

	# Create the feature gr
	feature.gr = GRanges(seqnames = feature_chr,
			     ranges = IRanges(start = feature_start, end = feature_end)
			    )

	# Find the overlaps
	overlaps.df = as.data.frame(as.matrix(findOverlaps(feature.gr, gene.gr)))

	if(any(nrow(overlaps.df))){

		# Enter in the genes
		for(i in 1:nrow(overlaps.df)){
			plot_gene(gene.df[overlaps.df$subjectHits[i],], y_low, y_high,
				  feature_start, feature_end, cex_title, geneName, x_pos_title)
		}

	}

}

# Plot gene
plot_gene = function(gene.v, ylow, yhigh, xstart, xend, cex_title, geneName, x_pos_title = 50){

	# Get y_mid
	ymid = (yhigh + ylow) / 2

	# Add in the text
	if(gene.v$strand == "+"){

		# Make the rectangle
		rect(gene.v$start, ymid + 0.1, gene.v$end, yhigh - 0.1, col = "gray")

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
		rect(gene.v$start, ylow + 0.1, gene.v$end, ymid - 0.1, col = "gray")

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





get_mnase_feature_density = function(bam_file_name, chr, start, end, fragL, fragH, bw){

	# Set the jar_file
	jar_file = "/data/illumina_pipeline/scripts/java_scripts/jar_files/get-mnase-density-R.jar"

	# Create the output string
	output_str = paste("java -jar", jar_file, bam_file_name, chr, start, end, fragL, fragH, bw)

	# Run the command
	mnase_dens.v = as.numeric(system(command = output_str, intern = T))

	# Return the mnase_dens.v
	return(mnase_dens.v)

}



























# Create the output typhoon plot
create_output_plot = function(plot_file_name, typhoon1.m, typhoon2.m,
			      typhoon1_title, typhoon2_title,
			      chr, start_coord, end_coord,
			      plot_width = 10, plot_height = 8, z_max = 5
			     ){


	# Open the plotting window
	png(file = plot_file_name, width = plot_width, height = plot_height, units = "in", res = 300)

	# Set up the screen
	close.screen(all.screens = T)

	# Each "screen" is like a mini plot
	# The syntax is left_coordinate,right_coordinate,bottom_coordinate,top_coordinate, and
	# each coordinate is scaled from 0 to 1
	scr.m = matrix(c(0.05, 1, 0.9, 1,
			 0.05, 1, 0.5, 0.9,
			 0.05, 1, 0.1, 0.5,
			 0.05, 1, 0, 0.1,
			 0, 0.05, 0.1, 0.9
			), ncol = 4, byrow = T
		      )

	# Split the screen
	# Each "screen" is given an index corresponding to the row in the scr.m matrix
	scr.s = split.screen(scr.m)

	# Set the left and right par
	left_par = 2
	right_par = 2

	# Make the gene schematic
	screen(scr.s[1])
	par(mar = c(0, left_par, 0, right_par))

	# The "set_chromatin_schematic() function just creates a blank plot
	set_chromatin_schematic()

	# The make_gene_schematic function pulls the genes from the following file:
	# /home/jab112/feature_files/yeast/gene_tables/sacCer2_ucsc_sgdGeneTable.csv
	make_gene_schematic(chr, start_coord, end_coord, y_low = 0, y_high = 1, proteinCoding = F)

	# Make the first typhoon plot
	screen(scr.s[2])
	par(mar = c(3, left_par, 3, right_par))
	dens_dot_plot(typhoon1.m, z_max = z_max, x_label = "", y_axt = "n", x_axt = "n", plot_title = typhoon1_title)
	abline(v = 465598, col = "red", lty = 2)

	# Enter in the axes
	make_typhoon_plot_axes(start_coord, end_coord)

	# Make the second typhoon plot
	screen(scr.s[3])
	par(mar = c(3, left_par, 3, right_par))
	dens_dot_plot(typhoon2.m, z_max = z_max, x_label = "", y_axt = "n", x_axt = "n", plot_title = typhoon2_title)
	abline(v = 465598, col = "red", lty = 2)

	# Enter in the axes
	make_typhoon_plot_axes(start_coord, end_coord)

	# Make the x-label
	screen(scr.s[4])
	par(mar = c(0, left_par, 0, right_par))
	set_chromatin_schematic()

	text(x = 0.5, y = 0.5, labels = paste("Chr ", chr, " (kb)", sep = ""), cex = 1.25)


	# Make the y-label
	screen(scr.s[5])
	par(mar = rep(0, 4))
	set_chromatin_schematic()

	text(x = 0.5, y = 0.5, labels = "Fragment Length (bp)", srt = 90, cex = 1.25)

	# Close the screen
	close.screen(all.screens = T)

	# Close the device
	dev.off()

}
