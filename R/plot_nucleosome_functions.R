# Get the 1-D feature density
GetMNaseFeatureDensity = function(bam_file_name, chr, start, end, fragL, fragH, bw){

  # Get the MNase Reads
  mnase_reads.gr = ConvertPairedReadBAMToGR(bam_file_name, chr, start_pos = start - 250, end_pos = end)

  # Subset on the reads in the range
  reads_subset.gr = mnase_reads.gr[BiocGenerics::width(mnase_reads.gr) >= fragL &
                                   BiocGenerics::width(mnase_reads.gr) <= fragH
                                  ]

  # Get the vector of reads
  read_pos.v = IRanges::mid(IRanges::ranges(reads_subset.gr))

  # Create the density (Note: not a true density since it doesn't sum to 1)
  read_pos.dens = suppressWarnings(
                    stats::density(read_pos.v, bw = bw, # weights = rep(1, length(read_pos.v)),
                                   from = start, to = end, n = end - start + 1
                                  )
                  )

  # Convert to a numeric vector
  read_pos.v = read_pos.dens$y * length(read_pos.v) / (end - start + 1)
  names(read_pos.v) = read_pos.dens$x

  # Scale to the number of reads
  return(read_pos.v)

}






GetDensityPeaks = function(cov.v, peak_width = 75, isPeakMax = TRUE){

	# Set the peak window
	peak_win = (peak_width - 1) / 2

	# Find the peaks
	cov_peaks.v = as.numeric(names(cov.v)[splus2R::peaks(x = cov.v, span = peak_width, strict = isPeakMax)])

	# Adjust the position
	peaks.df = data.frame(pos = cov_peaks.v,
	                      sig = cov.v[as.character(cov_peaks.v)]
	                     )

	return(peaks.df)

}

# Make the nucleosome
PlotNucleosome = function(nuc.df, y_max = 2, y0 = 0.5, yh = 0.2, nuc_col = "#FF0000"){

	# Set up the angle vector
	theta = seq(0, 2 * pi, length = 1000)

	# Set the y
	y = y0 + yh * sin(theta)

	for(i in 1:nrow(nuc.df)){

		# Get the position
		x0 = nuc.df$pos[i]

		# Find the coordinates for the nucleosome at each theta position
		x = x0 + 75 * cos(theta)

		# Find the signal color shading
		sig_shade = round(100 * nuc.df$sig[i] / y_max)

		if(sig_shade > 99){
			sig_shade = 99
		}

		sig_shade_str = formatC(sig_shade, flag = "0#", format = "d", width = 2)

		# Plot the nucleosome
		polygon(x, y, col = paste(nuc_col, sig_shade_str, sep = ""))

	}

}

PlotNucleosomeModel = function(bam_file_name, chr, start_coord, end_coord,
                               low_frag = 150, high_frag = 175,
                               min_peak_thresh = 0.2
                               ){

  # Get the nucleosome density in the region
  nuc_dens.v = GetMNaseFeatureDensity(bam_file_name, chr, start_coord, end_coord, low_frag, high_frag, bw = 20)

  # Find the signal peaks
  nuc_peaks.df = GetDensityPeaks(nuc_dens.v,
                                 x_mid = (start_coord + end_coord) / 2,
                                 min_thresh = min_peak_thresh
                                )

  # Set up the chromatin schematic
  SetChromatinSchematic(x_start = start_coord, x_end = end_coord, y_start = 0, y_end = 1)

  # Plot the nucleosomes
  PlotNucleosome(nuc_peaks.df)

}
