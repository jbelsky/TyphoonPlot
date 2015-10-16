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
                    stats::density(read_pos.v, bw = bw, weights = rep(1, length(read_pos.v)),
                                   from = start, to = end, n = end - start + 1
                                  )
                  )

  # Convert to a numeric vector
  read_pos.v = read_pos.dens$y
  names(read_pos.v) = read_pos.dens$x

  # Scale to the number of reads
  return(read_pos.v)

}






GetDensityPeaks = function(cov.v, x_mid, min_thresh = 1, peak_width = 75){

	# Set the peak window
	peak_win = (peak_width - 1) / 2

	# Perform an FFT on the coverage
	cov_fft.v = nucleR::filterFFT(cov.v, pcKeepComp = 0.02)

	# Find the peaks
	cov_peaks.v = nucleR::peakDetection(cov_fft.v, threshold = min_thresh, score = F, width = 1)

	if(!is.null(cov_peaks.v)){

		# Convert to GenomicRanges
		peaks.gr = GenomicRanges::GRanges(seqnames = 1,
				                              ranges = IRanges::IRanges(start = cov_peaks.v - 1, width = 1)
				                             )

		# Enter in the peak position and the signal
		GenomicRanges::values(peaks.gr)$peak = IRanges::mid(IRanges::ranges(peaks.gr))
		GenomicRanges::values(peaks.gr)$sig = cov_fft.v[GenomicRanges::values(peaks.gr)$peak]

		neg_idx = numeric()

		# Ensure that the peak is the highest signal in the range
		for(i in 1:length(peaks.gr)){

			# Get the maximum signal in the peak_width range
			start = GenomicRanges::values(peaks.gr)$peak[i] - peak_win
			end = GenomicRanges::values(peaks.gr)$peak[i] + peak_win

			if(start < 1){
				start = 1
			}else if(end > length(cov_fft.v)){
				end = length(cov_fft.v)
			}

			peak_range_sig = max(cov_fft.v[start:end])

			if(peak_range_sig > GenomicRanges::values(peaks.gr)$sig[i]){

				neg_idx = c(neg_idx, i)

			}

		}

		if(length(neg_idx) > 0){

			# Get the peak listing
			peaks.gr = peaks.gr[-neg_idx]

		}

		# Get the win
		win = (length(cov_fft.v) - 1)/2

		# Get the peaks
		peaks.df = data.frame(pos = GenomicRanges::values(peaks.gr)$peak - win,
				                  sig = GenomicRanges::values(peaks.gr)$sig
				                 )

		# Remove peaks within 40 bp of either end
		peaks.df = peaks.df[which(peaks.df$pos >= (-win + 40) & peaks.df$pos <= (win - 40)),]

	}else{

		peaks.df = data.frame(pos = numeric(), sig = numeric())

	}

	# Adjust the position
	peaks.df$pos = x_mid + peaks.df$pos

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
