# Function to convert a BAM file to a GR file
# By default reads in all the reads for the entire chromosome
# 	Otherwise, specify a particular start_pos and end_pos for just a specific region
ConvertPairedReadBAMToGR = function(bam_file_name, chr, start_pos = 1, end_pos = -1){

	# Create the BAM File object
	bf = Rsamtools::BamFile(bam_file_name, index = paste(bam_file_name, ".bai", sep = ""))

	# Get the chr list
	chr_length.v = Rsamtools::scanBamHeader(bf)$targets

	# Update the end_pos if necessary
	if(end_pos == -1){

		# Update the end_pos
		end_pos = chr_length.v[chr]

	}

	# Make a GR file for the chromosome
	chr.gr = GenomicRanges::GRanges(seqnames = chr,
			                            ranges = IRanges::IRanges(start = max(start_pos - 250, 1), end = end_pos),
			                            strand = "*"
			                           )

	# Specify the scan bam paramaeters
	p = Rsamtools::ScanBamParam(what = c("pos", "isize"),
                              which = chr.gr,
                              flag = Rsamtools::scanBamFlag(isMinusStrand = FALSE)
                             )

	# Get the reads that meet these conditions
	reads.l = Rsamtools::scanBam(bf, param = p)

	if(length(reads.l[[1]][["pos"]]) > 0){

		# Convert these reads to a GR object
		IP.gr = GenomicRanges::GRanges(seqnames = factor(chr, levels = names(chr_length.v)),
				                           ranges = IRanges::IRanges(start = reads.l[[1]][["pos"]],
                                                             width = reads.l[[1]][["isize"]]
                                                            ),
				                           strand = "*"
			                            )
		GenomeInfoDb::seqlengths(IP.gr) = chr_length.v

	}else{

		IP.gr = GenomicRanges::GRanges()

	}

	return(IP.gr)

}
