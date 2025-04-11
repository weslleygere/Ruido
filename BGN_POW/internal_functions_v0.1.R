process_channel <- function(channel_data, channel_name) {
  # Center the audio file around 0
  offset <- channel_data - mean(channel_data)
  
  # Calculate the spectrogram
  spect <- signal::specgram(
    x = offset,
    n = wl,
    Fs = audio@samp.rate,
    window = window,
    overlap = overlap
  )
  
  # Calculate the frame size
  frame <- ceiling(length(spect$t) / frame_bin)
  
  # Calculate the number of frequency bins
  freq_bins <- length(spect$f)
  
  # Calculate the absolute value of the coefficients in order to select only the real part (magnitude)
  spect_S <- abs(spect$S)
  
  # Calculate the dB value of the magnitude
  spect_S_log <- 10 * log10(spect_S / max(spect_S))
  
  # Create two matrices to store the BGN and POW values
  BGN <- POW <- matrix(NA, nrow = freq_bins, ncol = frame_bin)
  
  colnames(BGN) <- paste0(channel_name, "_", seq_along(1:frame_bin))
  colnames(POW) <- paste0(channel_name, "_", seq_along(1:frame_bin))
  rownames(BGN) <- rownames(POW) <- paste0("f", seq_along(1:freq_bins))
  
  for (i in 1:frame_bin) {
    # Define the start and end of the frame of each iteration, which corresponds to time_bin seconds of audio
    frame_start <- 1 + frame * (i - 1)
    frame_end <- ifelse(frame_start + frame - 1 > length(spect$t),
                        length(spect$t),
                        frame_start + frame - 1)
    
    # Select the amplitude data for the current frame
    ampdata <- spect_S_log[, frame_start:frame_end]
    # Set values below the minimum dB to the minimum dB
    ampdata[ampdata < db_threshold] <- db_threshold
    
    for (f in seq_len(nrow(ampdata))) {
      db_max <- max(ampdata[f, ])
      db_min <- min(ampdata[f, ])
      
      ###############################
      ######### Remover aqui ########
      # Applying Freedman-Diaconis rule to calculate the number of bins
      IQR <- quantile(ampdata[f, ], 0.75) - quantile(ampdata[f, ], 0.25)
      bin_width <- 2 * IQR / length(ampdata[f, ])^(1 / 3)
      num_bins <- ceiling((db_max - db_min) / bin_width)
      
      # Create the histogram with the calculated number of bins
      histo <- hist(
        x = ampdata[f, ],
        breaks = seq(db_min, db_max, length.out = num_bins + 1),
        plot = FALSE
      )
      
      # Find the mode of the histogram
      modal_intensity <- histo$mids[which.max(histo$counts)]
      
      # Calculate the BGN and POW values
      BGN[f, i] <- modal_intensity
      POW[f, i] <- db_max - modal_intensity
    }
  }
  
  return (list(BGN = BGN, POW = POW))
}
