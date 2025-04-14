# Auxiliary function to grab each time bin from a spectrogram
get_spect_bins <- function(frame_start, frame_end, full_spec) {
  full_spec[, frame_start:frame_end]
}

# Auxiliary function to process each channel
process_channel <- function(channel_data,
                            time_bins =  frame_bin,
                            wl = wl,
                            samp.rate = audio@samp.rate,
                            overlap = overlap,
                            db_threshold = db_threshold,
                            window = window) {
  # Center the audio file around 0
  offset <- channel_data - mean(channel_data)
  
  # Calculate the spectrogram
  spect <- signal::specgram(
    x = offset,
    n = wl,
    Fs = samp.rate,
    window = window,
    overlap = overlap
  )
  
  # Calculate the frame size
  frame <- ceiling(length(spect$t) / time_bins)
  
  # Calculate the number of frequency bins
  freq_bins <- length(spect$f)
  
  # Calculate the absolute value of the coefficients (magnitude)
  spect_S <- abs(spect$S)
  
  # Convert to decibels
  spect_S_log <- 10 * log10(spect_S / max(spect_S))
  
  # Converting values lower than the decibel threshold to the threshold value
  spect_S_log[spect_S_log < db_threshold] <- db_threshold
  
  # Creating a data frame containing the beginning and the end of each time bin in the spectrogram time unit
  time_frames <- data.frame(t(sapply(seq(time_bins), function(time,
                                                              frame_length = frame,
                                                              total_frames = length(spect$t)) {
    frame_start <- 1 + frame_length * (time - 1)
    frame_end <- min(frame_start + frame_length - 1, total_frames)
    return(data.frame(frame_start, frame_end))
  })))
  
  # R shenanigans, converting the values from character to numeric
  time_frames <- data.frame(sapply(time_frames, as.numeric))
  
  # Generating the BGN and POW data frame
  # Understanding the craziness:
  # First step: We isolate each time bin of our spectrogram into a list
  # Second step: We apply Towsey's metodology to each time bin
  ## Towsey's BGN: For each frequency window we calculate a histogram of values (Towsey said to compute a 100 bin, but we decided to adopt the Freedman-Diaconis rule) of decibels values.
  ## We also get max decibels values to calculate soundscape power
  ## Then we export a data.frame containing the values of BGN and POW for each frequency of each time bin
  # This proccess give us a matrix with cols = wl/2 and two rows (BGN and POW) so we transponse them and store them in a data frame
  # Next we bind all these data frames into a single one
  
  BGN_POW_df <- data.frame(do.call(cbind, (lapply(lapply(apply(time_frames, 1, function(x) {
    get_spect_bins(x[1], x[2], spect_S_log)
  }), function(spec) {
    apply(spec, 1, function(x) {
      db_max <- max(x)
      
      histo <- hist(x = x,
                    plot = FALSE,
                    breaks = 'FD')  # Freedman-Diaconis rule)
      
      modal_intensity <- histo$mids[which.max(histo$counts)]
      
      c(BGN = modal_intensity, POW = db_max - modal_intensity)
    })
  }), function(bin) {
    data.frame(t(bin))
  }))))
  
  # Separating BGN and POW in different data frames
  BGN <- BGN_POW_df[, grepl("BGN", colnames(BGN_POW_df))]
  POW <- BGN_POW_df[, grepl("POW", colnames(BGN_POW_df))]
  
  colnames(BGN) <- paste0("BGN", seq(time_bins))
  colnames(POW) <- paste0("POW", seq(time_bins))
  
  # Return them both again but no as a list
  return(list(BGN = BGN, POW = POW))
}

