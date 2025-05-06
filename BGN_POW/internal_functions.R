# Auxiliary function to grab each minute of a recording
get_sample_bins <- function(samples, samp.rate, time_bins, bin_size) {
  b <- seq(1, samples, by = samp.rate * bin_size)
  e <- pmin(b + samp.rate * bin_size - 1, samples)
  
  if (length(b) == length(e)) {
    data.frame(b, e)
  } else {
    data.frame(b, e = c(e, samples))
  }
  
}

process_channel <- function(channel_data,
                            time_bins,
                            ch,
                            bin_size,
                            wl,
                            samp.rate,
                            overlap,
                            db_threshold,
                            window,
                            histbreaks) {
  # Center the audio file around 0
  offset <- channel_data - mean(channel_data)
  
  all_samples <- get_sample_bins(length(offset), samp.rate, time_bins, bin_size)
  
  temp_holder <- apply(all_samples, 1, function(x) {
    list(signal::specgram(
      x = offset[x[1]:x[2]],
      n = wl,
      Fs = samp.rate,
      window = window,
      overlap = overlap
    )$S)
  })
  
  BGN_POW_df <- data.frame(do.call(cbind, lapply(lapply(temp_holder, function(single_bin) {
    
    spect_S <- abs(single_bin[[1]])
    
    # Convert to decibels
    spect_S <- 10 * log10(spect_S / max(spect_S))
    
    # Converting values lower than the decibel threshold to the threshold value
    spect_S[spect_S < db_threshold] <- db_threshold
    
    apply(spect_S, 1, function(x) {
      db_max <- max(x)
      db_min <- min(x)
      
      bin_width <- 2 * IQR(x) / length(x)^(1 / 3)
      
      num_bins <- ifelse(is.numeric(histbreaks),
                         histbreaks,
                         eval(parse(text = paste0("nclass.", histbreaks, "(x)"))))
      
      modal_intensity <- db_min + ((which.max(tabulate(findInterval(
        x = x,
        vec = seq(db_min, db_max, length.out = num_bins)
      )))) * bin_width)
      
      c(BGN = modal_intensity, POW = db_max - modal_intensity)
    })
    
  }), function(x)
    data.frame(t(x)))))
  
  # Separating BGN and POW in different data frames
  BGN <- setNames(data.frame(BGN_POW_df[, grepl("BGN", colnames(BGN_POW_df))]), paste0(paste0(toupper(ch), "_BGN", seq(time_bins))))
  POW <- setNames(data.frame(BGN_POW_df[, grepl("POW", colnames(BGN_POW_df))]), paste0(paste0(toupper(ch), "_POW", seq(time_bins))))
  
  # Return them both again but no as a list
  return(list(BGN = BGN, POW = POW))
  
}
