# Auxiliary function to grab each minute of a recording
get_sample_bins <- function(samples, samp.rate, bin_size) {
  b <- seq(1, samples, by = samp.rate * bin_size)
  e <- pmin(b + samp.rate * bin_size - 1, samples)
  
  if (length(b) == length(e)) {
    data.frame(b, e)
  } else {
    data.frame(b, e = c(e, samples))
  }
  
}

process_channel <- function(channel_data,
                            channel,
                            time_bin,
                            wl,
                            overlap,
                            db_threshold,
                            window,
                            histbreaks) {
  
  samp.rate <- channel_data@samp.rate
  
  frame_bin <- ceiling(length(channel_data) / samp.rate) / time_bin
  
  all_samples <- get_sample_bins(length(channel_data), samp.rate, time_bin)
  
  channel_data <- switch(
    channel,
    "stereo" = list("left" = channel_data@left, "right" = channel_data@right),
    "mono" = list(mono = tuneR::mono(channel_data, "both")@left),
    setNames(list(slot(channel_data, channel)), channel)
  )
  
  BGN_exp <- lapply(channel_data, function(x) {
    
    offset <- x - mean(x)
    
    temp_holder <- apply(all_samples, 1, function(y) {
      list(
        signal::specgram(
          x = offset[y[1]:y[2]],
          n = wl,
          Fs = samp.rate,
          window = window,
          overlap = overlap
        )$S
      )
    })
    
    BGN_POW_df <- data.frame(do.call(cbind, lapply(lapply(temp_holder, function(single_bin) {
      spect_S <- abs(single_bin[[1]])
      
      spect_S <- 10 * log10(spect_S / max(spect_S))
      
      spect_S[spect_S < db_threshold] <- db_threshold
      
      apply(spect_S, 1, function(x) {
        db_max <- max(x)
        db_min <- min(x)
        
        bin_width <- 2 * IQR(x) / length(x)^(1 / 3)
        
        num_bins <- ifelse(is.numeric(histbreaks), histbreaks, eval(parse(
          text = paste0("nclass.", histbreaks, "(x)")
        )))
        
        modal_intensity <- db_min + ((which.max(tabulate(
          findInterval(
            x = x,
            vec = seq(db_min, db_max, length.out = num_bins)
          )
        ))) * bin_width)
        
        c(BGN = modal_intensity, POW = db_max - modal_intensity)
      })
      
    }), function(x)
      data.frame(t(x)))))
    
    colnames(BGN_POW_df) <- paste0(rep(c("BGN", "POW"), frame_bin), rep(1:frame_bin, each = 2))
    
    # Separating BGN and POW in different data frames
    BGN <- data.frame(BGN_POW_df[, grepl("BGN", colnames(BGN_POW_df))])
    POW <- data.frame(BGN_POW_df[, grepl("POW", colnames(BGN_POW_df))])
    
    # Return them both again but no as a list
    return(list(BGN = BGN, POW = POW))
    
  })
  
  BGN_exp[["time_bins"]] <- setNames(round((all_samples$e - all_samples$b) / samp.rate),
                                 paste0("BIN", seq(frame_bin)))
  
  return(BGN_exp)
  
}
