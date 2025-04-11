library(tuneR)
library(signal)
library(tools)

bgnoise <- function(audio_file_path,
                    time_bin = 60,
                    db_threshold = -90,
                    downsample = FALSE,
                    target_samp_rate = NULL,
                    wl = 512,
                    window = signal::hanning(wl),
                    overlap = ceiling(length(window) / 2)) {
  
  # Check the audio extension
  if (tools::file_ext(audio_file_path) %in% c("mp3", "wav")) {
    audio <- if (tools::file_ext(audio_file_path) == "mp3") {
      tuneR::readMP3(audio_file_path)
    } else {
      tuneR::readWave(audio_file_path)
    }
  } else {
    stop("The audio file must be in MP3 or WAV format.")
  }

  # Downsample if requested
  if (downsample) {
    if (is.null(target_samp_rate)) {
      stop("Please provide a target sample rate for downsampling.")
    }
    audio <- tuneR::downsample(audio, target_samp_rate)
  }
  
  # Extract the base name of the audio file (without extension)
  audio_file_name <- tools::file_path_sans_ext(basename(audio_file_path))
  
  # Calculate the duration of the audio file in seconds
  audio_dur <- length(audio) / audio@samp.rate
  
  # Calculate the number of frames per time bin
  frame_bin <- ceiling(audio_dur / time_bin)
  
  # Auxiliary function to process each channel
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
    
    # Calculate the absolute value of the coefficients (magnitude)
    spect_S <- abs(spect$S)
    
    # Convert to decibels
    spect_S_log <- 10 * log10(spect_S / max(spect_S))

    # Create matrices for background noise (BGN) and power (POW)
    BGN <- POW <- matrix(NA, nrow = freq_bins, ncol = frame_bin)
    
    colnames(BGN) <- paste0(channel_name, "_", seq_len(frame_bin))
    colnames(POW) <- paste0(channel_name, "_", seq_len(frame_bin))
    rownames(BGN) <- rownames(POW) <- paste0("f", seq_len(freq_bins))
    
    for (i in 1:frame_bin) {
      frame_start <- 1 + frame * (i - 1)
      frame_end <- min(frame_start + frame - 1, length(spect$t))
      
      ampdata <- spect_S_log[, frame_start:frame_end]
      ampdata[ampdata < db_threshold] <- db_threshold
      
      for (f in seq_len(nrow(ampdata))) {
        db_max <- max(ampdata[f, ])
        db_min <- min(ampdata[f, ])
        
        histo <- hist(
          x = ampdata[f, ],
          plot = FALSE,
          breaks = 'FD'  # Freedman-Diaconis rule
        )
        
        modal_intensity <- histo$mids[which.max(histo$counts)]
        
        BGN[f, i] <- modal_intensity
        POW[f, i] <- db_max - modal_intensity
      }
    }

    return (list(BGN = BGN, POW = POW))
  }

  # Process left channel
  left_channel <- audio@left
  left_results <- process_channel(left_channel, "left")

  # Process right channel if stereo
  if (audio@stereo) {
    right_channel <- audio@right
    right_results <- process_channel(right_channel, "right")
    
    # Combine matrices
    BGN_combined <- cbind(left_results$BGN, right_results$BGN)
    POW_combined <- cbind(left_results$POW, right_results$POW)
    
  } else {
    BGN_combined <- left_results$BGN
    POW_combined <- left_results$POW
    
  }

  # Convert to data frames
  BGN_df <- as.data.frame(BGN_combined)
  POW_df <- as.data.frame(POW_combined)

  # Export to CSV
  write.csv(BGN_df, file = paste0(audio_file_name, "_BGN_results.csv"))
  write.csv(POW_df, file = paste0(audio_file_name, "_POW_results.csv"))

  # Return results
  return(list(
    BGN = BGN_df,
    POW = POW_df,
  ))
}