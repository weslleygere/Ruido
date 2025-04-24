library(tuneR)
library(signal)
library(tools)

bgnoise <- function(audiofile,
                    channel = "both",
                    time_bin = 60,
                    db_threshold = -90,
                    downsample = FALSE,
                    target_samp_rate = NULL,
                    wl = 512,
                    window = signal::hamming(wl),
                    overlap = ceiling(length(window) / 2),
                    histbreaks = "FD") {
  
  source("BGN_POW/internal_functions.R")
  
  if(!channel %in% c("left", "right", "both", "mono")) {
    stop("Please provide a valid channel: 'left', 'right', 'both', or 'mono'.")
  }
  
  # Check the audio extension
  audio <- if (is.character(audiofile)) {
    if (tools::file_ext(audiofile) %in% c("mp3", "wav")) {
      if (tools::file_ext(audiofile) == "mp3") {
        tuneR::readMP3(audiofile)
      } else {
        tuneR::readWave(audiofile)
      }
    } else {
      stop("The audio file must be in MP3 or WAV format.")
    }
  } else {
    audiofile
  }
  
  # Downsample if requested
  if (downsample) {
    if (is.null(target_samp_rate)) {
      stop("Please provide a target sample rate for downsampling.")
    }
    audio <- tuneR::downsample(audio, target_samp_rate)
  }
  
  samp.rate <- audio@samp.rate
  
  # Extract the base name of the audio file (without extension)
  # audio_file_name <- tools::file_path_sans_ext(basename(audio_file_path))
  
  # Calculate the duration of the audio file in seconds
  audio_dur <- length(audio) / samp.rate
  
  # Calculate the number of frames per time bin
  frame_bin <- ceiling(audio_dur / time_bin)
  
  # Process right channel if stereo
  if (channel == "both") {
    # Process left channel
    left_results <- process_channel(audio@left, ch = "left", time_bins = frame_bin, bin_size = time_bin, wl = wl, samp.rate = samp.rate, overlap = overlap, db_threshold = db_threshold, window = window, histbreaks = histbreaks)
    
    # Process right channel if stereo
    if (audio@stereo) {
      right_results <- process_channel(channel_data = audio@right, ch = "right", time_bins = frame_bin, bin_size = time_bin, wl = wl, samp.rate = samp.rate, overlap = overlap, db_threshold = db_threshold, window = window, histbreaks = histbreaks)
      
      # Combine results
      BGN_combined <- cbind(left_results$BGN, right_results$BGN)
      POW_combined <- cbind(left_results$POW, right_results$POW)
    } else {
      BGN_combined <- left_results$BGN
      POW_combined <- left_results$POW
    }
    
  } else if (channel == "mono") {
    main_results <- process_channel(channel_data = tuneR::mono(audio)@left, ch = "mono", time_bins = frame_bin,  bin_size = time_bin, wl = wl, samp.rate = samp.rate, overlap = overlap, db_threshold = db_threshold, window = window, histbreaks = histbreaks)
    BGN_combined <- main_results$BGN
    POW_combined <- main_results$POW
    
  } else {
    desired_channel <- slot(audio, channel)
    
    if (length(desired_channel) == 0) {
      stop("Provided audio channel is empty!")
    }
    
    main_results <- process_channel(channel_data = desired_channel, ch = channel, time_bins = frame_bin, wl = wl,  bin_size = time_bin, samp.rate = samp.rate, overlap = overlap, db_threshold = db_threshold, window = window, histbreaks = histbreaks)
    BGN_combined <- main_results$BGN
    POW_combined <- main_results$POW
  }
  return(list(BGN = as.data.frame(BGN_combined), POW = as.data.frame(POW_combined)))
}
