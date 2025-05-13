library(tuneR)
library(signal)
library(tools)

load_audio <- function(audiofile) {
  if (is.null(audiofile) || audiofile == "") {
    stop("The audiofile parameter cannot be NULL or empty.")
  }
  
  if (is.character(audiofile)) {
    if (!file.exists(audiofile)) {
      stop("The specified audio file does not exist.")
    }
    
    ext <- tolower(tools::file_ext(audiofile))
    if (ext == "mp3") {
      audio <- tuneR::readMP3(audiofile)
    } else if (ext == "wav") {
      audio <- tuneR::readWave(audiofile)
    } else {
      stop("Unsupported audio format. Supported formats are: .wav and .mp3.")
    }
    
    if (!inherits(audio, "Wave")) {
      stop("Failed to load the audio file as a valid 'Wave' object.")
    }
    return(audio)
    
  } else if (inherits(audiofile, "Wave")) {
    return(audiofile)
    
  } else {
    stop("Invalid input: provide a file path or a 'Wave' object.")
  }
}

get_audio_channel <- function(audio, channel = "stereo") {
  # Validate the channel parameter
  channel <- match.arg(channel, c("mono", "stereo", "left", "right"))
  
  # Mono channel
  if (channel == "mono") {
    if (!audio@stereo) {
      return(list(mono = audio@left)) # Already mono
    }
    return(list(mono = tuneR::mono(audio, which = "both")@left))
  }
  
  # Left channel
  if (channel == "left") {
    if (length(audio@left) == 0) {
      stop("The left channel is empty. Please provide a valid audio file with data in the left channel.")
    }
    return(list(left = audio@left))
  }
  
  # Right channel
  if (channel == "right") {
    if (!audio@stereo) {
      stop("The audio file is not stereo, so it does not have a right channel.")
    }
    if (length(audio@right) == 0) {
      stop("The right channel is empty. Please provide a valid stereo audio file.")
    }
    return(list(right = audio@right))
  }
  
  # Stereo channel
  if (channel == "stereo") {
    if (length(audio@left) == 0) {
      stop("The left channel is empty. Stereo audio requires a valid left channel.")
    }
    right <- if (audio@stereo && length(audio@right) > 0) audio@right else NULL
    return(list(left = audio@left, right = right))
  }
  
  # Fallback for unexpected cases
  stop("Unexpected error: invalid channel selection.")
}

bgnoise <- function(audiofile,
                    channel = "stereo",
                    time_bin = 60,
                    db_threshold = -90,
                    target_samp_rate = NULL,
                    wl = 512,
                    window = signal::hamming(wl),
                    overlap = ceiling(length(window) / 2),
                    histbreaks = "FD") {

  # Validate inputs
  if (!is.numeric(time_bin) || time_bin <= 0) stop("time_bin must be a positive number.")
  if (!is.numeric(db_threshold)) stop("db_threshold must be numeric.")
  if (!is.null(target_samp_rate) && (!is.numeric(target_samp_rate) || target_samp_rate <= 0)) {
    stop("target_samp_rate must be a positive number.")
  }
  if (!is.numeric(wl) || wl <= 0) stop("wl must be a positive number.")
  if (!is.numeric(overlap) || overlap < 0 || overlap > wl) stop("overlap must be between 0 and wl.")

  # Load audio
  audio <- load_audio(audiofile)

  # Downsample if needed
  if (!is.null(target_samp_rate)) {
    if (target_samp_rate >= audio@samp.rate) {
      stop("target_samp_rate must be less than the original sample rate.")
    }
    audio <- tuneR::downsample(audio, target_samp_rate)
  }

  samp.rate <- audio@samp.rate
  audio_dur <- length(audio) / samp.rate
  if (time_bin > audio_dur) stop("time_bin must be less than or equal to the audio duration.")
  frame_bin <- ceiling(audio_dur / time_bin)

  # Get audio channels
  channels <- get_audio_channel(audio, channel)

  # Load helper functions from internal_functions.R
  source("BGN_POW/internal_functions.R")

  # Process channels
  process <- function(data, ch) {
    process_channel(
      channel_data = data, ch = ch,
      time_bins = frame_bin, bin_size = time_bin,
      wl = wl, samp.rate = samp.rate, overlap = overlap,
      db_threshold = db_threshold, window = window, histbreaks = histbreaks
    )
  }

  if (channel == "mono") {
    results <- process(channels$mono, "mono")
    BGN_combined <- results$BGN
    POW_combined <- results$POW
  } else if (channel == "stereo") {
    left_results <- process(channels$left, "left")
    if (!is.null(channels$right)) {
      right_results <- process(channels$right, "right")
      BGN_combined <- cbind(left_results$BGN, right_results$BGN)
      POW_combined <- cbind(left_results$POW, right_results$POW)
    } else {
      BGN_combined <- left_results$BGN
      POW_combined <- left_results$POW
    }
  } else {
    results <- process(channels[[channel]], channel)
    BGN_combined <- results$BGN
    POW_combined <- results$POW
  }

  return(list(
    BGN = as.data.frame(BGN_combined),
    POW = as.data.frame(POW_combined)
  ))
}
