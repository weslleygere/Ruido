simple_sat <- function(soundfile,
                     time_bin = 60,
                     downsample = FALSE,
                     target_samp_rate = NULL,
                     wl = 512,
                     window = signal::hamming(wl),
                     overlap = ceiling(length(window) / 2),
                     channel = "both",
                     db_threshold = -90,
                     histbreaks = "FD",
                     powthr = 15,
                     bgnthr = 0.85) {
  source("BGN_POW/BGN_v0.1.R")
  
  BGN_POW <- bgnoise(
    soundfile,
    time_bin = time_bin,
    downsample = downsample,
    target_samp_rate = target_samp_rate,
    window = window,
    overlap = overlap,
    channel = channel,
    db_threshold = db_threshold,
    wl = wl,
    histbreaks = histbreaks
  )
  
  if (channel == "both") {
    BGN_left <- BGN_POW$BGN[grepl(toupper("left"), colnames(BGN_POW$BGN))]
    BGN_right <- BGN_POW$BGN[grepl(toupper("right"), colnames(BGN_POW$BGN))]
    
    BGN_Q_left <- setNames(quantile(unlist(BGN_left), probs = bgnthr), bgnthr)
    BGN_Q_right <- setNames(quantile(unlist(BGN_right), probs = bgnthr), bgnthr)
    
    BGN_saturation_left <- lapply(BGN_Q_left, function(Q)
      Q < BGN_left)
    BGN_saturation_right <- lapply(BGN_Q_right, function(Q)
      Q < BGN_right)
    
    BGN_saturation <- data.frame(Map(cbind, BGN_saturation_left, BGN_saturation_right))
    
  } else {
    BGN_Q <- setNames(quantile(unlist(BGN_POW$BGN), probs = bgnthr), bgnthr)
    
    BGN_saturation <- data.frame(lapply(BGN_Q, function(Q)
      Q < BGN_POW$BGN))
    
  }
  
  POW_saturation <- data.frame(lapply(powthr, function(Q)
    Q < BGN_POW$POW))
  
  sat_ex <- colSums(BGN_saturation |
                      POW_saturation) / (wl / 2)
  
  nbins <- length(sat_ex)
  
  names(sat_ex) <- if(channel == "both") {
    paste("SAT", rep(c("left", "right"), each = nbins / 2), seq(nbins / 2), sep = "_")
  } else {
    paste("SAT", channel, seq(nbins), sep = "_")
  }
  
  return(sat_ex)
  
}
