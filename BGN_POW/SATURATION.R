soundsat <- function(soundpath,
                     time_bin = 60,
                     downsample = FALSE,
                     target_samp_rate = NULL,
                     wl = 512,
                     window = signal::hamming(wl),
                     overlap = ceiling(length(window) / 2),
                     channel = "left",
                     db_threshold = -90,
                     histbreaks = "FD",
                     powthr = c(5.1, 20, 0.1),
                     bgnthr = c(0.51, 0.99, 0.02),
                     normality = "shapiro.test",
                     n_files = NULL) {
  # Loading necessary packages
  require(tuneR)
  require(signal)
  
  source("BGN_POW/BGN_v0.1.R")
  
  # Checking if the path actually exists and stopping the fuction in case it doesnt
  if (all(!dir.exists(soundpath)))
    stop("all provided soundpaths must be valid.")
  
  # Creating an object with the paths to all the recordings in the given folder
  soundfiles <- list.files(soundpath, full.names = TRUE, recursive = TRUE)
  soundfiles <- soundfiles[tools::file_ext(soundfiles) %in% c("mp3", "wav")]
  
  # For debugging purposes, only reads the first n files
  # Might become something official if I find a way to actually implement this in a better way
  if (!is.null(n_files)) {
    soundfiles <- soundfiles[1:n_files]
  }
  
  # Stop the function if the path does not contain at least 3 audio files
  # Using this function to calculate the saturation of only 1 or 2 audios does not make sense
  # I'm currently working on a fuction to calculate the saturation of only one audio
  if (length(soundfiles) < 3)
    stop("you must provide at least 3 recordings!")
  
  # Creating individual vectors for POW and BGN thresholds
  # POW thresholds are theorical values and BGN thresholds are quantile
  powthreshold <- seq(powthr[1], powthr[2], powthr[3])
  names(powthreshold) <- powthreshold
  bgnthreshold <- seq(bgnthr[1], bgnthr[2], bgnthr[3])
  
  # Creating a data.frame with each possible combination of the thresholds
  # Each row will contain a unique combination of the thresholds
  threshold_combinations <- setNames(expand.grid(powthreshold, bgnthreshold),
                                     c("powthreshold", "bgnthreshold"))
  
  # Creating a vector with the "name" of each combination
  # By "name" I mean a unique character string to each combination containing both of the threshold values
  combinations <- paste(threshold_combinations[, 1], threshold_combinations[, 2], sep = "/")

  # Making sure the user understand what he is about to get into.
  # Might add a timer in the future.
  print(
    paste(
      "Calculating saturation values for",
      length(soundfiles),
      "recordings using",
      length(combinations),
      "threshold combinations"
    )
  )
  
  # Creating a object to hold the normality values for the future
  # R is a weird language, so I did this to avoid errors in the future
  normal <- c()

# Function's main part.
  # First we need to calculate the BGN and POW values for each recoding on the given path, then we are going to create an "activity matrix" to determine the saturation.
  
  SAT_df <- data.frame(t(sapply(soundfiles, function(soundfile) {

    # Getting the recording BGN and POW values
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

    # Getting the quantile of background noise of the current recording to calculate the activity with BGN values
    BGN_Q <- setNames(quantile(unlist(BGN_POW$BGN), probs = seq(bgnthr[1], bgnthr[2], bgnthr[3])),
                      seq(bgnthr[1], bgnthr[2], bgnthr[3]))

    # Checking which values are above the threshold for BGN and POW
    BGN_saturation <- sapply(BGN_Q, function(Q)
      Q < BGN_POW$BGN)
    POW_saturation <- sapply(powthreshold, function(Q)
      Q < BGN_POW$POW)

    # Little message to tell the user the current progress of the function
    cat(
      "\r(",
      basename(soundfile),
      ") ",
      match(soundfile, soundfiles),
      " out of ",
      length(soundfiles),
      " recordinds concluded!",
      sep = ""
    )

        # This is the part where we calculate the saturation values for each threshold combination
        # The function will return a value between 0 and 1, where 0 means that no frame of the audio file was above the threshold and 1 means that all frames were above the threshold
    mapply(
      function(bgnthresh, powthresh) {
        sum(BGN_saturation[, paste(bgnthresh)] |
              POW_saturation[, paste(powthresh)]) /
          nrow(BGN_saturation)
        
      },
      threshold_combinations$bgnthreshold,
      threshold_combinations$powthreshold
    )
    
  })))
  
  colnames(SAT_df) <- combinations
  
  # Now, after we get the saturation values for every single threshold, we are going to the weird part of the function
  # In this index, we must calculate the normality of the data to define which threshold holds the correct information of saturation
  # We run a homosexuality test (I set shapiro.test as default since it generally yields "nicer looking" results, although using ks.test may yield more realistic results) to the values of saturation of each threshold and then we grab which threshold holds the more normal results
  # This is the quickiest part of the whole function.
    
  normal <- if (normality == "ks.test") {
    sapply(SAT_df, function(Q)
      ks.test(Q, pnorm)$p.value)
    
  } else if (normality == "shapiro.test") {
    as.numeric(sapply(SAT_df, function(uni)
      ifelse(length(unique(uni)) != 1, shapiro.test(uni)$p.value, 0)))
    
  } else {
    sapply(SAT_df, function(Q)
      eval(parse(text = paste0(
        normality, "(", Q, ")"
      )))$p.value)
    
  }
  
  # After running the normality tests, we name the strings to get which combination yields our results
  
  names(normal) <- combinations
  
  # Now we grab the threshold combination that yield the most normal results.
  # For shapiro.test, the higher the number, the higher the normality
  
  thresholds <- unlist(strsplit(names(which.max(normal)), split = "/"))
  
  # Here's a nice message telling the user the values of both thresholds and the result of the most normal test
  
  cat(
    "\n           Soundscape Saturation Results\n\n",
    "POW Threshold = ",
    as.numeric(thresholds[1]),
    " kHz        ",
    "BGN Threshold = ",
    as.numeric(thresholds[2]) * 100,
    "%\n",
    "            Normality test result = ",
    as.numeric(max(normal)),
    "\n ",
    sep = ""
  )
  
  # Here is the file we are going to export
  # In this list contains values of both thresholds, result of normality test and the values of saturation of the threshold that yield the most normal results
  
  export <- list(
    powthresh = numeric(0),
    bgntresh = numeric(0),
    normality = numeric(0),
    values = data.frame()
  )
  
  export["powthresh"] <- as.numeric(thresholds[1])
  export["bgntresh"] <- as.numeric(thresholds[2]) * 100
  export[normality] <- as.numeric(as.numeric(max(normal)))
  export[["values"]] <- data.frame(soundfile = basename(rownames(SAT_df)), SAT = SAT_df[[which.max(normal)]])
  
  return(export)
  
}
