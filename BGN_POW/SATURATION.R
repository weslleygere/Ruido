soundsat <- function(soundpath,
                     time_bin = 60,
                     downsample = FALSE,
                     target_samp_rate = NULL,
                     wl = 512,
                     window = signal::hamming(wl),
                     overlap = ceiling(length(window) / 2),
                     channel = "left",
                     powthr = c(5.1, 20, 0.1),
                     bgnthr = c(0.51, 0.99, 0.02),
                     normality = "shapiro.test",
                     n_files = NULL) {
  
  # Loading necessary packages
  require(tuneR)
  require(signal)

  # Checking if the path actually exists and stopping the fuction in case it doesnt
  if (!dir.exists(soundpath))
    stop("soundpath must contain the directory to your recordings!")
  
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
  bgnthreshold <- seq(bgnthr[1], bgnthr[2], bgnthr[3])
  
  # Creating a data.frame with each possible combination of the thresholds
  # Each row will contain a unique combination of the thresholds
  threshold_combinations <- expand.grid(powthreshold, bgnthreshold)
  colnames(threshold_combinations) <- c("powthreshold", "bgnthreshold")
  
  # Creating a vector with the "name" of each combination
  # By "name" I mean a unique character string to each combination containing both of the threshold values
  combinations <- paste(threshold_combinations[, 1], threshold_combinations[, 2], sep = "/")
  
  # Creating a object to hold the normality values for the future
  # R is a weird language, so I did this to avoid errors in the future
  normal <- c()
  
  # Creating the data.frame to store the values of saturation of each file and to each threshhold
  SAT_df <- data.frame(matrix(
    NA,
    nrow = length(soundfiles),
    ncol = length(combinations)
  ),
  row.names = basename(soundfiles))
  colnames(SAT_df) <- combinations
  
  # Begin the main loop
  # The loop will repeat to each provided file and will calculate the values of saturation for every threshold combinations
  for (soundfile in soundfiles) {
    
    # Storing the current file path
    # This will be useful in the future
    soundname <- soundfile
    
    # Reading the current file into the R environment
    c_soundfile <- if (tools::file_ext(soundfile) == "mp3") {
      readMP3(soundfile)
    } else {
      readWave(soundfile)
    }
    
    # Using the bgnoise function to extract values of background noise and soundscape power.
    # Currently I'm defaulting to use just one channel since this functions takes way too long to actually process everything. You would need a supercomputer if you wanted to extract saturation values from both left and right channels, but this can change in the future.
    
    BGN_POW <- bgnoise(
      c_soundfile,
      time_bin = time_bin,
      downsample = downsample,
      target_samp_rate = target_samp_rate,
      window = window,
      overlap = overlap,
      channel = channel
    )
    
    # Getting the quantile values of BGN necessary to calculate the index
    # I rename the values to access them more easily later on
    BGN_Q <- quantile(unlist(BGN_POW$BGN), probs = seq(bgnthr[1], bgnthr[2], bgnthr[3]))
    names(BGN_Q) <- seq(bgnthr[1], bgnthr[2], bgnthr[3])
    
    # Now we calculate the saturation values to each threshold
    # In this step we create two data.frames containing only TRUE or FALSE values
    # To simplify this, I'll give and example, if the current threshold for POW is 5.1 dB, this will check each value of POW given by the "bgnoise" function to each frame
    # If a value is greater than the current threshold, we will assign TRUE to its position in the data.frame "POW_saturation", if the value is not greater than the current threshold we assign the frame with FALSE
    # TRUE = 1 and FALSE = 0
    # I'm using sapply to apply the same function to every single element of the BGN_POW data.frame. Why sapply? Because it returns a vector and its simplier to use
    
    BGN_saturation <- sapply(BGN_Q, function(Q)
      Q < BGN_POW$BGN)
    POW_saturation <- sapply(powthreshold, function(Q)
      Q < BGN_POW$POW)
    colnames(POW_saturation) <- powthreshold
    
    # Loop inside the main loop
    # In this loop we are going to check the saturation of the current loaded file for every possible combination of our thresholds
    
    for (combination in 1:length(combinations)) {
      
      SAT_df[basename(soundname), combinations[combination]] <- 
        sum(BGN_saturation[, paste(threshold_combinations$bgnthreshold[combination])] |
        POW_saturation[, paste(threshold_combinations$powthreshold[combination])])/
        nrow(BGN_saturation)
      
    }
    
    # Since it takes so long to calculate, this is a nice message to tell the user that something is actually happening
    
    cat(
      "\r(",
      basename(soundname),
      ") ",
      match(soundname, soundfiles),
      " out of ",
      length(soundfiles),
      " recordinds concluded!",
      sep = ""
    )
    
  }
  
  # Now, after we get the saturation values for every single threshold, we are going to the weird part of the function
  # In this index, we must calculate the normality of the data to define which threshold holds the correct information of saturation
  # We run a homosexuality test (I set shapiro.test as default since it generally yields "nicer" results, although using ks.test may yield more realistic results) to the values of saturation of each threshold and then we grab which threshold holds the more normal results
  
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
    "      Normality test result = ",
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
  export[["values"]] <- data.frame(soundfile = rownames(SAT_df), SAT = SAT_df[[which.max(normal)]])
  
  return(export)
  
}
