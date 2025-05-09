sat_backup <- function(backup_path, 
                       od,
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
  
  source("BGN_POW/SATURATION.R")
  
  backup = backup_path
  
  backfiles <- list.files(backup_path, pattern = ".txt", full.names = TRUE)
  originalfiles <- list.files(od, full.names = TRUE)
  originalfiles <- originalfiles[tools::file_ext(originalfiles) %in% c("mp3", "wav")]
  
  remainingfiles <- originalfiles[!(basename(originalfiles) %in% basename(file_path_sans_ext(backfiles)))]
  
  powthreshold <- seq(powthr[1], powthr[2], powthr[3])
  names(powthreshold) <- powthreshold
  bgnthreshold <- seq(bgnthr[1], bgnthr[2], bgnthr[3])
  
  threshold_combinations <- setNames(expand.grid(powthreshold, bgnthreshold),
                                     c("powthreshold", "bgnthreshold"))
  
  combinations <- paste(threshold_combinations[, 1], threshold_combinations[, 2], sep = "/")
  
  half_wl <- wl / 2
  
  if (length(remainingfiles) == 0) {
    message("All files have already been processed.")
    
    SAT_df <- do.call(rbind,lapply(backfiles, read.table, header = TRUE))
    colnames(SAT_df) <- combinations
    
  } else {
    
    do.call(rbind, lapply(remainingfiles, function(soundfile) {
      BGN_POW <- tryCatch(
        bgnoise(
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
        ),
        error = function(e)
          e,
        warning = function(w)
          w
      )
      
      if (is(BGN_POW, "error") || is(BGN_POW, "warning")) {
        cat("\n", basename(soundfile), "is not valid!\nError:", BGN_POW$message, "\n")
      } else {
        BGN_Q <- setNames(quantile(unlist(BGN_POW$BGN), probs = seq(bgnthr[1], bgnthr[2], bgnthr[3])),
                          seq(bgnthr[1], bgnthr[2], bgnthr[3]))
        
        actv_size <- prod(dim(BGN_POW$BGN))
        
        BGN_saturation <- lapply(BGN_Q, function(Q)
          Q < BGN_POW$BGN)
        
        POW_saturation <- lapply(powthreshold, function(Q)
          Q < BGN_POW$POW)
        
        cat(
          "\r(",
          basename(soundfile),
          ") ",
          match(soundfile, originalfiles),
          " out of ",
          length(originalfiles),
          " remaining recordinds concluded!",
          sep = ""
        )
        
        sat_ex <- mapply(
          function(bgnthresh, powthresh) {
            colSums(BGN_saturation[[paste(bgnthresh)]] |
                      POW_saturation[[paste(powthresh)]]) / half_wl
          },
          threshold_combinations$bgnthreshold,
          threshold_combinations$powthreshold
        )
        
        rownames(sat_ex) <- paste(basename(soundfile), "SAT", 1:nrow(sat_ex), sep = "_")
        
        if (!is.null(backup)) {
          write.table(sat_ex, file = paste0(backup, "/", basename(soundfile), ".txt"), sep = "\t")
        }
        
        return(sat_ex)
      }
      
    }))
    
    SAT_df <- do.call(rbind,lapply(list.files(backup, pattern = ".txt",
                                              full.names = TRUE), read.table,
                                   header = TRUE))
    
    colnames(SAT_df) <- combinations
    
  }
  
  normal <- if (normality == "ks.test") {
    apply(SAT_df, 2, function(Q)
      ks.test(Q, pnorm)$p.value)
    
  } else if (normality == "shapiro.test") {
    apply(SAT_df, 2, function(x)
      ifelse(length(unique(x)) != 1, shapiro.test(x)$p.value, 0))
    
  } else {
    apply(SAT_df, 2, function(Q)
      eval(parse(text = paste0(
        normality, "(c(", paste((Q), collapse = ","), "))"
      )))$p.value)
    
  }
  
  thresholds <- unlist(strsplit(names(which.max(normal)), split = "/"))
  
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
  
  export <- list(
    powthresh = numeric(0),
    bgntresh = numeric(0),
    normality = numeric(0),
    values = data.frame()
  )
  
  export["powthresh"] <- as.numeric(thresholds[1])
  export["bgntresh"] <- as.numeric(thresholds[2]) * 100
  export[normality] <- as.numeric(as.numeric(max(normal)))
  export[["values"]] <- data.frame(SAT = SAT_df[, which.max(normal)],
                                   row.names = rownames(SAT_df))
  
  return(export)
  
}
