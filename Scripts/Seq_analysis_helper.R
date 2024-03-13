#function parser for cutadapt output

Cutadapt_parser <- 
  function(temp){
    
    i_1 <- which(grepl("Total read pairs processed", temp))
    i_2 <- which(grepl("Read 1 with adapter", temp))
    i_3 <- which(grepl("Read 2 with adapter", temp))
    i_4 <- which(grepl("Pairs that were too short", temp))
    i_5 <- which(grepl("Pairs written", temp))
    
    
      
    
      Total_reads = ifelse(!isEmpty(i_1), 
                           {str_split(temp[i_1], pattern = "\\s+")[[1]][5] %>% 
                               gsub(",", "", .) %>% 
                               as.numeric()},
                           NA)
      
      
      R1_with_adapter =   ifelse(!isEmpty(i_2), 
                                 {str_split(temp[i_2], pattern = "\\s+")[[1]][6] %>% 
                                     gsub(",", "", .) %>% 
                                     as.numeric()},
                                 NA)
      
      R1_with_adapter_prct = ifelse(!isEmpty(i_2), 
                                    {str_split(temp[i_2], pattern = "\\s+")[[1]][7] %>% 
                                        gsub("\\(|\\)|%", "", .) %>% 
                                        as.numeric()},
                                    NA)
      
      
      R2_with_adapter = ifelse(!isEmpty(i_3), 
                               {str_split(temp[i_3], pattern = "\\s+")[[1]][6] %>% 
                                   gsub(",", "", .) %>% 
                                   as.numeric()},
                               NA)
      
      R2_with_adapter_prct = ifelse(!isEmpty(i_3), 
                                    {str_split(temp[i_3], pattern = "\\s+")[[1]][7] %>% 
                                        gsub("\\(|\\)|%", "", .) %>% 
                                        as.numeric()},
                                    NA)
      
      
      Pairs_too_short = ifelse(!isEmpty(i_4), 
                               {str_split(temp[i_4], pattern = "\\s+")[[1]][6] %>% 
                                   gsub(",", "", .) %>% 
                                   as.numeric()},
                               NA)
      
      Pairs_too_short_prct = ifelse(!isEmpty(i_4),
                                    {str_split(temp[i_4], pattern = "\\s+")[[1]][7] %>% 
                                        gsub("\\(|\\)|%", "", .) %>% 
                                        as.numeric()},
                                    NA)
      
      
      Pairs_written = ifelse(!isEmpty(i_5),
                             {str_split(temp[i_5], pattern = "\\s+")[[1]][5] %>% 
                                 gsub(",", "", .) %>% 
                                 as.numeric()},
                             NA)
      
      Pairs_written_prct = ifelse(!isEmpty(i_5),
                                  {str_split(temp[i_5], pattern = "\\s+")[[1]][6] %>% 
                                      gsub("\\(|\\)|%", "", .) %>% 
                                      as.numeric()},
                                  NA)
      
      res <- list(Total_reads = Total_reads,
                  R1_with_adapter = R1_with_adapter,
                  R1_with_adapter_prct = R1_with_adapter_prct,
                  R2_with_adapter = R2_with_adapter,
                  R2_with_adapter_prct = R2_with_adapter_prct,
                  Pairs_too_short = Pairs_too_short,
                  Pairs_too_short_prct = Pairs_too_short_prct,
                  Pairs_written = Pairs_written,
                  Pairs_written_prct = Pairs_written_prct)
      
      
      res <- res[!is.na(res)]
      
      res <- as_tibble(res)
      
      return(res)
    
  }


# function to get all orientations of a primer
# from https://benjjneb.github.io/dada2/ITS_workflow.html

allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- Biostrings::DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = Biostrings::complement(dna), Reverse = Biostrings::reverse(dna), 
               RevComp = Biostrings::reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}

primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}


## create directory only if directory doesn't exist

dir.create_ifnot <- function(path){
  
  if(dir.exists(path)) message("directory exists")
  else dir.create(path) 
  
}


## function to find last common ancestor in taxonomy dataframe
coherent_taxonomy <- function(df) {
  # Check for coherence at each level for each sequence
  df_coherent <- df %>%
    summarise(across(everything(), ~ {
      # Ignore NA values in comparison
      unique_vals <- unique(na.omit(.))
      
      # If there is only one unique value, return it
      if (length(unique_vals) == 1) {
        return(unique_vals)
      }
      
      # Otherwise, set to NA
      return(NA)
    }))
  
  return(df_coherent)
}


#######
# alternative error model for DADA2 for binned Q-scores. 
# See https://github.com/benjjneb/dada2/issues/1307#issuecomment-957680971
#######

loessErrfun_mod4 <- function(trans) {
  qq <- as.numeric(colnames(trans))
  est <- matrix(0, nrow=0, ncol=length(qq))
  for(nti in c("A","C","G","T")) {
    for(ntj in c("A","C","G","T")) {
      if(nti != ntj) {
        errs <- trans[paste0(nti,"2",ntj),]
        tot <- colSums(trans[paste0(nti,"2",c("A","C","G","T")),])
        rlogp <- log10((errs+1)/tot)  # 1 psuedocount for each err, but if tot=0 will give NA
        rlogp[is.infinite(rlogp)] <- NA
        df <- data.frame(q=qq, errs=errs, tot=tot, rlogp=rlogp)
        
        # original
        # ###! mod.lo <- loess(rlogp ~ q, df, weights=errs) ###!
        # mod.lo <- loess(rlogp ~ q, df, weights=tot) ###!
        # #        mod.lo <- loess(rlogp ~ q, df)
        
        # jonalim's solution
        # https://github.com/benjjneb/dada2/issues/938
        mod.lo <- loess(rlogp ~ q, df, weights = log10(tot),degree = 1, span = 0.95)
        
        pred <- predict(mod.lo, qq)
        maxrli <- max(which(!is.na(pred)))
        minrli <- min(which(!is.na(pred)))
        pred[seq_along(pred)>maxrli] <- pred[[maxrli]]
        pred[seq_along(pred)<minrli] <- pred[[minrli]]
        est <- rbind(est, 10^pred)
      } # if(nti != ntj)
    } # for(ntj in c("A","C","G","T"))
  } # for(nti in c("A","C","G","T"))
  
  # HACKY
  MAX_ERROR_RATE <- 0.25
  MIN_ERROR_RATE <- 1e-7
  est[est>MAX_ERROR_RATE] <- MAX_ERROR_RATE
  est[est<MIN_ERROR_RATE] <- MIN_ERROR_RATE
  
  # enforce monotonicity
  # https://github.com/benjjneb/dada2/issues/791
  estorig <- est
  est <- est %>%
    data.frame() %>%
    mutate_all(funs(case_when(. < X40 ~ X40,
                              . >= X40 ~ .))) %>% as.matrix()
  rownames(est) <- rownames(estorig)
  colnames(est) <- colnames(estorig)
  
  # Expand the err matrix with the self-transition probs
  err <- rbind(1-colSums(est[1:3,]), est[1:3,],
               est[4,], 1-colSums(est[4:6,]), est[5:6,],
               est[7:8,], 1-colSums(est[7:9,]), est[9,],
               est[10:12,], 1-colSums(est[10:12,]))
  rownames(err) <- paste0(rep(c("A","C","G","T"), each=4), "2", c("A","C","G","T"))
  colnames(err) <- colnames(trans)
  # Return
  return(err)
}


