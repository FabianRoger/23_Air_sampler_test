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


### bbduk output parser

parse_bbduk_output <- function(output) {
  
  require(stringr)
  
  # Split the output into lines
  lines <- unlist(str_split(output, "\n"))
  
  # Extract the lines with relevant information
  input_line <- lines[grepl("Input:", lines)]
  contaminants_line <- lines[grepl("Contaminants:", lines)]
  low_quality_line <- lines[grepl("Low quality discards:", lines)]
  total_removed_line <- lines[grepl("Total Removed:", lines)]
  result_line <- lines[grepl("Result:", lines)]
  
  # Extract values from the lines
  input_reads <- as.numeric(str_extract(input_line, "\\d+"))
  contaminants_reads <- as.numeric(str_extract(contaminants_line, "\\d+"))
  contaminants_prct <- str_extract(contaminants_line, "\\d+\\.?\\d*%") %>% str_replace("%", "") %>% as.numeric()
  low_quality_discard_reads <- as.numeric(str_extract(low_quality_line, "\\d+"))
  low_quality_discard_prct <- str_extract(low_quality_line, "\\d+\\.?\\d*%") %>% str_replace("%", "") %>% as.numeric()
  total_removed_reads <- as.numeric(str_extract(total_removed_line, "\\d+"))
  total_removed_prct <- str_extract(total_removed_line, "\\d+\\.?\\d*%") %>% str_replace("%", "") %>% as.numeric()
  result_reads <- as.numeric(str_extract(result_line, "\\d+"))
  result_prct <- str_extract(result_line, "\\d+\\.?\\d*%") %>% str_replace("%", "") %>% as.numeric()
  
  # Create the tibble
  tibble(
    Input_reads = input_reads,
    PhiX_reads = contaminants_reads,
    PhiX_prct = contaminants_prct,
    N_reads = low_quality_discard_reads,
    N_prct = low_quality_discard_prct,
    Total_removed_reads = total_removed_reads,
    Total_removed_prct = total_removed_prct,
    Result_reads = result_reads,
    Result_prct = result_prct
  )
}


#parse ITSx output - function form Markus Schlegel @WSL

#' Reads a positions file from ITSx
read_itsx = function(prefix) {
  cols = c('OTU', 'len', 'SSU', 'ITS1', 'S58', 'ITS2', 'LSU', 'comment')
  domains = c('SSU', 'ITS1', 'S58', 'ITS2', 'LSU')
  x = read.delim(
    paste0(prefix, '.positions.txt'),
    header=F,
    col.names=cols,
    stringsAsFactors=F
  )
  x = cbind(
    x['OTU'], 
    apply(x[domains], 2, function(x) {
      out = sapply(strsplit(x, ': ', fixed=T), '[', 2)
      ifelse(out == 'Not found', NA, out)
    }),
    x['comment']
  )
  # remove --END--
  x = x[x$OTU != '--END--', ]
  x$ITSx_cat = with(x, ifelse(
    is.na(comment),
    'ok',
    gsub(' ', '_', gsub('! *$', '', gsub('only partial ', 'partial_', gsub('Broken or partial sequence, ', '', comment))))
  ))
  out_cols = c('OTU', 'SSU', 'ITS1', 'ITS2', 'LSU', 'ITSx_cat')
  nondetects = read.delim(
    paste0(prefix, '_no_detections.txt'),
    header=F,
    col.names=out_cols,
    stringsAsFactors=F
  )
  nondetects$ITSx_cat = 'no_ITS'
  rbind(
    x[, out_cols],
    nondetects
  )
}


#dataframe to match taxa codes
tax_match_df <- 
  tibble(taxonomy_level = c("d", "k", "p", "c", "o", "f", "g", "s", "SH"),
         taxon = c("devision","kingdom", "phylum","class","order","family","genus","species", "SH"))




# optimized lulu function
lulu_optimized <- function(
    otutable,
    matchlist,
    minimum_ratio_type = "min",
    minimum_ratio = 1,
    minimum_match = 84,
    minimum_relative_cooccurence = 0.95
) {
  start.time <- Sys.time()
  
  # Convert input to data.table
  otutable_dt <- as.data.table(otutable, keep.rownames = TRUE)
  setnames(otutable_dt, "rn", "OTUid")
  
  # Standardize matchlist columns
  colnames(matchlist) <- c("OTUid", "hit", "match")
  matchlist_dt <- as.data.table(matchlist)
  
  # Remove no-hits and self-hits
  matchlist_dt <- matchlist_dt[hit != "*" & OTUid != hit]
  
  # Filter out matches below minimum_match once
  matchlist_dt <- matchlist_dt[match > minimum_match]
  
  # Split the matchlist by OTUid to avoid repeated subsetting
  # result: a named list where each element is the subset of matchlist_dt for that OTUid
  # e.g. matchlist_by_OTUid[["OTU_1"]] is all hits for OTU_1
  matchlist_by_OTUid <- split(matchlist_dt, by = "OTUid", sorted = FALSE, keep.by = FALSE)
  # 'sorted=FALSE' ensures we do not reorder subsets by OTUid. This is a minor optimization.
  
  # Identify sample columns
  stats_cols <- setdiff(names(otutable_dt), "OTUid")
  
  # Create statistics table (total, spread)
  statistics_table <- otutable_dt[, .(
    OTUid = OTUid,
    total = rowSums(.SD),
    spread = rowSums(.SD > 0)
  ), .SDcols = stats_cols]
  
  # Order by spread and total (descending)
  setorder(statistics_table, -spread, -total)
  
  # Initialize parent_id column
  statistics_table[, parent_id := "NA"]
  
  # Reorder otutable_dt to match statistics_table order
  otutable_ordered <- otutable_dt[match(statistics_table$OTUid, otutable_dt$OTUid)]
  
  # Convert the counts (samples) into a matrix for fast numeric access
  samples_matrix <- as.matrix(otutable_ordered[, ..stats_cols])
  rownames(samples_matrix) <- statistics_table$OTUid
  
  # Precompute for each row which sample columns have nonzero counts
  # This helps speed up the co-occurrence intersection and ratio checks
  presence_indices_list <- lapply(seq_len(nrow(samples_matrix)), function(i) {
    which(samples_matrix[i, ] > 0)
  })
  
  # The main loop
  n_rows <- nrow(statistics_table)
  for (line in seq_len(n_rows)) {
    potential_child_id <- statistics_table$OTUid[line]
    # Indices of samples where child is present
    child_presence_cols <- presence_indices_list[[line]]
    
    # Retrieve potential hits from the pre-split list
    potential_hits_dt <- matchlist_by_OTUid[[potential_child_id]]
    # if there's nothing in the list, it means no matches above minimum_match
    if (is.null(potential_hits_dt) || nrow(potential_hits_dt) == 0) {
      # No potential parent found -> it remains its own parent
      statistics_table$parent_id[line] <- potential_child_id
      next
    }
    
    # last_relevant_entry: from original logic
    # we only consider potential parents among the top lines whose 'spread'
    # is >= child's spread. This is how the original code prunes the search space.
    last_relevant_entry <- sum(statistics_table$spread >= statistics_table$spread[line])
    
    # Among those last_relevant_entry lines, find which have OTUid in hits
    # This is a fast way to identify rows from the 'statistics_table' that match the hits
    # potential_hits_dt$hit is a character vector of parent OTU IDs
    hits_vec <- potential_hits_dt$hit
    if (!length(hits_vec)) {
      # If we have no hits, it's still self
      statistics_table$parent_id[line] <- potential_child_id
      next
    }
    
    # Which row(s) in the top 'last_relevant_entry' correspond to those hits
    # i.e. the row index for each "potential parent"
    potential_parents <- which(statistics_table$OTUid[1:last_relevant_entry] %in% hits_vec)
    
    success <- FALSE
    if (length(potential_parents) > 0) {
      # Evaluate each potential parent in ascending row order
      for (line2 in potential_parents) {
        if (!success) {
          # Parent's presence
          parent_presence_cols <- presence_indices_list[[line2]]
          
          # Relative co-occurrence:
          # intersection of child & parent's presence, divided by child's total presence
          cooccur_count <- length(intersect(child_presence_cols, parent_presence_cols))
          relative_cooccurence <- cooccur_count / length(child_presence_cols)
          
          if (relative_cooccurence >= minimum_relative_cooccurence) {
            # Calculate relative abundance ratio
            # Only consider the columns where the child is present
            parent_abundances_at_child_positions <- samples_matrix[line2, child_presence_cols]
            child_abundances_at_child_positions  <- samples_matrix[line,  child_presence_cols]
            
            if (minimum_ratio_type == "avg") {
              relative_abundance <-
                mean(parent_abundances_at_child_positions / child_abundances_at_child_positions)
            } else {
              relative_abundance <-
                min(parent_abundances_at_child_positions / child_abundances_at_child_positions)
            }
            
            if (relative_abundance > minimum_ratio) {
              # If the parent's row is above the child's row, copy parent's parent_id
              # else we consider that parent's own ID
              if (line2 < line) {
                statistics_table$parent_id[line] <- statistics_table$parent_id[line2]
              } else {
                statistics_table$parent_id[line] <- statistics_table$OTUid[line2]
              }
              success <- TRUE
            }
          }
        }
      }
    }
    
    # If no parent was found that satisfies all conditions, set parent_id to self
    if (!success) {
      statistics_table$parent_id[line] <- potential_child_id
    }
  }
  
  # Mark "parent" vs "merged"
  statistics_table[, curated := ifelse(OTUid == parent_id, "parent", "merged")]
  statistics_table[, rank := frank(-total, ties.method = "first")]
  
  # Combine nOTUid with sample counts
  curation_dt <- cbind(nOTUid = statistics_table$parent_id,
                       otutable_ordered[, ..stats_cols])
  
  # Sum rows by parent ID
  curation_table <- curation_dt[, lapply(.SD, sum), by = nOTUid, .SDcols = stats_cols]
  
  # Convert to data.frame
  curation_table_df <- as.data.frame(curation_table)
  row.names(curation_table_df) <- curation_table_df$nOTUid
  curation_table_df$nOTUid <- NULL
  
  # Prepare output
  curated_otus <- unique(statistics_table$parent_id)
  curated_count <- length(curated_otus)
  discarded_otus <- setdiff(statistics_table$OTUid, curated_otus)
  discarded_count <- length(discarded_otus)
  
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  
  # Convert to data.frame for consistency
  statistics_table_df <- as.data.frame(statistics_table)
  
  # Final result
  result <- list(
    curated_table = curation_table_df,
    curated_count = curated_count,
    curated_otus = curated_otus,
    discarded_count = discarded_count,
    discarded_otus = discarded_otus,
    runtime = time.taken,
    minimum_match = minimum_match,
    minimum_relative_cooccurence = minimum_relative_cooccurence,
    otu_map = statistics_table_df,
    original_table = otutable
  )
  
  return(result)
}
