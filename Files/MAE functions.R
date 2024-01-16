#### functions for reading byonic files####
# By Charley Cai for BIOL 448
# Supervised by Teesha Baker and Leonard Foster
# 2022-01-18

#' Function for parsing xlxs files. Starting and end positions are reported in 1-based, fully closed.
#' Keeps only distinct peptides in a fraction.
#' @param path The file path to an byonic xlsx file
#' @param exclude A string or vector to filter out proteins with strin in Protein name
#' @param min_pep Minimum peptide length
#' @param max_pep Maximum peptide length
#' @param min_z Minimum z score
#' @param cols Select columns
#' @return A dataframe

read_byonic <- function(path,
                        exclude = c("BOVIN", "Rev", "contam"),
                        min_pep = 9,
                        max_pep = 20,
                        min_z = 2,
                        cols = c("Accession", "Peptide.Sequence", "Description", "Start", "End", "Log.Prob")) {

  byonic.colnames <- c("Query", "Protein.Rank", "Peptide", "Glycans", "Modifications",
                       "Observed.m.z", "z", "Observed.M+H", "Calc.mass", "Off.by.x.error",
                       "Mass.error.ppm", "Start", "Cleavage", "Score", "Delta",
                       "Delta.mod", "Log.Prob", "#.of.unique.peptides", "Protein.Name",
                       "Protein.DB.number", "Mobility", "Scan.#", "Scan.Time")

  read_xlsx(path, sheet = "Spectra", col_names = byonic.colnames, skip = 1) %>%
    mutate(Peptide.Sequence = str_extract(str_remove_all(`Peptide`, "\\[.*?\\]"), "(?<=\\.).*(?=\\.)"),
           End = Start + str_length(Peptide.Sequence) - 1,
           Accession = str_extract(Protein.Name, "(?<=\\|)[^|]+(?=\\|)"),
           Description = Protein.Name)%>%
    filter(str_detect(Protein.Name, paste(exclude, collapse = '|'), negate = T)
           &
             str_length(Peptide.Sequence) >= min_pep
           &
             str_length(Peptide.Sequence) <= max_pep
           &
             Log.Prob >= min_z) %>%
    distinct(Peptide.Sequence, .keep_all = T) %>%
    select(one_of(cols))
}

#' takes a list of file paths to byonic files and merges them, taking distinct peptides
#'
#' @param fraction_paths A list of paths different fraction byonic files within the same replicate
#'
#' @param merge_fraction_dups A boolean indicating whether duplicate peptides between fractions should be merged
#'
read_and_merge_distinct_peptides_from_fractions <- function(fraction_paths, merge_fraction_dups = T, ...) {
  peptides <- map(fraction_paths, function(path) {
    tryCatch({
      read_byonic(path, ...)
    }, error = function(e) { })
  })

  peptides <- peptides %>% bind_rows()
  if (merge_fraction_dups) {
    return(peptides %>% distinct(Peptide.Sequence, .keep_all = T))
  } else {
    return(peptides)
  }
}

#' takes list of tibbles with replicates, each with fractions, calls read_and_merge for each fraction,
#' then rbinds replicates together
#'
#' @param condition A tibble containing a list of file paths as first column for a single byonic condition and a
#' column "rep" indicating distinct replicates

bind_replicates_countrep <- function(condition, min_percentage, out_path, ...) {
  reps <- unique(condition$rep)
  res <- vector("list", length(reps)) %>% setNames(reps)
  namesl = list()
  for (r in reps) {
    files <- filter(condition, rep == r)[[1]]
    res[[r]] <- read_and_merge_distinct_peptides_from_fractions(files)
    namesl[[r]] <- basename(files)[1]
  }
  
  total_reps <- length(reps)
  min_rep_count <- ceiling(total_reps * min_percentage)
  
  peptide_counts <- table(unlist(lapply(res, `[[`, "Peptide.Sequence")))
  valid_peptides <- names(peptide_counts[peptide_counts >= min_rep_count])
  
  peptide_counts2 <- data.frame(peptide_counts)
  peptide_counts2$"Peptide.Sequence" = peptide_counts2$Var1
  
  for (r in reps) {
    write.csv(res[[r]], paste0(out_path,"\\",gsub("\\..*","",namesl[[r]]),".csv"))
    res[[r]] = merge(x=res[[r]], y=peptide_counts2, by="Peptide.Sequence")
    res[[r]] = res[[r]][c("Accession","Peptide.Sequence","Description","Start","End","Log.Prob","Freq")]
  }
  
  return(bind_rows(res)%>% distinct(Peptide.Sequence, .keep_all = T) %>% filter(Peptide.Sequence %in% valid_peptides))
}



#' Find epitopes from list of peptides. Queenie's version
#'
#' For each Protein, group overlapping peptides reads within protein sequence and report the union
#' of the start and stop positions and count how many peptides contributed.
#'
#'The union of two sets is a new set that contains all of the elements that are in at least one of the two sets. 
#'The union is written as A∪B or “A or B”. 
#'The intersection of two sets is a new set that contains all of the elements that are in both sets. 
#'The intersection is written as A∩B or “A and B”.
#'
#' @param peptides A dataframe with columns of protein name, peptide, start, stop. Position indexing
#' is 1-based fully closed. Must contain columns names of "Accession", "Peptide.Sequence", "Start" "End"
#'
#' @return A dataframe of epitopes with protein name, start, stop, sequence, and number of occurances
find_epitopes_union <- function(peptides) {
  sorted <- arrange(peptides, Accession, Start, End)
  n_peptides <- nrow(sorted)
  result <- list(
    Accession = character(n_peptides),
    Description = character(n_peptides),
    Epitope = character(n_peptides),
    Start = double(n_peptides),
    End = double(n_peptides),
    Count = rep(1, n_peptides),
    Pep = character(n_peptides)
  )
  current_protein <- ""
  n_epitopes <- 0
  
  for (i in seq_len(n_peptides)) {
    # start new epitope if moved to a new protein or next peptide does not overlap with the previous epitope end
    if (sorted[[i, "Accession"]] != current_protein || sorted[[i, "Start"]] > result$End[n_epitopes]) {
      n_epitopes <- n_epitopes + 1
      current_protein <- sorted[[i, "Accession"]]
      
      result$Accession[[n_epitopes]] <- current_protein
      result$Description[[n_epitopes]] <- sorted[[i, "Description"]]
      result$Epitope[[n_epitopes]] <- sorted[[i, "Peptide.Sequence"]]
      result$Start[[n_epitopes]] <- sorted[[i, "Start"]]
      result$End[[n_epitopes]] <- sorted[[i, "End"]]
      result$Pep[[n_epitopes]] = sorted[[i, "Peptide.Sequence"]]
    } else {
      # set epitope end to the end of the current peptide if the current epitope ends further
      if (result$End[n_epitopes] < sorted[[i, "End"]]) {
        longer <- sorted[[i, "End"]] - result$End[n_epitopes]
        result$End[n_epitopes] <- sorted[[i, "End"]]
        result$Epitope[n_epitopes] <- paste0(result$Epitope[n_epitopes], str_sub(sorted[[i, "Peptide.Sequence"]], -longer))
        result$Pep[n_epitopes] <- paste0(result$Pep[n_epitopes],',' ,sorted[[i, "Peptide.Sequence"]])
        # increment epitope count
        result$Count[n_epitopes] <- result$Count[n_epitopes] + 1
      }
    }
  }
  
  result <- result %>%
    as_tibble() %>%
    slice_head(n = n_epitopes) %>%
    mutate(Length = End - Start + 1)
}

# Functions written by Teesha Baker to expand the script originally written by Charley Cai####
# Function to create FASTA protein file
create_fasta_file <- function(tibble, file_name) {
  # Open the file in write mode
  file <- file(file_name, "w")
  
  # Write each row as a FASTA entry
  for (i in 1:nrow(tibble)) {
    header <- paste0("Pep_", i)
    sequence <- tibble$Epitope[i]
    
    # Write the header line
    writeLines(paste0(">", header), file)
    
    # Write the sequence
    writeLines(sequence, file)
  }
  
  # Close the file
  close(file)
}


# Function to add empty column to a tibble
add_columns <- function(tibble_data) {
  tibble_data$Allele <- rep(NA, nrow(tibble_data))
  tibble_data$Matched.Binding.Core <- rep(NA, nrow(tibble_data))
  return(tibble_data)
}

# Add empty column to each tibble in the list using a for loop
#for (i in seq_along(cond_epitopes_union)) {
  #cond_epitopes_union[[i]] <- add_columns(cond_epitopes_union[[i]])
#}


#function to match binders to epitopes and add binding core to epitopes file 
add_matching_columns <- function(tibble_data, binders) {
  match_row <- NULL
  tibble_data$Allele <- NA
  tibble_data$Matched.Binding.Core <- NA
  for (i in 1:nrow(tibble_data)) {
    epitope <- tibble_data$Epitope[i]
    match_row <- binders[str_detect(epitope, binders$Peptide.Sequence), ]
    if (nrow(match_row) > 0) {
      match_row <- match_row[match_row$Rank == min(match_row$Rank), ]
      tibble_data$Allele[i] <- match_row$Allele
      tibble_data$Matched.Binding.Core[i] <- match_row$Peptide.Sequence
    }
  }
  return(tibble_data)
}




# Function to filter rows based on Count or Intensity columns
filter_rows <- function(tibble) {
  count_intensity_cols <- grep("Count|Intensity", names(tibble), value = TRUE)
  filtered_tibble <- tibble[rowSums(tibble[count_intensity_cols] != 0, na.rm = TRUE) > 0, ]
  return(filtered_tibble)
}

# Function to filter rows based on occurrence count in "Count" columns
filter_rows_occurrence <- function(tibble, min_occurrence) {
  count_cols <- grep("Spectral.Count", names(tibble), value = TRUE)
  occurrence_count <- rowSums(tibble[count_cols] > 0, na.rm = TRUE)
  threshold <- length(count_cols) * min_occurrence
  filtered_tibble <- tibble[occurrence_count > threshold, ]
  return(filtered_tibble)
}


# Function to remove unwanted columns based on prefixes
remove_unwanted_columns <- function(data, prefixes) {
  data %>%
    dplyr::select(-matches(paste0("^(", paste(prefixes, collapse = "|"), ")")))
}



# Create a function to filter tibbles based on Peptide.Sequence length
filter_tibble_by_length <- function(tibble, operator, threshold) {
  if (operator == "<") {
    filtered_tibble <- tibble[nchar(tibble$Peptide.Sequence) < threshold, ]
  } else if (operator == ">") {
    filtered_tibble <- tibble[nchar(tibble$Peptide.Sequence) > threshold, ]
  } else {
    stop("Invalid comparison operator. Use '<' or '>'.")
  }
  return(filtered_tibble)
}

convert_pdf_to_png <- function(pdf_file) {
  #PDF filename
  pdf_filename <- file.path(pdf_file)
  # Extract the PDF file name without the extension
  noext <- tools::file_path_sans_ext(pdf_file)
  
  # Read the PDF file and convert to images
  pdf_pages <- magick::image_read(pdf_file)
  
  # Save each page as a separate PNG file
  for (i in seq_along(pdf_pages)) {
    png_filename <- paste0(noext, "-001.png")
    magick::image_write(image = pdf_pages, path = png_filename, format = "png")
    cat("Converted:", pdf_file, "to", "png", "\n")
  }
}

# function to output csv files for further analysis (venn diagrams and distribution plots) for fragpipe
write_frag_pep_csvs = function(condition,out_path)
{
  
  file = condition
  num_cols = length(names(file))
  num_conds = (num_cols-16)/3
  
  nc = (num_cols - num_conds)+1
  test_sub = file[,c(1:16,nc:num_cols)]
  
  for(i in 1:num_conds){
    print(i)
    namep1 = names(test_sub[,17])
    namep1 = paste0(gsub("\\..*","",namep1),"Rep_",i)
    out_sub = test_sub[,c(1:16,16+i)]
    out_sub = subset(out_sub,out_sub[,17]>0)
    out_sub = out_sub[,c(1:16)]
    write.csv(out_sub,paste0(out_path,"\\",namep1,".csv"))
  }
}

count_byonic_spectra <- function(path,
                                 exclude = c( "Rev"),
                                 cols = c("Accession", "Peptide.Sequence", "Description", "Start", "End", "Log.Prob")) {
  
  byonic.colnames <- c("Query", "Protein.Rank", "Peptide", "Glycans", "Modifications",
                       "Observed.m.z", "z", "Observed.M+H", "Calc.mass", "Off.by.x.error",
                       "Mass.error.ppm", "Start", "Cleavage", "Score", "Delta",
                       "Delta.mod", "Log.Prob", "#.of.unique.peptides", "Protein.Name",
                       "Protein.DB.number", "Mobility", "Scan.#", "Scan.Time")
  
  dt = read_xlsx(path, sheet = "Spectra", col_names = byonic.colnames, skip = 1) %>%
    mutate(Peptide.Sequence = str_extract(str_remove_all(`Peptide`, "\\[.*?\\]"), "(?<=\\.).*(?=\\.)"),
           End = Start + str_length(Peptide.Sequence) - 1,
           Accession = str_extract(Protein.Name, "(?<=\\|)[^|]+(?=\\|)"),
           Description = Protein.Name)%>%
    filter(str_detect(Protein.Name, paste(exclude, collapse = '|'), negate = T)
    ) %>%
    select(one_of(cols))
  nrow(dt)
}


reorder_columns_alphabetically <- function(df) {
  # Extract the column names from the data frame
  col_names <- names(df)
  
  # Filter the column names to exclude 'position' and 'sequence' columns
  filtered_col_names <- col_names[!(col_names %in% c("position", "sequence"))]
  
  # Sort the remaining column names alphabetically
  sorted_col_names <- sort(filtered_col_names)
  
  # Reorganize the columns in the data frame
  df <- df[c("position", "sequence", sorted_col_names)]
  
  return(df)
}

# Define a function to read in FIMO files from MEMESuite
read_FIMO_files <- function(base_dir, file_name) {
  # Get a list of all subdirectories under the base directory
  subdirectories <- list.dirs(base_dir, full.names = TRUE, recursive = TRUE)
  # Filter the subdirectories to include only those containing "SARS-COV-2"
  subdirectories_subset <- subset(subdirectories, grepl(protein, subdirectories))
  # Initialize an empty list to store data frames
  data_list <- list()
  
  # Loop through the subdirectories
  for (subdir in subdirectories_subset) {
    # Construct the full file path
    file_path <- file.path(subdir, file_name)
    # Read the file into a data frame (adjust read.table() or read.csv() as needed)
    data <- readLines(file_path)  # Read the file as lines
    
    # Check if the file has at least 4 lines
    if (length(data) > 4) {
      data.table <- read.table(text = data, header = TRUE, sep = "\t")
      data_list[[file_path]] <- data.table
    } else {
      empty_data <- data.frame()
      data_list[[file_path]] <- empty_data
    }
  }
  
  return(data_list)
}

# Define a function to write a FASTA file (create_fasta_file but created a variable for the sequence column)
create_fasta_file_new <- function(tibble, column_header, file_name) {
  # Open the file in write mode
  file <- file(file_name, "w")
  
  # Write each row as a FASTA entry
  for (i in 1:nrow(tibble)) {
    header <- paste0("Pep_", i)
    sequence <- as.character(tibble[[column_header]][i])      
    # Write the header line
    writeLines(paste0(">", header), file)
    
    # Write the sequence
    writeLines(sequence, file)
  }
  
  # Close the file
  close(file)
}

# Define a function to find the union of peptides 
find_epitopes_union_new <- function(peptides_table, accession_col, start_col, end_col, peptide_sequence_col) {
  sorted <- arrange(peptides_table, !!sym(accession_col), !!sym(start_col), desc(!!sym(end_col)))
  n_peptides <- nrow(sorted)
  result <- list(
    Accession = character(n_peptides),
    Epitope = character(n_peptides),
    Start = double(n_peptides),
    End = double(n_peptides),
    Count = rep(1, n_peptides),
    Peptides = character(n_peptides)
  )
  current_protein <- ""
  n_epitopes=0
  i=1
  for (i in seq_len(n_peptides)) {
    # start a new epitope if moved to a new protein or the next peptide does not overlap with the previous epitope end
    if (sorted[[i, accession_col]] != current_protein || sorted[[i, start_col]] > result$End[n_epitopes]) {
      n_epitopes <- n_epitopes + 1
      current_protein <- sorted[i, accession_col]
      result$Accession[[n_epitopes]] <- current_protein
      result$Epitope[[n_epitopes]] <- sorted[[i, peptide_sequence_col]]
      result$Start[[n_epitopes]] <- sorted[[i, start_col]]
      result$End[[n_epitopes]] <- sorted[[i, end_col]]
      result$Peptides[[n_epitopes]] = sorted[[i, peptide_sequence_col]]
    } else {
      # set epitope end to the end of the current peptide if the current epitope ends further
      if (result$End[n_epitopes] < sorted[[i, end_col]]) {
        longer <- sorted[[i, end_col]] - result$End[n_epitopes]
        result$End[n_epitopes] <- sorted[[i, end_col]]
        result$Epitope[n_epitopes] <- paste0(result$Epitope[n_epitopes], str_sub(sorted[[i, peptide_sequence_col]], -longer))
        result$Peptides[n_epitopes] <- paste0(result$Peptides[n_epitopes],',' ,sorted[[i, peptide_sequence_col]])
        # increment epitope count
        result$Count[n_epitopes] <- result$Count[n_epitopes] + 1
      } else {
        # if both conditions are false, paste the current entry into the new result table
        n_epitopes <- n_epitopes + 1
        current_protein <- sorted[i, accession_col]
        result$Accession[[n_epitopes]] <- current_protein
        result$Epitope[[n_epitopes]] <- sorted[[i, peptide_sequence_col]]
        result$Start[[n_epitopes]] <- sorted[[i, start_col]]
        result$End[[n_epitopes]] <- sorted[[i, end_col]]
        result$Peptides[[n_epitopes]] = sorted[[i, peptide_sequence_col]]
      }
    }
  }
  
  result <- result %>%
    as_tibble() %>%
    mutate(Count = result$Count, Length = End - Start + 1) %>%
    slice_head(n = n_epitopes) %>%
    select(Accession, Epitope, Start, End, Count, Length, Peptides)
}

