#################################################################################
# switch_Lo2025.R
#
# Skeleton of function
# 1. data wrangling
# 2. identifying all switches that pass QC parameters
# 3. identifying cases of first switch and controls
#################################################################################

#' Antidepressant switching phenotyping (Lo et al, 2025)
#'
#' @description
#' This section summarises functions used in defining switching phenotypes described in:
#'
#' \emph{Antidepressant switching as a proxy phenotype for drug non-response: investigating
#' clinical, demographic and genetic characteristics}
#'
#' For details, please refer to the paper at \url{doi: https://doi.org/10.1016/j.bpsgos.2025.100502}
#'
#'
#' This function performs data wrangling on prescriptions at long format, identifies quality-controlled drug switches, and assembles a case-control dataframe.
#'
#' @param df A data frame containing prescription data. It should include patient identifiers, drug names, drug classes,
#'           and prescription dates.
#' @param id_col A string representing the column name for patient identifiers in the dataframe.
#' @param name_col A string representing the column name for drug names in the dataframe.
#' @param class_col A string representing the column name for drug classes in the dataframe.
#' @param date_col A string representing the column name for prescription dates in the dataframe.
#' @param time_switch Numeric. The maximum allowed time (in days) between prescription dates for a switch to be considered valid.
#' @param nudge_days Numeric. The buffer period (in days) added to remove cases of augmentation & prescription overlaps.
#' @param pre_switch_count Numeric. The maximum number of prescriptions for a drug allowed before the switch.
#' @param post_switch_count Numeric. The maximum number of prescriptions for a drug allowed after the switch.
#' @param total_count Numeric. The maximum total number of prescriptions for a drug across all prescribing journeys.
#' @param identify_controls Logical (TRUE / FALSE). Whether controls need to be isolated.
#' @param prescription_episode_days_control Numeric. The maximum number of days to consider prescriptions
#'   as part of the same prescription episode.
#' @param consecutive_prescriptions_days_control Numeric. The maximum number of days between consecutive prescriptions
#'   to qualify as consistent prescriptions.
#' @param class_filter Optional. A character vector specifying drug classes to filter switchers and non-switchers. Defaults to `NULL` (no filtering).
#'
#' @details
#' The function works in the following steps:
#' \itemize{
#'   \item \code{wrangle_switch_Lo2025}: Concatenating prescriptions from wide to long format.
#'   \item \code{switch_qc_Lo2025}: Identifying all switches, with quality control parameters.
#'   \item \code{case_control_assemble_Lo2025}: Assemble a case-control dataframe.
#' }
#'
#' @return A data frame with the following columns:
#' \itemize{
#'   \item \code{pt_id}: Unique patient identifier.
#'   \item \code{pre_drug}: Drug name prior to the switch or during the prescription episode for controls.
#'   \item \code{pre_class}: Drug class prior to the switch or during the prescription episode for controls.
#'   \item \code{switch_time}: Time (in days) between drug switch events (NA for controls).
#'   \item \code{post_drug}: Drug name after the switch (NA for controls).
#'   \item \code{post_class}: Drug class after the switch (NA for controls).
#'   \item \code{index_date}: Date of the drug switch or prescription episode initiation.
#'   \item \code{cc}: Case-control assignment (\code{case} for switchers, or \code{control} for non-switchers).
#' }
#'
#' @importFrom dplyr mutate ungroup filter rename group_by arrange
#' @importFrom stringr str_count
#' @importFrom tidyr unnest fill
#' @importFrom data.table as.IDate
#'
#' @seealso \code{wrangle_switch_Lo2025}, \code{switch_qc_Lo2025}, \code{case_control_assemble_Lo2025}
#'
#' @examples
#' # Example usage:
#' result <- switch_Lo2025(
#'   df = rx_demo_1,
#'   id_col = "ID", name_col = "drug", class_col = "class", date_col = "start_date",
#'   time_switch = 90, nudge_days = 5,
#'   pre_switch_count = 2, post_switch_count = 2, total_count = 3,
#'   identify_controls = TRUE,
#'   prescription_episode_days_control = 180,
#'   consecutive_prescriptions_days_control = 30)
#' @export

switch_Lo2025 = function(df, id_col, name_col, class_col, date_col,
                         time_switch, nudge_days, pre_switch_count, post_switch_count, total_count,
                         identify_controls,
                         prescription_episode_days_control = NULL,
                         consecutive_prescriptions_days_control = NULL,
                         class_filter = NULL
                         ) {

  # initiate dataframe
  rx_df = df

  # quick checking
  if (sum(is.na(rx_df[,date_col]) ) != 0) { warning("Warning: Please check for NAs on the column for prescription dates (date_col).")}

  # =========================================== #
  # Apply wrangling function
  # Assuming the prescription dataframe is in wide format, where:
  # every row is a prescription
  # =========================================== #
  wrangled_df = wrangle_switch_Lo2025(df = rx_df, id_col = id_col, name_col = name_col, class_col = class_col, date_col = date_col) %>% ungroup()

  # =========================================== #
  ### --- Get all QC-ed switches --- ###
  # =========================================== #
  switch_df = switch_qc_Lo2025(df = wrangled_df, time = time_switch, nudge_days = nudge_days, id_col = id_col,
                        name_col = "drug_names", class_col = "drug_classes", date_col = "drug_dates",
                      pre_switch_count = pre_switch_count, post_switch_count = post_switch_count,
                      total_count = total_count, pheno_col_only = F)
  # =========================================== #
  ### --- Take out cases and controls for first switch, and write as dataframe -- ###
  # =========================================== #
  phenotype = case_control_assemble_Lo2025(df = switch_df, id_col = id_col, name_col = "drug_names", class_col = "drug_classes",
                                    identify_controls = identify_controls,
                                    prescription_episode_days = prescription_episode_days_control,
                                    consecutive_prescriptions_days = consecutive_prescriptions_days_control,
                                    class_filter = class_filter)
  # return output
  return(phenotype)
}


#' Concatenating prescriptions from long to wide format
#'
#' This function concatenates prescription data (in long format) by grouping it by a specified ID column
#' and aggregating drug names, classes, and prescription dates into concatenated strings.
#' It outputs a dataset in wide format, with one row per participant by unique patient identifiers.
#'
#' @param df A data frame containing prescription data.
#' @param id_col A string specifying the column name for unique patient identifiers.
#' @param name_col A string specifying the column name for drug names.
#' @param class_col A string specifying the column name for drug classes.
#' @param date_col A string specifying the column name for prescription dates.
#'
#' @return A data frame grouped by \code{id_col}, with aggregated columns:
#' \itemize{
#'   \item \code{drug_names}: Concatenated drug names (separated by \code{;}).
#'   \item \code{drug_dates}: Concatenated prescription dates (separated by \code{;}).
#'   \item \code{drug_classes}: Concatenated drug classes (separated by \code{;}).
#' }
#'
#' @examples
#' # Example data
#' df <- data.frame(
#'   patient_id = c(1, 1, 2, 2),
#'   drug_name = c("Citalopram", "Citalopram", "Venlafaxine", "Venlafaxine"),
#'   drug_class = c("SSRI", "SSRI", "SNRI", "SNRI"),
#'   date = as.Date(c("2024-01-01", "2024-01-02", "2024-01-03", "2024-01-04"))
#' )
#'
#' # Wrangle data
#' result <- wrangle_switch_Lo2025(df, "patient_id", "drug_name", "drug_class", "date")
#'
#' @importFrom dplyr group_by arrange mutate distinct
#' @importFrom rlang sym
#' @export

wrangle_switch_Lo2025 = function(df, id_col, name_col, class_col, date_col) {
  data = df %>%
    dplyr::group_by(!!rlang::sym(id_col)) %>%
    dplyr::arrange(!!rlang::sym(id_col), !!rlang::sym(date_col), !!rlang::sym(name_col)) %>%
    dplyr::mutate(drug_names = paste(!!rlang::sym(name_col), collapse = ';')) %>%
    dplyr::mutate(drug_dates = paste(!!rlang::sym(date_col), collapse = ';')) %>%
    dplyr::mutate(drug_classes = paste(!!rlang::sym(class_col), collapse = ';')) %>%
    dplyr::distinct(!!rlang::sym(id_col), .keep_all = TRUE)
  return(data)
}

#' Identifying drug-switching events from prescription records
#'
#' This function identifies drug-switching events for a dataframe of prescriptions,
#' after running the wrangle_switch_Lo2025() function.
#'
#' It captures "switching events" based on specified time windows between prescriptions of different drugs.
#'
#' 'Nudge' days are allowed to avoid capturing prescriptions issued very closely
#' (overlapping prescriptions and treatment augmentation) as "switches".
#'
#' Quality control parameters (optional) are also provided on:
#' 1. Number of prescriptions of pre-switch drug before switching date (pre_switch_count)
#' 2. Number of prescriptions of pre-switch drug after switching date (post_switch_count)
#' 3. Total number of prescriptions of pre-switch drug across all prescribing journeys (total_count)
#'
#' These parameters can help to ensure transient exposure of pre-switch drug (with some cross-tapering allowed),
#' and to ensure the pre-switch drug is not prescribed in future episodes again.
#'
#' @importFrom dplyr rename rowwise mutate ungroup select
#' @importFrom rlang sym
#'
#' @param df A data frame containing patient prescription records.
#' @param time Numeric. The maximum allowed time (in days) between prescription dates for a switch to be considered valid.
#' @param nudge_days Numeric. The buffer period (in days) added to remove cases of augmentation & prescription overlaps.
#' @param id_col Character. The column name in `df` representing unique patient identifiers.
#' @param name_col Character. The column name in `df` representing drug names.
#' @param class_col Character. The column name in `df` representing drug classes.
#' @param date_col Character. The column name in `df` representing prescription dates.
#' @param pre_switch_count Numeric. The maximum number of prescriptions for a drug allowed before the switch.
#' @param post_switch_count Numeric. The maximum number of prescriptions for a drug allowed after the switch.
#' @param total_count Numeric. The maximum total number of prescriptions for a drug across all prescribing journeys.
#' @param pheno_col_only Logical. If `TRUE`, the output is limited to phenotype-related columns
#'   (switching details only). Default is `FALSE`.
#'
#' @return A data frame with columns containing original data and QC-passed switching details:
#' \itemize{
#'   \item \code{all_switch_names}: Non-QC drug switch names.
#'   \item \code{all_switch_classes}: Non-QC drug switch classes.
#'   \item \code{all_switch_dates}: Dates of non-QC drug switches.
#'   \item \code{time_bw_switches}: Time (in days) between non-QC drug switches.
#'   \item \code{switch_drug}: Drug names for QC-passed switches.
#'   \item \code{switch_class}: Drug classes for QC-passed switches.
#'   \item \code{switch_date}: Dates of QC-passed switches.
#'   \item \code{switch_time}: Time (in days) of QC-passed switches.
#' }
#'
#' @details
#' The function processes drug prescriptions and identifies switching events where the prescribed drug
#' changes between consecutive records.
#'
#' Switching details are extracted for non-QC events, and a secondary filtering step
#' applies QC criteria to determine eligible switches.
#'
#' The function internally uses the complementary `switch_details()` helper to compute switching parameters.
#' QC checks are applied based on prescription counts before, after, and in total, relative to specified thresholds.
#'
#' If `pheno_col_only` is set to `TRUE`, the output is simplified to include only switching-related columns.
#'
#' @examples
#' # Example dataset
#' dim(rx_demo_1)
#' #[1] 31  8
#'
#' # Apply wrangling function
#' rx_wrangle = wrangle_switch_Lo2025(rx_demo_1, "ID", "drug", "class", "start_date")
#'
#' # Apply switching function
#' rx_switch = switch_qc_Lo2025(
#'   rx_wrangle, time = 90, nudge_days = 5,
#'   id_col = "ID", name_col = "drug_names", class_col = "drug_classes",
#'   date_col = "drug_dates",
#'   pre_switch_count = 2, post_switch_count = 2, total_count = 3
#' )
#' @export

switch_qc_Lo2025 <- function(df, time, nudge_days,
                      id_col, name_col, class_col, date_col,
                      pre_switch_count, post_switch_count, total_count, pheno_col_only = F) {

  # Complimentary function: extracting details of switching
  switch_details <- function(df,
                             time_switch,
                             nudge,
                             pre,
                             post,
                             total,
                             parameter) {

    # time between non-qced drug switches
    time_switch = as.numeric(unlist(strsplit(as.character(df["time_bw_switches"]), ";")))
    # drug name for non-qced drug switches
    # NOTE: this drg_switch contains drug that the patient is taking "BEFORE" the drug switch
    drg_switch = unlist(strsplit(as.character(df["all_switch_names"]), ";"))
    ## drug class for non-qced drug switches
    class_switch = unlist(strsplit(as.character(df["all_switch_classes"]), ";"))
    # all drug prescriptions
    drg_all=unlist(strsplit(as.character(df["name_draft"]), ";"))
    # date of drug prescriptions
    kp_all = unlist(strsplit(as.character(df["date_draft"]), ";"))
    # date of switching
    kp_switch = unlist(strsplit(as.character(df["all_switch_dates"]), ";"))
    # indices of time switching within gap specified
    n= which(time_switch >= nudge_days & time_switch <= (time + nudge_days) )

    # if n contains indices, and so switching matching condition has occurred...
    if(length(n) > 0){
      # drug names for eligible switches
      drg_switch_n = drg_switch[n]
      # drug classes for eligible switches
      class_switch_n = class_switch[n]
      # dates for eligible switches
      kp_switch_n = kp_switch[n]
      # get the names of all drugs with eligible switches
      unique_drugs = unique(drg_switch_n)
      # effective switch information
      qced = NULL

      # ======================================== #
      # Quality control
      # Loop through each drug to check how many times have the drug been prescribed Before
      # ======================================== #

      for(i in 1:length(unique_drugs)){

        # --- Drug name, class and date --- ###
        # get the QC drug of interest
        switch_drug_qc = unique_drugs[i]
        # get drug class for QC drug of interest
        switch_class_qc = min(class_switch_n[which(drg_switch_n == switch_drug_qc)])
        # get switch date for QC drug of interest
        switch_date_qc= min(kp_switch_n[which(drg_switch_n == switch_drug_qc)])

        ### --- Switch time --- ###
        # get first date of switch for the QC drug of interest
        switch_date_qc= min(kp_switch_n[which(drg_switch_n == switch_drug_qc)])  #min of all dates, where corresponding drug is equivalent to for loop instance
        # taking overlapping index on "all_switch_dates" vector and "all_switch_name"
        time_index = which(kp_switch == switch_date_qc)  #index to which the switch occurred
        drug_index = which(drg_switch == switch_drug_qc)
        matched_index = intersect(time_index, drug_index)
        # get the correct time from time_bw_switches vector
        switch_time_qc = time_switch[matched_index]

        ### --- QC parameters --- ###
        # drug prescriptions after the switch date (including repeat prescriptions)
        drug_post= drg_all[which(as.Date(as.character(kp_all)) > as.Date(as.character(switch_date_qc)))]
        # number of times drug was prescribed after switch
        count_post= length(which(drug_post == switch_drug_qc))

        # drug prescriptions before the switch date (including repeat prescriptions)
        drug_pre= drg_all[which(as.Date(as.character(kp_all)) < as.Date(as.character(switch_date_qc)))]
        # number of times drug was prescribed before switch
        count_pre= length(which(drug_pre == switch_drug_qc))

        # total number of prescription for QC drug of interest
        count_all = length(which(drg_all == switch_drug_qc))

        # additional parameter to indicate what information should be extracted
        if (parameter == "drug") { output = switch_drug_qc }
        else if (parameter == "class") { output = switch_class_qc }
        else if (parameter == "date") { output = switch_date_qc }
        else if (parameter == "time") { output = switch_time_qc }

        # add qced information as vector if satisfy criteria
        if(count_pre <= pre_switch_count & count_post <= post_switch_count & count_all <= total_count){
          qced = c(qced, output) } else { NA }
      }

      # after looping on all qc drugs of interest, get that as summary for the particular patient
      paste(qced, collapse = ";") } else { NA }
  }


  data = df %>%
    dplyr::rename(name_draft = !!rlang::sym(name_col), class_draft = !!rlang::sym(class_col), date_draft = !!rlang::sym(date_col)) %>%
    dplyr::rowwise() %>%
    # Non QC-ed switches
    dplyr::mutate(all_switch_names = {drg = unlist(strsplit(as.character(name_draft), ";"))
    n= which(drg[-1] != drg[-length(drg)]) + 1    # identifies changes in drugs/ when different drugs are prescribed
    paste(c(drg[n[1] - 1], drg[n]), collapse = ";") }, # pulls out drug change journey
    all_switch_classes = {drg = unlist(strsplit(as.character(name_draft), ";"))      # pull out all prescriptions (by drug)
    class = unlist(strsplit(as.character(class_draft), ";"))  # pull out all prescriptions (by class)
    n= which(drg[-1] != drg[-length(drg)]) + 1    # identifies changes in drug when different drugs are prescribed
    paste(c(class[n[1] - 1], class[n]), collapse = ";") }, # pulls out respective classes in the drug change journey
    all_switch_dates = {kp = unlist(strsplit(as.character(date_draft), ";"))
    drg = unlist(strsplit(as.character(name_draft), ";"))
    n= which(drg[-1] != drg[-length(drg)]) + 1
    paste(kp[n], collapse = ";") }) %>%
    # Time between switches for all these non-qced switch
    dplyr::mutate(time_bw_switches = {kp_s = unlist(strsplit(as.character(all_switch_dates), ";"))  # switching dates
    kp = unlist(strsplit(as.character(date_draft), ";"))  # drug dates (used to obtain first prescription)
    kp_s= c(kp[1], kp_s)    # initial drug prescription + all switch dates
    paste(round(difftime(kp_s[-1], kp_s[-length(kp_s)], units = "days"), 1), collapse = ";") }) %>%
    dplyr::ungroup()

  # apply complementary functions to each row, to get the details of QC-ed switching
  # drug name
  data$switch_drug = apply(data, 1,
                           function(df) {switch_details(df, time_switch = time, nudge = nudge_days,
                                                        pre = pre_switch_count, post = post_switch_count, total = total_count,
                                                        parameter = "drug")})
  data$switch_class = apply(data, 1,
                            function(df) {switch_details(df, time_switch = time, nudge = nudge_days,
                                                         pre = pre_switch_count, post = post_switch_count, total = total_count,
                                                         parameter = "class")})
  data$switch_date = apply(data, 1,
                           function(df) {switch_details(df, time_switch = time, nudge = nudge_days,
                                                        pre = pre_switch_count, post = post_switch_count, total = total_count,
                                                        parameter = "date")})
  data$switch_time = apply(data, 1,
                           function(df) {switch_details(df, time_switch = time, nudge = nudge_days,
                                                        pre = pre_switch_count, post = post_switch_count, total = total_count,
                                                        parameter = "time")})
  # rename variables
  data = data %>% dplyr::ungroup() %>%
    dplyr::rename(!!rlang::sym(name_col) := name_draft, !!rlang::sym(class_col) := class_draft, !!rlang::sym(date_col) := date_draft)

  # convert "" to NA
  data <- data %>%
    mutate(
      all_switch_names = ifelse(all_switch_names == "", NA, all_switch_names),
      all_switch_classes = ifelse(all_switch_classes == "", NA, all_switch_classes),
      all_switch_dates = ifelse(all_switch_dates == "", NA, all_switch_dates),
      time_bw_switches = ifelse(time_bw_switches == "", NA, time_bw_switches),
      switch_drug = ifelse(switch_drug == "", NA, switch_drug),
      switch_class = ifelse(switch_class == "", NA, switch_class),
      switch_date = ifelse(switch_date == "", NA, switch_date),
      switch_time = ifelse(switch_time == "", NA, switch_time))

  if (pheno_col_only == T) {data = data %>% dplyr::select(!!rlang::sym(id_col), all_switch_names, all_switch_classes,
                                                          all_switch_dates, time_bw_switches,
                                                          switch_drug, switch_class, switch_date, switch_time) }
  return(data)
}

#' Assemble switchers and non-switchers (Lo et al, 2025)
#'
#' This function identifies "switchers" and "non-switchers" based on dataframes produced from \code{switch_qc_Lo2025()}.
#'
#' Details of "switchers" and "non-switchers" identification were described in the paper: \emph{Antidepressant switching as a proxy phenotype for drug non-response: investigating
#' clinical, demographic and genetic characteristics}
#'
#' For details, please refer to the paper at \url{doi: https://doi.org/10.1016/j.bpsgos.2025.100502}
#'
#' The function returns a phenotype dataset
#' with switchers and non-switcher assignments and related details.
#'
#' @param df A data frame containing the input dataset.
#' @param id_col Character. The name of the column representing unique patient identifiers.
#' @param name_col Character. The name of the column containing drug names (separated by ";").
#' @param class_col Character. The name of the column containing drug classes (separated by ";").
#' @param identify_controls Logical (TRUE / FALSE). Whether controls need to be isolated.
#' @param prescription_episode_days Numeric. The maximum number of days to consider prescriptions
#'   as part of the same prescription episode.
#' @param consecutive_prescriptions_days Numeric. The maximum number of days between consecutive prescriptions
#'   to qualify as consistent prescriptions.
#' @param class_filter Optional. A character vector specifying drug classes to filter switchers and non-switchers. Defaults to `NULL` (no filtering).
#'
#' @importFrom dplyr rename mutate filter group_by arrange lag ungroup distinct
#' @importFrom tidyr fill unnest
#' @importFrom stringr str_count str_replace
#' @importFrom data.table as.IDate
#'
#' @details
#' The function works in several steps:
#' \itemize{
#'   \item Case Identification: Identifies individuals who have drug switches meeting the specified criteria
#'   (e.g., first switch drug, class, and associated dates).
#'   \item Control Identification: Identifies individuals without drug switches but who meet prescription consistency
#'   requirements within a prescription episode.
#'   \item Case-Control Assembly: Combines eligible cases and controls into a single data frame, optionally filtered
#'   by specified drug classes.
#' }
#'
#' @return A data frame with the following columns:
#' \itemize{
#'   \item \code{pt_id}: Unique patient identifier.
#'   \item \code{pre_drug}: Drug name prior to the switch or during the prescription episode for controls.
#'   \item \code{pre_class}: Drug class prior to the switch or during the prescription episode for controls.
#'   \item \code{switch_time}: Time (in days) between drug switch events (NA for controls).
#'   \item \code{post_drug}: Drug name after the switch (NA for controls).
#'   \item \code{post_class}: Drug class after the switch (NA for controls).
#'   \item \code{index_date}: Date of the drug switch or prescription episode initiation.
#'   \item \code{cc}: Case-control assignment ("case" or "control").
#' }
#'
#' @examples
#' # Example dataset
#' dim(rx_demo_1)
#' #[1] 31  8
#'
#' # Apply wrangling function
#' rx_wrangle = wrangle_switch_Lo2025(rx_demo_1, "ID", "drug", "class", "start_date")
#'
#' # Apply switching function
#' rx_switch = switch_qc_Lo2025(
#'   rx_wrangle, time = 90, nudge_days = 5,
#'   id_col = "ID", name_col = "drug_names", class_col = "drug_classes",
#'   date_col = "drug_dates",
#'   pre_switch_count = 2, post_switch_count = 2, total_count = 3
#' )
#'
#' # Apply case-control assembling function
#' phenotype = case_control_assemble_Lo2025(df = rx_switch, id_col = "ID",
#'                                   name_col = "drug_names", class_col = "drug_classes",
#'                                   identify_controls = TRUE,
#'                                   prescription_episode_days = 180,
#'                                   consecutive_prescriptions_days = 90, class_filter = "SSRI")
#' @export

case_control_assemble_Lo2025 = function(df, id_col, name_col, class_col,
                                 prescription_episode_days = NULL, consecutive_prescriptions_days = NULL,
                                 identify_controls,
                                 class_filter = NULL) {

  # checking arguments
  if (identify_controls == TRUE & (is.null(prescription_episode_days) | is.null(consecutive_prescriptions_days)))
    {warning("Please specify prescription_episode_days and consecutive_prescriptions_days to identify controls.")}
  # next drug - get next drug for switchers
  # Function to get drug names and class after switch
  next_drug <- function(df, parameter) {

    # containing all switches (drug names, drug classes, dates of switch)
    # Note: These "drug names" and "drug classes" are drugs BEFORE switch
    rx_names = unlist(strsplit(df[["all_switch_names"]], ";"))
    rx_class = unlist(strsplit(df[["all_switch_classes"]], ";"))
    rx_dates = unlist(strsplit(df[["all_switch_dates"]], ";"))

    # identify drug for **first switch**
    first_switch_drug = unlist(strsplit(df[["switch_drug"]], ";"))[1]
    first_switch_date = unlist(strsplit(df[["switch_date"]], ";"))[1]

    if (!is.na(first_switch_drug)) {
      # take out the respective date from switches, but need to match based on correct drug name and date
      drug_index = which(rx_names == first_switch_drug)     # drug index match
      date_index = which(rx_dates == first_switch_date)     # date index match
      match_index = intersect(drug_index, date_index)
      # with the matched index, take out the next drug and class
      next_drug = rx_names[match_index + 1]
      next_class = rx_class[match_index + 1]
      # also get the index date
      # index for the first prescription of the drug that matches with first switch drug
      i = which(unlist(strsplit(df[["drug_names"]], ";")) == first_switch_drug)[1]
      # get the date for that prescription
      date = unlist(strsplit(df[["drug_dates"]], ";"))[i]
    } else {
      # if switch_drug is NA, just code them as NA (because there are no effective switch)
      next_drug = NA
      next_class = NA
      date = NA}
    # get outcomes as depending on parameter entered
    if (parameter == "name") {
      outcome = next_drug
    } else if (parameter == "class") {
      outcome = next_class
    } else if (parameter == "index_date") {
      outcome = date}
    return(outcome)}

  # isolate the switchers
  case <- df %>%
    dplyr::rename(pt_id = !!rlang::sym(id_col)) %>%
    dplyr::mutate(n_switch = ifelse(is.na(switch_drug), 0 ,stringr::str_count(switch_drug, ";")+1)) %>%
    dplyr::mutate(switch_date = as.character(switch_date)) %>%
    dplyr::mutate(cc = ifelse(n_switch == 0, "no_switch", "case")) %>%
    dplyr::filter(cc == "case")

  # checking to see if switchers identified
  if (nrow(case) == 0) {warning("No switchers identified.")}

  # details for switchers
  case_details <- data.frame(pt_id = case$pt_id,
                             # drug name before switch
                             pre_drug = apply(case, 1, function(df) {unlist(strsplit(df["switch_drug"], ";"))[1]}),
                             # drug class before switch
                             pre_class = apply(case, 1, function(df) {unlist(strsplit(df["switch_class"], ";"))[1]}),
                             # time between switch
                             switch_time = apply(case, 1, function(df) {unlist(strsplit(df["switch_time"], ";"))[1]}),
                             # drug name after switch
                             post_drug = unlist(lapply(apply(case, 1, function(df) {next_drug(df, parameter = "name")}), "[", 1)),
                             # drug class after switch
                             post_class = unlist(lapply(apply(case, 1, function(df) {next_drug(df, parameter = "class")}), "[", 1)),
                             # index date
                             index_date = unlist(lapply(apply(case, 1, function(df) {next_drug(df, parameter = "index_date")}), "[", 1))) %>%
    dplyr::mutate(index_date = data.table::as.IDate(index_date))

  # potential controls
  if (identify_controls == TRUE) {

    non_case = df %>% dplyr::rename(pt_id = !!rlang::sym(id_col)) %>%
      dplyr::filter(!pt_id %in% case_details$pt_id) %>%
      dplyr::group_by(pt_id) %>%
      dplyr::mutate(names = strsplit(drug_names, ";"), class = strsplit(drug_classes, ";"), dates = strsplit(drug_dates, ";")) %>%
      tidyr::unnest(cols = c(names,class,dates)) %>%
      dplyr::mutate(dates = data.table::as.IDate(dates))

    # Grouping of prescription episode
    non_case = non_case %>%
      dplyr::group_by(pt_id, !!rlang::sym(name_col)) %>%
      dplyr::arrange(dates) %>%
      dplyr::mutate(prev_drug_date = dplyr::lag(dates, n = 1, default = NA)) %>%
      dplyr::mutate(diff_days_drug = as.numeric(difftime(dates, prev_drug_date, units = "days"))) %>%
      dplyr::mutate(prescription_episode = ifelse(diff_days_drug > prescription_episode_days, seq_along(diff_days_drug), 1)) %>%
      tidyr::fill(prescription_episode) %>%
      dplyr::mutate(prescription_episode = ifelse(is.na(prescription_episode), 1, prescription_episode)) %>%
      dplyr::group_by(pt_id, !!rlang::sym(name_col), prescription_episode) %>%
      dplyr::mutate(n_drug_prescr = sum(!duplicated(dates))) %>% dplyr::ungroup()

    # Look at consecutive prescriptions
    non_case = non_case %>%
      dplyr::group_by(pt_id, !!rlang::sym(name_col), prescription_episode) %>%
      dplyr::arrange(dates) %>%
      dplyr::mutate(time_crit_n = ifelse(diff_days_drug < consecutive_prescriptions_days, 1, 0)) %>%
      dplyr::mutate(time_crit_tot = sum(time_crit_n, na.rm = TRUE)) %>%
      dplyr::ungroup()

    # Assess eligibility of controls
    controls = non_case %>%
      dplyr::group_by(pt_id) %>%
      dplyr::arrange(dates, !!rlang::sym(name_col)) %>%
      dplyr::mutate(no_switch = ifelse((n_drug_prescr >= 3 & time_crit_tot >= 1), 1, 0)) %>%
      dplyr::filter(no_switch == 1) %>%
      dplyr::group_by(pt_id, !!rlang::sym(name_col), prescription_episode) %>%
      dplyr::reframe(control_drug = unique(!!rlang::sym(name_col)),
                     control_class = unique(!!rlang::sym(class_col)),
                     date = min(as.Date(dates))) %>%
      dplyr::group_by(pt_id, control_class) %>%
      dplyr::filter(as.Date(date) == min(as.Date(date))) %>%
      dplyr::distinct(pt_id, control_class, date, .keep_all = TRUE) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(pre_drug = control_drug,
                    pre_class = control_class,
                    switch_time = NA,
                    post_drug = NA, post_class = NA,
                    index_date = date) %>%
      dplyr::select(pt_id, pre_drug, pre_class, switch_time, post_drug, post_class, index_date)


  }

  # assemble both cases and controls
  if (!is.null(class_filter)) {
    elig_case = case_details %>% dplyr::filter(grepl(paste(class_filter, collapse="|"), pre_class))
    if (identify_controls == TRUE) {elig_controls = controls %>% dplyr::filter(grepl(paste(class_filter, collapse="|"), pre_class)) }
  } else { elig_case = case_details
    if (identify_controls == TRUE) {elig_controls = controls} }
  # assign case control status
  if (identify_controls == TRUE) {pheno_df = rbind(elig_case %>% dplyr::mutate(cc = "case"), elig_controls %>% dplyr::mutate(cc = "control"))}
  if (identify_controls == FALSE) {pheno_df = elig_case %>% dplyr::mutate(cc = "case")}
  # reorder columns
  pheno_df = pheno_df %>% select(pt_id, cc, pre_drug, pre_class, switch_time, post_drug, post_class, index_date)
  return(pheno_df)
}

