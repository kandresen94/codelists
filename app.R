###############################################################################
# Shiny App Integration
# Dictionary: optimised for AURUM (Doesnt work for others yet)
#
# Definition: Shiny app to create codelists:
# This code 
# Creates a list of included codes based on the inclusion terms provided by the user.
# Removes unwanted from the included list based on the exclusion terms provided by the user.
# Compares the included list to a comparison codelist (if provided).
# Re-applies the exclusion terms in case the comparison codelist contains unwanted terms.
# Allows download of the final codelist and all of the exclusion terms for clinical review
#
# Author: Kirsty Andresen
# Date created: 01 November 2024
#
# UPDATE HISTORY:
# Date update 1: 
# Changes: 
###############################################################################

# Load required libraries
library(shiny)
library(dplyr)
library(stringr)
library(haven)

# Path for loading library
dict_path <- ("H:/")
# Path where saving codelists
save_path <- ("H:/")

# UI for Shiny App
ui <- fluidPage(
  titlePanel("Disease Codelist Generator for Aurum"),
  sidebarLayout(
    sidebarPanel(
      textInput("disease_name", "Enter Disease Name:"),
      textInput("inclusion_terms", "Enter Inclusion Terms (comma separated):"),
      textInput("exclusion_terms", "Enter Exclusion Terms (comma separated):"),
      fileInput("dictionary_file", "Upload Medical Dictionary (.txt):", accept = c(".txt")),
      fileInput("original_codelist", "Upload Comparison Codelist (.dta or .csv):", accept = c(".dta", ".csv")),
      actionButton("run_analysis", "Run"),
      downloadButton("download_final_list", "Download Final Codelist (.csv)"),
      downloadButton("download_exclusion_list", "Download Exclusion list (.csv)")
      
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Final List", DT::dataTableOutput("final_list_output")),
        tabPanel("All Excluded Codes", DT::dataTableOutput("all_excluded_codes_output"))
      )
    )
  )
)

# Server logic for Shiny App
server <- function(input, output) {
  final_list <- reactiveVal()
  all_excluded_codes <- reactiveVal()
  
  observeEvent(input$run_analysis, {
    # Load comparison codelist if provided
    original_codelist <- NULL
    if (!is.null(input$original_codelist)) {
      file_ext <- tools::file_ext(input$original_codelist$name)
      if (file_ext == 'dta') {
        original_codelist <- read_dta(input$original_codelist$datapath)
      } else if (file_ext == 'csv') {
        original_codelist <- read.csv(input$original_codelist$datapath, stringsAsFactors = FALSE)
      }
    }
    
    # Step 1: Load medical dictionary and prep for search
    medical <- NULL
    if (!is.null(input$dictionary_file)) {
      dict_ext <- tools::file_ext(input$dictionary_file$name)
      if (dict_ext == 'txt') {
        medical <- read.table(input$dictionary_file$datapath, header = TRUE, fill = TRUE, sep = "\t", encoding = "UTF-8", quote = "", stringsAsFactors = FALSE)
      }
    }
    if (!is.null(medical)) {
      medical <- medical %>% mutate(across(where(is.character), ~tolower(iconv(., from = "UTF-8", to = "ASCII//TRANSLIT", sub = ""))))
      medical <- medical %>% rename(medcodeid = MedCodeId, term = Term, readcode = CleansedReadCode, snomedcode = SnomedCTConceptId)
    }
    
    # Step 2: Define inclusion and exclusion terms from user input
    interms <- str_split(input$inclusion_terms, ",\\s*")[[1]]
    exterms <- str_split(input$exclusion_terms, ",\\s*")[[1]]
    incodes <- c() # SNOMED or Read
    excodes <- c() # SNOMED or Read
    
    # Step 3: Word search of CPRD medical code dictionary
    marker <- rep(0, nrow(medical))
    for (term in interms) {
      marker[str_detect(medical$term, regex(term, ignore_case = TRUE))] <- 1
    }
    
    for (code in incodes) {
      marker[medical$medcodeid %in% code] <- 1
      marker[medical$snomedcode %in% code] <- 1
    }
    
    medical <- medical %>% mutate(marker = marker)
    
    # Step 4: Sort and drop terms not captured by search terms
    medical <- medical %>% arrange(desc(marker), medcodeid)
    inclusion_list <- medical %>% filter(marker == 1) %>% select(-marker)
    
    # Step 5: Exclude unwanted terms to make search more specific
    first_draft <- inclusion_list
    for (word in exterms) {
      first_draft <- first_draft %>% filter(!str_detect(term, regex(word, ignore_case = TRUE)))
    }
    
    for (code in excodes) {
      first_draft <- first_draft %>% filter(!medcodeid %in% code, !snomedcode %in% code)
    }
    
    excluded_codes_first_draft <- inclusion_list %>% anti_join(first_draft, by = c("medcodeid", "term")) %>%
      mutate(source = "from dictionary") %>%
      select(medcodeid, term, source)
    
    if (!is.null(original_codelist)) {
      original_codelist <- original_codelist %>% mutate(medcodeid = as.numeric(medcodeid))
      
      not_included_in_update <- original_codelist %>% 
        select(medcodeid, term) %>%
        anti_join(first_draft, by = c("medcodeid", "term")) %>%
        arrange(medcodeid) %>%
        mutate(source = "original_codelist")
      
      new_codes_from_update <- first_draft %>% 
        select(medcodeid, term) %>%
        anti_join(original_codelist, by = c("medcodeid", "term")) %>%
        arrange(medcodeid) %>%
        mutate(source = "new")
      
      # append together
      comparison <- bind_rows(not_included_in_update, new_codes_from_update)
      
      all <- first_draft %>% select(medcodeid, term) %>% 
        full_join(comparison, by = c("medcodeid", "term")) %>%
        mutate(source = ifelse(is.na(source), "match", source)) %>%
        arrange(term) %>%
        select(medcodeid, term, source)
    } else {
      all <- first_draft %>% mutate(source = "match")
    }
    
    # Step 6: Re-exclude unwanted terms (removes codes we want to purposefully exclude from the original codelist)
    final_list <- all
    for (word in exterms) {
      final_list <- final_list %>% filter(!str_detect(term, regex(word, ignore_case = TRUE)))
    }
    
    for (code in excodes) {
      final_list <- final_list %>% filter(!medcodeid %in% code, !snomedcode %in% code)
    }
    
    excluded_codes_final_list <- all %>% anti_join(final_list, by = c("medcodeid", "term")) %>%
      mutate(source = "from codelist") %>%
      select(medcodeid, term, source)
    
    # Full list of excluded codes
    all_excluded_codes(all_excluded_codes <- bind_rows(excluded_codes_first_draft, excluded_codes_final_list))
    
    final_list(final_list %>% select(medcodeid, term, source))
    
    # Output final list and all excluded codes
    output$final_list_output <- DT::renderDataTable({
      DT::datatable(final_list(), options = list(pageLength = 30)) %>%
        DT::formatStyle(
          'source',
          target = 'row',
          backgroundColor = DT::styleEqual(c("match", "new"), c("darkgreen", "lightgreen"))
        )
    })
    
    output$all_excluded_codes_output <- DT::renderDataTable({
      DT::datatable(all_excluded_codes(), options = list(pageLength = 30)) %>%
        DT::formatStyle(
          'source',
          target = 'row',
          backgroundColor = DT::styleEqual(c("from codelist", "from dictionary"), c("darkred", "lightcoral"))
        )
    })
  })
  
  # Download final codelist
  output$download_final_list <- downloadHandler(
    filename = function() {
      paste0("cr_codelist_", input$disease_name, ".csv")
    },
    content = function(file) {
      write.csv(final_list(), file, row.names = FALSE)
    }
  )
  
  output$download_exclusion_list <- downloadHandler(
    filename = function() {
      paste0("cr_codelist_", input$disease_name, "_exclusion.csv")
    },
    content = function(file) {
      write.csv(all_excluded_codes(), file, row.names = FALSE)
    }
  )
}

# Run the application
options(shiny.maxRequestSize = 100*1024^2) # Increase maximum upload size to 100MB
shinyApp(ui = ui, server = server)
