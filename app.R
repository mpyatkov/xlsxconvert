library(shiny)
library(shinythemes)
library(readxl)
library(readr)
library(dplyr)
library(purrr)
library(stringr)

# --- App2 auxiliary validation function ---
file_validation <- function(file_path) {
  output_list <- list(data = NULL, header = NULL, status = FALSE, error_message = "")
  
  params <- tryCatch({
    read_xlsx(file_path, range = cell_rows(1:2), sheet = 1, col_names = TRUE)
  }, error = function(e) {
    # This will be caught and an error message generated
    NULL
  })
  
  if (is.null(params)) {
    output_list$error_message <- "Error reading the HEADER from the file. Ensure it has 2 rows."
    return(output_list)
  }
  
  data <- tryCatch({
    read_xlsx(file_path, skip = 3, sheet = 1)
  }, error = function(e) {
    # This will be caught and an error message generated
    NULL
  })
  
  if (is.null(data)) {
    output_list$error_message <- "Error reading the DATA from the file. Ensure data starts on row 4."
    return(output_list)
  }
  
  missing_columns_header <- setdiff(c("name", "description"), colnames(params))
  if (length(missing_columns_header) > 0) {
    output_list$error_message <- paste0("Missing HEADER columns: ", paste(missing_columns_header, collapse = ", "))
    return(output_list)
  }
  
  missing_columns_data <- setdiff(c("chr", "start", "end", "region_id", "strand", "color"), colnames(data))
  if (length(missing_columns_data) > 0) {
    output_list$error_message <- paste0("Missing DATA columns: ", paste(missing_columns_data, collapse = ", "))
    return(output_list)
  }
  
  params_list <- as.list(params)
  header <- str_glue('track name=\"{params_list$name}\" description=\"{params_list$description}\" visibility=3 itemRgb=On')
  
  data_processed <- data %>% 
    mutate(color = str_replace_all(color,"_",",")) %>% 
    mutate(score = 1000, thickStart = start, thickEnd = end) %>% 
    relocate(score, .after = region_id) %>% 
    relocate(all_of(c("thickStart", "thickEnd")), .before= color)
  
  output_list$header <- header
  output_list$data <- data_processed
  output_list$status <- TRUE
  output_list$error_message <- "OK"
  output_list
}

# ---- UI for App1 ----
app1_ui <- fluidPage(
  titlePanel("Convert XLSX to BED Files"),
  sidebarLayout(
    sidebarPanel(
      h4("Instructions"),
      helpText(
        "Upload one or more XLSX files for conversion. Each file must contain the following columns:",
        tags$ul(
          tags$li(strong("'chr'")),
          tags$li(strong("'start'")),
          tags$li(strong("'end'"))
        ),
        "Additional columns will be ignored."
      ),
      fileInput("app1_fileInput", "Select XLSX Files",
                multiple = TRUE, accept = c(".xlsx")),
      actionButton("app1_process", "Convert to BED", class = "btn btn-primary"),
      tags$hr(),
      h5("Example File"),
      tags$a(href = "example_xlsx2bed.xlsx", "Download Example XLSX File (Simple BED)", target = "_blank"),
      width = 4
    ),
    mainPanel(
      h3("Download Results"),
      uiOutput("app1_downloadUI"),
      tags$hr(),
      h5("Note: Ensure your files follow the specified format for successful conversion."),
      width = 8
    )
  )
)

# ---- UI for App2 ----
app2_ui <- fluidPage(
  titlePanel("Convert XLSX files to BED track files (UCSC Genome Browser format)"),
  sidebarLayout(
    sidebarPanel(
      h4("Instructions"),
      helpText(
        "Upload one or more XLSX files for conversion. Each file must follow format and instructions for UCSC browser upload in example.xlsx file:"
      ),
      fileInput("app2_fileInput", "Select XLSX Files",
                multiple = TRUE, accept = c(".xlsx")),
      actionButton("app2_process", "Convert to UCSC track", class = "btn btn-primary"),
      tags$hr(),
      h5("Example File"),
      tags$a(href = "example_xlsx2track.xlsx", "Download Example XLSX File (UCSC BED Track)", target = "_blank"),
      width = 4
    ),
    mainPanel(
      h3("Download Results"),
      uiOutput("app2_downloadUI"),
      tags$hr(),
      h5("Note: Ensure your files follow the specified format for successful conversion."),
      width = 8
    )
  )
)

# ---- Main UI ----
ui <- navbarPage(
  theme = shinythemes::shinytheme("cerulean"),
  title = "BED Conversion Tools",
  
  tags$head(
    tags$style(HTML("
            /* Target the specific download links by their IDs */
            #app1_downloadData, #app2_downloadData {
                color: red; /* Set the text color to red */
                font-weight: bold; /* Make the font bold to stand out */
            }
        "))
  ),
  
  tabPanel("XLSX to BED", app1_ui),
  tabPanel("XLSX to UCSC track", app2_ui)
)

# ---- Server for app1 ----
app1_server <- function(input, output, session) {
  processed_files <- reactiveVal(NULL)
  
  observeEvent(input$app1_process, {
    req(input$app1_fileInput)
    
    # Reset download link on new processing
    processed_files(NULL) 
    
    files <- input$app1_fileInput$datapath
    filenames <- input$app1_fileInput$name
    output_dir <- tempdir()
    
    # Use map2 to iterate and collect results
    results <- map2(files, filenames, function(file, filename) {
      tryCatch({
        data <- read_excel(file)
        required_cols <- c("chr", "start", "end")
        if (!all(required_cols %in% colnames(data))) {
          stop(paste("File is missing required BED columns (chr, start, end)."))
        }
        bed_data <- data %>%
          select(all_of(required_cols)) %>%
          arrange(chr, start, end)
        
        bed_file <- file.path(output_dir, paste0(tools::file_path_sans_ext(filename), ".bed"))
        write_delim(bed_data, bed_file, delim = "\t", col_names = FALSE)
        
        # Return the path on success
        return(list(status = "success", path = bed_file, filename = filename))
      }, error = function(e) {
        # Return the error on failure
        return(list(status = "error", message = e$message, filename = filename))
      })
    })
    
    # Separate successful and failed files
    successes <- keep(results, ~ .$status == "success")
    errors <- keep(results, ~ .$status == "error")
    
    # Show notifications for errors
    walk(errors, ~ showNotification(paste("Error processing", .$filename, ":", .$message), type = "error", duration = 10))
    
    # If there are any successful files, zip them
    if (length(successes) > 0) {
      output_files <- map_chr(successes, "path")
      zip_file <- file.path(output_dir, "converted_simple_bed_files.zip")
      zip(zip_file, output_files, flags = '-jr9X')
      processed_files(zip_file)
    }
  })
  
  output$app1_downloadUI <- renderUI({
    req(processed_files())
    downloadLink("app1_downloadData", "Download ZIP Archive of BED Files")
  })
  
  output$app1_downloadData <- downloadHandler(
    filename = function() {
      "converted_simple_bed_files.zip"
    },
    content = function(file) {
      file.copy(processed_files(), file)
    }
  )
}

# ---- Server for app2 ----
app2_server <- function(input, output, session) {
  processed_files <- reactiveVal(NULL)
  
  observeEvent(input$app2_process, {
    req(input$app2_fileInput)
    
    # Reset download link on new processing
    processed_files(NULL)
    
    files <- input$app2_fileInput$datapath
    filenames <- input$app2_fileInput$name
    output_dir <- tempdir()
    
    results <- map2(files, filenames, function(file, filename) {
      tryCatch({
        validated_data <- file_validation(file)
        if (isFALSE(validated_data$status)) {
          stop(validated_data$error_message)
        }
        bed_file <- file.path(output_dir, paste0(tools::file_path_sans_ext(filename), ".bed"))
        write_lines(validated_data$header, bed_file)
        write_tsv(validated_data$data, file = bed_file, append = TRUE, col_names = FALSE)
        
        return(list(status = "success", path = bed_file, filename = filename))
      }, error = function(e) {
        return(list(status = "error", message = e$message, filename = filename))
      })
    })
    
    successes <- keep(results, ~ .$status == "success")
    errors <- keep(results, ~ .$status == "error")
    
    walk(errors, ~ showNotification(paste("Error processing", .$filename, ":", .$message), type = "error", duration = 10))
    
    if (length(successes) > 0) {
      output_files <- map_chr(successes, "path")
      zip_file <- file.path(output_dir, "converted_ucsc_track_files.zip")
      zip(zip_file, output_files, flags = '-jr9X')
      processed_files(zip_file)
    }
  })
  
  output$app2_downloadUI <- renderUI({
    req(processed_files())
    downloadLink("app2_downloadData", "Download ZIP Archive of BED Files")
  })
  
  output$app2_downloadData <- downloadHandler(
    filename = function() {
      "converted_ucsc_track_files.zip"
    },
    content = function(file) {
      file.copy(processed_files(), file)
    }
  )
}

# ---- Main Server ----
server <- function(input, output, session) {
  # Call the server logic for each app directly
  app1_server(input, output, session)
  app2_server(input, output, session)
}

# ---- Run the app ----
shinyApp(ui = ui, server = server)