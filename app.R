###############################################################################
# app.R  –  RNAdegDB Explorer  v2
# - Graceful empty-result handling
# - Centre-aligned inputs → tabbed Results page
# - Interactive (searchable, sortable, downloadable) table via DT
# - Readable column names with units
# - Light aesthetic polish via {bslib}
###############################################################################

library(shiny)
library(httr)
library(jsonlite)
library(DT)     # interactive tables
library(bslib)  # Bootstrap 5 themes
library(magrittr)

## ---------------------------------------------------------------------------
## 1. Supabase credentials ----------------------------------------------------
## ---------------------------------------------------------------------------
SUPABASE_URL <- Sys.getenv("SUPABASE_URL")
SUPABASE_KEY <- Sys.getenv("SUPABASE_KEY")


auth_headers <- add_headers(
  apiKey       = SUPABASE_KEY,
  Authorization = paste("Bearer", SUPABASE_KEY)
)


# helper: pull ALL DISTINCT values of <column> from <table> -------------------
distinct_vals <- function(table, column, chunk = 1000) {
  out    <- character()    # will collect values here
  offset <- 0L
  
  repeat {
    # ask for a slice: limit <chunk> rows, starting at <offset>
    url <- sprintf(
      "%s/rest/v1/%s?select=%s&limit=%d&offset=%d",
      SUPABASE_URL, table, column, chunk, offset
    )
    
    txt <- content(GET(url, auth_headers), as = "text", encoding = "UTF-8")
    dat <- fromJSON(txt, simplifyDataFrame = TRUE)
    vals <- dat[[column]]
    
    # add & move to next slice
    out    <- c(out, vals)
    offset <- offset + chunk
    
    # stop once we got fewer rows than the slice size
    if (length(vals) < chunk) break
  }
  
  unique(out[!is.na(out)])
}


# pull choices once at app startup (fast & keeps UI code clean) ----------------
sample_choices   <- distinct_vals("sampledetails",  "sample")
cellline_choices <- distinct_vals("sampledetails",  "cell_line")
dataset_choices  <- distinct_vals("sampledetails",  "dataset")
gene_choices <- distinct_vals("featuredetails", "feature_id")

## ---------------------------------------------------------------------------
## 3.  THEME  -----------------------------------------------------------------
## ---------------------------------------------------------------------------
theme <- bs_theme(
  version   = 5,
  bootswatch = "flatly",
  base_font = font_google("Roboto")
)

## ---------------------------------------------------------------------------
## 4.  UI  --------------------------------------------------------------------
## ---------------------------------------------------------------------------
ui <- fluidPage(
  theme = theme,
  
  # -- app header -------------------------------------------------------------
  titlePanel(
    title = "RNAdecayCafe Explorer",
    icon("dna", class = "me-2") |> tagAppendAttributes(style = "color:#0d6efd;")
  ),
  
  # -- tabset -----------------------------------------------------------------
  tabsetPanel(id = "tabs",
              
              ## ---- 4.1  Query tab -----------------------------------------------------
              tabPanel(
                title = "Query",
                
                sidebarLayout(
                  
                  # left-hand filters ----------------------------------------------------
                  sidebarPanel(width = 4,
                               
                               selectizeInput(
                                 "genes",
                                 "Gene name(s)",
                                 choices = NULL,
                                 multiple = TRUE, 
                                 options = list(placeholder = "start typing...")
                               ),
                               
                               selectizeInput(
                                 "cell_lines", "Cell line",
                                 choices   = cellline_choices,
                                 multiple  = TRUE,
                                 options   = list(placeholder = "all")
                               ),
                               
                               selectizeInput(
                                 "datasets", "Dataset",
                                 choices   = dataset_choices,
                                 multiple  = TRUE,
                                 options   = list(placeholder = "all")
                               ),
                               
                               selectizeInput(
                                 "samples", "Sample",
                                 choices   = sample_choices,
                                 multiple  = TRUE,
                                 options   = list(placeholder = "all")
                               ),
                               
                               div(
                                 class = "text-center",
                                 actionButton("go", "Fetch data", class = "btn-primary btn-lg mt-2")
                               )
                  ),
                  
                  # right-hand documentation + messages ---------------------------------
                  mainPanel(width = 8,
                            
                            h4("How to use the filters"),
                            
                            tags$p(
                              "When you press ", tags$strong("Fetch data"), ", the RNA half-life ",
                              "data you requested will appear in the ",
                              tags$strong("Results"), " and ", tags$strong("Averages"), " tabs. ",
                              tags$strong("Results"), " contains every matching sample-level estimate; ",
                              tags$strong("Averages"), " summarises those values by cell line."
                            ),
                            
                            HTML("
            <ul>
              <li><b>Gene name(s)</b>: enter one or many gene symbols to retrieve their RNA-stability metrics.</li>
              <li><b>Sample</b>: choose specific SRA accessions / sample names. Leaving this set to <i>all</i> shows every sample in the Results tab. (It does not affect the Averages tab.)</li>
              <li><b>Cell line</b>: pick the cell line(s) you care about. The filter applies to both Results and Averages.</li>
              <li><b>Dataset</b>: restrict Results to one or more studies / citations. (Again, Averages always aggregate over all datasets.)</li>
            </ul>
          "),
                            
                            uiOutput("query_msg")  # “no rows matched…” alerts appear here
                  )
                )
              ),
              
              ## ---- 4.2  Results tab ---------------------------------------------------
              tabPanel(
                title = "Results",
                br(),
                DTOutput("result_tbl")
              ),
              
              ## ---- 4.3  Averages tab --------------------------------------------------
              tabPanel(
                title = "Averages",
                br(),
                DTOutput("avg_tbl")
              )
  )
)


## ---------------------------------------------------------------------------
## 5.  SERVER -----------------------------------------------------------------
## ---------------------------------------------------------------------------
server <- function(input, output, session) {
  
  updateSelectizeInput(
    session, "genes",
    choices = gene_choices,
    server = TRUE
  )
  
  
  build_filter <- function(col, vec, op = "in.") {
    if (length(vec) == 0) return(NULL)
    encoded <- URLencode(paste(vec, collapse = ","), reserved = TRUE)
    sprintf("%s=%s(%s)", col, op, encoded)
  }

  observeEvent(input$go, {
    
    ## --- 5.1 Build REST URL --------------------------------------------------
    filters <- c(
      build_filter("feature_id", input$genes),
      build_filter("sample",                      input$samples),
      build_filter("sampledetails.cell_line",     input$cell_lines),
      build_filter("sampledetails.dataset",       input$datasets)
    )
    where <- paste(filters, collapse = "&")
    
    base <- paste0(
      "rateconstants?",
      "select=",
      "sample,feature_id,kdeg,halflife,donorm_kdeg,donorm_halflife,reads,",
      "sampledetails!inner(cell_line,dataset,total_reads)"
    )
    
    
    url <- sprintf("%s/rest/v1/%s&%s", SUPABASE_URL, base, where)
    
    ## --- 5.2 Fetch -----------------------------------------------------------
    resp_txt <- content(GET(url, auth_headers), as = "text", encoding = "UTF-8")
    dat <- fromJSON(resp_txt, flatten = TRUE)
    
    ## --- 5.3 User feedback ---------------------------------------------------
    if (!is.data.frame(dat) || nrow(dat) == 0) {
      output$query_msg <- renderUI(
        div(class = "alert alert-warning text-center",
            "No rows matched your filter — check that the sample belongs to the dataset, ",
            "cell-line spelling, etc.")
      )
      updateTabsetPanel(session, "tabs", "Query")
      return()
    } else {
      output$query_msg <- renderUI(NULL)            # clear old messages
    }
    
    ## --- 5.4  Clean column names + order ------------------------------------
    names(dat) <- c(
      "sample", "gene_name",
      "kdeg (hr\u207B\u00B9)",          # hr⁻¹
      "halflife (hr)",
      "donorm kdeg (hr\u207B\u00B9)",
      "donorm halflife (hr)",
      "reads",
      "dataset",
      "cell_line", "total reads"
    )
    dat <- dat[, c("sample","gene_name","cell_line","dataset",
                   "kdeg (hr\u207B\u00B9)","donorm kdeg (hr\u207B\u00B9)",
                   "halflife (hr)","donorm halflife (hr)","reads")]
    
    ## --- 5.5  Render interactive table --------------------------------------
    output$result_tbl <- renderDT({
      num_cols <- names(dat)[sapply(dat, is.numeric)]
      datatable(
        dat,
        extensions = c("Buttons","FixedHeader"),
        options = list(
          pageLength   = 25,
          dom          = "Bfrtip",
          buttons      = c("copy", "csv", "excel"),
          fixedHeader  = TRUE,
          searchHighlight = TRUE
        ),
        class = "stripe hover compact",
        rownames = FALSE
      ) %>%
        formatRound(columns = num_cols, digits = 2)
    }) 
    
    ## -----------------------  NEW averages section  ----------------------------
    avg_filters <- c(
      build_filter("feature_id", input$genes),
      build_filter("cell_line", input$cell_lines)      # ← adds in.(...) if any
    )
    avg_where <- paste(avg_filters, collapse = "&")    # could be "" (empty)
    
    avg_url <- sprintf(
      "%s/rest/v1/averages?select=%s%s",
      SUPABASE_URL,
      paste(
        "feature_id,cell_line,avg_kdeg,avg_donorm_kdeg,avg_halflife,avg_donorm_halflife"
      ),
      if (nzchar(avg_where)) paste0("&", avg_where) else ""   # ← append only if present
    )
    
    avg_txt <- content(GET(avg_url, auth_headers), as = "text")
    avg_dat <- fromJSON(avg_txt, flatten = TRUE)
    
    if (is.data.frame(avg_dat) && nrow(avg_dat)) {
      names(avg_dat) <- c(
        "gene_name", "cell_line",
        "avg kdeg (hr\u207B\u00B9)", "avg donorm kdeg (hr\u207B\u00B9)",
        "avg halflife (hr)",        "avg donorm halflife (hr)"
      )
      output$avg_tbl <- renderDT({
        num_cols <- names(avg_dat)[sapply(avg_dat, is.numeric)]
        datatable(
          avg_dat,
          class = "stripe hover compact",
          rownames = FALSE,
          extensions = c("Buttons","FixedHeader"),
          options = list(
            dom            = "Bfrtip",
            buttons        = c("copy", "csv", "excel"),
            fixedHeader    = TRUE,
            pageLength     = 25,
            searchHighlight = TRUE
          )
        ) %>%
          formatRound(columns = num_cols, digits = 2)
      }) 
    } else {
      output$avg_tbl <- renderDT(NULL)   # or show “no rows” message
    }
    
    
    ## --- 5.7  Switch to Results tab -----------------------------------------
    updateTabsetPanel(session, "tabs", "Results")
  })
}

## ---------------------------------------------------------------------------
## 6.  RUN --------------------------------------------------------------------
## ---------------------------------------------------------------------------
shinyApp(ui, server)

