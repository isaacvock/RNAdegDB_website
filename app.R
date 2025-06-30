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
    title = "RNAdecayCafe Explorer"
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
                            
                            uiOutput("query_msg"),  # “no rows matched…” alerts appear here
                            
                            ### About block ####
                            hr(), # Thin horizontal rule
                            h4("About RNAdecayCafe"),
                            
                            tags$p(
                              "RNAdecayCafe is a curated set of high-quality, highly reproducible RNA ",
                              "half-life estimates derived from nucleotide-recoding RNA-seq (NR-seq) ",
                              "experiments—namely ", 
                              tags$a(
                                href = 'https://www.nature.com/articles/nmeth.4435',
                                "SLAM-seq"
                                )
                              , " and ", 
                              tags$a(
                                href = 'https://www.nature.com/articles/nmeth.4582',
                                "TimeLapse-seq"
                                ),
                              ". ",
                              "NR-seq currently represents the state-of-the-art in high-throughput RNA-",
                              "stability quantification."
                            ),
                            
                            tags$p(
                              "We re-processed all published NR-seq datasets using the ",
                              tags$a(href = 'https://www.biorxiv.org/content/10.1101/2024.10.14.617411v1',
                                     "EZbakR suite"), " (pipeline ",
                              tags$a(href = 'https://github.com/isaacvock/fastq2EZbakR', "fastq2EZbakR"),
                              ", analysis package ",
                              tags$a(href = 'https://github.com/isaacvock/EZbakR', "EZbakR"),
                              ") and found that pulse-label NR-seq data provide highly consistent half-life ",
                              "estimates across labs, protocols, and cell lines."
                            ),
                            
                            tags$p(
                              "The current release of RNAdecayCafe aggregates ",
                              tags$strong("66 samples"), ", ",
                              tags$strong("16 datasets"), ", ",
                              tags$strong("12 cell lines"), ", and ",
                              tags$strong("17 573 genes"), ". ",
                              "Below is a heatmap of the sample-to-sample log(half-life) Pearson correlations.",
                              " The main distinguishing factor between samples is the time for which cells were labeled,",
                              " which sets the dynamic range of the half-life estimates."
                            ),
                            
                            ### heatmap image
                            tags$img(
                              src = "heatmap.png",
                              alt = "Sample-to-sample correlation heatmap",
                              style = "max-width:80%; height:auto; margin-top:1rem;"
                            )
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
              ),
              
              ## ---- 4.4 Sample Details ---------------------------------------
              tabPanel(
                title = "Sample Details",
                br(),
                DTOutput("sample_tbl")
              ),
              
              ## ---- 4.5 Explaining tables ------------------------------------
              tabPanel(
                title = "Explaining the Tables",
                br(),
                
                h4("Results"),
                tags$ol(
                  tags$li(tags$b("sample"), ": SRA accession for the data that the estimate comes from."),
                  tags$li(tags$b("gene_name"), ": Gene symbol (hg38 annotation)."),
                  tags$li(tags$b("cell_line"), ": Cell line in which the RNA was measured."),
                  tags$li(tags$b("dataset"), ": Citation for the source data."),
                  tags$li(tags$b("kdeg"), ": Degradation-rate constant estimate (units = 1 / hr)."),
                  tags$li(
                    tags$b("donorm kdeg"),
                    ": Dropout-normalised ", tags$code("kdeg"),
                    " estimate. Normalizes out global differences in estimates seen between samples, often caused by dropout (see ",
                    tags$a(
                      href = "https://www.biorxiv.org/content/10.1101/2023.05.24.542133v1",
                      "our dropout preprint"
                    ),
                    " and ",
                    tags$a(
                      href="https://academic.oup.com/nar/article/52/7/e35/7612100",
                      "the Erhard lab's paper"
                    ),
                    " for details). We suggest using the dropout normalized values, especially when comparing across cell lines or samples,",
                    " but we provide the raw estimates for transparency. NOTE: dropout normalization does not impact",
                    " relative ordering of half-lives; thus, if RNA A is more stable than RNA B before normalization, it will still have a longer",
                    " halflife after normalization."
                  ),
                  tags$li(tags$b("halflife"), ": log(2) / kdeg; average lifetime of the RNA."),
                  tags$li(tags$b("donorm halflife"), ": Dropout–normalised halflife."),
                  tags$li(tags$b("reads"), ": Number of reads contributing to the estimate.")
                ),
                
                # ---------------- Averages table -----------------------------------------
                h4("Averages"),
                tags$ol(
                  tags$li(tags$b("gene_name"), ": as above."),
                  tags$li(tags$b("cell_line"), ": as above."),
                  tags$li(
                    tags$b("avg kdeg"),
                    ": Uncertainty-weighted mean kdeg across samples of that cell line."
                  ),
                  tags$li(
                    tags$b("avg donorm kdeg"),
                    ": Dropout–normalised uncertainty weighted mean kdeg."
                  ),
                  tags$li(
                    tags$b("avg halflife"), " and ", tags$b("avg donorm halflife"),
                    ": log(2) / the relevant avg kdeg."
                  )
                ),
                
                # ---------------- Sample-Details table -----------------------------------
                h4("Sample Details"),
                tags$ol(
                  tags$li(tags$b("sample"), ": SRA accession code."),
                  tags$li(tags$b("dataset"), ": Citation for the dataset."),
                  tags$li(
                    tags$b("pnew"),
                    ": T-to-C conversion rate in *new* (labelled) reads. Estimated by EZbakR."
                  ),
                  tags$li(
                    tags$b("pold"),
                    ": T-to-C conversion rate in *old* (unlabelled) reads. Estimated by EZbakR."
                  ),
                  tags$li(
                    tags$b("label time"),
                    ": Duration of s<sup>4</sup>U labelling; determines dynamic range of half-life estimates. Half-lives much shorter or longer than this are difficult to accurately estimate."
                  ),
                  tags$li(tags$b("cell line"), ": Cell line used."),
                  tags$li(tags$b("threePseq"), ": Indicates that the data comes from 3′-end sequencing (vs total RNA)."),
                  tags$li(
                    tags$b("total reads"),
                    ": Total number of reads after filtering intronic and multi-mapping reads."
                  ),
                  tags$li(
                    tags$b("median halflife"),
                    ": Median halflife for the sample. High values are suggestive of dropout; ",
                    "deep sequencing can lower medians by capturing many rapidly degraded, ",
                    "low-abundance transcripts (e.g., targets of NMD)."
                  )
                )
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
  
  
  ## ---- pull the full sampledetails table ------------------------------------
  sample_url <- sprintf("%s/rest/v1/sampledetails?select=*", SUPABASE_URL)
  sample_txt <- content(GET(sample_url, auth_headers), as = "text", encoding = "UTF-8")
  sample_dat <- fromJSON(sample_txt, flatten = TRUE)
  
  # optional: prettify column names
  names(sample_dat) <- gsub("_", " ", names(sample_dat), fixed = TRUE)
  
  output$sample_tbl <- renderDT({
    num_cols <- names(sample_dat)[sapply(sample_dat, is.numeric)]
    datatable(
      sample_dat,
      extensions = c("Buttons", "FixedHeader"),
      options = list(
        dom            = "Bfrtip",
        buttons        = c("copy", "csv", "excel"),
        fixedHeader    = TRUE,
        pageLength     = 25,
        searchHighlight = TRUE
      ),
      class = "stripe hover compact",
      rownames = FALSE
    ) %>%
      ## ---- precision tweaks --------------------------------------------------
      # pnew & pold  → 4 decimals
      formatRound(columns = c("pnew", "pold"), digits = 4) %>%
      
      formatRound(columns = c("median halflife"), digits = 2) %>%
      
      # label time   → 2 decimals       (name now has a space)
      formatRound(columns = "label time", digits = 2) %>%
      
      # total reads  → no decimals
      formatRound(columns = "total reads", digits = 0)
  })
  
  

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

