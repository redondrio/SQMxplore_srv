# Load packages
library(shiny)
library(shinyFiles)
library(SQMtools)
library(ggplot2)
library(DT)
library(data.table)

ui <- navbarPage("SQMxplore",
  # Page Input ----
  tabPanel("Dataset",
    fluidPage("",
      mainPanel(
        # If option added, include in switch at load section
        selectInput("type_load", "Select type of input",
          choices = c("Load SQM project",
            "Load from pre-saved RDS file")),
        h4("Select path"),
        shinyDirButton("samples_path", "Input directory", "Select"),
        h5(),
        verbatimTextOutput("out_samples_path", placeholder = TRUE),
        h1(),
        selectInput("project", "Select project",
        choices = c()), #updated
        h1(),
        actionButton("proj_load", "Load project"),
        h5("Currently loaded project:"),
        verbatimTextOutput("out_project", placeholder = TRUE)
      ) # Close panel
    ) # Close layout
  ), # Close Input page

  # Page Summary ----
  tabPanel("Summary",
    # Sidebar layout with input and output definitions
    fluidPage("", h1("Summary page"),
      # Display Table Output
      conditionalPanel(
        condition = "input.type_load=='Load SQM project'",
        h3("Reads summary"),
        DT::dataTableOutput("reads_sum"),
        fluidRow(
          column(6, plotOutput(outputId = "reads_readbar")),
          column(6, plotOutput(outputId = "reads_basebar"))
        ),

        h3("Contigs summary"),
        DT::dataTableOutput("contigs_sum"),

        h3("Taxa summary"),
        DT::dataTableOutput("taxa_sum"),

        h3("Orfs summary"),
        DT::dataTableOutput("orfs_sum"),
        fluidRow(
          column(6, selectInput(
            "orfs_row1", "Choose a value to display",
            choices = NA)),
          column(6, selectInput(
            "orfs_row2", "Choose a value to display",
            choices = NA))
        ),
        fluidRow(
          column(6, plotOutput(outputId = "orfs_bar1")),
          column(6, plotOutput(outputId = "orfs_bar2"))
        ),

        h3("Bins summary"),
        DT::dataTableOutput("bins_sum")
      ), # Close conditional panel

      conditionalPanel(
        condition = "input.type_load!='Load SQMlite from minimum tables'",
        h5("This object does not include summary tables")
      ) # Close conditional panel
    ) # Close layout
  ), # Close page

  # Page Tables ----
  tabPanel("Tables",
    # Sidebar layout with input and output definitions
    fluidPage("DataTable",
      # Row panel for inputs
      fluidRow(
        # Input:
        column(3,
          selectInput("lev1_tab",
            "Choose a section to display",
            choices = "", selected = ""), #updated
          selectInput("lev2_tab",
            "Choose a subsection to display",
            choices = ""), #updated
          selectInput("lev3_tab",
            "Choose a unit to display",
            choices = ""), #updated
          downloadButton("tabDown", "Download Table")
        ),
        column(8,
          selectizeInput("cols_tab",
            "Choose samples to display",
            choices = NA, multiple = TRUE) #updated
        )
      ), # Close row panel

      # Display Table Output
      DT::dataTableOutput("table")
    ) # Close layout
  ), # Close page

  # Page Taxonomy ----
  tabPanel("Plot Taxonomy",
    sidebarLayout(
      # Sidebar panel for inputs ----
      sidebarPanel(
        h4("Data"),
        # Input: Menu for the rank, units and samples
        selectInput("rank_tax", "Choose a rank to display",
          choices = "", selected = ""), #updated
        selectInput("count_tax", "Choose a unit to display",
          choices = "", selected = ""), #updated
        selectizeInput("samples_tax", "Select samples",
          choices = NULL,
          selected = NULL,
          multiple = TRUE), #updated

        h4("Taxa"),
        conditionalPanel(
          condition = "!input.sel_tax",
          # Input: Number for the number of taxa
          numericInput("n_tax", "Choose the number of taxa",
            value = 1, min = 0) #updated
        ),
        conditionalPanel(
          condition = "input.sel_tax",
          selectizeInput("tax_tax", "Selected taxa",
            choices = NULL,
            multiple = TRUE) #updated
        ),
        # Input: Write taxa names and override N
        checkboxInput("sel_tax", "Manually pick plotted taxa",
          value = FALSE),

        h4("Options"),
        checkboxInput("others_tax", "Show other reads",
          value = TRUE),
        conditionalPanel(
          condition = "!input.sel_tax",
          checkboxInput("unmapped_tax",
            "Ignore unmapped reads",
            value = FALSE),
          checkboxInput("unclass_tax",
            "Ignore unclassified reads",
            value = FALSE),
          checkboxInput("partial_tax",
            "Ignore partial classifications",
            value = FALSE),
        ),
        checkboxInput("rescale_tax", "Rescale to 100%",
                value = FALSE),
        numericInput("base_size_tax", "Font size", value = 11)
      ), # Close sidebar panel

      # Main panel for displaying outputs ----
      mainPanel(
        plotOutput(outputId = "taxPlot"),
          downloadButton("taxPlotDown", "Download Plot"),
          selectizeInput("dev_taxPlot", "Select device",
            choices = c("pdf", "jpeg", "png"),
            selected = "pdf"),
          numericInput("width_taxPlot", "Plot width",
            value = 20, min = 0),
          numericInput("height_taxPlot", "Plot height",
            value = 20, min = 0),
          selectizeInput("unit_taxPlot", "Units",
            choices = c("in", "cm", "mm", "px"),
            selected = "cm")
      ) # Close main panel
    ) # Close layout
  ), # Close Taxonomy page

  # Page Functions ----
  tabPanel("Plot Functions",
    sidebarLayout(
      # Sidebar panel for inputs ----
      sidebarPanel(
        h4("Data"),
        # Input: Menu for the rank, units and samples
        selectInput("fun_level_fun",
          "Choose an annotation to display",
          choices = "", selected = ""), #updated
        selectInput("count_fun",
          "Choose a unit to display",
          choices = "", selected = ""), #updated
        selectizeInput("samples_fun",
          "Select samples",
          choices = NULL,
          selected = NULL,
          multiple = TRUE), #updated

        h4("Functions"),
        conditionalPanel(
            condition = "!input.sel_fun && !input.load_fun",
          # Input: Number for the number of taxa
          numericInput("n_fun", "Choose the number of functions",
            value = 1, min = 0) #updated
        ),

          # Input: Write function names and override N
          checkboxInput("sel_fun", "Manually pick plotted functions",
            value = FALSE),
          conditionalPanel(
            condition = "input.sel_fun",
            selectizeInput("fun_fun", "Selected functions",
              choices = NULL,
              selected  = NULL,
              multiple = TRUE) #updated
          ),

          # Input: Load functions from ref file and override N
          checkboxInput("load_fun", "Load reference functions from file",
            value = FALSE),
          conditionalPanel(
            condition = "input.load_fun",
            selectizeInput("ref_fun", "Selected functions",
              choices = NULL,
              selected  = NULL,
              multiple = TRUE) #updated
          ),

        h4("Options"),
        conditionalPanel(
          condition = "!input.sel_fun",
          checkboxInput("unmapped_fun",
            "Ignore unmapped reads", value = FALSE),
          checkboxInput("unclass_fun",
            "Ignore unclassified reads", value = FALSE)
        ),
        numericInput("base_size_fun", "Change Font size", 11)
      ), # Close sidebar panel

      # Main panel for displaying outputs ----
      mainPanel(
        # Output: Histogram
          plotOutput(outputId = "funPlot"),
          downloadButton("funPlotDown", "Download Plot"),
          selectizeInput("dev_funPlot", "Select device",
            choices = c("pdf", "jpeg", "png"),
            selected = "pdf"),
          numericInput("width_funPlot", "Plot width",
            value = 20, min = 0),
          numericInput("height_funPlot", "Plot height",
            value = 20, min = 0),
          selectizeInput("unit_funPlot", "Units",
            choices = c("in", "cm", "mm", "px"),
            selected = "cm")
        ) # Close main panel
      ) # Close layout
  ) # Close Functions page
) # Close UI