# Load packages
library(shiny)
library(shinyFiles)
library(SQMtools)
library(ggplot2)
library(DT)
library(data.table)

# Server ----
server <- function(input, output, clientData, session) {
  # Initialise reactive values to store data
  reactiveData <- reactiveValues()

  # Set the paths ----
  roots <- c(root = "/srv/shiny-server/shiny_input")
  shinyDirChoose(input, "samples_path", roots = roots,
    filetypes = c("", "txt", "bigWig", "tsv", "csv", "bw"))

  # Update loading options ----
  observeEvent(input$samples_path, {
    updateSelectInput(session, "project", "Select project",
      choices = (list.files(parseDirPath(roots, input$samples_path))))
  })
  output$out_samples_path <- renderPrint(
    {parseDirPath(roots, input$samples_path)})

  # Tables dir to feed to loadSQMlite
  tab_dir <- reactive({
    paste0(parseDirPath(roots, input$samples_path), "/",
      input$project, "/results/tables/")
  })
  # Results dir to load parsed stat files
  res_dir <- reactive({
    paste0(parseDirPath(roots, input$samples_path), "/",
      input$project, "/results/")
  })

  # Load the SQM object and stats files ----
  observeEvent(input$proj_load, {
    tryCatch({
      showModal(modalDialog(title = "Loading", easyclose = TRUE))
      reactiveData$SQM <- switch(input$type_load,
        "Load SQM project" = {
          loadSQMlite(tab_dir())
        },
        "Load from pre-saved RDS file" = {
          readRDS(paste0(parseDirPath(roots, input$samples_path), "/",
            input$project))
        })
      # Load the stats file
      # add check.names=FALSE to prevent R from changing column names
      if (input$type_load == "Load SQM project") {
        reactiveData$reads_st <- read.csv(paste0(res_dir(),
          paste0("22.reads.tsv")), row.names = 1, header = TRUE, sep = "\t")
        reactiveData$contigs_st <- read.csv(paste0(res_dir(),
          paste0("22.contigs.tsv")), row.names = 1, header = TRUE, sep = "\t")
        reactiveData$taxa_st <- read.csv(paste0(res_dir(),
          paste0("22.taxa.tsv")), row.names = 1, header = TRUE, sep = "\t")
        reactiveData$orfs_st <- read.csv(paste0(res_dir(),
          paste0("22.orfs.tsv")), row.names = 1, header = TRUE, sep = "\t")
        reactiveData$orfs_st[is.na(reactiveData$orfs_st)] <- 0
        reactiveData$bins_st <- read.csv(paste0(res_dir(),
          paste0("22.bins.tsv")), row.names = 1, header = TRUE, sep = "\t")
      } else { # if no stat files, empty objects
        reactiveData$reads_st <- matrix()
        reactiveData$contigs_st <- matrix()
        reactiveData$taxa_st <- matrix()
        reactiveData$orfs_st <- matrix()
        reactiveData$bins_st <- matrix()
      }
      output$out_project <- renderText(isolate(input$project))
      showModal(modalDialog(title = "Loaded",
        "Your project is ready", easyclose = TRUE))
    }, warning = function(warn) {
      showModal(modalDialog(title = "Stopped",
        "There were warnings during loading", easyclose = TRUE))
    }, error = function(error) {
      showModal(modalDialog(title = "Loading error",
        "Please check project and stat files", easyclose = TRUE))
    }) # Close tryCatch
  }) # Close proj_load observer

  # Update Taxonomy Inputs ----
  observe({
    updateSelectInput(session, "rank_tax",
      choices = names(reactiveData$SQM[["taxa"]])
    )
  }) # Close observer

  observe({
    updateSelectInput(session, "count_tax",
      choices = names(reactiveData$SQM[["taxa"]][[input$rank_tax]])
    )
  }) # Close observer

  observe({
    updateSelectizeInput(session, "samples_tax",
      choices = reactiveData$SQM$misc$samples,
      selected = reactiveData$SQM$misc$samples
    )
  }) # Close observer

  observe({
    uniques <- unique(rownames(
      reactiveData$SQM[["taxa"]][[input$rank_tax]][[input$count_tax]]))
    if (length(uniques) > 10) {
      def_n_tax <- 10
    } else {
      def_n_tax <- length(uniques)
    }
    updateNumericInput(session, "n_tax",
      value = def_n_tax,
      max = length(uniques)
    )
  }) # Close observer

  observeEvent(c(input$proj_load, input$rank_tax), {
    updateSelectizeInput(session, "tax_tax",
      choices = unique(rownames(
        reactiveData$SQM[["taxa"]][[input$rank_tax]][[input$count_tax]]))
    )
  }) # Close observer

  # Output Taxonomy ----
  reactTaxPlot <- reactive({
    if (input$sel_tax) {
      # Plot by names of taxa
      plotTaxonomy(reactiveData$SQM,
                   rank = input$rank_tax,
                   count = input$count_tax,
                   tax = input$tax_tax,
                   others = input$others_tax,
                   samples = input$samples_tax,
                   ignore_unmapped = input$unmapped_tax,
                   ignore_unclassified = input$unclass_tax,
                   no_partial_classifications = input$partial_tax,
                   rescale = input$rescale_tax,
                   base_size = input$base_size_tax)
    } else {
      # Plot by number of taxa
      plotTaxonomy(reactiveData$SQM,
                   rank = input$rank_tax,
                   count = input$count_tax,
                   N = input$n_tax,
                   others = input$others_tax,
                   samples = input$samples_tax,
                   ignore_unmapped = input$unmapped_tax,
                   ignore_unclassified = input$unclass_tax,
                   no_partial_classifications = input$partial_tax,
                   rescale = input$rescale_tax,
                   base_size = input$base_size_tax)
    }
  })

  output$taxPlot <- renderPlot({
    print(reactTaxPlot())
  })

  output$taxPlotDown <- downloadHandler(
    filename = function() {
      paste("Tax_plot_", Sys.time(), ".pdf", sep = "")
    },
    content = function(file) {
      ggsave(file, reactTaxPlot(), device = "pdf")
    }
  )

  # Update Functions Inputs ----
  observe({
    updateSelectInput(session, "fun_level_fun",
      choices = names(reactiveData$SQM$functions),
    )
  })

  observeEvent(input$fun_level_fun, {
    # Add counts that exist in object and manually add percentage
    count_exist <- names(reactiveData$SQM[["functions"]][[input$fun_level_fun]])
    count_accepted <- c("abund", "percent", "bases", "tpm", "copy_number")
    updateSelectizeInput(session, "count_fun",
      choices = c("percent", intersect(count_exist, count_accepted))
    )
  })

  observe({
    updateSelectizeInput(session, "samples_fun",
      choices = reactiveData$SQM$misc$samples,
      selected = reactiveData$SQM$misc$samples
    )
  })

  observe({
    uniques <- unique(rownames(
      reactiveData$SQM[["functions"]][[input$fun_level_fun]][["abund"]]))
    if (length(uniques) > 10) {
      def_n_fun <- 10
    } else {
      def_n_fun <- length(uniques)
    }
    updateNumericInput(session, "n_fun",
      value = def_n_fun,
      max = length(uniques)
    )
  }) # Close observer

  observeEvent(c(input$proj_load, input$fun_level_fun), {
    updateSelectizeInput(session, "fun_fun",
      choices = rownames(
        reactiveData$SQM[["functions"]][[input$fun_level_fun]][["abund"]]),
      server = TRUE
    )
  }) # Close fun_level_fun observer

  # Output Functions ----
  reactFunPlot <- reactive({
    if (input$sel_fun) {
      plotFunctions(reactiveData$SQM,
                    fun_level = input$fun_level_fun,
                    count = input$count_fun,
                    fun = input$fun_fun,
                    samples = input$samples_fun,
                    ignore_unmapped = input$unmapped_fun,
                    ignore_unclassified = input$unclass_fun,
                    base_size = input$base_size_fun)
    } else {
      plotFunctions(reactiveData$SQM,
                    fun_level = input$fun_level_fun,
                    count = input$count_fun,
                    N = input$n_fun,
                    samples = input$samples_fun,
                    ignore_unmapped = input$unmapped_fun,
                    ignore_unclassified = input$unclass_fun,
                    base_size = input$base_size_fun)
    }
  })

  output$funPlot <- renderPlot({
    print(reactFunPlot())
  })

  output$funPlotDown <- downloadHandler(
    filename = function() {
      paste("Fun_plot_", Sys.time(), ".pdf", sep = "")
    },
    content = function(file) {
      ggsave(file, reactFunPlot(), device = "pdf")
    }
  )

  # Update Table Inputs ----
  observe({
    updateSelectInput(session, "lev1_tab",
      choices = switch(class(reactiveData$SQM),
        "SQM"     = c("orfs", "contigs", "bins", "taxa", "functions"),
        "SQMlite" = c("taxa", "functions"),
        "No project loaded"
      ),
      selected = switch(class(reactiveData$SQM),
        "SQM"     = "orfs",
        "SQMlite" = "taxa",
        "No project loaded"
      )
    )
  }) # Close observer

  observeEvent(input$lev1_tab, {
    updateSelectInput(session, "lev2_tab",
      choices = names(reactiveData$SQM[[input$lev1_tab]]),
      selected = names(reactiveData$SQM[[input$lev1_tab]])[1]
    )
  }) # Close lev1 observer

  observeEvent(c(input$lev1_tab, input$lev2_tab), {
    updateSelectInput(session, "lev3_tab",
      choices = switch(input$lev1_tab,
        "orfs"    = "This subsection does not have units",
        "contigs" = "This subsection does not have units",
        "bins"    = "This subsection does not have units",
        names(reactiveData$SQM[[input$lev1_tab]][[input$lev2_tab]])
      ),
      selected = switch(input$lev1_tab,
        "orfs"    = "This subsection does not have units",
        "contigs" = "This subsection does not have units",
        "bins"    = "This subsection does not have units",
        names(reactiveData$SQM[[input$lev1_tab]][[input$lev2_tab]][1])
      )
    )
  }) # Close lev1+2 observer

  observeEvent(reactiveData$SQM, {
    updateSelectizeInput(session, "cols_tab",
      choices  = reactiveData$SQM[["misc"]][["samples"]],
      selected = reactiveData$SQM[["misc"]][["samples"]]
    )
  }) # Close SQM observer

  # Output Table ----
  reactTable <- reactive({
    table <- switch(input$lev1_tab,
      "orfs"      = as.data.frame(
        reactiveData$SQM[[input$lev1_tab]][[input$lev2_tab]][, input$cols_tab]),
      "contigs"   = as.data.frame(
        reactiveData$SQM[[input$lev1_tab]][[input$lev2_tab]][, input$cols_tab]),
      "bins"      = as.data.frame(
        reactiveData$SQM[[input$lev1_tab]][[input$lev2_tab]][, input$cols_tab]),
      "taxa"      = as.data.frame(
        reactiveData$SQM[[input$lev1_tab]][[input$lev2_tab]][[input$lev3_tab]][, input$cols_tab]),
      "functions" = as.data.frame(
        reactiveData$SQM[[input$lev1_tab]][[input$lev2_tab]][[input$lev3_tab]][, input$cols_tab])
    )

    if (input$lev2_tab %in% c("KEGG", "COG")) {
      annot <- switch(input$lev2_tab,
        "KEGG" = "KEGG_names",
        "COG"  = "COG_names")
      fun_names <- reactiveData$SQM[["misc"]][[annot]][rownames(table)]
      fun_names[length(fun_names)] <- "Unmapped"
      names(fun_names)[length(fun_names)] <- "Unmapped"
      table <- cbind(fun_names, table)
      names(table)[1] <- c("Function name")
    }
    table
  })

  output$table <- DT::renderDataTable({
    DT::datatable(reactTable())
  })

  output$tabDown <- downloadHandler(
    filename = function() {
      paste("Table", Sys.time(), ".csv", sep = "")},
    content = function(file) {
      write.csv(reactTable(), file, row.names = TRUE)
      }
  )

  # Update Summary Inputs ----
  observe({
    updateSelectInput(session, "orfs_row1",
                      choices = reactiveData$orfs_st[, 1]
    )
    updateSelectInput(session, "orfs_row2",
                      choices = reactiveData$orfs_st[, 1]
    )
  }) # Close observer

  # Output Summary ----
  # Reads
  reactReadsSum <- reactive({ #create reactive function
    reactiveData$reads_st
  })
  output$reads_sum <- DT::renderDataTable({ #generate output
    DT::datatable(reactReadsSum(), options = list(
      paging = FALSE, searching = FALSE
    )
    )
  })

  # Barplots from Reads
  output$reads_readbar <- renderPlot({
    barplot(as.matrix(
      reactReadsSum()[1, ! names(reactReadsSum()) %in% c("Assembly")]),
      main = "Reads", col = "#008080", las = 2)
  })
  output$reads_basebar <- renderPlot({
    barplot(as.matrix(
      reactReadsSum()[2, ! names(reactReadsSum()) %in% c("Assembly")]),
      main = "Bases", col = "#00CED1", las = 2)
  })

  # Contigs
  reactContigsSum <- reactive({ #create reactive function
    reactiveData$contigs_st
  })
  output$contigs_sum <- DT::renderDataTable({ #generate output
    DT::datatable(reactContigsSum(), options = list(
      paging = FALSE, searching = FALSE
      )
    )
  })

  # Taxa
  reactTaxaSum <- reactive({ #create reactive function
    reactiveData$taxa_st
  })
  output$taxa_sum <- DT::renderDataTable({ #generate output
    DT::datatable(reactTaxaSum(), options = list(
      paging = FALSE, searching = FALSE
      )
    )
  })

  # Orfs
  reactOrfsSum <- reactive({ #create reactive function
    reactiveData$orfs_st
  })
  output$orfs_sum <- DT::renderDataTable({ #generate output
    DT::datatable(reactOrfsSum(), options = list(
      paging = FALSE, searching = FALSE
      )
    )
  })

  # Barplots from Orfs
  output$orfs_bar1 <- renderPlot({
    barplot(as.matrix(reactOrfsSum()[input$orfs_row1, ]),
      main = input$orfs_row1, col = "#FF7F50",
      ylim = c(0, max(reactOrfsSum()[input$orfs_row1, ])), las = 2)
  })
  output$orfs_bar2 <- renderPlot({
    barplot(as.matrix(reactOrfsSum()[input$orfs_row2, ]),
      main = input$orfs_row2, col = "#4169E1", las = 2)
  })

  # Bins
  reactBinsSum <- reactive({ #create reactive function
    reactiveData$bins_st
  })
  output$bins_sum <- DT::renderDataTable({ #generate output
    DT::datatable(reactBinsSum(), options = list(
      paging = FALSE, searching = FALSE
    )
    )
  })
} # Close server
