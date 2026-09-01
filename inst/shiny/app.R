library(shiny)
library(glycoPathDesign)

`%||%` <- function(x, y) if (is.null(x)) y else x

default_panel <- c("M9", "M8", "M5", "M5Gn", "M3Gn", "M3Gn2", "M3Gn3", "G1", "G2S1")

reaction_information <- data.frame(
  rate_class = c("MAN1", "MGAT1", "MAN2", "MGAT2", "MGAT4", "B4GALT", "ST6GAL", "FUT8"),
  plain_name = c(
    "Early mannose trimming", "Initiation of complex N-glycan processing",
    "Medial-Golgi mannose trimming", "Formation of the second GlcNAc antenna",
    "Addition of a further GlcNAc branch", "Terminal galactosylation",
    "Terminal alpha-2,6 sialylation", "Core fucosylation"
  ),
  stringsAsFactors = FALSE
)

format_number <- function(x, digits = 3) {
  if (!length(x) || is.na(x)) return("Not available")
  if (!is.finite(x)) return("Not finite")
  formatC(x, digits = digits, format = "fg", big.mark = ",")
}

rate_labels <- function(classes) {
  descriptions <- reaction_information$plain_name[match(classes, reaction_information$rate_class)]
  labels <- ifelse(is.na(descriptions), classes, paste0(classes, " — ", descriptions))
  setNames(classes, labels)
}

decision_row <- function(state, requirement, result, interpretation) {
  label <- c(pass = "Meets", fail = "Does not meet", pending = "Not assessed", info = "Report")[[state]]
  tags$tr(
    tags$td(tags$strong(requirement)),
    tags$td(class = paste("decision-state", state), label),
    tags$td(class = "decision-result", result),
    tags$td(interpretation)
  )
}

registered_pathway <- canonical_pathway()

registered_reactions <- data.frame(
  From = registered_pathway$edges$from,
  To = registered_pathway$edges$to,
  `Rate class` = registered_pathway$edges$rate_class,
  Weight = registered_pathway$edges$weight,
  check.names = FALSE, stringsAsFactors = FALSE
)

registered_reaction_choices <- setNames(
  seq_len(nrow(registered_reactions)),
  paste0(registered_reactions$From, " → ", registered_reactions$To,
         "  —  ", registered_reactions[["Rate class"]])
)

default_builder_reactions <- function() registered_reactions

default_builder_states <- function() {
  states <- registered_pathway$nodes$id
  measured <- registered_pathway$observations$glycoform[
    match(states, registered_pathway$observations$node)
  ]
  data.frame(
    State = states, `Measured class` = measured,
    Entry = unname(registered_pathway$entry[states]),
    Secretion = unname(registered_pathway$secretion[states]),
    check.names = FALSE, stringsAsFactors = FALSE
  )
}

ui <- fluidPage(
  tags$head(
    tags$meta(name = "viewport", content = "width=device-width, initial-scale=1"),
    tags$style(HTML("
      :root { --ink:#151515; --muted:#5f6264; --line:#cfd1d2; --soft:#f5f5f3;
              --paper:#ffffff; --pass:#355f4a; --warn:#7b5a20; --fail:#7b3732; }
      html, body { background:#fff; color:var(--ink); font-size:16px; }
      body { font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Arial,sans-serif; }
      .container-fluid { max-width:1500px; padding:0 32px 48px; }
      .app-header { margin:0 -32px; padding:22px 34px 18px; border-bottom:1px solid var(--line); }
      .app-header h1 { margin:0; font-size:28px; font-weight:700; letter-spacing:-.02em; }
      .app-header p { margin:5px 0 0; color:var(--muted); font-size:15px; }
      h2 { font-size:22px; margin:0 0 16px; } h3 { font-size:16px; margin:22px 0 10px; }
      p { line-height:1.48; }
      .nav-tabs { border-bottom:1px solid var(--line); margin:0 -32px 24px; padding:0 32px; }
      .nav-tabs>li>a { color:#353535; border:0!important; border-radius:0; padding:13px 18px 11px; }
      .nav-tabs>li.active>a,.nav-tabs>li.active>a:focus { color:#111; font-weight:650; border-bottom:3px solid #222!important; }
      .design-layout { display:grid; grid-template-columns:minmax(330px,420px) minmax(560px,1fr); gap:30px; }
      .input-pane { padding-right:28px; border-right:1px solid var(--line); }
      .result-pane { min-width:0; }
      .section-rule { border-top:1px solid var(--line); padding-top:18px; margin-top:20px; }
      .model-summary { display:grid; grid-template-columns:repeat(3,1fr); border:1px solid var(--line); margin:12px 0; }
      .model-summary div { padding:10px 8px; text-align:center; border-right:1px solid var(--line); }
      .model-summary div:last-child { border-right:0; }
      .model-summary strong { display:block; font-size:20px; }
      .model-summary span { color:var(--muted); font-size:12px; }
      .form-group { margin-bottom:13px; }
      .control-label { font-size:14px; }
      .selectize-input { border-color:#aeb1b2; border-radius:2px; box-shadow:none; }
      .selectize-control.multi .selectize-input>div { background:#ececeb; border:0; border-radius:2px; color:#111; }
      .btn { border-radius:2px; box-shadow:none!important; }
      .btn-primary { background:#222; border-color:#222; font-weight:650; }
      .btn-primary:hover,.btn-primary:focus { background:#000; border-color:#000; }
      .btn-default { border-color:#a8aaab; background:#fff; }
      .action-links { display:flex; gap:14px; margin:-4px 0 8px; font-size:13px; }
      .action-links a { color:#343434; text-decoration:underline; cursor:pointer; }
      .small-note { color:var(--muted); font-size:13px; line-height:1.45; }
      .builder-block { border-left:3px solid #777; padding-left:12px; margin:12px 0; }
      .builder-block h4 { font-size:14px; margin:14px 0 7px; }
      .builder-actions { display:flex; flex-wrap:wrap; gap:7px; margin:7px 0 12px; }
      .builder-actions .btn { font-size:12px; padding:5px 8px; }
      .dataTables_wrapper { font-size:12px; margin-bottom:5px; }
      table.dataTable { width:100%!important; table-layout:fixed; }
      table.dataTable tbody td { padding:6px 7px!important; white-space:normal!important; word-break:break-word; }
      table.dataTable thead th { padding:7px!important; white-space:normal!important; }
      details { border-top:1px solid var(--line); margin-top:14px; padding-top:11px; }
      summary { cursor:pointer; font-weight:600; font-size:14px; }
      .snfg-preview { display:flex; flex-wrap:wrap; gap:5px; margin:8px 0 2px; }
      .snfg-preview img { width:92px; height:auto; border:1px solid #d5d6d7; }
      .token { display:inline-block; border:1px solid #bbb; padding:4px 7px; margin:3px; font-size:13px; }
      .effective-cv { margin-top:-4px; color:var(--muted); font-size:13px; }
      .result-header { border-left:5px solid #555; padding:16px 19px; background:var(--soft); margin-bottom:20px; }
      .result-header.pass { border-color:var(--pass); }
      .result-header.warn { border-color:var(--warn); }
      .result-header.fail { border-color:var(--fail); }
      .result-header h2 { margin:0 0 5px; font-size:24px; }
      .result-header p { margin:0; color:#383838; }
      .result-context { display:flex; flex-wrap:wrap; gap:8px 18px; color:var(--muted); font-size:13px; margin-bottom:14px; }
      .decision-table { width:100%; border-collapse:collapse; background:#fff; }
      .decision-table th { text-align:left; font-size:12px; text-transform:uppercase; letter-spacing:.04em; color:var(--muted); border-bottom:2px solid #777; padding:8px 9px; }
      .decision-table td { padding:12px 9px; border-bottom:1px solid var(--line); vertical-align:top; line-height:1.4; }
      .decision-table td:nth-child(1) { width:19%; }
      .decision-table td:nth-child(2) { width:15%; font-weight:650; }
      .decision-table td:nth-child(3) { width:22%; }
      .decision-state.pass { color:var(--pass); }
      .decision-state.fail { color:var(--fail); }
      .decision-state.pending { color:var(--warn); }
      .decision-state.info { color:#333; }
      .next-action { border-top:2px solid #333; margin-top:22px; padding-top:15px; }
      .next-action h3 { margin:0 0 7px; font-size:17px; }
      .next-action p { margin:0 0 12px; }
      .candidate-result { margin-top:14px; border:1px solid var(--line); padding:15px; }
      .candidate-result h3 { margin:0 0 8px; }
      .candidate-panel { font-family:ui-monospace,SFMono-Regular,Menlo,monospace; line-height:1.6; }
      .technical-grid { display:grid; grid-template-columns:repeat(4,1fr); border:1px solid var(--line); margin:14px 0; }
      .technical-grid div { padding:12px; border-right:1px solid var(--line); }
      .technical-grid div:last-child { border-right:0; }
      .technical-grid strong { display:block; font-size:20px; }
      .technical-grid span { display:block; color:var(--muted); font-size:12px; margin-top:3px; }
      .conditional-note { color:var(--muted); font-size:13px; border-left:3px solid #999; padding-left:10px; margin-top:16px; }
      .advanced-layout { display:grid; grid-template-columns:320px minmax(560px,1fr); gap:30px; }
      .advanced-controls { border-right:1px solid var(--line); padding-right:26px; }
      .plain-section { border:1px solid var(--line); padding:18px; margin-bottom:18px; }
      .plain-section h2,.plain-section h3 { margin-top:0; }
      .definition-list dt { margin-top:13px; } .definition-list dd { margin-left:0; color:#3e3e3e; }
      .placeholder { border:1px solid var(--line); padding:30px; color:var(--muted); }
      .shiny-notification { border-radius:2px; }
      table { background:#fff; }
      @media (max-width:980px) {
        .design-layout,.advanced-layout { grid-template-columns:1fr; }
        .input-pane,.advanced-controls { border-right:0; padding-right:0; border-bottom:1px solid var(--line); padding-bottom:22px; }
        .decision-table { font-size:14px; }
      }
      @media (max-width:620px) {
        .container-fluid { padding-left:18px; padding-right:18px; }
        .app-header,.nav-tabs { margin-left:-18px; margin-right:-18px; padding-left:18px; padding-right:18px; }
        .model-summary,.technical-grid { grid-template-columns:1fr 1fr; }
        .decision-table th:nth-child(4),.decision-table td:nth-child(4) { display:none; }
      }
    "))
  ),
  div(class = "app-header",
      tags$h1("glycoPathDesign"),
      tags$p("Prospective design checks for glycoform biosynthetic inverse problems")),
  tabsetPanel(
    id = "app_mode",
    tabPanel("Design check", value = "design",
      div(class = "design-layout",
        div(class = "input-pane",
          tags$h2("Design inputs"),
          tags$h3("1  Pathway model"),
          radioButtons(
            "pathway_source", NULL,
            choices = c("Manuscript example" = "canonical", "Build a pathway" = "builder",
                        "Upload CSV" = "upload"),
            selected = "canonical", inline = FALSE
          ),
          conditionalPanel(
            "input.pathway_source == 'builder'",
            div(class = "builder-block",
              radioButtons(
                "builder_input_mode", "Reaction input",
                choices = c("Registered N-glycan reactions" = "guided",
                            "Custom definitions" = "custom"),
                selected = "guided", inline = FALSE
              ),
              conditionalPanel(
                "input.builder_input_mode == 'guided'",
                tags$p(class = "small-note",
                       "Choose a biochemically registered reaction. Source state, product state and effective rate class are assigned together."),
                selectInput("registered_reaction", "Registered reaction",
                            choices = registered_reaction_choices),
                actionButton("add_registered_reaction", "Add registered reaction",
                             class = "btn-default")
              ),
              conditionalPanel(
                "input.builder_input_mode == 'custom'",
                tags$p(class = "small-note",
                       "Enter identifiers for a reaction outside the registered model. Reuse a rate-class name when several edges share one effective parameter."),
                fluidRow(
                  column(6, textInput("custom_from", "From", "",
                                      placeholder = "e.g. precursor_A")),
                  column(6, textInput("custom_to", "To", "",
                                      placeholder = "e.g. product_A")),
                  column(8, textInput("custom_rate_class", "Effective rate class", "",
                                      placeholder = "e.g. branch_A")),
                  column(4, numericInput("custom_weight", "Weight", 1,
                                         min = 0.001, step = 0.1))
                ),
                actionButton("add_custom_reaction", "Add custom reaction",
                             class = "btn-default")
              ),
              tags$h4("Current reactions"),
              DT::DTOutput("reaction_editor"),
              div(class = "builder-actions",
                  actionButton("delete_reaction", "Delete selected", class = "btn-default"),
                  actionButton("clear_reactions", "Start blank", class = "btn-default")),
              tags$h4("States and measurement mapping"),
              DT::DTOutput("state_editor"),
              tags$p(class = "small-note", "Measured class may equal the state name or combine several states under one reported glycoform class."),
              div(class = "builder-actions",
                  actionButton("reset_builder", "Restore registered model", class = "btn-default"),
                  downloadButton("download_builder", "Export pathway files", class = "btn-default"))
            )
          ),
          conditionalPanel(
            "input.pathway_source == 'upload'",
            fileInput(
              "pathway_files", "Pathway CSV files", multiple = TRUE, accept = ".csv",
              placeholder = "Select nodes, edges and optional mapping/boundary files"
            ),
            downloadButton("download_templates", "Download CSV template", class = "btn-sm"),
            tags$details(
              tags$summary("Required columns"),
              tags$p(class = "small-note",
                     "nodes.csv: id, label; edges.csv: from, to, rate_class, weight. observations.csv maps node to glycoform. boundary.csv contains node, entry and secretion. The mapping and boundary files are optional."),
              textInput("custom_entry", "Entry node if boundary.csv is absent", ""),
              numericInput("custom_secretion", "Uniform secretion if boundary.csv is absent", 0.5, min = 0.001)
            )
          ),
          uiOutput("model_summary"),
          div(class = "section-rule",
            tags$h3("2  Reaction classes to distinguish"),
            uiOutput("rate_selector"),
            div(class = "action-links", actionLink("select_all_rates", "Select all"), actionLink("clear_rates", "Clear"))
          ),
          div(class = "section-rule",
            tags$h3("3  Glycoforms to measure together"),
            uiOutput("panel_selector"),
            div(class = "action-links", actionLink("select_all_glycoforms", "Select all"), actionLink("clear_glycoforms", "Clear")),
            uiOutput("selected_glycoforms")
          ),
          div(class = "section-rule",
            tags$h3("4  Expected analytical precision"),
            fluidRow(
              column(6, numericInput("assay_cv", "Single-measurement CV (%)", 5, min = 0.01, max = 100, step = 0.1)),
              column(6, numericInput("replicates", "Independent replicates", 1, min = 1, max = 100, step = 1))
            ),
            uiOutput("effective_cv_text"),
            tags$details(
              tags$summary("Rate-precision target"),
              numericInput("target_sd", "Target worst-direction log-rate SD", 0.25, min = 0.01, max = 2, step = 0.05),
              tags$p(class = "small-note", "A value of 0.25 corresponds to about 28% multiplicative uncertainty for one standard deviation.")
            )
          ),
          tags$p(class = "conditional-note",
                 "Results update automatically and are local to the supplied graph, reaction grouping, operating point, observation mapping, secretion model and error assumptions.")
        ),
        div(class = "result-pane",
          uiOutput("result_context"),
          uiOutput("result_header"),
          uiOutput("decision_table"),
          uiOutput("next_action"),
          uiOutput("candidate_result"),
          tags$details(
            tags$summary("Numerical diagnostics"),
            uiOutput("technical_metrics"),
            fluidRow(
              column(6, plotOutput("singular_plot", height = 300)),
              column(6, plotOutput("sensitivity_plot", height = 300))
            )
          ),
          br(), downloadButton("download_report", "Download design report")
        )
      )
    ),
    tabPanel("Robustness & fit", value = "advanced",
      div(class = "advanced-layout",
        div(class = "advanced-controls",
          tags$h2("Additional checks"),
          tags$h3("Operating-point robustness"),
          numericInput("n_points", "Operating points", 100, min = 10, max = 2000),
          numericInput("log_rate_sd", "Log-rate sampling SD", 0.7, min = 0.01),
          numericInput("robust_seed", "Random seed", 1, min = 1),
          actionButton("run_robustness", "Run robustness check", class = "btn-primary btn-block"),
          div(class = "section-rule",
            tags$h3("Observed composition"),
            fileInput("observed_file", "CSV with glycoform and proportion", accept = ".csv"),
            numericInput("fit_starts", "Optimisation starts", 8, min = 1, max = 100),
            actionButton("run_fit", "Fit pathway", class = "btn-primary btn-block")
          )
        ),
        div(
          div(class = "plain-section",
            tags$h2("Pathway used for these checks"),
            plotOutput("pathway_plot", height = 360)
          ),
          div(class = "plain-section",
            tags$h2("Robustness"),
            uiOutput("robustness_summary"),
            plotOutput("robustness_plot", height = 310)
          ),
          div(class = "plain-section",
            tags$h2("Forward-model fit"),
            uiOutput("fit_summary"),
            plotOutput("fit_plot", height = 310),
            tableOutput("rate_table")
          )
        )
      )
    ),
    tabPanel("Definitions", value = "definitions",
      div(class = "plain-section", style = "max-width:900px;margin:0 auto;",
        tags$h2("What the checks mean"),
        tags$dl(class = "definition-list",
          tags$dt("Independent composition dimensions"),
          tags$dd("A jointly normalised composition with m measured classes contains at most m - 1 independent directions."),
          tags$dt("Local rank"),
          tags$dd("The number of selected reaction directions that produce distinguishable local glycoform-response patterns."),
          tags$dt("Condition number"),
          tags$dd("The ratio between the most and least observable selected reaction directions. Larger values indicate greater sensitivity to noise."),
          tags$dt("Noise amplification g"),
          tags$dd("The worst-case conversion from relative glycoform-intensity noise to log-rate uncertainty. Smaller values are better."),
          tags$dt("Forward-model fit"),
          tags$dd("Whether the chosen pathway can reproduce an observed glycoform composition. Fit and rate recoverability are separate requirements.")
        ),
        tags$p(class = "conditional-note",
               "The tool evaluates identifiability conditional on a specified model. It does not prove that the model is biologically complete or that an inferred effective rate corresponds to a single enzyme.")
      )
    )
  )
)

server <- function(input, output, session) {
  reaction_rows <- reactiveVal(default_builder_reactions())
  state_rows <- reactiveVal(default_builder_states())
  state_proxy <- DT::dataTableProxy("state_editor")

  output$reaction_editor <- DT::renderDT({
    DT::datatable(
      reaction_rows(), rownames = FALSE, selection = "single",
      options = list(dom = "t", paging = FALSE, ordering = FALSE, autoWidth = FALSE,
                     scrollY = "220px", scrollCollapse = TRUE,
                     columnDefs = list(list(width = "26%", targets = c(0, 1)),
                                       list(width = "33%", targets = 2),
                                       list(width = "15%", targets = 3)))
    )
  }, server = FALSE)

  output$state_editor <- DT::renderDT({
    DT::datatable(
      state_rows(), rownames = FALSE, selection = "none",
      editable = list(target = "cell", disable = list(columns = c(0))),
      options = list(dom = "t", paging = FALSE, ordering = FALSE, autoWidth = FALSE,
                     columnDefs = list(list(width = "30%", targets = c(0, 1)),
                                       list(width = "20%", targets = c(2, 3))))
    )
  }, server = FALSE)

  observeEvent(input$state_editor_cell_edit, {
    updated <- DT::editData(state_rows(), input$state_editor_cell_edit,
                            proxy = state_proxy, rownames = FALSE, resetPaging = FALSE)
    updated[["Measured class"]] <- trimws(as.character(updated[["Measured class"]]))
    updated$Entry <- suppressWarnings(as.numeric(updated$Entry))
    updated$Secretion <- suppressWarnings(as.numeric(updated$Secretion))
    state_rows(updated)
  })

  observeEvent(reaction_rows(), {
    reactions <- reaction_rows()
    states <- unique(c(reactions$From[nzchar(reactions$From)], reactions$To[nzchar(reactions$To)]))
    old <- state_rows()
    if (!length(states)) {
      state_rows(old[0, , drop = FALSE])
      return()
    }
    updated <- data.frame(
      State = states, `Measured class` = states, Entry = 0, Secretion = 0.5,
      check.names = FALSE, stringsAsFactors = FALSE
    )
    old_index <- match(states, old$State)
    keep <- !is.na(old_index)
    if (any(keep)) {
      updated[keep, c("Measured class", "Entry", "Secretion")] <-
        old[old_index[keep], c("Measured class", "Entry", "Secretion"), drop = FALSE]
    }
    if (!any(is.finite(updated$Entry) & updated$Entry > 0)) updated$Entry[[1]] <- 1
    state_rows(updated)
  }, ignoreInit = TRUE)

  observeEvent(input$add_registered_reaction, {
    index <- suppressWarnings(as.integer(input$registered_reaction))
    if (!length(index) || is.na(index) || index < 1L || index > nrow(registered_reactions)) return()
    candidate <- registered_reactions[index, , drop = FALSE]
    current <- reaction_rows()
    duplicate <- nrow(current) > 0 && any(
      current$From == candidate$From & current$To == candidate$To &
        current[["Rate class"]] == candidate[["Rate class"]]
    )
    if (duplicate) {
      showNotification("That registered reaction is already in the pathway.", type = "message")
      return()
    }
    reaction_rows(rbind(current, candidate))
  })

  observeEvent(input$add_custom_reaction, {
    from <- trimws(input$custom_from %||% "")
    to <- trimws(input$custom_to %||% "")
    rate <- trimws(input$custom_rate_class %||% "")
    weight <- suppressWarnings(as.numeric(input$custom_weight %||% NA_real_))
    if (!nzchar(from) || !nzchar(to) || !nzchar(rate)) {
      showNotification("From, To and effective rate class are required.", type = "error")
      return()
    }
    if (identical(from, to)) {
      showNotification("From and To must be different states.", type = "error")
      return()
    }
    if (!is.finite(weight) || weight <= 0) {
      showNotification("Weight must be a positive number.", type = "error")
      return()
    }
    candidate <- data.frame(
      From = from, To = to, `Rate class` = rate, Weight = weight,
      check.names = FALSE, stringsAsFactors = FALSE
    )
    reaction_rows(rbind(reaction_rows(), candidate))
  })

  observeEvent(input$delete_reaction, {
    selected <- input$reaction_editor_rows_selected
    if (!length(selected)) {
      showNotification("Select one reaction row to delete.", type = "message")
      return()
    }
    current <- reaction_rows()
    if (nrow(current) <= 1) {
      showNotification("A pathway requires at least one reaction.", type = "error")
      return()
    }
    reaction_rows(current[-selected[[1]], , drop = FALSE])
  })

  observeEvent(input$clear_reactions, {
    reaction_rows(registered_reactions[0, , drop = FALSE])
  })

  observeEvent(input$reset_builder, {
    reaction_rows(default_builder_reactions())
    state_rows(default_builder_states())
  })

  builder_pathway <- reactive({
    reactions <- reaction_rows()
    validate(need(nrow(reactions) > 0, "Add at least one reaction."))
    validate(need(all(nzchar(reactions$From) & nzchar(reactions$To) & nzchar(reactions[["Rate class"]])),
                  "Every reaction requires source, product and rate-class names."))
    validate(need(all(is.finite(reactions$Weight) & reactions$Weight > 0),
                  "Every reaction weight must be a positive number."))

    settings <- state_rows()
    validate(need(nrow(settings) >= 2, "The reaction table must define at least two states."))
    validate(need(all(nzchar(settings[["Measured class"]])),
                  "Every state requires a measured-class name."))
    validate(need(all(is.finite(settings$Entry) & settings$Entry >= 0) && sum(settings$Entry) > 0,
                  "Entry values must be non-negative and at least one must be positive."))
    validate(need(all(is.finite(settings$Secretion) & settings$Secretion > 0),
                  "Secretion values must be positive."))

    nodes <- data.frame(id = settings$State, label = settings$State, stringsAsFactors = FALSE)
    edges <- data.frame(
      from = reactions$From, to = reactions$To,
      rate_class = reactions[["Rate class"]], weight = reactions$Weight,
      stringsAsFactors = FALSE
    )
    observations <- data.frame(
      node = settings$State, glycoform = settings[["Measured class"]], stringsAsFactors = FALSE
    )
    glyco_pathway(
      nodes, edges, observations,
      entry = setNames(settings$Entry, settings$State),
      secretion = setNames(settings$Secretion, settings$State),
      metadata = list(name = "Pathway built in glycoPathDesign", conditional = TRUE)
    )
  })

  pathway <- reactive({
    source <- input$pathway_source %||% "canonical"
    if (identical(source, "canonical")) return(canonical_pathway())
    if (identical(source, "builder")) return(builder_pathway())

    files <- input$pathway_files
    req(files)
    base_names <- tolower(basename(files$name))
    validate(need(all(c("nodes.csv", "edges.csv") %in% base_names),
                  "The custom pathway requires files named nodes.csv and edges.csv."))
    read_named <- function(name) {
      path <- files$datapath[match(name, base_names)]
      read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
    }
    nodes <- read_named("nodes.csv")
    edges <- read_named("edges.csv")
    observations <- if ("observations.csv" %in% base_names) read_named("observations.csv") else NULL

    if ("boundary.csv" %in% base_names) {
      boundary <- read_named("boundary.csv")
      validate(need(all(c("node", "entry", "secretion") %in% names(boundary)),
                    "boundary.csv must contain node, entry and secretion."))
      entry <- setNames(boundary$entry, boundary$node)
      secretion <- setNames(boundary$secretion, boundary$node)
    } else {
      declared_entry <- trimws(input$custom_entry %||% "")
      entry <- if (nzchar(declared_entry)) declared_entry else as.character(nodes$id[[1]])
      secretion <- input$custom_secretion %||% 0.5
    }
    glyco_pathway(nodes, edges, observations, entry = entry, secretion = secretion,
                  metadata = list(name = "User-supplied pathway", conditional = TRUE))
  })

  is_canonical <- reactive({
    p <- pathway()
    reference <- canonical_pathway()
    identical(p$nodes$id, reference$nodes$id) && identical(p$rate_classes, reference$rate_classes)
  })

  output$model_summary <- renderUI({
    p <- pathway()
    div(class = "model-summary",
        div(tags$strong(nrow(p$nodes)), tags$span("network states")),
        div(tags$strong(length(p$observed_classes)), tags$span("observable classes")),
        div(tags$strong(length(p$rate_classes)), tags$span("rate classes")))
  })

  output$rate_selector <- renderUI({
    p <- pathway()
    current <- intersect(isolate(input$rate_classes) %||% p$rate_classes, p$rate_classes)
    if (!length(current)) current <- p$rate_classes
    selectizeInput("rate_classes", NULL, choices = rate_labels(p$rate_classes), selected = current,
                   multiple = TRUE, options = list(plugins = list("remove_button"), closeAfterSelect = TRUE))
  })

  output$panel_selector <- renderUI({
    p <- pathway()
    fallback <- if (is_canonical()) intersect(default_panel, p$observed_classes) else p$observed_classes
    current <- intersect(isolate(input$glycoforms) %||% fallback, p$observed_classes)
    selectizeInput("glycoforms", NULL, choices = p$observed_classes, selected = current,
                   multiple = TRUE, options = list(plugins = list("remove_button"), closeAfterSelect = TRUE))
  })

  observeEvent(pathway(), {
    p <- pathway()
    panel <- if (is_canonical()) intersect(default_panel, p$observed_classes) else p$observed_classes
    updateSelectizeInput(session, "rate_classes", choices = rate_labels(p$rate_classes), selected = p$rate_classes, server = TRUE)
    updateSelectizeInput(session, "glycoforms", choices = p$observed_classes, selected = panel, server = TRUE)
  }, ignoreInit = TRUE)

  observeEvent(input$select_all_rates, {
    updateSelectizeInput(session, "rate_classes", selected = pathway()$rate_classes, server = TRUE)
  })
  observeEvent(input$clear_rates, updateSelectizeInput(session, "rate_classes", selected = character(), server = TRUE))
  observeEvent(input$select_all_glycoforms, {
    updateSelectizeInput(session, "glycoforms", selected = pathway()$observed_classes, server = TRUE)
  })
  observeEvent(input$clear_glycoforms, updateSelectizeInput(session, "glycoforms", selected = character(), server = TRUE))

  output$selected_glycoforms <- renderUI({
    selected <- intersect(input$glycoforms %||% character(), pathway()$observed_classes)
    if (!length(selected)) return(tags$p(class = "small-note", "Select at least two glycoform classes."))
    if (is_canonical()) {
      return(tags$details(
        tags$summary(paste("Show", length(selected), "selected SNFG structures")),
        div(class = "snfg-preview", lapply(selected, function(g) {
          tags$img(src = paste0("snfg_cards/", g, ".png"), alt = paste(g, "SNFG structure"), title = g)
        }))
      ))
    }
    div(lapply(selected, function(g) tags$span(class = "token", g)))
  })

  effective_cv <- reactive({
    cv <- as.numeric(input$assay_cv %||% 5) / 100
    n <- max(1, as.integer(input$replicates %||% 1))
    cv / sqrt(n)
  })

  output$effective_cv_text <- renderUI({
    tags$p(class = "effective-cv",
           paste0("Effective CV after averaging: ", format_number(100 * effective_cv()), "%"))
  })

  design_audit <- reactive({
    p <- pathway()
    panel <- intersect(input$glycoforms %||% character(), p$observed_classes)
    rates <- intersect(input$rate_classes %||% character(), p$rate_classes)
    req(length(panel) >= 2, length(rates) >= 1)
    check_panel(
      p, panel,
      target_sd = as.numeric(input$target_sd %||% 0.25),
      assay_cv = effective_cv(), rate_classes = rates
    )
  })

  output$result_context <- renderUI({
    p <- pathway()
    selected_glycoforms <- intersect(input$glycoforms %||% character(), p$observed_classes)
    selected_rates <- intersect(input$rate_classes %||% character(), p$rate_classes)
    div(class = "result-context",
        tags$span(tags$strong(if (is_canonical()) "Model: " else "Custom model: "),
                  p$metadata$name %||% if (is_canonical()) "Canonical N-glycan example" else "User-supplied pathway"),
        tags$span(tags$strong("Operating point: "), "reference rates"),
        tags$span(tags$strong("Selected: "),
                  paste(length(selected_glycoforms), "glycoforms /",
                        length(selected_rates), "rate classes")))
  })

  verdict <- reactive({
    x <- design_audit()
    if (!x$dimension_sufficient) {
      return(list(state = "fail", title = "Too few independent measurements",
                  text = paste0("At least ", x$n_rate_classes + 1,
                                " jointly normalised glycoform classes are required for ",
                                x$n_rate_classes, " selected rate classes.")))
    }
    if (!x$full_rank) {
      return(list(state = "fail", title = "Some selected reactions cannot be separated",
                  text = paste0("This panel distinguishes ", x$rank, " of ", x$n_rate_classes,
                                " local rate directions at the reference operating point.")))
    }
    if (!isTRUE(x$precision_target_met)) {
      return(list(state = "warn", title = "The panel separates the reactions, but precision is insufficient",
                  text = paste0("The effective CV is ", format_number(100 * effective_cv()),
                                "%; the selected target requires ", format_number(100 * x$max_assay_cv), "% or lower.")))
    }
    list(state = "pass", title = "The measurement design is locally sufficient",
         text = paste0("The panel separates all ", x$n_rate_classes,
                       " selected rate directions and meets the stated precision target at the reference operating point."))
  })

  output$result_header <- renderUI({
    p <- pathway()
    if (length(intersect(input$glycoforms %||% character(), p$observed_classes)) < 2 ||
        length(intersect(input$rate_classes %||% character(), p$rate_classes)) < 1) {
      return(div(class = "placeholder", "Select at least one rate class and two glycoform classes to evaluate the design."))
    }
    v <- verdict()
    div(class = paste("result-header", v$state), tags$h2(v$title), tags$p(v$text))
  })

  output$decision_table <- renderUI({
    req(length(input$glycoforms) >= 2, length(input$rate_classes) >= 1)
    x <- design_audit()
    dimension_result <- paste0(x$independent_dimensions, " available / ", x$n_rate_classes, " required")
    rank_result <- paste0(x$rank, " distinguishable / ", x$n_rate_classes, " selected")
    precision_result <- paste0(format_number(100 * effective_cv()), "% entered / ",
                               format_number(100 * x$max_assay_cv), "% maximum")
    condition_result <- if (x$full_rank) paste0(format_number(x$condition_number), "-fold sensitivity range") else "Not defined without full rank"
    tags$table(class = "decision-table",
      tags$thead(tags$tr(tags$th("Requirement"), tags$th("Status"), tags$th("Result"), tags$th("Meaning"))),
      tags$tbody(
        decision_row(if (x$dimension_sufficient) "pass" else "fail", "Panel dimension", dimension_result,
                     "Normalisation removes one independent direction."),
        decision_row(if (x$full_rank) "pass" else "fail", "Reaction separation", rank_result,
                     "Distinct local response patterns are required for rate attribution."),
        decision_row(if (isTRUE(x$precision_target_met)) "pass" else "fail", "Analytical precision", precision_result,
                     "The maximum CV is conditional on the chosen log-rate uncertainty target."),
        decision_row("info", "Conditioning", condition_result,
                     "Larger values indicate a greater difference between the best and worst observed rate directions."),
        decision_row("pending", "Forward-model fit", "Requires pilot or observed data",
                     "A full-rank panel is useful only if the specified pathway represents the measured composition.")
      )
    )
  })

  output$next_action <- renderUI({
    req(length(input$glycoforms) >= 2, length(input$rate_classes) >= 1)
    x <- design_audit()
    if (!x$dimension_sufficient) {
      deficit <- x$n_rate_classes + 1L - x$panel_size
      text <- paste0("Add at least ", deficit, " glycoform class", if (deficit == 1) "" else "es",
                     ", or reduce the number of rate classes to be estimated together.")
    } else if (!x$full_rank) {
      text <- "Change the measured glycoform set or combine reaction classes whose response patterns cannot be separated."
    } else if (!isTRUE(x$precision_target_met)) {
      required_n <- ceiling((as.numeric(input$assay_cv) / (100 * x$max_assay_cv))^2)
      replicate_text <- if (is.finite(required_n) && required_n > 1 && required_n <= 10000) {
        paste0(" Under independent analytical error, approximately ", required_n,
               " replicates would reach this effective CV; systematic error will not average away.")
      } else ""
      text <- paste0("Improve the effective CV to ", format_number(100 * x$max_assay_cv), "% or select a better-conditioned panel.", replicate_text)
    } else {
      text <- "Test the panel across plausible operating points, then use pilot data to confirm forward-model fit before interpreting individual rates."
    }
    div(class = "next-action", tags$h3("Recommended next step"), tags$p(text),
        actionButton("search_alternative", "Search alternative panels", class = "btn-default"))
  })

  candidate_visible <- reactiveVal(FALSE)
  observeEvent(input$search_alternative, candidate_visible(TRUE), ignoreInit = TRUE)

  candidate_search <- eventReactive(input$search_alternative, {
    x <- design_audit()
    size <- min(length(pathway()$observed_classes), max(x$panel_size, x$n_rate_classes + 1L))
    withProgress(message = "Evaluating candidate panels", value = 0.3, {
      panels <- enumerate_panels(
        pathway(), panel_size = size, max_panels = 5000, seed = 1,
        target_sd = x$target_sd, assay_cv = effective_cv(), rate_classes = x$rate_classes
      )
      incProgress(0.6)
      best <- optimise_panels(panels, objective = "amplification", n = 1L)
      if (!nrow(best)) return(list(best = NULL, panels = panels))
      list(best = best[1, , drop = FALSE], panels = panels,
           panel = strsplit(best$panel[[1]], " \\| ")[[1]])
    })
  }, ignoreInit = TRUE)

  output$candidate_result <- renderUI({
    req(candidate_visible())
    result <- candidate_search()
    if (is.null(result$best)) {
      return(div(class = "candidate-result", tags$h3("Alternative-panel search"),
                 tags$p("No locally full-rank panel was found among the evaluated candidates.")))
    }
    note <- if (result$panels$sampled) {
      paste0("A reproducible sample of ", format(result$panels$evaluated, big.mark = ","), " of ",
             format(result$panels$total_combinations, big.mark = ","), " possible panels was evaluated.")
    } else {
      paste0("All ", format(result$panels$total_combinations, big.mark = ","), " possible panels of this size were evaluated.")
    }
    div(class = "candidate-result",
        tags$h3("Best full-rank panel found"),
        tags$p(class = "candidate-panel", paste(result$panel, collapse = " · ")),
        tags$p(paste0("Maximum supported CV: ", format_number(100 * result$best$max_assay_cv),
                      "% · condition number: ", format_number(result$best$condition_number), ".")),
        tags$p(class = "small-note", note),
        actionButton("use_alternative", "Use this panel", class = "btn-default"))
  })

  observeEvent(input$use_alternative, {
    result <- candidate_search()
    if (length(result$panel)) {
      updateSelectizeInput(session, "glycoforms", choices = pathway()$observed_classes,
                           selected = result$panel, server = TRUE)
      candidate_visible(FALSE)
    }
  })

  output$technical_metrics <- renderUI({
    x <- design_audit()
    div(class = "technical-grid",
        div(tags$strong(paste0(x$rank, "/", x$n_rate_classes)), tags$span("local rank")),
        div(tags$strong(format_number(x$condition_number)), tags$span("condition number")),
        div(tags$strong(format_number(x$amplification)), tags$span("noise amplification g")),
        div(tags$strong(paste0(format_number(100 * x$max_assay_cv), "%")), tags$span("maximum CV")))
  })
  output$singular_plot <- renderPlot(plot(design_audit()))
  output$sensitivity_plot <- renderPlot(plot_sensitivity(design_audit(), main = "Rate-response patterns"), res = 120)
  output$download_report <- downloadHandler(
    filename = function() "glycoPathDesign_panel_audit.html",
    content = function(file) write_design_report(design_audit(), file)
  )

  output$download_templates <- downloadHandler(
    filename = function() "glycoPathDesign_pathway_template.zip",
    content = function(file) {
      directory <- tempfile("glycoPathDesign-template-")
      dir.create(directory)
      nodes <- data.frame(id = c("A", "B", "C"), label = c("precursor", "intermediate", "product"))
      edges <- data.frame(from = c("A", "B"), to = c("B", "C"), rate_class = c("r1", "r2"), weight = 1)
      observations <- data.frame(node = c("A", "B", "C"), glycoform = c("gA", "gB", "gC"))
      boundary <- data.frame(node = c("A", "B", "C"), entry = c(1, 0, 0), secretion = 0.5)
      write.csv(nodes, file.path(directory, "nodes.csv"), row.names = FALSE)
      write.csv(edges, file.path(directory, "edges.csv"), row.names = FALSE)
      write.csv(observations, file.path(directory, "observations.csv"), row.names = FALSE)
      write.csv(boundary, file.path(directory, "boundary.csv"), row.names = FALSE)
      old <- setwd(directory)
      on.exit(setwd(old), add = TRUE)
      utils::zip(file, c("nodes.csv", "edges.csv", "observations.csv", "boundary.csv"))
    }
  )

  output$download_builder <- downloadHandler(
    filename = function() "glycoPathDesign_built_pathway.zip",
    content = function(file) {
      built <- builder_pathway()
      directory <- tempfile("glycoPathDesign-built-")
      dir.create(directory)
      write.csv(built$nodes, file.path(directory, "nodes.csv"), row.names = FALSE)
      write.csv(built$edges, file.path(directory, "edges.csv"), row.names = FALSE)
      write.csv(built$observations, file.path(directory, "observations.csv"), row.names = FALSE)
      write.csv(data.frame(node = names(built$entry), entry = unname(built$entry),
                           secretion = unname(built$secretion)),
                file.path(directory, "boundary.csv"), row.names = FALSE)
      old <- setwd(directory)
      on.exit(setwd(old), add = TRUE)
      utils::zip(file, c("nodes.csv", "edges.csv", "observations.csv", "boundary.csv"))
    }
  )

  output$pathway_plot <- renderPlot(plot(pathway(), show_rate_labels = TRUE), res = 120)

  robustness <- eventReactive(input$run_robustness, {
    x <- design_audit()
    check_operating_points(
      pathway(), x$panel, n = input$n_points,
      log_rate_sd = input$log_rate_sd, seed = input$robust_seed,
      target_sd = input$target_sd, assay_cv = effective_cv(), rate_classes = x$rate_classes
    )
  }, ignoreInit = TRUE)

  output$robustness_summary <- renderUI({
    if ((input$run_robustness %||% 0) == 0) {
      return(tags$p(class = "small-note", "Run this check to assess whether the conclusion persists across plausible rate configurations."))
    }
    x <- robustness()
    div(class = "technical-grid",
        div(tags$strong(paste0(format_number(100 * x$full_rank_fraction), "%")), tags$span("operating points full rank")),
        div(tags$strong(format_number(x$condition_median)), tags$span("median condition number")),
        div(tags$strong(format_number(x$amplification_median)), tags$span("median noise amplification")),
        div(tags$strong(paste0(format_number(100 * x$max_cv_tenth_percentile), "%")), tags$span("10th-percentile maximum CV")))
  })
  output$robustness_plot <- renderPlot(plot(req(robustness())))

  fitted <- eventReactive(input$run_fit, {
    req(input$observed_file)
    p <- pathway()
    panel <- intersect(input$glycoforms %||% character(), p$observed_classes)
    req(length(panel) >= 2)
    observed <- read.csv(input$observed_file$datapath, stringsAsFactors = FALSE, check.names = FALSE)
    fit_pathway(p, observed, panel = panel, n_starts = input$fit_starts)
  }, ignoreInit = TRUE)

  output$fit_summary <- renderUI({
    if ((input$run_fit %||% 0) == 0) {
      return(tags$p(class = "small-note", "Upload a within-site composition to assess whether the pathway represents the observed glycoforms."))
    }
    x <- fitted()
    div(class = "technical-grid",
        div(tags$strong(format_number(x$fit_correlation)), tags$span("fit correlation")),
        div(tags$strong(paste0(x$audit$rank, "/", x$audit$n_rate_classes)), tags$span("local rank at fitted point")),
        div(tags$strong(format_number(x$audit$condition_number)), tags$span("condition number at fit")),
        div(tags$strong(format_number(x$objective)), tags$span("fit objective")))
  })
  output$fit_plot <- renderPlot(plot(req(fitted())))
  output$rate_table <- renderTable({
    x <- req(fitted())
    data.frame(rate_class = names(x$rates), effective_rate = unname(x$rates), log_rate = unname(x$log_rates))
  }, digits = 4, rownames = FALSE)
}

shinyApp(ui, server)
