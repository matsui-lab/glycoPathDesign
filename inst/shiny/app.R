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
    "Terminal sialylation (model class)", "Fucosylation (model class)"
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

# Layout depends on the pathway only; selecting rates or glycans never moves it.
pathway_map_coordinates <- function(p, canonical = FALSE) {
  if (canonical) {
    xy <- rbind(
      cbind(x = 90 + 160 * 0:7, y = 105),
      cbind(x = 1210 - 230 * 0:4, y = 350),
      cbind(x = 1210 - 230 * 0:4, y = 600),
      c(90,350), c(1210,850), c(980,850), c(750,850))
    rownames(xy) <- c("M9","M8","M7","M6","M5","M5Gn","M4Gn","M3Gn",
      "M3Gn2","G1","G2","G2S1","G2S2","M3Gn2F","G1F","G2F","G2S1F","G2S2F",
      "G1S1","M3Gn3","Gt2","Gt2S")
    return(xy[p$nodes$id, , drop = FALSE])
  }
  n <- nrow(p$nodes)
  if (requireNamespace("igraph", quietly = TRUE) && nrow(p$edges)) {
    g <- igraph::graph_from_data_frame(p$edges[, c("from","to")], directed = TRUE, vertices = p$nodes$id)
    xy <- tryCatch(igraph::layout_with_sugiyama(g)$layout,
      error = function(e) igraph::layout_in_circle(g))
    xy <- cbind(xy[,2], xy[,1])
  } else {
    theta <- seq(0, 2*pi, length.out = n+1)[seq_len(n)]
    xy <- cbind(cos(theta), sin(theta))
  }
  for (j in 1:2) xy[,j] <- 100 + (xy[,j] - min(xy[,j])) * 190
  rownames(xy) <- p$nodes$id
  xy
}

pathway_map_svg <- function(p, xy, rates, measured, canonical = FALSE, zoom = 100, fit_mode = FALSE, rate_values = NULL, proportions = NULL) {
  esc <- function(x) as.character(htmltools::htmlEscape(as.character(x), attribute = TRUE))
  width <- if (canonical) 1420 else max(700, max(xy[,1])+120)
  height <- if (canonical) 950 else max(300, max(xy[,2])+120)
  pieces <- c(sprintf('<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 %g %g" style="width:%g%%;height:auto;display:block" role="img" aria-label="Current pathway and selected reactions and measurements">', width,height,zoom),
    '<defs><marker id="map-arrow-gray" markerWidth="8" markerHeight="8" refX="7" refY="4" orient="auto"><path d="M0 0 L8 4 L0 8" fill="#bac0c4"/></marker><marker id="map-arrow-active" markerWidth="8" markerHeight="8" refX="7" refY="4" orient="auto"><path d="M0 0 L8 4 L0 8" fill="#853f53"/></marker></defs>')
  # Use one orientation for each unordered pair, including reverse reactions.
  edge_pairs <- vapply(seq_len(nrow(p$edges)), function(i)
    paste(sort(c(p$edges$from[i], p$edges$to[i])), collapse = "\r"), character(1))
  for (i in seq_len(nrow(p$edges))) {
    e <- p$edges[i,]; a <- xy[e$from,]; b <- xy[e$to,]
    x <- a[1]; y <- a[2]; u <- b[1]; v <- b[2]
    active <- e$rate_class %in% rates
    color <- if (active) "#853f53" else "#bac0c4"
    siblings <- which(edge_pairs == edge_pairs[i])
    if (length(siblings) > 1L) {
      delta <- b-a
      scale <- min(.4, 1/max(abs(delta)/c(55, if(canonical) 62 else 34)))
      start <- a+delta*scale; end <- b-delta*scale
      normal <- c(-delta[2], delta[1])/sqrt(sum(delta^2))
      if (e$from > e$to) normal <- -normal
      offset <- (match(i, siblings)-(length(siblings)+1)/2)*100
      control <- (a+b)/2 + normal*offset
      d <- sprintf("M %g %g Q %g %g %g %g", start[1],start[2],control[1],control[2],end[1],end[2])
      midpoint <- (start+2*control+end)/4
      tx <- midpoint[1]; ty <- midpoint[2]-10
    } else if (canonical && e$from == "M3Gn2" && e$to == "M3Gn3") {
      d <- sprintf("M %g %g H 1360 V %g H %g",x+49,y,v,u+55); tx <- 1350; ty <- 750
    } else if (canonical && e$from == "G1" && e$to == "G1S1") {
      d <- sprintf("M %g %g V 235 H %g V %g",x,y-58,u,v-64); tx <- 490; ty <- 226
    } else {
      delta <- b-a
      start_scale <- 1/max(abs(delta)/c(55,if(canonical) 62 else 34))
      start_scale <- min(.4,start_scale)
      start <- a+delta*start_scale; end <- b-delta*start_scale
      d <- sprintf("M %g %g L %g %g",start[1],start[2],end[1],end[2])
      tx <- (x+u)/2 + if (abs(u-x)<1) 38 else 0
      ty <- (y+v)/2 - if (abs(u-x)<1) 0 else 13
    }
    rate_label <- if (fit_mode) paste0(e$rate_class, " = ", if (is.null(rate_values)) "?" else format_number(rate_values[[e$rate_class]])) else e$rate_class
    label <- paste0(rate_label, if (e$weight != 1) paste0(" ×",e$weight) else "")
    label_markup <- esc(label)
    if (fit_mode) {
      ty <- ty - 13
      value <- if (is.null(rate_values)) "?" else format_number(rate_values[[e$rate_class]])
      label_markup <- sprintf('<tspan fill="#727b81">%s</tspan><tspan x="%g" dy="22" font-size="23" font-weight="700" fill="#9a3e51">%s</tspan>', esc(e$rate_class), tx, esc(value))
      if (e$weight != 1) label_markup <- paste0(label_markup, sprintf('<tspan font-size="11"> ×%s</tspan>', esc(e$weight)))
    }
    pieces <- c(pieces,sprintf('<g class="map-edge %s"><title>%s → %s: %s</title><path d="%s" fill="none" stroke="%s" stroke-width="%g" marker-end="url(#map-arrow-%s)"/><text x="%g" y="%g" text-anchor="middle" font-size="15" fill="%s" style="paint-order:stroke;stroke:white;stroke-width:5px">%s</text></g>',
      if(active) "selected" else "fixed",esc(e$from),esc(e$to),esc(e$rate_class),d,color,if(active) 2.5 else 1.5,if(active) "active" else "gray",tx,ty,color,label_markup))
  }
  for (node in p$nodes$id) {
    x <- xy[node,1]; y <- xy[node,2]
    g <- p$observations$glycoform[match(node,p$observations$node)]
    chosen <- !is.na(g) && g %in% measured
    h <- if(canonical) 114 else 60
    pieces <- c(pieces,sprintf('<g class="map-node %s"><title>%s</title><rect x="%g" y="%g" width="100" height="%g" rx="7" fill="%s" stroke="%s" stroke-width="%g"/>',
      if(chosen) "measured" else "unmeasured",esc(paste(node,if(!is.na(g)) paste("Measured as",g) else "Not observable")),x-50,y-h/2,h,
      if(chosen) "#edf5fc" else "#f7f8f8",if(chosen) "#377dac" else "#d9dfe2",if(chosen) 3 else 1))
    if (canonical) {
      pieces <- c(pieces,sprintf('<image href="snfg_cards/%s.png" x="%g" y="%g" width="86" height="%g"/>',esc(node),x-43,y-52,if (fit_mode) 73 else 104))
    } else {
      pieces <- c(pieces,sprintf('<text x="%g" y="%g" text-anchor="middle" font-size="13" fill="#35434b" textLength="%g" lengthAdjust="spacingAndGlyphs">%s</text>',x,y+if(fit_mode) -10 else 4,min(88,max(20,nchar(node)*8)),esc(node)))
    }
    if (fit_mode && chosen && !is.null(proportions) && g %in% names(proportions)) {
      fraction <- proportions[[g]]
      bar_y <- y + if (canonical) 24 else 0
      pooled <- sum(p$observations$glycoform == g) > 1
      if (pooled) {
        group_id <- match(g, unique(p$observations$glycoform))
        pieces <- c(pieces, sprintf('<text class="map-group-member" x="%g" y="%g" text-anchor="middle" font-size="12" fill="#315e7a">Group %d</text>', x, bar_y+16, group_id))
      } else pieces <- c(pieces, sprintf('<g class="map-proportion"><title>%s: %s%% of uploaded composition%s</title><rect x="%g" y="%g" width="78" height="9" rx="2" fill="#e1ebf2"/><rect x="%g" y="%g" width="%g" height="9" rx="2" fill="#4f8eb8"/><text x="%g" y="%g" text-anchor="middle" font-size="19" font-weight="700" fill="#315e7a">%s%%</text></g>',
        esc(g), format_number(100*fraction), if(pooled) " (shared measurement group)" else "",
        x-39,bar_y,x-39,bar_y,78*fraction,x,bar_y+28,format_number(100*fraction)))
    }
    pieces <- c(pieces,"</g>")
  }
  rendered <- paste0(c(pieces,"</svg>"),collapse="")
  if (zoom == 100) rendered <- sub("height:auto;display:block", if (fit_mode) "height:560px;display:block" else "height:520px;display:block", rendered, fixed = TRUE)
  groups <- unique(p$observations$glycoform)
  pooled <- groups[vapply(groups, function(g) sum(p$observations$glycoform == g) > 1, logical(1))]
  pooled <- intersect(pooled, intersect(measured, names(proportions)))
  if (fit_mode && length(pooled)) {
    rows <- vapply(pooled, function(g) sprintf('<li>Group %d: %s = %s%% (combined states: %s)</li>',
      match(g, groups), esc(g), trimws(format_number(100*proportions[[g]])),
      esc(paste(p$observations$node[p$observations$glycoform == g], collapse = ", "))), character(1))
    rendered <- paste0(rendered, '<div class="map-pooled-measurements"><strong>Shared measured classes</strong><ul>',
      paste0(rows, collapse=""), '</ul><p>Each percentage is the combined class total; individual state proportions are not measured separately.</p></div>')
  }
  HTML(rendered)
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
    tags$link(rel = "stylesheet", href = "design-map.css"),
    tags$script(HTML("document.addEventListener('change', function(e) { if(e.target.name === 'pathway_source' && e.target.value !== 'canonical') document.getElementById('pathway-editor').open = true; });")),
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



      /* Explanation: neutral = shared context; blue = supplied/selected input;
         burgundy = quantity sought; green = supported numerical result. */
      .how-page,.dc-demo { max-width:1240px; margin:28px auto; }
      .how-page h2,.dc-demo h2 { font-size:29px; margin:8px 0; letter-spacing:-.025em; }
      .how-page .eyebrow,.dc-demo .eyebrow { color:#626a70; font-size:11px; letter-spacing:.12em; font-weight:700; }
      .how-lead { color:#626a70; margin-bottom:24px; }
      .how-key { display:flex; gap:20px; flex-wrap:wrap; font-size:12px; color:#626a70; margin:0 0 24px; padding-bottom:12px; border-bottom:1px solid #e3e5e7; }
      .how-key span { display:flex; align-items:center; gap:7px; }
      .how-key i { width:9px; height:9px; display:inline-block; border-radius:2px; background:#a6abb0; }
      .how-key .input i { background:#356d98; } .how-key .question i { background:#b33135; }
      .how-key .result i { background:#376e57; }
      .how-flow { display:grid; grid-template-columns:1fr 28px 1.6fr 28px .9fr; align-items:center; }
      .how-card { align-self:stretch; border:1px solid #dadddf; border-radius:8px; padding:20px; background:#fff; }
      .how-card.unknown { background:#f5f6f7; border-color:#dadddf; }
      .how-card.output { background:#fff; border-color:#dadddf; }
      .how-card h3 { margin:8px 0 18px; font-size:20px; }
      .how-step { font-size:10px; letter-spacing:.08em; font-weight:700; color:#697178; }
      .how-card:first-child>.how-step { color:#356d98; }
      .how-card.unknown>.how-step { color:#b33135; }
      .how-card.output>.how-step { color:#376e57; }
      .how-join { text-align:center; font-size:25px; color:#a1a7ab; }
      .how-glycans { display:grid; grid-template-columns:repeat(3,1fr); gap:8px; align-items:end; }
      .how-glycans img { width:100%; height:auto; }
      .how-amount { height:85px; display:flex; align-items:end; justify-content:center; padding-bottom:5px; }
      .how-amount span { display:block; width:80%; text-align:center; padding-top:6px; font-weight:700; border-radius:3px 3px 0 0; background:#e2edf5!important; color:#295d84; }
      .how-note { font-size:12px; color:#697178; margin:14px 0 0; line-height:1.5; }
      .how-parameter { display:flex; align-items:center; justify-content:space-between; gap:10px; padding:15px 0; border-bottom:1px solid #e1e4e6; }
      .how-parameter strong { display:block; font-size:17px; color:#656e74; }
      .how-parameter small { color:#727a80; font-size:12px; }
      .how-number { font-size:36px; color:#b33135; font-weight:700; font-variant-numeric:tabular-nums; }
      .how-caption { color:#697178; font-size:12px; margin:14px 0 24px; }
      .how-check { display:grid; grid-template-columns:1fr 1fr; gap:25px; border-top:1px solid #dadddf; padding-top:19px; }
      .how-check h3 { font-size:17px; margin:0 0 6px; }
      .how-check p { font-size:13px; color:#697178; margin:0; }
      .how-check b { color:#b33135; }
      .dc-demo { margin-top:44px; padding-top:30px; border-top:2px solid #dadddf; }
      .dc-spread-band { margin:22px 0 28px; padding:20px; background:#fff; border:1px solid #dadddf; border-radius:8px; }
      .dc-spread-goal { font-size:18px; font-weight:650; margin-bottom:18px; }
      .dc-fixed-conditions { padding:12px 15px; background:#f3f4f5; color:#626a70; font-size:12px; border-radius:4px; }
      .dc-comparison { display:grid; grid-template-columns:repeat(2,minmax(0,1fr)); gap:20px; margin-top:24px; }
      .dc-plan { display:flex; flex-direction:column; border:1px solid #dadddf; border-radius:8px; overflow:hidden; }
      .dc-plan-heading { padding:18px 18px 12px; border-bottom:1px solid #e4e6e8; background:#fff; }
      .dc-plan-heading>.how-step { color:#356d98; }
      .dc-plan-heading h3 { font-size:20px; margin:8px 0 0; }
      .dc-change-note { display:inline-block; margin-top:10px; padding:4px 8px; border-radius:3px; background:#edf3f8; color:#356d98; font-size:12px; }
      .dc-aligned-outcome { padding:16px 12px 0; }
      .dc-conclusion { padding:14px 18px 18px; border-top:1px solid #e4e6e8; }
      .dc-conclusion>strong { display:block; font-size:16px; }
      .dc-audit-stats { display:flex; flex-wrap:wrap; gap:7px 16px; margin:12px 0; font-size:12px; font-variant-numeric:tabular-nums; color:#626a70; }
      .dc-conclusion details { margin-top:12px; }
      .dc-conclusion summary { color:#626a70; font-size:12px; }
      .dc-conclusion details p { font-size:13px; }
      .dc-error-aside { margin-top:24px; padding:15px 18px; border:1px solid #e1e4e6; border-radius:6px; background:#f7f8f8; }
      .dc-error-aside summary { font-size:13px; color:#626a70; }
      .dc-error-aside h3 { margin:14px 0; font-size:18px; }
      .dc-error-comparison { display:grid; grid-template-columns:1fr 1fr; gap:25px; }
      @media(max-width:850px) {
        .how-flow,.how-check,.dc-comparison,.dc-error-comparison { grid-template-columns:1fr; gap:14px; }
        .how-join { transform:rotate(90deg); }
        .how-glycans img { max-width:120px; display:block; margin:auto; }
        .how-key { gap:10px 16px; }
      }

      .how-workflow { background:#f6f7f8; border:1px solid #e0e3e5; border-radius:8px; padding:16px 20px; margin:0 0 28px; }
      .workflow-track { display:grid; grid-template-columns:1fr 36px 1fr 36px 1fr; gap:12px; align-items:start; margin-top:12px; }
      .workflow-stage { display:flex; flex-direction:column; gap:4px; padding:12px 16px; color:#626a70; border:1px solid #dadddf; border-radius:5px; }
      a.workflow-stage:hover { text-decoration:none; background:#fff; color:#343b40; }
      .workflow-stage small { font-size:10px; letter-spacing:.06em; }
      .workflow-stage strong { font-size:19px; color:#343b40; }
      .workflow-stage span { font-size:12px; }
      .workflow-stage.destination { background:#fff; border-color:#b99b9d; }
      .workflow-stage.destination em { color:#963439; font-size:12px; font-style:normal; font-weight:650; border-top:1px solid #eadfe0; padding-top:8px; margin-top:5px; }
      .workflow-arrow { text-align:center; padding-top:25px; color:#969fa5; font-size:23px; }
      .how-step-back { display:flex; gap:18px; align-items:center; margin-top:25px; padding:20px; background:#f3f5f7; border-left:3px solid #356d98; }
      .how-step-back>span { font-size:35px; color:#356d98; }
      .how-step-back strong { font-size:19px; }
      .how-step-back p { font-size:14px; color:#626a70; margin:5px 0 0; }
      .how-step-back>a { margin-left:auto; white-space:nowrap; font-size:14px; color:#356d98; }
      #how-fit-explanation,#how-design-demo { scroll-margin-top:20px; }
      @media(max-width:700px) { .workflow-track { grid-template-columns:1fr; gap:6px; } .workflow-arrow { padding:0; transform:rotate(90deg); } .how-step-back { flex-wrap:wrap; } }
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
      .pathway-example { margin:0 0 12px; padding:8px; border:1px solid var(--line); background:var(--soft); }
      .pathway-example figcaption { font-size:13px; font-weight:600; margin-bottom:5px; }
      .reaction-example { width:100%; border-collapse:collapse; font-size:13px; }
      .reaction-example th,.reaction-example td { padding:4px 5px; text-align:left; border-bottom:1px solid var(--line); }
      .reaction-example th { font-weight:600; vertical-align:top; }
      .reaction-example td { font-family:ui-monospace,SFMono-Regular,Menlo,monospace; }
      .pathway-example .small-note { margin:6px 0 0; font-size:12px; }
      .builder-block { margin:12px 0; }
      .builder-block h4 { font-size:14px; margin:14px 0 7px; }
      .builder-step { display:inline-block; margin-right:7px; padding:2px 5px; background:#ececeb; color:var(--ink); font-size:12px; }
      .reaction-composer { padding:12px; border:1px solid #b8c8cf; border-radius:5px; background:#f5f9fa; }
      .reaction-composer h4 { margin:0 0 10px; }
      .reaction-mode .shiny-options-group { display:flex; gap:0; border:1px solid #b8c8cf; border-radius:4px; overflow:hidden; background:#fff; }
      .reaction-mode .radio { flex:1; margin:0; }
      .reaction-mode .radio label { display:block; padding:8px 6px; text-align:center; font-size:13px; cursor:pointer; }
      .reaction-mode .radio input { position:absolute; opacity:0; }
      .reaction-mode .radio:has(input:checked) { background:#284c5b; color:#fff; }
      .reaction-mode .radio:focus-within { outline:2px solid #111; outline-offset:-2px; }
      .reaction-fields { display:grid; grid-template-columns:minmax(0,1fr) 18px minmax(0,1fr); align-items:end; gap:4px; }
      .reaction-fields .form-group { min-width:0; margin-bottom:8px; }
      .reaction-fields .control-label,.reaction-parameter .control-label { font-size:12px; }
      .reaction-fields input,.reaction-parameter input { font-family:ui-monospace,SFMono-Regular,Menlo,monospace; font-size:14px; }
      .reaction-fields input[readonly],.reaction-parameter input[readonly] { background:#e9f0f3; border-color:#c3d0d6; }
      .reaction-fields .conversion-arrow { text-align:center; font-size:21px; padding-bottom:11px; color:#284c5b; }
      .reaction-parameter .form-group { margin-bottom:8px; }
      .reaction-composer details { margin-top:4px; padding-top:5px; }
      .reaction-composer summary { font-size:12px; color:var(--muted); }
      .reaction-composer .add-row { width:100%; margin-top:10px; background:#284c5b; border-color:#284c5b; }
      .reaction-flow { text-align:center; color:#718b96; font-size:25px; line-height:1; padding:5px; }
      .pathway-table-header { display:flex; align-items:center; justify-content:space-between; gap:8px; }
      .pathway-table-header h4 { margin:0 0 8px; }
      .pathway-count { font-size:12px; color:var(--muted); white-space:nowrap; }
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
      .shared-model { margin:16px 0; padding:12px 16px; border:1px solid #dadddf; background:#f7f8f8; border-radius:6px; }
      .shared-model-body { max-width:650px; padding-top:18px; }
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
      tags$p("Design glycoform measurements for interpretable reaction-rate estimates")),
  tabsetPanel(
    id = "app_mode", selected = "how",
    tabPanel("Overview", value = "how",
      HTML('<section class="how-page" aria-label="Inputs, unknown reaction rates, and estimated values">
        <nav class="how-workflow" aria-label="Research workflow and explanation guide">
          <div class="how-step">THE RESEARCH WORKFLOW</div>
          <div class="workflow-track">
            <a class="workflow-stage" href="#how-design-demo"><small>BEFORE MEASUREMENT</small><strong>Design check</strong><span>Choose informative measurements</span></a>
            <span class="workflow-arrow" aria-hidden="true">→</span>
            <div class="workflow-stage measurement"><small>EXPERIMENT</small><strong>Measure</strong><span>Obtain glycan proportions</span></div>
            <span class="workflow-arrow" aria-hidden="true">→</span>
            <a class="workflow-stage destination" href="#how-fit-explanation"><small>AFTER MEASUREMENT</small><strong>Fit</strong><span>Estimate reaction rates</span><em>↓ First, understand this goal</em></a>
          </div>
        </nav>
        <div class="how-key" aria-label="Visual key"><span><i></i>Shared model &amp; assumptions</span><span class="input"><i></i>Input</span><span class="question"><i></i>What you want to know</span><span class="result"><i></i>Calculated precision</span></div>
        <div class="eyebrow" id="how-fit-explanation">THE GOAL · FIT</div>
        <h2>Put a number on each reaction’s rate.</h2>
        <p class="how-lead">The glycan proportions are measured. The reaction rates are unknown.</p>
        <div class="how-flow">
          <section class="how-card">
            <div class="how-step">1 · WHAT YOU PROVIDE</div>
            <h3>Measured proportions</h3>
            <div class="how-glycans">
              <div><div class="how-amount"><span style="height:80px;background:#9fc6db">40%</span></div><img src="snfg_cards/M5.png" alt="M5 glycan structure"></div>
              <div><div class="how-amount"><span style="height:40px;background:#e3c482">20%</span></div><img src="snfg_cards/M5Gn.png" alt="M5Gn glycan structure"></div>
              <div><div class="how-amount"><span style="height:80px;background:#a8cdb7">40%</span></div><img src="snfg_cards/M4Gn.png" alt="M4Gn glycan structure"></div>
            </div>

            <p class="how-note">All three glycans are measured here.<br>These proportions are the observed inputs.</p>
          </section>
          <div class="how-join" aria-hidden="true">→</div>
          <section class="how-card unknown">
            <div class="how-step">2 · THE NUMBERS WE WANT TO FIND</div>
            <h3>Known pathway. Unknown rates.</h3>
            <svg viewBox="0 0 530 250" style="width:100%;height:auto;margin-top:18px" role="img" aria-label="M5 converts to M5Gn through MGAT1, then to M4Gn through MAN2. The numerical rate of each reaction is unknown.">
              <defs><marker id="how-rate-arrow" markerWidth="7" markerHeight="7" refX="6" refY="3.5" orient="auto"><path d="M0 0L7 3.5L0 7" fill="#a0a7ad"/></marker></defs>
              <image href="snfg_cards/M5.png" x="0" y="30" width="110" height="110"/>
              <image href="snfg_cards/M5Gn.png" x="210" y="30" width="110" height="110"/>
              <image href="snfg_cards/M4Gn.png" x="420" y="30" width="110" height="110"/>
              <path d="M115 85H201" stroke="#a0a7ad" stroke-width="3" marker-end="url(#how-rate-arrow)"/><path d="M325 85H411" stroke="#a0a7ad" stroke-width="3" marker-end="url(#how-rate-arrow)"/>
              <g text-anchor="middle" fill="#687078" font-size="18" font-weight="650"><text x="158" y="58">MGAT1</text><text x="368" y="58">MAN2</text></g>
              <g text-anchor="middle" fill="#697178" font-size="15"><text x="158" y="124">Rate =</text><text x="368" y="124">Rate =</text></g>
              <rect x="128" y="137" width="60" height="65" rx="9" fill="#fce8e8"/><rect x="338" y="137" width="60" height="65" rx="9" fill="#fce8e8"/>
              <g text-anchor="middle" fill="#b33135" font-size="53" font-weight="750"><text x="158" y="189">?</text><text x="368" y="189">?</text></g>
              <text x="265" y="235" text-anchor="middle" font-size="15" fill="#697178">Find these two numerical values.</text>
            </svg>
            <p class="how-note">You supply the connections and reaction names.<br>Fit estimates the unknown rates from the measured proportions.</p>
          </section>
          <div class="how-join" aria-hidden="true">→</div>
          <section class="how-card output how-result">
            <div class="how-step">3 · FIT CALCULATES ESTIMATES</div>
            <h3>Estimated rate values</h3>
            <div class="how-parameter"><div><strong>MGAT1</strong><small>M5 → M5Gn</small></div><span class="how-number">1.5</span></div>
            <div class="how-parameter"><div><strong>MAN2</strong><small>M5Gn → M4Gn</small></div><span class="how-number">2.0</span></div>
            <p class="how-note">Effective rate parameters relative to a fixed secretion rate of 1. These are not enzyme amounts or directly measured enzyme activities.</p>
          </section>
        </div>
        <p class="how-caption">Illustrative numbers in a simplified three-glycan model: entry at M5, equal secretion, no further reactions. Not a fitted biological result.</p>
        <div class="how-step-back">
          <span aria-hidden="true">↶</span>
          <div><strong>Now step back to before measurement.</strong><p>Which glycans would you need to measure to make those rate estimates precise?</p></div>
          <a href="#how-design-demo">Design check ↓</a>
        </div>
      </section>'),
      div(id = "how-design-demo", class = "dc-demo",
        div(class = "eyebrow", "PLAN FOR THAT GOAL · DESIGN CHECK"),
        tags$h2("Which glycans should I measure?"),
        div(class = "dc-spread-band",
          tags$p(class = "dc-spread-goal", "You want to estimate MGAT1 and MAN2 reaction rates separately. You have not measured the glycans yet."),
          HTML('<svg viewBox="0 0 1120 365" style="width:100%;height:auto" role="img" aria-label="Choose which glycan proportions to measure: M5, M5Gn, or M4Gn. Red question marks indicate measurement choices. The goal, shown beneath each reaction, is a precise rate estimate with a narrow uncertainty range.">
<defs><marker id="choose-arrow" markerWidth="7" markerHeight="7" refX="6" refY="3.5" orient="auto"><path d="M0 0L7 3.5L0 7" fill="#a0a7ad"/></marker></defs>
<g font-family="Arial,sans-serif" text-anchor="middle">
<g fill="#b33135" font-size="24" font-weight="700"><text x="100" y="29">Measure ?</text><text x="560" y="29">Measure ?</text><text x="1020" y="29">Measure ?</text></g>
<image href="snfg_cards/M5.png" x="20" y="46" width="160" height="145"/>
<image href="snfg_cards/M5Gn.png" x="480" y="46" width="160" height="145"/>
<image href="snfg_cards/M4Gn.png" x="940" y="46" width="160" height="145"/>
<path d="M195 111H460" stroke="#a0a7ad" stroke-width="3" marker-end="url(#choose-arrow)"/>
<path d="M655 111H920" stroke="#a0a7ad" stroke-width="3" marker-end="url(#choose-arrow)"/>
<g fill="#687078" font-size="23" font-weight="700"><text x="330" y="89">MGAT1</text><text x="790" y="89">MAN2</text></g>
<g fill="#737b81" font-size="15"><text x="330" y="143">Reaction rate</text><text x="790" y="143">Reaction rate</text></g>
<path d="M330 154V190" stroke="#a0a7ad" stroke-width="2" marker-end="url(#choose-arrow)"/>
<path d="M790 154V190" stroke="#a0a7ad" stroke-width="2" marker-end="url(#choose-arrow)"/>
<rect x="193" y="202" width="274" height="142" rx="10" fill="#ffffff" stroke="#9eafa5" stroke-width="1.5" stroke-dasharray="5 4"/>
<rect x="653" y="202" width="274" height="142" rx="10" fill="#ffffff" stroke="#9eafa5" stroke-width="1.5" stroke-dasharray="5 4"/>
<g fill="#315e47" font-size="13" font-weight="700"><text x="330" y="225">GOAL · PRECISE RATE ESTIMATE</text><text x="790" y="225">GOAL · PRECISE RATE ESTIMATE</text></g>
<g stroke="#c2d2c8" stroke-width="2"><path d="M224 275H436 M684 275H896"/></g>
<g stroke="#376e57" stroke-width="3"><path d="M306 275H354 M306 265V285 M354 265V285 M766 275H814 M766 265V285 M814 265V285"/></g>
<g fill="#376e57"><circle cx="330" cy="275" r="7"/><circle cx="790" cy="275" r="7"/></g>
<g fill="#315e47" font-size="16" font-weight="700"><text x="330" y="250">MGAT1</text><text x="790" y="250">MAN2</text><text x="330" y="308">↔ Small estimation error</text><text x="790" y="308">↔ Small estimation error</text></g>
<g fill="#737b81" font-size="12"><text x="330" y="329">Target, not a calculated result</text><text x="790" y="329">Target, not a calculated result</text></g>
</g></svg>'),
          div(class = "dc-fixed-conditions", tags$strong("FIXED IN BOTH EXAMPLES · "),
            "MGAT1 rate = 1.5 · MAN2 rate = 2.0 · Measurement CV = 5%. Reference rates are hypothetical, not measured values.")
        ),
        tags$h3("Compare two choices of what to measure"),
        tags$p("For each choice, Design check calculates whether both rates can be estimated separately and with sufficient precision."),
        uiOutput("design_concept_result"),
        tags$p(class = "how-caption", "Simplified three-glycan model, not your current design. Local calculation at assumed MGAT1 = 1.5 and MAN2 = 2.0; secretion = 1, entry at M5. Target worst-direction log-rate SD = 0.25. These are planning assumptions, not fitted results.")
      )
    ),

    tabPanel("Design check", value = "design",
      tags$section(class = "design-plan",
        div(class = "plan-section-title", tags$h2("Inputs → pathway → assessment"),
            tags$a(href = "#design-results", "View results ↓")),
      div(class = "design-layout",
        div(class = "input-pane",
    div(class = "shared-model pathway-scope",
      tags$h2("1  Define the pathway"),
      radioButtons("pathway_source", NULL,
        choices = c("Manuscript example" = "canonical", "Build a pathway" = "builder", "Upload CSV" = "upload"),
        selected = "canonical", inline = FALSE),
      uiOutput("pathway_scope"),
      tags$details(id = "pathway-editor",
      tags$summary("Edit pathway"),
      div(class = "shared-model-body",
          tags$figure(class = "pathway-example",
            tags$figcaption("Example input · one reaction per row"),
            tags$table(class = "reaction-example",
              tags$thead(tags$tr(
                tags$th(scope = "col", "From glycan"),
                tags$th(scope = "col", "To glycan"),
                tags$th(scope = "col", "Effective rate class")
              )),
              tags$tbody(lapply(c(1L, 2L, 5L), function(i) {
                tags$tr(tags$td(registered_reactions$From[i]),
                        tags$td(registered_reactions$To[i]),
                        tags$td(registered_reactions[["Rate class"]][i]))
              }))
            ),
            tags$p(class = "small-note",
                   "Rate class = parameter name, not a speed value.")
          ),

          conditionalPanel(
            "input.pathway_source == 'builder'",
            div(class = "builder-block",
              div(class = "reaction-composer",
                tags$h4(tags$span(class = "builder-step", "1"), "Choose or write a reaction"),
                div(class = "reaction-mode",
                  radioButtons(
                    "builder_input_mode", NULL,
                    choices = c("Registered" = "guided", "Custom" = "custom"),
                    selected = "guided", inline = FALSE
                  )
                ),
                conditionalPanel(
                  "input.builder_input_mode == 'guided'",
                  selectInput("registered_reaction", "Registered N-glycan reactions",
                              choices = registered_reaction_choices)
                ),
                conditionalPanel(
                  "input.builder_input_mode == 'custom'",
                  tags$p(class = "small-note", "Custom definitions")
                ),
                tags$h4(tags$span(class = "builder-step", "2"), "Reaction to add"),
                conditionalPanel(
                  "input.builder_input_mode == 'guided'",
                  uiOutput("registered_fields")
                ),
                conditionalPanel(
                  "input.builder_input_mode == 'custom'",
                  div(class = "reaction-fields",
                    textInput("custom_from", "From glycan", "", placeholder = "e.g. M5"),
                    tags$span(class = "conversion-arrow", `aria-hidden` = "true", "→"),
                    textInput("custom_to", "To glycan", "", placeholder = "e.g. M5Gn")
                  ),
                  div(class = "reaction-parameter",
                    textInput("custom_rate_class", "Rate class · parameter name", "", placeholder = "e.g. MGAT1")
                  ),
                  tags$details(
                    tags$summary("Weight"),
                    numericInput("custom_weight", "Rate multiplier", 1, min = 0.001, step = 0.1)
                  )
                ),
                actionButton("add_reaction", "Add reaction to pathway", icon = icon("plus"),
                             class = "btn-primary add-row")
              ),
              div(class = "reaction-flow", `aria-hidden` = "true", "↓"),
              div(class = "pathway-table-header",
                tags$h4(tags$span(class = "builder-step", "3"), "Your pathway"),
                uiOutput("pathway_row_count", inline = TRUE)
              ),
              DT::DTOutput("reaction_editor"),
              div(class = "builder-actions",
                  actionButton("delete_reaction", "Delete selected", class = "btn-default"),
                  actionButton("clear_reactions", "Start blank", class = "btn-default")),
              tags$details(
                tags$summary("Optional: group indistinguishable glycans"),
                tags$p(class = "small-note",
                       "Only change this if your assay reports several glycans as one combined signal. Give those glycans the same measurement name; their amounts will be summed. Double-click a name to edit it."),
                DT::DTOutput("state_editor")
              ),
              tags$details(
                tags$summary("Advanced: inflow and secretion"),
                tags$p(class = "small-note",
                       "Entry sets the inflow into each glycan; Secretion sets its exit rate constant. Keep the example values unless you need different boundary assumptions."),
                DT::DTOutput("boundary_editor")
              ),
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
          uiOutput("model_summary")
      )
    ),
    ),

          div(class = "section-rule",
            tags$h3("2  Reaction rates to estimate"),
            uiOutput("rate_selection_note"),
            uiOutput("rate_selector"),
            div(class = "action-links", actionLink("select_all_rates", "Select all"), actionLink("clear_rates", "Clear"))
          ),
          div(class = "section-rule",
            tags$h3("3  Glycans to measure"),
            uiOutput("glycan_selection_note"),
            uiOutput("panel_selector"),
            div(class = "action-links", actionLink("select_all_glycoforms", "Select all"), actionLink("clear_glycoforms", "Clear")),
            uiOutput("selected_glycoforms")
          ),
          div(class = "section-rule",
            tags$h3("Measurement precision"),
            fluidRow(
              column(6, numericInput("assay_cv", "Single-measurement CV (%)", 5, min = 0.01, max = 100, step = 0.1)),
              column(6, numericInput("replicates", "Independent replicates", 1, min = 1, max = 100, step = 1))
            ),
            uiOutput("effective_cv_text"),
            tags$details(
              tags$summary("Rate-precision target"),
              numericInput("target_sd", "Target worst-direction log-rate SD", 0.25, min = 0.01, max = 2, step = 0.05),
              tags$p(class = "small-note", "The default 0.25 is an illustrative planning target, not a validated biological cutoff. It corresponds to about 28% multiplicative uncertainty for one standard deviation; choose a target appropriate to your question.")
            )
          ),
          tags$p(class = "conditional-note",
                 "Results update automatically and are local to the supplied graph, reaction grouping, operating point, observation mapping, secretion model and error assumptions.")
        ),
        div(class = "result-pane",
          div(class = "map-heading", tags$h2("Your pathway"),
            tags$p("Gray: pathway context · Burgundy arrows: rates to estimate · Blue boxes: glycans to measure")),
          div(class = "map-controls", sliderInput("pathway_zoom", "Diagram size (%)", min = 100, max = 250, value = 100, step = 25, width = "230px")),
          div(class = "pathway-map-viewport", uiOutput("design_pathway_map")),
          tags$p(class = "small-note", "Unselected reactions remain in the model with their rates fixed at the reference value (1). Unselected classes are excluded from the measured panel. Inflow and secretion are omitted from this diagram. Glycan structures depict model assumptions, not structures inferred from composition data."),
      div(class = "plan-to-result", `aria-hidden` = "true", "↓"),
      tags$section(class = "design-results", id = "design-results",
        tags$div(class = "plan-section-title", tags$h2("Results for this measurement plan"),
          tags$span("Updates with your selections")),
          uiOutput("result_header"),
          tags$details(tags$summary("Details & next steps"),
          uiOutput("result_context"), uiOutput("decision_table"),
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
        )
      )
      )
    ),
    tabPanel("Robustness & fit", value = "advanced",
      div(class = "advanced-layout fit-layout",
        div(class = "advanced-controls",
          tags$h2("Additional checks"),
          tags$h3("Operating-point robustness"),
          tags$p(class = "small-note", "Uses the pathway, target rates, measured panel and precision settings from Design check."),
          numericInput("n_points", "Operating points", 100, min = 10, max = 2000),
          numericInput("log_rate_sd", "Log-rate sampling SD", 0.7, min = 0.01),
          numericInput("robust_seed", "Random seed", 1, min = 1),
          actionButton("run_robustness", "Run robustness check", class = "btn-primary btn-block"),
          tags$hr(),
          tags$h3("Enter measured glycan proportions"),
          tags$p(class = "small-note", "One composition from one glycosylation site. One row per measured glycan."),
          uiOutput("fit_proportion_input"),
          actionButton("load_fit_example", "Load example CSV", class = "btn-block"),
          fileInput("observed_file", "Choose your CSV", accept = ".csv"),
          actionButton("run_fit", "Estimate reaction rates →", class = "btn-primary btn-block"),
          tags$details(tags$summary("CSV format & settings"),
            tableOutput("fit_csv_example"),
            tags$p(class = "small-note", "CSV uses decimals: 0.40 means 40%."),
            tags$p(class = "small-note", "Columns: glycoform, proportion. Supply at least two distinct glycans from this pathway. Values are normalized over uploaded glycans. Example values are simulated at the model’s reference rates, not experimental observations."),
            numericInput("fit_starts", "Optimisation starts", 8, min = 1, max = 100))
        ),
        div(class = "fit-main",
          div(class = "fit-map-heading", tags$h3("Estimate rates from these proportions"),
            uiOutput("fit_pathway_context"),
            actionLink("edit_pathway_from_fit", "Change pathway")),
          uiOutput("fit_diagram_key"),
          uiOutput("fit_state_message"),
          div(class = "fit-pathway", uiOutput("fit_pathway_diagram")),
          tags$p(class = "small-note", "Effective rate parameters under the pathway’s fixed secretion and inflow assumptions; not enzyme amounts or directly measured enzyme activities. B4GALT, ST6GAL and FUT8 denote the model classes GALT, SIAT and FUT; MGAT4 denotes late branching. These labels and glycan icons do not establish enzyme identity, linkage or fucose position from composition alone."),
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
  observeEvent(input$edit_pathway_from_fit, updateTabsetPanel(session, "app_mode", selected = "design"))

  output$design_concept_result <- renderUI({
    demo <- glyco_pathway(
      data.frame(id = c("M5", "M5Gn", "M4Gn")),
      data.frame(from = c("M5", "M5Gn"), to = c("M5Gn", "M4Gn"),
                 rate_class = c("MGAT1", "MAN2")), entry = "M5", secretion = 1)
    full_audit <- check_panel(demo, demo$nodes$id, log_rates = log(c(MGAT1 = 1.5, MAN2 = 2)), assay_cv = 0.05, target_sd = 0.25)
    noisy_audit <- check_panel(demo, demo$nodes$id, log_rates = log(c(MGAT1 = 1.5, MAN2 = 2)), assay_cv = 0.30, target_sd = 0.25)
    tagList(
      div(class = "dc-comparison dc-subsets", lapply(c("two", "three"), function(plan) {
        panel <- if (plan == "two") c("M5", "M4Gn") else demo$nodes$id
        audit <- check_panel(demo, panel, log_rates = log(c(MGAT1 = 1.5, MAN2 = 2)), assay_cv = 0.05, target_sd = 0.25)
        rate_sd <- c(MGAT1 = Inf, MAN2 = Inf)
        if (audit$full_rank) {
          y <- audit$composition
          noise_map <- diag(y) %*% (diag(length(y)) - outer(rep(1, length(y)), y))
          inverse <- solve(crossprod(audit$jacobian), t(audit$jacobian))
          covariance <- 0.05^2 * inverse %*% noise_map %*% t(noise_map) %*% t(inverse)
          rate_sd <- sqrt(diag(covariance))
        }
        outcome_boxes <- paste(vapply(seq_along(c("MGAT1", "MAN2")), function(i) {
          rate <- c("MGAT1", "MAN2")[[i]]
          cx <- c(160, 380)[[i]]
          col <- if (audit$full_rank) "#376e57" else "#b33135"
          background <- if (audit$full_rank) "#f0f6f2" else "#fcf4f3"
          interval <- if (audit$full_rank) {
            half_width <- 72 * rate_sd[[rate]] / 0.25
            sprintf('<line x1="%s" y1="268" x2="%s" y2="268" stroke="#bac6be"/><path d="M%s 264V273 M%s 257V273 M%s 264V273" stroke="#63756b"/><line x1="%s" y1="268" x2="%s" y2="268" stroke="%s" stroke-width="3"/><path d="M%s 260V276 M%s 260V276" stroke="%s" stroke-width="2"/><circle cx="%s" cy="268" r="4" fill="%s"/><g font-size="11" fill="#53645a"><text x="%s" y="291">−0.25</text><text x="%s" y="291">0</text><text x="%s" y="291">+0.25</text><text x="%s" y="308">log(rate / reference)</text></g>',
              cx-72,cx+72,cx-72,cx,cx+72,cx-half_width,cx+half_width,col,cx-half_width,cx+half_width,col,cx,col,cx-72,cx,cx+72,cx)
          } else {
            sprintf('<path d="M%s 268H%s M%s 259L%s 268L%s 277 M%s 259L%s 268L%s 277" fill="none" stroke="%s" stroke-width="2" stroke-dasharray="4 3"/><text x="%s" y="305" font-size="12" fill="%s">No finite marginal SD</text>',cx-72,cx+72,cx-62,cx-72,cx-62,cx+62,cx+72,cx+62,col,cx,col)
          }
          sprintf('<rect x="%s" y="197" width="190" height="173" rx="5" fill="%s" stroke="%s"/><text x="%s" y="219" fill="%s" font-size="10" font-weight="700">CALCULATED PRECISION</text><text x="%s" y="241" font-size="18" fill="%s" font-weight="700">%s</text>%s<text x="%s" y="335" font-size="16" fill="%s" font-weight="700">%s</text><text x="%s" y="355" font-size="11" fill="%s">%s</text>',
            cx-95,background,col,cx,col,cx,col,rate,interval,cx,col,
            if (audit$full_rank) sprintf("SD = %.3f", rate_sd[[rate]]) else "Not identifiable",cx,col,
            if (audit$full_rank) "Error bar: ±1 SD, not a 95% CI" else "Structural information limit")

        }, character(1)), collapse = "")
        outcome_svg <- sprintf('<svg viewBox="0 0 540 385" style="width:100%%;height:auto" role="img" aria-label="Selected glycans connected to predicted estimation precision for MGAT1 and MAN2"><g font-family="Arial,sans-serif" text-anchor="middle"><g font-size="13" font-weight="700" fill="#737b81"><text x="50" y="25">✓ Measure</text><text x="490" y="25">✓ Measure</text></g><text x="270" y="25" font-size="13" font-weight="700" fill="%s">%s</text><image href="snfg_cards/M5.png" x="10" y="44" width="80" height="80"/><rect x="224" y="38" width="92" height="92" rx="6" fill="%s" stroke="%s" stroke-width="2" stroke-dasharray="%s"/><image href="snfg_cards/M5Gn.png" x="230" y="44" width="80" height="80" opacity="%s"/><image href="snfg_cards/M4Gn.png" x="450" y="44" width="80" height="80"/><g stroke="#a0a7ad" stroke-width="2" fill="none"><path d="M98 81H219 M211 76L219 81L211 86 M318 81H439 M431 76L439 81L431 86 M160 121V183 M155 175L160 183L165 175 M380 121V183 M375 175L380 183L385 175"/></g><g fill="#687078" font-size="16" font-weight="700"><text x="160" y="67">MGAT1</text><text x="380" y="67">MAN2</text></g><g fill="#737b81" font-size="12"><text x="160" y="109">Reaction rate</text><text x="380" y="109">Reaction rate</text></g>%s</g></svg>',
          if (plan == "two") "#737b81" else "#356d98", if (plan == "two") "Not measured" else "✓ Measure", if (plan == "two") "#f4f5f6" else "#edf3f8", if (plan == "two") "#b5bbc0" else "#356d98", if (plan == "two") "4 3" else "none", if (plan == "two") "0.25" else "1", outcome_boxes)
        div(class = "dc-plan",
          div(class = "dc-plan-heading",
            div(class = "how-step", if (plan == "two") "EXAMPLE A · SELECT 2 OF 3" else "EXAMPLE B · SELECT ALL 3"),
            tags$h3(if (plan == "two") "Measure M5 and M4Gn" else "Also measure M5Gn"),
            div(class = "dc-change-note", if (plan == "two") "M5Gn excluded" else "Only change: include M5Gn")
          ),
          div(class = "dc-aligned-outcome", HTML(outcome_svg)),
          div(class = "dc-conclusion",
            tags$strong(if (audit$full_rank) "Adding M5Gn makes precise estimation possible here." else "M5 + M4Gn alone cannot resolve both rates."),
            tags$p(class = "how-note", "Local assessment under the shared assumptions above."),
            div(class = "dc-audit-stats",
              tags$span(paste0("Local rank: ", audit$rank, " / 2")),
              tags$span(if (audit$full_rank) paste0("Worst-direction SD: ", round(audit$achieved_sd, 3), " ≤ 0.25 target") else "Precision target: not assessable")
            ),
            tags$details(tags$summary("Statistical method and assumptions"),
              tags$p("Local first-order propagation of independent 5% relative intensity errors through compositional normalisation and the inverse log-rate Jacobian. Error bars show marginal standard deviations; the pass criterion uses the square root of the largest covariance eigenvalue (worst-direction SD). Neither the bars nor full local rank establish global identifiability or model validity."),
              if (audit$full_rank) tagList(
                tags$p(paste0("Worst-direction log-rate SD: ", round(audit$achieved_sd, 3), "; target ≤ 0.25.")),
                lapply(c("MGAT1", "MAN2"), function(rate) {
                  reference <- c(MGAT1 = 1.5, MAN2 = 2)[[rate]]
                  tags$p(sprintf("%s: predicted ±1 SD range %.2f–%.2f around reference %.1f (log-rate scale, transformed back).", rate, reference * exp(-rate_sd[[rate]]), reference * exp(rate_sd[[rate]]), reference))
                })
              ) else tags$p("Two normalised proportions give one independent quantity for two rate parameters. Different rate combinations can produce the same observed proportions.")
            )
          )
        )
      })),
      tags$details(class = "dc-error-aside",
        tags$summary("Next: what if the same three glycans are measured less precisely?"),
        tags$h3("Keep all three glycans selected. Change only the assumed error."),
        div(class = "dc-error-comparison",
          div(tags$strong("5% measurement CV"), tags$p(paste0("Expected log-rate SD: ", round(full_audit$achieved_sd, 2), " · meets the 0.25 target"))),
          div(tags$strong("30% measurement CV"), tags$p(paste0("Expected log-rate SD: ", round(noisy_audit$achieved_sd, 2), " · exceeds the 0.25 target")))
        )
      )
    )
  })

  reaction_rows <- reactiveVal(default_builder_reactions())
  state_rows <- reactiveVal(default_builder_states())

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
      state_rows()[, c("State", "Measured class"), drop = FALSE],
      colnames = c("Glycan", "Measurement name"), rownames = FALSE, selection = "none",
      editable = list(target = "cell", disable = list(columns = c(0))),
      options = list(dom = "t", paging = FALSE, ordering = FALSE, autoWidth = FALSE,
                     scrollY = "180px", scrollCollapse = TRUE,
                     columnDefs = list(list(width = "50%", targets = c(0, 1))))
    )
  }, server = FALSE)

  observeEvent(input$state_editor_cell_edit, {
    updated <- state_rows()
    edited <- DT::editData(updated[, c("State", "Measured class"), drop = FALSE],
                           input$state_editor_cell_edit, rownames = FALSE)
    updated[["Measured class"]] <- trimws(as.character(edited[["Measured class"]]))
    state_rows(updated)
  })

  output$boundary_editor <- DT::renderDT({
    DT::datatable(
      state_rows()[, c("State", "Entry", "Secretion"), drop = FALSE],
      colnames = c("Glycan", "Entry", "Secretion"), rownames = FALSE, selection = "none",
      editable = list(target = "cell", disable = list(columns = c(0))),
      options = list(dom = "t", paging = FALSE, ordering = FALSE, autoWidth = FALSE,
                     scrollY = "180px", scrollCollapse = TRUE)
    )
  }, server = FALSE)

  observeEvent(input$boundary_editor_cell_edit, {
    updated <- state_rows()
    edited <- DT::editData(updated[, c("State", "Entry", "Secretion"), drop = FALSE],
                           input$boundary_editor_cell_edit, rownames = FALSE)
    updated$Entry <- suppressWarnings(as.numeric(edited$Entry))
    updated$Secretion <- suppressWarnings(as.numeric(edited$Secretion))
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

  pending_reaction <- reactive({
    if (identical(input$builder_input_mode %||% "guided", "guided")) {
      index <- suppressWarnings(as.integer(input$registered_reaction))
      if (!length(index) || is.na(index) || index < 1L || index > nrow(registered_reactions)) return(NULL)
      return(registered_reactions[index, , drop = FALSE])
    }
    data.frame(
      From = trimws(input$custom_from %||% ""),
      To = trimws(input$custom_to %||% ""),
      `Rate class` = trimws(input$custom_rate_class %||% ""),
      Weight = suppressWarnings(as.numeric(input$custom_weight %||% NA_real_)),
      check.names = FALSE, stringsAsFactors = FALSE
    )
  })

  output$registered_fields <- renderUI({
    index <- suppressWarnings(as.integer(input$registered_reaction))
    req(length(index) == 1L, !is.na(index), index >= 1L, index <= nrow(registered_reactions))
    candidate <- registered_reactions[index, , drop = FALSE]
    fixed_field <- function(id, label, value) {
      div(class = "form-group",
          tags$label(class = "control-label", `for` = id, label),
          tags$input(id = id, type = "text", class = "form-control", readonly = "readonly", value = value))
    }
    tagList(
      div(class = "reaction-fields",
          fixed_field("registered_from", "From glycan", candidate$From[[1]]),
          tags$span(class = "conversion-arrow", `aria-hidden` = "true", "→"),
          fixed_field("registered_to", "To glycan", candidate$To[[1]])),
      div(class = "reaction-parameter",
          fixed_field("registered_rate", "Rate class · parameter name", candidate[["Rate class"]][[1]])),
      tags$details(tags$summary("Weight"),
                   fixed_field("registered_weight", "Rate multiplier", candidate$Weight[[1]]))
    )
  })

  output$pathway_row_count <- renderUI({
    n <- nrow(reaction_rows())
    tags$span(class = "pathway-count", paste(n, if (n == 1L) "reaction" else "reactions"))
  })

  observeEvent(input$add_reaction, {
    candidate <- pending_reaction()
    if (is.null(candidate)) return()
    from <- candidate$From[[1]]
    to <- candidate$To[[1]]
    rate <- candidate[["Rate class"]][[1]]
    weight <- candidate$Weight[[1]]
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
    current <- reaction_rows()
    duplicate <- nrow(current) > 0 && any(
      current$From == from & current$To == to & current[["Rate class"]] == rate
    )
    if (duplicate) {
      showNotification("That reaction is already in your pathway. Select a different conversion or use Start blank to begin again.", type = "message")
      return()
    }
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
    identical(p$nodes$id, reference$nodes$id) && identical(p$edges, reference$edges) &&
      identical(p$observations, reference$observations)
  })

  output$model_summary <- renderUI({
    p <- pathway()
    div(class = "model-summary",
        div(tags$strong(nrow(p$nodes)), tags$span("network states")),
        div(tags$strong(length(p$observed_classes)), tags$span("observable classes")),
        div(tags$strong(length(p$rate_classes)), tags$span("rate classes")))
  })

  output$pathway_scope <- renderUI({
    p <- pathway()
    tags$p(class = "scope-summary", tags$strong(p$metadata$name %||% "Custom pathway"),
      paste0(" · ", nrow(p$nodes), " glycans · ", nrow(p$edges), " reactions · ", length(p$rate_classes), " rate classes"))
  })
  output$rate_selection_note <- renderUI(tags$p(class = "selection-note",
    paste0(length(intersect(input$rate_classes, pathway()$rate_classes)), " selected · × to remove; type to add")))
  output$glycan_selection_note <- renderUI(tags$p(class = "selection-note",
    paste0(length(intersect(input$glycoforms, pathway()$observed_classes)), " selected · × to remove; type to add")))
  pathway_coordinates <- reactive(pathway_map_coordinates(pathway(), is_canonical()))
  output$design_pathway_map <- renderUI({
    pathway_map_svg(pathway(), pathway_coordinates(), input$rate_classes %||% character(),
      input$glycoforms %||% character(), is_canonical(), input$pathway_zoom %||% 100)
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

  outputOptions(output, "panel_selector", suspendWhenHidden = FALSE)
  outputOptions(output, "rate_selector", suspendWhenHidden = FALSE)

  observeEvent(pathway(), {
    p <- pathway()
    panel <- if (is_canonical()) intersect(default_panel, p$observed_classes) else p$observed_classes
    updateSelectizeInput(session, "rate_classes", choices = rate_labels(p$rate_classes), selected = p$rate_classes)
    updateSelectizeInput(session, "glycoforms", choices = p$observed_classes, selected = panel)
  }, ignoreInit = TRUE)

  observeEvent(input$select_all_rates, {
    updateSelectizeInput(session, "rate_classes", choices = rate_labels(pathway()$rate_classes), selected = pathway()$rate_classes)
  })
  observeEvent(input$clear_rates, updateSelectizeInput(session, "rate_classes", selected = character()))
  observeEvent(input$select_all_glycoforms, {
    updateSelectizeInput(session, "glycoforms", choices = pathway()$observed_classes, selected = pathway()$observed_classes)
  })
  observeEvent(input$clear_glycoforms, updateSelectizeInput(session, "glycoforms", selected = character()))

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
                           selected = result$panel)
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

  fit_sample_loaded <- reactiveVal(FALSE)
  observeEvent(input$load_fit_example, fit_sample_loaded(TRUE))
  observeEvent(input$observed_file, fit_sample_loaded(FALSE), ignoreInit = TRUE)
  fit_uses_example <- reactive(fit_sample_loaded() || is.null(input$observed_file))
  observed_data <- reactive({
    if (!fit_sample_loaded() && is.null(input$observed_file)) return(NULL)
    csv_path <- if (fit_sample_loaded()) system.file("extdata", "example_observed.csv", package = "glycoPathDesign") else input$observed_file$datapath
    observed <- tryCatch(read.csv(csv_path, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) NULL)
    validate(need(!is.null(observed) && all(c("glycoform", "proportion") %in% names(observed)),
                  "The CSV needs glycoform and proportion columns."))
    validate(need(nrow(observed) >= 2 && !anyNA(observed$glycoform) && !anyDuplicated(observed$glycoform),
                  "Supply at least two distinct glycoforms, with one row per glycoform."))
    validate(need(is.numeric(observed$proportion) && all(is.finite(observed$proportion)) &&
                    all(observed$proportion >= 0) && sum(observed$proportion) > 0,
                  "Proportions must be finite, non-negative numbers with a positive sum."))
    unknown <- setdiff(observed$glycoform, pathway()$observed_classes)
    validate(need(!length(unknown), paste("Glycoforms not in the pathway:", paste(unknown, collapse = ", "))))
    observed[, c("glycoform", "proportion"), drop = FALSE]
  })
  fitted <- eventReactive(input$run_fit, {
    req((input$run_fit %||% 0) > 0)
    observed <- fit_display_data()
    x <- withProgress(message = "Estimating reaction rates…", value = 0.3, {
      tryCatch(fit_pathway(pathway(), observed, n_starts = input$fit_starts %||% 8),
               error = function(e) validate(need(FALSE, conditionMessage(e))))
    })
    x$is_example <- fit_uses_example()
    x$input_data <- observed
    x
  }, ignoreInit = FALSE)
  current_fit <- reactive({
    if ((input$run_fit %||% 0) == 0) return(NULL)
    x <- fitted()
    if (!identical(x$pathway, pathway()) || !identical(x$input_data, fit_display_data()) ||
        !identical(x$is_example, fit_uses_example())) return(NULL)
    x
  })
  fit_example_data <- reactive({
    composition <- predict_composition(pathway())
    data.frame(glycoform = names(composition), proportion = unname(composition))
  })
  output$fit_csv_example <- renderTable({
    fit_example_data()
  }, digits = 4, spacing = "xs")
  fit_display_data <- reactive({
    d <- observed_data()
    if (!is.null(d)) return(d)
    fit_example_data()
  })
  output$fit_proportion_input <- renderUI({
    d <- fit_display_data()
    example <- fit_uses_example()
    percent <- 100*d$proportion/sum(d$proportion)
    div(class = "fit-proportion-input",
      tags$strong(if (fit_sample_loaded()) "Loaded: example_observed.csv" else if (example) "Example proportions · entire pathway" else paste(nrow(d), "measured glycans loaded")),
      div(class = "fit-proportion-rows", lapply(seq_len(nrow(d)), function(i) {
        div(class = "fit-proportion-row",
          tags$span(class = "glycan-name", d$glycoform[i]),
          div(class = "proportion-track", `aria-hidden` = "true",
            div(class = "proportion-fill", style = paste0("width:", percent[i], "%;"))),
          tags$strong(paste0(format_number(percent[i]), "%")))
      })),
      tags$small(if (example) "Simulated data, not measurements. Estimate these example rates, or upload your own values below."
                 else "Your values, normalized over the uploaded glycans.")
    )
  })
  output$fit_pathway_context <- renderUI({
    p <- pathway()
    tags$p(class = "small-note", tags$strong(p$metadata$name %||% "Custom pathway"),
      paste0(" · ", nrow(p$nodes), " glycans · ", length(p$rate_classes), " unknown rate parameters"))
  })
  output$fit_diagram_key <- renderUI({
    x <- current_fit()
    if (!is.null(x)) return(tags$p(class = "small-note fit-key",
      if (x$is_example) "Example fit · blue = simulated proportions · red numbers = estimated rates (not experimental results)"
      else "Fit completed · blue = your measured proportions · red numbers = estimated rates"))
    tags$p(class = "small-note fit-key",
      if (fit_uses_example()) "Simulated example for all glycans · Blue = proportions · Red ? = rates to estimate"
      else "Blue bars = your measured proportions · Red = reaction-rate estimates / unknowns")
  })
  output$fit_pathway_diagram <- renderUI({
    p <- pathway(); x <- current_fit(); d <- fit_display_data()
    pathway_map_svg(p, pathway_coordinates(), p$rate_classes,
      d$glycoform, is_canonical(), 100,
      fit_mode = TRUE, rate_values = if (is.null(x)) NULL else x$rates,
      proportions = setNames(d$proportion/sum(d$proportion), d$glycoform))
  })
  output$fit_state_message <- renderUI({
    x <- current_fit()
    if (is.null(x)) return(NULL)
    messages <- character()
    if (x$convergence != 0) messages <- c(messages,
      "The optimiser did not converge. The displayed rates are provisional; inspect the fit before interpreting them.")
    if (!x$audit$full_rank) messages <- c(messages,
      "These values are one fitted solution. This measurement panel cannot identify all rates separately; do not interpret the numbers as unique estimates.")
    if (!length(messages)) messages <- "Estimates calculated. Check agreement and precision below. Full local rank does not establish global uniqueness."
    tags$p(class = "conditional-note", paste(messages, collapse = " "))
  })

  output$fit_summary <- renderUI({
    if ((input$run_fit %||% 0) == 0) {
      return(tags$p(class = "small-note", "Upload a within-site composition to assess whether the pathway represents the observed glycoforms."))
    }
    x <- req(current_fit())
    div(class = "technical-grid",
        div(tags$strong(format_number(x$fit_correlation)), tags$span("fit correlation")),
        div(tags$strong(paste0(x$audit$rank, "/", x$audit$n_rate_classes)), tags$span("local rank at fitted point")),
        div(tags$strong(format_number(x$audit$condition_number)), tags$span("condition number at fit")),
        div(tags$strong(format_number(x$objective)), tags$span("fit objective")))
  })
  output$fit_plot <- renderPlot(plot(req(current_fit())))
  output$rate_table <- renderTable({
    x <- req(current_fit())
    data.frame(rate_class = names(x$rates), effective_rate = unname(x$rates), log_rate = unname(x$log_rates))
  }, digits = 4, rownames = FALSE)
}

shinyApp(ui, server)
