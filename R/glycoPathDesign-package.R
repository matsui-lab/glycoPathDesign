#' glycoPathDesign: measurement design for glycoform pathway inference
#'
#' `glycoPathDesign` asks whether a proposed site-specific glycoform panel can
#' distinguish the effective reaction-rate classes specified by a biosynthetic
#' pathway model, and what assay precision that distinction requires.
#'
#' The main workflow is:
#'
#' 1. Define a pathway with [glyco_pathway()].
#' 2. Inspect its forward prediction with [predict_composition()].
#' 3. Audit a proposed panel with [check_panel()].
#' 4. Test operating-point robustness with [check_operating_points()].
#' 5. Search alternatives with [enumerate_panels()] and [optimise_panels()].
#' 6. Check representation of observed data with [fit_pathway()].
#'
#' Outputs are local and conditional on the supplied biochemical graph,
#' grouping of edges into effective rate classes, entry and secretion model,
#' observation mapping, operating point, and error assumptions.
#'
#' @keywords internal
#' @importFrom stats setNames
"_PACKAGE"
