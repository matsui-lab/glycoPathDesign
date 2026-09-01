# glycoPathDesign

`glycoPathDesign` designs site-specific glycoform measurements for inverse
problems in biosynthetic pathways. Given a reaction graph and a set of effective
rate classes, it answers three practical questions before measurement:

1. Does the proposed glycoform panel contain enough independent composition
   dimensions?
2. Are the specified rate directions locally distinguishable, and how strongly
   will measurement noise be amplified?
3. Can the forward model represent an observed composition at the operating
   point where rates would be interpreted?

The package contains an R API and a Shiny application. It supports latent
biosynthetic states, multiple reactions sharing one effective rate class, and
aggregation of multiple latent states into a measured glycoform class.

## Install

Install the development version from GitHub:

```r
remotes::install_github("matsui-lab/glycoPathDesign")
```

Install a local source archive when working offline:

```r
install.packages("glycoPathDesign_0.5.0.tar.gz", repos = NULL, type = "source")
```

## Minimal workflow

```r
library(glycoPathDesign)

pathway <- canonical_pathway()

panel <- c(
  "M9", "M8", "M5", "M5Gn", "M3Gn",
  "M3Gn2", "M3Gn3", "G1", "G2S1"
)

audit <- check_panel(
  pathway,
  panel,
  target_sd = 0.25,
  assay_cv = 0.01
)

audit
plot(audit)
plot_sensitivity(audit)
supported_resolution(audit)
```

The maximum permitted coefficient of variation under this model is
`audit$max_assay_cv`. The condition number compares the most and least
observable local rate directions. The noise-amplification factor `g` converts
relative intensity noise into worst-direction log-rate uncertainty under the
implemented compositional error model.

## Define a pathway

```r
nodes <- data.frame(id = c("precursor", "product_a", "product_b"))
edges <- data.frame(
  from = c("precursor", "precursor"),
  to = c("product_a", "product_b"),
  rate_class = c("branch_a", "branch_b")
)

pathway <- glyco_pathway(
  nodes,
  edges,
  entry = "precursor",
  secretion = 0.5
)
```

Input tables use the following columns:

| table | required columns | optional columns |
|---|---|---|
| nodes | `id` | `label` |
| edges | `from`, `to`, `rate_class` | `weight` |
| observations | `node`, `glycoform` | — |

Rows of `observations` may assign several latent states to the same measured
glycoform. Edge weights scale reactions that share a rate-class parameter.

## Select a panel

```r
candidates <- enumerate_panels(
  pathway = canonical_pathway(),
  panel_size = 9,
  max_panels = 5000,
  seed = 1
)

optimise_panels(candidates, objective = "amplification", n = 10)
```

Enumeration is exhaustive when the number of combinations is no larger than
`max_panels`; otherwise the function reports that it evaluated a reproducible
random sample. It never presents a sampled search as a proof of minimality.

## Fit and then audit an observed composition

```r
observed <- read.csv(
  system.file("extdata", "example_observed.csv", package = "glycoPathDesign")
)

fit <- fit_pathway(canonical_pathway(), observed, n_starts = 8)
fit
plot(fit)
```

Fit quality and local recoverability are reported separately. A flexible model
may reproduce the composition while leaving individual rate directions weakly
resolved.

## Shiny application

```r
run_app()
```

The main screen accepts either the manuscript example or a user-supplied
pathway and places the design inputs beside one decision table. Users choose the
reaction classes to distinguish, the glycoforms measured as one normalised
composition, assay CV and replicate count. The result separates panel dimension,
local reaction separation, analytical precision, conditioning and the still
required forward-model-fit check, then recommends the next design change.

Custom models are uploaded as `nodes.csv`, `edges.csv` and optional
`observations.csv` and `boundary.csv` files; a generic template is downloadable
from the app. Alternatively, the in-browser builder has two explicit input
modes. The default guided mode adds a registered N-glycan reaction as one
validated source-state, product-state and effective-rate-class combination.
The custom mode accepts user-entered state and rate-class identifiers, including
reuse of one effective parameter across several edges. The
corresponding state table is generated automatically and exposes the measured-
class mapping, entry amount and secretion. Constructed pathways can be exported
as reusable CSV definitions. Alternative panels can be searched on demand.
Operating-point robustness, observed-composition fitting and diagnostic plots
remain available in a secondary tab. The canonical example retains exact SNFG
previews, while custom glycoform and rate-class identifiers are unrestricted.

## Scope

`glycoPathDesign` evaluates *local, model-conditional* identifiability and
precision. It does not establish that a pathway graph is biochemically complete,
convert effective rate classes into unique enzyme activities, or replace
experimental validation. Operating-point robustness should be checked whenever
rate values are uncertain.

## Author

Yusuke Matsui
