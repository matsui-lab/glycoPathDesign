# glycoPathDesign 0.6.0

- Adds the pathway input and visual guide.
- Separates fitted rates from fixed reference rates and displays pooled measurement totals once.
- Preserves distinct parallel and reverse reaction paths in the pathway diagram.
- Recognises starting solutions at floating-point precision, including CSV round trips, while retaining real nonconvergence warnings.
- Adds installation and Shiny launch instructions for the manuscript companion.

# glycoPathDesign 0.5.0

- Replaced free-form editing as the default pathway-building workflow with a
  registered N-glycan reaction selector that assigns source state, product
  state and effective rate class together.
- Separated novel states and rate classes into an explicit custom-definition
  mode with clearly labelled free-text identifiers.
- Changed the reaction table to a read-only summary with explicit add, delete,
  clear and restore actions, reducing accidental model-definition errors.

# glycoPathDesign 0.4.0

- Added an editable in-browser pathway builder so experimental users can define
  reactions without preparing CSV files.
- Reaction rows directly specify source state, product state, effective rate
  class and weight. State rows are generated automatically and expose measured-
  class mapping, entry amount and secretion.
- Added row creation and deletion, a resettable branched example and export of
  the constructed pathway as reusable package-compatible CSV files.
- Kept CSV upload as an advanced reuse route; built and uploaded pathways feed
  the same panel, precision, robustness and fit checks.

# glycoPathDesign 0.3.0

- Rebuilt the Shiny interface around one compact design-and-decision screen.
- Made user-supplied pathway models a primary workflow rather than a hidden
  technical option. A single multi-file input accepts nodes, edges, observation
  mappings and boundary conditions, with a downloadable generic template.
- Replaced explanatory cards with a decision table reporting panel dimension,
  reaction separation, analytical precision, conditioning and forward-model
  fit status.
- Added automatic next-step recommendations, on-demand alternative-panel
  search and direct transfer of the best evaluated panel into the design.
- Moved operating-point robustness, pathway graphics and observed-composition
  fitting into one secondary tab while retaining exact canonical SNFG previews.
- Restricted colour to small status accents and removed decorative surfaces,
  shadows and the four-step wizard.

# glycoPathDesign 0.2.0

- Replaced the analysis-tab-first Shiny interface with a four-step guided design
  workflow for experimental researchers.
- Added exact SNFG selection cards for all 22 canonical glycoform classes.
- Added plain-language Supported, Conditional, and Not supported verdicts with
  actionable recommendations and an explicit three-part design check.
- Added selection of the biosynthetic reaction classes that the experiment must
  distinguish; other pathway rates remain fixed in the local inverse problem.
- Added effective-CV calculation for independent technical replicates, a
  one-click complete example, candidate-panel suggestions, and a separate
  technical workspace retaining the full analysis controls.

# glycoPathDesign 0.1.0

- Added a validated pathway schema for nodes, reactions, effective rate classes,
  entry, secretion, and observation aggregation.
- Added steady-state forward prediction and numerical log-rate Jacobians.
- Added panel dimension, local rank, condition number, singular spectrum,
  compositional noise amplification, and assay-precision requirements.
- Added exhaustive or explicitly sampled panel search and optimisation.
- Added robustness audits across operating points.
- Added bounded multi-start fitting of pathway models to observed compositions.
- Added CSV, JSON, and HTML reports, base-R diagnostic graphics, an example
  canonical N-glycan network, tests, a vignette, and a Shiny application.
