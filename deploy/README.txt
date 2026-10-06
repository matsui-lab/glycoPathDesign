glycoPathDesign deployment entry point

This folder contains app.R, a thin launcher for the installed package.

1. Install R (>= 4.2.0), shiny, DT and igraph for the account running the app.
2. Install glycoPathDesign_0.6.0.tar.gz with R CMD INSTALL.
3. For a local check, run shiny::runApp("deploy", host="127.0.0.1", port=3838).
4. For Shiny Server, copy deploy/app.R into the chosen application directory.
   Ensure the Shiny Server service account can load the same package and dependencies.
5. Verify the deployed URL, model selection, example fit, pooled measurement
   display and asset loading before adding a hosted URL to the manuscript.

The launcher was tested locally on macOS with R 4.5.2 and Shiny 1.12.1.
Deployment on a production Linux server has not been tested. No hosted URL is
assigned here. A shinyapps.io deployment must resolve the package dependency
after its public release; this local-package launcher alone is not proof that
a cloud dependency build succeeds.

Official references:
https://shiny.posit.co/r/reference/shiny/latest/shinyapp.html
https://docs.posit.co/shiny-server/
