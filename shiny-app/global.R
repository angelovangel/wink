require(shiny)
require(shinydashboard)
require(shinyFiles)
require(shinypop) # remotes::install_github("dreamRs/shinypop")
require(shinyjs)
require(rmarkdown)
require(data.table)
require(DT)
require(dplyr)
require(sys) # https://github.com/jeroen/sys
require(here)

 si_fmt <- function(x) { 
 	system2(here( "bin/si-format.sh"), x, stdout = TRUE) 
 }

 