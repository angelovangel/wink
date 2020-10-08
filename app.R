#-----------------#
# Main WINK app
#-----------------#
# reactive read of stats and kraken output
# nextflow pipeline has to be started separately in the same folder, i.e. the shiny app reads results-wink

require(shiny)
require(shinydashboard)
require(shinyFiles)
require(shinypop)
require(shinyjs)
require(data.table)
require(DT)
require(dplyr)
require(sys) # https://github.com/jeroen/sys

si_fmt <- function(x) { system2("bin/si-format.sh", x, stdout = TRUE) }

ui <- dashboardPage(title = "WINK",
	dashboardHeader(title = "WINK - What's In my Nanopore reads, with Kraken2, in real-time", 
									titleWidth = '100%', 
									dropdownMenuOutput("notificationsMenu")
									), 
	dashboardSidebar(disable = TRUE),
	
	# 		shinyWidgets::sliderTextInput(
	# 			"abundfilter",
	# 			"Filter by abundance",
	# 			choices = c(0, 0.1, 1, 10, 20),
	# 			selected = 0,
	# 			post = "%",
	# 			grid = TRUE
	# 		)
	# 	)
	# ),
	dashboardBody(
		includeCSS("custom.css"),
		useShinyjs(),
		#useShinyalert(),
		use_notiflix_notify(position = "right-top", width = "380px"),
		
# 		tags$head(tags$style(HTML('
#         .skin-blue .main-header .logo {
#           background-color: light-blue;
#         }
#         .skin-blue .main-header .logo:hover {
#           background-color: light-blue;
#         }
#       '))),
		tabsetPanel(
			tabPanel("Sequencing overview", #--------------------------------------------------------------
							 fluidRow(
							 	box(width = 12, status = "warning", solidHeader = FALSE, collapsible = TRUE, 
							 			title = "Control panel",
							 		
							 		shinyDirButton(id = "fastq_pass_folder", 
							 									 label = "Select run folder", 
							 									 title = "Select the fastq_pass folder",
							 									 style = "color: #3498DB;",
							 									 #color = "primary", 
							 									 icon = icon("folder-open")),
							 		actionButton("reset", "Reset", style = "color: #3498DB;", icon = icon("sync"), onclick = "history.go(0);"),
							 		actionButton("run", "Start WINK", style = "color: #3498DB;", icon = icon("play")),
							 		actionButton("stop", "Stop WINK", style = "color: #3498DB;", icon = icon("power-off")),
							 		actionButton("more", "More options", style = "color: #3498DB;", icon = icon("cog"), class = "rightAlign"),
							 		tags$div(id = "optional_inputs",
							 						 checkboxInput("testrun", "Start test run"),
							 						 selectizeInput("nxf_profile", 
							 						 							 "Select nextflow profile", 
							 						 							 width = "40%",
							 						 							 choices = c("docker", "conda", "local"), 
							 						 							 selected = "docker"),
							 						 textInput("kraken_gz", 
							 						 					width = "40%",
							 						 					"Path to kraken2 database (gz file)", 
							 						 					value = "ftp://ftp.ccb.jhu.edu/pub/data/kraken2_dbs/minikraken_8GB_202003.tgz"),
							 						 
							 						 selectizeInput("taxlevel", 
							 						 							 "Taxonomic level for kraken2 abundance estimation", 
							 						 							 width = "40%",
							 						 							 choices = c("Domain" = "D", "Phylum" = "P", "Class" = "C", "Order" = "O", "Family" = "F", "Genus" = "G", "Species" = "S"), 
							 						 							 selected = "S")
							 						 
							 						 ),
							 		verbatimTextOutput("stdout")
							 	)
							 ),
							 fluidRow(
							 	box(width = 12, status = "warning", solidHeader = FALSE, collapsible = TRUE,
							 			title = "Run statistics",
							 	valueBoxOutput("nsamples", width = 3),
							 	valueBoxOutput("treads", width = 2),
							 	valueBoxOutput("tbases", width = 2),
							 	valueBoxOutput("n50", width = 2),
							 	valueBoxOutput("runtime", width = 3)
								)
							 ),
							 fluidRow(
							 	box(
							 		width = 12, status = "warning", solidHeader = FALSE, collapsible = TRUE,collapsed = TRUE,
							 		title = "Run statistics per sample",
							 		DT::dataTableOutput("stats", width = "100%", height = 500)
							 	)
							 )
			),#---------------------------------------------------------------
			tabPanel("Taxonomy and abundance", #-----------------------------
							 numericInput("maxrows", "Rows to show", 25),
							 tableOutput("ab_table"),
							 downloadButton("downloadCsv", "Download as CSV")
			), #---------------------------------------------------------------
			tabPanel("Nextflow output",
							 verbatimTextOutput("nxf_output")
			),
			tabPanel("Help")
		)
	)
)

server <- function(input, output, session) {
	
	options(shiny.launch.browser = TRUE)
	
	nx_notify_success(paste("Hello ", Sys.getenv("LOGNAME"))
	)
	
	# handling of shinyFiles
	volumes <- c(Home = fs::path_home(), getVolumes()() )
	shinyDirChoose(input, "fastq_pass_folder", 
								 roots = volumes, 
								 session = session)
	
	headers <- c("file", "format", "type", "num_seqs", 
							 "sum_len", "min_len", "avg_len", "max_len", 
							 "Q1", "Q2", "Q3", "sum_gap", "N50", "Q20_perc", "Q30_perc", 
							 "first_write", "last_write")
	
	# REACTIVES ------------------------------------------------------------------
	# reactive to get actual list of stats files
	statsFiles <- reactivePoll(
		1000, session,
		checkFunc = function() {
			if( file.exists("results-wink/latest-stats") ) {
				file.info("results-wink/latest-stats")$mtime
			} else {
				""
				}
		},
		valueFunc = function() {
			list.files("results-wink/latest-stats", pattern = "*stats.txt", full.names = TRUE)
		}
		)
	
	#reactive to fread stats files
	statsData <- reactivePoll(
		1000, session, 
		checkFunc = function() {
			if ( any(file.exists(statsFiles() )) )
				file.info( statsFiles() )$mtime
			else
				""
		}, 
		valueFunc = function() {
			if ( any(file.exists(statsFiles() )) )
				rbindlist( lapply(statsFiles(), fread, sep = " ", col.names = headers))
			else setNames(data.frame(matrix(ncol = 17, nrow = 0, 0)), headers)
				
		}
		)
	pxout <- reactiveFileReader(1000, session, ".pxout", readLines)

	# reactive vals for storing total, mapped reads, nxf process info...
	seqData <- reactiveValues(nsamples = 0, treads = 0, tbases = 0, n50 = 0, runtime = 0)
	nxf <- reactiveValues(pid = 0)
	
	# OBSERVERS------------------------------------------------------------------
	# observer for optional inputs
	hide("optional_inputs")
	observeEvent(input$more, {
		shinyjs::toggle("optional_inputs")
	})
	shinyjs::disable("stop")
	
	observe({
		seqData$nsamples <- nrow( statsData() )
		seqData$treads <- si_fmt( sum(statsData()$num_seqs, na.rm = TRUE) )
		seqData$tbases <- si_fmt( sum(statsData()$sum_len, na.rm = TRUE) )
		seqData$n50 <-  mean(statsData()$N50, na.rm = TRUE)
		seqData$runtime <- difftime( max(statsData()$last_write), min(statsData()$first_write) ) 
	})
	
	
	# RENDERS ------------------------------------------------------------------
	
	output$stdout <- renderPrint({
		
		# build nxf call and print to stdout
		if (is.integer(input$fastq_pass_folder)) {
			cat("No fastq_pass folder selected\n")
		} else {
			# hard set fastq folder
		selectedFolder <<- parseDirPath(volumes, input$fastq_pass_folder)
		
		nxf_args <<- c("run" ,
									 "main.nf",
									 "--fastq_pass", selectedFolder, 
									 "--kraken_gz", input$kraken_gz)
		cat("nextflow", nxf_args)
		}
	})
	# CALLS TO NEXTFLOW PIPELINE ------------------------------------------------------------------
	
	# start
	observeEvent(input$run, {
		if(is.integer(input$fastq_pass_folder)) {
			nx_notify_error("Select run folder first!")
		} else {
			
			pid <- sys::exec_background("nextflow", 
																	args = nxf_args, 
																	std_out = ".pxout")
			nxf$pid <- pid
			nx_notify_success(paste("Nextflow pipeline started with pid:", pid))
			
			shinyjs::enable("stop")
			shinyjs::disable(id = "run")
			shinyjs::toggleCssClass("stop", "yellow")
			
			shinyjs::html(selector = ".logo", 
										html = paste("<p style='background-color:#E67E22;'>Nextflow pipeline running, pid: ", 
																 nxf$pid, "</p>")
										)
			
			writeLines("", ".pxout")
			
		}
	})
	
	# and kill 
	observeEvent(input$stop, {
		if (nxf$pid != 0) {
			tools::pskill(nxf$pid)
			nx_notify_warning(paste("Nextflow pipeline with pid:", nxf$pid, "was stopped!"))
			shinyjs::enable("run")
			shinyjs::toggleCssClass("stop", "yellow")
			shinyjs::disable("stop")
			shinyjs::html(selector = ".logo", 
										html = "WINK - What's In my Nanopore reads, with Kraken2, in real-time")
		}
	})
	
	
	# ------------------------------------------------------------------
	output$ab_table <- renderTable({
		
	})
	
	 output$nxf_output <- renderPrint({
	 	pxout()
	 	#runjs("document.getElementById('nxf_output').scrollTo(0,1e9);") # scroll the page to bottom with each message, 1e9 is just a big number
	 })
	#
	output$stats <- renderDataTable({
		df <- statsData() %>% 
			dplyr::mutate(bases = sum_len,
										bases_human = si_fmt(bases)) %>%
			dplyr::select(file, num_seqs, bases, bases_human, max_len, N50, Q20_perc, last_write)
		
			datatable(df, filter = 'top',
								extensions = 'Buttons', 
								options = list(dom = 'Btp', 
															 buttons = c('copy', 'csv', 'excel')
															 ), 
								rownames = FALSE, class = 'hover row-border') %>%
				DT::formatStyle('num_seqs',
												background = styleColorBar(c(0, max(df$num_seqs)), 'skyblue'),
												backgroundSize = '95% 70%',
												backgroundRepeat = 'no-repeat',
												backgroundPosition = 'right') %>%
				DT::formatStyle('bases',
												background = styleColorBar(c(0, max(df$bases)), 'skyblue'),
												backgroundSize = '95% 70%',
												backgroundRepeat = 'no-repeat',
												backgroundPosition = 'right') %>%
				DT::formatStyle('max_len',
												background = styleColorBar(c(0, max(df$max_len)), 'skyblue'),
												backgroundSize = '95% 70%',
												backgroundRepeat = 'no-repeat',
												backgroundPosition = 'right') %>%
				DT::formatStyle('N50',
												background = styleColorBar(c(0, max(df$N50)), 'skyblue'),
												backgroundSize = '95% 70%',
												backgroundRepeat = 'no-repeat',
												backgroundPosition = 'right') 
	})
	
	output$nsamples <- renderValueBox({
		valueBox(
			value = seqData$nsamples,
			subtitle = "Detected barcodes", 
			color = 'light-blue'
		)
	})
	
	output$treads <- renderValueBox({
		valueBox(
			value = seqData$treads, 
			subtitle = "Total reads",
			color = 'light-blue'
		)
	})
	
	output$tbases <- renderValueBox({
		valueBox(
			value = seqData$tbases, 
			subtitle = "Total bases",
			color = 'light-blue'
		)
	})
	
	output$n50 <- renderValueBox({
		valueBox(
			value = paste( round(seqData$n50, digits = 0), "bp"), 
			subtitle = "Mean N50",
			color = 'light-blue'
		)
	})
	
	output$runtime <- renderValueBox({
		valueBox(
			value = paste( round(as.numeric(seqData$runtime, units = 'hours'), digits = 2), "hours" ),
			subtitle = "Running time",
			color = 'light-blue'
		)
	})
	
	session$onSessionEnded(function() {
		#system2("rm", args = c("-rf", "work")) # comment out on deploy
	})
	
}
shiny::shinyApp(ui, server)
