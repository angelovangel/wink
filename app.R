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
							 						 column(width = 12,
							 							checkboxInput("testrun", "Simulate run with test data", width = "100%"),
							 							checkboxInput("skip_kraken", "Skip kraken2, show only run statistics", width = "100%"),
							 							checkboxInput("weakmem", "Do not load kraken2 database in RAM (use on weak machines)")
							 							),
							 						 column(width = 6,
							 						 selectizeInput("nxf_profile", 
							 						 							 "Select nextflow profile", 
							 						 							 width = "100%",
							 						 							 choices = c("docker", "conda", "local"), 
							 						 							 selected = "docker"),
							 						 ),
							 						 column(width = 6,
							 						 		selectizeInput("taxlevel", 
							 						 									 "Taxonomic level for kraken2 abundance", 
							 						 									 width = "100%",
							 						 									 choices = c("Domain" = "D", "Phylum" = "P", "Class" = "C", "Order" = "O", "Family" = "F", "Genus" = "G", "Species" = "S"), 
							 						 									 selected = "S")
							 						 ),
							 						 column(width = 6,
							 						 textInput("kraken_gz", 
							 						 					width = "100%",
							 						 					"Path to kraken2 database (gz file)", 
							 						 					value = "ftp://ftp.ccb.jhu.edu/pub/data/kraken2_dbs/minikraken_8GB_202003.tgz")
							 						 )
							 		),
							 		tags$hr(),
							 		column(widt = 12,
							 			verbatimTextOutput("stdout"),
							 			verbatimTextOutput("log_output")
							 		)
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
							 	box(width = 12, 
							 			status = "warning", solidHeader = FALSE, collapsible = TRUE,collapsed = FALSE,
							 			title = "Run statistics per sample",
							 			DT::dataTableOutput("stats", width = "100%", height = 500)
							 	)
							 )
			),#---------------------------------------------------------------
			tabPanel("Taxonomy and abundance", #-----------------------------
							 fluidRow(
							 	box(width = 12, status = "warning", solidHeader = FALSE, collapsible = TRUE,
							 			title = "Kraken read assignment per sample",
							 			infoBoxOutput("current_barcode", width = 3),
							 			infoBoxOutput("all_reads", width = 3),
							 			infoBoxOutput("ass_reads", width = 3),
							 			infoBoxOutput("unass_reads", width = 3)
							 	)
							 ),
							 fluidRow(
							 	box(width = 6, 
							 			status = "warning", solidHeader = FALSE, collapsible = TRUE,collapsed = FALSE,
							 			title = "Bracken results per sample", 
							 			sliderInput("topn", "Top N species to show", value = 3, min = 1, max = 5, step = 1),
							 			DT::dataTableOutput("ab_table", width = "100%", height = 500)
							 			),
							 	box(width = 6,
							 			status = "warning", solidHeader = FALSE, collapsible = TRUE, collapsed = FALSE,
							 			title = "Detailed bracken results",
							 			shinyWidgets::sliderTextInput("filterFreq", 
							 																		"Filter by abundance", 
							 																		choices = c(0, 0.1, 1, 5, 10, 20), 
							 																		selected = 0,
							 																		post = "%", grid = TRUE),
							 			#sliderInput("filterFreq", "Filter by abundance", min = 0, max = 1, value = 0, step = 0.1),
							 			DT::dataTableOutput("ab_table_detail", width = "100%", height = 500))
							 )
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
	
	#---create and read watch file that accumulates ps monitor of the nxf process
	nxf_logfile <- tempfile() 
	system2( "touch", c(nxf_logfile, "&& echo 'ps monitor for the wink nextflow process' >", nxf_logfile) )
	nxf_logfile_data <- reactiveFileReader( 1000, session, nxf_logfile, readLines )
	
	#---the same for the nxf output
	nxf_outfile <- tempfile()
	system2( "touch", c(nxf_outfile, "&& echo 'wink nextflow pipeline output' >", nxf_outfile) )
	nxf_outfile_data <- reactiveFileReader( 1000, session, nxf_outfile, readLines )
	
	
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
	brackenFiles <- reactivePoll(
		1000, session, 
		checkFunc = function() {
			if( file.exists("results-wink/latest-bracken") ) {
				file.info("results-wink/latest-bracken")$mtime
			} else {
				""
			}
			}, 
		valueFunc = function() {
			l <- list.files("results-wink/latest-bracken", pattern = "*.tsv", full.names = TRUE)
			setNames(l, basename(l)) # also returns l
		}
	)
	
	brackenData <- reactivePoll(
		1000, session,
		checkFunc = function() {
			if ( any( file.exists(brackenFiles()) ) ) 
				file.info( brackenFiles() )$mtime
			else
				""
			}, 
		valueFunc =  function(){
			if ( any(file.exists(brackenFiles() )) ) {
				
				rbindlist( lapply(brackenFiles(), fread), idcol = "file" )
			} else {
				brackenColNames <- c("file", "name", "taxonomy_id", "kraken_assigned_reads", "new_est_reads", "freq")
				setNames(data.frame(matrix(ncol = 6, nrow = 0, 0)), brackenColNames)
			}
		})
	
	
	# this is the summary table, one row per barcode
	brackenDataLeft <- reactive({
		brackenData() %>% 
		dplyr::filter(name != "unclassified") %>% # remove these here, they are used for the krakenData reactives
		mutate(freq = round(freq, 3)) %>%
		group_by(file) %>%
		slice_head(n = input$topn) %>%
		summarize_at( c("name"), .funs = "paste", collapse = " | " )
	})
	
	# reactive vals for storing total, mapped reads, nxf process info...
	seqData <- reactiveValues(nsamples = 0, treads = 0, tbases = 0, n50 = 0, runtime = 0)
	krakenData <- reactiveValues(all_reads = 0, assigned_reads = 0, unassigned_reads = 0)
	nxf <- reactiveValues(pid = NULL, watch = NULL)
	
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
			cat("wink command preview, select a fastq_pass folder to start\n")
		} else {
			# hard set fastq folder
		selectedFolder <<- parseDirPath(volumes, input$fastq_pass_folder)
		skip_kraken <<- ifelse(input$skip_kraken, "--skip_kraken", "")
		
		nxf_args <<- c("run" ,
									 "main.nf",
									 "--fastq_pass", selectedFolder,
									 skip_kraken,
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
			nxf$pid <- sys::exec_background("nextflow", 
																	args = nxf_args, 
																	std_out = nxf_outfile
																	)
			
			nx_notify_success(paste("Nextflow pipeline started with pid", nxf$pid))
			
			nxf$watch <- sys::exec_background("bin/watch-pid.sh", 
										 args = nxf$pid, 
										 std_out = nxf_logfile
										 )
			
			shinyjs::enable("stop")
			shinyjs::disable(id = "run")
			shinyjs::toggleCssClass("stop", "yellow")
			
			shinyjs::html(selector = ".logo", 
										html = paste("<p style='background-color:#E67E22;'>Nextflow pipeline running with pid ", 
																 nxf$pid, "</p>")
										)
			
			writeLines("", ".pxout", )
			
		}
	})
	
	# and kill 
	observeEvent(input$stop, {
		
		if (nxf$pid != 1 && nxf$watch != 1) { # don't kill pideinz please:)
			tools::pskill(nxf$pid)
			#tools::pskill(nxf$watch)
			
			nx_notify_warning(paste("Nextflow pipeline with pid", nxf$pid, "was stopped!"))
			shinyjs::enable("run")
			shinyjs::toggleCssClass("stop", "yellow")
			shinyjs::disable("stop")
			shinyjs::html(selector = ".logo", 
										html = "WINK - What's In my Nanopore reads, with Kraken2, in real-time")
		}
	})
	
	
	# ------------------------------------------------------------------
	
	 output$log_output <- renderPrint({
	 	cat({ nxf_logfile_data() %>% tail(2)}, sep = "\n" )
	 	})
	
	output$nxf_output <- renderPrint({
		cat(nxf_outfile_data(), sep = "\n")
		#runjs("document.getElementById('nxf_output').scrollTo(0,1e9);") # scroll the page to bottom with each message, 1e9 is just a big number
	})
	
	# value boxes outputs for stats tab--------------------------
	output$stats <- DT::renderDataTable({
		df <- statsData() %>% 
			dplyr::mutate(file = basename(tools::file_path_sans_ext(file)),
										bases = sum_len,
										bases_h = si_fmt(bases)) %>%
			dplyr::select(file, num_seqs, bases, bases_h, max_len, N50, Q20_perc, last_write)
		
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
			subtitle = "Nanopore running time",
			color = 'light-blue'
		)
	})
	
	# value boxes outputs for kraken tab---------------------------
	
	output$current_barcode <- renderInfoBox({
		infoBox(
			 value = ifelse(is.null(last_selection$row_value), 
			 							 "Select a sample", 
			 							 last_selection$row_value
			 							 ),
			 title = "Sample", 
			icon = icon("vial"),
			color = 'light-blue'
		)
	})
	
	output$all_reads <- renderInfoBox({
		infoBox(
			value = ifelse(is.null(last_selection$row_value), 
										 0, 
										 krakenData$all_reads),
			title = "All reads",
			subtitle = last_selection$row_value,
			icon = icon("bars"),
			color = 'light-blue'
		)
	})
	
	output$ass_reads <- renderInfoBox({
		infoBox(
			value = ifelse(is.null(last_selection$row_value), 
										 0, 
										 paste(krakenData$assigned_reads, 
										"(", 
										round(krakenData$assigned_reads/sum(krakenData$assigned_reads, krakenData$unassigned_reads)*100, 0),
										"%)"
										)
							),
			title = "Assigned reads", 
			subtitle = last_selection$row_value,
			icon = icon("database"),
			color = 'light-blue'
		)
	})
	
	output$unass_reads <- renderInfoBox({
		infoBox(
				value = ifelse(is.null(last_selection$row_value), 
											 0, 
											 paste(krakenData$unassigned_reads, 
											 			"(", 
											 			round(krakenData$unassigned_reads/sum(krakenData$assigned_reads, krakenData$unassigned_reads)*100, 0),
											 			"%)"
											 )
				),
				title = "Unassigned reads",
				subtitle = last_selection$row_value,
				icon = icon("question"),
				color = 'light-blue'
		)
	})
	
	
	# retain selected row in a temp reactive variable to enable persistent selection
	# without this, the row is de-selected after each referesh of the table
	last_selection <- reactiveValues(row_value = NULL)
	
	# store the values each time a row is selected
	observeEvent(input$ab_table_rows_selected,{ 
		
		last_selection$row_value = brackenDataLeft()$file[input$ab_table_rows_selected]
		#print(last_selection$row_value)
		
	})
	
	#------------------------------------------------
	output$ab_table <- DT::renderDataTable({
		
		# which is the number of the previous row now?
		newselection <- which(brackenDataLeft()$file %in% last_selection$row_value)
		
		DT::datatable(brackenDataLeft(), 
									selection = list(mode = "single", selected = newselection), #keep persistant row selection!
									caption = HTML(paste("Top" ,tags$b(input$topn), "kraken-assigned species per sample (click on a row for more info)")),
									#filter = 'top',
									extensions = 'Buttons', 
									options = list(dom = 'Btp', 
																 buttons = c('copy', 'csv', 'excel')
									), 
									rownames = FALSE, class = 'hover row-border')
	})
	
	taxdb_left = "https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id="
	taxdb_right = "&lvl=3&lin=f&keep=1&srchmode=1&unlock"
	
	output$ab_table_detail <- DT::renderDataTable({
		# validate ab_table is not empty
		validate(
			need(!is.null(last_selection$row_value), "No data here"),
			need( nrow(brackenData() ) != 0, "No data here" )
		)
		
		# use the row value reactive to filter
		df_assigned <- brackenData() %>% 
			dplyr::filter(file == last_selection$row_value, name != "unclassified") %>%
			
			dplyr::mutate(
				taxonomyID = paste0(taxdb_left, taxonomy_id, taxdb_right), #two-step mutate to build href
				kraken_reads = kraken_assigned_reads,
				bracken_corr_reads = new_est_reads, 
				freq = round(freq*100, 2)
				) %>% 
			dplyr::mutate(
				taxonomyID = paste0("<a target = '_blank' href='", taxonomyID, "'>", taxonomy_id, "</a>")
				) %>%
			
			dplyr::select(name, taxonomyID, kraken_reads, bracken_corr_reads, freq) %>%
			dplyr::filter(freq >= input$filterFreq)
		
		df_unassigned <- brackenData() %>% 
			dplyr::filter(file == last_selection$row_value, name == "unclassified")
		
		# set kraken statistics for value box
		krakenData$assigned_reads = sum( df_assigned$kraken_reads )
		krakenData$unassigned_reads = sum( df_unassigned$kraken_assigned_reads)
		krakenData$all_reads = sum(df_assigned$kraken_reads, df_unassigned$kraken_assigned_reads)
		
		caption <- if_else(is.na(last_selection$row_value), "Select a sample from the table on the left", 
											 paste("Bracken abundance table for", tags$b(last_selection$row_value)))
		DT::datatable(df_assigned, 
									escape = FALSE,
									selection = "single",
									caption = HTML(caption),
									#caption = HTML(paste("Bracken abundance table", tags$b(rowselected))),
									#filter = 'top',
									extensions = 'Buttons', 
									options = list(dom = 'Btp', 
																 buttons = c('copy', 'csv', 'excel')
									), 
									rownames = FALSE, class = 'hover row-border') %>%
		DT::formatStyle('freq',
										background = styleColorBar(c(0, 100), 'skyblue'),
										backgroundSize = '95% 70%',
										backgroundRepeat = 'no-repeat',
										backgroundPosition = 'right') %>%
		DT::formatStyle('bracken_corr_reads',
											background = styleColorBar(c(0, max(df_assigned$bracken_corr_reads, na.rm = T)), 'skyblue'),
											backgroundSize = '95% 70%',
											backgroundRepeat = 'no-repeat',
											backgroundPosition = 'right') %>%
		DT::formatStyle('kraken_reads',
											background = styleColorBar(c(0,max(df_assigned$kraken_reads, na.rm = T)), 'skyblue'),
											backgroundSize = '95% 70%',
											backgroundRepeat = 'no-repeat',
											backgroundPosition = 'right')
	})
	
	session$onSessionEnded(function() {
		isolate( tools::pskill(nxf$watch) )
		unlink(nxf_logfile)
		unlink(nxf_outfile)
		#system2("rm", args = c("-rf", "work")) # comment out on deploy
	})
	
}
shiny::shinyApp(ui, server)
