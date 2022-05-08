

server <- function(input, output, session) {
	
	options(shiny.launch.browser = TRUE)
	
	nx_notify_success(paste("Hello ", Sys.getenv("LOGNAME"))
	)
	
	#---create and read watch file that accumulates ps monitor of the nxf process
	nxf_logfile <- tempfile() 
	system2( "touch", c(nxf_logfile, "&& echo 'ps monitor for the wink nextflow process' >", nxf_logfile) )
	nxf_logfile_data <- reactiveFileReader( 1000, session, nxf_logfile, readLines )
	
	#---the same for the nxf output
	nxf_outfile <- tempfile()
	system2( "touch", c(nxf_outfile, "&& echo 'wink nextflow pipeline output' >", nxf_outfile) )
	nxf_outfile_data <- reactiveFileReader( 1000, session, nxf_outfile, readLines )
	
	# BEGIN handling of shinyFiles - fastq_pass, kraken_db and results===================================
	volumes <- c(Home = fs::path_home(), 
							 wink = fs::path_real("../"), 
							 getVolumes()() )
	
	shinyDirChoose(input, "fastq_pass_folder", 
								 roots = volumes, 
								 session = session)
	
	shinyDirChoose(input, "kraken_db_folder", 
								 roots = volumes, 
								 session = session)
	
	res_folder <- reactive({
		shinyDirChoose(input, "results_folder", 
									 roots = volumes, 
									 session = session)
		return( parseDirPath(volumes, input$results_folder) )
	})
	# and render what has been selected
	output$fastq_folder_selected <- renderText(paste("<b>Selected fastq folder: </b>", parseDirPath(volumes, input$fastq_pass_folder)))
	output$kraken_db_selected <- renderText( paste("<b>Selected kraken db folder: </b>", parseDirPath(roots = volumes, input$kraken_db_folder)) )
	output$results_folder_selected <- renderText(paste( "<b>Selected results folder: </b>", res_folder() ))
	output$db_used <- renderText( paste("<b>kraken db used: </b>", parseDirPath(roots = volumes, input$kraken_db_folder)) )
	
	#-----------------------------------------------------------------------
	stats_folder <- reactive({
		file.path(res_folder(), "latest-stats") 
	})
	
	bracken_folder <- reactive({
		file.path(res_folder(), "latest-bracken")
	})
	#
	# END handling of shinyFiles===================================
	
	headers <- c("file", "format", "type", "num_seqs", 
							 "sum_len", "min_len", "avg_len", "max_len", 
							 "Q1", "Q2", "Q3", "sum_gap", "N50", "Q20_perc", "Q30_perc", 
							 "first_write", "last_write")
	
	# REACTIVES ------------------------------------------------------------------
	# reactive to get actual list of stats files
	
	statsFiles <- reactivePoll(
		1000, session,
		checkFunc = function() {
			if( !is.null(stats_folder() ) ) {
				file.info( stats_folder() )$mtime
			} else {
				""
			}
		},
		valueFunc = function() {
			list.files(stats_folder(), pattern = "*stats.txt", full.names = TRUE)
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
			if( !is.null(bracken_folder()) ) {
				file.info( bracken_folder() )$mtime
			} else {
				""
			}
		}, 
		valueFunc = function() {
			l <- list.files(bracken_folder(), pattern = "*.tsv", full.names = TRUE)
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
			#dplyr::filter(name != "unclassified") %>% # remove these here, they are used for the krakenData reactives
			#mutate(freq = round(freq, 3)) %>%
			group_by(file) %>%
			slice_head(n = input$topn) %>%
			summarize( across("name", .fns = paste, collapse = " | "), 
								 across("kraken_assigned_reads", .fns = sum, na.rm = TRUE), 
								 .groups = "keep" )
	})
	
	# reactive vals for storing total, mapped reads, nxf process info...
	seqData <- reactiveValues(nsamples = 0, treads = 0, tbases = 0, n50 = 0, runtime = 0)
	krakenData <- reactiveValues(all_reads = 0, assigned_reads = 0, unassigned_reads = 0)
	nxf <- reactiveValues(pid = NULL, watch = NULL, status = NULL)
	
	# OBSERVERS------------------------------------------------------------------
	# observer for optional inputs
	hide("optional_inputs")
	observeEvent(input$more, {
		shinyjs::toggle("optional_inputs")
	})
	shinyjs::disable("stop")
	
	observe({
		seqData$nsamples <- nrow( statsData() )
		seqData$treads <- si_fmt( as.character(sum(statsData()$num_seqs, na.rm = TRUE)) )
		seqData$tbases <- si_fmt( as.character(sum(statsData()$sum_len, na.rm = TRUE)) )
		seqData$n50 <-  mean(statsData()$N50, na.rm = TRUE) 
		seqData$runtime <- difftime( max(statsData()$last_write), min(statsData()$first_write) ) # use asPOSIXct here!
		nxflog <- nxf_logfile_data()
		nxf$status <- ifelse(length( nxflog ) == 1, 
												 "Uknown", 
												 # shorten this?
												 gsub(pattern = "\\s+", " ", tail(nxflog, 1)) %>% trimws() %>% strsplit(split = " ")  %>% unlist() %>% tail(1)
		)
	})
	
	
	# RENDERS ------------------------------------------------------------------
	
	# Build parameters for nextflow run =========================================================
	output$stdout <- renderPrint({
		# build nxf call and print to stdout
		# force all three folders to be defined before allowing run
		if ( any(is.integer(input$fastq_pass_folder), is.integer(input$kraken_db_folder), is.integer(input$results_folder)) ) {
			cat("wink command preview, select a fastq_pass folder and a kraken2 database folder to start\n")
			shinyjs::disable("run")
		} else {
			shinyjs::enable("run")
			# hard set fastq folder
			fastq_folder <<- parseDirPath(volumes, input$fastq_pass_folder)
			kraken_db_folder <<- parseDirPath(volumes, input$kraken_db_folder)
			
			# by default results are in the app execution directory
			wink_results_folder <<- case_when( is.integer(input$results_folder) ~ c("--results", here("results-wink")), # 
																				 TRUE ~ c("--results", res_folder() ) 
			)
			
			skip_kraken <<- ifelse(input$skip_kraken, "--skip_kraken", "") # this works because both T and F are length 1, does not work for nxf_profile
			weakmem <<- ifelse(input$weakmem, "--weakmem", "")
			nxf_profile <<- case_when( input$nxf_profile == "local" ~ "", 
																 TRUE ~ c("-profile", input$nxf_profile) )
			
			nxf_args <<- c("run" ,
										 here("main.nf"),
										 "--fastq_pass", fastq_folder,
										 wink_results_folder,
										 skip_kraken,
										 weakmem,
										 "--kraken_db", kraken_db_folder, 
										 "--taxlevel", input$taxlevel,
										 nxf_profile)
			cat("nextflow", nxf_args)
		}
	})
	# Build parameters for nextflow run =========================================================
	
	# CALLS TO NEXTFLOW PIPELINE ----------------------------------------------------------------
	
	# start
	observeEvent(input$run, {
			nxf$pid <- sys::exec_background("nextflow", 
																			args = nxf_args, 
																			std_out = nxf_outfile
			)
			#nx_notify_success(paste("Nextflow pipeline started with pid", nxf$pid))
			nxf$watch <- sys::exec_background(here("bin/watch-pid.sh"), 
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
		
	})
	
	# stop app if the nxf process fails
	observeEvent(nxf$status, {
		if (nxf$status == "Z") {
			nx_report_warning("Nextflow exited!", "If you pressed Stop WINK then this is OK, check the logfiles otherwise. Press Reset to start over.")
			shinyjs::disable("stop")
			shinyjs::toggleCssClass("stop")
			shinyjs::html(selector = ".logo", 
										html = "WINK - What's In my Nanopore reads, with Kraken2, in real-time")
		}
	})
	
	# and kill 
	observeEvent(input$stop, {
		
		if (nxf$pid != 1 && nxf$watch != 1) { # don't kill pideinz please:)
			tools::pskill(nxf$pid)
			#tools::pskill(nxf$watch)
			
			#nx_notify_warning(paste("Nextflow pipeline with pid", nxf$pid, "was stopped!"))
			#shinyjs::enable("run")
			shinyjs::toggleCssClass("stop")
			shinyjs::disable("stop")
			shinyjs::html(selector = ".logo", 
										html = "WINK - What's In my Nanopore reads, with Kraken2, in real-time")
		}
	})
	
	
	# ------------------------------------------------------------------
	
	output$log_output <- renderPrint({
		cat({ nxf_logfile_data() %>% tail(2)}, sep = "\n" )
		#cat(file = stderr(), "nxf_status: ", nxf$status, "\n")
	})
	
	output$nxf_output <- renderPrint({
		cat(nxf_outfile_data(), sep = "\n")
		#runjs("document.getElementById('nxf_output').scrollTo(0,1e9);") # scroll the page to bottom with each message, 1e9 is just a big number
	})
	
	# OUTPUTS FOR STATS TAB--------------------------
	output$stats <- DT::renderDataTable({
		df <- statsData() %>% 
			dplyr::mutate(file = basename(tools::file_path_sans_ext(file)),
										bases = sum_len,
										bases_h = si_fmt( as.character(bases)) ) %>% # issues with bit64 integers? integer64-> as.character -> si_fmt works
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
	
	# OUTPUTS FOR KRAKEN TAB---------------------------
	
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
										 			round(krakenData$assigned_reads/krakenData$all_reads * 100, 0),
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
										 			round(krakenData$unassigned_reads/krakenData$all_reads * 100, 0),
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
									caption = HTML(paste("Top" ,tags$b(input$topn), "kraken-assigned species per sample and # of reads assigned to them (click on a row for more info)")),
									#filter = 'top',
									extensions = 'Buttons', 
									options = list(dom = 'Btp', 
																 buttons = c('copy', 'csv', 'excel')
									), 
									colnames = c("barcode", paste("top", input$topn, "hits"), "assigned reads"),
									rownames = FALSE, class = 'hover row-border')
	})
	
	taxdb_left = "https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id="
	taxdb_right = "&lvl=1&lin=s"
	
	output$ab_table_detail <- DT::renderDataTable({
		# validate ab_table is not empty
		validate(
			need(!is.null(last_selection$row_value), "No data here"),
			need( nrow(brackenData() ) != 0, "No data here" )
		)
		
		# use the row value reactive to filter
		df_assigned <- brackenData() %>% 
			dplyr::filter(file == last_selection$row_value) %>% #name != "unclassified") %>%
			
			dplyr::mutate(
				taxonomyID = paste0(taxdb_left, taxonomy_id, taxdb_right), #two-step mutate to build href
				kraken_reads = kraken_assigned_reads,
				bracken_corr_reads = new_est_reads, 
				bracken_freq = round(freq*100, 2)
			) %>% 
			dplyr::mutate(
				taxonomyID = paste0("<a target = '_blank' href='", taxonomyID, "'>", taxonomy_id, "</a>")
			) %>%
			
			dplyr::select(name, taxonomyID, kraken_reads, bracken_corr_reads, bracken_freq) %>%
			dplyr::filter(bracken_freq >= input$filterFreq)
		
		df_unassigned <- brackenData() %>% 
			dplyr::filter(file == last_selection$row_value, name == "unclassified")
		
		# set kraken statistics for value box
		krakenData$assigned_reads = sum(df_assigned$kraken_reads, na.rm = TRUE)
		krakenData$unassigned_reads = sum(df_unassigned$kraken_assigned_reads, na.rm = TRUE)
		krakenData$all_reads = sum(df_assigned$kraken_reads, df_unassigned$kraken_assigned_reads, na.rm = TRUE)
		
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
			DT::formatStyle('bracken_freq',
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
	
	
	#
	# generate and download rmarkdown report ------------------------------
	#
	output$download_rmarkdown <- downloadHandler(
		
		filename = paste0("wink-report-", format(Sys.time(), "%Y%m%d_%H%M%S"), ".html"),
		content = function(file) {
			# Copy the report file to a temporary directory before processing it, in
			# case we don't have write permissions to the current working dir 
			tempReport <- file.path(tempdir(), "report.Rmd")
			file.copy("report.Rmd", tempReport, overwrite = TRUE)
			
			# Set up parameters to pass to Rmd document
			params <- list(
				n50 = seqData$n50,
				statsData = statsData(),
				brackenData = brackenData() %>% dplyr::filter(freq >= input$filterFreq/100), #filtered data goes in the report, here it is still as fraction!!
				filter_used = input$filterFreq,
				db_used = parseDirPath(roots = volumes, input$kraken_db_folder),
				total_barcodes = seqData$nsamples,
				total_bases = seqData$tbases,
				total_reads = seqData$treads, # sum( statsData()$num_seqs, na.rm = T),
				ass_reads = si_fmt( as.character(sum( brackenData()$kraken_assigned_reads, na.rm = T)) ),
				run_time = paste( round(as.numeric(seqData$runtime, units = 'hours'), digits = 2), "hours" )
			)
			
			# Knit the document, passing in the `params` list, and eval it in a
			# child of the global environment (this isolates the code in the document
			# from the code in this app).
			rmarkdown::render(tempReport, output_file = file,
												params = params,
												envir = new.env(parent = globalenv())
			)
		}
	)
	#
	# generate and download rmarkdown report ------------------------------
	#
	
	# auto-generate report?
	
	observe({
		invalidateLater(90000, session) # 15*60*1000 msec
		reportname <- paste0("wink-autoreport-", format(Sys.time(), "%Y%m%d_%H%M%S"), ".html")
		tempReport <- file.path(tempdir(), "report.Rmd")
		file.copy("report.Rmd", tempReport, overwrite = TRUE)
		
		# Set up parameters to pass to Rmd document
		isolate(
			params <- list(
				n50 = seqData$n50,
				statsData = statsData(),
				brackenData = brackenData() %>% dplyr::filter(freq >= input$filterFreq /100), #filtered data goes in the report, here it is still as fraction!!
				filter_used =  input$filterFreq,
				db_used = parseDirPath(roots = volumes, input$kraken_db_folder),
				total_barcodes = seqData$nsamples,
				total_bases = seqData$tbases,
				total_reads = seqData$treads, # sum( statsData()$num_seqs, na.rm = T),
				ass_reads = si_fmt( as.character(sum( brackenData()$kraken_assigned_reads, na.rm = T)) ),
				run_time = paste( round(as.numeric(seqData$runtime, units = 'hours'), digits = 2), "hours" )
			)
		)
		# generate report only if checkbox and results folder is selected
		if( input$autoreport & !is.integer(input$results_folder) ) {
			rmarkdown::render(tempReport, 
											output_file = reportname,
											output_dir = file.path(getwd(), "auto-reports"),
											params = params,
											envir = new.env(parent = globalenv())
		)
			print( paste("auto report generated on ", Sys.time()) )
			nx_notify_success("Report saved!")
		}
		
		
	})
	# end auto-generate report

	session$onSessionEnded(function() {
		isolate( tools::pskill(nxf$watch) )
		unlink(nxf_logfile)
		unlink(nxf_outfile)
		#system2("rm", args = c("-rf", "work")) # comment out on deploy
	})
	
}
