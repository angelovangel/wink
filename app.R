#-----------------#
# Main WINK app
#-----------------#
# reactive read of stats and kraken output
# nextflow pipeline has to be started separately in the same folder, i.e. the shiny app reads results-wink
require(shiny)
require(shinydashboard)
require(shinyWidgets)
require(data.table)
require(DT)
require(dplyr)
require(stringr)
#require(prettyunits)
require(shinyjs)

si_fmt <- function(x) { system2("bin/si-format.sh", x, stdout = TRUE) }

ui <- dashboardPage(
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
		tags$head(tags$style(HTML('
        .skin-blue .main-header .logo {
          background-color: green;
        }
        .skin-blue .main-header .logo:hover {
          background-color: green;
        }
      '))),
		tabsetPanel(
			tabPanel("Sequencing overview", #--------------------------------------------------------------
							 fluidRow(
							 	box(width = 12, status = "warning", solidHeader = FALSE, collapsible = TRUE, 
							 			title = "Control panel",
							 		actionBttn("fastq_pass_folder", "Select Run folder", style = 'minimal', color = "primary", icon = icon("folder")),
							 		actionButton("run", "Start WINK", icon = icon("play"))
							 		
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
							 tableOutput("rawtable"),
							 downloadButton("downloadCsv", "Download as CSV")
			), #---------------------------------------------------------------
			tabPanel("Help")
		)
	)
)

server <- function(input, output, session) {
	
	options(shiny.launch.browser = TRUE)
	session$onSessionEnded(function() {
		stopApp() # comment out on deploy
	})
	
	headers <- c("file", "format", "type", "num_seqs", 
							 "sum_len", "min_len", "avg_len", "max_len", 
							 "Q1", "Q2", "Q3", "sum_gap", "N50", "Q20_perc", "Q30_perc", 
							 "first_write", "last_write")
	
	# REACTIVES ------------------------------------------------------------------
	# reactive to get actual list of stats files
	statsFiles <- reactivePoll(
		1000, session,
		checkFunc = function() {
			if( file.exists("results-wink/latest-stats") )
				file.info("results-wink/latest-stats")$mtime
			else
				""
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

	# reactive vals for storing total and mapped reads
	seqData <- reactiveValues(nsamples = 0, treads = 0, tbases = 0, n50 = 0, runtime = 0)
	
	# OBSERVERS------------------------------------------------------------------
	observe({
		seqData$nsamples <- nrow( statsData() )
		seqData$treads <- si_fmt( sum(statsData()$num_seqs, na.rm = TRUE) )
		seqData$tbases <- si_fmt( sum(statsData()$sum_len, na.rm = TRUE) )
		seqData$n50 <-  mean(statsData()$N50, na.rm = TRUE)
		seqData$runtime <- difftime( max(statsData()$last_write), min(statsData()$first_write) ) 
	})
	
	# observeEvent(input$simulate, {
	# 	shinyjs::disable("simulate")
	# 	# system call here, so that I can use wait=F
	# 	p <- system2("./simulateST.R", 
	# 							 args = c("-f" ,"testdata/HMW_Zymo.tsv",
	# 							 				 "-o", "stream.tsv",
	# 							 				 "-s", "1"), 
	# 							 wait = FALSE)
	# 	
	# })
	# 
	
	# RENDERS ------------------------------------------------------------------
	
	output$rawtable <- renderTable({
		statsData()
	})
	#
	output$stats <- renderDataTable({
		df <- statsData() %>% 
			dplyr::mutate(bases = sum_len,
										bases_human = si_fmt(bases)) %>%
			dplyr::select(file, num_seqs, bases, bases_human, max_len, N50, Q20_perc, last_write)
		
			datatable(df, filter = 'top', options = list(dom = 'tp'), rownames = FALSE, class = 'hover row-border') %>%
				DT::formatStyle('num_seqs',
												background = styleColorBar(c(0, max(df$num_seqs)), 'skyblue'),
												backgroundSize = '95% 70%',
												backgroundRepeat = 'no-repeat',
												backgroundPosition = 'right') %>%
				DT::formatStyle('bases',
												background = styleColorBar(c(0, max(df$bases)), 'lightgreen'),
												backgroundSize = '95% 70%',
												backgroundRepeat = 'no-repeat',
												backgroundPosition = 'right') %>%
				DT::formatStyle('max_len',
												background = styleColorBar(c(0, max(df$max_len)), 'lightgreen'),
												backgroundSize = '95% 70%',
												backgroundRepeat = 'no-repeat',
												backgroundPosition = 'right') %>%
				DT::formatStyle('N50',
												background = styleColorBar(c(0, max(df$N50)), 'lightgreen'),
												backgroundSize = '95% 70%',
												backgroundRepeat = 'no-repeat',
												backgroundPosition = 'right') 
	})
	
	output$nsamples <- renderValueBox({
		valueBox(
			value = seqData$nsamples,
			subtitle = "Detected barcodes", 
			color = 'blue'
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
			color = 'olive'
		)
	})
	
	output$n50 <- renderValueBox({
		valueBox(
			value = paste( round(seqData$n50, digits = 0), "bp"), 
			subtitle = "Mean N50",
			color = 'green'
		)
	})
	
	output$runtime <- renderValueBox({
		valueBox(
			value = paste( round(as.numeric(seqData$runtime, units = 'hours'), digits = 2), "hours" ),
			subtitle = "Running time",
			color = 'green'
		)
	})
	# output$treads <- renderValueBox({
	# 	readsData$total <- fileData()[1,3] %>% as.numeric()
	# 	
	# 	valueBox(
	# 		value = prettyNum(readsData$total, big.mark = ","),
	# 		subtitle = HTML("<b>Total reads</b> | <b>Step:</b> ", 
	# 										last(fileData()$step),
	# 										" | <b>Timestamp:</b>",
	# 										as.character(parse_date_time(last(fileData()$time), orders = "a b! d! HMS Y")) 
	# 		),
	# 		icon = icon("dna")
	# 	)
	# })
	# 
	# output$abundanceTable <- renderFormattable({
	# 	if (nrow(filteredData()) == 0)
	# 		return()
	# 	
	# 	df <- filteredData() %>%
	# 		dplyr::filter(step == max(step)) %>%
	# 		dplyr::select(species, prob, sAligned, Abundance_m) %>%
	# 		dplyr::mutate(Abundance_m = formattable::percent(Abundance_m)) %>%
	# 		dplyr::arrange(desc(sAligned)) %>%
	# 		dplyr::mutate(prob = round(prob, 2)) %>%
	# 		# just the top 20?
	# 		head(20)
	# 	
	# 	formattable(df, list(sAligned = formattable::color_bar("lightgreen"),
	# 											 Abundance_m = formattable::color_bar("lightgreen"),
	# 											 prob = formattable::color_bar("lightgreen")))
	# 	
	# })
	# 
	# output$abundancePlot <- renderBubbles({
	# 	if (nrow(filteredData()) == 0)
	# 		return()
	# 	df <- filteredData() %>%
	# 		dplyr::filter(step == max(step)) %>%
	# 		dplyr::arrange(desc(Abundance_m)) %>%
	# 		# assign colors according to prob
	# 		dplyr::mutate(prob_color = scales::col_numeric("YlGn", domain = c(0,1))(prob)) %>%
	# 		# just the top 20?
	# 		head(20)
	# 	
	# 	bubbles(value = df$sAligned, 
	# 					label = df$species, 
	# 					key = df$species, 
	# 					color = df$prob_color,
	# 					tooltip = paste("Prob = ", df$prob, "\n", "Reads = ", df$sAligned))
	# 	
	# })
}
shiny::shinyApp(ui, server)
