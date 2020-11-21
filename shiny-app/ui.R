ui <- dashboardPage(
	title = "WINK",
	dashboardHeader(
		title = "WINK - What's In my Nanopore reads, with Kraken2, in real-time",
		titleWidth = '100%',
		dropdownMenuOutput("notificationsMenu")
	),
	dashboardSidebar(disable = TRUE),
	
	dashboardBody(
		tags$head(tags$style(
			HTML("
			.yellow {
				background-color: #F4D03F;
				}
			.red {
			color: red;
			}
						 ")
		)),
		#includeCSS("custom.css"),
		useShinyjs(),
		#useShinyalert(),
		use_notiflix_notify(position = "right-top", width = "380px"),
		
		tabsetPanel(
			tabPanel(
				"Sequencing overview",
				#--------------------------------------------------------------
				fluidRow(
					box(
						width = 12,
						status = "warning",
						solidHeader = FALSE,
						collapsible = TRUE,
						title = "Control panel",
						
						shinyDirButton(
							id = "fastq_pass_folder",
							label = "Select fastq_pass folder",
							title = "Select fastq_pass folder",
							style = "color: #3498DB;",
							#color = "primary",
							icon = icon("folder-open")
						),
						shinyDirButton(
							id = "kraken_db_folder",
							label = "Select kraken database folder", 
							title = "Select kraken database folder",
							style = "color: #3498DB;",
							icon = icon("folder-open")
						),
						actionButton(
							"reset",
							"Reset",
							style = "color: #3498DB;",
							icon = icon("sync"),
							onclick = "history.go(0);"
						),
						actionButton("run", "Start WINK", style = "color: #3498DB;", icon = icon("play")),
						actionButton(
							"stop",
							"Stop WINK",
							style = "color: #3498DB;",
							icon = icon("power-off")
						),
						actionButton(
							"more",
							"More options",
							style = "color: #3498DB;",
							icon = icon("cog"),
							class = "rightAlign"
						),
						tags$div(
							id = "optional_inputs",
							column(
								width = 12,
								#checkboxInput("testrun", "Simulate run with test data", width = "100%"),
								checkboxInput("skip_kraken", "Skip kraken2, show only run statistics", width = "100%"),
								checkboxInput(
									"weakmem",
									"Do not load kraken2 database in RAM (use on weak machines)"
								),
								shinyDirButton(
									id = "results_folder",
									label = "Select results folder",
									title = "Select results folder, default is results-wink",
									#style = "color: #3498DB;",
									#color = "primary",
									icon = icon("folder-open")
								), 
								tags$hr(),
							),
							column(
								width = 4,
								selectizeInput(
									"nxf_profile",
									"nextflow profile",
									width = "100%",
									choices = c("docker", "conda", "local"),
									selected = "docker"
								)
							),
							column(
								width = 4,
								selectizeInput(
									"taxlevel",
									"Taxonomic level for kraken2 abundance",
									width = "100%",
									choices = c(
										"Domain" = "D",
										"Phylum" = "P",
										"Class" = "C",
										"Order" = "O",
										"Family" = "F",
										"Genus" = "G",
										"Species" = "S"
									),
									selected = "S"
								)
							)
						),
						tags$hr(),
						column(
							widt = 12,
							verbatimTextOutput("stdout"),
							verbatimTextOutput("log_output")
						)
					)
				),
				fluidRow(
					box(
						width = 12,
						status = "warning",
						solidHeader = FALSE,
						collapsible = TRUE,
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
						width = 12,
						status = "warning",
						solidHeader = FALSE,
						collapsible = TRUE,
						collapsed = FALSE,
						title = "Run statistics per sample",
						DT::dataTableOutput("stats", width = "100%", height = 500)
					)
				)
			),
			#---------------------------------------------------------------
			tabPanel("Taxonomy and abundance", #-----------------------------
							 fluidRow(
							 	box(
							 		width = 12,
							 		status = "warning",
							 		solidHeader = FALSE,
							 		collapsible = TRUE,
							 		title = "Kraken read assignment per sample",
							 		infoBoxOutput("current_barcode", width = 3),
							 		infoBoxOutput("all_reads", width = 3),
							 		infoBoxOutput("ass_reads", width = 3),
							 		infoBoxOutput("unass_reads", width = 3),
							 		downloadButton("download_rmarkdown", "Save report", style = "color: #3498DB;"),
							 		tags$a("Only data filtered by abundance will be exported")
							 	)
							 ), 
							 # 
							 textOutput(outputId = "db_used"),
							 
							 fluidRow(
							 	box(
							 		width = 6,
							 		status = "warning",
							 		solidHeader = FALSE,
							 		collapsible = TRUE,
							 		collapsed = FALSE,
							 		title = "Bracken results per sample",
							 		sliderInput(
							 			"topn",
							 			"Top N species to show",
							 			value = 3,
							 			min = 1,
							 			max = 5,
							 			step = 1
							 		),
							 		DT::dataTableOutput("ab_table", width = "100%", height = 500)
							 	),
							 	box(
							 		width = 6,
							 		status = "warning",
							 		solidHeader = FALSE,
							 		collapsible = TRUE,
							 		collapsed = FALSE,
							 		title = "Detailed bracken results",
							 		shinyWidgets::sliderTextInput(
							 			"filterFreq",
							 			"Filter by abundance",
							 			choices = c(0, 0.1, 1, 5, 10, 20),
							 			selected = 0,
							 			post = "%",
							 			grid = TRUE
							 		),
							 		#sliderInput("filterFreq", "Filter by abundance", min = 0, max = 1, value = 0, step = 0.1),
							 		DT::dataTableOutput("ab_table_detail", width = "100%", height = 500)
							 	)
							 )),
			#---------------------------------------------------------------
			tabPanel("Nextflow output",
							 verbatimTextOutput("nxf_output")
							 ),
			tabPanel("Help", includeMarkdown(here( "README.md" ))
							 )
		)
	)
)
