library(markdown)

navbarPage("Analysis tool",
	tags$head(
		tags$link(rel = "stylesheet", type = "text/css", href = "bootstrap.css")
	),
	tabPanel("Load dada",
		sidebarLayout(
			sidebarPanel(
				#checkboxInput('prev','Last config',FALSE),
				
				#conditionalPanel(condition = "!input.prev",
					fluidRow(
						fileInput('datafile', 'Choose data CSV file',accept=c('text/csv', 'text/comma-separated-values,text/plain')),
						selectInput('clt', 'Select cluster', choices = NULL),
						selectInput('strata', 'Select strata', choices = NULL, multiple = T),
						selectInput('weight', 'Select weights', choices = NULL),
						selectInput('fpop', 'Select population', choices = NULL),
						actionButton("desButton", "Apply",style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
						verbatimTextOutput("design")
					)
				#),
				#conditionalPanel(condition = "input.prev",
				#	fluidRow(
				#		actionButton("prevButton", "Apply",style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
				#		verbatimTextOutput("design")
				#	)
				#)
			),
			mainPanel(
				DT::dataTableOutput("out")
			)
		)
	),
	tabPanel("Plot",
        sidebarLayout(
			sidebarPanel(
				fluidRow(
					column(10,
						h4("DATA")
					),
					column(2,
						actionButton("plotButton", "Plot", style="color: #fff; background-color: #337ab7; border-color: #2e6da4; margin-top: 25px, align: right")
					)
				),
				fluidRow(
					column(12,
					conditionalPanel(condition = "input.smult",
						selectInput('sm_xi', 'Target variable (select mutl)', NULL,multiple=T)
					),
					conditionalPanel(condition = "!input.smult",
						selectInput('xi', 'Target variable', NULL)
					)
					)
				),
				fluidRow(
					column(2,
						h6("Mutl target"),
						checkboxInput('smult','S-mult Q',FALSE)
						#)
					),
					
					column(2,
						selectInput('nomb', 'Numeric?', c("yes","no"), selected="no")
					),
					
					column(2,
						h6("Confidence"),
						checkboxInput("Y_cl",NULL, FALSE)
					),
					column(2,
						conditionalPanel(condition = "input.Y_cl",
							numericInput('CL', 'Level',0.95, min=0.80,	max=0.99, step =0.005)
						)
					)
				),

				fluidRow(
					column(12,
						selectInput('yi', 'Disagregation', NULL, multiple=T)
					)
				),
				fluidRow(
					column(8,
						selectInput('D1filter', 'Filter variable 1', choices = "None", selected="None" ,multiple=F)
					),
					column(4,
							selectInput('T1filter', 'Filter out choices',choices = "None", selected="None" ,multiple=T)
					)
				),
				fluidRow(
					column(8,
						selectInput('D2filter', 'Filter varaible 2', choices = "None", selected="None" ,multiple=F)
					),
					column(4,
							selectInput('T2filter', 'Filter out choices',choices = "None", selected="None" ,multiple=T)
					)
				),
				
				
				### here check to hide the options
				fluidRow(
					column(9,
						h4("GRAPHS option")
					),
					column(3,
							actionButton("logButton", "Log analysis", style="color: #fff; background-color: #F15B55; border-color: #2e6da4")
					)
				),
				fluidRow(
					column(4,
						selectInput('gptyp', 'Type graph', c("Line","Bar chart","Heat map","Pie Chart"),selected="Bar chart")
					),
					column(2,
						conditionalPanel(condition = "input.gptyp=='Pie Chart'",
							numericInput('sdonut', '% donut',0, min=0,	max=100, step =10)
						)
					),
					column(1,
						h6("Flip"),
						checkboxInput("flip", NULL, FALSE)
					),
					column(1,
						h6("Sort"),
						checkboxInput("sort", NULL, FALSE)
					),
					column(2,
						numericInput('hgt', 'Height',400,min=100,	max=3000, step =50)
					),
					column(2,
						numericInput('wdh', 'Width',750,min=100,max=3000, step =50)
					)
				),
				fluidRow(
					column(2,
					     colourpicker::colourInput("H_col", "Color",  palette = "limited",
								allowedCols = c("#F15B55","#F2A0A1","#D3CAB7","#E3E4E5","#C7C8CA","#A7A9AC","#505758","#fff67a","#f69e61","#a5c9a1","#667A95","#89A5C9")
							)
					),
					column(3,
						conditionalPanel(condition = "input.gptyp=='Pie Chart'|input.gptyp=='Bar chart'|input.gptyp=='Line'",
							selectInput('col_ramp', 'Color ramp', c("REACH red","REACH blue","Reds","Purples","Oranges","Greys","Greens","Blues","Spectral"))
						)
					),					
					column(1,
						conditionalPanel(condition = "input.gptyp=='Pie Chart'|input.gptyp=='Bar chart'|input.gptyp=='Line'",
							h6("revert"),
							checkboxInput("col_revert",NULL, FALSE)
						)
					),
					column(3,
						selectInput('ffont', 'Font', font_family, selected="League Gothic")
					),
					column(2,
						numericInput('fsize', 'Font size',6,min=0.5, max=40, step =0.5)					
					)
				),
				fluidRow(
					column(6,
						selectInput('tabxfilter', 'Hide target categories', NULL,multiple=T),
						textInput('laby', 'Target label')
					),
				  
					column(6,
						selectInput('tabyfilter', 'Hide disagregation categories', NULL,multiple=T),
						textInput('labx', 'Disagregation label')
					)
				),width = 5
			),
			mainPanel(
			
			
				plotOutput('plot1',width = "100%"),
				width = 7
			)
		)
	),

	tabPanel("Table",
		textOutput("textest"),
		tableOutput("tableOut1"),
		conditionalPanel(condition = "input.nomb=='yes'",
					tableOutput("statest")
				)
	),
	tabPanel("Analysis Log", 
		sidebarLayout(
			sidebarPanel(
				fluidRow(
					actionButton("refresh", "Refresh", style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
					br(),
					downloadButton("downloadBtn", "Download log"),
					br(),
					actionButton("delete", "Delete older operation", style="color: #fff; background-color: #F15B55; border-color: #2e6da4")
				),
				width = 2
			),
			mainPanel(
				DT::dataTableOutput("test_log")
			)
		)
	),
	tabPanel("Recreate log", 
		sidebarLayout(
			sidebarPanel(
				fluidRow(
					fileInput('logfile', 'Choose CSV data file',accept=c('text/csv', 'text/comma-separated-values,text/plain')),
					br(),
					fileInput('data_update', 'Choose analysis log CSV file',accept=c('text/csv', 'text/comma-separated-values,text/plain')),
					br(),
					downloadButton("report", "Download report"),
					br(),
					downloadButton("report_doc", "Download report word"),
					br()
				),
				width = 3
			),
			mainPanel(
				DT::dataTableOutput("log_new")
			)
		)
	)
)

