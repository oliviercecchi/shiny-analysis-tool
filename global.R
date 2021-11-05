require(devtools)
require(plyr)
require(dplyr)
require(car)
require(reshape2)
require(stringr)
require(survey)
require(ggplot2)
require(ggthemes)
require(pander)
require(knitr)
require(extrafont)
require(shiny)
require(shinyjs)
require(gtools)
require(RCurl)
require(sysfonts)
require(showtext)
require(shinythemes)
require(pander)
require(RColorBrewer)
require(markdown)
require(colourpicker)
require(DT)
require(readr)


nam<-c(				
	"datafile",
	"clt",
	"strata",
	"weight",
	"fpop",
	"xi",
	"sm_xi",
	"yi",
	"D1filter",
	"T1filter",
	"D2filter",
	"T2filter",
	"nomb",
	"gptyp",
	"CL",
	"Y_cl",
	"fsize",
	"ffont",
	"hgt",
	"wdh",
	"tabxfilter",
	"tabyfilter",
	"labx",
	"laby",
	"H_col",
	"sdonut",
	"flip",
	"col_ramp",
	"col_revert",
	"smult",
	"sort"
)					


logf<-read.table(text = "", col.names = nam)			

responsesDir<-"analysis_log"

fonttable<-fonttable()
font_family<-mixedsort(unique(as.character(fonttable[["FamilyName"]])))

source("function_v1.0.r")
options(shiny.maxRequestSize=1000*1024^2)
options("survey.lonely.psu"="adjust")





