function(input, output, session) {
    
   db <- reactive({
		inFile <- input$datafile
		if (is.null(inFile)){
		  return(NULL)
		 }else{
			dt<-read.csv(inFile$datapath,check.names=F)
  	dt<-sanit_data(dt)
		no_disagregation<-rep("all",nrow(dt))
		dt<-cbind(no_disagregation,dt)
		names(dt)[1]<-"No_disagregation"
		dt<-class_assess(dt)
		return(dt)
		}
	})
	
  observe({
  	updateSelectInput(session, "xi", choices = colnames(db()),selected =NULL)
  	updateSelectInput(session, "sm_xi", choices = colnames(db()),selected =NULL)
  	updateSelectInput(session, "yi", choices = colnames(db()),selected ="No_disagregation")
  	updateSelectInput(session, "clt", choices = c("None",colnames(db())),selected="None")
  	updateSelectInput(session, "strata", choices = c("None",colnames(db())),selected="None")
  	updateSelectInput(session, "weight", choices = c("None",colnames(db())),selected="None")
  	updateSelectInput(session, "fpop", choices = c("None",colnames(db())),selected="None")
  })
  
  observe({
    theclass<-ifelse(lapply(db(),class)[which(colnames(db()) %in% input$xi)]%in%c("numeric","decimal"),"yes","no")
    updateSelectInput(session, "nomb", selected=theclass)
  })
  
  #log parameter for faster reloads
  observe({
	  updateSelectInput(session, "D1filter", choices = c("None",colnames(db())),selected="None")
  })
  
  observe({
	  updateSelectInput(session, "D2filter", choices = c("None",colnames(db())),selected="None")
  })
  
  observe({
	  updateSelectInput(session, "T1filter", choices = c("None",levels(as.factor(db()[[input$D1filter]]))),selected="None")
	})
  
  observe({
	  updateSelectInput(session, "T2filter", choices = c("None",levels(as.factor(db()[[input$D2filter]]))),selected="None")
  })
  
  dfr<-eventReactive(input$plotButton,{
  #   
  # if(is.null(input$xi)input$xi==""){
  #   stop(showNotification("target variable missing",type="message"))
  #   
  # }  else if(is.null(input$yi)|input$yi==""){
  #   stop(message("disagregation missing",type="message"))
  # }
    
		if(input$Y_cl){
			CL<-input$CL
		} else {
			CL<- NA
		}
		
		if(input$smult){
			index<-which(names(db()) %in% c(input$sm_xi,input$yi,input$strata,input$clt,input$weight,input$fpop,input$D2filter,input$D1filter))
		} else {
			index<-which(names(db()) %in% c(input$xi,input$yi,input$strata,input$clt,input$weight,input$fpop,input$D2filter,input$D1filter))
		}
		
		de<-db()[,index]
		
		if(input$D1filter!="None"){de<-de[de[[input$D1filter]]%in%input$T1filter,]}
		if(input$D2filter!="None"){de<-de[de[[input$D2filter]]%in%input$T2filter,]}
		
		if(!input$smult){
			de<-na.omit(de)
		}
		
		if(input$smult){
		  de[,input$sm_xi]<- lapply(de[,input$sm_xi],function(x) {car::recode(x,"c('Yes','yes','TRUE',TRUE,'true','1',1)='yes'; c(0,'0','No','no',FALSE,'FALSE','false')='no'")}) %>% data.frame
		  de<-de[de[,input$sm_xi] %>% apply(.,1,function(x){!all(is.na(x))}) %>% which,
		         de %>%  apply(.,2,function(x){!all(is.na(x))}) %>% which]
		  # sm_xi<-sm_xi[sm_xi%in%names(de)]
		}
		
		if(input$smult){
		  de[,input$sm_xi]<-lapply(de[,input$sm_xi],as.factor)
		} else  if(input$nomb=="yes"){
		  if(length(input$xi)>1){de[,input$xi]<-lapply(de[,input$xi],coerc)}else{de[[input$xi]] <- de[[input$xi]] %>% coerc}
		} else{
		  if(length(input$xi)>1){de[,input$xi]<-lapply(de[,input$xi],as.factor)}else{de[[input$xi]] <- de[[input$xi]] %>% as.factor}
		}
		
		if(length(input$yi)>1){de[,input$yi]<-lapply(de[,input$yi],as.factor)}else{de[[input$yi]] <- de[[input$yi]] %>% as.factor}
		

		if(input$clt=="None"){
			clust<-formula(paste0("~",1))
		}else{
			clust<-formula(paste0("~",input$clt))
		}

		if(all(input$strata=="None")){
			strat<-NULL
		}else if(length(input$strata) >1){
			strat<-formula(paste0("~interaction(",paste(input$strata,collapse=","),")"))
		}else{
			strat<-formula(paste0("~",input$strata))
		}
		 
		if(input$weight=="None"){
			wgh=NULL
		}else{
			wgh<-formula(paste0("~",input$weight))
		}
		
		if(input$fpop=="None"){
			fpop=NULL
		}else{
			fpop<-formula(paste0("~",input$fpop))
		}
		
		de<- de %>% type_convert_silent()
		
		dsamp<-svydesign(
		  id=clust,
		  strata=strat,
		  weights=wgh,
		  fpc=fpop, 
		  data=de,
		  nest=TRUE
		)	
		
		# do the aggregation
		if(nrow(de)==0){
		  showModal(modalDialog(
		    title = "Warning",
		    "No valid survey",
		    easyClose = TRUE
		  ))
		  dbgr<-NULL
		} else if(input$smult){
		  dbgr<-process_smultiple(input$sm_xi,input$yi,dsamp,CL,stat.test=F)
		}else  if(input$nomb=="yes" & length(input$yi)==1 ){
		  type<-as.character(input$type)
		  dbgr<-process_num(input$xi,input$yi,dsamp,CL, type=input$type,stat.test=F)
		# }else if(length(input$yi)>1){
		#   dbgr<-process_mutl_disa(input$xi,input$yi,dsamp,CL,input$nomb,stat.test=F)
		}else{
		  dbgr<-process_categories(input$yi,input$xi,dsamp,CL,stat.test=F)
		}
	  dbgr
	 
   })
  
  	
	
	s_design<-eventReactive(input$desButton,{
		if(input$clt=="None"){
			clust<-formula(paste0("~",1))
		}else{
			clust<-formula(paste0("~",input$clt))
		}
		
		if(all(input$strata=="None")){
			strat<-NULL
		}else if(length(input$strata) >1){
			strat<-formula(paste0("~interaction(",paste(input$strata,collapse=","),")"))
		}else{
			strat<-formula(paste0("~",input$strata))
		}
		 
		if(input$weight=="None"){
			wgh=NULL
		}else {
			wgh<-formula(paste0("~",input$weight))
		}
		
		if(input$fpop=="None"){
			fpop=NULL
		}else {
			fpop<-formula(paste0("~",input$fpop))
		}
		 
		 desi<-svydesign(
			id=clust,
			strata=strat,
			weights=wgh,
			fpc=fpop, 
			data=db(),
			nest=TRUE
		)
		return(desi)
	})
	
	
	
	output$design<-renderPrint({s_design()})
	
	labx<-renderText({
		 if(input$nomb=="yes"){
			 if(input$labx==""){
				 labx<-gsub("_"," ",input$xi)
			 }else{
				 labx<-input$labx
			 }
		}else{
			if(input$labx==""){
				labx<-gsub("_"," ",input$yi)
			}else{
				labx<-input$labx
			}
		}
	})

	laby<-renderText({
		if(input$nomb=="yes"){
			if(input$laby==""){
				laby<-gsub("_"," ",input$yi)
			}else{
				laby<-input$laby
			}
		}else{
			if(input$laby==""){
				laby<-gsub("_"," ",input$xi)
			}else{
				laby<-input$laby
			}
		}
	})

  observe({
  	  updateSelectInput(session, "tabxfilter", choices = c("None",levels(as.factor(dfr()[[1]]$rowname))),selected="None")
  	  updateSelectInput(session, "tabyfilter", choices = c("None",levels(as.factor(dfr()[[1]]$variable))),selected="None")
  })
  
	output$tableOut1 <- DT::renderDT({
	  df<-dfr()[[1]]
	  gh <- db()[[input$yi]] %>% table %>% data.frame
	  names(gh)<-c("rowname","valid_n")
	  df<-merge(df,gh,by="rowname",all.x=T)
	  df$type_data<-input$nomb
	  df$aggregation_type<-if(input$nomb=="yes"){if(input$type=="sum"){"sum"}else{"average"}}else{"percentage"}
	  if(all((c("ci_lw","ci_up") %in% names(df)  ))){
	    df$ci_lw<-round(as.numeric(df$ci_lw),2)
	    df$ci_up<-round(as.numeric(df$ci_up),2)
	  }
	  df$value<-round(as.numeric(df$value),2)
    as.data.frame(df)
	})

	# show the table on the right of the loading area
  	output$out <- renderPrint({summary(db())})
	
    output$plot1 <- renderPlot({
      if(!is.null(dfr())){
    		df<-dfr()[[1]][!dfr()[[1]]$rowname%in%input$tabxfilter,]
    		df<-df[!df$variable%in%input$tabyfilter,]
    		df<-sanit_dt(df)
    
    		if(input$smult){
    				spl<-str_split(input$sm_xi[1],"_._")
    				lab<-spl[[1]][length(spl[[1]])-1]
    					df$xi<-rep(lab,nrow(df))
    			} else {
    					df$xi<-rep(input$xi,nrow(df))
    			}
    				
    			df$yi<-rep(input$yi,nrow(df))
    	
    			df<-df[!is.na(df$value),]
    
    		
    		if(!input$smult & input$nomb=="yes"){
    			cdata <- df
    			ordeR<-cdata[with(cdata, order(-value)), ]$rowname
    		}else{	
    			if(length(unique(df$variable))==1){
    				cdata <- df
    				ordeR<-cdata[with(cdata, order(-value, variable)), ]$variable
    			}else{
    				cdata <- ddply(df, "variable", summarise, mean = mean(value,na.rm=TRUE))
    				ordeR<-cdata[with(cdata, order(-mean, variable)), ]$variable
    			}
    		}
    		
    		par(mar = c(0, 4.1, 0, 1))
    		colorsc<-color_ramp(df,input$nomb,input$col_ramp,input$col_revert,input$H_col,input$smult,NULL)
    	
        	 if(input$gptyp=="Line"){
    			    graph_line(df,laby(),input$fsize,input$ffont,input$H_col,input$nomb,input$col_ramp,input$col_revert,input$Y_cl)
        	 }else if(input$gptyp=="Pie Chart"){		 
    			    piechar(df,laby(),input$sdonut,input$fsize,input$ffont,colorsc,input$Y_cl)
        	 }else{ 
        		  graph_crosst_sc(df,labx(),laby(),input$nomb,input$fsize,input$sort,ordeR,input$ffont,input$gptyp,input$flip,input$H_col,input$smult,input$col_ramp,input$col_revert,input$Y_cl)
        	 }
      }
    },height = reactive(input$hgt), width = reactive(input$wdh))
	
		
	# title of data
	output$textest<-renderText({
	  if(!is.null(dfr())){
		  paste0("x= ",labx()," ; y = ",laby()," ; ",dfr()[[3]])
	  }
	})
	
		# output$statest<-renderTable({
		# 	dfr()[[2]]
		# },include.rownames=T)
		# 
		
	logt<-reactive({
		as.data.frame(
			t(sapply(nam, function(x) paste(as.character(input[[x]]),collapse=" -/- ")))
			)
	})
	
	humanTime <- function() format(Sys.time(), "%Y%m%d-%H%M%OS")
		
	saveData <- function(data) {
	  fileName <- sprintf("%s_%s.csv",
						  humanTime(),
						  digest::digest(data))
	  
	  write.csv(x = data, file = file.path(responsesDir, fileName),
				row.names = FALSE, quote = TRUE)
	}

	# action to take when submit button is pressed
	observeEvent(input$logButton, {
	  saveData(logt())
	  showNotification("Operation saved",type="message")
	})
	
	loadData <- function() {
	  files <- list.files(file.path(responsesDir), full.names = TRUE)
	  data <- lapply(files, read.csv, stringsAsFactors = FALSE)
	  data <- do.call(rbind, data)
	  data
	}
	
	output$test_log <- DT::renderDataTable(
		  loadData(),
		  rownames = FALSE,
		  options = list(searching = FALSE, lengthChange = FALSE)
	)
	
	observeEvent(input$refresh, {
			output$test_log<-DT::renderDataTable(
		  loadData(),
		  rownames = FALSE,
		  options = list(searching = FALSE, lengthChange = FALSE)
		  )
	})

	# get the log cleaned
	observeEvent(input$delete, {
		files<-list.files(path = paste0(getwd(),"/",responsesDir), 
			pattern = NULL, all.files = FALSE,
           full.names = FALSE, recursive = FALSE,
           ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
		   
		if(length(files)>0){
			file<-sapply(files,function(x){paste0(responsesDir,"//",x)})
			file.remove(file)
		}
	})
	
	output$downloadBtn <- downloadHandler(
	  filename = function() { 
		sprintf("analysis_plan_%s.csv", humanTime())
	},
	  content = function(file) {
		write.csv(loadData(), file, row.names = FALSE)
	}
	)
	
	# redoing the analysis from the cleaning log
	 log_analysis <- reactive({
		inFile <- input$data_update
		if (is.null(inFile)){
		  return(NULL)
		 }else{
		dt<-read.csv(inFile$datapath)
		return(dt)
		}
	})
	
	db_up <- reactive({
		inLog <- input$logfile
		if (is.null(inLog)){
		  return(NULL)
		 }else{
		dt<-read.csv(inLog$datapath,check.names=F)
		dt<-sanit_data(dt)
		no_disagregation<-rep("all",nrow(dt))
		dt<-cbind(no_disagregation,dt)
		names(dt)[1]<-"No_disagregation"
		dt<-class_assess(dt)
		return(dt)
		}
	})
	
	
	questions <- reactive({
	  inFile <- input$questions
	  if (is.null(inFile)){
	    return(NULL)
	  }else{
	    dt<-read.csv(inFile$datapath)
	    return(dt)
	  }
	})	

	choices <- reactive({
	  inFile <- input$choices
	  if (is.null(inFile)){
	    return(NULL)
	  }else{
	    dt<-read.csv(inFile$datapath)
	    return(dt)
	  }
	})
	
	
  	# output$log_new <- DT::renderDataTable({log_analysis()})
	
	output$log_new<-DT::renderDataTable(
		  as.data.frame(log_analysis()),
		  rownames = FALSE,
		  options = list(searching = FALSE, lengthChange = FALSE)
		  )
	
	
  	#output$log_new <- renderTable({log_analysis()})
	#output$db_new <- renderTable({db_up()})
	
	
	 output$report <- downloadHandler(
      # For PDF output, change this to "report.pdf"
      filename = "analysis_report.html",
      content = function(file) {
        # Copy the report file to a temporary directory before processing it, in
        # case we don't have write permissions to the current working dir (which
        # can happen when deployed).
        # tempReport <- file.path(tempdir(), "report.rmd")
        # file.copy("report.rmd", tempReport, overwrite = TRUE)

        # Set up parameters to pass to rmd document

 params <- list(		
				db_up = db_up(),
				clt = as.character(log_analysis()[["clt"]]),
				strata = strsplit(as.character(log_analysis()[["strata"]])," -/- "),
				weight = as.character(log_analysis()[["weight"]]),
				fpop = as.character(log_analysis()[["fpop"]]),
				xi = as.character(log_analysis()[["xi"]]),
				sm_xi =  strsplit(as.character(log_analysis()[["sm_xi"]])," -/- "),
				T1filter = strsplit(as.character(log_analysis()[["T1filter"]])," -/- "),
				T2filter = strsplit(as.character(log_analysis()[["T2filter"]])," -/- "),
				yi = strsplit(as.character(log_analysis()[["yi"]])," -/- "),
				D1filter =  as.character(log_analysis()[["D1filter"]]),
				D2filter =  as.character(log_analysis()[["D2filter"]]),
				nomb = as.character(log_analysis()[["nomb"]]),
				gptyp = as.character(log_analysis()[["gptyp"]]),
				CL = as.numeric(log_analysis()[["CL"]]),
				Y_cl = as.logical(log_analysis()[["Y_cl"]]),
				fsize = as.numeric(log_analysis()[["fsize"]]),
				ffont = as.character(log_analysis()[["ffont"]]),
				hgt = as.numeric(log_analysis()[["hgt"]]),
				wdh = as.numeric(log_analysis()[["wdh"]]),
				tabxfilter =  strsplit(as.character(log_analysis()[["tabxfilter"]])," -/- "),
				tabyfilter = strsplit(as.character(log_analysis()[["tabyfilter"]])," -/- "),
				labx = as.character(log_analysis()[["labx"]]),
				laby = as.character(log_analysis()[["laby"]]),
				H_col = as.character(log_analysis()[["H_col"]]),
				sdonut = as.numeric(log_analysis()[["sdonut"]]),
				flip = as.logical(log_analysis()[["flip"]]),
				col_ramp = as.character(log_analysis()[["col_ramp"]]),
				smult = as.logical(log_analysis()[["smult"]]),
				col_revert = as.logical(log_analysis()[["col_revert"]]),
				sort = as.logical(log_analysis()[["sort"]])
		)		


        # Knit the document, passing in the `params` list, and eval it in a
        # child of the global environment (this isolates the code in the document
        # from the code in this app).
        rmarkdown::render("report.rmd", output_file = file,
          params = params,
          envir = new.env(parent = globalenv())
        )
      }
    )	
	
	
	
	 output$report_doc <- downloadHandler(
      # For PDF output, change this to "report.pdf"
      filename = "analysis_report.docx",
      content = function(file) {
        # Copy the report file to a temporary directory before processing it, in
        # case we don't have write permissions to the current working dir (which
        # can happen when deployed).
        # tempReport <- file.path(tempdir(), "report_doc.rmd")
        # file.copy("report_doc.rmd", tempReport, overwrite = TRUE)

        # Set up parameters to pass to rmd document
  

      # Set up parameters to pass to rmd document

        # Set up parameters to pass to rmd document

 params <- list(		
				db_up = db_up(),
				clt = as.character(log_analysis()[["clt"]]),
				strata = strsplit(as.character(log_analysis()[["strata"]])," -/- "),
				weight = as.character(log_analysis()[["weight"]]),
				fpop = as.character(log_analysis()[["fpop"]]),
				xi = as.character(log_analysis()[["xi"]]),
				sm_xi =  strsplit(as.character(log_analysis()[["sm_xi"]])," -/- "),
				T1filter = strsplit(as.character(log_analysis()[["T1filter"]])," -/- "),
				T2filter = strsplit(as.character(log_analysis()[["T2filter"]])," -/- "),
				yi = strsplit(as.character(log_analysis()[["yi"]])," -/- "),
				D1filter =  as.character(log_analysis()[["D1filter"]]),
				D2filter =  as.character(log_analysis()[["D2filter"]]),
				nomb = as.character(log_analysis()[["nomb"]]),
				gptyp = as.character(log_analysis()[["gptyp"]]),
				CL = as.numeric(log_analysis()[["CL"]]),
				Y_cl = as.logical(log_analysis()[["Y_cl"]]),
				fsize = as.numeric(log_analysis()[["fsize"]]),
				ffont = as.character(log_analysis()[["ffont"]]),
				hgt = as.numeric(log_analysis()[["hgt"]]),
				wdh = as.numeric(log_analysis()[["wdh"]]),
				tabxfilter =  strsplit(as.character(log_analysis()[["tabxfilter"]])," -/- "),
				tabyfilter = strsplit(as.character(log_analysis()[["tabyfilter"]])," -/- "),
				labx = as.character(log_analysis()[["labx"]]),
				laby = as.character(log_analysis()[["laby"]]),
				H_col = as.character(log_analysis()[["H_col"]]),
				sdonut = as.numeric(log_analysis()[["sdonut"]]),
				flip = as.logical(log_analysis()[["flip"]]),
				col_ramp = as.character(log_analysis()[["col_ramp"]]),
				smult = as.logical(log_analysis()[["smult"]]),
				col_revert = as.logical(log_analysis()[["col_revert"]]),
				sort = as.logical(log_analysis()[["sort"]])
		)		

		
        # Knit the document, passing in the `params` list, and eval it in a
        # child of the global environment (this isolates the code in the document
        # from the code in this app).
        rmarkdown::render("report_doc.rmd", output_file = file,
          params = params,
          envir = new.env(parent = globalenv())
        )
      }
    )
}
