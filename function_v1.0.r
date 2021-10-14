proper<-function(s) sub("(.)", ("\\U\\1"), tolower(s), pe=TRUE)

wrapit <- function(text) {
  wtext <- paste(strwrap(text,width=40),collapse=" \n ")
  return(wtext)
}

coerc<-function(x){x<-as.numeric(as.character(x))}

sanit_dt<-function(x){
	for (j in 1:ncol(x)){
		if(class(x[,j])!="numeric"){
			x[,j]<-iconv(x[,j], to="ASCII//TRANSLIT")
			x[,j]<-gsub("_"," ",x[,j])
			x[,j]<-proper(x[,j])
		}
	}
	return(x)
}

sanit_data<-function(x){
	names(x)<-iconv(names(x), to="ASCII//TRANSLIT")
	names(x)<-gsub("[()'?{}:&!@*%]","",names(x))
	names(x)<-gsub("[,:]"," ",names(x))
	names(x)<-gsub("[/-]","_._",names(x))
	names(x)<-gsub("  "," ",names(x))
	names(x)<-gsub("  "," ",names(x))
	names(x)<-gsub(" ",".",names(x))
	
	for (j in 1:ncol(x)){
			x[,j]<-iconv(x[,j], to="ASCII//TRANSLIT")
			x[,j]<-car::recode(x[,j],"'FALSE'='no';'TRUE'='yes'")# ;TRUE='Yes';FALSE='No'")
			x[,j]<-car::recode(x[,j],"' '=NA;''=NA;'na'=NA;'n/a'=NA;NULL=NA;'NULL'=NA;'FALSE'='no';'TRUE'='yes'")# ;TRUE='Yes';FALSE='No'")
			x[,j]<-gsub("[:,]","_",x[,j])
			x[,j]<-gsub("[()'?{}:&!@*]","",x[,j])
			x[,j]<-gsub("[éèê]","e",x[,j])
			x[,j]<-gsub("[äàâ]","a",x[,j])
			x[,j]<-gsub("[öô]","o",x[,j])
			x[,j]<-tolower(x[,j])
	}
	return(x)
}

sanit_name<-function(x){
	names(x)<-iconv(names(x), to="ASCII//TRANSLIT")
	names(x)<-gsub("[()'?{}:&!@*%]","",names(x))
	names(x)<-gsub("[,:]"," ",names(x))
	names(x)<-gsub("[/-]","_._",names(x))
	names(x)<-gsub("  "," ",names(x))
	names(x)<-gsub("  "," ",names(x))
	names(x)<-gsub(" ",".",names(x))
	
	return(x)
}

class_assess<-function(df){
  require(readr)
  df %>% type_convert
} 


# small function to merge CI and data used into process numbers
merge_cfint<-function(a,dbgr,y){
  
  index<-do.call(paste, c(dbgr[,1:2], sep = ":"))
  dbgr<-cbind(dbgr,index) 
  
  rown<-data.frame(do.call('rbind', strsplit(as.character(row.names(a)),"\\:",fixed=F)))
  rown[,2]<-gsub(paste0("^",y),"",rown[,2])
  cfint<-data.frame(a,rown)
  
  index1<-do.call(paste, c(cfint[,4:3], sep = "::"))
  index2<-do.call(paste, c(cfint[,3:4], sep = "::"))
  if(length(which(index%in%index1))>length(which(index%in%index2))){
    index<-index1	
  } else {
    index<-index2
  }
  
  cfint<-cbind(cfint,index)[,c(1,2,5)]
  names(cfint)[1:2]<-c("ci_lw","ci_up")
  dbgr<-merge(dbgr,cfint,by.x="id",by.y="index")[,-1]
  return(dbgr)
}



# create mean aggregatation and t-test
process_num<-function(x,y,dstrat,CL,type="avg",stat.test=FALSE){
  for (dis in 1:length(x)) 
  {
    f<-formula(paste0("~",y,"+",x[dis]))
    ft<-formula(paste0(x[dis],"~",y))
    fx<-formula(paste0("~",x[dis]))
    fy<-formula(paste0("~",y))
    
    if(type=="avg"){myfuny<-svymean}else if(type=="sum"){myfuny<-svytotal}
    
    mytable<- mytable<-svyby(fx,fy,dstrat,myfuny,na.rm=T)
    mytable<-mytable[,-3]
    
    dstrat_confint<-dstrat
    a<-NA
    try(a<-confint(svyby(fx,fy, dstrat_confint,myfuny,na.rm.by=T,na.rm.all=T,drop.empty.groups=T),level=CL) %>% as.data.frame(),silent=TRUE)
    
    if(!all(is.na(a))){
      names(a)<-c("ci_lw","ci_up")
      dbgr<-cbind(mytable,a)
    } else {
      dbgr<-mytable
    }
    names(dbgr)[2]<-"value"
    
    dbgr<-data.frame(append(dbgr,"variable", after =  1),check.names=FALSE,stringsAsFactors=F)
    
    names(dbgr)[2]<-"variable"
    dbgr$variable<-rep(x[dis],nrow(dbgr))
    
    stest<-NA
    s_dstrat<-NA
    
    lev<-unique(dstrat$variable[[y]])
    sq<-length(lev)
    cr<-as.data.frame(matrix(,nrow=sq,ncol=sq))
    row.names(cr)<-lev
    names(cr)<-lev
    
    # create a matrix with p-values
    
    if(stat.test&(length(x)==1)){	
      for (i in lev)
      {
        for(j in lev)
        {
          if(i!=j){
            s_dstrat<-subset(dstrat,(!is.na(dstrat$variable[[x[dis]]])&!is.na(dstrat$variable[[y]]))&(dstrat$variable[[y]]==i|dstrat$variable[[y]]==j))
            try(stest<-svyttest(ft, design=s_dstrat),silent=TRUE)
            if(!all(is.na(stest))){
              cr[which(row.names(cr)==i),which(names(cr)==j)]<-round(as.numeric(stest$p.value),3)
            }	
          }
        }
      }
    }  else {
      cr<-NA
    }
    
    if (dis==1){
      out<-dbgr
    }else{
      out<-rbind.fill(out,dbgr)
    }
  }
  names(out)[1]<-"rowname"
  dbgr<-out
  
  if(stat.test&(length(x)==1)){		
    tname<-paste0("Pairwise design based t-test; Valid n=", nrow(dstrat$variable)," p-values = ")
  }else{
    tname<-paste0("Valid n=", nrow(dstrat$variable))
  }
  
  return(list(data=dbgr, pvalue = cr , test.name = tname))	
}


process_categories <- function (x,y,dstrat,CL,stat.test=FALSE){

  f<-formula(paste0("~",x,"+",y))
  fx<-formula(paste0("~",x))
  fy<-formula(paste0("~",y))
  tname<-NA;stat<-NA;stdf<-NA;pval<-NA;test<-NA
  
  if(length(na.omit(unique(dstrat$variable[[y]])))==1)
  {
    dbgr<-data.frame(
      rowname=unique(dstrat$variable[[x]]),
      variable=na.omit(unique(dstrat$variable[[y]])),
      value=1
    )
  }else{
    mytable<-svytable(f, dstrat, Ntotal = NULL, round = F)
    
    dbgr<-
      mytable %>% 
      apply(.,1,function(x){x/sum(x)}) %>% 
      melt
    
    dbgr$id<-paste0(dbgr[[x]],"::",dbgr[[y]])
    
    dbcr<- mytable %>% melt
    
    dbcr$id<-paste0(dbcr[[x]],"::",dbcr[[y]])
    dbcr$count<-dbcr$value
    dbcr<-dbcr[,c("id","count")]
    
    dbgr<-merge(dbgr,dbcr,by="id")
    
    if(ncol(dbgr)==1){
      dbgr$variable<-unique(dstrat$variable[[y]])
      dbgr$rowname<-rownames(dbgr)
    }else{
      names(dbgr)[names(dbgr)%in%y]<-c("variable")
      names(dbgr)[names(dbgr)%in%x]<-c("rowname")
    }
    
    # try chi square
    if(stat.test){
      try({
        test<-svychisq(f,dstrat,statistic = "Chisq")
        tname<-test$method
        stat<-round(test$statistic,3)
        stdf<-test$parameter
        pval<-round(test$p.value,3)
      },silent=TRUE)
    }
    
    a<-NA
    try(a<-confint(svyby(fy,fx, dstrat,svymean,na.rm.all=T,na.rm.by=F,drop.empty.groups=F),level=CL),silent=TRUE)
    
    if (ncol(dbgr)==1){
      dbgr<-cbind(dbgr,a)
      dbgr$rownname<-rownames(dbgr)
    } else if(!all(is.na(a))){
      dbgr<-merge_cfint(a,dbgr,y)
    }
  }
  
  dbgr<-dbgr[!is.na(dbgr$value),]
  dbgr[,1]<-factor(dbgr[,1], levels=sort(unique(dbgr[,1]), decreasing=TRUE))
  dbgr[,2]<-factor(dbgr[,2], levels=sort(unique(dbgr[,2]), decreasing=TRUE))
  
  if(stat.test){
    tname<-paste0(tname,", p-value=",pval,"; Valid n: ",nrow(dstrat$allprob),"\n\n\n")
    statis<-pval
  } else {
    tname<-paste0("Valid n: ",nrow(dstrat$allprob),"\n\n\n")
    statis<-NA
  }
  return(list(data=dbgr, pvalue = statis , test.name = tname))	
}

# aggregation select multiple
process_smultiple<-function (x,y,dstrat,CL,stat.test=FALSE){
  for (dis in 1:length(x)) 
  {
    level<-tail(strsplit(x[dis],"_._")[[1]],1)
    
    factlev<-c('1',"yes","TRUE",TRUE)
    output<-process_categories(y,x[dis],dstrat,CL,stat.test=stat.test)
    dbgr<-output$data
    dbgr<-dbgr[dbgr$variable%in%factlev,]
    dbgr$variable<-rep(level,nrow(dbgr))
    
    if (dis==1){
      out<-dbgr
    }else{
      out<-rbind.fill(out,dbgr)
    }
  }
  test="None"
  pval="None"
  dbgr<-out
  return(list(data=dbgr, pvalue = pval , test.name = test))	
}


# aggregation with multiple disaggregation
process_mutl_disa<-function(x,y,dstrat,CL,nb,stat.test=FALSE){
  
  f<-formula(paste0("~",x,"+",paste0("interaction(",paste(y,collapse=","),")")))
  ft<-formula(paste0(x,"~",paste0("interaction(",paste(y,collapse=","),")")))
  fx<-formula(paste0("~",x))
  fy<-formula(paste0("~",paste0("interaction(",paste(y,collapse=","),")")))
  
  if(nb=="no"){
    process_categories(y,x,dstrat,CL,stat.test=stat.test)
  } else if(nb=="yes"){
    process_num(x,y,dstrat,CL,stat.test=stat.test)	
  }
}	





color_ramp<-function(x,nb,col_ramp_typ,col_revert,def_col,smult,var_col){
	
	if(is.null(var_col)){
		if((!smult & nb=="yes")|is.na(col_ramp_typ)|col_ramp_typ=="NA"){
			nb_col<-1
		} else if((smult & nb=="no")|any(names(x)%in%"facet")){
			k<-1
			nb_col<-length(unique(x$rowname))
		}else{
			k<-2
			nb_col<-length(unique(x$variable))
		}
		
	} else {
		nb_col<-length(unique(x[[var_col]]))
	}
	
	if(nb_col==1){
		color<-def_col
	} else {
			
			if(col_ramp_typ == "REACH red"){
				if(nb_col==2){color<-c("#A7A9AC","#D3CAB7")
				} else if (nb_col==1){color<-c("#F15B55")
				} else if (nb_col==3){color<-c("#F15B55","#D3CAB7","#A7A9AC")
				} else if (nb_col==4){color<-c("#F15B55","#D3CAB7","#C7C8CA","#A7A9AC")
				} else if (nb_col==5){color<-c("#F15B55","#F2A0A1","#D3CAB7","#C7C8CA","#A7A9AC")
				} else if (nb_col==6){color<-c("#F15B55","#F2A0A1","#D3CAB7","#E3E4E5","#C7C8CA","#A7A9AC","#505758")
				} else {color<-rep(c("#F15B55","#F2A0A1","#D3CAB7","#E3E4E5","#C7C8CA","#A7A9AC"),50)}
			
			
			} else if(col_ramp_typ == "REACH blue"){
				if(nb_col==2){color<-c("#667A95","#D3CAB7")
				} else if (nb_col==1){color<-c("#667A95")
				} else if (nb_col==3){color<-c("#667A95","#D3CAB7","#A7A9AC")
				} else if (nb_col==4){color<-c("#667A95","#D3CAB7","#C7C8CA","#A7A9AC")
				} else if (nb_col==5){color<-c("#667A95","#89A5C9","#D3CAB7","#C7C8CA","#A7A9AC")
				} else if (nb_col==6){color<-c("#667A95","#89A5C9","#D3CAB7","#E3E4E5","#C7C8CA","#A7A9AC","#505758")
				} else {color<-rep(c("#667A95","#89A5C9","#D3CAB7","#E3E4E5","#C7C8CA","#A7A9AC"),50)}
			
			} else {
				if(nb_col<=8){
					color<-brewer.pal(nb_col,col_ramp_typ)
				}else{
					color<-rep(brewer.pal(8,col_ramp_typ),50)
				}
			}
		
			if(col_revert==TRUE){
				color<-rev(color)
			}
	}
	color
}


hyph<-function(mot,limit){
  nb<-0
  b<-0
  a<-strsplit(as.character(mot),"(\\.|\\_|\\:\\ )")
  for (ct in 1:length(a[[1]])){
    if(ct==1)
    {
      lab<-a[[1]][ct]
    }else{
      if(!is.na(nchar(a[[1]][ct+1]))&b>=limit+nchar(a[[1]][ct+1]))
      {
        lab<-paste0(lab,"\n",a[[1]][ct])
        nb<-nb+limit+2+nchar(a[[1]][ct]) 
      }else{
        lab<-paste0(lab," ",a[[1]][ct])
        nb<-nb+1
      }
    }
    b<-nchar(lab)-nb
  }
  lab
}


graph_crosst_sc<-function(datag,labx,laby,nb,fsize,is_ordered,ordeR,font_family,type_graph,flip,H_col,smult,col_ramp,col_revert,cl_yn){
	
	# here can percentage
	fct<-100
	
	#names(datag)[1]<-"rowname"
	
	if(cl_yn){
		if(nb=="yes"){
			datag$ci_lw<-coerc(datag$ci_lw)
			datag$ci_up<-coerc(datag$ci_up)
		}else{
			datag$ci_lw<-coerc(datag$ci_lw)*fct
			datag$ci_up<-coerc(datag$ci_up)*fct
		}
	}
	
	if(nb=="yes"){
		 datag$value<-coerc(datag$value)
	}else{
		 datag$value<-coerc(datag$value)*fct
	}
		
  colorsc<-color_ramp(datag,nb,col_ramp,col_revert,H_col,smult,NULL)


  if(type_graph=="Heat map"){
		typ<-"heat"
		
		datag$variable<-factor(datag$variable,levels=rev(levels(factor(datag$variable))))
		p <- ggplot(datag, aes(x=rowname, y=variable)) 
		p<-p+geom_tile(aes(fill = value),colour = "white")
		if(cl_yn){
			if(nb=="yes"){
				p<-p+geom_text(aes(label=paste(round(value,2)," (",round(ci_lw,2)," , ",round(ci_up,2),")",sep="")), family=font_family,size = fsize,color="#42423E")
			}else{
				p<-p+geom_text(aes(label=paste(round(value,1),"% (",round(ci_lw,1)," , ",round(ci_up,1),")",sep="")), family=font_family,size = fsize,color="#42423E")
			}
		}else{
			if(nb=="yes"){
				p<-p+geom_text(aes(label=paste(round(value,2),sep="")), family=font_family,size = fsize,color="#42423E")
			} else {
				p<-p+geom_text(aes(label=paste(round(value,1),"%",sep="")), family=font_family,size = fsize,color="#42423E")
			}
		}
		if(is_ordered==TRUE){
			p<-p+scale_y_discrete(limits=rev(ordeR))
		}
	
		p<-p+scale_fill_gradient(low = "white",high = H_col )+
		 scale_x_discrete(expand = c(0, 0)) +
		 xlab("") + 
		 ylab("") 
		
  }else{
		
		if(smult & nb=="No"){
			typ<-"barVmult"
			if(flip){
				datag$rowname<-factor(datag$rowname,levels=rev(levels(factor(datag$rowname))))
			}
			p<-ggplot(data=datag, aes(x=variable, y=value, fill=rowname))
			tp<-labx;labx<-laby;laby<-tp
			p<-p+geom_bar(stat="identity",position = position_dodge(width = 0.9), alpha=.9)
			p<-p+scale_x_discrete(limits=rev(ordeR))
			
		} else if(any(names(datag)%in%"facet")){
			typ<-"facet_barV"
			datag$rowname<-as.character(datag$rowname)
			datag$facet<-as.character(datag$facet)
			datag$variable<-as.character(datag$variable)
			if(flip){
				datag$variable<-factor(datag$variable,levels=rev(levels(factor(datag$variable))))
			}
			
			tp<-labx;labx<-laby;laby<-tp
			
			if(length(unique(datag$rowname))==1){
				p<-ggplot(data=datag, aes(x=variable, y=value, fill=facet))#,group=facet))
				p<-p+geom_bar(stat="identity",position = position_dodge(width = 0.9), alpha=.9)
			} else {
				p<-ggplot(data=datag, aes(x=variable, y=value, fill=rowname))#,group=facet))
				p<-p+geom_bar(stat="identity",position = position_dodge(width = 0.9), alpha=.9)
			}
			
			if(cl_yn){
					p<-p+geom_errorbar(aes(ymin=ci_lw, ymax=ci_up),position = position_dodge(width = 0.9), alpha=.7,width=.5,color="#42423E")
			}
			
			if(length(unique(datag$rowname))!=1){
				p<-p+facet_wrap(~ facet)
			}
			
		} else if(!smult & nb=="yes" & length(unique(datag$variable))==1 ){
			typ<-"int"
			if(flip){
				datag$rowname<-factor(datag$rowname,levels=rev(levels(factor(datag$rowname))))
			}
			p<-ggplot(data=datag, aes(x=rowname, y=value))
			p<-p+geom_bar(fill=colorsc,stat="identity",position = position_dodge(width = 0.9), alpha=.9)
			
		}else if(!smult & nb=="yes" & length(unique(datag$variable))>1 ){
			typ<-"barHy"
			 colorsc<-color_ramp(datag,"no",col_ramp,col_revert,H_col,smult,NULL)
			if(flip){
				datag$rowname<-factor(datag$rowname,levels=rev(levels(factor(datag$rowname))))
			}
			p<-ggplot(data=datag, aes(x=rowname, y=value, fill=variable))
			p<-p+geom_bar(stat="identity",position = position_dodge(width = 0.9), alpha=.9)
				
		}else if(!smult & !(any(levels(as.factor(datag[,2]))=="Oui")) & !(any(levels(as.factor(datag[,2]))=="yes")) & !(any(levels(as.factor(datag[,2]))=="oui")) ){
			typ<-"barHy"
			if(flip){
				datag$rowname<-factor(datag$rowname,levels=rev(levels(factor(datag$rowname))))
				datag$variable<-factor(datag$variable,levels=rev(levels(factor(datag$variable))))
			}
			p<-ggplot(data=datag, aes(x=rowname, y=value, fill=variable))
			p<-p+geom_bar(stat="identity",position = position_dodge(width = 0.9), alpha=.9)
				
		}else{
				typ<-"barV"
				p<-ggplot(data=datag, aes(x=variable, y=value, fill=rowname))
				tp<-labx;labx<-laby;laby<-tp
				p<-p+geom_bar(stat="identity",position = position_dodge(width = 0.9), alpha=.9)
				if(is_ordered==TRUE){
					p<-p+scale_x_discrete(limits=rev(ordeR))
				}else{
					p<-p+scale_x_discrete()
				}
		}
	
	
		if(cl_yn & typ!="facet_barV"){
			p<-p+geom_errorbar(aes(ymin=ci_lw, ymax=ci_up),position = position_dodge(width = 0.9), alpha=.7,width=.5,color="#42423E")
		}
		
		# if(!smult & nb=="yes" & length(unique(datag$variable))==1 ){
			# p<-p+geom_errorbar(aes(ymin=datag$ci_lw, ymax=datag$ci_up),position = position_dodge(width = 0.9), alpha=.7,width=.5,color="#42423E")
		# }
		
		if(typ=="facet_barV"){
			p<-p+geom_text(aes(label = paste0(formatC(round(value,1), digits = 1, format = "f"),"%")),position = position_dodge(width = 0.75),family=font_family,size = round(fsize*1.2,0),color="#42423E")+
			scale_fill_manual(guide = guide_legend(reverse=TRUE),hyph(laby,10),values=rev(colorsc))			
		}else if(typ=="int"){
			p<-p+geom_text(aes(label = formatC(round(value,2),digits = 2,format = "f")),position = position_dodge(width = 0.75),family=font_family,size = round(fsize*1.2,0),color="#42423E")
		}else if(typ=="barV"){
			if(nb=="yes"){
				p<-p+geom_text(aes(label = formatC(round(value,2), digits = 2, format = "f")),position = position_dodge(width = 0.75),family=font_family,size = round(fsize*1.2,0),color="#42423E")+
				scale_fill_manual(guide = guide_legend(reverse=FALSE),hyph(laby,10),values=rev(colorsc))
			}else{
				p<-p+geom_text(aes(label = paste0(formatC(round(value,1), digits = 1, format = "f"),"%")),position = position_dodge(width = 0.75),family=font_family,size = round(fsize*1.2,0),color="#42423E")+
				scale_fill_manual(guide = guide_legend(reverse=FALSE),hyph(laby,10),values=rev(colorsc))
			}
		}else{
			if(nb=="yes"){
				p<-p+geom_text(aes(label = formatC(round(value,2), digits = 2, format = "f")),position = position_dodge(width = 0.75),family=font_family,size = round(fsize*1.2,0),color="#42423E")+
				scale_fill_manual(guide = guide_legend(reverse=FALSE),hyph(laby,10),values=rev(colorsc))
			}else{
				p<-p+geom_text(aes(label = paste0(formatC(round(value,1), digits = 1, format = "f"),"%")),position = position_dodge(width = 0.75),family=font_family,size = round(fsize*1.2,0),color="#42423E")+
				scale_fill_manual(guide = guide_legend(reverse=FALSE),hyph(laby,10),values=rev(colorsc))
			}
		}
		
		if(is_ordered==TRUE){
			p<-p+scale_x_discrete(limits=rev(ordeR))
			p<-p+scale_y_continuous(breaks=NULL)
		}else{
			p<-p+scale_y_continuous(breaks=NULL)
		}
		
  }
  if(typ=="barVmult"){
    p<-p+theme_pander()+
		theme(
			legend.position = c(0.8,0.2),
			legend.key.width = unit(2, "cm"),
			legend.key.height = unit(1, "cm"),
			legend.title=element_text(family=font_family,size=round(fsize*3.4,0),face="plain",color="#42423E"),
            legend.text=element_text(family=font_family,size=round(fsize*3.4,0),color="#42423E"),
			axis.text.x=element_text(family=font_family,size=round(fsize*3.4,0),color="#42423E"),
            axis.text.y=element_text(family=font_family,size=round(fsize*3.4,0),color="#42423E"),
			axis.title.x = element_blank(),
			axis.title.y = element_blank()
      )+labs(x = hyph(labx,10))
	  
	}else if(typ=="facet_barV"){
    p<-p+theme_pander()+
		theme(
			legend.key.width = unit(2, "cm"),
			legend.key.height = unit(1, "cm"),
			legend.title=element_text(family=font_family,size=round(fsize*3.4,0),face="plain",color="#42423E"),
            legend.text=element_text(family=font_family,size=round(fsize*3.4,0),color="#42423E"),
			axis.text.x=element_text(family=font_family,size=round(fsize*3.4,0),color="#42423E"),
            axis.text.y=element_text(family=font_family,size=round(fsize*3.4,0),color="#42423E"),
			axis.title.x = element_blank(),axis.title.y = element_blank(),
			strip.text.x = element_text(family=font_family,size=round(fsize*3.4,0),face="plain",color="#42423E"),
			strip.text.y = element_text(family=font_family,size=round(fsize*3.4,0),face="plain",color="#42423E")
    )+labs(x = hyph(labx,10))
	  
	}else if(typ!="heat"){
    p<-p+theme_pander()+
		theme(
			legend.key.width = unit(2, "cm"),
			legend.key.height = unit(1, "cm"),
			legend.title=element_text(family=font_family,size=round(fsize*3.4,0),face="plain",color="#42423E"),
            legend.text=element_text(family=font_family,size=round(fsize*3.4,0),color="#42423E"),
			axis.text.x=element_text(family=font_family,size=round(fsize*3.4,0),color="#42423E"),
            axis.text.y=element_text(family=font_family,size=round(fsize*3.4,0),color="#42423E"),
			axis.title.x = element_blank(),axis.title.y = element_blank()
      )+labs(x = hyph(labx,10))
	  
  }else{
    p.bot<-p+ theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA,color="white", size=0.5, linetype="solid"),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_rect(fill="white"),
        plot.background = element_rect(fill="white"),
        legend.position = "none", 
        axis.text.x = element_blank(),
        plot.margin = unit(c(1,0,0,0), "cm"),
        axis.text.y=element_text(family=font_family,size=round(fsize*3.4,0),color="#42423E"))

	p.top   <-  p.bot + theme(
		axis.text.x = element_text(family=font_family,size=round(fsize*3.4,0)),
		axis.text.y = element_text(family=font_family,size=round(fsize*3.4,0),color="white")
		)  + coord_cartesian(ylim = c(0,0))
		
	require(gtable)
	#Extract Grobs
	g1<-ggplotGrob(p.top)
	g2<-ggplotGrob(p.bot)
	#Bind the tables
	g<-gtable:::rbind_gtable(g1, g2, "first")
	#Remove a row between the plots
	g <- gtable_add_rows(g, unit(-1.25,"cm"), pos=nrow(g1))
	#draw
	panels <- g$layout$t[grep("panel", g$layout$name)]
	g$heights[panels] <- unit.c(unit(0,"null"),unit(2,"null"))
	grid.newpage()
	grid.draw(g)
	
  
  }
  if(typ!="heat"){
	if(flip){
		p<-p+coord_flip()
	}
  }
  if(typ!="heat"){
   print(p)
  }
}


graph_line<-function(datag,laby,fsize,font_family,H_col,nb,col_ramp,col_revert,cl_yn){
	
	
	#names(datag)[1]<-"rowname"
	
	if(cl_yn){
		#names(datag)[4:5]<-c("ci_lw","ci_up")
		datag$ci_lw<-as.numeric(as.character(datag$ci_lw))
		datag$ci_up<-as.numeric(as.character(datag$ci_up))
	}
	
		datag$value<-as.numeric(as.character(datag$value))
		
	
	if(nb=="yes"){
	
		colorsc<-color_ramp(datag,nb,col_ramp,col_revert,H_col,FALSE,NULL)
		datag$rowname<-factor(datag$rowname,levels=levels(factor(datag$rowname)))
		p<-ggplot(data=datag, aes(x=rowname, y=value,group=1))
		p<-p+geom_ribbon(aes(ymin = ci_lw, ymax = ci_up), fill = H_col,alpha=0.05) 
		p<-p+geom_path(color=H_col,stat="identity",position = position_dodge(width = 0.9), alpha=.9,size=1.5)	
			
		p<-p+geom_errorbar(aes(ymin=ci_lw, ymax=ci_up),position = position_dodge(width = 0.9), alpha=.7,width=.5,color="#42423E")
		p<-p+geom_text(aes(label = paste0(formatC(round(value,2), digits = 1, format = "f"))),position = position_dodge(width = 0.75),family=font_family,size = round(fsize*0.8,0),color="#42423E")+
				scale_fill_manual(guide = guide_legend(reverse=TRUE),hyph(laby,10),values=colorsc)+
				scale_y_continuous(breaks=NULL)
				
		
		p<-p+theme_pander()+
		theme(
			 legend.key.width = unit(2, "cm"),
			 legend.key.height = unit(1, "cm"),
			 legend.title=element_text(family=font_family,size=round(fsize*3.4,0),face="plain",color="#42423E"),
			 legend.text=element_text(family=font_family,size=round(fsize*3.4,0),color="#42423E"),
			 axis.text.x=element_text(family=font_family,size=round(fsize*3,0),color="#42423E"),
			 axis.text.y=element_text(family=font_family,size=round(fsize*3,0),color="#42423E"),
			 axis.title.x = element_blank(),
			 axis.title.y = element_blank()
		)
				
		p<-p+labs(x = hyph(laby,30))
	
	} else {
	
		

		if(any(names(datag)%in%"facet")){
			typ<-"facet_barV"
			datag<-class_assess(datag)
			datag$variable<-as.factor(as.character(datag$variable))
			datag$facet<-as.factor(as.character(datag$facet))
			
		
			if(length(unique(datag$rowname))!=1){
				p<-ggplot(data=datag, aes(x=rowname, y=value, colour=variable, group=variable))
				colorsc<-color_ramp(datag,"no",col_ramp,col_revert,H_col,TRUE,"variable")
			} else {
				p<-ggplot(data=datag, aes(x=variable, y=value, colour=facet,  group=facet))
				colorsc<-color_ramp(datag,"no",col_ramp,col_revert,H_col,TRUE,"facet")
			}
			
			p<-p+geom_line(size=0.8)
						
			if(cl_yn){
				p<-p+geom_ribbon(aes(ymin = ci_lw, ymax = ci_up) ,alpha=0.1,fill = "grey70",colour="white") 
			}
			
			p<-p+scale_color_manual(guide = guide_legend(reverse=FALSE),hyph(laby,10),values=colorsc)
			
			if(length(unique(datag$rowname))!=1){
				p<-p+facet_wrap(~ facet)
			}
			
		} else {
		
			datag$variable<-as.numeric(as.character(datag$variable))
			colorsc<-color_ramp(datag,nb,col_ramp,col_revert,H_col,FALSE,"rowname")
		
			p<-ggplot(data=datag, aes(x=variable, y=value, colour=rowname,group=rowname))
			p<-p+geom_line(size=0.8)
			if(cl_yn){
				p<-p+geom_ribbon(aes(ymin = ci_lw, ymax = ci_up) ,alpha=0.1,fill = "grey70",colour="white") 
			}
			
			p<-p+scale_color_manual(guide = guide_legend(reverse=FALSE),hyph(laby,10),values=colorsc)
		}
		
		p<-p+theme_bw()+
		theme(
			 legend.key.width = unit(2, "cm"),
			 legend.key.height = unit(1, "cm"),
			 legend.title=element_text(family=font_family,size=round(fsize*3.4,0),face="plain",color="#42423E"),
			 legend.text=element_text(family=font_family,size=round(fsize*3.4,0),color="#42423E"),
			 axis.text.x=element_text(family=font_family,size=round(fsize*3,0),color="#42423E"),
			 axis.text.y=element_text(family=font_family,size=round(fsize*3,0),color="#42423E"),
			 panel.border = element_blank(),
			 panel.grid.major.x = element_blank() ,
			 # panel.grid.minor = element_blank(),
			 #panel.grid.major = element_blank(),
			 axis.title.x = element_blank(),
			 axis.title.y = element_blank()
		)
	
	}
	
	  
	print(p)
	  
}


piechar<-function(datag,laby,sdonut,fsize,font_family,colorsc,cl_yn){
	
	# Add addition columns, needed for drawing with geom_rect.
	#names(datag)[1]<-"rowname"
	datag = datag[order(datag$rowname), ]
	tf_cl<-FALSE
	
	if(cl_yn){
		# names(datag)[4:5]<-c("ci_lw","ci_up")
		tf_cl<-TRUE
	}
	
	
	
	datag$fraction = do.call(c,tapply(datag$value, datag$rowname, function(x) {x / sum(x)},simplify = TRUE))
	
	datag$ymax = do.call(c,tapply(datag$fraction, datag$rowname, function(x) {cumsum(x)},simplify = TRUE))
	datag$ymin = do.call(c,tapply(datag$ymax, datag$rowname, function(x) {c(0, head(x, n=-1))},simplify = TRUE))
	datag$pos = do.call(c,tapply(datag$fraction, datag$rowname, function(x) {cumsum(x) -  x/2},simplify = TRUE))
	datag<-as.data.frame(datag)
	

	
	p = ggplot(datag, aes(fill=as.factor(variable), ymax=ymax, ymin=ymin, xmax=100, xmin=sdonut)) 
	p = p + geom_rect(colour="white",size=0.5) 
	
	p<-p+ facet_wrap( ~ rowname, ncol=2)
	
	if(tf_cl==TRUE){
		p = p + geom_text(aes(x=mean(c(100,sdonut)), y=pos, label = paste(round(fraction*100,1),"% (",round(ci_lw,1)," , ",round(ci_up,1),")",sep="")), size=round(fsize*0.8,0), family=font_family) 
	}else{
		p = p + geom_text(aes(x=mean(c(100,sdonut)), y=pos, label = paste(round(fraction*100,1),"%",sep="")), size=round(fsize*0.8,0), family=font_family) 
	}
	p = p + coord_polar(theta="y") 
	p = p + xlim(c(0, 100)) 
	p = p + scale_fill_manual(guide = guide_legend(reverse=FALSE),hyph(laby,20),values=rev(colorsc)) 
	p<-p+theme_pander()+
	theme(
		 legend.key.width = unit(2, "cm"),
		 legend.key.height = unit(1, "cm"),
		 legend.title=element_text(family=font_family,size=round(fsize*2,0),face="plain",color="#42423E"),
         legend.text=element_text(family=font_family,size=round(fsize*2,0),color="#42423E"),
		 axis.text.x=element_text(family=font_family,size=round(fsize*2,0),color="#808080"),
         axis.text.y=element_text(family=font_family,size=round(fsize*2,0),color="#42423E"),
		 axis.title.x = element_blank(),
		 panel.grid=element_blank(),
		 axis.text=element_blank(),
		 axis.ticks=element_blank(),
		 axis.title.y = element_blank()
    )
	print(p)
	
}
