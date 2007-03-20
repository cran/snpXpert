
###################################
# main function
###################################
snpXpert<-function()
{  
     library(tcltk)
#     library(dgc.genetics)
#     library(genetics)
     library(gap)
     library(tcltk2)
     library(tkrplot)
     library(ORMDR)
     tclRequire("Iwidgets")

     Result.tt <<- tktoplevel()
     tkwm.title(Result.tt,"BIBS SNP Analyzer")

     topMenu<-tkmenu(Result.tt)
     tkconfigure(Result.tt,menu=topMenu)
     File.Menu <- tkmenu(topMenu,tearoff=FALSE)
     Init.Menu <- tkmenu(topMenu,tearoff=FALSE)
     Anal.Menu <- tkmenu(topMenu,tearoff=FALSE)

     tkadd(File.Menu,"command",label="Read Data",command=function() readF())
     tkadd(File.Menu,"command",label="Save Result",command=function() print("Save Result"))
     tkadd(File.Menu,"command",label="Quit",command=function() tkdestroy(Result.tt))

     tkadd(Init.Menu,"command",label="Set Class variable",command=function() selectClass())
     tkadd(Init.Menu,"command",label="Set SNPs",command=function() selectSNPs())
     tkadd(Init.Menu,"command",label="Check Genotype",command=function() makeForm.data(SNPlist.id))
     tkadd(Init.Menu,"command",label="Set Case/Control",command=function() CaseControl.set())

     tkadd(Anal.Menu,"command",label="HWE",command=function() HWE(SNPlist.id,Class.data))
     tkadd(Anal.Menu,"command",label="allele X2",command=function() chiX2.allele.test(SNPlist.id,Class.data))
     tkadd(Anal.Menu,"command",label="genotype X2",command=function() chiX2.genotype.test(SNPlist.id,Class.data))
#     tkadd(Anal.Menu,"command",label="Cochran-Armitage test",command=function() CA.test(SNPlist.id,Class.data))
     tkadd(Anal.Menu,"command",label="Regression",command=function() regression.test(SNPlist.id,Class.data))
     tkadd(Anal.Menu,"command",label="Logistic Regression",command=function() logistic.test(SNPlist.id,Class.data))
     tkadd(Anal.Menu,"command",label="LD",command=function() LD.plot(SNPlist.id,Class.data))
     tkadd(Anal.Menu,"command",label="HaploScore",command=function() haploScore(SNPlist.id,Class.data))
     tkadd(Anal.Menu,"command",label="MDR",command=function() MDR(SNPlist.id,Class.data))

     tkadd(topMenu,"cascade",label="File",menu=File.Menu)
     tkadd(topMenu,"cascade",label="Initiation",menu=Init.Menu)
     tkadd(topMenu,"cascade",label="Analysis",menu=Anal.Menu)
 
     tkfocus(Result.tt)

     ###############
     tkpack(tn <<- tkwidget(Result.tt, "iwidgets::tabnotebook"))
     tkconfigure(tn,                         # prettyfication taken from incrtcl
            tabpos="n",
            width=1000,
            height=350,
            angle=0,
            bevelamount=2,
            gap=2,
            margin=2,
            tabborders=0,
            tabbackground="white",
            background="lightgray",
            backdrop="lightgray")



###################################
# read data file(csv file) and put it into RAW.data 
###################################

readF<-function()
{    fileName<-tclvalue(tkgetOpenFile())
      if (!nchar(fileName))
     {    tkmessageBox(message="No file was selected!")
     } else
     {   RAW.data<<-read.csv(fileName,header=T)
     }
}


###################################
# set Class variable
###################################

selectClass<-function()
{
   list.name <<- colnames(RAW.data)
   rb <- NULL
   tt <- tktoplevel()
   rbValue <- tclVar("0")
   tkgrid(tklabel(tt,text="Select Class variable!"))
   for (i in 1:length(list.name) ) {
        rbtemp <- tkradiobutton(tt)
        tkconfigure(rbtemp, variable=rbValue, value=list.name[i])
        tkgrid(tklabel(tt,text=list.name[i]),rbtemp)
   }

   OnOK <- function()
   {
       tkdestroy(tt)
       classID<<-which(list.name==as.character(tclvalue(rbValue)))
    }
   Selected <- 0 
   OK.but <- tkbutton(tt,text="OK",command=function() Selected<-OnOK())
   tkgrid(OK.but)
   tkfocus(tt)
}

###################################
# set SNPs
###################################

selectSNPs<-function()
{
   list.name <<- colnames(RAW.data)
   cb <- NULL
   initialValues <- NULL
   cbValue <- NULL
   tt <- tktoplevel()
   for (i in 1:length(list.name)) {
        initTemp <- 0
        initialValues <- c(initialValues, initTemp)
   }

   tkgrid(tklabel(tt,text="Select SNPs:"))
   cbArray <- NULL
   for (i in 1:length(list.name) ) {
        cbArray[[i]] <- tkcheckbutton(tt)
        cbValue[[i]] <- tclVar(initialValues[i])
        tkconfigure(cbArray[[i]], variable=cbValue[[i]])
        tkgrid(tklabel(tt,text=list.name[i]),cbArray[[i]],sticky="e")
   }
   OnOK <- function()
   {
	for (i in 1:length(list.name)) {
		cbResult <<- c(cbResult, as.character(tclvalue(cbValue[[i]])))
	}
             SNPlist.id<<-which(cbResult=="1")
       tkdestroy(tt)
    }
   cbResult <- NULL
   OK.but <- tkbutton(tt,text="OK",command=OnOK)
   tkgrid(OK.but)
   tkfocus(tt)

}


###################################
# check Genotype
###################################

genotypeDialog <- function(SNP.id,entryInit="",entryWidth=20,returnValOnCancel="ID_CANCEL")
{   title<-names(RAW.data)[SNP.id]
    SNP.table<<-table(RAW.data[,SNP.id])
    modelD <- tktoplevel() ;   tkfocus(modelD);   tkwm.title(modelD,title)
    VarTcl<-list();Widget<-list()
    for(i in 1:length(SNP.table))
    { VarTcl[[i]]<-tclVar(paste(as.character(names(SNP.table)[i])))
      Widget[[i]]<-tkentry(modelD,width=paste(entryWidth),textvariable=VarTcl[[i]])
    }
   for(i in 1:length(SNP.table))
    {   tkgrid(tklabel(modelD,text="  "));   
        tkgrid(tklabel(modelD,text=as.character(names(SNP.table)[i])),Widget[[i]])
   }
    ReturnVal<<-names(SNP.table)
    onOK <- function()
    {  
         for(i in 1:length(SNP.table))
          if(tclvalue(VarTcl[[i]])!="") ReturnVal[i]<<- tclvalue(VarTcl[[i]])
         ReturnVal<<-matrix(ReturnVal)
         rownames(ReturnVal)<<-names(SNP.table)
          tkdestroy(modelD)
    }
    onCancel <- function()
    {     ReturnVal<<-returnValOnCancel
         tkdestroy(modelD)
    }
    OK.but     <-tkbutton(modelD,text="   OK   ",command=onOK)
    Cancel.but <-tkbutton(modelD,text=" Cancel ",command=onCancel)
    tkgrid(OK.but,Cancel.but) ;    tkgrid(tklabel(modelD,text="    "))
    tkfocus(modelD) ;    tkwait.window(modelD)
    return(ReturnVal)
} 

dataType<-function(SNP.id,SNP.match)
{ data<-as.character(RAW.data[,SNP.id])
  for(i in 1:nrow(SNP.match))
  data[data==rownames(SNP.match)[i]]<-SNP.match[i]
  return(data)
}

makeForm.data<-function(SNPlist.id)
{   Form.data<<-RAW.data
    for(i in SNPlist.id)
   {   a<-genotypeDialog(SNP.id=i)
        temp<-dataType(SNP.id=i,a) 
        Form.data[,i]<<-temp
   }
}


###################################
# set Case/Control
# output : Class.data, control.id, case.id
###################################

CaseControl.set<-function()
{  tt <- tktoplevel()
   cb.case<-list()
   cbValue.case<-list()
   cb.control<-list()
   cbValue.control<-list()
   class.table<-table(Form.data[,classID])
   tkwm.title(tt,"Select Case(Disease) and Control(Normal) group")
   tkgrid(tklabel(tt,text="Select Case(Disease) and Control(Normal) group"))
   tkgrid(tklabel(tt,text=""))
   tkgrid(tklabel(tt,text=""),tklabel(tt,text="Case(Disease)"),tklabel(tt,text="Control(Normal)"))

   for(i in 1:length(class.table))
   {   cb.case[[i]]<- tkcheckbutton(tt)
       cbValue.case[[i]] <- tclVar("0")
       tkconfigure(cb.case[[i]],variable=cbValue.case[[i]])
        cb.control[[i]]<- tkcheckbutton(tt)
       cbValue.control[[i]] <- tclVar("0")
       tkconfigure(cb.control[[i]],variable=cbValue.control[[i]])
        tkgrid(tklabel(tt,text=names(class.table)[i]),cb.case[[i]],cb.control[[i]]) 
   }
   OnOK <- function()
   {
       cbVal.case<<-NULL
       cbVal.control<<-NULL
       for(i in 1:length(table(Form.data[,classID])))
       {   cbVal.case <<- c(cbVal.case,as.character(tclvalue(cbValue.case[[i]])))
            cbVal.control <<- c(cbVal.control,as.character(tclvalue(cbValue.control[[i]])))
       }
      tkdestroy(tt)
   control.id<<-rep(FALSE,length(Form.data[,classID]))
   temp.id<-which(cbVal.control=="1")
   for(i in 1:length(temp.id))
     control.id[Form.data[,classID]==names(class.table)[temp.id[i]]]<<-TRUE
  case.id<<-rep(FALSE,length(Form.data[,classID]))
   temp.id<-which(cbVal.case=="1")
   for(i in 1:length(temp.id))
     case.id[Form.data[,classID]==names(class.table)[temp.id[i]]]<<-TRUE
  
   test.selected.id<<-which(control.id | case.id)
   Class.data<<-rep(NA,length(Form.data[,classID]))
   Class.data[control.id]<<-"Control"
   Class.data[case.id]<<-"Case"

   }
   OK.but <- tkbutton(tt,text="OK",command=OnOK)
   tkgrid(OK.but)
   tkfocus(tt)


}


###################################
# HWE
###################################

HWE<-function(SNPlist.id,Class.data)
{
     Control.id<-which(Class.data=="Control")
     HWE.table<-NULL
     for(i in SNPlist.id)
     {   temp<-as.matrix(Form.data[Control.id,i])
         temp[temp=="."]<-NA
         hap1<-apply(temp,1,function(x){strsplit(x,split="")[[1]][1]})
         hap2<-apply(temp,1,function(x){strsplit(x,split="")[[1]][2]})
         hap<-cbind(hap1,hap2)
         temp<-apply(hap,1,function(x){paste(x[1],x[2],sep="/")})
         temp[temp=="NA/NA"]<-NA
         temp.result<-HWE.chisq(genotype(temp))
         temp.table<-table(hap)
         sort.id<-sort.list(temp.table/sum(temp.table),decreasing=T)
         temp.table2<-table( factor(apply(hap,1,function(x){sum(x==names(temp.table)[sort.id[2]])})))
         temp.table2<-as.character(temp.table2[c("0","1","2")])
         temp.table2[is.na(temp.table2)]<-"0"
         temp.freq<-paste(round(temp.table/sum(temp.table),3)[sort.id],"(",temp.table[sort.id] ,")",sep="")
         temp.label<-c(names(temp.table)[sort.id[1]],temp.freq[1],names(temp.table)[sort.id[2]],temp.freq[2])
         HWE.table<-rbind(HWE.table,c(length(which(!is.na(temp))),temp.label, temp.table2,
                                 round(temp.result$statistic,3),round(temp.result$p.value,3),
                                 round(HWE.exact(genotype(temp))$p.value,3)))
    }  
    
    HWE.table<-cbind(colnames(Form.data)[SNPlist.id],HWE.table)

    HWE.table<-rbind(c("name","n","Major(C)","C.freq","Minor(R)","R.freq","CC","CR","RR","X2-Statistic","p-value","exact test p-value"),HWE.table)

    colnames(HWE.table)<-NULL
    rownames(HWE.table)<-NULL    
    HWE.result<<-HWE.table
   tclarray1 <- tclArray()
#    result.display(HWE.table,"Hardy-Weinberg Equilibrium test result (Control group)",10,tclarray)
    t.title<-"Hardy-Weinberg Equilibrium test result : Control ( "
#    t.title<-paste(t.title,")  Control(")
    t.temp<-table(Form.data[Control.id,classID]);t.temp<-t.temp[t.temp!=0]
    for(i in 1:length(t.temp))
    {  t.title<-paste(t.title,as.character(names(t.temp)[i]))
        if(i != length(t.temp)) t.title<-paste(t.title,",")
    }
    t.title<-paste(t.title,")")

  add.tab.list("HWE",t.title,HWE.result,10,tclarray1)
    
}

add.tab.list<-function(result.label1,result.label2,result.table,cw,tclarray)
{
    tbn <- tclvalue(tkadd(tn, label=result.label1))
    tkpack(tbw <- .Tk.newwin(tbn))
    tkpack(fr <- tkframe(tbw))
    tkpack(lb <- tklabel(fr, text=result.label2))
    ID <- paste(tn$ID, evalq(num.subwin<-num.subwin+1, tn$env), sep=".")
    win <- .Tk.newwin(ID)
    assign(ID, tbw, envir = tn$env)
    assign("parent", tn, envir = tbw$env)
   tclRequire("Tktable")
    for (i in (0:(nrow(result.table)-1)))
      for (j in (0:(ncol(result.table)-1)))
         tclarray[[i,j]] <- result.table[i+1,j+1]
    table1 <- tkwidget(tbw,"table",variable=tclarray,rows=nrow(result.table),cols=ncol(result.table),
               titlerows=1,selectmode="extended",colwidth=cw,background="white")

    tkpack(table1)
    #list(tbw,fr,lb) # return all three in case you need them later
}

result.display<-function(result.table,title,cw,tclarray)
{

   tclRequire("Tktable")
#   tclarray <- tclArray()
    for (i in (0:(nrow(result.table)-1)))
      for (j in (0:(ncol(result.table)-1)))
         tclarray[[i,j]] <- result.table[i+1,j+1]

    tt<-tktoplevel()
    tkwm.title(tt,title)
    table1 <- tkwidget(tt,"table",variable=tclarray,rows=nrow(result.table),cols=ncol(result.table),
               titlerows=1,selectmode="extended",colwidth=cw,background="white")
    tkpack(table1)
}


###################################
# allele X2
###################################

chiX2.allele.test<-function(SNPlist.id,Class.data)
{
     Control.id<-which(Class.data=="Control")
     Case.id<-which(Class.data=="Case")
     Class<-c(rep(0,length(Control.id)),rep(1,length(Case.id)))
     chiX2.table<-NULL
     for(i in SNPlist.id)
     {   temp<-as.matrix(Form.data[c(Control.id,Case.id),i])
         temp[temp=="."]<-NA
         hap1<-apply(temp,1,function(x){strsplit(x,split="")[[1]][1]})
         hap2<-apply(temp,1,function(x){strsplit(x,split="")[[1]][2]})
         hap<-c(hap1,hap2)
         temp.table<-table(rep(Class,2),hap)
         X2.result<-chisq.test(temp.table,correct=F)
         temp.freq<-temp.table/sum(temp.table)
         temp.table.hap<-table(hap)
         sort.id<-sort.list(temp.table.hap,decreasing=T)
         major.allele<-names(temp.table.hap)[sort.id[1]]
         minor.allele<-names(temp.table.hap)[sort.id[2]]
         temp.label<-NULL
          temp.label<-c(major.allele,
                               paste(round(temp.table[c(2,1),major.allele]/sum(temp.table),3),
                                          "(",temp.table[c(2,1),major.allele],")",sep=""),
                               minor.allele,
                               paste(round(temp.table[c(2,1),minor.allele]/sum(temp.table),3),
                                          "(",temp.table[c(2,1),minor.allele],")",sep=""))
                n<-length(Class[!is.na(temp)])
         chiX2.table.temp<-c(table(Class[!is.na(temp)])[c("1","0")],temp.label,round(X2.result$statistic,3),round(X2.result$p.value,4),
                                     round( fisher.test(temp.table)$p.value,4))
          chiX2.table<-rbind(chiX2.table,chiX2.table.temp)
    }  

    chiX2.table<-cbind(colnames(Form.data)[SNPlist.id],chiX2.table)

    chiX2.table<-rbind(c("name","n.case","n.cont.","Major(C)","Case","Cont.","Minor(R)","Case","Cont.",
                                  "X2-Statistic","p-value","exact test p-value"),chiX2.table)

    colnames(chiX2.table)<-NULL
    rownames(chiX2.table)<-NULL    
    chiX2.allele.result<<-chiX2.table
    t.title<-"Allele : X2 test result : Case("
    t.temp<-table(Form.data[Case.id,classID]);t.temp<-t.temp[t.temp!=0]
    for(i in 1:length(t.temp))
    { t.title<-paste(t.title,as.character(names(t.temp)[i]))
      if(i != length(t.temp)) t.title<-paste(t.title,",")
    }
    t.title<-paste(t.title,")  Control(")
    t.temp<-table(Form.data[Control.id,classID]);t.temp<-t.temp[t.temp!=0]
    for(i in 1:length(t.temp))
    {  t.title<-paste(t.title,as.character(names(t.temp)[i]))
        if(i != length(t.temp)) t.title<-paste(t.title,",")
    }
    t.title<-paste(t.title,")")
    chiX2.allele.result<<-chiX2.table
   tclarray2 <- tclArray()
#    result.display(chiX2.table,t.title,20,tclarray)
add.tab.list("allele:X2",t.title,chiX2.allele.result,10,tclarray2)

}

###################################
# genotype X2
###################################

chiX2.genotype.test<-function(SNPlist.id,Class.data)
{
     Control.id<-which(Class.data=="Control")
     Case.id<-which(Class.data=="Case")
     Class<-c(rep(0,length(Control.id)),rep(1,length(Case.id)))
     chiX2.table<-NULL
     for(i in SNPlist.id)
     {   temp<-as.matrix(Form.data[c(Control.id,Case.id),i])
         temp[temp=="."]<-NA
         hap1<-apply(temp,1,function(x){strsplit(x,split="")[[1]][1]})
         hap2<-apply(temp,1,function(x){strsplit(x,split="")[[1]][2]})
         hap<-cbind(hap1,hap2)
         temp.table1<-table(hap)    
         sort.id<-sort.list(temp.table1,decreasing=T)
         major.allele<-names(temp.table1)[sort.id[1]]
         minor.allele<-names(temp.table1)[sort.id[2]]
         temp<-factor(apply(hap,1,function(x){sum(x==names(temp.table1)[sort.id[2]])}),level=c(0,1,2))
         temp.table<-table(Class,temp);temp.table.test<-temp.table
         if( sum(apply(temp.table.test,2,sum)==0)!=0)
         {     id<-which(apply(temp.table.test,2,sum)==0)
               temp.table.test<-temp.table.test[,-id]
          }
         X2.result<-chisq.test(temp.table.test,correct=F)
         n<-length(Class[!is.na(temp)])
         temp.label<-c(paste(major.allele,major.allele,sep=""),
                               paste(round(temp.table[c(2,1),1]/n,3),"(",temp.table[c(2,1),1],")",sep=""),
                               paste(major.allele,minor.allele,sep=""),
                               paste(round(temp.table[c(2,1),2]/n,3),"(",temp.table[c(2,1),2],")",sep=""),
                               paste(minor.allele,minor.allele,sep=""),
                               paste(round(temp.table[c(2,1),3]/n,3),"(",temp.table[c(2,1),3],")",sep=""))
          chiX2.table.temp<-c(table(Class[!is.na(temp)])[c("1","0")],temp.label,round(X2.result$statistic,3),round(X2.result$p.value,4),
                                     round( fisher.test(temp.table)$p.value,4))
          chiX2.table<-rbind(chiX2.table,chiX2.table.temp)
   }  

    chiX2.table<-cbind(colnames(Form.data)[SNPlist.id],chiX2.table)

    chiX2.table<-rbind(c("name","n.case","n.cont.","CC","Case","Cont.","CR","Case","Cont.","RR","Case","Cont.",
                                       "X2-Statistic","p-value","exact test p-value"),chiX2.table)
    colnames(chiX2.table)<-NULL
    rownames(chiX2.table)<-NULL    
    chiX2.genotype.result<<-chiX2.table
    t.title<-"Genotype : X2 test result : Case("
    t.temp<-table(Form.data[Case.id,classID]);t.temp<-t.temp[t.temp!=0]
    for(i in 1:length(t.temp))
    { t.title<-paste(t.title,as.character(names(t.temp)[i]))
      if(i != length(t.temp)) t.title<-paste(t.title,",")
    }
    t.title<-paste(t.title,")  Control(")
    t.temp<-table(Form.data[Control.id,classID]);t.temp<-t.temp[t.temp!=0]
    for(i in 1:length(t.temp))
    {  t.title<-paste(t.title,as.character(names(t.temp)[i]))
        if(i != length(t.temp)) t.title<-paste(t.title,",")
    }
    t.title<-paste(t.title,")")
    chiX2.genotype.result<<-chiX2.table 
   tclarray3 <- tclArray()
#    result.display(chiX2.table,t.title,18,tclarray)
add.tab.list("genotype:X2",t.title,chiX2.genotype.result,10,tclarray3)

}

###################################
# Cochran-Armitage test
###################################
#
#CA.test<-function(SNPlist.id,Class.data)
#{
#     Control.id<-which(Class.data=="Control")
#     Case.id<-which(Class.data=="Case")
#     Class<-c(rep(0,length(Control.id)),rep(1,length(Case.id)))
#     CA.table<-NULL
#     for(i in SNPlist.id)
#     {   temp<-as.matrix(Form.data[c(Control.id,Case.id),i])
#         temp[temp=="."]<-NA
#        hap1<-apply(temp,1,function(x){strsplit(x,split="")[[1]][1]})
#         hap2<-apply(temp,1,function(x){strsplit(x,split="")[[1]][2]})
#         hap<-cbind(hap1,hap2)
#         temp.table<-sort(table(c(hap)),decreasing=T)
#         major.allele<-names(temp.table)[1]
#         temp.n<-apply(hap,1,function(x){length(which(x!=major.allele))})
#         temp.result<-Cochran.Armitage.test(temp.n, Class)
#        temp.table1<-table(hap)    
#         sort.id<-sort.list(temp.table1,decreasing=T)
#         major.allele<-names(temp.table1)[sort.id[1]]
#         minor.allele<-names(temp.table1)[sort.id[2]]
#         temp<-factor(apply(hap,1,function(x){sum(x==names(temp.table1)[sort.id[2]])}),level=c(0,1,2))
#         temp.table<-table(Class,temp);temp.table.test<-temp.table
#         if( sum(apply(temp.table.test,2,sum)==0)!=0)
#         {     id<-which(apply(temp.table.test,2,sum)==0)
#               temp.table.test<-temp.table.test[,-id]
#          }
#         n<-length(Class[!is.na(temp)])
#         temp.label<-c(paste(major.allele,major.allele,sep=""),
#                               paste(round(temp.table[c(2,1),1]/n,3),"(",temp.table[c(2,1),1],")",sep=""),
#                               paste(major.allele,minor.allele,sep=""),
#                               paste(round(temp.table[c(2,1),2]/n,3),"(",temp.table[c(2,1),2],")",sep=""),
#                               paste(minor.allele,minor.allele,sep=""),
#                               paste(round(temp.table[c(2,1),3]/n,3),"(",temp.table[c(2,1),3],")",sep=""))
#          CA.table<-rbind(CA.table,c(table(Class[!is.na(temp)])[c("1","0")],temp.label,round#(temp.result$statistic,3),round(temp.result$p.value,4)))
#    }  
#
#    CA.table<-cbind(colnames(Form.data)[SNPlist.id],CA.table)
#
#   CA.table<-rbind(c("name","n.case","n.cont.","CC","Case","Cont.","CR","Case","Cont.","RR","Case","Cont.",
#                                       "X2-Statistic","p-value"),CA.table)
#    colnames(CA.table)<-NULL
#    rownames(CA.table)<-NULL    
#    CA.result<<-CA.table
#    t.title<-"Cochran-Armitage test result : Case("
#    t.temp<-table(Form.data[Case.id,classID]);t.temp<-t.temp[t.temp!=0]
#    for(i in 1:length(t.temp))
#    { t.title<-paste(t.title,as.character(names(t.temp)[i]))
#      if(i != length(t.temp)) t.title<-paste(t.title,",")
#    }
#    t.title<-paste(t.title,")  Control(")
#    t.temp<-table(Form.data[Control.id,classID]);t.temp<-t.temp[t.temp!=0]
#    for(i in 1:length(t.temp))
#    {  t.title<-paste(t.title,as.character(names(t.temp)[i]))
#        if(i != length(t.temp)) t.title<-paste(t.title,",")
#    }
#    t.title<-paste(t.title,")")
#   tclarray4 <- tclArray()
##    result.display(chiX2.table,t.title,18,tclarray)
#   add.tab.list("Cochran-Armitage",t.title,CA.table,10,tclarray4)
#}
#

###################################
# regression
###################################
           

regression.test<-function(SNPlist.id,Class.data)
{    Control.id<-which(Class.data=="Control")
      Case.id<-which(Class.data=="Case")
      Class<-c(rep(0,length(Control.id)),rep(1,length(Case.id)))
      Cov.id<<-NULL
###select Y.test()

     list.name <<- colnames(RAW.data)[-SNPlist.id]
     cb <- NULL
     initialValues <- NULL
     cbValue <- NULL
     tt1 <- tktoplevel()
     for (i in 1:length(list.name)) {
         initTemp <- 0
         initialValues <- c(initialValues, initTemp)
     }
     tkgrid(tklabel(tt1,text="Select continuous response variable (Y):"))
     cbArray <- NULL
     rbValue<-tclVar("0")
     for (i in 1:length(list.name) ) {
        cbArray <- tkradiobutton(tt1)
        tkconfigure(cbArray, variable=rbValue,value=list.name[i])
        tkgrid(tklabel(tt1,text=list.name[i]),cbArray)
    }
    OnOK1 <- function()
    {
             Y.id<<-which(colnames(RAW.data)== tclvalue(rbValue))
             tkdestroy(tt1)
####
     list.name <<- colnames(RAW.data)[SNPlist.id]
     cb <- NULL
     initialValues <- NULL
     cbValue <- NULL
     tt2 <- tktoplevel()
     for (i in 1:length(list.name)) {
         initTemp <- 0
         initialValues <- c(initialValues, initTemp)
     }
     tkgrid(tklabel(tt2,text="Select SNPs:"))
     cbArray <- NULL
     for (i in 1:length(list.name) ) {
        cbArray[[i]] <- tkcheckbutton(tt2)
        cbValue[[i]] <- tclVar(initialValues[i])
        tkconfigure(cbArray[[i]], variable=cbValue[[i]])
        tkgrid(tklabel(tt2,text=list.name[i]),cbArray[[i]],sticky="e")
    }
    OnOK2 <- function()
    {
	for (i in 1:length(list.name)) {
		cbResult <<- c(cbResult, as.character(tclvalue(cbValue[[i]])))
	}
             SNPlist.id.test<<-SNPlist.id[which(cbResult=="1")]
             tkdestroy(tt2)

######### select SNPMode()
             temp.list.name <- colnames(RAW.data)[SNPlist.id.test]
             list.name<-NULL
             for(i in 1:length(temp.list.name))
             {  list.name<-c(list.name,paste(temp.list.name[i],"Additive",sep=" : "))
                temp<-substring("                                            ",1,nchar(temp.list.name[i]))
                list.name<-c(list.name,paste(temp,c("Dominant", "Recessive"),sep=" : "))
             }
             rb <- NULL
             tt3 <- tktoplevel()
             rbValue<-NULL
             tkgrid(tklabel(tt3,text="Select Class variable!"))
             for(ii in 1:length(temp.list.name))
             {   rbValue[[ii]] <- tclVar("0")
                 for (i in 1:3 ) {
                     rbtemp <- tkradiobutton(tt3)
                     tkconfigure(rbtemp, variable=rbValue[[ii]], value=list.name[3*(ii-1)+i])
                     tkgrid(tklabel(tt3,text=list.name[i]),rbtemp)
                }
            }
            OnOK3 <- function()
            {
               tkdestroy(tt3)
               SNPMode<<-NULL
               for(ii in 1:length(temp.list.name))
                   SNPMode<<-c(SNPMode,which(list.name[(3*(ii-1)+1):(3*ii)]==as.character(tclvalue(rbValue[[ii]])))) 
#####

              list.name <<- colnames(RAW.data)[-c(SNPlist.id,classID,Y.id)]
              cb <- NULL
              initialValues <- NULL
              cbValue <- NULL
              tt4 <- tktoplevel()
              for (i in 1:length(list.name)) {
                   initTemp <- 0
                   initialValues <- c(initialValues, initTemp)
              }
              tkgrid(tklabel(tt4,text="Select Covariates:"))
              cbArray <- NULL
              for (i in 1:length(list.name) ) {
                   cbArray[[i]] <- tkcheckbutton(tt4)
                   cbValue[[i]] <- tclVar(initialValues[i])
                   tkconfigure(cbArray[[i]], variable=cbValue[[i]])
                   tkgrid(tklabel(tt4,text=list.name[i]),cbArray[[i]],sticky="e")
              }
              OnOK4 <- function()
              {    for (i in 1:length(list.name)) {
	            cbResult <<- c(cbResult, as.character(tclvalue(cbValue[[i]])))
	      }
                   Cov.id<<-c(1:ncol(RAW.data))[-c(SNPlist.id,classID,Y.id)][which(cbResult=="1")] 
                   tkdestroy(tt4)
 
                   test.data<-as.matrix(log(Form.data[c(Control.id,Case.id),Y.id]))
                   temp.SNP<-as.matrix(Form.data[c(Control.id,Case.id),SNPlist.id.test])
                   for(i in 1:ncol(temp.SNP))
                  {    temp<-as.matrix(temp.SNP[,i])
                       temp[temp=="."]<-NA
                       hap1<-apply(temp,1,function(x){strsplit(x,split="")[[1]][1]})
                       hap2<-apply(temp,1,function(x){strsplit(x,split="")[[1]][2]})
                       hap<-cbind(hap1,hap2)
                       MA<-names(sort(table(hap),decreasing=T))[1]
                       if(SNPMode[i]==1)
                      {     temp.MODE<-apply(hap,1,function(x){2-sum(x==MA)})
                      } else if(SNPMode[i]==2) 
                      {     temp.MODE<-apply(hap,1,function(x){ifelse(sum(x==MA)==2,0,1)})
                      } else if(SNPMode[i]==3) 
                      {     temp.MODE<-apply(hap,1,function(x){ifelse(sum(x==MA)==0,1,2)})
                      } 
                      test.data<-cbind(test.data,as.matrix(temp.MODE))
                  }
                  colnames(test.data)<-c("Class",colnames(RAW.data)[SNPlist.id.test])
                  test.data<-cbind(test.data,Form.data[c(Control.id,Case.id),Cov.id])
                  colnames(test.data)<-c(colnames(RAW.data)[c(Y.id,SNPlist.id.test,Cov.id)])
                  test.data<-data.frame(test.data)
                  model.id<-c(SNPlist.id.test,Cov.id)
                  model<-paste(colnames(RAW.data)[Y.id]," ~ ",sep=" ")
                  for(i in 1:length(model.id))
                 {     model<-paste(model,colnames(RAW.data)[model.id[i]],sep="")
                       if(i!=length(model.id)) 
                             model<-paste(model,"+",sep=" ")
                 }
                 fit <- lm(as.formula(model), data=test.data)
                 REG.summary<<- matrix(round(summary(fit)$coefficients[-1,c(1,2,4)],3),ncol=3)
                 REG.summary<<-cbind(colnames(test.data)[-1],REG.summary)
                 REG.summary<<-rbind(c("name ","Estimate","Std.Error","P-value"),REG.summary)
                 model.t<-paste("log(",colnames(RAW.data)[Y.id],")"," ~ ",sep=" ")
                  for(i in 1:length(SNPlist.id.test))
                 {     model.t<-paste(model.t,colnames(RAW.data)[SNPlist.id.test[i]],
                                             "(",c("Add","Dom","Rec")[SNPMode[i]],")",sep="")
                       if(i!=length(SNPlist.id.test)) 
                             model.t<-paste(model.t,"+",sep=" ")
                 }
                 if(length(Cov.id)!=0)
                 {   for(i in 1:length(Cov.id))
                    {     model.t<-paste(model.t,"+",sep=" ")
                          model.t<-paste(model.t,colnames(RAW.data)[Cov.id[i]],sep=" ")
                    }
                 }
                 t.title<-"Regression result : Case("
                 t.temp<-table(Form.data[Case.id,classID]);t.temp<-t.temp[t.temp!=0]
                 for(i in 1:length(t.temp))
                 {   t.title<-paste(t.title,as.character(names(t.temp)[i]))
                    if(i != length(t.temp)) t.title<-paste(t.title,",")
                 }
                 t.title<-paste(t.title,")  Control(")
                 t.temp<-table(Form.data[Control.id,classID]);t.temp<-t.temp[t.temp!=0]
                 for(i in 1:length(t.temp))
                {     t.title<-paste(t.title,as.character(names(t.temp)[i]))
                      if(i != length(t.temp)) t.title<-paste(t.title,",")
                }
                t.title<-paste(t.title,") : Model  ",model.t,sep="")

                 tclarray.Reg<-tclArray()
                 add.tab.list("Reg.",t.title,REG.summary,15,tclarray.Reg)
             }
             cbResult <- NULL
             OK.but4<- tkbutton(tt4,text="OK",command=OnOK4)
             tkgrid(OK.but4)
             tkfocus(tt4)
            }
              OK.but3 <- tkbutton(tt3,text="OK",command=function() Selected<-OnOK3())
              tkgrid(OK.but3)
              tkfocus(tt3)
    }
   cbResult <- NULL
   OK.but2 <- tkbutton(tt2,text="OK",command=OnOK2)
   tkgrid(OK.but2)
   tkfocus(tt2)
    }
  cbResult <- NULL
   OK.but1 <- tkbutton(tt1,text="OK",command=OnOK1)
   tkgrid(OK.but1)
   tkfocus(tt1)

 } # end of regression function


###################################
# logistic regression
###################################
 

logistic.test<-function(SNPlist.id,Class.data)
{    Control.id<-which(Class.data=="Control")
      Case.id<-which(Class.data=="Case")
      Class<-c(rep(0,length(Control.id)),rep(1,length(Case.id)))
      Cov.id<<-NULL
###selectSNPs.test()

     list.name <<- colnames(RAW.data)[SNPlist.id]
     cb <- NULL
     initialValues <- NULL
     cbValue <- NULL
     tt1 <- tktoplevel()
     for (i in 1:length(list.name)) {
         initTemp <- 0
         initialValues <- c(initialValues, initTemp)
     }
     tkgrid(tklabel(tt1,text="Select SNPs:"))
     cbArray <- NULL
     for (i in 1:length(list.name) ) {
        cbArray[[i]] <- tkcheckbutton(tt1)
        cbValue[[i]] <- tclVar(initialValues[i])
        tkconfigure(cbArray[[i]], variable=cbValue[[i]])
        tkgrid(tklabel(tt1,text=list.name[i]),cbArray[[i]],sticky="e")
    }
    OnOK1 <- function()
    {
	for (i in 1:length(list.name)) {
		cbResult <<- c(cbResult, as.character(tclvalue(cbValue[[i]])))
	}
            SNPlist.id.test<<-SNPlist.id[which(cbResult=="1")]
             tkdestroy(tt1)

#######  SNPMode(SNPlist.id.test)
             temp.list.name <- colnames(RAW.data)[SNPlist.id.test]
             list.name<<-NULL
             for(i in 1:length(temp.list.name))
             {  list.name<<-c(list.name,paste(temp.list.name[i],"Additive",sep=" : "))
                temp<-substring("                                            ",1,nchar(temp.list.name[i]))
                list.name<<-c(list.name,paste(temp,c("Dominant", "Recessive"),sep=" : "))
             }
             rb <- NULL
             tt2 <- tktoplevel()
             rbValue<-NULL
             tkgrid(tklabel(tt2,text="Select Class variable!"))
             for(ii in 1:length(temp.list.name))
             {   rbValue[[ii]] <- tclVar("0")
                 for (i in 1:3 ) {
                     rbtemp <- tkradiobutton(tt2)
                     tkconfigure(rbtemp, variable=rbValue[[ii]], value=list.name[3*(ii-1)+i])
                     tkgrid(tklabel(tt2,text=list.name[i]),rbtemp)
                }
            }
            OnOK2 <- function()
            {
               tkdestroy(tt2)
               SNPMode<<-NULL
               temp.list.name <- colnames(RAW.data)[SNPlist.id.test]
               for(ii in 1:length(temp.list.name))
                   SNPMode<<-c(SNPMode,which(list.name[(3*(ii-1)+1):(3*ii)]==as.character(tclvalue(rbValue[[ii]]))))

######### selectCovariate()

              list.name <<- colnames(RAW.data)[-c(SNPlist.id,classID)]
              cb <- NULL
              initialValues <- NULL
              cbValue <- NULL
              tt3 <- tktoplevel()
              for (i in 1:length(list.name)) {
                   initTemp <- 0
                   initialValues <- c(initialValues, initTemp)
              }
              tkgrid(tklabel(tt3,text="Select Covariates:"))
              cbArray <- NULL
              for (i in 1:length(list.name) ) {
                   cbArray[[i]] <- tkcheckbutton(tt3)
                   cbValue[[i]] <- tclVar(initialValues[i])
                   tkconfigure(cbArray[[i]], variable=cbValue[[i]])
                   tkgrid(tklabel(tt3,text=list.name[i]),cbArray[[i]],sticky="e")
              }
              OnOK3 <- function()
              {    for (i in 1:length(list.name)) {
	            cbResult <<- c(cbResult, as.character(tclvalue(cbValue[[i]])))
	      }
                  Cov.id<<-c(1:ncol(RAW.data))[-c(SNPlist.id,classID)][which(cbResult=="1")] 
  
                   tkdestroy(tt3)
                   test.data<-as.matrix(Class)
                   temp.SNP<-as.matrix(Form.data[c(Control.id,Case.id),SNPlist.id.test])
                   for(i in 1:ncol(temp.SNP))
                  {    temp<-as.matrix(temp.SNP[,i])
                       temp[temp=="."]<-NA
                       hap1<-apply(temp,1,function(x){strsplit(x,split="")[[1]][1]})
                       hap2<-apply(temp,1,function(x){strsplit(x,split="")[[1]][2]})
                       hap<-cbind(hap1,hap2)
                       MA<-names(sort(table(hap),decreasing=T))[1]
                       if(SNPMode[i]==1)
                      {     temp.MODE<-apply(hap,1,function(x){2-sum(x==MA)})
                      } else if(SNPMode[i]==2) 
                      {     temp.MODE<-apply(hap,1,function(x){ifelse(sum(x==MA)==2,0,1)})
                      } else if(SNPMode[i]==3) 
                      {     temp.MODE<-apply(hap,1,function(x){ifelse(sum(x==MA)==0,1,2)})
                      } 
                      test.data<-cbind(test.data,as.matrix(temp.MODE))
                  }
                  colnames(test.data)<-c("Class",colnames(RAW.data)[SNPlist.id.test])
                  test.data<-cbind(test.data,Form.data[c(Control.id,Case.id),Cov.id])
                  colnames(test.data)<-c("Class",colnames(RAW.data)[c(SNPlist.id.test,Cov.id)])
                  test.data<-data.frame(test.data)
                  model.id<-c(SNPlist.id.test,Cov.id)
                  model<-"Class ~"
                  for(i in 1:length(model.id))
                 {     model<-paste(model,colnames(RAW.data)[model.id[i]],sep=" ")
                       if(i!=length(model.id)) 
                             model<-paste(model,"+",sep=" ")
                 }
                 fit <- glm(as.formula(model), family=binomial,data=test.data)
                 summary(fit)$coefficients
                 model.summary<- matrix(round(summary(fit)$coefficients[-1,c(1,2,4)],3),ncol=3)
                 OR.temp<-cbind(exp(model.summary[,1]),exp(model.summary[,1]-model.summary[,2]*1.96),
                                                                           exp(model.summary[,1]+model.summary[,2]*1.96))
                 model.summary<-round(cbind(model.summary,as.matrix(OR.temp)),3)
                 model.summary<-cbind(colnames(test.data)[-1],model.summary)
                 model.summary<-cbind(rownames(model.summary),model.summary)
                 colnames(model.summary)<-c("name ","Estimate","Std.Error","P-value","Odds Ratio","OR:Lower","OR:Upper")
                 model.summary<-rbind(colnames(model.summary),model.summary)
               model.summary[seq(from=nrow(model.summary),by=-1,length=length(Cov.id)),5:7]<-" "
                Log.result<<-model.summary
                t.title<-"Logistic regression result : Case("
                t.temp<-table(Form.data[Case.id,classID]);t.temp<-t.temp[t.temp!=0]
                for(i in 1:length(t.temp))
               {    t.title<-paste(t.title,as.character(names(t.temp)[i]))
                    if(i != length(t.temp)) t.title<-paste(t.title,",")
               }
               t.title<-paste(t.title,")  Control(")
               t.temp<-table(Form.data[Control.id,classID]);t.temp<-t.temp[t.temp!=0]
               for(i in 1:length(t.temp))
               {    t.title<-paste(t.title,as.character(names(t.temp)[i]))
                   if(i != length(t.temp)) t.title<-paste(t.title,",")
               }
                t.title<-paste(t.title,") : Model  ",model,sep="")
                tclarray5<-tclArray() 
                add.tab.list("Logistic Reg.",paste("Logistic regression result ",model,sep=" : "),Log.result,10,tclarray5)
             }
             cbResult <- NULL
             OK.but3<- tkbutton(tt3,text="OK",command=OnOK3)
             tkgrid(OK.but3)
             tkfocus(tt3)
        }
        Selected <- 0 
       OK.but2 <- tkbutton(tt2,text="OK",command=function() Selected<-OnOK2())
       tkgrid(OK.but2)
       tkfocus(tt2)
   }
   cbResult <- NULL
   OK.but1 <- tkbutton(tt1,text="OK",command=OnOK1)
   tkgrid(OK.but1)
   tkfocus(tt1)
}


###################################
# LD plot
###################################
 
LD.plot<-function(SNPlist.id,Class.data)
{     Control.id<-which(Class.data=="Control")
      Case.id<-which(Class.data=="Case")
      Class<-c(rep(0,length(Control.id)),rep(1,length(Case.id)))

   list.name <<- colnames(RAW.data)[SNPlist.id]
   cb <- NULL
   initialValues <- NULL
   cbValue <- NULL
   tt <- tktoplevel()
   for (i in 1:length(list.name)) {
        initTemp <- 0
        initialValues <- c(initialValues, initTemp)
   }

   tkgrid(tklabel(tt,text="Select SNPs:"))
   cbArray <- NULL
   for (i in 1:length(list.name) ) {
        cbArray[[i]] <- tkcheckbutton(tt)
        cbValue[[i]] <- tclVar(initialValues[i])
        tkconfigure(cbArray[[i]], variable=cbValue[[i]])
        tkgrid(tklabel(tt,text=list.name[i]),cbArray[[i]],sticky="e")
   }
   OnOK <- function()
   {    for (i in 1:length(list.name)) {
	cbResult <<- c(cbResult, as.character(tclvalue(cbValue[[i]])))
        }
        SNPlist.id.test<<-SNPlist.id[which(cbResult=="1")]
        tkdestroy(tt)
        hap.data<-NULL
        for(i in 1:length(SNPlist.id.test))
        {  temp<-as.matrix(Form.data[c(Control.id,Case.id),SNPlist.id.test[i]])
            temp[temp=="."]<-NA          
            hap1<-apply(temp,1,function(x){strsplit(x,split="")[[1]][1]})
            hap2<-apply(temp,1,function(x){strsplit(x,split="")[[1]][2]})
            hap<-cbind(hap1,hap2)
            temp<-apply(hap,1,function(x){paste(x[1],x[2],sep="/")})
            temp[temp=="NA/NA"]<-NA
            hap.data<-cbind(hap.data,temp)
        }
        colnames(hap.data)<-colnames(Form.data)[SNPlist.id.test]
        data<-makeGenotypes(data.frame(hap.data))
        LD.result<-LD(data)
        tclarray<-tclArray()        
        add.plot.list("LD plot","What is it",LD.result,tclarray)
    }
   cbResult <- NULL
   OK.but <- tkbutton(tt,text="OK",command=OnOK)
   tkgrid(OK.but)
   tkfocus(tt)
}

add.plot.list<-function(result.label1,result.label2,LD.result,tclarray)
{
        result.label1<-"LD plot"
        result.label2<-"What is it"
        cw<-15
        tclarray<-tclArray()

        tbn <- tclvalue(tkadd(tn, label=result.label1))
        tkpack(tbw <- .Tk.newwin(tbn))
        tkpack(fr <- tkframe(tbw))
        tkpack(lb <- tklabel(fr, text=result.label2))
        ID <- paste(tn$ID, evalq(num.subwin<-num.subwin+1, tn$env), sep=".")
        win <- .Tk.newwin(ID)
        assign(ID, tbw, envir = tn$env)
        assign("parent", tn, envir = tbw$env)
        tkpack(tkrplot(fr,function(){LDtable(LD.result)}))
}
 

###################################
#haploscore
###################################
 

haploScore<-function(SNPlist.id,Class.data)
{
      Control.id<-which(Class.data=="Control")
      Case.id<-which(Class.data=="Case")
      Class<-c(rep(0,length(Control.id)),rep(1,length(Case.id)))

      list.name <<- colnames(RAW.data)[SNPlist.id]
      cb <- NULL
      initialValues <- NULL
      cbValue <- NULL
      tt <- tktoplevel()
      for (i in 1:length(list.name)) {
           initTemp <- 0
           initialValues <- c(initialValues, initTemp)
      }

      tkgrid(tklabel(tt,text="Select SNPs:"))
      cbArray <- NULL
      for (i in 1:length(list.name) ) {
           cbArray[[i]] <- tkcheckbutton(tt)
           cbValue[[i]] <- tclVar(initialValues[i])
           tkconfigure(cbArray[[i]], variable=cbValue[[i]])
           tkgrid(tklabel(tt,text=list.name[i]),cbArray[[i]],sticky="e")
      }
      OnOK <- function()
      {     for (i in 1:length(list.name)) {
	    cbResult <<- c(cbResult, as.character(tclvalue(cbValue[[i]])))
            }
            SNPlist.id.test<<-SNPlist.id[which(cbResult=="1")]
            tkdestroy(tt)
            hap.data<-NULL
            for(i in 1:length(SNPlist.id.test))
            {    temp<-as.matrix(Form.data[c(Control.id,Case.id),SNPlist.id.test[i]])
                 temp[temp=="."]<-NA          
                 hap1<-apply(temp,1,function(x){strsplit(x,split="")[[1]][1]})
                 hap2<-apply(temp,1,function(x){strsplit(x,split="")[[1]][2]})
                 hap<-cbind(hap1,hap2)
                 temp<-apply(hap,1,function(x){paste(x[1],x[2],sep="/")})
                 temp[temp=="NA/NA"]<-NA
                 hap.data<-cbind(hap.data,temp)
            }
            colnames(hap.data)<-colnames(Form.data)[SNPlist.id.test]

            HAP1<-NULL
            for(i in 1:ncol(hap.data))
            {    HAP1<-cbind(HAP1,genetics::hap(hap.data[,i])$hap1,genetics::hap(hap.data[,i])$hap2)
            }
            HAP1[is.na(HAP1)]<-"Z"
            HAP1<-data.frame(id=1:nrow(HAP1),class=Class,HAP1)
            HAP1.result<-hap.em(id=HAP1$id,data=HAP1[,-c(1:2)],locus.label=colnames(RAW.data)[SNPlist.id.test])

            score.id<-factor(Class)
            temp.result<-hap.score(y=score.id,geno=HAP1[,-c(1:2)],locus.label=colnames(RAW.data)[SNPlist.id.test])

            haplotype<-apply(temp.result$haplotype,1,
                       function(x){   h.name<-NULL  ;   for(i in 1:length(x))  h.name<-paste(h.name,x[i],sep="");
                                          return(h.name) })


            sort.id<-sort.list(temp.result$hap.prob,decreasing=T)

            haplo.result<-cbind(haplotype[sort.id],round(cbind(temp.result$hap.prob[sort.id], 
            temp.result$score.haplo[sort.id], temp.result$score.haplo.p[sort.id]),4))
            haplo.result<-rbind(c("haplotype","frequency","score","p-value"),haplo.result)

            title2<-"haplotype \n"
            for(i in 1:length(SNPlist.id.test))
            {    title2<-paste(title2,colnames(RAW.data)[SNPlist.id.test[i]],sep=" ")
                 if(i%%5==0)   
                 {    title2<-paste(title2,"\n",sep="")
                 } else if(i!=length(SNPlist.id.test))
                 {    title2<-paste(title2,",",sep="")
                 }
            }
               
            tclarray<-tclArray()
            add.tab.list("hap.score",title2,haplo.result,15,tclarray)
            
      }
      cbResult <- NULL
      OK.but <- tkbutton(tt,text="OK",command=OnOK)
      tkgrid(OK.but)
      tkfocus(tt)
}
  

###################################
# MDR
###################################
 
MDR<-function(SNPlist.id,Class.data)
{  
      Control.id<-which(Class.data=="Control")
      Case.id<-which(Class.data=="Case")
      Class<-c(rep(0,length(Control.id)),rep(1,length(Case.id)))



      list.name <<- colnames(RAW.data)[SNPlist.id]
      cb <- NULL
      initialValues <- NULL
      cbValue <- NULL
      tt <- tktoplevel()
      for (i in 1:length(list.name)) {
           initTemp <- 0
           initialValues <- c(initialValues, initTemp)
      }

      tkgrid(tklabel(tt,text="Select SNPs:"))
      cbArray <- NULL
      for (i in 1:length(list.name) ) {
           cbArray[[i]] <- tkcheckbutton(tt)
           cbValue[[i]] <- tclVar(initialValues[i])
           tkconfigure(cbArray[[i]], variable=cbValue[[i]])
           tkgrid(tklabel(tt,text=list.name[i]),cbArray[[i]],sticky="e")
      }
      textEntryVarTcl <- tclVar(paste("  "))
      textEntryWidget <- tkentry(tt,width=paste(20),textvariable=textEntryVarTcl)
      tkgrid(tklabel(tt,text="       "))
      tkgrid(tklabel(tt,text="number of SNP combinations"),textEntryWidget)
      tkgrid(tklabel(tt,text="       "))
      ReturnVal <- "ID_CANCEL"
      OnOK <- function()
      {     for (i in 1:length(list.name)) {
	    cbResult <<- c(cbResult, as.character(tclvalue(cbValue[[i]])))
            }
            nbr<<- as.numeric(tclvalue(textEntryVarTcl))

            SNPlist.id.test<<-SNPlist.id[which(cbResult=="1")]
#  SNPlist.id.test<-SNPlist.id[1:6]
            tkdestroy(tt)
            mdr.data<-as.matrix(Class)
            temp.SNP<-as.matrix(Form.data[c(Control.id,Case.id),SNPlist.id.test])
            for(i in 1:ncol(temp.SNP))
            {    temp<-as.matrix(temp.SNP[,i])
                 temp[temp=="."]<-NA
                 hap1<-apply(temp,1,function(x){strsplit(x,split="")[[1]][1]})
                 hap2<-apply(temp,1,function(x){strsplit(x,split="")[[1]][2]})
                 hap<-cbind(hap1,hap2)
                 MA<-names(sort(table(hap),decreasing=T))[1]
                 temp.MODE<-apply(hap,1,function(x){2-sum(x==MA)})
                 mdr.data<-cbind(mdr.data,as.matrix(temp.MODE))
            }
            colnames(mdr.data)<-c("class",colnames(RAW.data)[SNPlist.id.test])
            for(i in 1:ncol(mdr.data))
                 if(length(which(is.na(mdr.data[,i])))!=0)  mdr.data<-mdr.data[-c(which(is.na(mdr.data[,i]))),]
            mdr.result<-mdr.c(mdr.data,colresp=1,cs=1,combi=nbr,cv.fold=10)
 # mdr.result1<-rmdr.fast(mdr.data,cvk=10,nbr)
            MDR.result<-NULL
            for(i in 1:nbr)
            {  if(i!=1)  
               {    MDR.result<-paste(MDR.result, colnames(RAW.data)[SNPlist.id.test][mdr.result$min.comb[i,]],sep=" : ")
               } else
               { MDR.result<-paste(MDR.result, colnames(RAW.data)[SNPlist.id.test][mdr.result$min.comb[i,]],sep="  ")
               }
            }
            MDR.result<-cbind(MDR.result,round(mdr.result$train.erate,3),round(mdr.result$test.erate,3))
            MDR.result<-rbind(c("SNPs","Miss. Class. Error","Predict Error"),MDR.result)      
            title2<-"Multi-Dimensional Reduction \n"
            for(i in 1:length(SNPlist.id.test))
            {    title2<-paste(title2,colnames(RAW.data)[SNPlist.id.test[i]],sep=" ")
                 if(i%%5==0)   
                 {    title2<-paste(title2,"\n",sep="")
                 } else if(i!=length(SNPlist.id.test))
                 {    title2<-paste(title2,",",sep="")
                 }
            }
            tclarray.MDR<-tclArray()
            add.tab.list("MDR",title2,MDR.result,15,tclarray.MDR)  
## ORMDR

            ormdr.result<-ormdr(mdr.data,bestcombi=as.numeric(mdr.result$best.combi),cs=1,colresp=1,CI.Asy=TRUE,CI.Boot=TRUE,B=5000)
            ormdr.result<-rbind(colnames(ormdr.result),ormdr.result)
            title3<-"Multi-Dimensional Reduction \n Best combination : "
            for(i in 1:length(mdr.result$best.combi))
            {    title3<-paste(title3,colnames(RAW.data)[SNPlist.id.test[as.numeric(mdr.result$best.combi[i])]],sep=" ")
                 if(i%%5==0)   
                 {    title3<-paste(title3,"\n",sep="")
                 } else if(i!=length(SNPlist.id.test))
                 {    title3<-paste(title3,",",sep="")
                 }
            }

              tclarray.ORMDR<-tclArray()
            add.tab.list("ORMDR",title3,ormdr.result,15,tclarray.ORMDR)  
     
      }
      cbResult <- NULL
      OK.but <- tkbutton(tt,text="OK",command=OnOK)
      tkgrid(OK.but)
      tkfocus(tt)
}

}
 