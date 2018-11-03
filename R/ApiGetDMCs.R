#' Get TOPK=500 DMCs and non-DMCs using moderated-t test
#' @param betaValue    A matrix from TCGA array data
#' @param TOPK         An integer number, default 500. Number of DMCs/non-DMCs.
#' @param tumorNum     A postive number, First tumorNum columns in betaValue
#'                     are tumor samples. If tumorNum is NULL, first half of
#'                     columns are considered as tumor samples,
#' @param filterProbes Logistic. defalut is FALSE. The code use all probes in betaValue.
#'                     If TRUE,
#'                     you can use default good probes provided in our code.
#'                     you can also provide your good probes in userProbes.
#' @param userProbes   A number list. The row numbers in betaValue.
#'                     These rows are considered as good probes.
#' return DMCs        (TOPK DMCs and TOPK non-DMCs row index in betaValue)
#' @note              User can provide the good probes indexes (row number)
#'                    to filter the probes.
#'                    A global variable goodProbes are used in this function.
#'                    goodProbes: probes with SNPs at the CpG or single base
#'                    extension sites, and corss-reative probes are removed.
#'                    More details see the reference paper.
#' @export
#'
ApiGetDMCs <-function(betaValue,TOPK=500,tumorNum=NULL,
                      filterProbes=FALSE,userProbes=NULL){
  ###  get the TOPK=500 MDCs using moderated-t test
  anrow<-dim(betaValue)[1] #row number
  ancol<-dim(betaValue)[2] #col number
  #goodProbes<-NULL
  #data(goodProbes,)


  if (filterProbes){
    ### filter probes
    if (is.null(userProbes)){
      #using good probes same as paper
      tmpidx<-goodProbes<=anrow
      #only keep probs < row number FALSE/TRUE
      xbetaSub<-betaValue[goodProbes[tmpidx],] #
      #xannoSub<-annot[goodProbes[tmpidx],]
      xannoSub<-annotGeneNames[goodProbes[tmpidx]]

    }else{
      #using probes provided by user
      tmpidx<-userProbes
      xbetaSub<-betaValue[tmpidx,]
      #xannoSub<-annot[tmpidx,]
      xannoSub<-annotGeneNames[tmpidx]
    }

  }else{
    ### not filter probes
    xbetaSub=betaValue
    #xannoSub=annot
    xannoSub=annotGeneNames
  }
  ### remove NA
  myidx<-which(!is.na(rowSums(xbetaSub)))
  betaSub<-xbetaSub[myidx,]  ## remove 'NA'
  #annoSub<-xannoSub[myidx,] ## remove id contain 'NA'
  annoSub<-xannoSub[myidx]
  ### moderated-t test of logit -beta value to find DMCs and non-DMCs
  Lbeta=.logit(betaSub)

  if (is.null(tumorNum)){
    tumorNum=ancol/2 ## T_1-T_n tumor/cancer;   N_1-N_n normal
  }
  #n=ncol(Lbeta)/2
  #
  mixedT=Lbeta[,1:tumorNum]  ## tumor / cancer
  pureN=Lbeta[,(tumorNum+1):ancol] ## normal
  controlNum=ancol-tumorNum
  #n<-tumorNum
  #n<-controlNum
  mv =cbind(pureN,mixedT) ## Lbeta #cbind(pureN, mixedT)
  ### mv=cbind(mixedT,pureN)
  #pd = c(paste("N",1:n, sep = ""), paste("C",1:n, sep = "")) #normal cancer
  #colnames(mv) = pd
  design <- stats::model.matrix(~0 + factor(c(rep("Control",controlNum), rep("Tumor",tumorNum)))) #control /tumor
  colnames(design) <- c("Control", "Tumor")  # control, tumor
  ### moderated t-test   functions in limma packages
  #library(limma)
  contrast.matrix <- limma::makeContrasts("Tumor-Control", levels=design)
  fit0 <- limma::lmFit(mv, design)
  fit1 <- limma::contrasts.fit(fit0, contrast.matrix)
  fit2 <- limma::eBayes(fit1)
  pv = fit2$p.value
  #
  ior = order(pv) ## increasing p-value
  ior2 = order(pv, decreasing = T) ## decreasing p-values
  # # # top 500
  #print(length(ior))
  if(length(ior)<TOPK ){
    myDMCs<-NULL
    stop("TOPK length is larger than dataset")
  }
  else{
    #DMCs=annoSub[ior[1:TOPK],]$`Composite Element REF`#Chromosome
    #nonDMCs=annoSub[ior2[1:TOPK],]$`Composite Element REF`
    #DMCs=ior[1:TOPK]
    #nonDMCs=ior2[1:TOPK]
    DMCs=annoSub[ior[1:TOPK]]
    nonDMCs=annoSub[ior2[1:TOPK]]
    myDMCs=cbind(DMCs,nonDMCs)
    #print(nonDMCs)
    #print(DMCs)
  }

   return(myDMCs)
}

##if(FALSE){


### ------ logit
.logit <-function (a) {
  #logit transform
  x<-a
  return(log2(x/(1-x)))
}

#}
