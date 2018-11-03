#' Bayesian Purity Model (BPM) Main functions.
#'
#' @param  betaValue   A matrix,TCGA methlation array data. Each row: loci,
#'                     Tumor1,Tumor2,...,Normal1,Nomral2,...
#'
#' @param  TOPK        A number. Number of DMCs/nonDMCs selected
#' @param  tumorNum    The number of tumor samples.
#'                     if NULL, the default number is half of column number of dataset.
#' @param filterProbes Logistic. defalut is FALSE. The code use all probes in betaValue.
#'                     If TRUE,
#'                     you can use default good probes provided in our code.
#'                     you can also provide your good probes in userProbes.
#' @param userProbes   A number list. The row numbers in betaValue.
#'                     These rows are considered as good probes.

#' @return tumor purity estimation of tumor samples
#' @export
#' @examples
#' ### need to install package "limma"
#' ### source("https://bioconductor.org/biocLite.R");biocLite("limma");
#' BayPM(simUCEC,20,2);
#'
BayPM <-function(betaValue,TOPK=500,tumorNum=NULL, filterProbes=FALSE,userProbes=NULL){
  ###dmcf: DMCs sites fiel  (differentially methylated CpG sites)
  ### filter the dataset
  #anrow<-dim(betaValue)[1] #row number
  ancol<-dim(betaValue)[2] #col number



  #############============================
  #cat("Find DMCs  (in ApiGetDMCs.R).. \n")
  ### get DMCs from  ApiGetDMCs.R
  DMCs_nonDMCs=ApiGetDMCs(betaValue,TOPK,tumorNum,filterProbes,userProbes)
  #write.csv(DMCs_nonDMCs, file=dmcfn) ## DMCs-non-DMCs
  cat("DMCs and non-DMCs were identified by moderated-t statistics.. \n")
  if(is.null(DMCs_nonDMCs)){
    stop("not find DMCs or non-DMCs.")
  }else{
    DMCs<-DMCs_nonDMCs[,1] # DMCs
    nonDMCs<-DMCs_nonDMCs[,2] # non-DMCs

  }

  ###find DMCs's index in annot
  #DMCsidx=match(DMCs,annot$`Composite Element REF`)
  #find nonDMCs's index in annot
  #nonDMCsidx=match(nonDMCs,annot$`Composite Element REF`)
  #annotGeneNames<-NULL
  #data(annotGeneNames,envir=enviroment())

  DMCsidx=match(DMCs,annotGeneNames)
  nonDMCsidx=match(nonDMCs,annotGeneNames)

  #print(dim(betaValue))
  #print("=======")
  #print(DMCsidx)
  #print("==========")
  DMCsBeta=betaValue[DMCsidx,] # beta value of DMCs
  nonDMCsBeta=betaValue[nonDMCsidx,] # beta value of non-DMCs

  ##ancol=ncol(betaValue) #column number of Lbeta

  if (is.null(tumorNum)){
    tumorNum=ancol/2 # cancer,normal.
    print("TumorNum not provided. First half columns used as tumor samples.")
  }
  ### rt=50 #repeat 50 times
  ### top500 DMCs, top500 non-DMCs
  m=TOPK #
  #n=tumorNum #sample number
  alp_est=rep(NA,tumorNum) #matrix(NA,n,rt) #store alpha:tumor purity
  xbar_est=rep(NA,m) #matrix(NA,m,rt)
  #store xi: mode of beta-value of each row of tumor samples x
  #z=nonDMCsBeta[,1:sampleN] #mix cancer   loci*sampleN
  #y=nonDMCsBeta[,(sampleN+1):(sampleN+sampleN)] #normal


  # ===============Estimate phi
  # Step1:  estimate phi  using  mean/var of y (normal)
  cat("Step1 estimate phi ...\n")
  y_all=betaValue[,(tumorNum+1):(ancol)] #normal samples
  # mean/var of y: m1/v1
  m1 = rowMeans(y_all, na.rm = T)
  v1 = apply(y_all, 1, stats::var)
  # estimate phi: mode of beta-value of each row of control samples y
  phi = m1 + (2*m1 - 1)*v1 / (m1*(1-m1) - 3*v1)
  # phi=y_bar if phi outside (0,1)
  # phi[phi >= 1 | phi <= 0] = m1[phi >= 1 | phi <= 0]
  phi[phi>=1 &!is.na(phi)]=m1[phi>=1 & !is.na(phi)]
  phi[phi<=0 &!is.na(phi)]=m1[phi<=0 & !is.na(phi)]


  # ================Estmate v_z noise intensity
  # Step2/3: estimate v_z using non-DMCs by maximum likelihood estimation
  cat("Step2/3 estimate v_z ...\n")
  nonDMCs_z<-nonDMCsBeta[,1:tumorNum]  # z is tumor
  #print(c(dim(phi),length(phi),length(nonDMCsidx)))
  nonDMCs_phi<-phi[nonDMCsidx]
  nu_z <- estimateNu(nonDMCs_z,nonDMCs_phi)


  # ===============Sampleing xi and alpha(tumor purity)
  # STep4: sampel xi alpha form DMCs
  cat("Step4 sample alpha and xi ...\n")
  DMCs_z=DMCsBeta[,1:tumorNum]  # tumor/cancer  DMCs
  DMCs_y=DMCsBeta[,(tumorNum+1):ancol] #normal in DMCs
  nm = rowMeans(DMCs_y, na.rm = T) #nomral mean
  tm = rowMeans(DMCs_z, na.rm = T) #tumor mean
  # hypermethylation 1; hypomethylation 2;
  mstates = rep(NA, m)
  mstates[which(tm > nm)] = 1 #hyper
  mstates[which(tm < nm)] = 2 #hypo

  set.seed(0)
  # fullSampler from purityGS_mode.R;
  # output: x_bar,x_last,x_sample,xpar,nab(nv),alp(alpha)
  res = fullSampler(DMCs_y, DMCs_z, mstates, matrix(c(0.5,0.5,0.5,0.5),2), maxit = 500, burnin = 3000, n_ab0 = nu_z)

  eryn = 10
  alp_est = colMeans(res$alp[(1:(length(res$nab)/eryn))*eryn,])

  #xbar_est= res$x_bar
  return(alp_est)

}
