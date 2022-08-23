
#function

#' @title MVMR Mode: A plurality valid sensitivity analysis for multivariable MR
#' @description This function implements a plurality valid (e.g. mode based) estimator for MVMR
#' @param Bexp SNP-exposure associations.
#' @param Bout SNP-outcome associations.
#' @param SEexp Standard error of the SNP-exposure associations.
#' @param SEout Standard error of the SNP-outcome associations
#' @param Psi tuning parameter for the contamination mixture (CM) method (see the original method's documentation for more details)
#' @param CIMin The lower bound on the 95% CI to search for when using CM
#' @param CIMax The upper bound of the 95% CI to search for when using CM
#' @param CIStep The space between each step in CM (see the original method's documentation for more details)
#' @param alpha The p-value threshold
#' @param residual the method used to create residuals. Can be either "IVW", "Egger", or  "Robust". N.B. the method has only been validated for the defult (IVW)
#' @param  Mode The type of plurality valid  estimator to use. Can be "CM" for contamination mixture, or "MBE" for a traditional mode based estimator. n.b. the default (CM) performed best in the validation
#' @param  weighting For the MBE, if to use a weighed MBE or weighted.
#' @param  stderror For MBE, the type of standard error. Please see MendelianRandomization R package documentation for more detail
#' @param  phi tuning parameter for the MBE. Please see MendelianRandomization R package documentation for more detail
#'  @param distribution distribution assumed for the MBE. Please see MendelianRandomization R package documentation for more detail
#'  @param iterations the number of iterations to use for the SE of the MBE
#' @keywords Two-sample Mendelian Randomization, 2SMR, MVMR, Multivaribile MR
#' @export




mv_mrmode<-function(Bout, Bexp, SEout, SEexp, Psi=0, CIMin = NA, CIMax = NA, CIStep = 0.001, alpha = 0.05, residual="IVW",Mode="CM",weighting = "weighted",stderror = "simple",phi = 1,distribution = "normal",iterations = 10000) {

  #Function to run CM using residuals for two exposure
  bvmodeR<-function(Bout, Bexp1, Bexp2, SEout, SEexp1, SEexp2,residual){
    if (residual=="IVW"){
      e1<- summary( lm(Bexp1~Bexp2+0,weights = 1/SEout^2))$res
      o1<-  summary(lm(Bout~Bexp2+0,weights = 1/SEout^2))$res
    }else if(residual=="Egger"){
      e1<- summary( lm(Bexp1~Bexp2,weights = 1/SEout^2))$res
      o1<-  summary(lm(Bout~Bexp2,weights = 1/SEout^2))$res
    }else if(residual=="Robust"){
      e1<- summary( robustbase::lmrob(Bexp1~Bexp2-1,weights = 1/SEout^2,k.max=500, maxit.scale = 500))$res
      o1<- summary(robustbase::lmrob(Bout~Bexp2-1,weights = 1/SEout^2,k.max=500, maxit.scale = 500))$res
      #}else if(residual=="Median"){
      #  e1<-as.numeric( summary( quantreg::rq(Bexp1~Bexp2-1,weights = 1/SEout^2))$residuals)
      #  o1<- as.numeric(  summary(quantreg::rq(Bout~Bexp2-1,weights = 1/SEout^2))$residuals)
    }else {return(warnings("residual must be set to: IVW, Egger, Robust"))} #, or Median

    if (Mode == "CM"){
      out<-MendelianRandomization::mr_conmix(MendelianRandomization::mr_input(bx = e1, bxse = SEexp1, by = o1,byse = SEout/SEout), psi = Psi, CIMin = CIMin, CIMax = CIMax, CIStep = CIStep, alpha = alpha)
    } else if (Mode == "MBE"){
      out<-MendelianRandomization::mr_mbe(MendelianRandomization::mr_input(bx = e1, bxse = SEexp1, by = o1,byse = SEout/SEout),  weighting = weighting, stderror = stderror,  phi = phi, seed = 314159265, iterations = iterations, distribution = distribution,  alpha = alpha)
    } else {return(warnings("Mode must take the value of CM or MBE"))}

    return(out)}

  #code for running if there is only one exposure
  if (length(Bout)==length(Bexp)){
    if (Mode == "CM"){
      out_t<-MendelianRandomization::mr_conmix(MendelianRandomization::mr_input(bx = Bexp, bxse = SEexp, by = Bout,byse = SEout), psi = Psi, CIMin = CIMin, CIMax = CIMax, CIStep = CIStep, alpha = alpha)
    } else if (Mode == "MBE"){
      out<-MendelianRandomization::mr_mbe(MendelianRandomization::mr_input(bx = Bexp, bxse = SEexp, by = Bout,byse = SEout),  weighting = weighting, stderror = stderror,  phi = phi, seed = 314159265, iterations = iterations, distribution = distribution,  alpha = alpha)
    } else {warnings("Mode must take the value of CM or MBE")}



    #code for running mode function with two or more exposures
  }else{
    Exposure<-1:ncol(Bexp)
    out_t<-data.frame(Exposure)
    for (i in 1:ncol(Bexp)){
      mr1<-bvmodeR(Bout,Bexp1=Bexp[,i],Bexp2=Bexp[,-i],SEout,SEexp1=SEexp[,i],SEexp2=SEexp[,-i],residual=residual )
      out_t$Estimate[out_t$Exposure==i]<-mr1$Estimate
      if (Mode == "CM"){
        out_t$CILower[out_t$Exposure==i]<-min(mr1$CILower)
        out_t$CIUpper[out_t$Exposure==i]<-max(mr1$CIUpper )
      } else if (Mode == "MBE"){
        out_t$StdError[out_t$Exposure==i]<-mr1$StdError
        out_t$CILower[out_t$Exposure==i]<-mr1$CILower
        out_t$CIUpper[out_t$Exposure==i]<-mr1$CIUpper
      }
      out_t$Pvalue[out_t$Exposure==i]<-mr1$Pvalue
      if (Mode == "CM"){
        out_t$NvalidSNPs[out_t$Exposure==i]<-length(mr1$Valid)
        out_t$Nrange[out_t$Exposure==i]<-length(mr1$CIUpper)}}}

  return(out_t)     }

