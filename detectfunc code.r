library(bbmle)
library(MASS)
library(corpcor) #for make.positive.definite

setClass("dfmod", representation(type="character", model="mle2", formulae="list"))

#Trapezoidal rule of integration
integrate2 <- function(f,lower,upper,...,subs=100)
{ h <- (upper-lower)/subs
  sq <- lower+(1:(subs-1))*h
  (upper-lower) * (f(lower,...)+2*sum(f(sq,...))+f(upper,...)) / (2*subs)
}
#Cosine expansion terms
#x: vector of data 
#prm: vector of coefficients
expan <- function(x, prm)
{  n <- length(x)
   nterms <- length(prm)
   m <- cos(pi*rep(1:nterms,each=n)*rep(x,nterms))
   m <- matrix(m*rep(prm,each=n), nrow=n)
   1+apply(m,1,sum)
}
#Key function, either (if b parameter present in prm) hazard-rate or (otherwise) half-normal
#x: vector of data 
#prm: named list of coefficients, must contain s, plus b if hazard-rate intended
keyfunc <- function(x, prm)
{	if("b" %in% names(prm))
    if(!is.na(prm$b)) 1-exp(-(prm$s/x)^prm$b) else
      exp(-x^2/(2*prm$s^2)) else
      exp(-x^2/(2*prm$s^2)) 
    
}
#Detection function
#x: vector of data 
#type: one of "angle" (line-type) or "radius" (point-type)
#kprm: named list of key function coefficients (see keyfunc)
#eprm: vector of expansion term coefficients (see expan)
#mprm: named list of parameters defining initial increase in detection function, 
#      must contain d (rate) and e (position)
DF <- function(x,mx,type,kprm,eprm=NULL,mprm=NULL)
{  res <- keyfunc(x,kprm)
   if(!is.null(eprm))
     if(sum(is.na(eprm))<length(eprm))
   {  eprm <- eprm[!is.na(eprm)]
      res <- res * expan(x/mx,eprm)
      res[res<0] <- 0
      res <- res / (keyfunc(0,kprm) * expan(0,eprm))
   }
   if(!is.null(mprm))
     if(!is.na(mprm[[1]])) res <- res/(1 + exp(mprm$d*(mprm$e-x)))
   if(type=="radius") res <- res*x
   res
}

#Detection probability density function
#x: vector of data 
#type: one of "angle" (line-type DF) or "radius" (point-type DF)
#lns: log s, width parameter
#lnb: log b, hazard rate parameter
#lnd: log d, increase rate parameter
#e: increase location parameter
#c1,c2: expansion term parameters
#log: logical, return PDF or log(PDF)
PDF <- function(x, type, lns, lnb=NULL, lnd=NULL, e=NULL, c1=NULL, c2=NULL, mx=max(x), log=FALSE)
{
  n <- length(x)
  ns <- length(lns)
  ne <- if(is.null(e)) 1 else length(e)
  np <- ns*ne
  if(is.null(lnb))
  { f <- function(i,s,e) 
      integrate2(DF, 0, mx, mx, type, 
                 kprm=list(s=s[i]), 
                 eprm=c(c1,c2), 
                 mprm=if(is.null(e)) NULL else list(d=exp(lnd), e=e[i]))
    kprm <- list(s=rep(exp(lns),len=n))
  }else
  { f <- function(i,s,e) 
      integrate2(DF, 0, mx, mx, type, 
                 kprm=list(s=s[i], b=exp(lnb)), 
                 eprm=c(c1,c2), 
                 mprm=if(is.null(e)) NULL else list(d=exp(lnd), e=e[i]))
    kprm <- list(s=rep(exp(lns),len=n), b=exp(lnb))
  }
  intgrl <- sapply(1:np, f, rep(exp(lns),ne), rep(e,ns))
  mprm <- if(is.null(e)) NULL else list(d=exp(lnd),e=rep(e,len=n))
  res <- DF(x, mx, type, kprm, c(c1,c2), mprm=mprm) / rep(intgrl,len=n) 
  res[res<1e-322] <- 1e-322
  if(log) res <- log(res)
  res
}
#Fit a detection function
# type: whether to analyse angle or radius data
# key: key function to use (normal or hazard rate)
# form: function form (monotonic or increasing)
# order: integer number of cosine expansion terms in 0:2
# f: list of formulae defining covariates for lns and/or e, eg list(lns~mass,e~mass)
# data: obligatory dataframe containing at least a columns with the same name as type
#       plus columns for any covariates named in f
#...: additional argumants to pass to mle2
# (NB minuslogl, start and data defined internally and cannot be reset; parameters can be only if f NULL)
fitdf <- function(type=c("angle","radius"), key=c("normal","hazard"), 
                  form=c("monotonic","increasing"), order=0, f=NULL, data, ...)
{ type <- match.arg(type)
  key <- match.arg(key)
  form <-match.arg(form)
  if(!order %in% 0:2) stop("order must be integer in 0:2")
  if(!type %in% names(data)) stop("data must contain a column named same as type")
  if(!is.null(f)) 
  { if(!is.list(f) | (is.list(f) &
      sum(unlist(lapply(lapply(f, class), function(x) x=="formula"))) < length(f)))
        stop("when provided, f must be a list of formulae")
    depvars <- unlist(lapply(f, function(x) all.vars(x)[1]))
    covars <- unlist(lapply(f, function(x) all.vars(x)[-1]))
    if(sum(depvars %in% c("lns","e")) < length(depvars))
      stop("formula left-hand side(s) must be lns and/or e")
    if(sum(covars %in% names(data)) < length(covars))
      stop("formula right-hand sides must all be named in data")
    if("e" %in% depvars & form=="monotonic")
      stop("attempting to provide a linear model for e with monotonic form")
  }

  varnms <- c(type, unique(unlist(lapply(f, function(x) attr(terms(x), "term.labels")))))
  data <- data.frame(data[,varnms])
  data <- data.frame(data[apply(!is.na(data),1,prod)==1,])
  names(data) <- varnms
  
  PDFprms <- paste("'",type,"',lns", sep="")
  if(key=="hazard") PDFprms <- paste(PDFprms,"lnb",sep=",")
  if(form=="increasing") PDFprms <- paste(PDFprms,"lnd=lnd","e=e",sep=",")
  if(order==1) PDFprms <- paste(PDFprms,"c1=c1",sep=",")
  if(order==2) PDFprms <- paste(PDFprms,"c1=c1","c2=c2", sep=",")
  formula <- formula(paste(type," ~ PDF(",PDFprms,")"))

  inits <- if(type=="angle") log(2*sd(data$angle)) else log(2*sd(data$radius))
  start <- list(lns=inits)
  if(key=="hazard") start <- c(start, lnb=1)
  if(form=="increasing") start <- c(start, lnd=2, e=0)
  if(order==1) start <- c(start, c1=1)
  if(order==2) start <- c(start, c1=1, c2=-1)

  mod <- mle2(formula, start=start, data=data, parameters=f, ...)
  if(is.null(f) | length(unlist(lapply(f, function(x) attr(terms(x), "term.labels"))))==0)
    f <- list(NULL)
  new("dfmod", type=type, model=mod, formulae=f)
}
#Estimate effective detection parameter an standard error
#If lns or e modelled with covariates, a named list of those covariates 
#must be supplied at which to estimate ED
edest <- function(mod, covars=NULL, reps=1000)
{
  #Function generates a predictor-by-parameter matrix of parameters from a fitted DF object mod
  #If model includes covariates, a named list of covars is required at which to predict
  predict.prms <- function(cfs, covars=NULL)
  { 
    #Function sums coefficients for parameters with linear models, 
    # based on formula list f as originally supplied to fitdF
    #cf: coefficients for the dependent variable
    #cvnms / cvnm: predictor variable names
    #nm: predictor variable name plus "."
    #cv: covariate values for the given predictor variable
    sumcoefs <- function(f)
    { cf <- cfs[grep(paste(f[[2]],".",sep=""),cfnms, fixed=T)]
      res <- cf[grep("Intercept", names(cf))]
      cvnms <- attr(terms(f), "term.labels")
      for(cvnm in cvnms)
      { nm <- paste(".",cvnm, sep="")
        cv <- covars[[match(cvnm,names(covars))]]
        if(is.factor(mod@model@data[[cvnm]]))
        { if(sum(cv %in% levels(mod@model@data[[cvnm]]))<length(cv))
          (stop("factor level(s) provided in covars not found in data"))
          nm <- paste(nm, cv, sep="")
          j <- sapply(nm, grep, names(cf))
          res <- res + unlist(lapply(j, function(j) if(length(j)==0) 0 else cf[j]))
        }else
          res <- res + cv*cf[grep(nm,names(cf))]
      }
      res
    }
    
    cfnms <- names(cfs)
    ncfs <- cfs[match(c("lns","lnb","lnd","e","c1","c2"),cfnms)]
    if(!is.null(mod@formulae[[1]]))
    { ncfs <- matrix(rep(ncfs, each=length(covars[[1]])), ncol=6)
      fcfs <- lapply(mod@formulae, sumcoefs)
      depvars <- lapply(mod@formulae, function(f) f[[2]])
      smatch <- match("lns",depvars)
      ematch <- match("e",depvars)
      if(!is.na(smatch) & length(fcfs[[smatch]])>0) ncfs[,1] <- fcfs[[smatch]]
      if(!is.na(ematch) & length(fcfs[[ematch]])>0) ncfs[,4] <- fcfs[[ematch]]
    } else
    ncfs <- matrix(ncfs,nrow=1)  
    dimnames(ncfs)[[2]] <- c("lns","lnb","lnd","e","c1","c2")
    ncfs
  }
  
  eddcalc <- function(prm)
  { minf <- function(EDD)
    {  if(type=="angle") x <- EDD else x <- EDD^2/2
       sumIN <- x - integrate2(DF,0,EDD, mx,type,kprm,eprm,mprm)
       sumOUT <- integrate2(DF,EDD,mx, mx,type,kprm,eprm,mprm)
       (sumIN - sumOUT)^2
    }
    kprm <- list(s=exp(prm[1]), b=exp(prm[2]))
    mprm <- list(d=exp(prm[3]), e=prm[4])
    eprm <- c(c1=prm[5], c2=prm[6])
    optimise(minf, c(0,mx))$minimum
  }
  
  if(!is.null(mod@formulae[[1]]))
  { if(is.null(covars)) stop("model requires covars (named list of covariate values at which to predict)")
    cvnms <- unlist(lapply(mod@formulae, function(x) attr(terms(x), "term.labels")))
    if(sum(cvnms %in% names(covars)) < length(cvnms)) stop("covars doesn't contain values for all the necessary covariates")
    if(sum(names(covars) %in% cvnms) < length(covars)) warning("some covars not among model covariates so ignored")
    if(length(covars)>1)
    { ns <- unlist(lapply(covars, length))
      if(min(ns)!=max(ns)) stop("elements of covars must have the same length")
    }
  } else
  if(!is.null(covars))
  { warning("covars provided to model without covariates so ignored")
    covars <- NULL
  }

  type <- mod@type
  mx <- max(mod@model@data[[type]])
  ed <- apply(predict.prms(coef(mod@model),covars), 1, eddcalc) #effective detection estimate(s)
  vc <- vcov(mod@model)
  if(det(vc)<=0) vc <- make.positive.definite(vc)
  prmsamp <- mvrnorm(reps, coef(mod@model), vc) #bootstrapped ed estimates
  pmat <- apply(prmsamp, 1, predict.prms, covars)
  if(is.null(covars))
  { edsamp <- apply(pmat, 2, eddcalc)
    cls <- quantile(edsamp, c(0.025,0.975))
    names(cls) <- c("lowerCL","upperCL")
    if(type=="angle")
      2*c(angle=ed, se=sd(edsamp), cls) else
      c(radius=ed, se=sd(edsamp), cls) 
  }else
  { dim(pmat) <- c(length(covars[[1]]),6,reps)
    edsamp <- apply(pmat, 3, function(smat) apply(smat, 1, eddcalc))
    if(is.matrix(edsamp))
    { ses <- apply(edsamp,1,sd)
      cls <- apply(edsamp,1,quantile,c(0.025,0.975))
      dimnames(cls)[[1]] <- c("lowerCL","upperCL")
    } else
    { ses <- sd(edsamp)
      cls <- quantile(edsamp, c(0.025,0.975))
      names(cls) <- c("lowerCL","upperCL")
    }
    if(type=="angle")
      data.frame(covars, angle=2*ed, se=2*ses, t(2*cls)) else
      data.frame(covars, radius=ed, se=ses, t(cls))
  }
}
setMethod("plot", "dfmod",
          function(x, hcol=1, bars=10, add=FALSE, ...)
          {	if(!is.null(x@formulae[[1]])) stop("Cannot (currently) plot covariate models")
            type <- x@type
            dat <- x@model@data[[type]]
            cfs <- coef(x@model)
            nms <- names(cfs)
            xx <- seq(0,max(dat),max(dat)/1000)
            brks <- seq(0, max(dat), length.out=bars+1)
            PDF <- PDF(xx, type, cfs[match("lns", nms)], cfs[match("lnb", nms)],
                      cfs[match("lnd", nms)], cfs[match("e", nms)],
                      cfs[match("c1", nms)], cfs[match("c2", nms)])
            if(add) lines(xx,pdf,...) else 
		if(type=="angle"){
			Angle <- xx
			plot(Angle, PDF, type="l", ...)
		} else 
		{
			Radius <- xx
			plot(Radius, PDF, type="l", ...)
		}			
            if(bars>0 & add==FALSE)
            { hh <- hist(dat,breaks=brks,plot=F)
              stps <- c(hh$density,0)
              lines(brks, stps, type="s", col=hcol)
            }
          }
)

