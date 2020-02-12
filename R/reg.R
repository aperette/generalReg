#' @title  Normal regression model with general parametrization
#'
#' @description Determine the estimates of the parameters of a
#' general parametrization model with a Newthon Raphson method
#'
#'
#' @param formula a nonlinear model formula including variables and parameters
#' @param formula_var a nonlinear model formula for the diagonal of covariance matrix
#' @param data A data frame in which to evaluate the variables in \code{formula} and \code{formula_var}.
#' Can also be a list or an environment, but not a matrix
#' @param start A named list of starting estimates. When \code{start} is missing, a very cheap guess for \code{start} is tried and parameters names are automatically identified
#' @param control A list of control parameters. See 'Details'
#' @param bias_correction Logical. Should a bias correction be estimated for the parameters?
#'
#'
#' @details The control argument is a list that can supply any of the following components:
#'  \itemize{
#'  \item control$max_it Maximum number of iterations
#'  \item control$reltol Maximum difference between iterations to consider convergence
#'}
#'
#'  @examples
#' \dontrun{
#'
#' fit <- reg_general(y~alfa+X1^(gama)+beta*X2,~100+sigma*X3,data=data,
#' start=list(alfa=100,beta=1,gama=0.5,sigma=0.1),bias_correction = T)
#'
#' }
#'
#' @export


reg_general=function(formula=NULL,
                   formula_var=NULL,
                   data,
                   start=NULL,
                   method="NR",
                   control=list(reltol=0.0001,max_it=500,kappa="AUTO",verbose=0),
                   bias_correction=T){

  if(is.null(formula)){
    stop('arguments "formula" is missing, with no default',call. = F)}
  if(is.null(formula_var)){
    stop('arguments "formula_var" is missing, with no default',call. = F)}
  if(!inherits(formula,"formula"))
      stop("formula has a indefinite format",call. = F)
  if(!inherits(formula_var,"formula"))
      stop("formula_var has a indefinite format",call. = F)
  if(!inherits(control,"list"))
    stop("control must be a list",call. = F)
  if(!is.null(start)){
    if(!is.list(start))
      stop("start must be a list",call. = F)
    theta=names(start)}
  if(is.numeric(control$kappa))
    if(control$kappa<=0 | control$kappa>1)
      stop("kappa must be in interval (0,1]",call. = F)
  if(is.character(control$kappa))
    if(control$kappa!="AUTO")
      stop('kappa must be numeric or "AUTO"',call. = F)
  if(!is.na(method))
    if(!method %in% c("pso","optim","NR","gamlss"))
      stop('method should be one of "NR", "pso", "optim","gamlss"',call. = F)
  class0=class(data)
  if(class0!="data.frame"){
    data = tryCatch(expr=as.data.frame(data), error=function(e) NA)
    if(class(data)!="data.frame")
      stop("Cannot coerce class '", class0,"' to a data.frame",call. = F)}
  if(is.na(method)) method="NONE"
  if(is.null(control$reltol)) control$reltol=0.0001
  if(is.null(control$max_it)) control$max_it=500
  if(is.null(control$kappa)) control$kappa=1
  if(is.null(control$verbose)) control$verbose=0

  inputs = list(formula=formula,formula_var=formula_var,control=control,start=start,bias_correction=bias_correction)

  resposta=data[,as.character(formula)[2]]
  mu=as.character(formula)[3]
  S=tail(as.character(formula_var),1)

  cl <- match.call()
  n=nrow(data)
  q=1
  par=definir_parametros(c(mu,S),data)
  covar=par[[1]]
  data = data %>% dplyr::select_(.dots=c(as.character(formula)[2],"covar"))
  if(is.null(start)){
    theta=par[[2]]
  }else{
    out = par[[2]][!par[[2]] %in% theta]
    if(length(out)>0)
      stop(paste0(out[1]," not found in start"),call. = F)
    par[[2]]=theta
  }

  if(control$verbose!=0) cat("Parameters:",paste0(theta,collapse=", "),"\n")
  if(control$verbose!=0) cat("Predictors:",paste0(covar,collapse=","),"\n")

  l=paste0("-0.5*log(",S,") -0.5*(",mu,"-resposta)^2/(",S,")")

  Di=array()
  Vi=array()
  Pd=array()
  ld=array()
  a=matrix(NA,nrow=length(theta),ncol=length(theta))
  C=matrix(NA,nrow=length(theta),ncol=length(theta))
  for(k in 1:length(theta)){
    Di[k]=as.character(Deriv::Deriv(mu,theta[k]))
    Vi[k]=as.character(Deriv::Deriv(S,theta[k]))
    ld[k]=as.character(Deriv::Deriv(l,theta[k]))
    Pd[k]=as.character(Deriv::Deriv(paste0("((",S,")^0.5)"),theta[k]))
    for(j in 1:length(theta)){
      a[k,j]=as.character(Deriv::Deriv(Di[k],theta[j]))
      C[k,j]=as.character(Deriv::Deriv(Vi[k],theta[j]))}
  }

  Fi=rbind(Di,Vi)
  aC=rbind(c(a),c(C))

  eval(parse(text=paste0(
    "gerar_l = function(theta, data,resposta){
    N=nrow(data)
    -0.5*sum(log(rep(1,N)*",estruturar(S,par,asArray=T),") +
    ((rep(1,N)*",estruturar(mu,par,asArray=T),")-resposta)^2/(rep(1,N)*",estruturar(S,par,asArray=T),"))}"
    )))

  eval(parse(text=paste0(
    "log_like=function(theta, data){
    sum(dnorm(c(",paste(resposta,collapse=","),"),",
    paste("rep(1,",length(resposta),")*",estruturar(mu,par,asArray=T),collapse = ","),
    ",sqrt(",paste("rep(1,",length(resposta),")*",estruturar(S,par,asArray=T),collapse = ","),"),log = T))
    }"
  )))


  eval(parse(text=paste0(
    "gerar_mu = function(theta, data){
    N=nrow(data)
    c(",
    paste("rep(1,N)*",estruturar(mu,par),collapse = ",")
    ,")}"
  )))

  eval(parse(text=paste0(
    "gerar_sigma = function(theta, data){
    N=nrow(data)
    Matrix::.sparseDiagonal(x=c(",
    paste("rep(1,N)*",estruturar(S,par),collapse = ","),
    "),n=N)}"
  )))

  eval(parse(text=paste0(
    "gerar_D = function(theta, data){
    N=nrow(data)
    Matrix::Matrix(c(",
    paste("rep(1,N)*",estruturar(t(Di),par),collapse = ","),
    "), ncol=N,byrow=T) %>% as.vector %>% Matrix::Matrix(ncol=",
    length(theta),", sparse=T,byrow=T)}"
  )))

  eval(parse(text=paste0(
    "gerar_V = function(theta, data){
    N=nrow(data)
    Matrix::Matrix(c(",
    paste("rep(1,N)*",estruturar(t(Vi),par),collapse = ","),
    "), ncol=N,byrow=T) %>% as.vector %>% Matrix::Matrix(ncol=",
    length(theta),", sparse=T,byrow=T)}"
  )))

  gerar_J = function(D){
    dt=dplyr::tibble(valor=as.vector(D),
                     id=rep(1:nrow(D),ncol(D)),
                     par=rep(1:ncol(D),each=nrow(D)))
    dt2=dplyr::left_join(dt,dt,by="id") %>%
      dplyr::mutate(total=2*valor.x*valor.y,
                    aux=1) %>%
      dplyr::select(-valor.x,-valor.y)
    dplyr::bind_rows(dt2,dt2 %>% dplyr::mutate(aux=0, total=0)) %>%
      dplyr::arrange(par.x,par.y,id,aux) %>%
      .$total %>%
      Matrix::Matrix(nrow=2*nrow(D),sparse=T)
  }

  eval(parse(text=paste0(
    "gerar_G = function(theta, data){
    N=nrow(data)
    Matrix::Matrix(c(",
    paste("rep(1,N)*",estruturar(t(aC),par),collapse = ","),
    "), ncol=N,byrow=T) %>% as.vector %>% Matrix::Matrix(ncol=",
    ncol(aC),", sparse=T,byrow=T)}"
  )))

  eval(parse(text=paste0(
    "gerar_F = function(theta, data){
    N=nrow(data)
    Matrix::Matrix(c(",
    paste("rep(1,N)*",estruturar(Matrix::t(Fi),par),collapse = ","),
    "), ncol=N,byrow=T) %>% as.vector %>% Matrix::Matrix(ncol=",
    ncol(Fi),", sparse=T,byrow=T)}"
  )))

  correcao_vies=function(){

    D=gerar_D(theta_val,data)
    G=gerar_G(theta_val,data)
    J=gerar_J(D)

    phi=-0.5*(G+J)

    K=solve(Matrix::t(Fn)%*%Hn%*%Fn,tol=1e-2000)
    e=phi%*%as.vector(K)

    K%*%Matrix::t(Fn)%*%Hn%*%e}

  gerar_H = function(sigma){
    Matrix::diag(sigma) %>%
    {rbind(1/.,0.5*(.^-2))} %>%
      as.vector() %>%
      {Matrix::.sparseDiagonal(x=.,n=length(.))}
  }

  gerar_s=function(media,sigma,theta,p=control$kappa,all=T){
    res=resposta-media
    s2=-(Matrix::diag(sigma)-res^2)
    s=c(rbind(res,s2))
    if(all==F | p=="AUTO") return(s)
    S=Fn %*% unlist(theta) + s*p
    return(S)}

  #skovgaard:
  eval(parse(text=paste0(
    "gerar_Pd = function(theta, data){
    N=nrow(data)
    Matrix::Matrix(c(",
    paste("rep(1,N)*",estruturar(Pd,par),collapse = ",")
    ,"), nrow=N, sparse=T)}"
  )))

  eval(parse(text=paste0(
    "gerar_C = function(theta, data){
    N=nrow(data)
    Matrix::Matrix(c(",
    paste("rep(1,N)*",estruturar(c(C),par),collapse = ",")
    ,"), nrow=N, sparse=T)}"
  )))

  eval(parse(text=paste0(
    "gerar_D2 = function(theta, data){
    N=nrow(data)
    Matrix::Matrix(c(",
    paste("rep(1,N)*",estruturar(c(a),par),collapse = ",")
    ,"), nrow=N, sparse=T)}"
  )))

  gerar_M2=function(media,sigma){
    z=resposta-media
    c(-Matrix::diag(sigma),
           -2*z*Matrix::diag(sigma),
           -2*z*Matrix::diag(sigma),
           -4*(z^2)*Matrix::diag(sigma)+2*(Matrix::diag(sigma))^2) %>%
      matrix(ncol=length(z),byrow = TRUE) %>%
      apply(2,function(x) Matrix::Matrix(x,ncol=2,byrow=TRUE)) %>%
      Matrix::bdiag()
  }


  loglike_sim=array()
  erros=array()
  interacao=list()
  dif=NA
  if(!is.null(start))
    start= unlist(start,use.names = F)
  if(is.null(start))
    start= rep(1,length(theta))
  theta_val = data.frame(matrix(start,nrow=1,dimnames = list(NULL,theta)))

  if(method=="pso"){
    if(control$verbose!=0) cat("In pso, start values are not used\n")
    logit=function(p){
      p=log(p/(1-p))
      ifelse(is.infinite(p),NA,p)}

    opt = pso::psoptim(rep(NA,length(theta)),
                        fn=function(par){
                          aux=-log_like(logit(par),data)
                          if(is.na(aux)) aux=Inf
                          return(aux)},
                       upper=1,lower=0,
                        control=list(maxit=control$max_it,max.restart=1,reltol=control$reltol,
                                     trace=control$verbose,REPORT=control$verbose))
    k=opt$counts[2]
    if(opt$convergence==2)
      stop("Maximal number of iterations reached",call. = F)
    theta_val_novo = logit(opt$par)
    names(theta_val_novo)=theta
    theta_val=data.frame(valor=logit(opt$par),nome=theta) %>% tidyr::spread(nome,valor)
  }
  if(method=="gamlss"){

    mu_names = theta[stringr::str_detect(mu,theta)]
    sigma_names = theta[stringr::str_detect(S,theta)]
    suppressWarnings(expr={
      require(gamlss.dist,quietly = T)
      fit = gamlss.nl::nlgamlss(y=eval(parse(text=as.character(formula)[2])),
                                mu.formula=formula,
                                sigma.formula=formula_var,
                                mu.start = start[stringr::str_detect(mu,theta)],
                                sigma.start = start[stringr::str_detect(S,theta)],
                                control=gamlss.nl::NL.control(iterlim=control$max_it,steptol=control$reltol,
                                                             gradtol = control$reltol),
                                llik.output = control$verbose>0,
                                family=gamlss.dist::gamlss.family(NO2(sigma.link="identity")),
                                data=data)
    })
    k=fit$iter
    theta_val_novo = fit$coefficients
    names(theta_val_novo)=theta
    theta_val=data.frame(valor=fit$coefficients,nome=theta) %>% tidyr::spread(nome,valor)
  }
  if(method=="optim"){
    opt = optim(unlist(theta_val),
                fn=function(par){
                  aux=-log_like(par,data)
                  if(is.na(aux)) aux=Inf
                  return(aux)},
                control=list(maxit=control$max_it,reltol=control$reltol,
                             trace=1*(control$verbose>0),REPORT=control$verbose))
    k=opt$counts[1]
    if(opt$convergence==1)
      stop("Maximal number of iterations reached",call. = F)
    theta_val_novo = opt$par
    names(theta_val_novo)=theta
    theta_val=data.frame(valor=opt$par,nome=theta) %>% tidyr::spread(nome,valor)
  }
  if(method=="NR"){
    k=1
    dif=1
  while(dif>control$reltol & k <= control$max_it){

    sigma=gerar_sigma(theta_val,data)
    media=gerar_mu(theta_val,data)
    Fn=gerar_F(theta_val,data)
    if(sum(is.na(Fn))>0) stop("Missing values in matrix F",call. = F)
    Hn=gerar_H(sigma)
    score = Matrix::t(Fn)%*%Hn%*%gerar_s(media,sigma,theta_val,all=F)
    fisher = Matrix::t(Fn)%*%Hn%*%Fn
    inv_fisher= tryCatch(expr={Matrix::solve(fisher,tol=1e-2000)},
                         error=function(e)NA)
    loglike_sim[k]=sum(dnorm(resposta,media,sqrt(Matrix::diag(sigma)),log = T))

    if(is.numeric(control$kappa)){
      kappa=control$kappa
      Sn=gerar_s(media,sigma,theta_val)
      # theta_val_novo = tryCatch(expr={(Matrix::solve(Matrix::t(Fn)%*%Hn%*%Fn,tol=1e-2000) %*% Matrix::t(Fn)%*%Hn%*%Sn)[,1]},
      #                         error=function(e)NA)
      theta_val_novo = tryCatch(expr={unlist(theta_val) + control$kappa*(inv_fisher %*% score)[,1]},
                              error=function(e)NA)}
    if(control$kappa=="AUTO"){
      theta_atual=theta_val
      suppressWarnings(expr={
        kappa=optim(0.5,function(par){
          # teste = tryCatch(expr={
          #   (Matrix::solve(Matrix::t(Fn)%*%Hn%*%Fn,tol=1e-2000) %*% Matrix::t(Fn)%*%Hn%*%gerar_s(media,sigma,theta_val,p = par[1]))[,1]},
          teste = tryCatch(expr={
            (unlist(theta_atual) + par[1]*(inv_fisher %*% score))[,1]},error=function(e)NA)
          theta_val[1,]=teste
          -sum(dnorm(resposta,gerar_mu(theta_val,data),sqrt(Matrix::diag(gerar_sigma(theta_val,data))),log = T))},
          method="Brent",lower=0,upper=1)$par
        })
      theta_val_novo = tryCatch(expr={
        (Matrix::solve(Matrix::t(Fn)%*%Hn%*%Fn,tol=1e-2000) %*% Matrix::t(Fn)%*%Hn%*%gerar_s(media,sigma,theta_val,p = kappa))[,1]},
        error=function(e)NA)
      }


    if(is.na(theta_val_novo)[1]) stop("Non convergence",call. = F)
    names(theta_val_novo) = theta

    interacao[[k]]=theta_val_novo
    k=k+1
    dif = max(abs(theta_val-theta_val_novo)/abs(theta_val))
    dif = ifelse(is.infinite(dif),1,dif)
    dif = min(dif,kappa)
    erros[k-1] = mean((media-resposta)^2)
    theta_val[1,] = theta_val_novo

    if(control$verbose>0 & (k %% control$verbose)==0 & control$kappa=="AUTO")
      cat(paste0("It ",k,
                 ": log-likelihood=",round(loglike_sim[k-1],2),
                 " error_max=",formatC(dif,format = "e", digits = 2),
                 ", kappa=",round(kappa,3),"\n"))

    if(control$verbose>0 & (k %% control$verbose)==0 & is.numeric(control$kappa))
      cat(paste0("It ",k,
                 ": log-likelihood=",round(loglike_sim[k-1],2),
                 " error_max=",formatC(dif,format = "e", digits = 2),"\n"))


    if(is.na(dif) | dif>10000000) stop("Non convergence",call. = F)

  }
  if(k >= control$max_it | is.na(dif)) stop("Maximal number of iterations reached",call. = F)
  k=k-1
  # if(which.max(loglike_sim)!=length(loglike_sim) &
  #    lm(erros~index,data=data.frame(erros,index=1:length(erros)))$coef[2]>0)
  #   warning("Estimate may not be a Maximum Likelihood Estimation",call. = F)
  }

  sigma=gerar_sigma(theta_val,data)
  media=gerar_mu(theta_val,data)
  Fn=gerar_F(theta_val,data)
  Hn=gerar_H(sigma)
  Sn=gerar_s(media,sigma,theta_val)

  #Fisher obs
  G=gerar_G(theta_val,data)
  M2 = gerar_M2(media,sigma)
  obs = Matrix::t(Fn) %*% Hn %*% M2 %*% Hn %*% Fn +
    matrix(t(gerar_s(media,sigma,theta_val,all = F)) %*% Hn %*% G,ncol=length(theta))

  corrigido=NA
  if(bias_correction){
    bias=correcao_vies()
    corrigido=theta_val-bias[,1]
    corrigido=unlist(corrigido)
  }

if(is.null(names(theta_val))) names(theta_val)= theta

  out <- list(
    parameters=unlist(theta_val),
    corrected_parameters=corrigido,
    fitted.values = media,
    residuals = resposta-media,
    var = sigma,
    it = k-1,
    score = Matrix::t(Fn)%*%Hn%*%gerar_s(media,sigma,theta_val,all=F),
    fisher = Matrix::t(Fn)%*%Hn%*%Fn,
    fisher_obs = -1*obs,
    loglike = -sum(dnorm(resposta,media,sqrt(Matrix::diag(sigma)),log = T)),
    target=resposta,
    epsilon=dif,
    theta=theta_val[1,],
    functions=list(
    function_loglike=log_like,
    function_mu = gerar_mu,
    function_sigma = gerar_sigma,
    function_l=gerar_l,
    function_F=gerar_F,
    function_H=gerar_H,
    function_s=gerar_s,
    function_D=gerar_D,
    function_V=gerar_V,
    function_G=gerar_G,
    function_C=gerar_C,
    function_D2=gerar_D2,
    function_M2=gerar_M2),
    inputs = inputs,
    data = data,
    call = cl
  )
  class(out) <- "genReg"

  return(out)
}

#' @export
setClass("genReg")


#' @method print genReg
#' @export
print.genReg <- function(x, ...){
  cat("Call:\n")
  print(x$call)

  ndec=stringr::str_length(stringr::str_extract(as.character(x$parameters), "\\.[0-9]*")) - 1

  nround=(min(5,max(ndec,na.rm=T)))
  cat("\n\nCoefficients:\n")
  print(round(x$parameters,nround))
  if(x$inputs$bias_correction){
  cat("\nCorrected coefficients:\n")
  print(round(x$corrected_parameters,nround))}
  cat("\n")
}

#' @method summary genReg
#' @export
summary.genReg <- function(x, ...){
  cat("\nCall:\n")
  print(x$call)

  ndec=stringr::str_length(stringr::str_extract(as.character(x$parameters), "\\.[0-9]*")) - 1
  nround=(min(5,max(ndec,na.rm=T)))
  cat("\n\nResiduals:\n")
  print(summary(x$res)[-4])
  cat("\n\nCoefficients:\n")

  df=length(x$target)-length(x$parameters)+1
  coef=x$parameters
  se=x$fisher_obs %>% solve %>% diag %>% sqrt
  tv=coef/se
  pv=pt(2*abs(tv),df,lower.tail=F)
  #Pr(>|t|)
  print(data.frame(Estimate=round(coef,nround),Std.Error=round(se,nround), t.value=tv, Pv=pv,
             row.names = names(x$parameters)))
  cat("\n")
}



#' @method coef genReg
#' @export
coef.genReg <- function(x, ...){
  x$parameters
}

#' @method confint genReg
#' @export
confint.genReg <- function(x, level=0.95, bias_correction = "both",...){
  if(level<=0 | level>=1 | !is.numeric(level))
    stop("level must be in interval (0,1)",call. = F)
  if(bias_correction=="both")
    bias_correction=NULL
  flag_null=is.null(bias_correction)
  if(is.null(bias_correction))
    bias_correction=x$inputs$bias_correction

  normal_int=conf_interval(x,x$parameters, level)
  if(bias_correction){
    if(is.na(x$corrected_parameters[1])) stop("Corrected parameters not founded",call. = F)
    bias_int=conf_interval(x,x$corrected_parameters, level)}

  if(bias_correction==FALSE){
    return(normal_int)
  }
  if(bias_correction==TRUE & flag_null==FALSE){
    return(bias_int)
  }
  return(list(EMV_intervals=normal_int,corrected_intervals=bias_int))
}

#' @method vcov genReg
#' @export
vcov.genReg <- function(x, bias_correction="both", tolerance=1e-20,...){
  if(bias_correction=="both")
    bias_correction=NULL
  flag_null=is.null(bias_correction)
  if(is.null(bias_correction))
    bias_correction=x$inputs$bias_correction

  normal_fisher=Matrix::solve(fisher_inf(x,x$parameters),tol = tolerance)
  if(bias_correction){
    if(is.na(x$corrected_parameters[1])) stop("Corrected parameters not founded",call. = F)
    bias_fisher=Matrix::solve(fisher_inf(x,x$corrected_parameters),tol=tolerance)}

  if(bias_correction==FALSE){
    return(normal_fisher)
  }
  if(bias_correction==TRUE & flag_null==FALSE){
    return(bias_fisher)
  }
  return(list(vcov=normal_fisher,corrected_vcov=bias_fisher))

}

#' @method logLik genReg
#' @export
logLik.genReg <- function(x, bias_correction=FALSE, ...){
  if(bias_correction){
    if(is.na(x$corrected_parameters[1])) stop("Corrected parameters not founded",call. = F)
    return(-x$functions$function_loglike(x$corrected_parameters,x$data))}
  return(x$loglike)
}

#' @method predict genReg
#' @export
predict.genReg <- function(x, newdata=NULL, type = "mean", bias_correction = NULL, ...){
  if(is.null(bias_correction)) bias_correction=x$inputs$bias_correction
  if(!type %in% c("mean","var")) stop("type must be 'mean' or 'var'",call. = FALSE)
  parameters=sapply(x$parameters,list)
  if(bias_correction){
    if(is.na(x$corrected_parameters[1])) stop("Corrected parameters not founded",call. = F)
    parameters=sapply(x$corrected_parameters,list)}
  if(type == "mean"){
    if(is.null(newdata)) newdata=x$data
    return(x$functions$function_mu(parameters,newdata))}
  if(type == "var")
    if(is.null(newdata)) newdata=x$data
  return(Matrix::diag(x$functions$function_sigma(parameters,newdata)))
}

#' @export
likelihood_ratio <- function(x, parameters,correction=FALSE,control=NULL,start=NULL){
  if(class(x)!="genReg")
    stop("x must be a genReg class",call. = F)

  par_teste=parameters

  formula=Reduce(paste, deparse(x$inputs$formula))
  formula_var=Reduce(paste, deparse(x$inputs$formula_var))

  theta_names = names(x$parameters)
  equals_par = names(parameters)[names(parameters) %in% names(x$parameters)]
  dif_par = names(x$parameters)[!names(x$parameters) %in% equals_par]
  wrong_par = names(parameters)[!names(parameters) %in% names(x$parameters)]
  if(length(wrong_par)>1) stop(paste0("Parameter ",dif_par[1]," not found"),call. = F)

  if(length(equals_par)<length(x$parameters)){
    for(i in 1:length(equals_par)){
      formula=stringr::str_replace_all(formula,equals_par[i],as.character(parameters[[equals_par[i]]]))
      formula_var=stringr::str_replace_all(formula_var,equals_par[i],as.character(parameters[[equals_par[i]]]))
    }

    if(!is.null(control$reltol)) x$inputs$control$reltol=control$reltol
    if(!is.null(control$max_it)) x$inputs$control$max_it=control$max_it
    if(!is.null(control$kappa)) x$inputs$control$kappa=control$kappa
    x$inputs$start = x$inputs$start[dif_par]
    if(is.null(start)) start=x$inputs$start
    control=x$inputs$control
    control$verbose=0
    cat("MLE under H0:\n")

    x2=reg_general(
      formula=as.formula(formula),
      formula_var=as.formula(formula_var),
      start=x$inputs$start,
      data=x$data,
      control=control)
    cat(paste(names(coef(x2)),round(coef(x2),3),sep=":",collapse="  "),"\n")
    par_teste = coef(x2) %>% data.frame(nome=names(.),valor=.) %>% tidyr::spread(nome,valor) %>%
      merge(parameters) %>%
      dplyr::select(theta_names)
  }
  loglike=function(y,media,var){
    sum(dnorm(y,mean=media,sqrt(var),log=T))}
  theta=coef(x) %>% data.frame(nome=names(.),valor=.) %>% tidyr::spread(nome,valor)
  p=length(theta)
  q=length(parameters)
  data=x$data

  mu=x$functions$function_mu(theta,x$data)
  mu0=x$functions$function_mu(par_teste,x$data)
  S=x$functions$function_sigma(theta,x$data)
  S0=x$functions$function_sigma(par_teste,x$data)
  var=Matrix::diag(S)
  var0=Matrix::diag(S0)
  y=x$target

  l1=loglike(y,mu,var)
  l0=loglike(y,mu0,var0)

  LR=2*(l1-l0)
  if(LR<0) stop("Likelihood ratio test statistic is negative",call. = F)

  pv=pchisq(LR,df = q,lower.tail = F)

  if(correction==FALSE){return(list(LR=LR,p_value=pv))}

  ####Skovgaard correction:
  z=y-mu
  z0=y-mu0
  u0=(z0^2)/var0
  V=x$functions$function_V(theta,data)
  V0=x$functions$function_V(par_teste,data)
  D=x$functions$function_D(theta,data)
  D0=x$functions$function_D(par_teste,data)

  P=sqrt(var)
  P0=var0^0.5
  Pd=0.5*diag(1/sqrt(var))%*%V
  a=z/P
  u00=(a^2)*(P0^2)/var0


  gerar_J=function(theta){
    S_aux=x$functions$function_sigma(theta,x$data)
    var=Matrix::diag(S_aux)
    V_aux=x$functions$function_V(theta,data)
    D_aux=x$functions$function_D(theta,data)
    z=y-x$functions$function_mu(theta,x$data)
    T=D_aux+z*V_aux/var
    B=-z*D_aux -0.5*V_aux
    A=-V_aux/var^2
    aux1=rep(1:length(theta),length(theta))
    aux2=rep(1:length(theta),each=length(theta))
    Ad = -2*A[,aux1]*V_aux[,aux2]/var
    Ad = Ad - x$functions$function_C(theta,data)/var^2
    E = -0.5*(Ad*(var-z^2)) - x$functions$function_D2(theta,data)*z/var
    G = apply(B[,aux1]*A[,aux2] + E,2,sum) %>% matrix(nrow=length(theta))
    Matrix::t(T)%*%solve(S_aux)%*%D_aux + G
  }


  J = gerar_J(theta)
  J0 = gerar_J(par_teste)

  T00=D0+a*P0*V0/var0

  B00=-a*P0*D0 -0.5*V0
  A0=-V0/var0^2

  aux1=rep(1:length(theta),length(theta))
  aux2=rep(1:length(theta),each=length(theta))
  Ad0 = -2*A0[,aux1]*V0[,aux2]/var0
  Ad0 = Ad0 - x$functions$function_C(par_teste,data)/var0^2

  E00 = -0.5*(Ad0*(var0-(a*P0)^2)) - x$functions$function_D2(par_teste,data)*a*P0/var0
  G00 = apply(B00[,aux1]*A0[,aux2] + E00,2,sum) %>% matrix(nrow=length(theta))

  J00 = Matrix::t(T00)%*%solve(S0)%*%D0 + G00
  J00=as.matrix(J00)

  R = Pd*a + D

  Q = D+z*V/var
  Q00 = D0+a*P0*V0/var0

  Ud0 = Matrix::t(Q00)%*%solve(S0)%*%R
  U0= Matrix::t(x$functions$function_F(par_teste,data)) %*% x$functions$function_H(S0) %*% x$functions$function_s(mu0,S0,par_teste,all = F)

  ld = Matrix::t(R)%*%solve(x$var)%*%(-z)
  ld0 = Matrix::t(R)%*%solve(S0)%*%(-z0)

  w=which(!names(theta) %in% names(parameters))

  p1=Matrix::det(J)^0.5*Matrix::det(Ud0)^-1*Matrix::det(J0[w,w])^0.5*Matrix::det(J00[w,w])^-0.5*Matrix::det(J00)^0.5
  p2=(Matrix::t(U0)%*%solve(J00,tol=1e-200) %*% U0)^(p/2)
  p3=LR^(q/2-1)
  p4=Matrix::t(ld-ld0)%*%Matrix::solve(Ud0,tol=1e-20)%*%U0
  rho=p1*p2[1,1]/(p3*p4[1,1])
  LR2=LR-2*log(rho)

  #if((0 + (p1<0) + (p2[1,1]<0) + (p3<0) + (p4[1,1]<0))>0)
  if(rho<0){
    stop("Negative values are founded in Skovgaard's correction",call. = F)}


  if(LR2<0){
    LR2=LR*(1-1/LR*log(rho))^2
    warning("Correção quadrática utilizada",call. = F)
  }

  pv2=tryCatch(pchisq(LR2,df = q,lower.tail = F), error=function(e) NA)


  out = list(LR=LR,p_value=pv, LR_correction=LR2, p_value_correction=pv2)
  class(out) = "lr_ratio"
  return(out)
}

#' @method print lr_ratio
#' @export
print.lr_ratio <- function(x, ...){
  cat("Likelihood ratio test:\n")
  x$LR=paste0("LR:",round(x$LR,3))
  x$p_value=paste0("p-value:",formatC(x$p_value,format = "e", digits = 2))
  cat(c(x$LR,x$p_value),"\n")
  if(length(x)>2){
    cat("Skovgaard correction:\n")
  x$LR_correction=paste0("LR*:",round(x$LR_correction,3))
  x$p_value_correction=paste0("p-value:",formatC(x$p_value_correction,format = "e", digits = 2))
  cat(c(x$LR_correction,x$p_value_correction),"\n")

  }
}
