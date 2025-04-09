expectile_peri_ls <- function(y, f, tau, n.cores=1) {
  ep.parallel1<-function(yy,ff,tau,tt) {
    ntau = length(tau)
    coef = matrix(0, 2, ntau)
    for(i in 1:ntau){
      taui = tau[i]
      x1 = cos(2*pi*ff*tt)
      x2 = sin(2*pi*ff*tt)
      o = rep(1,length(x1))
      X = cbind(o,x1,x2)
      fit <- expectreg.ls(yy ~ X[,2]+ X[,3], mstop=200, expectiles=taui, quietly=TRUE)
      coef[,i] = c(fit$coefficients[[1]],fit$coefficients[[2]])
    }
    return(coef)
  }
  
  ns<-length(y)
  nf<-length(f)
  tt<-c(1:ns)
  ntau<-length(tau)
  
  yy=y
  tmp=list()
  for(i in 1:length(f)){
    tmp[[i]] = ep.parallel1(yy,f[i],tau,tt)
  }
  
  out<-lapply(tmp,FUN=function(x) {apply(x^2,2,sum)})
  out<-matrix(unlist(out),ncol=ntau,byrow=T)
  out<-out*ns/4
  
  out[out<0]<-0
  if(ntau==1) out<-c(out)
  out
}

spec.normalize<-function(qper) {
  if(is.matrix(qper)) {
    return( apply(qper,2,FUN=function(x) { x/sum(x) }) )
  } else {
    return( qper/sum(qper) )
  }
}

