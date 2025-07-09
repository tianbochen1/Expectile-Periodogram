# simulation
source('ep.R')
library(fGarch)
library(quantreg)
library(expectreg)
ncores = 2
tau=(5:95)/100



#######EEG example
tau=(5:95)/100
library(R.matlab)
data = readMat('Patient-3_seizure-2_Q5-Q6_17-21.mat')
data = data$tmpdata
data = data[1,1000:4000]

n <- length(data)
freq <- c(0:(n-1))/n
freq <- freq[freq>0 & freq <0.1]
ep = expectile_peri_ls(data, freq, tau, n.cores=8)
pg = (1/n)*(abs(fft(data)[2:301]))^2

#1
par(mar=c(4.1,4,1.5,1.5),mgp=c(3,0.5,0)) 
plot((1:3001)/1000, data ,type='l',xaxt="n",yaxt="n",xlab='',ylab='',ylim=c(-150,520))
axis(side = 1, tck = -0.02) ;axis(side = 2, tck = -0.02)
title(xlab="Time (Second)",ylab=bquote('Voltage: '~mu~'V'), line=2, cex.lab=1.2)
grid()
lines((1:3001)/1000, rep(expectile(data,0.1),3001),col='blue')
lines((1:3001)/1000, rep(expectile(data,0.9),3001),col='red')
legend("topright", c('Expectile: 0.1','Expectile: 0.9'),inset = 0.01,
       lty=c(1, 1, 1),  col=c( "blue", "red"),
       seg.len=2)
lines(c(1.4,1.85),c(-40,-40),lty=2)
lines(c(1.4,1.85),c(90,90),lty=2)
lines(c(1.4,1.4),c(-40,90),lty=2)
lines(c(1.85,1.85),c(-40,90),lty=2)
text(1.625,-100,'Main spike')

lines(c(1.45,1.65),c(160,160),lty=2)
lines(c(1.45,1.65),c(400,400),lty=2)
lines(c(1.45,1.45),c(160,400),lty=2)
lines(c(1.65,1.65),c(160,400),lty=2)
text(1.55,450,'Bursts')


plot((1:150)/3000,(pg/sum(pg))[1:150],ylim=c(0, 0.38) ,type='l',xaxt="n",yaxt="n",xlab='',ylab='')
lines((1:150)/3000, (ep[,86]/sum(ep[,86]))[1:150], col='blue')

axis(side = 1, tck = -0.02) ;axis(side = 2, tck = -0.02)
title(xlab="Freq (cycles per ms)",ylab="Periodograms", line=2, cex.lab=1.2)
grid()
text(10/3000, 0.35 ,'Main spike')
text(75/3000, 0.05 ,'Bursting')
legend("topright", c(expression('Ordinary','Expectile: 0.9')),inset = 0.01,
       lty=c(1, 1),  col=c("black",'blue'),
       seg.len=2)


#3
par(mar=c(4.1,4,1.5,1.5),mgp=c(3,0.5,0)) 
plot((1:150)/3000, (ep[,6]/sum(ep[,6]))[1:150],ylim=c(0, 0.45), type='l', col='black', xlab='', ylab='',xaxt='n',yaxt='n')
axis(side = 1, tck = -0.02) ;axis(side = 2, tck = -0.02)
grid()
title(xlab="Freq (cycles per ms)",ylab="Expectile periodogram", line=2, cex.lab=1.2)
text(10/3000, 0.42 ,'Main spike')
text(75/3000, 0.05 ,'No bursting')

legend("topright", c(expression('Expectile: 0.1')),inset = 0.01,
       lty=c(1, 1, 1),  col=c("black"),
       seg.len=2)

#4  
ep = spec.normalize(ep)
par(mar=c(4.1,4,1.5,1.5),mgp=c(3,0.5,0)) 
image.plot((30:120)/3000, tau ,ep[30:120,] ,xaxt="n",yaxt="n",xlab='',ylab='')
axis(side = 1, tck = -0.02) ;axis(side = 2, tck = -0.02)
title(xlab="Freq (cycles per ms)",ylab="Expectile", line=2, cex.lab=1.2)

#####################sp500 example
tau=(5:95)/100
data = read.csv('sp.csv')
data = data$Close[14540:22103]
data = diff(log(data))
ep1 = EP_est(data, 0.025, tau, 8)

pg = (1/length(data))*(abs(fft(data)[2:190]))^2
pgsmo = smooth.spline(pg,spar=0.8)$y
pgsmo = pgsmo/sum(pgsmo)
ep = spec.normalize(ep1$qp)
epsmo = ep
for(i in 1:91){
  epsmo[,i] = smooth.spline(ep[,i],spar=0.8)$y
}
epsmo = spec.normalize(epsmo)

#1
par(mar=c(4.1,4,1.5,1.5),mgp=c(3,0.5,0)) 
plot((1:7563)/252+1986, data ,type='l',xaxt="n",yaxt="n",xlab='',ylab='',ylim=c(-0.1,0.113), xlim = c(1986,2016))
axis(side = 1, tck = -0.02) ;axis(side = 2, tck = -0.02)
title(xlab="Year",'Daily log return', line=2, cex.lab=1.2)
grid()
lines((1:7563)/252+1986, 0.014+0.007*sin((1:7563)/410 -0.5),col='blue')
lines((1:7563)/252+1986, -0.014-0.007*sin((1:7563)/410 -0.5),col='blue')
legend("top", c('Daily log return','10-year cycle'),inset = 0.01,
       lty=c(1, 1),  col=c( "black", "blue"),
       seg.len=2)

#2
par(mar=c(4.1,4,1.5,1.5),mgp=c(3,0.5,0)) 
plot((1:189)/7560, pg/sum(pg),lty = 1, type='l', col='black', xlab='', ylab='',xaxt='n',yaxt='n',ylim=c(0, 0.13))
lines((1:189)/7560,ep[,86],lty=1,col='blue')
grid()
axis(side = 1, tck = -0.02) ;axis(side = 2, tck = -0.02)
title(xlab="Freq (cycles per day)",ylab="Periodograms", line=2, cex.lab=1.2)
legend("topright", c(expression('Ordinary','Expectile: 0.9')),inset = 0.01,
       lty=c(1, 1),  col=c("black","blue"),
       seg.len=2)
text(8/3000, 0.12 ,'10-year cycle')

#3
image.plot((1:189)/7560, tau ,ep ,xaxt="n",yaxt="n",xlab='',ylab='',xlim=c(0,0.025))
axis(side = 1, tck = -0.02) ;axis(side = 2, tck = -0.02)
title(xlab="Freq (cycles per day)",ylab="Expectile", line=2, cex.lab=1.2)


#4
par(mar=c(4.1,4,1.5,1.5),mgp=c(3,0.5,0)) 
plot((1:189)/7560, pgsmo/sum(pgsmo),lty = 1, type='l', col='black', xlab='', ylab='',xaxt='n',yaxt='n',ylim=c(0, 0.08))
lines((1:189)/7560,epsmo[,86],lty=1,col='blue')
grid()
axis(side = 1, tck = -0.02) ;axis(side = 2, tck = -0.02)
title(xlab="Freq (cycles per day)",ylab="Periodograms", line=2, cex.lab=1.2)
legend("topright", c(expression('Ordinary','Expectile: 0.9')),inset = 0.01,
       lty=c(1, 1),  col=c("black","blue"),
       seg.len=2)

#5
image.plot((1:189)/7560, tau ,epsmo ,xaxt="n",yaxt="n",xlab='',ylab='',xlim=c(0,0.025))
axis(side = 1, tck = -0.02) ;axis(side = 2, tck = -0.02)
title(xlab="Freq (cycles per day)",ylab="Expectile", line=2, cex.lab=1.2)



#6
load('garch_sdf.rdata')
par(mar=c(4.1,4,1.5,1.5),mgp=c(3,0.5,0)) 
plot(1:199/398, pgram[2:200], type='l', col='black',ylim=c(0.004,0.01), xlab='', ylab='',xaxt='n',yaxt='n')
lines(1:199/398, garch_sdf[,86],col='blue')
axis(side = 1, tck = -0.02) ;axis(side = 2, tck = -0.02)
grid()
title(xlab="Freq",ylab="Periodograms", line=2, cex.lab=1.2)
legend(x = 'topright', 
       legend = c('Ordinary','Expectile: 0.9'),
       lty = 1,
       col = c( 'black','blue'),inset = 0.01)

#7
image.plot(1:199/398, y = tau, z=garch_sdf,xaxt='n',yaxt='n',xlab='',ylab='',xlim=c(0,0.5))
axis(side = 1, tck = -0.02) ;axis(side = 2, tck = -0.02)
title(xlab="Freq",ylab="Expectile", line=2, cex.lab=1.2)


##mixture example
s = 5000
n = 200
tau = (5:95)/100
mix_ep = matrix(0,100-1,91)
mix_pg = rep(0,99)
mix_qp = rep(0,99)
set.seed(1)
for(i in 1:s){   ##GARCH
  data1 = mymix(200)
  n <- length(data1)
  freq <- c(0:(n-1))/n
  freq <- freq[freq>0 & freq <0.5]
  q1 <- expectile_peri_ls(data1, freq, tau, n.cores=8)
  q1 = spec.normalize(q1)
  mix_ep = mix_ep + q1/s
  
  mix_pg = mix_pg + smooth.spline((1/n)*(abs(fft(data1)[1:99]))^2,cv=T)$y/s
  mix_qp = mix_qp + smooth.spline(lap.spec.new2(data1,freq,0.5,intercept=T,type=1,n.cores=8),cv=T)$y/s
  print(i)}

for(i in 1:91){
  mix_ep[,i] = smooth.spline(mix_ep[,i], cv=T)$y
}


par(mar=c(4.1,4.5,1.5,1.5),mgp=c(3,0.5,0)) 
image.plot(x = (1:99)/198,y = (5:95)/100, z=spec.normalize(mix_ep),xaxt='n',yaxt='n',xlab='',ylab='')
axis(side = 1, tck = -0.02) ;axis(side = 2, tck = -0.02)
title(xlab="Freq",ylab="Expectile", line=2, cex.lab=1.2)

par(mar=c(4.1,4.5,1.5,1.5),mgp=c(3,0.5,0)) 
plot((1:99)/198, mix_pg/sum(mix_pg),type='l',xlab='',ylab='',col='blue',xaxt='n',yaxt='n',ylim=c(0,0.041))
lines((1:99)/198, mix_ep[,6]/sum(mix_ep[,6]), type='l',col='red')
lines((1:99)/198, mix_ep[,86]/sum(mix_ep[,86]), type='l',col='black')
title(xlab="Freq",ylab="Periodograms", line=2, cex.lab=1.2)
axis(side = 1, tck = -0.02) ;axis(side = 2, tck = -0.02)
grid()
legend(x = 'topright', 
       legend = c('Ordinary','Expectile: 0.1','Expectile: 0.9'),
       lty = 1,
       col = c( 'blue','red','black') ,inset = 0.01)



################################################ simulation 1#################################
#######################omega = .25################################
b0 = 1
b1 = 0.9
b2 = 1
omega0 = 2*pi*0.09
omega1 = 2*pi*0.12
omegac = 2*pi*0.25
r = 0.6   
n = 200
tau = 0.8
ncores = 2

ep_smo25 = 0
qp_smo25 = 0
pg_smo25 = 0
ep_res25 = 0
qp_res25 = 0
pg_res25 = 0
k = 5000
set.seed(1)
for (i in 1:k){
  xt =  as.numeric(arima.sim(list(order=c(2,0,0),ar=c(2*r*cos(omegac),-r^2),sd=1),n=n))
  at = b0 + b1*cos(omega0*(1:n)) + b2*sin(omega1*(1:n))
  yt = xt*at
  freq <- c(0:(n-1))/n
  freq <- freq[freq>0 & freq <0.5]
  
  ep25 = expectile_peri_ls(yt, freq, tau, n.cores=1)
  ep_res25 = ep_res25 + ep25/k
  ep_no25 = expectile_peri_ls(xt, freq, tau, n.cores=1)
  ep_smo25 = ep_smo25 + smooth.spline(ep_no25, cv=T)$y/k
  pg25 = spec.pgram(yt, plot=F)$spec
  pg_res25 = pg_res25 + pg25/k
  pg_no25 = spec.pgram(xt, plot=F)$spec
  pg_smo25 = pg_smo25 + smooth.spline(pg_no25, cv=T)$y/k
  
  
  qp25 = lap.spec.new2(yt,freq,0.5,intercept=T,type=1,n.cores=1)
  qp_res25 = qp_res25 + qp25/k
  qp_no25 = lap.spec.new2(xt,freq,0.5,intercept=T,type=1,n.cores=1)
  qp_smo25 = qp_smo25 + smooth.spline(qp_no25, cv=T)$y/k
  print(i)
}


#######################omega = .30###############################  
omegac = 2*pi*0.3
tau = 0.8
ep_smo30 = 0
qp_smo30 = 0
pg_smo30 = 0
ep_res30 = 0
qp_res30 = 0
pg_res30 = 0
k = 5000
set.seed(1)
for (i in 1:k){
  xt =  as.numeric(arima.sim(list(order=c(2,0,0),ar=c(2*r*cos(omegac),-r^2),sd=1),n=n))
  at = b0 + b1*cos(omega0*(1:n)) + b2*sin(omega1*(1:n))
  yt = xt*at
  freq <- c(0:(n-1))/n
  freq <- freq[freq>0 & freq <0.5]
  
  ep30 = expectile_peri_ls(yt, freq, tau, n.cores=1)
  ep_res30 = ep_res30 + ep30/k
  ep_no30 = expectile_peri_ls(xt, freq, tau, n.cores=1)
  ep_smo30 = ep_smo30 + smooth.spline(ep_no30, cv=T)$y/k
  
  pg30 = spec.pgram(yt, plot=F)$spec
  pg_res30 = pg_res30 + pg30/k
  pg_no30 = spec.pgram(xt, plot=F)$spec
  pg_smo30 = pg_smo30 + smooth.spline(pg_no30, cv=T)$y/k
  
  
  qp30 = lap.spec.new2(yt,freq,0.5,intercept=T,type=1,n.cores=1)
  qp_res30 = qp_res30 + qp30/k
  qp_no30 = lap.spec.new2(xt,freq,0.5,intercept=T,type=1,n.cores=1)
  qp_smo30 = qp_smo30 + smooth.spline(qp_no30, cv=T)$y/k
  print(i)
}

par(mar=c(4.1,4.5,1.5,1.5),mgp=c(3,0.5,0)) 
plot((1:99)/198, pg_smo25[2:100]/sum(pg_smo25[2:100]),type='l',xlab='',ylab='',col='blue',xaxt='n',yaxt='n')
lines((1:99)/198, qp_smo25/sum(qp_smo25), type='l',col='red', xlab='Freq',ylab='')
lines((1:99)/198, ep_smo25/sum(ep_smo25), type='l',col='black', xlab='Freq',ylab='')
title(xlab="Freq",ylab="Periodograms", line=2, cex.lab=1.2)
axis(side = 1, tck = -0.02) ;axis(side = 2, tck = -0.02)
grid()
legend(x = 'topleft', 
       legend = c('Ordinary','Laplace','Expectile: 0.8'),
       lty = 1,
       col = c( 'blue','red','black') ,inset = 0.01)

par(mar=c(4.1,4.5,1.5,1.5),mgp=c(3,0.5,0)) 
plot((1:99)/198, pg_res25[2:100]/sum(pg_res25[2:100]),type='l',xlab='',ylab='',col='blue',ylim=c(0,0.045),xaxt='n',yaxt='n')
lines((1:99)/198, qp_res25/sum(qp_res25), type='l',col='red', xlab='Freq',ylab='')
lines((1:99)/198, ep_res25/sum(ep_res25), type='l',col='black', xlab='Freq',ylab='')
title(xlab="Freq",ylab="Periodograms", line=2, cex.lab=1.2)
axis(side = 1, tck = -0.02) ;axis(side = 2, tck = -0.02)
grid()
legend(x = 'topright', 
       legend = c('Ordinary','Laplace','Expectile: 0.8'),
       lty = 1,
       col = c( 'blue','red','black') ,inset = 0.01)

plot((1:99)/198, pg_smo30[2:100]/sum(pg_smo30[2:100]),type='l',xlab='',ylab='',col='blue',xaxt='n',yaxt='n')
lines((1:99)/198, qp_smo30/sum(qp_smo30), type='l',col='red', xlab='Freq',ylab='')
lines((1:99)/198, ep_smo30/sum(ep_smo30), type='l',col='black', xlab='Freq',ylab='')
title(xlab="Freq",ylab="Periodograms", line=2, cex.lab=1.2)
axis(side = 1, tck = -0.02) ;axis(side = 2, tck = -0.02)
panel.first = grid(6, 6, col = 'grey',lty = 3)
legend(x = 'topleft', 
       legend = c('Ordinary','Laplace','Expectile: 0.8'),
       lty = 1,
       col = c( 'blue','red','black'),inset = 0.01 )

par(mar=c(4.1,4.5,1.5,1.5),mgp=c(3,0.5,0)) 
plot((1:99)/198, pg_res30[2:100]/sum(pg_res30[2:100]),type='l',xlab='',ylab='',col='blue',ylim=c(0,0.045),xaxt='n',yaxt='n')
lines((1:99)/198, qp_res30/sum(qp_res30), type='l',col='red', xlab='Freq',ylab='')
lines((1:99)/198, ep_res30/sum(ep_res30), type='l',col='black', xlab='Freq',ylab='')
title(xlab="Freq",ylab="Periodograms", line=2, cex.lab=1.2)
axis(side = 1, tck = -0.02) ;axis(side = 2, tck = -0.02)
grid()
legend(x = 'topright', 
       legend = c('Ordinary','Laplace','Expectile: 0.8'),
       lty = 1,
       col = c( 'blue','red','black') ,inset = 0.01)

#########################image
omegac = 2*pi*0.25
k = 1000
set.seed(1)
ep_ima25 = 0
for (i in 1:k){
  xt =  as.numeric(arima.sim(list(order=c(2,0,0),ar=c(2*r*cos(omegac),-r^2),sd=1),n=n))
  at = b0 + b1*cos(omega0*(1:n)) + b2*sin(omega1*(1:n))
  yt = xt*at
  freq <- c(0:(n-1))/n
  freq <- freq[freq>0 & freq <0.5]
  
  ep_no25 = expectile_peri_ls(xt, freq, (5:95)/100, n.cores=16)
  ep_no25 = spec.normalize(ep_no25)
  ep_ima25 =  ep_ima25 + ep_no25
  print(i)
}
for(i in 1:91){
  ep_ima25[,i] = smooth.spline(ep_ima25[,i], cv=T)$y
}
par(mar=c(4.1,4.5,1.5,1.5),mgp=c(3,0.5,0)) 
image.plot(x = (1:99)/198,y = (5:95)/100, z=spec.normalize(ep_ima25),xaxt='n',yaxt='n',xlab='',ylab='')
axis(side = 1, tck = -0.02) ;axis(side = 2, tck = -0.02)
title(xlab="Freq",ylab="Expectile", line=2, cex.lab=1.2)


omegac = 2*pi*0.3
k = 500
set.seed(1)
ep_ima30 = 0
for (i in 1:k){
  xt =  as.numeric(arima.sim(list(order=c(2,0,0),ar=c(2*r*cos(omegac),-r^2),sd=1),n=n))
  at = b0 + b1*cos(omega0*(1:n)) + b2*sin(omega1*(1:n))
  yt = xt*at
  freq <- c(0:(n-1))/n
  freq <- freq[freq>0 & freq <0.5]
  
  ep_no30 = expectile_peri_ls(xt, freq, (5:95)/100, n.cores=16)
  ep_no30 = spec.normalize(ep_no30)
  ep_ima30 =  ep_ima30 + ep_no30
  print(i)
}
for(i in 1:91){
  ep_ima30[,i] = smooth.spline(ep_ima30[,i], cv=T)$y
}
par(mar=c(4.1,4.5,1.5,1.5),mgp=c(3,0.5,0)) 
image.plot(x = (1:99)/198,y = (5:95)/100, z=spec.normalize(ep_ima30),xaxt='n',yaxt='n',xlab='',ylab='')
axis(side = 1, tck = -0.02) ;axis(side = 2, tck = -0.02)
title(xlab="Freq",ylab="Expectile", line=2, cex.lab=1.2)


################################Fisher's test############################################
pfish = function(x){
  q = length(x)
  g = q * max(x)/sum(x)
  js = 0:q
  s = 1-js*g/q
  for(i in 1:(q+1)){
    if (s[i] <0) s[i] = 0
  }
  pp = 1 - sum((-1)^js * choose(q,js) * s^(q-1))
  return(pp)
}
b0 = 1
b1 = 0.9
b2 = 0
omega0 = 2*pi*0.1
omega1 = 2*pi*0.12
omegac = 2*pi*0.3
r = 0.6   
n = 200
tau = c(0.85,0.9,0.95)
ncores = 2
g = rep(0,24)
k = 1000
set.seed(2)
for (i in 1:k){
  xt =  as.numeric(arima.sim(list(order=c(2,0,0),ar=c(2*r*cos(omegac),-r^2),sd=1),n=n))
  at = b0 + b1*cos(omega0*(1:n)) + b2*sin(omega1*(1:n))
  yt = xt*at
  freq <- c(0:(n-1))/n
  freq <- freq[freq>0 & freq <0.5]
  
  ep = expectile_peri(yt, freq, tau, n.cores=ncores)[1:50, ]
  ep7 = ep[,1]
  ep8 = ep[,2]
  ep9 = ep[,3]
  
  g[1] = g[1] + (pfish(ep7) < 0.01)
  g[2] = g[2] + (pfish(ep7) < 0.05)
  g[3] = g[3] + (pfish(ep7) < 0.1)
  g[4] = g[4] + (pfish(ep8) < 0.01)
  g[5] = g[5] + (pfish(ep8) < 0.05)
  g[6] = g[6] + (pfish(ep8) < 0.1)
  g[7] = g[7] + (pfish(ep9) < 0.01)
  g[8] = g[8] + (pfish(ep9) < 0.05)
  g[9] = g[9] + (pfish(ep9) < 0.1)
  
  qp = lap.spec.new2(yt,freq,tau,intercept=T,type=1,n.cores=ncores)[1:50, ]
  qp7 = qp[,1]
  qp8 = qp[,2]
  qp9 = qp[,3]
  
  g[10] = g[10] + (pfish(qp7) < 0.01)
  g[11] = g[11] + (pfish(qp7) < 0.05)
  g[12] = g[12] + (pfish(qp7) < 0.1)
  g[13] = g[13] + (pfish(qp8) < 0.01)
  g[14] = g[14] + (pfish(qp8) < 0.05)
  g[15] = g[15] + (pfish(qp8) < 0.1)
  g[16] = g[16] + (pfish(qp9) < 0.01)
  g[17] = g[17] + (pfish(qp9) < 0.05)
  g[18] = g[18] + (pfish(qp9) < 0.1)
  
  qp5 = lap.spec.new2(yt,freq,0.5,intercept=T,type=1,n.cores=ncores)[1:50]
  g[19] = g[19] + (pfish(qp5) < 0.01)
  g[20] = g[20] + (pfish(qp5) < 0.05)
  g[21] = g[21] + (pfish(qp5) < 0.1)
  
  pg = spec.pgram(yt, plot=F)$spec[1:50]
  g[22] = g[22] + (pfish(pg) < 0.01)
  g[23] = g[23] + (pfish(pg) < 0.05)
  g[24] = g[24] + (pfish(pg) < 0.1)
  print(i)
}

##########################smoothing
ncores = 12
set.seed(1)
GT=list() 

GT[[1]] = rep(0,99)
GT[[2]] = rep(0,199)
GT[[3]] = rep(0,399)
GT[[4]] = rep(0,799)

set.seed(1)
omega0 = 2*pi*0.09
omega1 = 2*pi*0.12
omegac = 2*pi*0.25
r = 0.6   
tau = 0.9
for(i in 1:4000){  ##AR2
  data1 = as.numeric(arima.sim(list(order=c(2,0,0),ar=c(2*r*cos(omegac),-r^2),sd=1),n=200))
  n <- length(data1)
  freq <- c(0:(n-1))/n
  freq <- freq[freq>0 & freq<0.5]
  q1 <- expectile_peri_ls(data1,freq,tau,n.cores=ncores)
  q1 = spec.normalize(q1)
  
  data2 = as.numeric(arima.sim(list(order=c(2,0,0),ar=c(2*r*cos(omegac),-r^2),sd=1),n=400))
  n <- length(data2)
  freq <- c(0:(n-1))/n
  freq <- freq[freq>0 & freq<0.5]
  q2 <- expectile_peri_ls(data2,freq,tau,n.cores=ncores)
  q2 = spec.normalize(q2)
  
  data3 = as.numeric(arima.sim(list(order=c(2,0,0),ar=c(2*r*cos(omegac),-r^2),sd=1),n=800))
  n <- length(data3)
  freq <- c(0:(n-1))/n
  freq <- freq[freq>0 & freq<0.5]
  q3 <- expectile_peri_ls(data3,freq,tau,n.cores=ncores)
  q3 = spec.normalize(q3)
  
  data4 = as.numeric(arima.sim(list(order=c(2,0,0),ar=c(2*r*cos(omegac),-r^2),sd=1),n=1600))
  n <- length(data4)
  freq <- c(0:(n-1))/n
  freq <- freq[freq>0 & freq<0.5]
  q4 <- expectile_peri_ls(data4,freq,tau,n.cores=ncores)
  q4 = spec.normalize(q4)
  
  GT[[1]] = GT[[1]] + q1/4000
  GT[[2]] = GT[[2]] + q2/4000
  GT[[3]] = GT[[3]] + q3/4000
  GT[[4]] = GT[[4]] + q4/4000
  if(i %% 20 == 0){
    print(c('ar',i))}
}

GT[[5]] = rep(0,99)
GT[[6]] = rep(0,199)
GT[[7]] = rep(0,399)
GT[[8]] = rep(0,799)
for(i in 1:4000){  ##mix
  data1 = mymix(200)
  n <- length(data1)
  freq <- c(0:(n-1))/n
  freq <- freq[freq>0 & freq<0.5]
  q1 <- expectile_peri_ls(data1,freq,tau,n.cores=ncores)
  q1 = spec.normalize(q1)
  
  data2 = mymix(400)
  n <- length(data2)
  freq <- c(0:(n-1))/n
  freq <- freq[freq>0 & freq<0.5]
  q2 <- expectile_peri_ls(data2,freq,tau,n.cores=ncores)
  q2 = spec.normalize(q2)
  
  data3 = mymix(800)
  n <- length(data3)
  freq <- c(0:(n-1))/n
  freq <- freq[freq>0 & freq<0.5]
  q3 <- expectile_peri_ls(data3,freq,tau,n.cores=ncores)
  q3 = spec.normalize(q3)
  
  data4 = mymix(1600)
  n <- length(data4)
  freq <- c(0:(n-1))/n
  freq <- freq[freq>0 & freq<0.5]
  q4 <- expectile_peri_ls(data4,freq,tau,n.cores=ncores)
  q4 = spec.normalize(q4)
  
  GT[[5]] = GT[[5]] + q1/4000
  GT[[6]] = GT[[6]] + q2/4000
  GT[[7]] = GT[[7]] + q3/4000
  GT[[8]] = GT[[8]] + q4/4000
  if(i %% 20 == 0){
    print(c('mix',i))}
}

for(i in 1:8){
  GT[[i]] = smooth.spline(GT[[i]],cv=T)$y
}

set.seed(1) #####AR2
nsim = 5000
mse1 = 0; kl1 = 0
mse2 = 0; kl2 = 0
mse3 = 0; kl3 = 0
mse4 = 0; kl4 = 0
for(j in 1:nsim){
  data1 = as.numeric(arima.sim(list(order=c(2,0,0),ar=c(2*r*cos(omegac),-r^2),sd=1),n=200))
  n <- length(data1)
  freq <- c(0:(n - 1)) / n 
  freq <- freq[freq>0 & freq <0.5]
  tmp = expectile_peri_ls(data1, freq, tau,8)
  res1 = smooth.spline(tmp, cv=T)$y
  res1 = res1/sum(res1)
  mse1 = mse1 + (sum((res1 - GT[[1]])^2)/length(res1))/nsim
  kl1 = kl1 + KL.plugin(res1, GT[[1]])/nsim
  
  
  data2 = as.numeric(arima.sim(list(order=c(2,0,0),ar=c(2*r*cos(omegac),-r^2),sd=1),n=400))
  n <- length(data2)
  freq <- c(0:(n - 1)) / n 
  freq <- freq[freq>0 & freq <0.5]
  tmp = expectile_peri_ls(data2, freq, tau,8)
  res2 = smooth.spline(tmp, cv=T)$y
  res2 = res2/sum(res2)
  mse2 = mse2 + (sum((res2 - GT[[2]])^2)/length(res2))/nsim
  kl2 = kl2 + KL.plugin(res2, GT[[2]])/nsim
  
  data3 = as.numeric(arima.sim(list(order=c(2,0,0),ar=c(2*r*cos(omegac),-r^2),sd=1),n=800))
  n <- length(data3)
  freq <- c(0:(n - 1)) / n 
  freq <- freq[freq>0 & freq <0.5]
  tmp = expectile_peri_ls(data3, freq, tau,8)  
  res3 = smooth.spline(tmp, cv=T)$y
  res3 = res3/sum(res3)
  mse3 = mse3 + (sum((res3 - GT[[3]])^2)/length(res3))/nsim
  kl3 = kl3 + KL.plugin(res3, GT[[3]])/nsim
  
  data4 = as.numeric(arima.sim(list(order=c(2,0,0),ar=c(2*r*cos(omegac),-r^2),sd=1),n=1600))
  n <- length(data4)
  freq <- c(0:(n - 1)) / n 
  freq <- freq[freq>0 & freq <0.5]
  tmp = expectile_peri_ls(data4, freq, tau,8)
  res4 = smooth.spline(tmp, cv=T)$y
  res4 = res4/sum(res4)
  mse4 = mse4 + (sum((res4 - GT[[4]])^2)/length(res4))/nsim
  kl4 = kl4 + KL.plugin(res4, GT[[4]])/nsim
  print(c('ar',j))
}


set.seed(1) #####mix
nsim = 5000
mse11 = 0; kl11 = 0
mse22 = 0; kl22 = 0
mse33 = 0; kl33 = 0
mse44 = 0; kl44 = 0
for(j in 1:nsim){
  data1 = mymix(200)
  n <- length(data1)
  freq <- c(0:(n - 1)) / n 
  freq <- freq[freq>0 & freq <0.5]
  tmp = expectile_peri_ls(data1, freq, tau,8)
  res11 = smooth.spline(tmp, cv=T)$y
  res11 = res11/sum(res11)
  mse11 = mse11 + (sum((res11 - GT[[5]])^2)/length(res11))/nsim
  kl11 = kl11 + KL.plugin(res11, GT[[5]])/nsim
  
  
  data2 = mymix(400)
  n <- length(data2)
  freq <- c(0:(n - 1)) / n 
  freq <- freq[freq>0 & freq <0.5]
  tmp = expectile_peri_ls(data2, freq, tau,8)
  res22 = smooth.spline(tmp, cv=T)$y
  res22 = res22/sum(res22)
  mse22 = mse22 + (sum((res22 - GT[[6]])^2)/length(res22))/nsim
  kl22 = kl22 + KL.plugin(res22, GT[[6]])/nsim
  
  data3 = mymix(800)
  n <- length(data3)
  freq <- c(0:(n - 1)) / n 
  freq <- freq[freq>0 & freq <0.5]
  tmp = expectile_peri_ls(data3, freq, tau,8)
  res33 = smooth.spline(tmp, cv=T)$y
  res33 = res33/sum(res33)
  mse33 = mse33 + (sum((res33 - GT[[7]])^2)/length(res33))/nsim
  kl33 = kl33 + KL.plugin(res33, GT[[7]])/nsim
  
  data4 = mymix(1600)
  n <- length(data4)
  freq <- c(0:(n - 1)) / n 
  freq <- freq[freq>0 & freq <0.5]
  tmp = expectile_peri_ls(data4, freq, tau,8)
  res44 = smooth.spline(tmp, cv=T)$y
  res44 = res44/sum(res44)
  mse44 = mse44 + (sum((res44 - GT[[8]])^2)/length(res44))/nsim
  kl44 = kl44 + KL.plugin(res44, GT[[8]])/nsim
  print(c('mix',j))
}

par(mar=c(4.1,4.5,1.5,1.5),mgp=c(3,0.5,0)) 
plot(c(200,400,800,1600),c(mse1,mse2,mse3,mse4)/mse1,xlab='',ylab='',xaxt='n',yaxt='n')
lines(c(200,400,800,1600),c(mse1,mse2,mse3,mse4)/mse1)
points(c(200,400,800,1600),c(kl1,kl2,kl3,kl4)/kl1,col='blue',pch=5)
lines(c(200,400,800,1600),c(kl1,kl2,kl3,kl4)/kl1,col='blue')

axis(side = 1, tck = -0.02) ;axis(side = 2, tck = -0.02)
title(xlab="Sample size",ylab="Error", line=2, cex.lab=1.2)
grid()
legend("topright", c(expression('MSE ('%*%'5.57e-06)'),expression('KL    ('%*%'2.50e-02)')),inset = 0.01,lty=c(1,1),pch=c(1,5),
       col=c('black','blue'))


par(mar=c(4.1,4.5,1.5,1.5),mgp=c(3,0.5,0)) 
plot(c(200,400,800,1600),c(mse11,mse22,mse33,mse44)/mse11,xlab='',ylab='',xaxt='n',yaxt='n')
lines(c(200,400,800,1600),c(mse11,mse22,mse33,mse44)/mse11)
points(c(200,400,800,1600),c(kl11,kl22,kl33,kl44)/kl11,col='blue',pch=5)
lines(c(200,400,800,1600),c(kl11,kl22,kl33,kl44)/kl11,col='blue')
axis(side = 1, tck = -0.02) ;axis(side = 2, tck = -0.02)
title(xlab="Sample size",ylab="Error", line=2, cex.lab=1.2)
grid()
legend("topright", c(expression('MSE ('%*%'1.48e-05)'),expression('KL    ('%*%'5.32e-02)')),inset = 0.01,lty=c(1,1),pch=c(1,5),
       col=c('black','blue'))



