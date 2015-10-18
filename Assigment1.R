
# Assignment 1

# function description

	PDEsolution <- function(a,b,N,scheme){
		# input parameters
		# a,b - coefficients in: -au"+bu'=x^2
		# N - discratization parameter
		# difference scheme ('central' or 'bakwards')
		
		
		#1 numerical solution________________________________________
		
		#computed parameters
		#step
			h<-1/N
		#boundary conditions
		 	b1<- 1; b2<-0
		# Difference scheme's coefficients
			if (scheme =="central"){
				k1 <- -a/(h^2) - b/(2*h); k2 <- 2*a/(h^2); k3<- -a/(h^2)+b/(2*h)	
			}	
			if (scheme =="backward"){
				k1 <- -a/(h^2) - b/(2*h); k2 <- 2*a/(h^2)+b/h; k3<- -a/(h^2)	
			}	
		# creating tridiagonal matrix mA
			d1 <- c(rep(k1,N-1),0)		# subdiagonal vector with length = N-1
			d2 <- c(1,rep(k2,N-1),1)	# diagonal vector with length = N
			d3 <- c(0,rep(k3,N-1))		# superdiagonal vector with length  = N-1

			diag.num <- -outer(seq(d2),seq(d2),"-")
			mA <- diag(d2)
			mA[diag.num == 1] <- d3
			mA[diag.num == -1] <- d1	

		# creating right-hand side vector vF
			xi <- seq(0, 1, by = h)			# vector with xi, i = 1...N+1 
			vF <- c(b1, xi[2:N]^2 ,b2)		# right-hand side vector f(x)
	
		# solving linear system of equations mA %*% vU =vF
			vU <-  solve(mA, vF)

		#2 exact solution___________________________________________________
		
			if (b == 0) {
				u <- -(xi^4)/(12*a) - (1-1/(12*a))*xi + 1 			
			}
			else{
				#c1 <- (-1/(3*b) - a/(b^2) - (2*a^2)/(b^3) - 1)/(exp(b/a)-1)
				#c2 <- 1 - c1
				#u <- c1*exp(xi*b/a) + c2 + (xi^3)/(3*b) + (a*xi^2)/(b^2) + (2*xi*a^2)/(b^3)
				c1 <- (-1/(3*b) - a/(b^2) - (2*a^2)/(b^3) - 1)/(1-exp(-b/a))
				c2 <- 1 - (-1/(3*b) - a/(b^2) - (2*a^2)/(b^3) - 1)*exp(-b/a)/(1-exp(-b/a))
				c3 <- 1/(3*b)
				c4 <- a/b^2
				c5 <- 2*a^2/b^3
				
				u <- c1*exp((xi-1)*b/a) + c2 + c3*xi^3 + c4*xi^2 + c5*xi
			}
		
		result <- data.frame(u,vU,xi,c1,c2,c3,c4,c5)
		names(result) <- c("exact","numeric","mesh","c1","c2","c3","c4","c5")
		
		return(result)
	}
# parameters variation__________________

a <- 1e-3; b <- 1; N <- 128
#scheme <-"central"
scheme <-"backward"
	
result <- PDEsolution(a,b,N,scheme)
result
err <- max(abs(result$exact - result$numeric))
err

r11_16 <- PDEsolution(1e-1,1,16,"central")
r11_32 <- PDEsolution(1e-1,1,32,"central")

r10_16 <- PDEsolution(1,0,16,scheme)
r10_32 <- PDEsolution(1,0,32,scheme)

r1e3_16<-PDEsolution(1e-3,1,16,"central")
r1e3_32<-PDEsolution(1e-3,1,32,"central")
r1e3_128<-PDEsolution(1e-3,1,128,"central")
	

#plots__________________
xi_16 <- seq(0, 1, by = 1/16)			# vector with xi, i = 1...N+1
xi_32 <- seq(0, 1, by = 1/32)
xi_128 <- seq(0, 1, by = 1/128)
colors<-c("blue","red","black","green")

setwd("/Users/daxelka/Education/Numerical solution of PDE")
getwd()

#png("result.png")
par(mfrow=c(1,2))
plot( xi_16, abs(r11_16$exact - r11_16$numeric), type='l', col = colors[1], xlab="xi", ylab="error", main = "a=1, b=1")
lines( xi_32, abs(r11_32$exact - r11_32$numeric), col=colors[2])
#lines( w$Sig2,w$Zcen, col=colors[3])
#lines( w$Sig3,w$Zcen, col=colors[4])
legend("topleft", legend=c("N=16", "N=32"), col = colors, lwd=1)

plot( xi_16, abs(r10_16$exact - r10_16$numeric), type='l', col = colors[1], xlab="xi", ylab="error", main = "a=1, b=0")
lines( xi_32, abs(r10_32$exact - r10_32$numeric), col=colors[2])
legend("topleft", legend=c("N=16", "N=32"), col = colors, lwd=1)
#dev.off()

par(mfrow=c(2,2))
plot(r1e3_16$mesh, r1e3_16$exact, type="l", col=colors[1])
lines(r1e3_16$mesh, r1e3_16$numeric, col = colors[2])
plot(r1e3_32$mesh, r1e3_32$exact, type="l", col=colors[1])
lines(r1e3_32$mesh, r1e3_32$numeric, col = colors[2])
plot(r1e3_128$mesh, r1e3_128$exact, type="l", col=colors[1])
lines(r1e3_128$mesh, r1e3_128$numeric, col = colors[2])

par(mfrow=c(1,2))
plot(r11_16$mesh, r11_16$exact, type="l", col=colors[1])
lines(r11_16$mesh, r11_16$numeric, col = colors[2])
plot(r11_32$mesh, r11_32$exact, type="l", col=colors[1])
lines(r11_32$mesh, r11_32$numeric, col = colors[2])