#install.packages("GetoptLong")
library(GetoptLong)

# Assignment 1

# input parameters
# a,b - coefficients in: -au"+bu'=x^2
# N - discratization parameter
# difference scheme ('central' or 'bakwards')

Equation <- function(a,b) {
		
	#boundary conditions
	b1 <- 1; b2 <- 0
	
	#righthand side function
	RightHand <- function(xi) { xi^2 }
	
	ExactSolution <- function(xi) {
		if (b == 0) {
			u <- -(xi^4)/(12*a) - (1-1/(12*a))*xi + 1 			
		} 
		else {
			c1 <- (-1/(3*b) - a/(b^2) - (2*a^2)/(b^3) - 1)/(1-exp(-b/a))
			c2 <- 1 - (-1/(3*b) - a/(b^2) - (2*a^2)/(b^3) - 1)*exp(-b/a)/(1-exp(-b/a))
			c3 <- 1/(3*b)
			c4 <- a/b^2
			c5 <- 2*a^2/b^3
			
			u <- c1*exp((xi-1)*b/a) + c2 + c3*xi^3 + c4*xi^2 + c5*xi
		}
	}
	
	# Difference scheme's coefficients
	# k1*U_{i-1} +k2*U_i + k3*U_{i+1} = f(x_i)
	DiffSchemeCoeff <- function(h, scheme) {
		if (scheme =="central") {
			list(k1 = -a/(h^2) - b/(2*h),
				 k2 = 2*a/(h^2),
				 k3 = -a/(h^2) + b/(2*h))
		}	
		else if (scheme =="backward"){
			list(k1 = -a/(h^2) - b/h,
				 k2 = 2*a/(h^2)+b/h,
				 k3 = -a/(h^2))
		}		
	}
	
		
	out <- list(boundary1 = b1, 
				boundary2 = b2, 
				right_hand = RightHand,
				exact_solution = ExactSolution,
				diff_scheme_coeff = DiffSchemeCoeff)
	
	return(out)
}


PDEsolution <- function(a, b, N, scheme = "central") { 
	equation <- Equation(a, b) 
	#1 numerical solution________________________________________
	
	#computed parameters
	#step
	h <- 1/N
	
	# creating right-hand side vector vF
	xi <- seq(0, 1, by = h)			# vector with xi, i = 1...N+1 
    vF <- c(equation$boundary1, equation$right_hand(xi[2:N]), equation$boundary2)		# right-hand side vector f(x)
	
	# coefficients of difference scheme
	coefficients <- equation$diff_scheme_coeff(h, scheme)	
		
	# creating tridiagonal matrix mA
		d1 <- c(   rep(coefficients$k1, N-1), 0)		# subdiagonal vector with length = N-1
		d2 <- c(1, rep(coefficients$k2, N-1), 1)	    # diagonal vector with length = N
		d3 <- c(0, rep(coefficients$k3, N-1)   )		# superdiagonal vector with length  = N-1

		diag.num <- -outer(seq(d2), seq(d2),"-")
		mA <- diag(d2)
		mA[diag.num == 1] <- d3
		mA[diag.num == -1] <- d1	


	# solving linear system of equations mA %*% vU =vF
		vU <-  solve(mA, vF)

	result <- data.frame(equation$exact_solution(xi), vU, xi)
	names(result) <- c("exact","numeric","mesh")
	
	error = max(abs(result$exact - result$numeric))
	out <- list(data=result, error=error, a=a, b=b, N=N, scheme=scheme)
	
	class(out) <- "PDEsolution"
	
	return(out)
}

print.PDEsolution <- function(object) {
	cat("Max. error: ", object$error, "\n\n")
	print(object$data)
}

plot.PDEsolution <- function(object) {
	data <- object$data
	title <- qq("a=@{object$a}, b=@{object$b}, N=@{object$N}, \n @{object$scheme}")
	plot(data$mesh, data$exact, type="b", col="blue", xlab ="x", ylab="u(x)", main = title)
	lines(data$mesh, data$numeric, col = "red")
	legend("bottomleft", legend=c("exact", "numerical"), col = c("blue", "red"), lwd=1)
}

# parameters variation__________________
PDEsolution(1,0,16,"central")

par(mfrow=c(2,2))
plot(PDEsolution(1e-3,1,16,"central"))
plot(PDEsolution(1e-3,1,128,"central"))
plot(PDEsolution(1e-3,1,32,"backward"))
plot(PDEsolution(1e-3,1,128,"backward"))


