# Assignment 1
#
# Compute numerical solution of boundary value problem:
# -au"+bu'= x^2,  x in (0;1); u(0) = 1, u(1) = 0
# N - discratization parameter
# difference scheme 'central' or 'bakwards'
# Compare with exact solution.

install.packages("GetoptLong")
library(GetoptLong)

Equation <- function(a,b) {
  # Computes the exact solution of the problem and parameters 
  # are needed for numerical solution.
  #
  # Args:
  #   a: coefficient at -u"
  #   b: coefficient at u'
  #  	
  # Returns:
  #   list with parameters of b.v. problem:
  #     boundary.left: first boundary condition
  #     boundary.right: second boundary condition
  #  	right.hand: righthand side function
  #     exact.solution: function calculating exact solution
  #     diff.scheme.coeff: vector of difference scheme coefficients

  boundary.left <- 1; boundary.right <- 0
	
  RightHand <- function(xi) { xi^2 } 
    
  ExactSolution <- function(xi) {
  	# Computes exact solution of b.v. problem.
  	# 
  	# Args:
  	#   xi: mesh vector
  	#
  	# Returns:
  	#   solution vector u(xi).
		if (b == 0) {
			u <- -(xi^4)/(12*a) - (1-1/(12*a))*xi + 1 			
		} else {
			c1 <- (-1/(3*b) - a/(b^2) - (2*a^2)/(b^3) - 1)/(1-exp(-b/a))
			c2 <- 1 - (-1/(3*b) - a/(b^2) - (2*a^2)/(b^3) - 1)*exp(-b/a)/(1-exp(-b/a))
			c3 <- 1/(3*b)
			c4 <- a/b^2
			c5 <- 2*a^2/b^3
			
			u <- c1*exp((xi-1)*b/a) + c2 + c3*xi^3 + c4*xi^2 + c5*xi
		}
	}

	DiffSchemeCoeff <- function(h, scheme) {
	  # Computes the difference scheme's coefficients
	  #
	  # Args:
	  #   h: step of discretisation
	  #   scheme: difference scheme "central" or "backward"
	  #
	  # Returns:
	  #   vector(k1,k2,k3): k1,k2,k3 - coefficients in diff. scheme 
	  #                     k1*U_{i-1} +k2*U_i + k3*U_{i+1}
		if (scheme == "central") {
			list(k1 = -a/(h^2) - b/(2*h),
				 k2 = 2*a/(h^2),
				 k3 = -a/(h^2) + b/(2*h))
		}	
		else if (scheme == "backward") {
			list(k1 = -a/(h^2) - b/h,
				 k2 = 2*a/(h^2)+b/h,
				 k3 = -a/(h^2))
		}		
	 }
	
		
	list(boundary.left = boundary.left, 
		 boundary.right = boundary.right, 
		 right.hand = RightHand,
		 exact.solution = ExactSolution,
		 diff.scheme.coeff = DiffSchemeCoeff)
	
}

NumSolution <- function(a, b, N, scheme = "central") {
  # Computes numerical solution of b.v. problem and compares it with exact solution
  #
  # Args:
  #   a: coefficient at -u"
  #   b: coefficient at u'
  #   N: parameter of discretisation
  #   scheme: difference scheme "central" or "backward"
  #
  # Returns:
  #   list of:
  #     data: data.frame of 3 vectors: exact result, numerical result and mesh.
  #     error: max error of numerical solution, max|Ui - u(xi)|
  #     a: coefficient at -u"
  #     b: coefficient at u'
  #     N: parameter of discretisation
  #     scheme: difference scheme "central" or "backward"   
  equation <- Equation(a, b)  # calculating all needed parameters of the problem
	
  h <- 1/N  # step of discretisation

  xi <- seq(0, 1, by = h)	 # mesh-vector of xi, i = 1...N+1 
  vF <- c(equation$boundary.left, 
          equation$right.hand(xi[2:N]), 
          equation$boundary.right)  # right-hand side vector f(xi)
  	
  # Creating tridiagonal matrix mA
  coefficients <- equation$diff.scheme.coeff(h, scheme)  # coefficients of difference scheme
		
  d1 <- c(   rep(coefficients$k1, N-1), 0)		# subdiagonal vector with length = N-1
  d2 <- c(1, rep(coefficients$k2, N-1), 1)	    # diagonal vector with length = N
  d3 <- c(0, rep(coefficients$k3, N-1)   )		# superdiagonal vector with length  = N-1

  diag.num <- -outer(seq(d2), seq(d2),"-")
  mA <- diag(d2)
  mA[diag.num == 1] <- d3
  mA[diag.num == -1] <- d1	


  # Solving system of linear equations mA %*% vU =vF
  vU <-  solve(mA, vF)

  result <- data.frame(exact = equation$exact.solution(xi), numeric = vU, mesh = xi)
	
  error = max(abs(result$exact - result$numeric))
  out <- list(data=result, error=error, a=a, b=b, N=N, scheme=scheme)
	
  class(out) <- "NumSolution"
	
  return(out)
}

print.NumSolution <- function(object) {
	# Prints output of NumSolution
	cat("Max. error: ", object$error, "\n\n")
	print(object$data)
}

plot.NumSolution <- function(object) {
	# Creates plot with numerical and exact solutions vs. xi
	data <- object$data
	title <- qq("N=@{object$N}, @{object$scheme} scheme, \n max. error=@{round(object$error,3)}")
	plot(data$mesh, data$numeric, type="b", col="blue", xlab ="x", ylab="u(x)", main = title, cex = 0.6, cex.main = 1, font.main = 2)
	lines(data$mesh, data$exact, col = "red")
	legend("bottomleft", legend=c("numerical", "exact"), col = c("blue", "red"), lwd=1, pch = c(1, NA))
}
# _____________________________________________________________________________

NumSolution(1e-3,1,32,"central")$error
NumSolution(1e-3,1,128,"central")$error

NumSolution(1e-3,1,32,"backward")$error
NumSolution(1e-3,1,128,"backward")$error

#png("case3.png")
par(mfrow = c(2,2), oma  = c(0,0,2,0))
plot(NumSolution(1e-3,1,32,"central"))
plot(NumSolution(1e-3,1,128,"central"))
plot(NumSolution(1e-3,1,32,"backward"))
plot(NumSolution(1e-3,1,128,"backward"))
title("Solution for a = 0.001, b = 1", outer = TRUE)
#dev.off()
#_____________________________________________






