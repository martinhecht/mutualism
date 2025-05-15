## Changelog:
# MH 0.0.20 2025-05-16: based on c:\Users\martin\Dropbox\151_goodcontrols\01_Posteingang\20250228_Steffen_Formel\bias.pdf

## Documentation
#' @title
#' @description
#' @param A autoregression matrix
#' @param F Number of focal processes
#' @param c1 Number of X-related confounders
#' @param c2 Number of Y-related confounders
#' @return

## Function definition
bias.formula <- function( A, F, c1, c2 ){
		
		# packages
		# require( reshape2 ) # melt()
		
		# total number of confounders
		C <- c1 + c2
		
		# total number of variables
		V <- C + F
		
		if ( F > 0 ){
		
			if( c1 > 0 ) c1t <- (F+1):(F+c1) else c1t <- NULL
			if( c2 > 0 ) c2t <- (F+1 +c1):V  else c2t <- NULL
			
			ct <- list( c1t, c2t )
			
			for( f in 1:F ){
				A[ct[[f]], f] <- beta_covfocal
				A[f, ct[[f]]] <- beta_covfocal
			}
			
			A[c1t, c1t][lower.tri(A[c1t, c1t], diag = FALSE)] <- beta_covfocal
			A[c1t, c1t][upper.tri(A[c1t, c1t], diag = FALSE)] <- beta_covfocal
			
			A[c2t, c2t][lower.tri(A[c2t, c2t], diag = FALSE)] <- beta_covfocal
			A[c2t, c2t][upper.tri(A[c2t, c2t], diag = FALSE)] <- beta_covfocal
			
			A[1:F,1:F][lower.tri(A[1:F,1:F])] <- beta_focal
			A[1:F,1:F][upper.tri(A[1:F,1:F])] <- beta_focal
			
			## Check the max eigenvalue (should be less than 1 for stationarity)
			cat( paste0( "max eigenvalue: ", max(eigen(A)$values), "\n" ) ); flush.console()
			
			## terms for bias formulas
			BY <- A[1:F, 1:F, drop=FALSE]
			if( C>0 ) BYZ <- A[1:F, (F+1):(F+C), drop=FALSE] else BYZ <- NULL
			if( C>0 ) BZY <- A[(F+1):(F+C), 1:F, drop=FALSE] else BZY <- NULL
			if( C>0 ) BZZ <- A[(F+1):(F+C), (F+1):(F+C), drop=FALSE] else BZZ <- NULL
			
			## Var-cov matrix of time-specific residuals
			Sigma <- diag(rep(sigma, V))
			
			## stationary covariance matrix for first time point
			Sigma1 <- irow( solve( diag(dim(A)[1]^2) - t(A) %x% t(A) ) %*% row(Sigma) )
			
			## terms for bias formulas
			if( C>0 ) SigmaZY <- Sigma1[(F+1):(F+C), 1:F] else SigmaZY <- NULL
			SigmaY <- Sigma1[1:F, 1:F]
			
			# identity matrices
			IF <- diag( F )
			
			
			####### results #######
			
			## Eq. 12, Biases of trait matrix
			if( C > 0 ){
				deltaTY <- solve( IF - BY ) %*% BYZ %*% ( BZZ %*% SigmaZY - SigmaZY %*% solve( SigmaY ) %*% BY %*% SigmaY + BZY %*% SigmaY ) %*% solve( SigmaY ) %*% solve( IF - BY ) %*% SigmaY
		
				## Eq. 11, Biases of autoregression matrix (RI-CLPM as analysis model)
				deltaBY <- ( BYZ %*% SigmaZY - ( IF - BY ) %*% ( deltaTY ) ) %*% solve( SigmaY )
				# relative Bias (in percent)
				deltaBYrel <- deltaBY / BY * 100
				
				## Biases of autoregression matrix (CLPM as analysis model), simplified Eq. 11 
				deltaBYstar <- ( BYZ %*% SigmaZY ) %*% solve( SigmaY )
				# relative Bias (in percent)
				deltaBYstarrel <- deltaBYstar / BY * 100
				
			} else {
				deltaTY <- deltaBY <- deltaBYstar <- matrix( 0, F, F )
			}
			## Compensation rate
			CR <- ( deltaBYstar - deltaBY ) / deltaBYstar * 100
		
		} else {
			deltaTY <- deltaBY <- deltaBYrel <- deltaBYstar <- deltaBYstarrel <- CR <- NULL
		}
		
		# return object
		ret <- list( deltaTY, deltaBY, deltaBYrel, deltaBYstar, deltaBYstarrel, CR )
		names( ret ) <- c( "deltaTY", "deltaBY", "deltaBYrel", "deltaBYstar", "deltaBYstarrel", "CR" )		
		
		# return
		return( ret )
}

### development
# user.profile <- shell( "echo %USERPROFILE%", intern=TRUE )
# Rfiles.folder <- file.path( user.profile,
                                    # "Dropbox/96_mutualism/mutualism/R" )
# Rfiles <- list.files( Rfiles.folder , pattern="*.R" )
# Rfiles <- Rfiles[ !Rfiles %in% c("bias.formula.R") ]
# for( Rfile in Rfiles ){
	# source( file.path( Rfiles.folder, Rfile ) )
# }

## Autoregressive effects (for all variables)
# alpha <- 0.3

## Cross-lagged effects of the focal effect (X and Y)
# beta_focal <- 0.1

## Main confounding effects: Cross-lagged effects 
## (1) between X-related confounders and X
## (2) between Y-related confounders and Y
## (3) within X-related confounders
## (4) within Y-related confounders
## They are all set to be equal in this simulation but can be easily changed
# beta_covfocal <- 0.1

## Trivial confounding effects: Cross-lagged effects 
## (1) between X-related confounders and Y
## (2) between Y-related confoudners and X
## (3) between X-related confounder and Y-related confounders
##  They are all set to be equal in this simulation but can be easily changed
# beta_covnofocal <- 0.01

## Time-specific residual variance
# sigma <- 0.5

## Number of focal processes
# F <- 1

## Number of X-related confounders
# c1 <- 1

## Number of Y-related confounders
# c2 <- 0

## total number of confounders
# C <- c1 + c2

## total number of variables
# V <- C + F

# Autoregression matrix
# A <- matrix(NA, V, V)
# diag(A) <- alpha

# A[lower.tri(A, diag = FALSE)] <- beta_covnofocal 
# A[upper.tri(A, diag = FALSE)] <- beta_covnofocal

# bias.formula( A=A, F=F, c1=c1, c2=c2 )


### test
# require( testthat )
# test_file("../tests/testthat/XXXXX.R")
