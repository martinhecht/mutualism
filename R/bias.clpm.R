## Changelog:
# MH 0.0.3 2022-06-06: minor aesthetic fixes, nothing core functional
# MH 0.0.2 2022-06-02: updated based on clpm_bias_v2.pdf 2022-06-02
# MH 0.0.1 2022-05-31: based on clpm_bias_v1.pdf 2022-05-29

## Documentation
#' @title
#' @description
#' @param Afull autoregression matrix of the full CLPM
#' @param CovMatrixfull either first time point covariance matrix (stationary=FALSE) or process error covariance matrix (stationary=TRUE) of the full CLPM
#' @param F1 number of processes of reduced model
#' @param F2 number of omitted processes in the reduced model
#' @param stationary (logical), FALSE (default): CovMatrixfull is first time point covariance matrix, TRUE: CovMatrixfull is process error covariance matrix
#' @param extended.results (logical), FALSE (default): bias as matrix is returned, TRUE: additional results (including bias in long format) is returned
#' @return

## Function definition
bias.clpm <- function( Afull, CovMatrixfull, F1, F2, stationary=TRUE, extended.results=FALSE ){
		
		# F
		F <- F1 + F2
		
		# Afull submatrices
		A11 <- Afull[1:F1,1:F1,drop=FALSE]
		A12 <- Afull[1:F1,(F1+1):F,drop=FALSE]

		## CovMatrix
		# if stationary=FALSE: covariance matrix of first time point
		#  if stationary=TRUE: covariance matrix of process errors
		if (!stationary){
			Gfull <- CovMatrixfull
		} else {
			Gfull <- irow( solve( diag(F^2) - Afull %x% Afull ) %*% row( CovMatrixfull ) )
		}
		
		# Gfull submatrices
		G11 <- Gfull[1:F1,1:F1,drop=FALSE]
		G21 <- Gfull[(F1+1):F,1:F1,drop=FALSE]		
		
		## bias
		# c:\Users\martin\Dropbox\96_mutualism\03_models\07_clpm_matr3T\01_bias\clpm_bias_v1.pdf (2022-05-29)
		#bias.first.term <- A11 %*% ( G11 %*% - diag( F1 ) )
		#bias.second.term <- A12 %*% G21 %*% solve( G11 )
		#bias <- bias.first.term + bias.second.term
		# 0.0.2 2022-06-02: updated based on clpm_bias_v2.pdf 2022-06-02
		# c:\Users\martin\Dropbox\96_mutualism\03_models\07_clpm_matr3T\01_bias\clpm_bias_v2.pdf (2022-06-02)
		bias <- A12 %*% G21 %*% solve( G11 )
		
		# return object
		if( !extended.results ) {
			ret <- bias
		} else {
			# R package
			require( reshape2 ) # melt()
			
			## long format
			# bias
			bias.long <- melt( bias )
			# true
			true <- melt( A11 )
			# merge
			bias.long <- merge( bias.long, true, by=c("Var1","Var2") )
			colnames( bias.long ) <- c("row", "col", "bias", "true")
			# value
			bias.long$value <- bias.long$true + bias.long$bias
			# type
			bias.long$type <- "cross-lagged effect"
			bias.long$type[ bias.long$row==bias.long$col ] <- "autoregressive effect"
			# parameter label
			bias.long$par <- NA
			bias.long$par[ bias.long$type %in% "autoregressive effect" ] <- paste0( "ar", bias.long$row[ bias.long$type %in% "autoregressive effect" ] )
			bias.long$par[ bias.long$type %in% "cross-lagged effect" ] <- paste0( "cl", bias.long$col[ bias.long$type %in% "cross-lagged effect" ], "}", bias.long$row[ bias.long$type %in% "cross-lagged effect" ] )

			bias.long <- bias.long[,c("type","par","row", "col", "bias", "true", "value")]
			bias.long <- bias.long[ order( bias.long$type ) , ]
			rownames( bias.long ) <- seq( along=rownames( bias.long ) )
			
			# return list
			# ret <- list( bias, bias.first.term, bias.second.term, bias.long )
			ret <- list( bias, bias.long )
			names( ret ) <- c( "bias", "bias.long" )
		}
		
		# return
		return( ret )
}

### development
# user.profile <- shell( "echo %USERPROFILE%", intern=TRUE )
# Rfiles.folder <- file.path( user.profile,
                                    # "Dropbox/96_mutualism/mutualism/R" )
# Rfiles <- list.files( Rfiles.folder , pattern="*.R" )
# Rfiles <- Rfiles[ !Rfiles %in% c("bias.clpm.R") ]
# for( Rfile in Rfiles ){
	# source( file.path( Rfiles.folder, Rfile ) )
# }

### test1
# F <- 3
# A <- matrix( c( 0.5,0,0, 0,0.5,0, 0.5,0.5,0.5 ), nrow=F, ncol=F )
# G <- matrix( c( 1,0.5,0.5, 0.5,1,0.5, 0.5,0.5,1 ), nrow=F, ncol=F )
# bias.clpm( A, G, 2, 1 )
# bias.clpm( A, G, 2, 1, extended.results=TRUE )
# bias.clpm( A, G, 2, 1, stationary=FALSE )

### test2
# F <- 4
# A <- matrix( c( 0.5,0,0,0, 0,0.5,0,0, 0.25,0.25,0.5,0, 0.25,0.25,0.5,0.5 ), nrow=F, ncol=F )
# G <- matrix( c( 1,0.5,0.5,0.5, 0.5,1,0.5,0.5, 0.5,0.5,1,0.5, 0.5,0.5,0.5,1 ), nrow=F, ncol=F )
# bias.clpm( A, G, 2, 2, extended.results=TRUE )

### test
# require( testthat )
# test_file("../tests/testthat/XXXXX.R")
