## Changelog:
# MH 0.0.17 2025-04-11: bugfix var.names (for >26 processes)
# MH 0.0.2 2022-09-22
# MH 0.0.1 2022-06-01

## Documentation
#' @title
#' @description
#' @param A autoregression matrix
#' @param Q process error covariance matrix
#' @param G first time point covariance matrix covariance matrix
#' @param stationary (logical), TRUE (default): G is constrained to process covariance matrix, FALSE: G needs to be specified
#' @param N number of persons
#' @param T number of time points
#' @param seed "random" or number
#' @return

## Function definition
data.gen.clpm.riclpm <- function( A, Q, G=NULL, B=NULL, stationary=TRUE, N=1, T=2, seed="random" ){
		
		# packages
		require( MASS ) # mvrnorm()
		require( reshape2 ) # melt()
		
		# F (number of processes)
		F <- dim( A )[1]
		
		# constrain G
		if( stationary ){
			G <- irow( solve( diag(F^2) - A %x% A ) %*% row( Q ) )
		}
		
		# seed
		if( seed %in% "random" ){
			seed <- sample( 1:999999999, 1 )
		}
		set.seed( seed )
		
		## initialization of output array
		x <- array( NA, dim=c(F,T,N) )
		# indiviudal means (for RI-CLPM)
		if( !is.null( B ) ){
			muj <- array( NA, dim=c(F,N) )

			# term used for t >= 2 data generation with individual means
			# (compute that once here to avoid computing it in the heave loop)
			IminusA <- (diag(F) - A)
		}

		## data generation
		
		# indiviudal means (for RI-CLPM)
		if( !is.null( B ) ){
			for( n in 1:N ){
				muj[,n] <- mvrnorm( 1, mu=rep(0,F), Sigma=B )
			}
		}
		
		# time series
		for( n in 1:N ){
			for( t in 1:T ){
				# first time point
				if( t %in% 1 ){
					x[,1,n] <- mvrnorm( 1, mu=rep(0,F), Sigma=G )
					# indiviudal means (for RI-CLPM)
					if( !is.null( B ) ){
						x[,1,n] <- x[,1,n] + muj[,n]
					}
				# t >= 2
				} else {
					x[,t,n] <- A %*% x[,t-1,n] + mvrnorm( 1, mu=rep(0,F), Sigma=Q )
					# indiviudal means (for RI-CLPM)
					if( !is.null( B ) ){
						x[,t,n] <- x[,t,n] + IminusA %*% muj[,n]
					}
				}
			}
		}		

		## reshapen
		# very long ("long long")
		dll <- melt(x)
		colnames( dll ) <- c( "var", "t", "n", "value" )
		# long
		dl <- dcast( dll, n + t ~ var, value.var="value" )
		let <- c( letters[24:26], letters[1:23] )
		# MH 0.0.17 2025-04-11: bugfix var.names (for >26 processes)
		#var.names <- paste0( "x", ncol( dl ) - 2 ) changed to:
		#var.names <- c( "x", "y", paste0( "z", 1:(ncol( dl ) - 4) ) )
		if( ncol( dl ) <= 28 ) var.names <- let[ 1:(ncol( dl ) - 2) ] else var.names <- c( "x", "y", paste0( "z", 1:(ncol( dl ) - 4) ) )
		colnames( dl )[ 3:ncol( dl ) ] <- var.names
		# wide
		dw.l <- sapply( var.names, function( vn ) {d<-dcast( dl, n ~ t, value.var=vn ); colnames(d)[-1] <- paste0(vn,".",colnames(d)[-1]); return(d)}, simplify=FALSE )
		dw <- Reduce(function(x, y) merge(x, y,all=TRUE,by="n"),dw.l,accumulate=FALSE )

		# return object
		ret <- list( dll, dl, dw, seed )
		names( ret ) <- c( "dll", "dl", "dw", "seed" )
		
		# return
		return( ret )
}

### development
# user.profile <- shell( "echo %USERPROFILE%", intern=TRUE )
# Rfiles.folder <- file.path( user.profile,
                                    # "Dropbox/96_mutualism/mutualism/R" )
# Rfiles <- list.files( Rfiles.folder , pattern="*.R" )
# Rfiles <- Rfiles[ !Rfiles %in% c("data.gen.clpm.R") ]
# for( Rfile in Rfiles ){
	# source( file.path( Rfiles.folder, Rfile ) )
# }

### test1a
# F <- 3
# A <- matrix( c( 0.5,0,0, 0,0.5,0, 0.5,0.5,0.5 ), nrow=F, ncol=F )
# Q <- matrix( c( 1,0.5,0.5, 0.5,1,0.5, 0.5,0.5,1 ), nrow=F, ncol=F )
# ( dl <- data.gen.clpm.riclpm( A, Q, N=5, T=10 )$dl )

### test1b
# G <- matrix( c( 2,1,1, 1,2,1, 1,1,2 ), nrow=F, ncol=F )
# ( dl <- data.gen.clpm.riclpm( A, Q, G, stationary=FALSE, N=5, T=10 )$dl )

### test1c
# B <- matrix( c( 2,1,1, 1,2,1, 1,1,2 ), nrow=F, ncol=F )
# ( dl <- data.gen.clpm.riclpm( A, Q, G, B, stationary=FALSE, N=5, T=10 )$dl )


### test
# require( testthat )
# test_file("../tests/testthat/XXXXX.R")
