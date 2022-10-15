## Changelog:
# MH 0.0.9 2022-10-15: added F parameter, F<=4 recommended (then with max.tries=1000 seems to always work)
# MH 0.0.1 2022-06-02

## Documentation
#' @title
#' @description
#' @param
#' @param
#' @param
#' @return


## Function definition
gen.matrices <- function(F=3,verbose=FALSE){
		
		# R package
		require( matrixcalc ) # is.positive.definite()
		require( R.utils ) # withTimeout()
		
		# ret <- structure( "message", class = c("try-error", "character") )
		# while( inherits( ret, "try-error" ) ){
		# ret <- NULL
		# while( is.null( ret ) ){
			# ret <- try( withTimeout( .gen.matrices(F=F,verbose=verbose), timeout = 10, onTimeout = "silent" ) )
		# }
		
		ret <- .gen.matrices(F=F,verbose=verbose)
		
		return( ret )
}

.gen.matrices <- function(F,verbose){

	# start time                                                         
	# start <- Sys.time()
	
	# number of processes
	# F <- 3
	
	G <- matrix(0,nrow=F,ncol=F)
	asG <- matrix(0,nrow=F,ncol=F)
	# while( any( diag( asG ) < 0.1 ) | any( diag( asG ) > 10 ) | is.positive.definite( asG ) ){
	# while( !is.positive.definite( G ) || !is.positive.definite( asG ) ){
	tries <- 1
	max.tries <- 1000
	while( !is.positive.definite( G ) || !is.positive.definite( asG ) && (tries <= max.tries) ){

		if( verbose ) { cat( paste0( tries, "\n" ) ); flush.console() }

		# A <- matrix( c( runif(1,0.25,0.75),runif(1,0,0.5),runif(1,0,0.5),
						# runif(1,0,0.5),runif(1,0.25,0.75),runif(1,0,0.5),
						# runif(1,0,0.5),runif(1,0,0.5),runif(1,0.25,0.75) ), nrow=F, ncol=F )
		# off-diagonals
		A <- matrix( sapply( 1:(F^2), function(x) runif(1,0,0.5) ), nrow=F,ncol=F )
		# diagonals
		for( f in 1:F ) A[f,f] <- runif(1,0.25,0.75)
			
		G <- matrix(0,nrow=F,ncol=F)
		asG <- matrix(0,nrow=F,ncol=F)

		# while( ( !is.positive.definite( G ) | !is.positive.definite( asG ) ) & !(tries < max.tries) ){
		# while( ( !is.positive.definite( G ) || !is.positive.definite( asG ) ) && (tries <= max.tries) ){
				

				
				# G <- matrix( c( runif(1,0.5,1.5),NA,NA,
								# runif(1,-0.5,0.5),runif(1,0.5,1.5),NA,
								# runif(1,-0.5,0.5),runif(1,-0.5,0.5),runif(1,0.5,1.5) ), nrow=F, ncol=F )
				
				# off-diagonals
				G <- matrix( sapply( 1:(F^2), function(x) runif(1,-0.5,0.5) ), nrow=F,ncol=F )
				# diagonals
				for( f in 1:F ) G[f,f] <- runif(1,0.5,1.5)
				# mirror
				G[ upper.tri(G) ] <- t(G)[ upper.tri(G) ]
				
				asG <- irow( solve( diag(F^2) - A %x% A ) %*% row( G ) )
				# warum auch immer, asG wird nicht als symmetrische Matrix von is.positive.definite erkannt,
				# obwohl natürlich quasi symmetrisch
				asG[ upper.tri(asG) ] <- t(asG)[ upper.tri(asG) ]
				
				tries <- tries + 1

		# }
		# if( tries >= max.tries ) asG <- matrix(0,nrow=F,ncol=F)

	}
	# run time                                                           
	# print( runtime <- Sys.time() - start )

	ret <- list( A, G, asG )
	names( ret ) <- c( "A", "G", "asG" )

    return( ret )
}

### development

# user.profile <- shell( "echo %USERPROFILE%", intern=TRUE )
# Rfiles.folder <- file.path( user.profile,
                                    # "Dropbox/96_mutualism/mutualism/R" )
# Rfiles <- list.files( Rfiles.folder , pattern="*.R" )
# Rfiles <- Rfiles[ !Rfiles %in% c("gen.matrices.R", ".gen.matrices.R") ]
# for( Rfile in Rfiles ){
	# source( file.path( Rfiles.folder, Rfile ) )
# }

# print( m <- gen.matrices(F=4,verbose=TRUE) )
# is.positive.definite( m$G )
# is.positive.definite( m$asG )

# while( TRUE ) { print( m <- gen.matrices(F=4,verbose=TRUE) ); print( is.positive.definite( m$G ) ) ; print( is.positive.definite( m$asG ) ); flush.console() }


### test
# require( testthat )
# test_file("../tests/testthat/XXXXX.R")
