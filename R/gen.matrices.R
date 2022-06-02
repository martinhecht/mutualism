## Changelog:
# MH 0.0.1 2022-06-02

## Documentation
#' @title
#' @description right now only for F=3
#' @param
#' @param
#' @param
#' @return


## Function definition
gen.matrices <- function(){
		
		# R package
		require( matrixcalc ) # is.positive.definite()
		require( R.utils ) # withTimeout()
		
		# ret <- structure( "message", class = c("try-error", "character") )
		# while( inherits( ret, "try-error" ) ){
		ret <- NULL
		while( is.null( ret ) ){
			ret <- try( withTimeout( .gen.matrices(), timeout = 1, onTimeout = "silent" ) )
		}
		
		return( ret )
}

.gen.matrices <- function(){

	# start time                                                         
	# start <- Sys.time()
	
	# number of processes
	F <- 3
	
	asG <- matrix(0,nrow=F,ncol=F)
	while( any( diag( asG ) < 0.1 ) | any( diag( asG ) > 10 ) ){

		A <- matrix( c( runif(1,0.25,0.75),runif(1,0,0.5),runif(1,0,0.5),
						runif(1,0,0.5),runif(1,0.25,0.75),runif(1,0,0.5),
						runif(1,0,0.5),runif(1,0,0.5),runif(1,0.25,0.75) ), nrow=F, ncol=F )
		
		G <- matrix(0,nrow=F,ncol=F)
		tries <- 1
		max.tries <- 100
		while( !is.positive.definite( G ) & !is.positive.definite( asG ) & tries < max.tries ){
				
				# cat( paste0( tries, " " ) ); flush.console()
				
				G <- matrix( c( runif(1,0.5,1.5),NA,NA,
								runif(1,-0.5,0.5),runif(1,0.5,1.5),NA,
								runif(1,-0.5,0.5),runif(1,-0.5,0.5),runif(1,0.5,1.5) ), nrow=F, ncol=F )
				G[ lower.tri(G) ] <- G[ upper.tri(G) ]
				
				asG <- irow( solve( diag(F^2) - A %x% A ) %*% row( G ) )
				# warum auch immer, asG wird nicht als symmetrische Matrix von is.positive.definite erkannt,
				# obwohl natürlich quasi symmetrisch
				asG[ lower.tri(asG) ] <- asG[ upper.tri(asG) ]
				tries <- tries + 1
		}
		if( tries >= max.tries ) asG <- matrix(0,nrow=F,ncol=F)

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

# while( TRUE ) { print( gen.matrices() ); flush.console() }


### test
# require( testthat )
# test_file("../tests/testthat/XXXXX.R")
