## Changelog:
# MH 0.0.8 2022-10-14: initial programming

## Function definition
#### only F (number of processes) = 2 !!!
#### t1.stationary: only FALSE !!!
gen.lavaan.syntax <- function( model=c("CLPM","RI-CLPM"), T=2, t1.stationary=FALSE ){
		
		syn <- NULL
		
		if( model %in% "RI-CLPM" ){
			if( !t1.stationary ){
				# RICLPM lavaan syntax
				# https://jeroendmulder.github.io/RI-CLPM/lavaan.html
				syn <- "# Create between components (random intercepts)"
				syn <- c(syn,	paste0( "RIx =~", paste(paste0( "1*x.",1:T ), collapse=" + " ), "\n",
										"RIy =~", paste(paste0( "1*y.",1:T ), collapse=" + " )  )
							)
				syn <- c(syn, "# Create within-person centered variables" )
				for( t in 1:T ){
					syn <- c(syn,	# latent
									paste0( "fx.",t," =~ 1*x.",t, "\n",
											"fy.",t," =~ 1*y.",t )
							)
				}
				syn <- c(syn, "# Estimate the lagged effects between the within-person centered variables" )
				for( t in 2:T ){
					syn <- c(syn,	paste0("fx.",t," ~ a1*fx.",t-1," + a12*fy.",t-1,"\n",
										   "fy.",t," ~ a2*fy.",t-1," + a21*fx.",t-1   ),
									paste0("fx.",t," ~~ evar1*fx.",t,"\n",
										   "fy.",t," ~~ evar2*fy.",t,"\n",
										   "fx.",t," ~~ ecov*fy.",t)
									
										   )
				}
				syn <- c(syn, paste0("# first time point\n",
									 "fx.1 ~~ fvar1*x.1\n",
									 "fy.1 ~~ fvar2*y.1\n",
									 "fx.1 ~~ fcov*y.1\n",
									 "# Estimate the variance and covariance of the random intercepts\n",
									 "RIx ~~ RIxVar*RIx\n",
									 "RIy ~~ RIyVar*RIy\n",
									 "RIx ~~ RIxyCov*RIy" )	)
				syn <- paste( syn, collapse="\n" )
			
			}
		}
		
		if( model %in% "CLPM" ){
			if( !t1.stationary ){		
				# data analysis
				syn <- "x.1 ~~ fvar1*x.1
						y.1 ~~ fvar2*y.1
						x.1 ~~ fcov*y.1"
				for( t in 2:T ){
					syn <- c(syn,	paste0("x.",t," ~ a1*x.",t-1," + a12*y.",t-1,"\n",
										   "y.",t," ~ a2*y.",t-1," + a21*x.",t-1   ),
									paste0("x.",t," ~~ evar1*x.",t,"\n",
										   "y.",t," ~~ evar2*y.",t,"\n",
										   "x.",t," ~~ ecov*y.",t)
									
										   )
				}
				syn <- paste( syn, collapse="\n" )
			}
		}
		
		# return
		return( syn )
}

### development
# user.profile <- shell( "echo %USERPROFILE%", intern=TRUE )
# Rfiles.folder <- file.path( user.profile,
                                    # "Dropbox/96_mutualism/mutualism/R" )
# Rfiles <- list.files( Rfiles.folder , pattern="*.R" )
# Rfiles <- Rfiles[ !Rfiles %in% c("gen.lavaan.syntax.R") ]
# for( Rfile in Rfiles ){
	# source( file.path( Rfiles.folder, Rfile ) )
# }

### 
# cat( gen.lavaan.syntax( "RI-CLPM", T=4 ) )
# cat( gen.lavaan.syntax( "CLPM", T=4 ) )


### test
# require( testthat )
# test_file("../tests/testthat/XXXXX.R")
