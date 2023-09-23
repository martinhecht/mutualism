## Changelog:
# MH 0.0.15 2023-09-23: RI-CLPM, now with constraints for random process means, RI-CLPM2/3 discarded!!
# MH 0.0.14 2023-09-22: RI-CLPM3, intercept set to 1 for t >= 2, freely estimated parameter for intercept (=mean) at t=1
# MH 0.0.11 2023-05-08: RI-CLPM2
# MH 0.0.10 2023-05-08: bug fixes in RI-CLPM (first time point variance now free, b now not into first time point) and removed latent
# MH 0.0.8 2022-10-14: initial programming

## Function definition
#### only F (number of processes) = 2 !!!
#### t1.stationary: only FALSE !!!
gen.lavaan.syntax <- function( model=c("CLPM","RI-CLPM"), T=2, t1.stationary=FALSE ){
		
		# cautionary note
		if( model %in% "RI-CLPM" ) cat( "\nNote: Starting from version 0.0.15, the RI-CLPM model has changed (now with constraints for random process means)!\n      Be cautious when using old scripts that depend on versions <= 0.0.14\n\n" )
		
		syn <- NULL

		if( model %in% c("RI-CLPM") ){
			if( !t1.stationary ){
				# RICLPM lavaan syntax
				syn <- "### Random intercepts cross-lagged panel model ###"
				syn <- c( syn, "#   (random process means implementation)\n" )
				# MH 0.0.15 2023-09-23: new RI-CLPM (with constraints for random process means)
				syn <- c(syn, "# Random process means" )

				syn <- c(syn, paste( c("mu1 =~ 1*x.1", paste0( "con1a*x.",2:T ), paste0( "con1b*y.",2:T )), collapse=" + " ),
							  paste( c("mu2 =~ 1*y.1", paste0( "con2a*y.",2:T ), paste0( "con2b*x.",2:T )), collapse=" + " ) )
				syn <- c(syn, c( "con1a == (1-a1)",
								 "con2a == (1-a2)",
								 "con1b == -1*a21",
								 "con2b == -1*a12" ) )
				
				syn <- c(syn, paste0("# Random process means (co)variances\n",
									 "mu1 ~~ mu1var*mu1\n",
									 "mu2 ~~ mu2var*mu2\n",
									 "mu1 ~~ mucov*mu2" )	)				
				syn <- c(syn, paste0("# First time point (co)variances\n",
									 "x.1 ~~ f1var*x.1\n",
									 "y.1 ~~ f2var*y.1\n",
									 "x.1 ~~ fcov*y.1" ) )
				syn <- c(syn, "# Autoregressive and cross-lagged effects" )
				for( t in 2:T ){
					syn <- c(syn,	paste0("x.",t," ~ a1*x.",t-1," + a12*y.",t-1,"\n",
										   "y.",t," ~ a2*y.",t-1," + a21*x.",t-1 ) )
				}
				syn <- c(syn, "# Process error (co)variances" )
				for( t in 2:T ){
					syn <- c(syn,	paste0("x.",t," ~~ e1var*x.",t,"\n",
										   "y.",t," ~~ e2var*y.",t,"\n",
                                           "x.",t," ~~ ecov*y.",t ) )
				}				

				syn <- paste( syn, collapse="\n" )
			}
		}
		
		if( model %in% "CLPM" ){
			if( !t1.stationary ){		
				# data analysis
				syn <- paste0("# First time point (co)variances\n",
							  "x.1 ~~ fvar1*x.1\n",
							  "y.1 ~~ fvar2*y.1\n",
							  "x.1 ~~ fcov*y.1")
				syn <- c(syn, "# Autoregressive and cross-lagged effects" )
				for( t in 2:T ){
					syn <- c(syn,	paste0("x.",t," ~ a1*x.",t-1," + a12*y.",t-1,"\n",
										   "y.",t," ~ a2*y.",t-1," + a21*x.",t-1   ) )
				}
				syn <- c(syn, "# Process error (co)variances" )
				for( t in 2:T ){
					syn <- c(syn,	paste0("x.",t," ~~ evar1*x.",t,"\n",
										   "y.",t," ~~ evar2*y.",t,"\n",
										   "x.",t," ~~ ecov*y.",t) )
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
# syn <- gen.lavaan.syntax( "RI-CLPM", T=3 )
# cat(syn)
# cat( gen.lavaan.syntax( "RI-CLPM3", T=4 ) )
# cat( gen.lavaan.syntax( "CLPM", T=4 ) )


### test
# require( testthat )
# test_file("../tests/testthat/XXXXX.R")
