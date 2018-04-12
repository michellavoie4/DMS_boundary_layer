#---------------------------------------------------------------------------#
# Modeling cell surface DMS enrichment                                      #
# Fig. 2. 3 and 4 of Lavoie et al. (2018)                                   #
#                                                                           #
# Lavoie M. 2018                                                            #
#---------------------------------------------------------------------------#

#---------------------------------------------------------------------------#
# Fig. 2 Modeling cell surface DMS enrichment in low DMS producers          #
#                                                                           #
# Lavoie M. 2018                                                            #
#---------------------------------------------------------------------------#

graphics.off()
rm(list=ls())

#------------#
# Parameters #
#------------#

J1    <- 2.8e-16               # DMS flux for I. galbana (mol DMS cm-2 s-1) (Spiese et al. 2009)
r1    <- 2.4e-4                # Cell radius of I. galbana (cm) (Spiese et al. 2009)
Area  <- 4 * pi * r1^2         # Area in cm^2
Vol   <- 4/3 * pi * (1e4*r1)^3 # Volume in um^3
Vol_L <- 1e-15 * Vol           # Volume in Liter of biovolume
R     <- J1 * Area / Vol_L     # production rate in mol DMS Lcell-1 s-1

#------------------#
# Model parameters #
#------------------#

D  <- 1.19e-5         # DMS diffusion coefficient at 20 degC (cm^2 s^-1) (Saltzman et al. 1993)
E  <- c( 0.1, 1, 10 ) # Shear rate (s^-1)
v  <- 1e-2            # Seawater kinematic viscosity at 20 degC (cm^2 s^-1) (Karp-Boss et al. 1996) 


#-----------#
# Variables #
#-----------#

x <- seq( 2e-4, 5e-2, length=1000 ) # cell radii (cm)
J <- R/1e3 * x/3                   # DMS flux (mol cm^-2 s^-1); (1/1e3) is converting L of the volume in cm^3


#---------------#
# Model Results #
#---------------#

# Reynolds number (no units) (Karp-Boss et al 1996)
Re  <- outer( x^2, E/v )

# Peclet number (no units) (Karp-Boss et al 1996)
Pe  <- outer( x^2, E/D)

# Sherwood number (no units) (Karp-Boss et al 1996)
Sh <- Re * 0

i     <- Re<1e-1 & Pe<1e-2
j     <- Re<3e-1 & Pe>=1e-2 & Pe<1e2
k     <- Re<1    & Pe>=1e2
l     <- Re>=1

Sh[i] <- 1     + 0.29 * Pe[i]^(1/2)
Sh[j] <- 1.014 + 0.15 * Pe[j]^(1/2)
Sh[k] <-         0.55 * Pe[k]^(1/3)
Sh[l] <- NA

# Boundary layer thickness (cm) (Spiese et al 2016)
delta  <- x / Sh

# DMS cell surface enrichment (mol L^-1)
y <- (1e3 * J * x / D) * (1 - (x / (x + delta)))

#------------------#
# Plotting output  #
#------------------#

x     <- x  * 1e4 # Convert cell radius in um

ri    <- which.min( abs( x - r1*1e4 ) )
r1_um <- x[ri]    # Radius closest to r1


#--- Graph log-log

#tiff("Fig2_LowDMS1.tiff")                                        

par( mar = c(5, 5, 2, 2) )

# Solid, dashed and dotted lines correspond to the three shear rates from vector E, respectively.

plot(  x, y[,1], xlim=c(1, 1000), ylim=c(1e-12, 1e-6), xlab=expression(paste("cell radius ", mu, "m")), ylab=expression("DMS enrichment " ~  (mol ~ L^ ~ {-1})), log="xy", type="l", xaxp=c(1, 1000, 1), xaxt="n", yaxp=c(1e-12, 1e-6, 1), yaxt="n")
lines( x, y[,2], lty=2 )
lines( x, y[,3], lty=3 )
points( r1_um, y[ri,1], pch=19, cex=1)
axis( 1, at=c(1, 1e1, 1e2, 1e3) )
axis( 2, at=c(1e-12, 1e-11, 1e-10, 1e-9, 1e-8, 1e-7, 1e-6) )
axis( 1, at=c(2:9/1, 2:9/0.1, 2:9/0.01, 2:9/0.001), labels=FALSE, tcl=-.2) # tcl=-0.2 for minor ticks
axis( 2, at=c(2:9/1e12, 2:9/1e11, 2:9/1e10, 2:9/1e9, 2:9/1e8, 2:9/1e7), labels=FALSE, tcl=-.2)

# Median DMS concentration in the surface ocean and interquartile range (Lana et al 2011)
abline(h=1.9E-09, lty=4)  
rect(0.1,1.14E-09,2000,2.85E-09, col = rgb(0.5,0.5,0.5,0.1), lty=0)

#dev.off()

#-------------#
# References  #  
#-------------#

# Karp-Boss L, Boss E, Jumars PA (1996). Nutrient fluxes to planktonic osmotrophs in the presence of fluid motion. Oceanography and Marine Biology : an Annual Review 34: 71-107.
# Lana A, Bell TG, Simó R, Vallina SM, Ballabrera-Poy J, Kettle AJ, Dachs J, Bopp L, Saltzman ES, Stefels J, Johnson JE, Liss PS (2011). An updated climatology of surface dimethlysulfide concentrations and emission fluxes in the global ocean. Global Biogeochem Cycles 25: GB1004
# Saltzman ES, King DB, Holmen K, Leck C (1993). Experimental determination of the diffusion coefficient of dimethylsulfide in water. J Geophys Res, C 98: 16481-16486.
# Spiese CE, Kieber DJ, Nomura CT, Kiene RP (2009). Reduction of dimethylsulfoxide to dimethylsulfide by marine phytoplankton. Limnol Oceanogr 54: 560-570.
# Spiese CE, Le T, Zimmer RL, Kieber DJ (2016). Dimethylsulfide membrane permeability, cellular concentrations and implications for physiological functions in marine algae. J Plankton Res 38: 41-54.


#------------------------------------------------------------------------------------------#
# Fig. 3 Modeling cell surface DMS enrichment in High N replete/N-limited DMS producers    #
#                                                                                          #
# Lavoie M. 2018                                                                           #
#------------------------------------------------------------------------------------------#

#---------------------------------------#
# Scenario 3 a) Assuming J scales with r #
#---------------------------------------#

#------------------#
# Model parameters #
#------------------#

J1 <- 8.6E-16                    # DMS flux for E.huxleyi (N replete) (mol DMS cm-2 s-1) (Sunda et al 2007)
r1 <- 2E-04                     # cell radius of E. huxleyi (cm) (Sunda et al 2007)
Area <- 4*pi*(r1^2)              # Area in cm^2
Vol <- (4/3)*pi*((10000*r1)^3)   # Volume in um^3
Vol_L <- 1E-15*Vol               # Volume in Liter of biovolume
R1 <- J1 * (Area/Vol_L)           # production rate in mol DMS Lcell-1 s-1

J2 <- 4.12E-15                    # DMS flux for E huxleyi N -limited (mol DMS cm-2 s-1) (Sunda et al 2007)
r2 <- 2.18E-04                     # cell radius of E. huxleyi (cm) (Sunda et al 2007)
Area2 <- 4*pi*(r2^2)               # Area in cm^2
Vol2 <- (4/3)*pi*((10000*r2)^3)    # Volume in um^3
Vol_L2 <- 1E-15*Vol2               # Volume in Liter of biovolume
R2 <- J2 * (Area2/Vol_L2)          # production rate in mol DMS Lcell-1 s-1

R <- c ( R1, R2 )

D <- 1.19E-05                    # DMS diffusion coefficient (cm^2 /s) (Saltzman et al. 1993)
E <- 0.1                       # Shear rates (s-1)
v <- 1E-02                       # Seawater kinematic viscosity (cm^2 s-1) (Karp-Boss et al. 1996)

#-----------#
# Variables #
#-----------#

x <- seq(2E-04, 5E-02, length=1000)  # cell radii (cm)
J <- outer (x/3, R/1e3)               # DMS flux (mol cm-2 s-1)  We assume that J scales with r

#-----------------#
# Model Functions #
#-----------------#

# Reynolds number (no units) (Karp-Boss et al 1996)
Re  <- x^2 * E/v

# Peclet number (no units) (Karp-Boss et al 1996)
Pe  <- x^2 * E/D

# Sherwood number (no units) (Karp-Boss et al 1996)
Sh <- Re * 0

i     <- Re<1e-1 & Pe<1e-2
j     <- Re<3e-1 & Pe>=1e-2 & Pe<1e2
k     <- Re<1    & Pe>=1e2
l     <- Re>=1

Sh[i] <- 1     + 0.29 * Pe[i]^(1/2)
Sh[j] <- 1.014 + 0.15 * Pe[j]^(1/2)
Sh[k] <-         0.55 * Pe[k]^(1/3)
Sh[l] <- NA

# Boundary layer thickness (cm) (Spiese et al 2016)
delta  <-  x / Sh 

# DMS cell surface enrichment (mol L^-1)
y <- (1e3 * J * x / D) * (1 - (x / (x + delta)))

#------------------#
# Plotting output  #
#------------------#

x     <- x  * 1e4 # Convert cell radius in um

ri    <- which.min( abs( x - r1*1e4 ) )
r1_um <- x[ri]    # Radius closest to r1

ri_2    <- which.min( abs( x - r2*1e4 ) )
r2_um <- x[ri_2]    # Radius closest to r2

#------------
# theoretical r^2 relationship
y1 <- 1e-10 * x^2

# theoretical linear relationship
y2 <- 3e-9 * x

#--- Graph log-log

#tiff("Fig3a_HighDMS.tiff")                         

par( mar = c(5, 5, 2, 2) )

# Solid, dashed and dotted lines correspond to the three shear rates from vector E, respectively.

plot(  x, y[, 1], xlim=c(1, 1000), ylim=c(1e-12, 1e-6), xlab=expression(paste("cell radius ", mu, "m")), ylab=expression("DMS enrichment " ~  (mol ~ L^ ~ {-1})), log="xy", type="l", xaxp=c(1, 1000, 1), xaxt="n", yaxp=c(1e-17, 1e-6, 1), yaxt="n")
lines(  x, y[, 2], lty=2)
points( r1_um, y[ri, 1], pch=19, cex=1)
points( r2_um, y[ri_2, 2], pch=1, cex=1)
lines(x[1:10], y1[1:10], lty=3)   # Add a theoretical r^2 relationship
text(x=2, y=7E-10, labels=expression(r^{2}), srt = 0, cex=2)
lines(x[1:10], y2[1:10], lty=3)   # Theoretical linear relationship
text(x=2, y=1E-08, labels=expression(r^{1}), srt = 0, cex=2)
axis( 1, at=c(1, 1e1, 1e2, 1e3) )
axis( 2, at=c(1e-12, 1e-11, 1e-10, 1e-9, 1e-8, 1e-7, 1e-6) )
axis( 1, at=c(2:9/1, 2:9/0.1, 2:9/0.01, 2:9/0.001), labels=FALSE, tcl=-.2) # tcl=-0.2 for minor ticks
axis( 2, at=c(2:9/1e12, 2:9/1e11, 2:9/1e10, 2:9/1e9, 2:9/1e8, 2:9/1e7), labels=FALSE, tcl=-.2)

# Median DMS concentration in the surface ocean and interquartile range (Lana et al 2011)
abline(h=1.9E-09, lty=4)
rect(0.1,1.14E-09,2000,2.85E-09, col = rgb(0.5,0.5,0.5,0.1), lty=0)

#dev.off()

#--------------------------------------------------------#
# Scenario 3b) Assuming a constant J for all cell radius #
#--------------------------------------------------------#

#------------------#
# Model parameters #
#------------------#

J1 <- 8.6E-16                    # DMS flux for E.huxleyi (N replete) (mol DMS cm-2 s-1) (Sunda et al)
r1 <- 2E-04                     # cell radius of E. huxleyi (cm) (Sunda et al 2007)

J2 <- 4.12E-15                    # DMS flux for E huxleyi N -limited (mol DMS cm-2 s-1) (Dacey and Wakeham 1986)
r2 <- 2.18E-04                     # cell radius of E. huxleyi (cm) (Sunda)

D <- 1.19E-05          # DMS diffusion coefficient (cm^2 /s)
E <- 0.1               # Shear rate (s-1)
v <- 1E-02             # Seawater kinematic viscosity (cm^2 s-1) (Karp-Boss et al 1996)

#-----------#
# Variables #
#-----------#

x <- seq(2E-04, 5E-02, length=1000)  # cell radii (cm)
J <- c(J1,J2)                              # DMS flux (mol cm-2 s-1)  We assume that J is constant

#-----------------#
# Model Functions #
#-----------------#

# Reynolds number (no units) (Karp-Boss et al 1996)
Re  <- x^2 * E/v

# Peclet number (no units) (Karp-Boss et al 1996)
Pe  <- x^2 * E/D

# Sherwood number (no units) (Karp-Boss et al 1996)
Sh <- Re * 0

i     <- Re<1e-1 & Pe<1e-2
j     <- Re<3e-1 & Pe>=1e-2 & Pe<1e2
k     <- Re<1    & Pe>=1e2
l     <- Re>=1

Sh[i] <- 1     + 0.29 * Pe[i]^(1/2)
Sh[j] <- 1.014 + 0.15 * Pe[j]^(1/2)
Sh[k] <-         0.55 * Pe[k]^(1/3)
Sh[l] <- NA

# Boundary layer thickness (cm) (Spiese et al 2016)
delta  <- x / Sh

# DMS cell surface enrichment (mol L^-1)
y <- outer( x /D * (1 - (x / (x + delta))) , 1e3 * J )

#------------------#
# Plotting output  #
#------------------#

x     <- x  * 1e4 # Convert cell radius in um

ri    <- which.min( abs( x - r1*1e4 ) )
r1_um <- x[ri]    # Radius closest to r1

ri_2    <- which.min( abs( x - r2*1e4 ) )
r2_um <- x[ri_2]    # Radius closest to r2

#------------
# theoretical r^2 relationship
y1 <- 1e-10 * x^2

# theoretical linear relationship
y2 <- 3e-9 * x

#--- Graph log-log

#tiff("Fig3b_cstJ_HighDMS.tiff")                              

par( mar = c(5, 5, 2, 2) )

# Solid, dashed and dotted lines correspond to the three shear rates from vector E, respectively.

plot(  x, y[,1], xlim=c(1, 1000), ylim=c(1e-12, 1e-6), xlab=expression(paste("cell radius ", mu, "m")), ylab=expression("DMS enrichment " ~  (mol ~ L^ ~ {-1})), log="xy", type="l", xaxp=c(1, 1000, 1), xaxt="n", yaxp=c(1e-12, 1e-6, 1), yaxt="n")

lines( x, y[,2], lty=2 )
points( r1_um, y[ri,1], pch=19, cex=1)
points( r2_um, y[ri_2, 2], pch=1, cex=1)
lines(x[1:10], y1[1:10], lty=3)   # Add a theoretical r^2 relationship
text(x=2, y=7E-10, labels=expression(r^{2}), srt = 0, cex=2)
lines(x[1:10], y2[1:10], lty=3)   # Theoretical linear relationship
text(x=2, y=1E-08, labels=expression(r^{1}), srt = 0, cex=2)

axis( 1, at=c(1, 1e1, 1e2, 1e3) )
axis( 2, at=c(1e-12, 1e-11, 1e-10, 1e-9, 1e-8, 1e-7, 1e-6) )
axis( 1, at=c(2:9/1, 2:9/0.1, 2:9/0.01, 2:9/0.001), labels=FALSE, tcl=-.2) # tcl=-0.2 for minor ticks
axis( 2, at=c(2:9/1e12, 2:9/1e11, 2:9/1e10, 2:9/1e9, 2:9/1e8, 2:9/1e7), labels=FALSE, tcl=-.2)

# Median DMS concentration in the surface ocean and interquartile range (Lana et al 2011)
abline(h=1.9E-09, lty=4)
rect(0.1,1.14E-09,2000,2.85E-09, col = rgb(0.5,0.5,0.5,0.1), lty=0)

#dev.off()

#-------------#
# References  #  
#-------------#

# Karp-Boss L, Boss E, Jumars PA (1996). Nutrient fluxes to planktonic osmotrophs in the presence of fluid motion. Oceanography and Marine Biology : an Annual Review 34: 71-107.
# Lana A, Bell TG, Simó R, Vallina SM, Ballabrera-Poy J, Kettle AJ, Dachs J, Bopp L, Saltzman ES, Stefels J, Johnson JE, Liss PS (2011). An updated climatology of surface dimethlysulfide concentrations and emission fluxes in the global ocean. Global Biogeochem Cycles 25: GB1004
# Saltzman ES, King DB, Holmen K, Leck C (1993). Experimental determination of the diffusion coefficient of dimethylsulfide in water. J Geophys Res, C 98: 16481-16486.
# Spiese CE, Le T, Zimmer RL, Kieber DJ (2016). Dimethylsulfide membrane permeability, cellular concentrations and implications for physiological functions in marine algae. J Plankton Res 38: 41-54.
# Sunda WG, Hardison R, Kiene RP, Bucciarelli E, Harada H (2007). The effect of nitrogen limitation on cellular DMSP and DMS release in marine phytoplankton: climate feedback implications. Aquat. Sci. 69:341-51.

#------------------------------------------------------------------------------------------#
# Fig.4 Modeling cell surface DMS enrichment in High DMS producers under mechanic stress   #
#                                                                                          #
# Lavoie M. 2018                                                                           #
#------------------------------------------------------------------------------------------#

#---------------------------------------#
# Scenaio 4 a) Assuming J scales with r #
#---------------------------------------#

#------------------#
# Model parameters #
#------------------#

J1 <- 3.4E-13                    # DMS flux for Heterocapsa (mol DMS cm-2 s-1) (Niki et al 2000)
r1 <- 7.7E-04                    # cell radius of Heterocapsa (cm) (Niki et al 2000)
Area <- 4*pi*(r1^2)              # Area in cm^2
Vol <- (4/3)*pi*((10000*r1)^3)   # Volume in um^3
Vol_L <- 1E-15*Vol               # Volume in Liter of biovolume
R <- J1 * (Area/Vol_L)           # production rate in mol DMS Lcell-1 s-1

D <- 1.19E-05         # DMS diffusion coefficient (cm^2 /s) (Saltzman et al 1993)
E <- c( 0.1, 1, 10 )  # Shear rate (s-1)
v <- 1E-02            # Seawater kinematic viscosity (cm^2 s-1) (Karp-Boss et al 1996)

#-----------#
# Variables #
#-----------#

x <- seq(2E-04, 5E-02, length=1000)  # cell radii (cm)
J <- R/1e3 * (x/3)                  # DMS flux (mol cm-2 s-1)  We assume that J scales with r

#-----------------#
# Model Functions #
#-----------------#

# Reynolds number (no units) (Karp-Boss et al 1996)
Re  <- outer( x^2, E/v )

# Peclet number (no units) (Karp-Boss et al 1996)
Pe  <- outer( x^2, E/D)

# Sherwood number (no units) (Karp-Boss et al 1996)
Sh <- Re * 0

i     <- Re<1e-1 & Pe<1e-2
j     <- Re<3e-1 & Pe>=1e-2 & Pe<1e2
k     <- Re<1    & Pe>=1e2
l     <- Re>=1

Sh[i] <- 1     + 0.29 * Pe[i]^(1/2)
Sh[j] <- 1.014 + 0.15 * Pe[j]^(1/2)
Sh[k] <-         0.55 * Pe[k]^(1/3)
Sh[l] <- NA

# Boundary layer thickness (cm) (Spiese et al 2016)
delta  <- x / Sh

# DMS cell surface enrichment (mol L^-1)
y <- (1e3 * J * x / D) * (1 - (x / (x + delta))) 

#------------------#
# Plotting output  #
#------------------#

x     <- x  * 1e4 # Convert cell radius in um

ri    <- which.min( abs( x - r1*1e4 ) )
r1_um <- x[ri]    # Radius closest to r1

#--- Graph log-log

#tiff("Fig4a_HighDMS_Stressed.tiff")                                      

par( mar = c(5, 5, 2, 2) )

# Solid, dashed and dotted lines correspond to the three shear rates from vector E, respectively.

plot(  x, y[,1], xlim=c(1, 1000), ylim=c(1e-9, 1e-4), xlab=expression(paste("cell radius ", mu, "m")), ylab=expression("DMS enrichment " ~  (mol ~ L^ ~ {-1})), log="xy", type="l", xaxp=c(1, 1000, 1), xaxt="n", yaxp=c(1e-12, 1e-6, 1), yaxt="n")
lines( x, y[,2], lty=2 )
lines( x, y[,3], lty=3 )
points( r1_um, y[ri,1], pch=19, cex=1)
axis( 1, at=c(1, 1e1, 1e2, 1e3) )
axis( 2, at=c(1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4) )
axis( 1, at=c(2:9/1, 2:9/0.1, 2:9/0.01, 2:9/0.001), labels=FALSE, tcl=-.2) # tcl=-0.2 for minor ticks
axis( 2, at=c(2:9/1e9, 2:9/1e8, 2:9/1e7, 2:9/1e6, 2:9/1e5, 2:9/1e4), labels=FALSE, tcl=-.2)

# Median DMS concentration in the surface ocean and interquartile range (Lana et al 2011)
abline(h=1.9E-09, lty=4)
rect(0.1,1.14E-09,2000,2.85E-09, col = rgb(0.5,0.5,0.5,0.1), lty=0)

#dev.off()

#--------------------------------------------------------#
# Scenario 4b) Assuming a constant J for all cell radius #
#--------------------------------------------------------#

#------------------#
# Model parameters #
#------------------#

J1 <- 3.4E-13     # flux DMS pour Heterocapsa (Niki et al 2000)
r1 <- 7.7E-04        # cell radius of Heterocapsa (Niki et al 2000)

D <- 1.19E-05         # DMS diffusion coefficient (cm^2 /s) (Saltzman et al 1993)
E <- c( 0.1, 1, 10 )  # Shear rate (s-1)
v <- 1E-02           # Seawater kinematic viscosity (cm^2 s-1) (Karp-Boss et al 1996)

#-----------#
# Variables #
#-----------#

x <- seq(2E-04, 5E-02, length=1000)  # cell radii (cm)
J <- J1                # DMS flux (mol cm-2 s-1)  We assume that J is constant

#-----------------#
# Model Functions #
#-----------------#

# Reynolds number (no units) (Karp-Boss et al 1996)
Re  <- outer( x^2, E/v )

# Peclet number (no units) (Karp-Boss et al 1996)
Pe  <- outer( x^2, E/D)

# Sherwood number (no units) (Karp-Boss et al 1996)
Sh <- Re * 0

i     <- Re<1e-1 & Pe<1e-2
j     <- Re<3e-1 & Pe>=1e-2 & Pe<1e2
k     <- Re<1    & Pe>=1e2
l     <- Re>=1

Sh[i] <- 1     + 0.29 * Pe[i]^(1/2)
Sh[j] <- 1.014 + 0.15 * Pe[j]^(1/2)
Sh[k] <-         0.55 * Pe[k]^(1/3)
Sh[l] <- NA

# Boundary layer thickness (cm) (Spiese et al 2016)
delta  <- x / Sh

# DMS cell surface enrichment (mol L^-1)
y <- (1e3 * J * x / D) * (1 - (x / (x + delta)))

#------------------#
# Plotting output  #
#------------------#

x     <- x  * 1e4 # Convert cell radius in um

ri    <- which.min( abs( x - r1*1e4 ) )
r1_um <- x[ri]    # Radius closest to r1

#--- Graph log-log

#tiff("Fig4b_cstJ_HighDMS_stressed.tiff")                                        

par( mar = c(5, 5, 2, 2) )

# Solid, dashed and dotted lines correspond to the three shear rates from vector E, respectively.

plot(  x, y[,1], xlim=c(1, 1000), ylim=c(1e-9, 1e-5), xlab=expression(paste("cell radius ", mu, "m")), ylab=expression("DMS enrichment " ~  (mol ~ L^ ~ {-1})), log="xy", type="l", xaxp=c(1, 1000, 1), xaxt="n", yaxp=c(1e-12, 1e-6, 1), yaxt="n")
lines( x, y[,2], lty=2 )
lines( x, y[,3], lty=3 )
points( r1_um, y[ri,1], pch=19, cex=1)
axis( 1, at=c(1, 1e1, 1e2, 1e3) )
axis( 2, at=c(1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4) )
axis( 1, at=c(2:9/1, 2:9/0.1, 2:9/0.01, 2:9/0.001), labels=FALSE, tcl=-.2) # tcl=-0.2 for minor ticks
axis( 2, at=c(2:9/1e9, 2:9/1e8, 2:9/1e7, 2:9/1e6, 2:9/1e5, 2:9/1e4), labels=FALSE, tcl=-.2)

# Median DMS concentration in the surface ocean and interquartile range (Lana et al 2011)
abline(h=1.9E-09, lty=4)
rect(0.1,1.14E-09,2000,2.85E-09, col = rgb(0.5,0.5,0.5,0.1), lty=0)

#dev.off()

#-------------#
# References  #  
#-------------#

# Niki T, Kunugi M, Otsuki A (2000). DMSP-lyase activity in five marine phytoplankton species: it's potential importance in DMS production. Mar Biol 136: 759-764.
# Karp-Boss L, Boss E, Jumars PA (1996). Nutrient fluxes to planktonic osmotrophs in the presence of fluid motion. Oceanography and Marine Biology : an Annual Review 34: 71-107.
# Lana A, Bell TG, Simó R, Vallina SM, Ballabrera-Poy J, Kettle AJ, Dachs J, Bopp L, Saltzman ES, Stefels J, Johnson JE, Liss PS (2011). An updated climatology of surface dimethlysulfide concentrations and emission fluxes in the global ocean. Global Biogeochem Cycles 25: GB1004
# Saltzman ES, King DB, Holmen K, Leck C (1993). Experimental determination of the diffusion coefficient of dimethylsulfide in water. J Geophys Res, C 98: 16481-16486.
# Spiese CE, Le T, Zimmer RL, Kieber DJ (2016). Dimethylsulfide membrane permeability, cellular concentrations and implications for physiological functions in marine algae. J Plankton Res 38: 41-54.
