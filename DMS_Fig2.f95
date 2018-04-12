!---------------------------------------------------------------------------!
! Modeling cell surface DMS enrichment                                       !
! This is the Fortran 95 code for the results of Fig 2                      !
!                                                                           !
! Lavoie M. 2018                                                            !
!---------------------------------------------------------------------------!

program DMS_Fig2

implicit none
integer, parameter :: num = 5000, inc = 1
integer, parameter :: ikind = selected_real_kind(p=18)
real (kind = ikind), dimension (num) :: x2
real (kind = ikind), dimension (num, 3) :: E
real (kind = ikind) :: J1, R1, pi, area, vol, vol_l, start, finish, i2, step
real (kind = ikind), dimension (num, 3) :: R, Diff, visc, x, DMSflux, Re, Pe, Sh, delta, Pext, y
integer :: j, i         ! i and j are control variables or counters

call cpu_time(start)

!------------------!
! Model parameters !
!------------------!
J1 = 2.8E-16                                    ! DMS flux for I galbana (mol DMS cm-2 s-1)
r1 = 2.4E-4                                     ! Cell radius of I galbana (cm)
pi = 4. * atan(1.0)
area = 4. * pi * r1**2.0                        ! Cell surface area in cm^2

vol = (4.0/3.0) * pi * ((10000 * r1)**3)        ! Cell volume in um^3
vol_l = 1E-15 * vol                             ! Cell volume in Liter of biovolume

print*, area, vol, vol_l

do j = 1, num, inc
  	do i = 1, 3
	E(j, i) = (10.**i) / 100.            ! Fill the R array of shear rates (s-1)
	R(j, i) = J1 * (area/vol_l)          ! Fill the R array of DMS production rate in mol DMS Lcell-1 s-1. R is repeated for each x for 3 E (so there are 300 values)
	Diff(j, i) = 1.19E-05                ! Fill the Diff array diffusion coefficients (cm^2 s-1). Diff. is repeated for each x and for each of the three E
	visc(j, i) = 1.0E-2                  ! Seawater kinematic viscosity (cm^2 s-1) is repeated for each x and for each of the three E
	end do
end do
! Print the matrix of each parameter
print*
print*, "Here is the E matrix", E
print*
print*, "Here is the R matrix", R
print*
print*, "Here is the Diff matrix", Diff
print*
print*, "Here is the visc matrix", visc

!-----------------!
! Model variables !
!-----------------!

i2 = 0.95
step =  0.05
do j = 1, num, inc
i2 = step + i2
x2(i2) = i2 / 5000          ! Create a numeric vector inside a loop since the fortran loop (j = 1, num, inc) requires integers.
  	do i = 1, 3
	x(j, i) = x2(i2)      ! cell radii (cm) from 2E-04 to 5E-02 cm or 2 um to 500 um and repeated at each of the three E.
	end do
end do
print*
print*, "Here is the x matrix", x 

DMSflux = R/1000. * (x/3.)                          ! DMS flux (mol cm-2 s-1) repeated at all x and for three E ; (1/1000) is converting L of cell volume in cm^3
print*
print*, "Here is the matrix DMSflux", DMSflux

!------------------------------------------------!
! Model equations
!------------------------------------------------!

! Calculate Reynolds and Peclet numbers
Re = (x*x * E / visc)     ! Reynold number (no units)
Pe = (x*x * E / Diff)     ! Peclet number (no units)
print*
print*, "matrix Re :", Re

! Calculating the Sherwood number (Sh, no units) at three shear rate (E)
Sh = 0             

do i = 1, num, inc
  do j = 1, 3 
  if (Re(i,j) < 0.1 .and. Pe(i,j) < 0.01) then
    Sh(i,j) = 1 + 0.29 * Pe(i,j)**0.5
  else if (Re(i,j) < 0.3 .and. Pe(i,j) < 100. .and. Pe(i,j) >= 0.01) then
    Sh(i,j) = 1.014 + 0.15 * Pe(i,j)**0.5
  else if (Re(i,j) < 1. .and. Pe(i,j) >= 100.) then
    Sh(i,j) = 0.55 * Pe(i,j)**(1./3.)
  else if (Re(i,j) >= 1.) then
    Sh(i,j) = 1.0E+15      ! Inputtiny a very large number, so that the final enrichments at Re >= 1 are near 0.
  end if
  end do
end do

! Calculating the cell surface DMS enrichment
delta = x / Sh                 ! Boundary layer thickness (cm) at the three E
y = (1000. * DMSflux * x / Diff) * (1. - (x / (x + delta)))          ! DMS cell surface enrichmeent (mol/L) at the three E

! printing the final cell surface DMS enrichment at different E
print*
print*, "DMS cell surface enrichment", y
print*
print*, "DMS cell surface enrichment at E = 0.1 s-1 : "
print*, y(:,1)
print*
print*, "DMS cell surface enrichment at E = 1 s-1 : "
print*, y(:,2)
print*
print*, "DMS cell surface enrichment at E = 10 s-1 : "
print*, y(:,3)

print*, x(1,1), x(1,2), x(1,3)
print*, x(5000,1)

! Printing in a text file the cell surface DMS enrichment at different shear rates
open(unit=10, file = "DMS_enrich.txt", status = "unknown")
write(10, *) "cell radius =               ", x(:,1), "Enrichment at 0.1 s-1  =        ", &
y(:,1), "Enrichment at 1 s-1  =       ", y(:,2), "Enrichment at 10 s-1 =          ", y(:,3)
close(10)

!Calculate the CPU execution time
call cpu_time(finish)
print*
print*, finish
print*, start
print*, "Time = ", finish-start, "seconds"

end program
