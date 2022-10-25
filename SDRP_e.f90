double precision function SDRP (frq_in, f0_in,Y_in,gamma0_in,gamma2_in,delta0_in,delta2_in,strength_in)

use MSFLIB
use MSIMSL
USE PORTLIB

interface

double precision function FMB(vel)
double precision vel
end function FMB

double precision function SDRP_UNDERINT(vel)
double precision vel
end function SDRP_UNDERINT

end interface

double precision frq_in, f0_in,Y_in,gamma0_in,gamma2_in,delta0_in,delta2_in,strength_in

double precision frq, f0,Y,gamma0,gamma2,delta0,delta2,strength
common           frq, f0,Y,gamma0,gamma2,delta0,delta2,strength

double precision kB, c, h, m_aem, pi, k_vs_aem
common /CB1/ kB, c, h, m_aem, pi, k_vs_aem

DATA kB /1.38064852D-23/, c /2.997925D08/, h /6.62607004D-34/, m_aem /1.6605402D-27/
DATA pi /3.14159265359D0/, k_vs_aem /0.120272478907D0/
!all units for SI system except k_vs_aem which is m_aem in grams divided by kB!

double precision absorp, errest

frq=frq_in
f0=f0_in
Y=Y_in
gamma0=gamma0_in
delta0=delta0_in
gamma2=gamma2_in
delta2=delta2_in
strength=strength_in

CALL DQDAGI(SDRP_UNDERINT, 0.D0, 1, 1.D-12, 1.D-12, absorp, errest)

SDRP = strength*absorp*(frq/f0)**2/pi

end function SDRP

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
double precision function FMB(vel)
!Maxwell-Boltzmann distribution
!for the velocity absolute value(?)
!mass - molar mass in g/mol
!tmpr - T in K
!vel - in m/s (???)
double precision vel
!double precision M1,vel2
double precision vel2

double precision kB, c, h, m_aem, pi, k_vs_aem
common /CB1/kB, c, h, m_aem, pi, k_vs_aem

!M1 = mass*k_vs_aem/tmpr !optimization, eh?
vel2 = vel*vel			!

!FMB = 0.797884560803*vel2*M1*DSQRT(M1)
!FMB = FMB*DEXP(-vel2*0.5*M1)

FMB = 2.256758D0*vel2*DEXP(-vel2)

end function FMB
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
double precision function PAR_VS_V(par0, par2, vel)
!v-dependence for gamma and delta (see paper)
!par0 == gamma0 or delta0 (in freq. units)
!par2 == gamma2 or delta2 (in freq. units)
!vel - velocity (in the same units as everywhere
!tmpt - T in Kelvin
!mass - molar mass in g/mol
double precision par0, par2, vel
!double precision vel_prob !most probable velocity
double precision kB, c, h, m_aem, pi, k_vs_aem
common /CB1/kB, c, h, m_aem, pi, k_vs_aem

!vel_prob=DSQRT(2*tmpr/(mass*k_vs_aem))
!PAR_VS_V=par0+par2*((vel/vel_prob)**2-1.5D0)
PAR_VS_V=par0+par2*(vel**2-1.5D0)

end function PAR_VS_V
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
double precision function SDRP_UNDERINT(vel)!,f0,gamma0,gamma2,delta0,delta2,y,strength,tmpr,mass)
!function which is integrated to obtain SDRP profile
!SDRP is Speed Dep-t Rosenkranz Profile ()
!f0 - central frequency, gamma0/2 - width (collisional and wind)
!delta0/2 - shift (collisional and wind), Y - 1-st order mixing
!strendth - integral intensity, tmpr - T in kelvin, mass - molar mass in g/mol
interface

double precision function FMB(vel)
double precision vel !,tmpr,mass
end function FMB

double precision function PAR_VS_V(par0, par2, vel)
double precision par0, par2, vel
end function PAR_VS_V

end interface

double precision vel

double precision kB, c, h, m_aem, pi, k_vs_aem
common      /CB1/kB, c, h, m_aem, pi, k_vs_aem

double precision frq, f0,Y,gamma0,gamma2,delta0,delta2,strength
common           frq, f0,Y,gamma0,gamma2,delta0,delta2,strength

double precision gamma_tot, delta_tot, dif_plus, dif_minus

gamma_tot = PAR_VS_V(gamma0,gamma2,vel)
delta_tot = PAR_VS_V(delta0,delta2,vel)
dif_plus = frq - f0 - delta_tot
dif_minus = frq + f0 + delta_tot

SDRP_UNDERINT = (gamma_tot+Y*dif_plus)/(gamma_tot**2+dif_plus**2)
SDRP_UNDERINT =	SDRP_UNDERINT + (gamma_tot-Y*dif_minus)/(gamma_tot**2+dif_minus**2)
SDRP_UNDERINT =	SDRP_UNDERINT*FMB(vel)

end function SDRP_UNDERINT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!