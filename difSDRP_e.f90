double precision function difSDRP(frq_in, f0_in,Y_in,gamma0_in,gamma2_in,delta0_in,delta2_in,strength_in,dev_in)

use MSFLIB
use MSIMSL
USE PORTLIB

interface

double precision function FMB(vel)
double precision vel
end function FMB

double precision function difSDRP_UNDERINT(vel)
double precision vel
end function difSDRP_UNDERINT

end interface

double precision frq_in, f0_in,Y_in,gamma0_in,gamma2_in,delta0_in,delta2_in,strength_in,dev_in

double precision frq, f0,Y,gamma0,gamma2,delta0,delta2,strength,dev
common           frq, f0,Y,gamma0,gamma2,delta0,delta2,strength,dev

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
dev = dev_in

CALL DQDAGI(difSDRP_UNDERINT, 0.D0, 1, 1.D-12, 1.D-12, absorp, errest)

difSDRP = strength*absorp*(frq/f0)**2/pi

end function difSDRP



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
double precision function difSDRP_UNDERINT(vel)!,f0,gamma0,gamma2,delta0,delta2,y,strength,tmpr,mass)
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

double precision frq, f0,Y,gamma0,gamma2,delta0,delta2,strength,dev
common           frq, f0,Y,gamma0,gamma2,delta0,delta2,strength,dev

double precision gamma_tot, delta_tot, dif_plus, dif_minus

gamma_tot = PAR_VS_V(gamma0,gamma2,vel)
delta_tot = PAR_VS_V(delta0,delta2,vel)
dif_plus = frq + dev - f0 - delta_tot
dif_minus = frq - dev - f0 - delta_tot

difSDRP_UNDERINT = (gamma_tot+Y*dif_plus)/(gamma_tot**2+dif_plus**2)
difSDRP_UNDERINT =	difSDRP_UNDERINT - (gamma_tot+Y*dif_minus)/(gamma_tot**2+dif_minus**2)
difSDRP_UNDERINT =	difSDRP_UNDERINT*FMB(vel)

end function difSDRP_UNDERINT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!