subroutine RUNFIT(npts, npar, frq, signal, plco, plar, tl, dl, tpl, LL, params0, f_sd, params1, dparams1, resid1, JAC_FLAG, a_fit_log, a_rms_log, a_params, a_out)

implicit integer (i-k)
implicit logical (f) 
implicit character*32 (a)

interface

function DIAG(A)
double precision, DIMENSION(:,:), INTENT(IN) :: A
double precision, DIMENSION(size(A,1),size(A,2)) :: DIAG
end function DIAG

function IDENT(n)
integer n
double precision IDENT(n,n)
end function IDENT

end interface

integer, intent(in):: npts, npar, f_sd
!npts - number of total points
!npar - number of total parameters

!f_sd is flag: 1 -> g2!=0, 0 -> g2==0

double precision, intent(in):: frq(npts), signal(npts), plco, plar, tl, dl, tpl, LL, params0(npar)
double precision, intent(in):: JAC_FLAG(npar)
! frq - frequencies 
! signal - absorption values
! plco - CO pressure
! plar - Ar pressure
! tl - temperature
! dl - deviation
! tpl - record type
! LL - cell length
! params0 - initial parameters values

character*255, intent(in)::a_fit_log, a_rms_log, a_params, a_out

double precision, intent(out):: resid1(npts), params1(npar), dparams1(npar)

double precision ap, rms1, rms2, model1(npts), params(npar), steparams(npar), JAC(npts, npar), resid(npts)
double precision AM1(npar,npar), AM2(npar,npar), AV(npar)
double precision ps(npar), dp(npar), rel_thres

!double precision  k_st, COstr, pi, frabi, ngamma, gdop, NVC
!integer nconstpar !number of parameters common for all recordings (ones we need to obtain)
!COMMON /CO_const/ k_st, COstr, pi, frabi, nconstpar, ngamma, gdop, NVC !constants set in the main.f90 unit

double precision  k_st, COstr, pi, frabi, ngamma, gdop, Lcell, S_Tdep, NVC
integer nconstpar !number of parameters common for all recordings (ones we need to obtain)
COMMON /CO_const/ k_st, COstr, pi, frabi, nconstpar, ngamma, gdop, Lcell, S_Tdep, NVC !constants set in the main.f90 unit

params = params0

rel_thres = 1.D-4

!params(6) = 1.00023

if (f_sd.eq.0) then
	params(3) = 1.D-20  ! accounting for no sd-effect	
	params(4) = 0.d0  ! means that g2 = d2 = 0 
	!params(6) = 0.995066 ! SD and noSD intensities are different
	params(13) = 0.D0
end if


open(10,file=a_fit_log,status='UNKNOWN') !,position='APPEND')
open(11,file=a_rms_log,status='UNKNOWN') !,position='APPEND')
1100 format(I5,I5, D15.9)
1101 format(11(D15.9,2X))	 !number in front of D15.9 should be NCONSTPAR

write(*,*)'fit start, initial residual calc'
call MDL(npts, npar, frq, signal, model1, resid, plco, tl, dl, tpl, LL, params)
write(*,*)'fit start, initial residual calc done'

rms1 = dsqrt(sum(resid(:)**2))
write(*,*)'fit start, initial RMS done'

WRITE(*,*) 'starting RMS: ', rms1
if (f_sd.eq.1) WRITE(*,*) '!!!!START SD PROFILE FIT!!!!'
if (f_sd.eq.0) WRITE(*,*) '!!!!START CLASSIC FIT!!!!'
call time(a_time)
WRITE(10,*) a_time, ' !!!!START!!!!'
WRITE(11,1100) -1,0,rms1
write(11,1101) params0(1:11)

k=0
f_end=.false.

do while (.not.(f_end))!exit condition is gradient changing slowly or fit is too long
	ap=1d-20 !ap is lambda in Levenverg-Marquardt
	i=0
	call JACOBI(npts, npar, frq, signal, plco, tl, dl, tpl, LL, params,JAC,JAC_FLAG)
	if (f_sd.eq.0) then
		JAC(:,3)  = 0.D0
		JAC(:,4)  = 0.D0
		JAC(:,13) = 0.D0		
	end if
	if (f_sd.eq.1) write(*,*)'SD FIT step = ',k, ',   substep = ', i, '  jacobi done!'
	if (f_sd.eq.0) write(*,*)'CLASSIC FIT step = ',k, ',   substep = ', i, '  jacobi done!'
	f_p=.false.
	do while (.not.f_p)

		AM1=ap*ident(npar)+matmul(transpose(JAC),JAC)
		resid(:)=(signal(:)-model1(:))
		AV=matmul(transpose(JAC),resid)
		CALL DLINDS(npar, AM1, npar, AM2, npar)
		steparams=matmul(AM2,AV)
		call time(a_time)
		write(*,*)a_time,'step = ',k,'   substep = ',i
		write(10,*)a_time,'step = ',k,'   substep = ',i
		write(*,*) 'f0 step = ', steparams(1), ' MHz'
		write(*,*) 'calculating the new residual'
		call MDL(npts, npar, frq, signal, model1, resid, plco, tl, dl, tpl, LL, params+steparams)
		resid(:)=resid(:)
		rms2 = dsqrt(sum(resid(:)**2))
		write(*,*) 'old rms = ', rms1, '   new rms = ', rms2, '   lambda = ', ap
		write(*,*) 'relative difference = ', (rms2-rms1)/rms1
		write(10,*) 'old rms = ', rms1, '   new rms = ', rms2, '   lambda = ', ap
		WRITE(11,1100) k,i,rms2
		WRITE(11,*) 'Dopler width = ', gdop
		WRITE(11,1101) params(1:11)

	    if (rms2.le.rms1)then                                    !if next-step residual less than current (and the relative 
	        if (abs((rms2-rms1)/rms1).gt.(rel_thres)) then            !difference is greater than 0.0005), we make a step 
				write(*,*) 'next rms is smaller, STEP FORWARD'
				write(10,*) 'next rms is smaller, STEP FORWARD'
				f_p=.true.                                      !in the direction defined above
				params=params+steparams                               !
				rms1=rms2                                        ! 
			else if (abs((rms2-rms1)/rms1).le. (rel_thres)) then      !if the difference between current and next residuals is too small,
				write(*,*) 'rms changes too slow, STOP!'
				write(10,*) 'rms changes too slow, STOP!'
				f_p=.true.                                      !we force the algorithm to finish
				f_end=.true.                                    !
			end if
		else if (rms2.gt.rms1) then !if next-step residual greater than current, we "overstepped" minimal rms
			ap=ap*1.D4             !so we increase lambda, which makes the step less
			write(*,*) 'next rms is greater, change LAMBDA and repeat'                       
			write(10,*) 'next rms is greater, change LAMBDA and repeat'                       
			if (i.eq.10) then
				f_end=.true.
				f_p=.true.
				write(*,*)'can"t converge to a solution, BREAK' 
				write(10,*)'can"t converge to a solution, BREAK' 
			end if
		end if
		i=i+1
		
	end do
	k=k+1
end do

close(10)
close(11)

!!!!!!!!!!!!CALCULATION OF THE UNCERTAINTIES!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call MDL(npts, npar, frq, signal, model1, resid, plco, tl, dl, tpl, LL, params)

rms2 = dsqrt(sum(resid(:)**2)/(npts-1)) !stdev of the residual (hope to be noise-like to treat it as noise sigma)

AM1 = matmul(transpose(JAC),JAC) ! evaluation according to Rosenkranz, step 1
!write(*,*) 'jacobi*jacobi'
! the following is protection against zero diagonal values in AM1 
! it may happen if k-th parameter is not applicable to fitted recordings, in this case we will have zero everywhere
! for this partial parameter, and zero diag.element makes matrix inversion impossible
do k=1,npar
	if (AM1(k,k).eq.0.D0) AM1(k,k)=1.D-6  
end do
! end of protection
!write(*,*) 'protection'


CALL DLINDS(npar, AM1, npar, AM2, npar) ! evaluation step 2 (see the origins in the memo related to the MPM/ECS revision paper!)

!write(*,*) 'inversion'

do ipar=1,npar
	dp(ipar)=rms2*dsqrt(AM2(ipar,ipar))
end do

ps = params

write(*,*) 'uncertainty done'
! stat.uncert. of the i-th parameter is sqrt of diag.element of the inverted matmul(transpose(JAC),JAC)
! multiplied to the residual's sigma
! matrix AM2 multiplied to rms2**2 is the covariance matrix
! normalized to sqrt(AM2(i,i)*AM2(j,j)) will be the correlation matrix

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do ipar = 1,npar
	do jpar = ipar+1,npar
		AM2(ipar,jpar) = AM2(ipar,jpar)/dsqrt(AM2(ipar,ipar)*AM2(jpar,jpar))
		AM2(jpar,ipar) = AM2(ipar,jpar)
	end do
	AM2(ipar,ipar)=1.D0
end do

write(*,*) 'correlation done'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
open(12,file=a_params, status = 'UNKNOWN', POSITION = 'APPEND')

write(*,*) a_params

1212 format(23(F15.6,2X))
1313 format(11(F15.6,2X))
!						       frq	  f_err	  G0	G0_err  G2    G2_err  D2    D2_err Y     Y_err   Scale  Scale_err	pow.f.  BL0   BL1    BL2      BL3  2nd harmnc	Nu_vc	Nu_err
WRITE(12,1212) plco, plar, tl, ps(1), dp(1), ps(2), dp(2), ps(3), dp(3), ps(4), dp(4), ps(5), dp(5), ps(6), dp(6),      ps(7), ps(8), ps(9), ps(10), ps(11), ps(12),    ps(13), dp(13)
write(*,*) 'params written'
close(12)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

write(*,*) a_out

!open(13, file = a_out, status = 'UNKNOWN')

!do i_str=1, npar
!   write(13,1313) AM2(i_str,:)
!end do

!close(13) 

!write(*,*) 'correlation written'

params1 = params
dparams1 = dp
!write(*,*) '!!'

resid1 = resid

!write(*,*) '!!!'

write(*,*) 'return from fit'
!write(*,*) 'PRESS ENTER'
!read(*,*)
return
end subroutine RUNFIT



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine init_params(npts, npar, frq, signal, plco, plar, tl, dl, tpl, LL, params0)

integer, intent(in):: tpl, npts, npar
!npts - number of total points
!npar - number of total parameters

double precision, intent(in):: frq(npts), signal(npts), plco, plar, tl, dl, LL 
double precision, intent(out):: params0(npar)

double precision tmpS, tmpScl, tmpG0, tmpG2, tmpctr, tmpBL0, tmpBL1 !, tmpBL2, tmpPow
double precision tmpxmax, tmpxmin, tmpdelta, tmpamp, tmpdop, tmpIm

double precision  k_st, COstr, pi, frabi, ngamma, gdop, Lcell, S_Tdep, NVC
integer nconstpar !number of parameters common for all recordings (ones we need to obtain)
COMMON /CO_const/ k_st, COstr, pi, frabi, nconstpar, ngamma, gdop, Lcell, S_Tdep, NVC !constants set in the main.f90 unit

!write(*,*) 'type ', tpl
params0 = 0.D0

tmpS = COstr*K_st*plco*LL/tl
write(*,*) tmpS/LL
!read(*,*)
tmpDop = gdop*dsqrt(tl)

tmpxmax = frq(sum(maxloc(signal))) ! location of the maximum signal

tmpBL0 = 0.D0

do ibl = 1,10
	tmpBL0 = tmpBL0+(signal(ibl)+signal(npts+1-ibl))	
end do

tmpBL0 = 5.D-2*tmpBL0
tmpBL1 = (signal(npts)-signal(1))/(frq(npts)-frq(1))



if ((tpl.eq.1).or.(tpl.eq.2)) then !for the differential recording with freq. deviation

	tmpxmin = frq(sum(minloc(signal)))
	tmpctr = 0.5D0*(tmpxmax+tmpxmin)
	!tmpctr = 115271.18
	!write(*,*) 'center ', tmpctr
	tmpdelta = abs(tmpxmax - tmpxmin)
	!write(*,*) 'min-max distance is ', tmpdelta
	tmpamp = 0.5D0*abs(maxval(signal)-minval(signal))
	!write(*,*) 'amplitude ', tmpamp
	!tmpG0 = 0.25D0*dsqrt(3.D0)*tmpdelta+0.5D0*dsqrt(3.D0*abs(0.25D0*tmpdelta**2-dl**2))	
	tmpG0 = (0.16041/0.045)*plco + 2.4698D0*plar
	!tmpG0 = dsqrt(tmpG0**2+0.093D0**2)
	!write(*,*) 'HWHM ', tmpG0
	tmpG2 = 0.1*tmpG0
	!tmpG2 = 0.343*plco + 0.322*plar !0.1*tmpG0
	!tmpS = tmpS*LL
	!write(*,*) 'line intensity ', tmpS
	!write(*,*) 'test ', difSDVe(0.D0,tmpctr,tmpdop,tmpG0,0.D0,0.D0,0.D0,tmpS,dl,1.D0,0.D0)	
!                   difSDVe(sg,sg0,GamD,Gam0,Gam2,Shift0,Shift2,strength_in,dev,scale,pow)	
	tmpScl = 1.0D0*tmpamp/difSDVe(tmpxmax,tmpctr,tmpdop,tmpG0,tmpG2,0.D0,0.D0,tmpS,dl,1.D0,0.D0,tpl)
	!difSDV(tmpxmax,tmpctr,tmpdop,tmpG0,tmpG2,0.D0,0.D0,tmpS,dl)
	!difSDV(sg,sg0,GamD,Gam0,Gam2,Shift0,Shift2,strength_in,dev)
	write(*,*) 'scale ', tmpScl
	!read(*,*)
end if

if (tpl.eq.3) then
tmpamp = maxval(signal)-tmpBL0
tmpctr = tmpxmax
tmpG0 = 3.4D0*plco + 2.5D0*plar
tmpG2 = 0.1*tmpG0
call HTP(tmpctr,tmpDop,tmpG0,tmpG2,0.D0,0.D0,NVC, 0.D0,tmpctr,tmpScl,tmpIm)
!tmpScl = 1-dexp(-tmpScl*tmpS)
tmpScl = 1.D0*tmpamp/tmpScl !SDRP(tmpctr, tmpctr,0.D0,tmpG0,tmpG2,0.D0,0.D0,tmpS)
end if

if (tpl.eq.0) then	!for direct profile observation as is
	tmpS = tmpS*1.D5  !1/cm to 1/km
	tmpctr = tmpxmax
	tmpamp = maxval(signal)-tmpBL0
	imax = sum(maxloc(signal))
	do while (abs(signal(imax)-maxval(signal)).le.0.5D0*tmpamp)
		imax = imax+1
	end do
	!tmpG0 = abs(frq(imax)-tmpxmax)
	tmpG0 = 3.5D0*plco + 2.8D0*plar
	tmpG2 = 0.1*tmpG0	
	tmpScl = 1.D0 !*tmpamp/SDRP(tmpctr, tmpctr,0.D0,tmpG0,tmpG2,0.D0,0.D0,tmpS)
	!tmpScl = 1.D0
	tmpBL1 = 0.D0
end if

params0(1) = tmpctr	!center
params0(2) = tmpG0  !*0.5	! G or G0
params0(3) = tmpG2  ! G2 for SD profile
params0(4) = 0.D0	!D2	for SD profile
params0(5) = 0.D0	!Y	mixing
params0(6) = tmpScl*1.0	!scaling intensity parameter
params0(7) = 0.D0	!power correction factor for RAD and VID recordings
params0(8) = tmpBL0	!constant baseline
params0(9) = tmpBL1	!linear baseline
if (tpl.eq.0) params0(9) = 0.D0
params0(10)= 0.D0	!square baseline
params0(11)= 0.D0	!cubic baseline

end subroutine init_params