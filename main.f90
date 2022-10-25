program CO_SINGLEFIT

!I added this comment for test to track changes

use MSFLIB
use MSIMSL

implicit integer (i-k)
implicit logical (f) 
implicit character*32 (a)

interface

double precision function SDRP (frq_in, f0_in,Y_in,gamma0_in,gamma2_in,delta0_in,delta2_in,strength_in)
double precision frq_in, f0_in,Y_in,gamma0_in,gamma2_in,delta0_in,delta2_in,strength_in
end function SDRP

double precision function difSDRP(frq_in, f0_in,Y_in,gamma0_in,gamma2_in,delta0_in,delta2_in,strength_in,dev)
double precision frq_in, f0_in,Y_in,gamma0_in,gamma2_in,delta0_in,delta2_in,strength_in,dev
end function difSDRP

double precision function difSDV(sg,sg0,GamD,Gam0,Gam2,Shift0,Shift2,strength_in,dev)
double precision sg,sg0,GamD,Gam0,Gam2,Shift0,Shift2,strength_in,dev
end function difSDV

double precision function difSDVe(sg,sg0,GamD,Gam0,Gam2,Shift0,Shift2,strength_in,dev, scale,pow,tpl)
double precision sg,sg0,GamD,Gam0,Gam2,Shift0,Shift2,strength_in,dev, scale,pow
integer tpl
end function difSDVe

double precision function difSDVe_add(sg,sg0,GamD,Gam0,Gam2,Shift0,Shift2,strength_in,dev, scale,pow,amp,NVC,tpl)
double precision sg,sg0,GamD,Gam0,Gam2,Shift0,Shift2,strength_in,dev, scale,pow,amp,NVC
integer tpl
end function difSDVe_add

function DIAG(A)
double precision, DIMENSION(:,:), INTENT(IN) :: A
double precision, DIMENSION(size(A,1),size(A,2)) :: DIAG
end function DIAG

function IDENT(n)
integer n
double precision IDENT(n,n)
end function IDENT

double precision function Q_PART(Tq,Tref,npart,tpart,qpart)
double precision Tq, Tref
integer npart
double precision tpart(npart),qpart(npart)
end function Q_PART


end interface

integer npts, npar, nfil,dirlength

!double precision k_st, tmpW, tmpxmax, tmpxctr, tmpdel,tmpA,tmpBL0,tmpBL1, COstr,pi, Lcell
!double precision tmpS, tmpDike !, tmpPow, tmpG0, tmpG2, tmpD0, tmpD2, tmpY,tmpCself, tmpCar
!double precision rms1, frabi, ngamma, gdop, k_vs_aem, clight, molm
!arrays for fit
double precision, allocatable:: frq(:), signal(:), resid(:), resid_sd(:) 
double precision  k_vs_aem, plco, plar, tl, dl, LL 
double precision, allocatable:: model1(:),model2(:), JAC_FLAG(:) !, JAC(:,:)
integer tpl 
!frq - list of f-s, signal - list of signals, plco - CO pressure for each point, plar - Ar pressure f.e.p.
!tl - temperature f.e.p.,
! dl - deviation (if needed), tpl - type (0 - resonator, profile as is, absolute absorption value (intensity not fitted),-->
!--> 1 - RAD, differential profile, relative absorption (intensity fitted), ind - record index

!double precision tmpS, tmpScl, tmpG0, tmpG2, tmpctr, tmpBL0, tmpBL1 !, tmpBL2, tmpPow
!double precision tmpxmax, tmpxmin, tmpdelta, tmpamp, tmpdop

double precision tpart(400),qpart(400)	!change size here if needed!!!
integer npart


!aux array used in the fitting routine
!double precision, allocatable:: AM1(:,:),AM2(:,:),AV(:) !, DER(:,:), DERINV(:,:)
!double complex, allocatable:: AM3(:,:)

!aux arrays									 
double precision, allocatable:: params0(:), params(:),params_sd(:), dp(:), dp_sd(:)
double precision, allocatable:: par_start_self(:), par_start_foreign(:)

double precision qual, qual_sd, molm
character($MAXPATH) a_data,a_data_f,a_data_out, a_params_sd, a_params_cl, a_corr_cl, a_corr_sd


character($MAXPATH) a_list

character*128 a_header, a_folder,a_folder_out, a_partition, a_aux

CHARACTER($MAXPATH) dir, a_script

double precision  k_st, COstr, pi, frabi, ngamma, gdop, Lcell, t_pl, Elow, c2, frq0, S_Tdep, NVC
integer nconstpar
COMMON /CO_const/ k_st, COstr, pi, frabi, nconstpar, ngamma, gdop, Lcell, S_Tdep, NVC !Elow, c2, frq0

!!!!!!!!constants!!!!!!!!!!!!!
k_st = 2894.947919D20 !recalc coef-t for the intensity
Lcell = 10.D0         !cm, cell length
frabi = 20.D-2
pi = 3.14159265359D0
clight = 2.997925D08 ! speed of light, m/s
k_vs_aem = 0.120272478907D-3
c2 = 1.4388 !cm*K c2 = hc/k_B for intensity calcs
!k_vs_aem = atomic mass unit [kg] divided by k_B [J/K]
NVC = 0.0D0
!write(*,*) clight, k_vs_aem

COstr = 3.3D-24       !cm/molec HITRAN intensity J=0-1/1-2 line  (3.3e-24/2.566e-23) 2.566d-23 !
ngamma = 0.75D0 ! -0.1
molm = 28.01D0        !molar mass
gdop = 115271.197D0

open(999, file = 'mol_params.txt', STATUS='OLD')
read(999,*) COstr !intensity
write(*,*) COstr, '  ***'
read(999,*) molm  !molar mass
write(*,*) molm, '  ***'
read(999,*) ngamma !T-dependence power
write(*,*) ngamma, '  ***'
read(999,*) frq0	!central frequency
write(*,*) frq0, '  ***'
read(999,*) a_partition	!file with the partition sums
write(*,*) a_partition, '***'
read(999,*) Elow !lower level energy, cm-1
write(*,*) Elow, '***'

read(999,*) a_subdir !lower level energy, cm-1
write(*,*) a_subdir, '***'
close(999)


a_partition=trim(a_partition)
open(998, file=a_partition, STATUS='OLD')

npart=400
!here we read partition sums from the file with the name specified in the mol_params.txt
!for we do not work with extra-high T-s, we limit maximum T as 400K, but if necessary one can change this limit
do ipart=1,npart
	read(998,*) tpart(ipart), qpart(ipart)
end do

close(998)

if (frq0.lt.1000.D0) then
	frq0 = frq0*1.D3
end if
										   
!COstr = COstr*Q_PART(tl,296.D0,npart,tpart,qpart)*

gdop = (frq0/clight)*dsqrt(2.D0*dlog(2.D0)/(molm*k_vs_aem))
!part of the Dopler width

dir = FILE$CURDRIVE
dirlength = GETDRIVEDIRQQ(dir)
!nconstpar = 13 !11 without CO line strength and f_rabi
!D2 for CO and Ar both were BANNED due to uselessness
!nconstpar = 17 ! see comments below
npar = 13
allocate(params_sd(npar), params(npar), params0(npar), dp(npar), dp_sd(npar))
allocate(JAC_FLAG(npar))
allocate(par_start_self(npar), par_start_foreign(npar))

!a_subdir='CO_230'

if (a_subdir.eq.'') then 
	a_folder = '\data\'
else
a_folder = '\data\'//trim(a_subdir)//'\'
end if

a_folder_out = '\data_out\'
a_script = trim(dir)//trim(a_folder_out)//'view_resid.gpl'
!a_list = trim(dir)//trim(a_folder)//'CO_list.txt'
a_list = 'CO_list.txt'
a_params = 'params_out.txt'

a_params_sd = trim(dir)//'\params_out_SD.txt'
a_params_cl = trim(dir)//'\params_out_CL.txt'


open(1,FILE=a_list, STATUS='OLD')
!open(10, FILE=a_params, STATUS = 'UNKNOWN', POSITION='APPEND')

do i=1,3  !first 3 lines in the list file are reserved for the header, need to skip
	read(1,*)a_header
	!write(*,*)trim(a_header)
end do

read(1,*) nfil !first line in CO_LIST is the number of files to process
write(*,*) nfil, ' files in the list'

1414 format(3(F9.3,2X),F15.6,2X,5(F11.6,2X),F15.6,2X,13(F11.6,2X))
1313 format(F15.6,2X,3(F9.6,2X))



open(2, file='params_start.txt', STATUS='OLD')
!params(:) = 0.D0
JAC_FLAG(:) = 1.D0
222 FORMAT(2(D12.3,2X))
do ipar=1,npar
	READ(2,*)par_start_self(ipar), par_start_foreign(ipar), JAC_FLAG(ipar), a_aux
	!write(*,*) par_start_self(ipar), par_start_foreign(ipar), JAC_FLAG(ipar)
end do
!READ(*,*)
close(2)


do i=1,nfil !just read line by line
	
	read(1,*)a_data,npts,plco,plar,tl,dl,tpl,LL	
	allocate(frq(npts), signal(npts), resid(npts), resid_sd(npts))
	allocate(model1(npts), model2(npts))
	if (tl.le.100.D0) tl = tl + 273.15

	write(*,*) 'temperature ', tl
	write(*,*) 'Dopler width', gdop*dsqrt(tl)
	!read(*,*)

	a_data_f = trim(dir)//trim(a_folder)//a_data
	write(*,*) 'READ FROM', a_data_f
	a_data_out = trim(dir)//trim(a_folder_out)//'r_'//a_data	
	a_corr_cl = trim(dir)//trim(a_folder_out)//'cor_cl_'//a_data
	a_corr_sd = trim(dir)//trim(a_folder_out)//'cor_sd_'//a_data
	!a_data_out = trim(a_folder_out)//'r_'//a_data	

	!!!! read the data from the file a_data_f
	open(2,file = a_data_f ,STATUS='OLD')
	do ipt = 1, npts
		read(2,*) frq(ipt), signal(ipt)
	end do
	close(2)

	if (tpl.eq.0) then
		frq = frq*1.D3 !resonator recordings frq in GHz, recalc to MHz
		signal = signal*1.D5 !1/cm to 1/km
	end if

	write(*,*) 'data read success from file ', trim(a_data)
	!write(*,*) 'type', tpl

	S_Tdep = Q_PART(tl,296.D0,npart,tpart,qpart)*exp(-c2*Elow/tl)*(1-exp(-c2*frq0/(tl*29979.25D0)))
	S_Tdep = S_Tdep/exp(-c2*Elow/296.D0)
    S_Tdep = S_Tdep/(1-exp(-c2*frq0/(296.D0*29979.25D0)))
	NVC = plco*0.3D0

	call init_params(npts, npar, frq, signal, plco, plar, tl, dl, tpl, LL, params0)
!        init_params(npts, npar, frq, signal, plco, tl, dl, tpl, LL, params0)

	params0(2) = (plco*par_start_self(2)+plar*par_start_foreign(2))
	params0(3) = (plco*par_start_self(3)+plar*par_start_foreign(3))
	params0(5) = (plco*par_start_self(5)+plar*par_start_foreign(5))
	params0(7) = par_start_self(7)
	params0(12) = par_start_self(12)
    params0(13) = (plco*par_start_self(13)+plar*par_start_foreign(13))
	!write(*,*) params0(13)	
	!read(*,*)
	call MDL(npts, npar, frq, signal, model1, resid, plco, tl, dl, tpl, LL, params0)
	!write(*,*) '!'
	call MDL(npts, npar, frq, signal, model1, resid_sd, plco, tl, dl, tpl, LL, params0)
	!write(*,*) '!!'

	open(13, file=a_data_out, status = 'UNKNOWN')

	do i1 = 1, npts
	    if (tpl.eq.0) then
			write(13,1313) (frq(i1)-params0(1))*1.D-3, signal(i1), resid_sd(i1), resid(i1)
		else
			write(13,1313) frq(i1)-params0(1), signal(i1), resid_sd(i1), resid(i1)
		end if
	end do
	close(13)
	!goto 99
	!read(*,*)
	!write(*,*) '!!!'
	call RUNFIT(npts, npar, frq, signal, plco, plar, tl, dl, tpl, LL, params0, 0, params, dp, resid, JAC_FLAG, 'fit_log.txt', 'rms_log.txt', a_params_cl, a_corr_cl)

	!params0 = params
	!params0(3) = (plco*par_start_self(3)+plar*par_start_foreign(3))

	call RUNFIT(npts, npar, frq, signal, plco, plar, tl, dl, tpl, LL, params0, 1, params_sd, dp_sd, resid_sd, JAC_FLAG, 'fit_log.txt', 'rms_log.txt', a_params_sd, a_corr_sd)
		
	!write(10,*) plco, plar, tl, params0(1), params0(2)

	write(*,*) 'write residual to file ', a_data_out

	open(13, file=a_data_out, status = 'UNKNOWN')
	
	call MDL(npts, npar, frq, signal, model1, resid, plco, tl, dl, tpl, LL, params)
	call MDL(npts, npar, frq, signal, model1, resid_sd, plco, tl, dl, tpl, LL, params_sd)

	do i1 = 1, npts
	    if (tpl.eq.0) then
			write(13,1313) (frq(i1)-params(1))*1.D-3, signal(i1)/maxval(signal), resid_sd(i1)/maxval(signal), resid(i1)/maxval(signal)
		else
			write(13,1313) frq(i1)-params(1), signal(i1)/maxval(signal), resid_sd(i1)/maxval(signal), resid(i1)/maxval(signal)
		end if
	end do

	close(13)

99	open(14, file = 'params_out_essential.txt', status='UNKNOWN', position = 'APPEND')

	t_pl = (tl/296.)**ngamma

	!if (tpl.eq.0) then
	!	params(2)=params(2)*1.D-3
	!	dp(2) = dp(2)*1.D-3 
	!	params_sd(2)=params_sd(2)*1.D-3
	!	dp_sd(2) = dp_sd(2)*1.D-3
	!	params_sd(3)=params_sd(3)*1.D-3
	!	dp_sd(3) = dp_sd(3)*1.D-3
	!	params_sd(4)=params_sd(4)*1.D-3
	!	dp_sd(4) = dp_sd(2)*1.D-3
	!end if
	!								 frq		f_err		G		   G_err     Y          Y_err 	  frq		   f_err		G0		      G0_err      G2              G2_err     D2                D2_err     Y                Y_err	I SI	  err    I SD          err 
	write(14,1414)	plco, plar, tl, params(1), dp(1), params(2)*t_pl, dp(2), params(5)*t_pl, dp(5), params_sd(1), dp_sd(1),params_sd(2)*t_pl, dp_sd(2),params_sd(3)*t_pl, dp_sd(3),params_sd(4)*t_pl, dp_sd(4),params_sd(5)*t_pl, dp_sd(5), params(6), dp(6),params_sd(6), dp_sd(6)

	close(14)

	qual    = maxval(signal)/dsqrt(sum(resid(:)**2)/(npts-1))
	qual_sd = maxval(signal)/dsqrt(sum(resid_sd(:)**2)/(npts-1))
	a_data_out = 'r_'//trim(a_data)

	if (i.eq.1) then
		open(15, file=a_script, status = 'UNKNOWN')
		write(15,*) 'load "header.gpl"'
		write(15,*) 'K = 50.'		
		write(15,*) 'S_v = 1.'
		write(15,*) 'S_sd = 0.'
		write(15,*) 'S_data = 3.'
		write(15,*) 'set multiplot layout 1,', nfil
		write(15,*) 'xpco = 0.05' !universal positions of labels showing pressure and fit quality
		write(15,*) 'ypco = 0.9'
		write(15,*) 'xpar = 0.05'
		write(15,*) 'ypar = 0.85'
		write(15,*) 'xqlu = 0.05'
		write(15,*) 'yqlu = 0.45'
		write(15,*) 'xqld = 0.05'
		write(15,*) 'yqld = 0.15'
	else
		open(15, file=a_script, status = 'UNKNOWN', position='APPEND')
	end if

	write(15,*) 'unset label'
	apc=''
	apa=''
	aq=''
	aqs=''
	write(apc,'(f9.3)') plco
	write(apa,'(f9.3)') plar
	write(aq,'(I8)') int(qual)
	write(aqs,'(I8)') int(qual_sd)
	if (i.eq.1) then
		write(15,*) 'set label "P_{CO} =', trim(apc), '" at graph xpco, graph ypco'
		write(15,*) 'set label "P_{Ar} =', trim(apa), '" at graph xpar, graph ypar'
		write(15,*) 'set label "Q =', trim(aq), '" at graph xqlu, graph yqlu'
		write(15,*) 'set label "Q =', trim(aqs), '" at graph xqld, graph yqld'
	else
		write(15,*) 'set label "', trim(apc), '" at graph xpco, graph ypco'
		write(15,*) 'set label "', trim(apa), '" at graph xpar, graph ypar'
		write(15,*) 'set label "', trim(aq), '" at graph xqlu, graph yqlu'
		write(15,*) 'set label "', trim(aqs), '" at graph xqld, graph yqld'
	end if
	if (tpl.eq.1) then
		ixtics = int(abs(dl)*5.D0) !int((params(2)*pco(i)+params(4)*par(i))*10)
	elseif (tpl.eq.2) then
		ixtics = int(abs(dl)*7.D0)
	elseif (tpl.eq.3) then
		ixtics = int(abs(maxval(frq)-params(1))/3.) !0.1*floor(abs(maxval(frq)-params(1))*10/3.) !
	elseif (tpl.eq.0) then
		ixtics = 10
	end if
	if (ixtics.eq.0) then
		write(15,*) 'set xtics 0.4'
	else 
		write(15,*) 'set xtics ', ixtics
	end if
	ixr = int(2.5*float(ixtics))
	if (ixr.eq.0) ixr = 1
	aixr=''
	write(aixr,'(I3)') ixr
	write(15,*) 'set xrange [-', trim(aixr), ':', trim(aixr), ']'
	write(15,'(A128)') 'plot "' //trim(a_data_out)// '" using 1:($2+S_data) w p ls 1 ti "",\' 
	write(15,'(A128)') '     "'// trim(a_data_out)// '" using 1:($4*K+S_v) w l ls 3 ti "",\'
	write(15,'(A128)') '     "'// trim(a_data_out)// '" using 1:($3*K+S_sd) w l ls 3 ti ""'
	write(15,*) ''

	close(15)

	deallocate(frq, signal, resid, resid_sd)
	deallocate(model1, model2)
end do

close(1)
!close(10)
deallocate(params_sd, params, params0, dp, dp_sd)

write(*,*) 'Dopler width ', gdop
write(*,*) 'PRESS ENTER'
read(*,*)

end	program CO_SINGLEFIT


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine MDL(npts, npar, frq, signal, model1, resid, plco, tl, dl, tpl, LL, params )

interface

double precision function SDRP (frq_in, f0_in,Y_in,gamma0_in,gamma2_in,delta0_in,delta2_in,strength_in)
double precision frq_in, f0_in,Y_in,gamma0_in,gamma2_in,delta0_in,delta2_in,strength_in
end function SDRP

double precision function difSDRP(frq_in, f0_in,Y_in,gamma0_in,gamma2_in,delta0_in,delta2_in,strength_in,dev)
double precision frq_in, f0_in,Y_in,gamma0_in,gamma2_in,delta0_in,delta2_in,strength_in,dev
end function difSDRP

double precision function difSDV(sg,sg0,GamD,Gam0,Gam2,Shift0,Shift2,strength_in,dev)
double precision sg,sg0,GamD,Gam0,Gam2,Shift0,Shift2,strength_in,dev
end function difSDV

double precision function difSDVe(sg,sg0,GamD,Gam0,Gam2,Shift0,Shift2,strength_in,dev, scale,pow,tpl)
double precision sg,sg0,GamD,Gam0,Gam2,Shift0,Shift2,strength_in,dev, scale,pow
integer tpl
end function difSDVe

double precision function difSDVe_add(sg,sg0,GamD,Gam0,Gam2,Shift0,Shift2,strength_in,dev, scale,pow,amp,NVC,tpl)
double precision sg,sg0,GamD,Gam0,Gam2,Shift0,Shift2,strength_in,dev, scale,pow,amp,NVC
integer tpl
end function difSDVe_add

subroutine HTP(sg0,GamD,Gam0,Gam2,Shift0,Shift2,anuVC,eta,sg,LS_pCqSDHC_R,LS_pCqSDHC_I)
double precision sg0,GamD,Gam0,Gam2,Shift0,Shift2,anuVC,eta,sg,LS_pCqSDHC_R,LS_pCqSDHC_I
end subroutine HTP


end interface

!===============================!

integer, intent(in):: npts, npar, tpl

double precision, intent(in):: params(npar) 
double precision, intent(in):: plco, tl, dl, LL
double precision, intent(in):: frq(npts), signal(npts)

double precision, intent(out)::	resid(npts), model1(npts)

double precision tmpS, tmpPow, tmpDop, tmpG0, tmpG2, tmpD0, tmpD2, tmpY, tmpBL0, tmpBL1, tmpBL2
double precision tmpScl

double precision  k_st, COstr, pi, frabi, ngamma, gdop, Lcell, S_Tdep, NVC
integer nconstpar
COMMON /CO_const/ k_st, COstr, pi, frabi, nconstpar, ngamma, gdop, Lcell, S_Tdep, NVC ! Elow, c2, frq0, npart, tpart, qpart

tmpD0 = 0.D0
tmpD2 = 0.D0
tmpY = 0.D0
tmpBL0 = 0.D0
tmpBL1 = 0.D0
tmpBL2 = 0.D0
tmpScl = 1.D0
tmpPow = 0.D0

tmpDop = gdop*dsqrt(tl)	
tmpG0 = params(2)   !HWHM (classic Gamma or SD Gamma0)
tmpG2 = params(3)  	!SD-component of HWHM Gamma2 for QSD... profiles
tmpD2 = params(4)	!SD-component of center Delta2 for QSD... profiles
tmpScl = params(6)	!Scaling factor for line area
tmpBL0 = params(8)	!baseline constant term
tmpBL1 = params(9)	!baseline linear term
tmpBL2 = params(10)	!baseline square term

tmpS = COstr*S_Tdep  !Q_PART(tl,296.D0,npart,tpart,qpart)*exp(-c2*Elow/tl)*(1-exp(-c2*frq0/(tl*29979.25D0)))	
!tmpS = tmpS/exp(-c2*Elow/296.D0)
!tmpS = tmpS/(1-exp(-c2*frq0/(296.D0*29979.25D0)))

if (tpl.eq.0) then
	tmpS = tmpS*K_st*plco*1.D5/(tl) !*4.343
	tmpY = params(5) ! 1-st order mixing patameter
end if

if (tpl.ge.1) then
	tmpS = tmpS*K_st*plco*LL/(tl)			
	tmpPow = params(7) !correction due to radiation source power depending on frequency
	tmpY = 0.D0	
end if



do i = 1, npts

if (tpl.eq.0) then	
	model1(i)=SDRP(frq(i), params(1),tmpY,tmpG0,tmpG2,0.D0,tmpD2,tmpS)*tmpScl 
	model1(i) = model1(i)+tmpBL0
	if (tmpBL1.ne.0.D0) model1(i) = model1(i)+tmpBL1*(frq(i)-params(1))
	if (tmpBL2.ne.0.D0) model1(i) = model1(i)+tmpBL2*frq(i)**2 !(frq(i)-params(1))**2
end if


if (tpl.eq.1) then		
	model1(i) = difSDVe_add(frq(i), params(1),tmpDop,tmpG0,tmpG2,tmpD0,tmpD2,tmpS,dl,tmpScl,tmpPow, params(12),params(13),tpl)
	!if (tmpBL0.ne.0.D0) BL constant
	model1(i) = model1(i)+tmpBL0
	!if (tmpBL1.ne.0.D0) BL linear
	model1(i) = model1(i)+tmpBL1*(frq(i)-params(1))
	!if (tmpBL2.ne.0.D0) BL square
	model1(i) = model1(i)+tmpBL2*(frq(i)-params(1))**2
	! BL cube
	!model1(i) = model1(i)+params(11)*(frq(i)-params(1))**3
end if

if (tpl.eq.2) then
	model1(i) = difSDVe_add(frq(i), params(1),tmpDop,tmpG0,tmpG2,tmpD0,tmpD2,tmpS,dl,tmpScl,tmpPow, params(12),params(13),tpl)
	if (tmpBL0.ne.0.D0) model1(i) = model1(i)+tmpBL0
	if (tmpBL1.ne.0.D0) model1(i) = model1(i)+tmpBL1*(frq(i)-params(1))
	if (tmpBL2.ne.0.D0) model1(i) = model1(i)+tmpBL2*(frq(i)-params(1))**2
	model1(i) = model1(i)+params(11)*(frq(i)-params(1))**3
end if

!qRe1 = (1.D0-dexp(-(qRe1+amp*qRe3)*strength_in))*(1+pow0*(sg+dev-sg0))
if (tpl.eq.3) then
	!call qSDV(params(1),tmpDop,tmpG0,tmpG2,tmpD0,tmpD2,frq(i),model1(i),tmpY)
	!write(*,*) i
	call HTP(params(1),tmpDop,tmpG0,tmpG2,tmpD0,tmpD2,params(13), 0.D0,frq(i),model1(i),tmpY)
	!write(*,*) i, model1(i)
	!model1(i) = tmpScl*(1-dexp(-model1(i)*tmpS))*(1.D0+tmpPow*(frq(i)-params(1)))
	model1(i) = model1(i)*tmpS*tmpScl*(1.D0+tmpPow*(frq(i)-params(1)))
	if (tmpBL0.ne.0.D0) model1(i) = model1(i)+tmpBL0
	if (tmpBL1.ne.0.D0) model1(i) = model1(i)+tmpBL1*(frq(i)-params(1))
	if (tmpBL2.ne.0.D0) model1(i) = model1(i)+tmpBL2*(frq(i)-params(1))**2
	model1(i) = model1(i)+params(11)*(frq(i)-params(1))**3
end if

end do

resid(:) = signal(:)-model1(:)

end subroutine MDL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine JACOBI(npts, npar, frq, signal, plco, tl, dl, tpl, LL, params,JAC,JAC_FLAG)

interface

double precision function SDRP (frq_in, f0_in,Y_in,gamma0_in,gamma2_in,delta0_in,delta2_in,strength_in)
double precision frq_in, f0_in,Y_in,gamma0_in,gamma2_in,delta0_in,delta2_in,strength_in
end function SDRP

double precision function difSDRP(frq_in, f0_in,Y_in,gamma0_in,gamma2_in,delta0_in,delta2_in,strength_in,dev)
double precision frq_in, f0_in,Y_in,gamma0_in,gamma2_in,delta0_in,delta2_in,strength_in,dev
end function difSDRP

double precision function difSDV(sg,sg0,GamD,Gam0,Gam2,Shift0,Shift2,strength_in,dev)
double precision sg,sg0,GamD,Gam0,Gam2,Shift0,Shift2,strength_in,dev
end function difSDV

double precision function difSDVe(sg,sg0,GamD,Gam0,Gam2,Shift0,Shift2,strength_in,dev,scale,pow,tpl)
double precision sg,sg0,GamD,Gam0,Gam2,Shift0,Shift2,strength_in,dev,scale,pow
integer tpl
end function difSDVe

double precision function difSDVe_add(sg,sg0,GamD,Gam0,Gam2,Shift0,Shift2,strength_in,dev, scale,pow,amp,NVC,tpl)
double precision sg,sg0,GamD,Gam0,Gam2,Shift0,Shift2,strength_in,dev, scale,pow,amp,NVC
integer tpl
end function difSDVe_add

double precision function Q_PART(Tq,Tref,npart,tpart,qpart)
double precision Tq, Tref
integer npart
double precision tpart(npart),qpart(npart)
end function Q_PART

end interface

integer, intent(in):: npts,npar
integer, intent(in):: tpl 

double precision, intent(in):: params(npar) 
double precision, intent(in):: tl, dl, LL

double precision, intent(in):: plco
double precision, intent(in):: frq(npts), signal(npts),JAC_FLAG(npar)

double precision, intent(out)::JAC(npts,npar)

double precision tmppar(npar), noblpar(npar), dpar, tmpmodel1(npts),tmpmodel2(npts), resid(npts)
double precision tmpS, tmpScl, tmpDop, tmpPow, tmpG0, tmpG2, tmpD0, tmpD2, tmpY

double precision  k_st, COstr, pi, frabi, ngamma, gdop, Lcell, S_Tdep, NVC
integer nconstpar
COMMON /CO_const/ k_st, COstr, pi, frabi, nconstpar, ngamma, gdop, Lcell, S_Tdep, NVC ! Elow, c2, frq0



3333 format(I10, F15.3,2X,F12.6,2X,F12.6)

JAC = 0.D0
noblpar(:) = params(:)
dpar = 1.D-5

do ipar = 8,11 
	noblpar(ipar)=0.D0
end do

!!!params 1 to 11 (f0, G0-2, D0, Y)
call MDL(npts, npar, frq, signal, tmpmodel1, resid, plco, tl, dl, tpl, LL, noblpar)

do ipar = 1,5
    if (JAC_FLAG(ipar).eq.1) then
		tmppar(:) = noblpar(:)
		tmppar(ipar) = tmppar(ipar)+dpar	
		call MDL(npts, npar, frq, signal, tmpmodel2, resid, plco, tl, dl, tpl,LL, tmppar)	     
		JAC(:,ipar)=(tmpmodel2(:)-tmpmodel1(:))/dpar  
		!write(*,*) 'jacobi, parameter number', ipar, sum(JAC(:,ipar)**2)
	end if
end do

!derivative by the power factor parameter for RAD recordings

if ((tpl.ge.1).and.(JAC_FLAG(7).eq.1)) then
	tmppar(:) = noblpar(:)
	tmppar(7) = tmppar(7)+dpar	
	call MDL(npts, npar, frq, signal, tmpmodel2, resid, plco, tl, dl, tpl,LL, tmppar)
	JAC(:,7)=(tmpmodel2(:)-tmpmodel1(:))/dpar
end if

!derivative by the amplitude of the additional line from BWO 2-nd harmonic
if ((tpl.ge.1).and.(JAC_FLAG(12).eq.1)) then
	tmppar(:) = noblpar(:)
	tmppar(12) = tmppar(12)+dpar	
	call MDL(npts, npar, frq, signal, tmpmodel2, resid, plco, tl, dl, tpl,LL, tmppar)
	JAC(:,12)=(tmpmodel2(:)-tmpmodel1(:))/dpar
end if
if (abs(frq(1)-params(1)).le.14.) JAC(:,12) = 0.D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!derivative by the Nu_vel_chng parameter (Dicke effect)
!correlates with params(3), fit only with JAC(:,3)=0.D0
if ((tpl.ge.1).and.(JAC_FLAG(13).eq.1)) then
	tmppar(:) = noblpar(:)
	tmppar(13) = tmppar(13)+dpar	
	call MDL(npts, npar, frq, signal, tmpmodel2, resid, plco, tl, dl, tpl,LL, tmppar)
	JAC(:,13)=(tmpmodel2(:)-tmpmodel1(:))/dpar
end if

tmpS = COstr*S_Tdep

!tmpS = COstr*Q_PART(tl,296.D0,npart,tpart,qpart)*exp(-c2*Elow/tl)*(1-exp(-c2*frq0/(tl*29979.25D0)))
!tmpS = tmpS/exp(-c2*Elow/296.D0)
!tmpS = tmpS/(1-exp(-c2*frq0/(296.D0*29979.25D0)))

if (tpl.ge.1) then
	tmpS = tmpS*K_st*plco*LL/(tl)				
else
	tmpS = tmpS*K_st*plco*1.D5/(tl)				
end if

tmpDop = gdop*dsqrt(tl)
tmpD0 = 0.D0
tmpG0 = params(2)
tmpG2 = params(3)
tmpD2 = params(4)
tmpY = params(5)

do i=1,npts
if (JAC_FLAG(8).eq.1) JAC(i,8)=1.D0   !constant baseline part derivative = 1.
if (tpl.eq.1) then
	if (JAC_FLAG(6).eq.1) JAC(i,6) = difSDVe_add(frq(i),params(1),tmpDop,tmpG0,tmpG2,tmpD0,tmpD2,tmpS,dl,1.D0,tmpPow,params(12),params(13),tpl)
	if (JAC_FLAG(9).eq.1) JAC(i,9) =(frq(i)-params(1))
	if (JAC_FLAG(10).eq.1) JAC(i,10)=(frq(i)-params(1))**2
	if (JAC_FLAG(11).eq.1) JAC(i,11)=(frq(i)-params(1))**3
else if (tpl.eq.2) then
    if (JAC_FLAG(6).eq.1) then
	   JAC(i,6) = difSDVe_add(frq(i),params(1),tmpDop,tmpG0,tmpG2,tmpD0,tmpD2,tmpS,dl,params(6)+dpar,tmpPow,params(12),params(13),tpl)
	   JAC(i,6) = JAC(i,6) - difSDVe_add(frq(i),params(1),tmpDop,tmpG0,tmpG2,tmpD0,tmpD2,tmpS,dl,params(6),tmpPow,params(12),params(13),tpl)
	   JAC(i,6) = JAC(i,6)/dpar
	end if
	if (JAC_FLAG(9).eq.1) JAC(i,9) =(frq(i)-params(1))
	if (JAC_FLAG(10).eq.1) JAC(i,10)=(frq(i)-params(1))**2
	if (JAC_FLAG(11).eq.1) JAC(i,11)=(frq(i)-params(1))**3
else if (tpl.eq.3) then
	!!call qSDV(params(1),tmpDop,tmpG0,tmpG2,tmpD0,tmpD2,frq(i),JAC(i,6),tmpY)
	
	if (JAC_FLAG(6).eq.1) then
		call HTP(params(1),tmpDop,tmpG0,tmpG2,tmpD0,tmpD2,params(13), 0.D0,frq(i),JAC(i,6),tmpY)	
		JAC(i,6) = JAC(i,6)*tmpS
	end if
	!!JAC(i,6) = (1-dexp(-JAC(i,6)*tmpS))*(1.D0+tmpPow*(frq(i)-params(1)))
	if (JAC_FLAG(9).eq.1) JAC(i,9) =(frq(i)-params(1))
	if (JAC_FLAG(10).eq.1) JAC(i,10)=(frq(i)-params(1))**2
	if (JAC_FLAG(11).eq.1) JAC(i,11)=(frq(i)-params(1))**3
else if (tpl.eq.0) then
	if (JAC_FLAG(6).eq.1) JAC(i,6) = SDRP(frq(i),params(1),tmpY,tmpG0,tmpG2,tmpD0,tmpD2,tmpS) !*0.D0 !constant intensity
	if (JAC_FLAG(9).eq.1) JAC(i,9) =(frq(i)-params(1))
	if (JAC_FLAG(10).eq.1) JAC(i,10)=frq(i)**2 !(frq(i)-params(1))**2
	if (JAC_FLAG(11).eq.1) JAC(i,11)=(frq(i)-params(1))**3
endif
end do

end subroutine JACOBI

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function DIAG(A)
double precision, DIMENSION(:,:), INTENT(IN) :: A
double precision, DIMENSION(size(A,1),size(A,2)) :: DIAG
integer i

DIAG=0.d0
do i=1,min(size(A,1),size(A,2))
DIAG(i,i)=A(i,i)
end do

end function DIAG

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function IDENT(n)
integer n,i
double precision IDENT(n,n)

IDENT=0.d0

do i=1,n
IDENT(i,i)=1.d0
end do

end function IDENT

double precision function Q_PART(Tq,Tref,npart,tpart,qpart)
!returns ratio Q(Tref)/Q(Tq) for the intensity
double precision Tq, Tref
integer npart
double precision tpart(npart),qpart(npart)

double precision QTq, Qref, Tleft, Tright, Qleft, Qright

!first find the location of the Tref in the list
do ipart=1,npart
	if ((tpart(ipart).le.Tref).and.(tpart(ipart+1).gt.Tref)) then
	irleft = ipart
	exit !as soon as we found an element of tpart which is less or equal to Tref, we stop
	endif
end do

if (tpart(irleft).eq.Tref) then
	Qref = qpart(irleft) !if Tref coincides with the element of tpart, no calculations needed
else
	Tleft = tpart(irleft)
	Tright = tpart(irleft+1)
	Qleft = qpart(irleft)
	Qright = qpart(irleft+1)
	Qref = Qleft+(Qright-Qleft)*(Tref-Tleft)/(Tright-Tleft) 
	!if Tref between two elements of tpart, linear interpolation is used
endif

!now the same for Tq, the temperature to which intensity is recalculated
do ipart=1,npart
	if ((tpart(ipart).le.Tq).and.(tpart(ipart+1).gt.Tq)) then
	irleft = ipart
	exit !as soon as we found an element of tpart which is less or equal to Tq, we stop
	endif
end do

if (tpart(irleft).eq.Tq) then
	QTq = qpart(irleft) !if Tq coincides with the element of tpart, no calculations needed
else
	Tleft = tpart(irleft)
	Tright = tpart(irleft+1)
	Qleft = qpart(irleft)
	Qright = qpart(irleft+1)
	QTq = Qleft+(Qright-Qleft)*(Tq-Tleft)/(Tright-Tleft) 
	!if Tref between two elements of tpart, linear interpolation is used
endif


Q_PART = Qref/Qtq

end function
















