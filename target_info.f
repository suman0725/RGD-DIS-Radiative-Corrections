	subroutine target_info(theta_tar,theta,tb,ta,z,a,tgt_len,tgt_thick,aux,teff)
C+______________________________________________________________________________
!
! TARGET_INFO - Calculate material thickness before and after scattering.
!
! ARGUMENTS:
!
!   THETA_TAR:  (R*8):  Target angle - 0 deg is normal to beam, increase CCW.
!			(i.e. increasing angle faces SOS).
!   THETA:	(R*8):	Scattering angle in Degrees (Input).
!   TB:		(R*8):	Thickness BEFORE scatter, in radiation lengths.
!   TA:		(R*8):	Thickness AFTER  scatter, in radiation lengths.
!   Z:		(R*8):	'Z' of target nucleus.
!   A:		(R*8):	'A' of target nucleus.
!   tgt_len:    (R*8):  target length in cm.
!   M_TGT:	(R*8):	Mass in GeV/c^2 of target nucleus.
!   AUX:	(R*8):	Array of length 10 containing:
!			aux(1) = RESOL parameter for DEEPSIG smearing (not used).
!			aux(2) = Separation energy in GeV of target nucleus.
!			aux(3) = F(0)  parameter for F(y) function
!			aux(4) = BigB  parameter for F(y)
!			aux(5) = a     parameter for F(y)
!			aux(6) = b     parameter for F(y)
!			aux(7) = alpha parameter for F(y)
!
! AUX definitions for old f(y) model:
!		!	aux(1) = Mass of recoiling nucleon in GeV/c^2.
!		!	aux(2) = Separation energy in GeV of target nucleus.
!		!	aux(3) = Y0 parameter for F(y) function
!		!	aux(4) = B  parameter for F(y)
!		!	aux(5) = C  parameter for F(y)
!		!	aux(6) = D  parameter for F(y)
!		!	aux(7) = RESOL parameter for DEEPSIG smearing.
!
C-______________________________________________________________________________

	implicit none

C Get various constants.

	include 'math_physics.inc'
	include 'spect.inc'

C Declare arguments.

	real*8 theta,tb,ta,z,a,aux(7),theta_tar

C Local declarations.

	integer*4 i,j,itar
	logical found

	real*8 rho,tgt_len,tgt_thick,tgt_rl
	real*8 thrad,thrad_tar
	real*8 t1,t2,teff	!effective thickness for e+/e- pair production.

! Parameters for NE18 target calculations

	real*8 x0_al,rho_al,x0_cm_al
	real*8 X0_air,rho_air,x0_cm_air
	real*8 x0_mylar,rho_mylar,x0_cm_mylar
	real*8 x0_cm_kevlar
	parameter (X0_Al	= 24.01)
	parameter (rho_Al	= 2.70)
	parameter (X0_cm_Al	= X0_Al/rho_Al)
	parameter (X0_air	= 36.66)
	parameter (rho_air	= .001205)
	parameter (X0_cm_air	= X0_air/rho_air)
	parameter (X0_mylar	= 39.95)
	parameter (rho_mylar	= 1.39)
	parameter (X0_cm_mylar	= x0_mylar/rho_mylar)
	parameter (X0_cm_kevlar	= 55.2)		!Calculation from C.Keppel.

C Lookup table from which we find the radiation length and AUX parameters.
C DAVE G NOTE: Stuff for He3,Be,Cu mostly placeholders.
C I'm going to try and use the parameters from Atti, Faralli, and West 
C (nucl-th/9812078)

	real*8		lookup(11,10)/			!data vs. Z.

!  From Nadia's iteration: August 2007
!  From Suman's iteration: March 2026
!------------------------------------------------------------------------------------------
!        H     2H       3He     4He      Be      C       Al     Fe     Cu       Au      Sn
!-------------------------------------------------------------------------------------------
!
     > 1.,        1.,         2.,         2.,        4., 
     > 6.,        13.,    26.,    29.,        79.,    50.,  !Z's
     > 1.,        2.,         3.,         4.,        9., 
     > 12.,       27.,    56.,    64.,        197.,  118.,  !A's
     > 61.28,       122.6,      71.07,      94.32,     65.19,
     > 42.7,      24.01,  13.84,  12.86,      6.461,  8.820, !R.L.s in g/cm^2.
     > 0.00,        0.14,       0.10,       0.16,      0.10,  
     > 0.25,      0.25,   0.25,   0.10,       0.25,  0.25, !RESOL's
     > 0.00,        0.00225,    0.00549,    0.02020,   0.00928,  
     > 0.01727,   0.0099, 0.01060,0.00855,    0.00693, 0.00788, !E_SEP's in GeV.
     > 0.0056000,  0.00904157, 0.0053094,  0.0040197,  0.0034814,
     > 0.0031822,  0.0032783,  0.0028900,  0.00287403,  0.0026424,   0.0027788,  !f0
     > 0.0000000,  0.000852159, 0.0021843,  0.0013449,  0.0011608,
     > 0.0013591,  0.0013474,  0.0014016,  0.0008866,  0.0007632,   0.0008359,   !bigB
     > 0.0000000,  0.0077272, 0.0028864,  0.0026986,  0.0031195,
     > 0.0030265,  0.0029698,  0.0031802,  0.0030959,  0.0030654,    0.0030834,   !a
     > 0.0000000,  0.0093937, 0.0103492,  0.0074941,  0.0078398,
     > 0.0070505,  0.0065760,  0.0072635,  0.0070945,  0.0067678,      0.0069602, !b
     > 20.000,    45.384,    64.247,   100.256,   110.967,  
     > 137.285,   131.845,   165.700,   132.458,   132.452,   132.456/  !alpha


C ============================ Executable Code =================================

C Get target specific stuff from lookup table.

	rho = tgt_thick/tgt_len
	found = .false.
	do i = 1,11				!loop over known targets.
	  if (lookup(i,1).eq.z.and.float(int(a)).eq.lookup(i,2)) then
	    found = .true.
	    itar=i
	    tgt_rl = tgt_thick/lookup(i,3)	!Divide by radiation length.
	    do j = 1,7
	      aux(j) = lookup(i,j+3)
	    enddo
	  endif
	enddo
	if (.not.found) then
	  write(6,*) 'cant find target in lookup table!'
          return			!Quit if couldn't find info.
	endif

C Compute TA and TB for target.  Target Chamber windows added at end.

	thrad = theta*d_r		!Angle in radians.
	thrad_tar = theta_tar*d_r

	if (z.ge.4) then		!Solid target.
	  tb = tgt_rl/2./abs(cos(thrad_tar))
	  ta = tgt_rl/2./abs(cos(thrad+spect_sign*thrad_tar))

	elseif (z.eq.1..or.z.eq.2) then		!Liq. 1H or 2H target.
	  tb = tgt_rl/2.
     >       + t_can_front/X0_cm_Al

	  ta = tgt_rl/2. + t_can_back/X0_cm_Al  ! this is for tuna can targets
c	  ta = tgt_rl/2.16 + t_can_back/X0_cm_Al  ! this is for tuna can targets
cdg	  if (tgt_len.le.7) then		! 4 cm target
cdg	    if (theta.lt.32.) then		! Out the end.
cdg	      ta = tgt_rl/2.
cdg     >           + t_can_back/X0_cm_Al/cos(thrad) !NOT CORRECTED FOR CURVED ENDCAP
cdg	    else				! Out the side.
cdg	      ta = can_radius*rho/lookup(itar,3)/sin(thrad)
cdg     >           + t_can_side/X0_cm_Al/sin(thrad)
cdg	    endif
cdg	  else                                  !15 cm target, always side.
cdg	      ta =  can_radius*rho/lookup(itar,3)/sin(thrad)
cdg     >           + t_can_side/X0_cm_Al/sin(thrad)
cdg	  endif
	endif

! The parameterization for the pair production cross section is cross
! section per nucleon, per effective tgt len.  teff is (5+t1)*(5+t2)/10
! where t1,t2 are twice the thicknesses in radiation lengths (%) before and
! after the center of the target. (i.e. t1=tgt_rl/cos(theta_tar))

	t1 = 100.*(2.*tb)	!double and convert to %
	t2 = 100.*(2.*ta)	!	"	"
	teff = (5.+t1)*(5+t2)/10.

! Contribution due to vacuum chamber walls, air, and entrance to spectrometer.

	ta = ta + t_exitwin/X0_cm_Al		!Al exit window.
     >          + t_air/X0_cm_air 		!air before magnets 
     >          + t_mylar/X0_cm_mylar		!mylar entrance windows.
     >          + t_kevlar/X0_cm_kevlar		!kevlar entrance windows.


	return	

1001	format(a)
	end
