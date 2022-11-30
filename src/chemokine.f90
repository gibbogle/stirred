! Chemokine data

module chemokine

use global

implicit none

!type receptor_type
!	character(8) :: name
!	logical :: used
!	integer :: chemokine
!	integer :: sign
!	real(REAL_KIND) :: strength
!	real(REAL_KIND) :: level(5)
!	real(REAL_KIND) :: saturation_threshold
!	real(REAL_KIND) :: refractory_time
!end type
!type(receptor_type), target :: receptor(MAX_RECEPTOR)

integer, parameter :: N1D = 1   !20

type chemokine_type
	character(24) :: name
	logical :: used
	logical :: present
	logical :: use_secretion
	logical :: constant
	logical :: controls_growth
	logical :: controls_death
	real(REAL_KIND) :: bdry_conc
	real(REAL_KIND) :: diff_coef
	real(REAL_KIND) :: membrane_diff_in
	real(REAL_KIND) :: membrane_diff_out
	logical :: decay
	real(REAL_KIND) :: halflife
	real(REAL_KIND) :: decay_rate
	real(REAL_KIND) :: max_cell_rate	! Vmax
	real(REAL_KIND) :: MM_C0			! Km
	real(REAL_KIND) :: Hill_N
	real(REAL_KIND) :: medium_diff_coef	! diffusion coefficient in the medium
!	real(REAL_KIND) :: medium_dlayer	! unstirred layer thickness
!	real(REAL_KIND) :: medium_M			! mass of constituent
!	real(REAL_KIND) :: medium_U			! total blob uptake rate
!	real(REAL_KIND) :: medium_Cext		! far-field concentration
!	real(REAL_KIND) :: medium_Cbnd		! boundary concentration
!	real(REAL_KIND) :: total_flux_prev
!	real(REAL_KIND) :: medium_Cbnd_prev

!	real(REAL_KIND) :: diff_reduction_factor
!	real(REAL_KIND), allocatable :: coef(:,:)
!	real(REAL_KIND), allocatable :: conc
!	real(REAL_KIND), allocatable :: grad(:,:,:,:)
!	real(REAL_KIND), allocatable :: Cave_b(:,:,:)
!	real(REAL_KIND), allocatable :: Cprev_b(:,:,:)
!	real(REAL_KIND), allocatable :: Fprev_b(:,:,:)
!	real(REAL_KIND), allocatable :: Fcurr_b(:,:,:)
	real(REAL_KIND) :: Cmedium(N1D)     ! the vertical conc distribution (index 1 is the well bottom)
end type
type(chemokine_type), target :: chemo(MAX_CHEMO)

type ODEdiff_type
!	integer :: ichemo
!	integer :: nextra
!	integer :: nintra
!	integer :: nvars
!	integer, allocatable :: ivar(:,:,:)
!	integer, allocatable :: varsite(:,:)
!	integer, allocatable :: icoef(:,:)
!	integer, allocatable :: iexcoef(:,:)
!	integer, allocatable :: vartype(:)
!	integer, allocatable :: cell_index(:)
!	integer, allocatable :: isite_extra(:)
!	integer, allocatable :: isite_intra(:)
!	integer, allocatable :: extra_isite(:)
!!	integer, allocatable :: ncoef(:)
	real(REAL_KIND) :: deltaC_soft
	real(REAL_KIND) :: k_soft
	real(REAL_KIND) :: C1_soft
end type
type(ODEdiff_type) :: ODEdiff

integer :: nchemo, chemomap(MAX_CHEMO)
real(REAL_KIND) :: Conc(2*MAX_CHEMO)    ! average IC and EC concs (where in monolayer model EC is at the cell layer)
real(REAL_KIND) :: Cmediumave(MAX_CHEMO)    ! average medium concentrations
!real(REAL_KIND) :: CglucoseMedium(N1D)
!real(REAL_KIND) :: CdrugMedium(2,0:2,N1D)

contains

!----------------------------------------------------------------------------------------
! Convert halflife in hours to a decay rate /sec
!----------------------------------------------------------------------------------------
real(REAL_KIND) function DecayRate(halflife)
real(REAL_KIND) :: halflife

if (halflife == 0) then		! No decay
	DecayRate = 0
else
	DecayRate = log(2.0)/(halflife*60*60)
endif
end function

!----------------------------------------------------------------------------------------
! Note units:
! distance		cm
! volume		cm^3
! time			s
! diff coeff	cm^2.s^-1
! mass			mol
! concentration	mM where 1 mM = 1.0e-6 mol.cm^-3
! consumption	mol.cell^-1.s^-1
! production	mol.cell^-1.s^-1
! Changed 28/07/2015 to use mumol for mass, mumol/s for flux
!-----------------------------------------------------------
! mass			mumol = 10^-6 mol
! flux			mumol/s
! concentration	mM where 1 mM = 1 mumol.cm^-3
! Note: concentration*volume = mM*cm3 = mumol
!----------------------------------------------------------------------------------------
subroutine SetupChemo
integer :: ichemo

chemo(OXYGEN)%name = 'Oxygen'
chemo(GLUCOSE)%name = 'Glucose'
chemo(TRACER)%name = 'Tracer'
chemo(OXYGEN)%decay_rate = 0
chemo(GLUCOSE)%decay_rate = 0
chemo(TRACER)%decay_rate = 0

Conc = 0

do ichemo = 1,MAX_CHEMO
	chemo(ichemo)%present = .false.
	if (chemo(ichemo)%used) then
!		chemo(ichemo)%medium_dlayer = d_layer	! for now, using the same unstirred layer thickness for all constituents
		if (ichemo == OXYGEN .or. ichemo == GLUCOSE .or. ichemo == TRACER) then
			chemo(ichemo)%present = .true.
		endif
	endif
enddo
do ichemo = 1,2
	Conc(ichemo) = chemo(ichemo)%bdry_conc
	Conc(MAX_CHEMO + ichemo) = chemo(ichemo)%bdry_conc
	chemo(ichemo)%Cmedium = chemo(ichemo)%bdry_conc
enddo
!call AllocateConcArrays
!if (use_FD) then
!	do ichemo = 1,MAX_CHEMO
!		if (.not.chemo(ichemo)%used) cycle
!		if (allocated(chemo(ichemo)%Cave_b)) deallocate(chemo(ichemo)%Cave_b)
!		if (allocated(chemo(ichemo)%Cprev_b)) deallocate(chemo(ichemo)%Cprev_b)
!		if (allocated(chemo(ichemo)%Fprev_b)) deallocate(chemo(ichemo)%Fprev_b)
!		if (allocated(chemo(ichemo)%Fcurr_b)) deallocate(chemo(ichemo)%Fcurr_b)
!		allocate(chemo(ichemo)%Cave_b(NXB,NYB,NZB))
!		allocate(chemo(ichemo)%Cprev_b(NXB,NYB,NZB))
!		allocate(chemo(ichemo)%Fprev_b(NXB,NYB,NZB))
!		allocate(chemo(ichemo)%Fcurr_b(NXB,NYB,NZB))
!		chemo(ichemo)%diff_reduction_factor = 0.5               ! default value, may need to be adjusted
!	enddo
!endif
!chemo(OXYGEN)%diff_reduction_factor = 0.3
end subroutine

!----------------------------------------------------------------------------------------
! Consumption rate = Rmax*C/(MM_C0 + C)  (Michaelis-Menten)
! where Km = MM_C0 >= Rmax*T*10^6/Vextra_cm3 is the O2 concentration at which cell uptake is halved,
! and T is the maximum expected time step
! This is to ensure that all the O2 in the site is not depleted in a time step, making
! O2 conc go negative.  For now it seems reasonable to assume that glucose uptake
! varies in proportion to O2 uptake, i.e. we only need Oxygen M-M
!----------------------------------------------------------------------------------------
subroutine SetMMParameters

!chemo(OXYGEN)%MM_C0 = 0.00133		! 1 mmHg = 1.33 uM (Kevin suggests 1 - 2 uM)
!chemo(GLUCOSE)%MM_C0 = chemo(OXYGEN)%MM_C0*chemo(GLUCOSE)%max_cell_rate/chemo(OXYGEN)%max_cell_rate
write(logmsg,'(a,e12.4)') 'Oxygen MM_C0: ',chemo(OXYGEN)%MM_C0
call logger(logmsg)
write(logmsg,'(a,e12.4)') 'Glucose MM_C0: ',chemo(GLUCOSE)%MM_C0
call logger(logmsg)
end subroutine

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
!subroutine AllocateConcArrays
!integer :: ic
!
!write(logmsg,*) 'AllocateConcArrays'
!call logger(logmsg)
!do ic = 1,MAX_CHEMO
!	if (chemo(ic)%used) then
!		if (allocated(chemo(ic)%conc)) then
!			call logger("chemo(ic)%conc already allocated")
!			deallocate(chemo(ic)%conc)
!!			stop
!		endif
!		allocate(chemo(ic)%conc(NX,NY,NZ))
!		if (allocated(chemo(ic)%grad)) then
!			call logger("chemo(ic)%grad already allocated")
!			deallocate(chemo(ic)%grad)
!!			stop
!		endif
!		allocate(chemo(ic)%grad(3,NX,NY,NZ))
!		chemo(ic)%conc = chemo(ic)%bdry_conc	
!	endif
!enddo
!write(logmsg,*) 'did AllocateConcArrays'
!call logger(logmsg)
!end subroutine

!----------------------------------------------------------------------------------------
! Set the concentrations at a site when it is on the boundary
!----------------------------------------------------------------------------------------
!subroutine SetBdryConcs(site)
!integer :: site(3)
!integer :: ichemo
!
!do ichemo = 1,MAX_CHEMO
!	if (chemo(ichemo)%used) then
!		chemo(ichemo)%conc(site(1),site(2),site(3)) = chemo(ichemo)%bdry_conc
!	endif
!enddo
!end subroutine

!----------------------------------------------------------------------------------------
! Reset the concentrations at a site that is no longer on the boundary
!----------------------------------------------------------------------------------------
!subroutine ResetConcs(site0)
!integer :: site0(3)
!integer :: site(3), ichemo, k, n
!real(REAL_KIND) :: csum(MAX_CHEMO)
!
!csum = 0
!n = 0
!do k = 1,27
!	if (k == 14) cycle
!	site = site0 + jumpvec(:,k)
!	if (occupancy(site(1),site(2),site(3))%indx(1) == OUTSIDE_TAG) cycle
!!	if (associated(occupancy(site(1),site(2),site(3))%bdry)) cycle
!	n = n+1
!	do ichemo = 1,MAX_CHEMO
!		if (chemo(ichemo)%used) then
!			csum(ichemo) = csum(ichemo) + chemo(ichemo)%conc(site(1),site(2),site(3))
!		endif
!	enddo
!enddo
!
!do ichemo = 1,MAX_CHEMO
!	if (chemo(ichemo)%used) then
!		chemo(ichemo)%conc(site0(1),site0(2),site0(3)) = csum(ichemo)/n
!	endif
!enddo
!end subroutine

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
!subroutine ShowConcs
!integer :: x, y, z, k, x1, site(3), ichemo
!real(REAL_KIND) :: C(MAX_CHEMO,100)
!logical :: start
!
!y = blob_centre(2) 
!z = blob_centre(3)
!start = .false.
!k = 0
!do x = blobrange(1,1),blobrange(1,2)
!	if (occupancy(x,y,z)%indx(1) == OUTSIDE_TAG) then
!		if (start) then
!			exit
!		else
!			cycle
!		endif
!	endif
!	if (.not.start) x1 = x
!	start = .true.
!	k = k+1
!	do ichemo = 1,MAX_CHEMO
!		if (.not.chemo(ichemo)%used) cycle
!	    C(ichemo,k) = chemo(ichemo)%conc(x,y,z)
!	enddo
!enddo
!do ichemo = 1,MAX_CHEMO
!	if (.not.chemo(ichemo)%used) cycle
!enddo
!!if (C(1) < 100) then
!!	write(logmsg,*) 'First conc < 100: istep: ',istep, C(1),occupancy(x1,y,z)%indx(1)
!!	call logger(logmsg)
!!	site = (/x1,y,z/)
!!	write(logmsg,'(a,L)') 'isbdry: ',isbdry(site)
!!	call logger(logmsg)
!!	stop
!!endif
!end subroutine

!----------------------------------------------------------------------------------------
! A drug and its metabolites continue to be solved for if any one has a medium boundary 
! concentration that exceeds Cthreshold
!----------------------------------------------------------------------------------------
subroutine UpdateChemomap
integer :: ichemo, iparent, idrug, im
logical :: present

!write(*,*) 'UpdateChemomap'
nchemo = 0
do ichemo = 1,MAX_CHEMO
	if (chemo(ichemo)%present) then
		nchemo = nchemo + 1
		chemomap(nchemo) = ichemo
	endif
enddo
end subroutine

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
subroutine CheckDrugPresence
integer :: idrug, iparent, im, ichemo, nmet

! Check drugs and their metabolites
do idrug = 1,ndrugs_used
	nmet = drug(idrug)%nmetabolites
	iparent = DRUG_A + (MAX_METAB+1)*(idrug-1)
	if (chemo(iparent)%present) then		! simulation with this drug has started
		if (.not.drug_gt_cthreshold(idrug)) then
			write(logmsg,'(a,i2,a,a)') 'Removing drug and metabolites, concentrations below threshold: ',idrug,' ',chemo(iparent)%name
			call logger(logmsg)
!			write(logmsg,'(a,4e12.3)') 'concs: ',(chemo(iparent+im)%medium_Cbnd,im=0,nmet)
!			call logger(logmsg)
			write(logmsg,'(a,e12.3)') 'threshold concentration: ',Cthreshold
			call logger(logmsg)
			do im = 0,nmet
				ichemo = iparent + im
				chemo(ichemo)%present = .false.
			enddo
			call RemoveDrug(idrug)
		endif
	endif
enddo
end subroutine

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
subroutine RemoveDrug(idrug)
integer :: idrug
integer :: iparent, ichemo, im, kcell, nmet
type(cell_type), pointer :: cp

iparent = DRUG_A + (MAX_METAB+1)*(idrug-1)
nmet = drug(idrug)%nmetabolites
do im = 0,nmet
	ichemo = iparent + im
!	chemo(ichemo)%medium_Cbnd = 0
	chemo(ichemo)%Cmedium = 0
	Conc(ichemo) = 0
	Conc(MAX_CHEMO+ichemo) = 0
enddo
do kcell = 1,nlist
	if (cell_list(kcell)%state == DEAD) cycle
    cp => cell_list(kcell)
	do im = 0,nmet
		ichemo = iparent + im
		cp%Cin(ichemo) = 0
	enddo
enddo
end subroutine

!----------------------------------------------------------------------------------
! Computes metabolism rate as a fraction of the maximum cell rate
! Use the "soft landing" option for Hill_N = 1 if MM_threshold = 0
!----------------------------------------------------------------------------------
real(REAL_KIND) function O2_metab(C)
integer :: ichemo
real(REAL_KIND) :: C
real(REAL_KIND) :: metab

ichemo = OXYGEN
if (ichemo == OXYGEN) then
	if (chemo(ichemo)%Hill_N == 2) then
		if (C > 0) then
			metab = C*C/(chemo(ichemo)%MM_C0*chemo(ichemo)%MM_C0 + C*C)
		else
			metab = 0
		endif
	else
		if (MM_THRESHOLD > 0) then
			if (C > ODEdiff%C1_soft) then
				metab = (C-ODEdiff%deltaC_soft)/(chemo(ichemo)%MM_C0 + C - ODEdiff%deltaC_soft)
			elseif (C > 0) then
				metab = ODEdiff%k_soft*C*C
			else
				metab = 0
			endif
		else
			if (C > 0) then
				metab = C/(chemo(ichemo)%MM_C0 + C)
			else
				metab = 0
			endif
		endif
	endif
endif
O2_metab = metab
end function

!----------------------------------------------------------------------------------
! Computes metabolism rate as a fraction of the maximum cell rate
!----------------------------------------------------------------------------------
real(REAL_KIND) function glucose_metab(C)
real(REAL_KIND) :: C
integer :: N
real(REAL_KIND) :: MM_C0

N = chemo(GLUCOSE)%Hill_N
MM_C0 = chemo(GLUCOSE)%MM_C0
if (C > 0) then
	glucose_metab = C**N/(MM_C0**N + C**N)
else
	glucose_metab = 0
endif
end function



end module