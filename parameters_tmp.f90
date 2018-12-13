!==============================================================================
! Module parameters defines all global parameters for our problem.
!==============================================================================
module parameters
 integer Nx 
 integer Ny 
 integer, parameter :: Norb = 3          !number of orbitals in the basis
 integer N                               !total number of orbitals
 integer Nbi                             !total number of Bi orbitals
 integer nclass                          !number of distance class
 integer num_beta_steps
 integer nwarms, ninv, nmeas, nbin
 integer if_X_displace                   ! No displacement for checking code
 double precision mu
 double precision tsp                    !O-Bi overlap integral
 double precision tpp                    !O-O overlap integral
 double precision ep                     !Oxygen 2p site energy
 double precision es                     !Bi 6s site energy
 double precision beta
 double precision beta_max      
 double precision beta_min
 double precision spring_const           !harmonic potential
 double precision alpha                  !Anharmonic potential
 double precision d0                     !equilibrium bond distance
 character(len=80) :: fname
 character(len=10) :: s1,s2,s3,s4,s5,s6,s7,s8
contains
 subroutine init_parameters()
 implicit none
 Nx = Nxval
 Ny = Nyval
 N = Nx*Ny*Norb
 NBi = Nx*Ny
 nclass = (Nx/2+1)*(Ny/2+1)
 mu = muval
 es = esval
 ep = epval
 tsp = tspval
 tpp = tppval
 d0 = 1.0d0
 alpha = alphaval
 spring_const = springval
 if_X_displace = 1 

 beta_max = bemaxval
 beta_min = beminval
 num_beta_steps = nbetaval

 nwarms = nwarmval
 ninv = ninvval    ! print warmup progress per ninv steps
 nmeas = nmeasval 
 nbin = nbinval

 !======================================================
 ! Initialize files for recording results VS temperature
 500 format(a80)
 open(unit=11,file='dirval'//'/data_'//'fnameval'//'.txt',status="replace")
 write(unit=11,fmt=500) 'beta    n_avg    n_err    E_avg    E_err    sp_avg    sp_err    bp_avg    bp_err'

 return
 end subroutine init_parameters
end module parameters
