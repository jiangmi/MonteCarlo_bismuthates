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
 Nx = 4
 Ny = 4
 N = Nx*Ny*Norb
 NBi = Nx*Ny
 nclass = (Nx/2+1)*(Ny/2+1)
 mu = 9.0
 es = 6.42
 ep = 2.42
 tsp = 2.08
 tpp = 0.056
 d0 = 1.0d0
 alpha = 300
 spring_const = 10
 if_X_displace = 1

 beta_max = 10.0d0
 beta_min = 10.0d0
 num_beta_steps = 1

 nwarms = 1000
 ninv = 10    ! print warmup progress per ninv steps
 nmeas = 5000 
 nbin = 10

 !======================================================
 ! Initialize files for recording results VS temperature
 open(unit=11,file='E_'//'Nx4_mu9.0_es6.42_ep2.42_tsp2.08_tpp0.056_alpha300_c10'//'.txt',status="replace")
 open(unit=12,file='sp_'//'Nx4_mu9.0_es6.42_ep2.42_tsp2.08_tpp0.056_alpha300_c10'//'.txt',status="replace")
 open(unit=13,file='bp_'//'Nx4_mu9.0_es6.42_ep2.42_tsp2.08_tpp0.056_alpha300_c10'//'.txt',status="replace")

 return
 end subroutine init_parameters
end module parameters
