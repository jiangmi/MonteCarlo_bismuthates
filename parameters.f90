!==============================================================================
! Module parameters defines all global parameters for our problem.
!==============================================================================
module parameters
 integer Nx 
 integer Ny 
 integer, parameter :: Norb = 3          !number of orbitals in the basis
 integer N                               !total number of orbitals
 integer Nbi                             !total number of Bi orbitals
 integer num_beta_steps
 integer nwarms, ninv, nmeas, nbin
 integer if_X_displace                   ! No displacement for checking code
 double precision mu
 double precision tps                    !O-Bi overlap integral
 double precision ep                     !Oxygen 2p site energy
 double precision es                     !Bi 6s site energy
 double precision beta
 double precision beta_max      
 double precision beta_min
 double precision spring_const           !harmonic potential
 double precision alpha                  !Anharmonic potential
 double precision d0                     !equilibrium bond distance
contains
 subroutine init_parameters()
 implicit none
 Nx = 8
 Ny = 8
 N = Nx*Ny*Norb
 NBi = Nx*Ny
 mu = -3.0d0
 es = -2.0d0
 ep = -6.0d0
 tps = 1.6d0!0.0d0
 d0 = 1.8d0
 alpha = 0.d0!3000.0d0
 spring_const = 0.d0!100.d0
 if_X_displace = 0

 beta_max = 1.0d0
 beta_min = 1.0d0
 num_beta_steps = 1

 nwarms = 1
 ninv = 1    ! print warmup progress per ninv steps
 nmeas = 4 
 nbin = 2
return
 end subroutine init_parameters
end module parameters
