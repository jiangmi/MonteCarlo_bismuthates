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
 character(len=80) :: fname
 character(len=10) :: s1,s2,s3,s4,s5,s6,s7,s8
contains
 subroutine init_parameters()
 implicit none
 Nx = 2
 Ny = 2
 N = Nx*Ny*Norb
 NBi = Nx*Ny
 mu = -3.0d0
 es = -2.0d0
 ep = -6.0d0
 tps = 1.6d0!0.0d0
 d0 = 1.8d0
 alpha = 3000.0d0
 spring_const = 100.d0
 if_X_displace = 1

 beta_max = 10.0d0
 beta_min = 1.0d0
 num_beta_steps = 9

 nwarms = 10
 ninv = 10    ! print warmup progress per ninv steps
 nmeas = 10 
 nbin = 2

 !======================================================
 ! Initialize files for recording results VS temperature
 write(s1,'(i1)') Nx
 write(s2,'(f4.1)') mu
 write(s3,'(f4.1)') es
 write(s4,'(f4.1)') ep
 write(s5,'(f3.1)') tps
 write(s6,'(f3.1)') d0
 write(s7,'(i2)') nwarms
 write(s8,'(i2)') nmeas
 fname = 'Nx'//adjustl(trim(s1))//'_mu'//adjustl(trim(s2))//'_es'//adjustl(trim(s3))//'_ep'//adjustl(trim(s4)) &
         //'_tps'//adjustl(trim(s5))//'_d0'//adjustl(trim(s6))//'_nwarm'//adjustl(trim(s7))//'_nmeas'//adjustl(trim(s8))//'.txt'
 !write(*,*) fname
 open(unit=11,file='E_'//fname,status="replace")
 open(unit=12,file='sp_'//fname,status="replace")
 open(unit=13,file='bp_'//fname,status="replace")

 return
 end subroutine init_parameters
end module parameters
