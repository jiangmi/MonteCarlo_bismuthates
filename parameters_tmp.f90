!==============================================================================
! Module parameters defines all global parameters for our problem.
!==============================================================================
module parameters
 integer Nx, Ncx
 integer Ny, Ncy 
 integer, parameter :: Norb = 3          !number of orbitals in the basis
 real*8,  parameter :: pi = 3.1415926
 integer N, Nc                           !total number of orbitals
 integer Nbi                             !total number of Bi orbitals
 integer nclass                          !number of distance class
 integer num_beta_steps
 integer nwarms, ninv, nmeas, nbin, nbeta
 integer if_X_displace                   ! No displacement for checking code
 integer if_print_MC
 integer travel_cluster                  ! if using travelling cluster sampling
 double precision mu
 double precision tsp                    !O-Bi overlap integral
 double precision tpp                    !O-O overlap integral
 double precision ep                     !Oxygen 2p site energy
 double precision es                     !Bi 6s site energy
 double precision dXamp                  !dX amplitude
 double precision beta
 double precision, dimension(:),   allocatable :: betas
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
 N  = Nx*Ny*Norb
 Ncx = Ncxval
 Ncy = Ncyval
 Nc  = Ncx*Ncy*Norb
 NBi = Nx*Ny
 nclass = (Nx+1)*(Ny+1)
 mu = muval
 es = esval
 ep = epval
 tsp = tspval
 tpp = tppval
 d0 = 1.0d0
 alpha = alphaval
 spring_const = springval
 if_X_displace = 1
 if_print_MC = 1
 dXamp = dXampval
 travel_cluster = 1

 nbeta = nbetaval
 allocate(betas(1:nbeta))
 betas = betasval

 nwarms = nwarmval
 ninv = ninvval    ! print warmup progress per ninv steps
 nmeas = nmeasval 
 nbin = nbinval

 !======================================================
 ! Initialize files for recording results VS temperature
 510 format(a170)
 open(unit=11,file='dirval'//'/data_'//'fnameval'//'.txt',status="replace")
 write(unit=11,fmt=510) 'beta     n_avg      n_err      E_avg      E_err      |X|      sp_avg      sp_err     bp_avg     bp_err    chi_s     chi_s_err    chi_px    chi_px_err   chi_py   chi_py_err'

 520 format(a142)
 open(unit=12,file='dirval'//'/data_sublat_'//'fnameval'//'.txt',status="replace")
 write(unit=12,fmt=520) 'beta     n_avg      n_Bi1      n_Bi2      n_px1      n_px2     n_py1      n_py2      n_A1g1     n_A1g2      Sp1        Sp2       Bp1       Bp2'

 550 format(a41)
 open(unit=7,file='dirval'//'/SpBp_r_'//'fnameval'//'.txt',status="replace")
 write(unit=7,fmt=550) 'dx    dy     Sp      err      Bp      err'

 open(unit=30,file='dirval'//'/Ek_'//'fnameval'//'.txt',status="replace")

 ! Print various quantities along the MC updates
 if (if_print_MC==1) then
   open(unit=26,file='dirval'//'/MC_Xp_'//'fnameval'//'.txt',status="replace")
   open(unit=13,file='dirval'//'/MC_XLs_'//'fnameval'//'.txt',status="replace")
   open(unit=14,file='dirval'//'/MC_ns_'//'fnameval'//'.txt',status="replace")
   open(unit=15,file='dirval'//'/MC_npx_'//'fnameval'//'.txt',status="replace")
   open(unit=16,file='dirval'//'/MC_nLs_'//'fnameval'//'.txt',status="replace")
   open(unit=17,file='dirval'//'/MC_nLd_'//'fnameval'//'.txt',status="replace")
   open(unit=18,file='dirval'//'/MC_nLx_'//'fnameval'//'.txt',status="replace")
   open(unit=19,file='dirval'//'/MC_nLy_'//'fnameval'//'.txt',status="replace")
   open(unit=20,file='dirval'//'/MC_Sp_XLs_'//'fnameval'//'.txt',status="replace")
   open(unit=21,file='dirval'//'/MC_Bp_XLs_'//'fnameval'//'.txt',status="replace")
   open(unit=22,file='dirval'//'/MC_Sp_Ls_'//'fnameval'//'.txt',status="replace")
   open(unit=23,file='dirval'//'/MC_Bp_Ls_'//'fnameval'//'.txt',status="replace")
   open(unit=24,file='dirval'//'/MC_Sp_'//'fnameval'//'.txt',status="replace")
   open(unit=25,file='dirval'//'/MC_Bp_'//'fnameval'//'.txt',status="replace")
 endif

 return
 end subroutine init_parameters
end module parameters
