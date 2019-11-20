module measurements
 use parameters, only: nbin
 integer cnt  ! count the measurements in each bin
 integer o1, o2

 ! accumulated quantites
 double precision aenergy, ax, ax2
 double precision antot, anbis_avg, anox, anoy
 double precision anpLs_avg  !Number of oxygen holes with Ls symmetry
 double precision anpLd_avg, anpLx_avg, anpLy_avg
 double precision asp_site_avg
 double precision abp_site_avg
 double precision anbi1, anbi2, anpx1, anpx2, anpy1, anpy2
 double precision anpLs1, anpLs2, asp1, asp2, abp1, abp2
 double precision, dimension(:), allocatable  :: aEk
 double precision, dimension(:), allocatable  :: anbis
 double precision, dimension(:), allocatable  :: anpLs
 double precision, dimension(:), allocatable  :: anpLd
 double precision, dimension(:), allocatable  :: anpLx
 double precision, dimension(:), allocatable  :: anpLy
 double precision, dimension(:), allocatable  :: aspolaron
 double precision, dimension(:), allocatable  :: abpolaron
 double precision, dimension(:), allocatable  :: aspolaron_ij
 double precision, dimension(:), allocatable  :: abpolaron_ij
 double precision, dimension(:), allocatable  :: achi_sc
 double precision, dimension(:,:,:), allocatable  :: achi_ch

 !utility arrays
 double precision, dimension(:), allocatable  :: fermi
 double precision, dimension(:,:), allocatable  :: ns
 double precision, dimension(:,:), allocatable  :: nLs
 double precision, dimension(:,:), allocatable  :: nr

 !Below site-dependent quantities along MC updates
 double precision, dimension(:), allocatable  :: ns_mc
 double precision, dimension(:), allocatable  :: npx_mc
 double precision, dimension(:), allocatable  :: nLs_mc
 double precision, dimension(:), allocatable  :: nLd_mc
 double precision, dimension(:), allocatable  :: nLx_mc
 double precision, dimension(:), allocatable  :: nLy_mc
 double precision, dimension(:), allocatable  :: Sp_mc
 double precision, dimension(:), allocatable  :: Bp_mc
 double precision, dimension(:), allocatable  :: Sp1_mc
 double precision, dimension(:), allocatable  :: Bp1_mc

 !Below are for binned quantites
 double precision, dimension(:), allocatable  :: benergy
 double precision, dimension(:), allocatable  :: bx
 double precision, dimension(:), allocatable  :: bx2
 double precision, dimension(:), allocatable  :: bnbis_avg
 double precision, dimension(:), allocatable  :: bnox
 double precision, dimension(:), allocatable  :: bnoy
 double precision, dimension(:), allocatable  :: bntot
 double precision, dimension(:), allocatable  :: bnpLs_avg
 double precision, dimension(:), allocatable  :: bnpLd_avg
 double precision, dimension(:), allocatable  :: bnpLx_avg
 double precision, dimension(:), allocatable  :: bnpLy_avg
 double precision, dimension(:), allocatable  :: bsp_site_avg
 double precision, dimension(:), allocatable  :: bbp_site_avg

 double precision, dimension(:,:), allocatable  :: bEk
 double precision, dimension(:,:), allocatable  :: bnbis
 double precision, dimension(:,:), allocatable  :: bnpLs
 double precision, dimension(:,:), allocatable  :: bnpLd
 double precision, dimension(:,:), allocatable  :: bnpLx
 double precision, dimension(:,:), allocatable  :: bnpLy
 double precision, dimension(:,:), allocatable :: bspolaron
 double precision, dimension(:,:), allocatable :: bbpolaron
 double precision, dimension(:,:), allocatable :: bspolaron_ij
 double precision, dimension(:,:), allocatable :: bbpolaron_ij
 double precision, dimension(:,:), allocatable  :: bchi_sc
 double precision, dimension(:,:,:,:), allocatable  :: bchi_ch

 !sublattice dependent quantities
 double precision, dimension(:), allocatable  :: bnbi1
 double precision, dimension(:), allocatable  :: bnbi2
 double precision, dimension(:), allocatable  :: bnpx1
 double precision, dimension(:), allocatable  :: bnpx2
 double precision, dimension(:), allocatable  :: bnpy1
 double precision, dimension(:), allocatable  :: bnpy2
 double precision, dimension(:), allocatable  :: bnpLs1
 double precision, dimension(:), allocatable  :: bnpLs2
 double precision, dimension(:), allocatable  :: bsp1
 double precision, dimension(:), allocatable  :: bsp2
 double precision, dimension(:), allocatable  :: bbp1
 double precision, dimension(:), allocatable  :: bbp2
contains
 !=============================================================================
 subroutine get_err(bins,mean,std)
 use parameters, only: nbin
 implicit none
 double precision, dimension(1:nbin) :: bins
 double precision mean, std
 mean = sum(bins(:))/dfloat(nbin)
 std = sqrt(sum((bins(:)-mean)**2.0d0)/dfloat(nbin-1))
 return
 end subroutine get_err
 !=============================================================================
 subroutine allocate_quantities()
 use parameters
 implicit none
 allocate(fermi(0:N-1))
 allocate(ns(0:Nbi-1, 0:N-1))
 allocate(nLs(0:Nbi-1, 0:N-1))
 allocate(nr(0:Nbi-1, 0:N-1))
 allocate(aEk(0:N-1))
 allocate(anbis(0:Nbi-1))
 allocate(anpLs(0:Nbi-1))
 allocate(anpLd(0:Nbi-1))
 allocate(anpLx(0:Nbi-1))
 allocate(anpLy(0:Nbi-1))
 allocate(aspolaron(0:Nbi-1))
 allocate(abpolaron(0:Nbi-1))
 allocate(aspolaron_ij(0:nclass-1))
 allocate(abpolaron_ij(0:nclass-1))
 allocate(achi_sc(0:Norb-1))
 allocate(achi_ch(0:nclass-1, 0:Norb-1, 0:Norb-1))
 allocate(benergy(1:nbin))
 allocate(bx(nbin))
 allocate(bx2(nbin))
 allocate(bnbis_avg(nbin))
 allocate(bnox(nbin))
 allocate(bnoy(nbin))
 allocate(bntot(nbin))
 allocate(bnpLs_avg(nbin))
 allocate(bnpLd_avg(nbin))
 allocate(bnpLx_avg(nbin))
 allocate(bnpLy_avg(nbin))
 allocate(bnbi1(nbin))
 allocate(bnbi2(nbin))
 allocate(bnpx1(nbin))
 allocate(bnpx2(nbin))
 allocate(bnpy1(nbin))
 allocate(bnpy2(nbin))
 allocate(bnpLs1(nbin))
 allocate(bnpLs2(nbin))
 allocate(bsp1(nbin))
 allocate(bsp2(nbin))
 allocate(bbp1(nbin))
 allocate(bbp2(nbin))
 allocate(bEk(0:N-1, nbin))
 allocate(bnbis(0:Nbi-1, nbin))
 allocate(bnpLs(0:Nbi-1, nbin))
 allocate(bnpLd(0:Nbi-1, nbin))
 allocate(bnpLx(0:Nbi-1, nbin))
 allocate(bnpLy(0:Nbi-1, nbin))
 allocate(ns_mc(0:Nbi-1))
 allocate(npx_mc(0:Nbi-1))
 allocate(nLs_mc(0:Nbi-1)) 
 allocate(nLd_mc(0:Nbi-1))
 allocate(nLx_mc(0:Nbi-1))
 allocate(nLy_mc(0:Nbi-1))
 allocate(Sp_mc(0:Nbi-1))
 allocate(Bp_mc(0:Nbi-1))
 allocate(Sp1_mc(0:Nbi-1))
 allocate(Bp1_mc(0:Nbi-1))
 allocate(bsp_site_avg(nbin))
 allocate(bbp_site_avg(nbin))
 allocate(bspolaron(0:Nbi-1, nbin))
 allocate(bbpolaron(0:Nbi-1, nbin))
 allocate(bspolaron_ij(0:nclass-1, nbin))
 allocate(bbpolaron_ij(0:nclass-1, nbin))
 allocate(bchi_sc(0:Norb-1, nbin))
 allocate(bchi_ch(0:nclass-1, 0:Norb-1, 0:Norb-1, nbin))
end subroutine allocate_quantities
 !=============================================================================
 subroutine deallocate_quantities()
 implicit none
 deallocate(fermi)
 deallocate(ns, nLs, nr)
 deallocate(aEk, bEk)
 deallocate(anbis, anpLs, anpLd, anpLx, anpLy)
 deallocate(aspolaron, abpolaron)
 deallocate(aspolaron_ij, abpolaron_ij)
 deallocate(achi_sc, achi_ch)
 deallocate(bnbis_avg, bnpLs_avg, bnpLd_avg, bnpLx_avg, bnpLy_avg)
 deallocate(bnbis, bnpLs, bnpLd, bnpLx, bnpLy)
 deallocate(bsp_site_avg, bbp_site_avg)
 deallocate(bnbi1, bnbi2)
 deallocate(bnpx1, bnpx2)
 deallocate(bnpy1, bnpy2)
 deallocate(bnpLs1, bnpLs2)
 deallocate(bsp1, bsp2)
 deallocate(bbp1, bbp2)
 deallocate(ns_mc, npx_mc, nLs_mc, Sp_mc, Bp_mc, Sp1_mc, Bp1_mc)
 deallocate(bspolaron, bbpolaron)
 deallocate(bspolaron_ij, bbpolaron_ij)
 deallocate(bchi_sc, bchi_ch)
 end subroutine deallocate_quantities
 !=============================================================================
 subroutine zero_accumulators()
 use parameters
 implicit none
 cnt = 0
 aenergy = 0.0d0
 aX = 0.0d0
 ax2 = 0.0d0
 anox = 0.0d0
 anoy = 0.0d0
 aEk = 0.0d0
 anbis = 0.0d0
 antot = 0.0d0
 anpLs = 0.0d0
 anpLd = 0.0d0
 anpLx = 0.0d0
 anpLy = 0.0d0
 anbis_avg = 0.0d0
 anpLs_avg = 0.0d0
 anpLd_avg = 0.0d0
 anpLx_avg = 0.0d0
 anpLy_avg = 0.0d0
 asp_site_avg = 0.d0
 abp_site_avg = 0.d0
 aspolaron = 0.0d0
 abpolaron = 0.0d0
 anbi1 = 0.0d0
 anbi2 = 0.0d0
 anpx1 = 0.0d0
 anpx2 = 0.0d0
 anpy1 = 0.0d0
 anpy2 = 0.0d0
 anpLs1 = 0.0d0
 anpLs2 = 0.0d0
 asp1 = 0.0d0
 asp2 = 0.0d0
 abp1 = 0.0d0
 abp2 = 0.0d0
 aspolaron_ij = 0.0d0
 abpolaron_ij = 0.0d0
 achi_sc = 0.0d0
 achi_ch = 0.0d0
 end subroutine zero_accumulators
 !=============================================================================
 subroutine populate_bins(bin)
 implicit none
 integer bin
 !print*, 'Populating bin ', bin
 benergy(bin) = aenergy/dfloat(cnt)
 bX(bin) = aX/dfloat(cnt)
 bX2(bin) = aX2/dfloat(cnt)
 bntot(bin) = antot/dfloat(cnt)
 bnbis_avg(bin) = anbis_avg/dfloat(cnt)
 bnox(bin) = anox/dfloat(cnt)
 bnoy(bin) = anoy/dfloat(cnt)
 bnpLs_avg(bin) = anpLs_avg/dfloat(cnt)
 bnpLd_avg(bin) = anpLd_avg/dfloat(cnt)
 bnpLx_avg(bin) = anpLx_avg/dfloat(cnt)
 bnpLy_avg(bin) = anpLy_avg/dfloat(cnt)
 bEk(:,bin) = aEk/dfloat(cnt)
 bnbis(:,bin) = anbis/dfloat(cnt)
 bnpLs(:,bin) = anpLs/dfloat(cnt)
 bnpLd(:,bin) = anpLd/dfloat(cnt)
 bnpLx(:,bin) = anpLx/dfloat(cnt)
 bnpLy(:,bin) = anpLy/dfloat(cnt)
 bsp_site_avg(bin) = asp_site_avg/dfloat(cnt)
 bbp_site_avg(bin) = abp_site_avg/dfloat(cnt)

 bnbi1(bin) = anbi1/dfloat(cnt)
 bnbi2(bin) = anbi2/dfloat(cnt)
 bnpx1(bin) = anpx1/dfloat(cnt)
 bnpx2(bin) = anpx2/dfloat(cnt)
 bnpy1(bin) = anpy1/dfloat(cnt)
 bnpy2(bin) = anpy2/dfloat(cnt)
 bnpLs1(bin) = anpLs1/dfloat(cnt)
 bnpLs2(bin) = anpLs2/dfloat(cnt)
 bsp1(bin) = asp1/dfloat(cnt)
 bsp2(bin) = asp2/dfloat(cnt)
 bbp1(bin) = abp1/dfloat(cnt)
 bbp2(bin) = abp2/dfloat(cnt)

 bspolaron(:,bin) = aspolaron/dfloat(cnt)
 bbpolaron(:,bin) = abpolaron/dfloat(cnt)
 bspolaron_ij(:,bin) = aspolaron_ij/dfloat(cnt)
 bbpolaron_ij(:,bin) = abpolaron_ij/dfloat(cnt)
 bchi_sc(:,bin) = achi_sc/dfloat(cnt)
 bchi_ch(:,:,:,bin) = achi_ch/dfloat(cnt)
 return
 end subroutine populate_bins

 !============================================================================
 ! do_measurements below
 !============================================================================
 subroutine do_measurements(X)
 use parameters
 use cluster, only: return_index_for_coordinates, dclass, dclass_F, expqr, phase
 use monte_carlo, only: compute_total_E, get_H
 implicit none
 integer ix, iy, jx, jy, n1, n2, n3, n4
 integer ixp, iyp, ixm, iym
 integer jxp, jyp, jxm, jym
 integer i, j, k, info, i1, i2, j1, j2
 integer lwork
 double precision, dimension(0:N-1,0:N-1) :: U
 double precision, allocatable, dimension(:) :: work
 double precision Xi_Ls, Xi_Ld, Xi_Lx, Xi_Ly, Xj_Ls
 double precision a_innLd, a_innLx, a_innLy
 double precision a_innLd2, a_innLx2, a_innLy2
 double precision tmp1, tmp2, tmp3
 double precision energy, fac, fac1, factor
 double precision, dimension(0:N-1) :: X
 double precision, dimension(0:N-1) :: Ek

 !Zero the record for each MC measurement
 ns_mc = 0.0d0
 npx_mc = 0.0d0
 nLs_mc = 0.0d0
 nLd_mc = 0.0d0 
 nLx_mc = 0.0d0
 nLy_mc = 0.0d0
 Sp_mc = 0.0d0
 Bp_mc = 0.0d0
 Sp1_mc = 0.0d0
 Bp1_mc = 0.0d0

 fermi = 0.0d0
 ns = 0.0d0
 nLs = 0.d0
 nr = 0.d0

 cnt = cnt + 1
 call compute_total_E(energy,X)
 aenergy = aenergy + energy/N
 do i = Nbi,N-1
  aX  = aX  + X(i)/(N-NBi)
  aX2 = aX2 + X(i)*X(i)/(N-NBi)
 enddo 
 
 !Measure the electronic quantities
 call get_H(U,X)

 lwork = 3*N-1
 allocate(work(1:lwork))
 call dsyev('V','U',N,U,N,Ek,work,lwork,info)
 if (info.ne.0) then
    print *,"dsyev error", N, lwork, info
 endif

 aEk = aEk + Ek

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!  Obtain useful arrays  !!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 do n1 = 0,N-1
   fermi(n1) = 1.0d0/(exp(beta*(Ek(n1)-mu))+1.0d0)
 enddo
 
 !loop over Bi sites
 do ix = 0,Nx-1
  do iy = 0,Ny-1
   i = return_index_for_coordinates(ix,iy,0)
   ixp = return_index_for_coordinates(ix,iy,1)
   iyp = return_index_for_coordinates(ix,iy,2)
   ixm = return_index_for_coordinates(ix-1,iy,1)
   iym = return_index_for_coordinates(ix,iy-1,2)

   !sum over eigenstates
   do n1 = 0,N-1
     ns(i,n1)  = U(i,n1)*U(i,n1)
     nLs(i,n1) = (0.5d0*(U(ixp,n1) - U(ixm,n1) + U(iyp,n1) - U(iym,n1)))**2.0
     nr(i,n1)  = ns(i,n1) + nLs(i,n1)
   enddo
  enddo
 enddo

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! start computing quantities
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 do ix = 0,Nx-1
  do iy = 0,Ny-1
   i = return_index_for_coordinates(ix,iy,0)
   ixp = return_index_for_coordinates(ix,iy,1) 
   iyp = return_index_for_coordinates(ix,iy,2) 
   ixm = return_index_for_coordinates(ix-1,iy,1)
   iym = return_index_for_coordinates(ix,iy-1,2)

   Xi_Ls = 0.5d0*(X(ixm) - X(ixp) + X(iym) - X(iyp))
   Xi_Ld = 0.5d0*(X(ixm) - X(ixp) - X(iym) + X(iyp))
   Xi_Lx = (X(ixm) + X(ixp))/sqrt(2.0)
   Xi_Ly = (X(iym) + X(iyp))/sqrt(2.0)

   !sum over eigenstates
   do n1 = 0,N-1
     fac = 2.0d0*fermi(n1)/Nbi
     tmp1 = fac*ns(i,n1)
     tmp2 = fac*U(ixp,n1)*U(ixp,n1)
     tmp3 = fac*U(iyp,n1)*U(iyp,n1)
     anbis(i) = anbis(i) + tmp1*Nbi
     anbis_avg = anbis_avg + tmp1
     anox  = anox  + tmp2
     anoy  = anoy  + tmp3
     antot = antot + tmp1+tmp2+tmp3

     a_innLd = 0.5d0*(U(ixp,n1) - U(ixm,n1) - U(iyp,n1) + U(iym,n1))
     a_innLx = (U(ixp,n1) + U(ixm,n1))/sqrt(2.0)
     a_innLy = (U(iyp,n1) + U(iym,n1))/sqrt(2.0)
     a_innLd2 = a_innLd*a_innLd
     a_innLx2 = a_innLx*a_innLx
     a_innLy2 = a_innLy*a_innLy

     anpLs(i)  = anpLs(i) + fac*nLs(i,n1)*Nbi
     anpLd(i)  = anpLd(i) + fac*a_innLd2*Nbi
     anpLx(i)  = anpLx(i) + fac*a_innLx2*Nbi
     anpLy(i)  = anpLy(i) + fac*a_innLy2*Nbi
     anpLs_avg = anpLs_avg + fac*nLs(i,n1)
     anpLd_avg = anpLd_avg + fac*a_innLd2
     anpLx_avg = anpLx_avg + fac*a_innLx2
     anpLy_avg = anpLy_avg + fac*a_innLy2

     !record quantities for each MC measurement
     if (if_print_MC==1) then
       ns_mc(i) = ns_mc(i) + tmp1*Nbi
       npx_mc(i) = npx_mc(i) + tmp2*Nbi
       nLs_mc(i) = nLs_mc(i) + fac*nLs(i,n1)*Nbi
       nLd_mc(i) = nLd_mc(i) + fac*a_innLd2*Nbi
       nLx_mc(i) = nLx_mc(i) + fac*a_innLx2*Nbi
       nLy_mc(i) = nLy_mc(i) + fac*a_innLy2*Nbi
     endif

     !Also compute orbital occupancy for two sublattices
     ! /0.5 below accounts for Nbi/2 sites of sublattice
     if (mod(ix+iy,2)==0) then
       anbi1 = anbi1 + tmp1/0.5
       anpx1 = anpx1 + tmp2/0.5
       anpy1 = anpy1 + tmp3/0.5
       anpLs1 = anpLs1 + fac*nLs(i,n1)/0.5
     else
       anbi2 = anbi2 + tmp1/0.5
       anpx2 = anpx2 + tmp2/0.5
       anpy2 = anpy2 + tmp3/0.5
       anpLs2 = anpLs2 + fac*nLs(i,n1)/0.5
     endif

     !!!!!!!!!!!!!!!!!!!!!!!!!!
     !!   single polaron     !!
     !!!!!!!!!!!!!!!!!!!!!!!!!!
     ! aspolaron does not need Nbi because of (i) index
     tmp1 = 2.0d0*fermi(n1)*Xi_Ls*nr(i,n1)
     aspolaron(i) = aspolaron(i) + tmp1
     asp_site_avg = asp_site_avg + tmp1/Nbi

     !record quantities for each MC measurement
     if (if_print_MC==1) then
       Sp_mc(i)  = Sp_mc(i) + tmp1
       Sp1_mc(i) = Sp1_mc(i) + 2.0d0*fermi(n1)*nLs(i,n1)
     endif

     !Also compute single polaron density for two sublattices
     ! /0.5 below accounts for Nbi/2 sites of sublattice
     if (mod(ix+iy,2)==0) then
       asp1 = asp1 + tmp1/(Nbi/2.0)
     else
       asp2 = asp2 + tmp1/(Nbi/2.0)
     endif

     ! second sum over eigenstates for bipolaron (Bp) and Sp correlation
     do n2 = 0,N-1
       !!!!!!!!!!!!!!!!!!!!!!!!!!
       !!       bipolaron      !!
       !!!!!!!!!!!!!!!!!!!!!!!!!!
       fac1 = 2.0d0*fermi(n2)/Nbi
       tmp1 = fermi(n1)*fermi(n2)* Xi_Ls*nr(i,n1)*nr(i,n2)
       abpolaron(i) = abpolaron(i) + tmp1
       abp_site_avg = abp_site_avg + tmp1/Nbi

       !record quantities for each MC measurement
       if (if_print_MC==1) then
         Bp_mc(i) = Bp_mc(i) + tmp1
         Bp1_mc(i) = Bp1_mc(i) + fermi(n1)*fermi(n2)* nLs(i,n1)*nLs(i,n2)
       endif

       !Also compute bipolaron density for two sublattices
       ! /0.5 below accounts for Nbi/2 sites of sublattice
       if (mod(ix+iy,2)==0) then
         abp1 = abp1 + tmp1/(Nbi/2.0)
       else
         abp2 = abp2 + tmp1/(Nbi/2.0)
       endif

       do jx = 0,Nx-1
         do jy = 0,Ny-1
           j = return_index_for_coordinates(jx,jy,0)
           jxp = return_index_for_coordinates(jx  ,jy  ,1)
           jxm = return_index_for_coordinates(jx-1,jy  ,1)
           jyp = return_index_for_coordinates(jx  ,jy  ,1)
           jym = return_index_for_coordinates(jx  ,jy-1,1)

           k = dclass(i,j)
          ! print*, 'meas ', k,i,j

           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           !!  single polaron staggered correlation  !!
           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           Xj_Ls = -0.5d0*(X(jxp) - X(jxm) + X(jyp) - X(jym))
           tmp1 = fermi(n1)*fermi(n2)*Xi_Ls*Xj_Ls*phase(i,j)/dclass_F(k)
           factor = 4.0*tmp1
           aspolaron_ij(k) = aspolaron_ij(k) + factor*nr(i,n1)*nr(j,n2)

           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           !!  bipolaron staggered correlation  !!
           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           ! need additional two sum over eigenstates for bipolaron etc.
           ! X_iLs*X_jLs * n_iup*n_idn*n_jup*n_jdn
           tmp2 = 0.0
           do n3 = 0,N-1
             do n4 = 0,N-1
               factor = tmp1*fermi(n3)*fermi(n4)
               tmp2 = tmp2 + factor* nr(i,n1)*nr(i,n2)*nr(j,n3)*nr(j,n4)
             enddo
           enddo

           abpolaron_ij(k) = abpolaron_ij(k) + tmp2

           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           !!  correct Sp staggered correlation by Bp  !!
           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           ! need additional two sum over eigenstates for bipolaron etc.
           ! X_iLs*X_jLs * n_iup*n_jup*n_jdn temporary term
           tmp3 = 0.0
           do n3 = 0,N-1
             factor = 2.0*tmp1*fermi(n3)
             tmp3 = tmp3 - factor* (nr(i,n1)*nr(j,n2)*nr(j,n3)  &
                                   +nr(j,n1)*nr(i,n2)*nr(i,n3))
           enddo

           ! correction
           aspolaron_ij(k) = aspolaron_ij(k) + 4.0*tmp2 + tmp3

           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                                                 
           !!      s-wave susceptibilities     !!                                                 
           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
           do o1 = 0,Norb-1
             i1 = return_index_for_coordinates(ix,iy,o1)
             i2 = return_index_for_coordinates(jx,jy,o1)
             achi_sc(o1) = achi_sc(o1) + fermi(n1)*fermi(n2)/Nbi   &
                 *U(i1,n1)*U(i1,n1)*U(i2,n2)*U(i2,n2)
           enddo

           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                                                 
           !!      charge susceptibilities     !!                                                 
           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
           do o1 = 0,Norb-1
             do o2 = 0,Norb-1
               j1 = return_index_for_coordinates(ix,iy,0)
               j2 = return_index_for_coordinates(jx,jy,0)
               i1 = return_index_for_coordinates(ix,iy,o1)
               i2 = return_index_for_coordinates(jx,jy,o2)
               achi_ch(k,o1,o2) = achi_ch(k,o1,o2) + 4.0d0*fermi(n1)*fermi(n2)/Nbi   &
                   *expqr(k,j1,j2)*U(i1,n1)*U(i1,n1)*U(i2,n2)*U(i2,n2)
                  ! *U(i1,n1)*U(i1,n1)*U(i2,n2)*U(i2,n2)
             enddo
           enddo
         enddo
       enddo  !end loop jx, jy
     enddo !end loop n2
    enddo  !end loop n1
    
    !Print quantities for each MC measurement
    if (if_print_MC==1) then
      530 format(f13.5)
      540 format(f13.5,f13.5,f13.5,f13.5)
      write(unit=26,fmt=540)  X(ixp), X(ixm), X(iyp), X(iym)
      write(unit=13,fmt=530)  Xi_Ls
      write(unit=14,fmt=530)  ns_mc(i)
      write(unit=15,fmt=530)  npx_mc(i)
      write(unit=16,fmt=530)  nLs_mc(i)
      write(unit=17,fmt=530)  nLd_mc(i)
      write(unit=18,fmt=530)  nLx_mc(i)
      write(unit=19,fmt=530)  nLy_mc(i)
      write(unit=20,fmt=530)  Sp_mc(i)/Xi_Ls
      write(unit=21,fmt=530)  Bp_mc(i)/Xi_Ls
      write(unit=22,fmt=530)  Sp1_mc(i)
      write(unit=23,fmt=530)  Bp1_mc(i)
      write(unit=24,fmt=530)  Sp_mc(i)
      write(unit=25,fmt=530)  Bp_mc(i)
    endif

   enddo
 enddo   !end loop ix, iy

 deallocate(work)

 return
 end subroutine do_measurements


end module measurements
