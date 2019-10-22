module measurements
 use parameters, only: nbin
 integer cnt  ! count the measurements in each bin
 integer o1, o2

 ! accumulated quantites
 double precision aenergy, ax, ax2
 double precision antot, anbis_avg, anox, anoy
 double precision anpa1g_avg  !Number of oxygen holes with A1g symmetry
 double precision asp_site_avg
 double precision abp_site_avg
 double precision anbi1, anbi2, anpx1, anpx2, anpy1, anpy2
 double precision anpA1g1, anpA1g2, asp1, asp2, abp1, abp2
 double precision, dimension(:), allocatable  :: anbis
 double precision, dimension(:), allocatable  :: anpa1g
 double precision, dimension(:), allocatable  :: aspolaron
 double precision, dimension(:), allocatable  :: abpolaron
 double precision, dimension(:), allocatable  :: aspolaron_ij
 double precision, dimension(:), allocatable  :: achi_sc
 double precision, dimension(:,:,:), allocatable  :: achi_ch

 !Below site-dependent quantities along MC updates
 double precision, dimension(:), allocatable  :: ns_mc
 double precision, dimension(:), allocatable  :: npx_mc
 double precision, dimension(:), allocatable  :: nLs_mc
 double precision, dimension(:), allocatable  :: Sp_mc
 double precision, dimension(:), allocatable  :: Bp_mc

 !Below are for binned quantites
 double precision, dimension(:), allocatable  :: benergy
 double precision, dimension(:), allocatable  :: bx
 double precision, dimension(:), allocatable  :: bx2
 double precision, dimension(:), allocatable  :: bnbis_avg
 double precision, dimension(:), allocatable  :: bnox
 double precision, dimension(:), allocatable  :: bnoy
 double precision, dimension(:), allocatable  :: bntot
 double precision, dimension(:), allocatable  :: bnpa1g_avg
 double precision, dimension(:), allocatable  :: bsp_site_avg
 double precision, dimension(:), allocatable  :: bbp_site_avg

 !sublattice dependent quantities
 double precision, dimension(:), allocatable  :: bnbi1
 double precision, dimension(:), allocatable  :: bnbi2
 double precision, dimension(:), allocatable  :: bnpx1
 double precision, dimension(:), allocatable  :: bnpx2
 double precision, dimension(:), allocatable  :: bnpy1
 double precision, dimension(:), allocatable  :: bnpy2
 double precision, dimension(:), allocatable  :: bnpA1g1
 double precision, dimension(:), allocatable  :: bnpA1g2
 double precision, dimension(:), allocatable  :: bsp1
 double precision, dimension(:), allocatable  :: bsp2
 double precision, dimension(:), allocatable  :: bbp1
 double precision, dimension(:), allocatable  :: bbp2

 double precision, dimension(:,:), allocatable  :: bnbis
 double precision, dimension(:,:), allocatable  :: bnpa1g
 double precision, dimension(:,:), allocatable :: bspolaron
 double precision, dimension(:,:), allocatable :: bbpolaron
 double precision, dimension(:,:), allocatable :: bspolaron_ij
 double precision, dimension(:,:), allocatable  :: bchi_sc
 double precision, dimension(:,:,:,:), allocatable  :: bchi_ch

contains
 !=============================================================================
 subroutine get_err(bins,mean,std)
 use parameters, only: nbin
 implicit none
 integer i
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
 allocate(anbis(0:Nbi-1))
 allocate(anpa1g(0:Nbi-1))
 allocate(aspolaron(0:Nbi-1))
 allocate(abpolaron(0:Nbi-1))
 allocate(aspolaron_ij(0:nclass-1))
 allocate(achi_sc(0:Norb-1))
 allocate(achi_ch(0:nclass-1, 0:Norb-1, 0:Norb-1))
 allocate(benergy(1:nbin))
 allocate(bx(nbin))
 allocate(bx2(nbin))
 allocate(bnbis_avg(nbin))
 allocate(bnox(nbin))
 allocate(bnoy(nbin))
 allocate(bntot(nbin))
 allocate(bnpa1g_avg(nbin))
 allocate(bnbi1(nbin))
 allocate(bnbi2(nbin))
 allocate(bnpx1(nbin))
 allocate(bnpx2(nbin))
 allocate(bnpy1(nbin))
 allocate(bnpy2(nbin))
 allocate(bnpA1g1(nbin))
 allocate(bnpA1g2(nbin))
 allocate(bsp1(nbin))
 allocate(bsp2(nbin))
 allocate(bbp1(nbin))
 allocate(bbp2(nbin))
 allocate(bnbis(0:Nbi-1, nbin))
 allocate(bnpa1g(0:Nbi-1, nbin))
 allocate(ns_mc(0:Nbi-1))
 allocate(npx_mc(0:Nbi-1))
 allocate(nLs_mc(0:Nbi-1)) 
 allocate(Sp_mc(0:Nbi-1))
 allocate(Bp_mc(0:Nbi-1))
 allocate(bsp_site_avg(nbin))
 allocate(bbp_site_avg(nbin))
 allocate(bspolaron(0:Nbi-1, nbin))
 allocate(bbpolaron(0:Nbi-1, nbin))
 allocate(bspolaron_ij(0:nclass-1, nbin))
 allocate(bchi_sc(0:Norb-1, nbin))
 allocate(bchi_ch(0:nclass-1, 0:Norb-1, 0:Norb-1, nbin))
end subroutine allocate_quantities
 !=============================================================================
 subroutine deallocate_quantities()
 implicit none
 deallocate(anbis)
 deallocate(anpa1g)
 deallocate(aspolaron)
 deallocate(abpolaron)
 deallocate(aspolaron_ij)
 deallocate(achi_sc)
 deallocate(achi_ch)
 deallocate(bnbis_avg)
 deallocate(bnpa1g_avg)
 deallocate(bnbis)
 deallocate(bnpa1g)
 deallocate(bsp_site_avg)
 deallocate(bbp_site_avg)
 deallocate(bnbi1)
 deallocate(bnbi2)
 deallocate(bnpx1)
 deallocate(bnpx2)
 deallocate(bnpy1)
 deallocate(bnpy2)
 deallocate(bnpA1g1)
 deallocate(bnpA1g2)
 deallocate(bsp1)
 deallocate(bsp2)
 deallocate(bbp1)
 deallocate(bbp2)
 deallocate(ns_mc)
 deallocate(npx_mc)
 deallocate(nLs_mc)
 deallocate(Sp_mc)
 deallocate(Bp_mc)
 deallocate(bspolaron)
 deallocate(bbpolaron)
 deallocate(bspolaron_ij)
 deallocate(bchi_sc)
 deallocate(bchi_ch)
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
 anbis = 0.0d0
 antot = 0.0d0
 anpa1g = 0.0d0
 anbis_avg = 0.0d0
 anpa1g_avg = 0.0d0
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
 anpA1g1 = 0.0d0
 anpA1g2 = 0.0d0
 asp1 = 0.0d0
 asp2 = 0.0d0
 abp1 = 0.0d0
 abp2 = 0.0d0
 aspolaron_ij = 0.0d0
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
 bnpa1g_avg(bin) = anpa1g_avg/dfloat(cnt)
 bnbis(:,bin) = anbis/dfloat(cnt)
 bnpa1g(:,bin) = anpa1g/dfloat(cnt)
 bsp_site_avg(bin) = asp_site_avg/dfloat(cnt)
 bbp_site_avg(bin) = abp_site_avg/dfloat(cnt)

 bnbi1(bin) = anbi1/dfloat(cnt)
 bnbi2(bin) = anbi2/dfloat(cnt)
 bnpx1(bin) = anpx1/dfloat(cnt)
 bnpx2(bin) = anpx2/dfloat(cnt)
 bnpy1(bin) = anpy1/dfloat(cnt)
 bnpy2(bin) = anpy2/dfloat(cnt)
 bnpA1g1(bin) = anpA1g1/dfloat(cnt)
 bnpA1g2(bin) = anpA1g2/dfloat(cnt)
 bsp1(bin) = asp1/dfloat(cnt)
 bsp2(bin) = asp2/dfloat(cnt)
 bbp1(bin) = abp1/dfloat(cnt)
 bbp2(bin) = abp2/dfloat(cnt)

 bspolaron(:,bin) = aspolaron/dfloat(cnt)
 bbpolaron(:,bin) = abpolaron/dfloat(cnt)
 bspolaron_ij(:,bin) = aspolaron_ij/dfloat(cnt)
 bchi_sc(:,bin) = achi_sc/dfloat(cnt)
 bchi_ch(:,:,:,bin) = achi_ch/dfloat(cnt)
 return
 end subroutine populate_bins

 !============================================================================
 ! do_measurements below
 !============================================================================
 subroutine do_measurements(X)
 use parameters
 use cluster, only: return_index_for_coordinates, dclass, dclass_F, expqr
 use monte_carlo, only: compute_total_E, get_H
 implicit none
 integer ix, iy, jx, jy, tau, nn, nnp, mm, nnn
 integer ixp, iyp, ixm, iym
 integer jxp, jyp, jxm, jym
 integer i, j, k, info, i1, i2, j1, j2
 integer lwork
 double precision, dimension(0:N-1,0:N-1) :: U
 double precision, allocatable, dimension(:) :: work
 double precision Xi_A1g, Xj_A1g
 double precision a_inn, a_inn2, a_jnn, a_jnn2
 double precision a_innp, a_innp2, a_jnnp, a_jnnp2
 double precision s_inn, s_inn2, s_innp, s_innp2
 double precision s_jnn, s_jnn2, s_jnnp, s_jnnp2 
 double precision tmp1, tmp2, tmp3
 double precision energy, fermi, fac, fermi1, fac1, factor
 double precision, dimension(0:N-1) :: X
 double precision, dimension(0:N-1) :: Ek

 !Zero the record for each MC measurement
 ns_mc = 0.0d0
 npx_mc = 0.0d0
 nLs_mc = 0.0d0
 Sp_mc = 0.0d0
 Bp_mc = 0.0d0

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

 !loop over Bi sites
 do ix = 0,Nx-1
  do iy = 0,Ny-1
   i = return_index_for_coordinates(ix,iy,0)
   ixp = return_index_for_coordinates(ix,iy,1) 
   iyp = return_index_for_coordinates(ix,iy,2) 
   ixm = return_index_for_coordinates(ix-1,iy,1)
   iym = return_index_for_coordinates(ix,iy-1,2)
   Xi_A1g = -0.5d0*(X(ixp) - X(ixm) + X(iyp) - X(iym))

   !sum over eigenstates
   do nn = 0,N-1  
     fermi = 1.0d0/(exp(beta*(Ek(nn)-mu))+1.0d0)
     fac = 2.0d0*fermi/Nbi
     tmp1 = fac*U(i,nn)*U(i,nn)
     tmp2 = fac*U(ixp,nn)*U(ixp,nn)
     tmp3 = fac*U(iyp,nn)*U(iyp,nn)
     anbis(i) = anbis(i) + tmp1*Nbi
     anbis_avg = anbis_avg + tmp1
     anox  = anox  + tmp2
     anoy  = anoy  + tmp3
     antot = antot + tmp1+tmp2+tmp3

     a_inn = 0.5d0*(U(ixp,nn) - U(ixm,nn) + U(iyp,nn) - U(iym,nn))
     s_inn = U(i,nn)
     a_inn2 = a_inn*a_inn
     s_inn2 = s_inn*s_inn

     anpa1g(i)  = anpa1g(i) + fac*a_inn2*Nbi
     anpa1g_avg = anpa1g_avg + fac*a_inn2

     !record quantities for each MC measurement
     if (if_print_MC==1) then
       ns_mc(i) = ns_mc(i) + tmp1*Nbi
       npx_mc(i) = npx_mc(i) + tmp2*Nbi
       nLs_mc(i) = nLs_mc(i) + fac*a_inn2*Nbi
     endif

     !Also compute orbital occupancy for two sublattices
     ! /0.5 below accounts for Nbi/2 sites of sublattice
     if (mod(ix+iy,2)==0) then
       anbi1 = anbi1 + tmp1/0.5
       anpx1 = anpx1 + tmp2/0.5
       anpy1 = anpy1 + tmp3/0.5
       anpA1g1 = anpA1g1 + fac*a_inn2/0.5
     else
       anbi2 = anbi2 + tmp1/0.5
       anpx2 = anpx2 + tmp2/0.5
       anpy2 = anpy2 + tmp3/0.5
       anpA1g2 = anpA1g2 + fac*a_inn2/0.5
     endif

     ! aspolaron does not need Nbi because of (i) index
     tmp1 = 2.0d0*fermi*Xi_A1g*(s_inn2 + a_inn2)
    ! aspolaron(i) = aspolaron(i) + Xi_A1g*tmp1
     aspolaron(i) = aspolaron(i) + tmp1
     asp_site_avg = asp_site_avg + tmp1/Nbi

     !record quantities for each MC measurement
     if (if_print_MC==1) then
       Sp_mc(i) = Sp_mc(i) + tmp1
     endif

     !Also compute single polaron density for two sublattices
     ! /0.5 below accounts for Nbi/2 sites of sublattice
     if (mod(ix+iy,2)==0) then
       asp1 = asp1 + tmp1/(Nbi/2.0)
     else
       asp2 = asp2 + tmp1/(Nbi/2.0)
     endif

     ! second sum over eigenstates for bipolaron etc.
     do nnp = 0,N-1     
       !!!!!!!!!!!!!!!!!!!!!!!!!!
       !!       bipolaron      !!
       !!!!!!!!!!!!!!!!!!!!!!!!!!
       fermi1 = 1.0d0/(exp(beta*(Ek(nnp)-mu))+1.0d0)
       fac1 = 2.0d0*fermi1/Nbi

       a_innp = 0.5d0*(U(ixp,nnp) - U(ixm,nnp) + U(iyp,nnp) - U(iym,nnp))
       s_innp = U(i,nnp)
       a_innp2 = a_innp*a_innp
       s_innp2 = s_innp*s_innp
       
       tmp1 = fermi*fermi1* Xi_A1g*   &
              (s_inn2 + a_inn2)*(s_innp2 + a_innp2)
       abpolaron(i) = abpolaron(i) + tmp1
       abp_site_avg = abp_site_avg + tmp1/Nbi

       !record quantities for each MC measurement
       if (if_print_MC==1) then
         Bp_mc(i) = Bp_mc(i) + tmp1
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

           !!!!!!!!!!!!!!!!!!!!!!!!
           !!       <L_i*L_j>    !!
           !!!!!!!!!!!!!!!!!!!!!!!!
           Xj_A1g = -0.5d0*(X(jxp) - X(jxm) + X(jyp) - X(jym))
           a_jnnp = 0.5d0*(U(jxp,nnp) - U(jxm,nnp) + U(jyp,nnp) - U(jym,nnp))
           a_jnn  = 0.5d0*(U(jxp,nn) - U(jxm,nn) + U(jyp,nn) - U(jym,nn))
           a_jnn2  = a_jnn *a_jnn
           a_jnnp2 = a_jnnp*a_jnnp

           factor = fac*fac1*Xi_A1g*Xj_A1g/dclass_F(k)
           aspolaron_ij(k) = aspolaron_ij(k) + factor*(          &
                                      s_inn2*s_jnnp2             &
                                    - s_inn*s_jnn*s_innp*s_jnnp  &
                                    + s_inn2*a_jnn2              &
                                    - s_inn*a_jnn*s_innp*a_jnnp  &
                                    + a_inn2*s_jnn2              &
                                    - a_inn*s_jnn*a_innp*s_jnnp  &
                                    + a_inn2*a_jnnp2             &
                                    - a_inn*a_jnn*a_innp*a_jnnp )

           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                                                 
           !!      s-wave susceptibilities     !!                                                 
           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
           do o1 = 0,Norb-1
             i1 = return_index_for_coordinates(ix,iy,o1)
             i2 = return_index_for_coordinates(jx,jy,o1)
             achi_sc(o1) = achi_sc(o1) + fermi*fermi1/Nbi   &
                 *U(i1,nn)*U(i1,nn)*U(i2,nnp)*U(i2,nnp)
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
               achi_ch(k,o1,o2) = achi_ch(k,o1,o2) + 4.0d0*fermi*fermi1/Nbi   &
                   *expqr(k,j1,j2)*U(i1,nn)*U(i1,nn)*U(i2,nnp)*U(i2,nnp)
                  ! *U(i1,nn)*U(i1,nn)*U(i2,nnp)*U(i2,nnp)
             enddo
           enddo
         enddo
       enddo  !end loop jx, jy
     enddo !end loop nnp
    enddo  !end loop nn
    
    !Print quantities for each MC measurement
    if (if_print_MC==1) then
      530 format(f13.5)
      write(unit=13,fmt=530)  Xi_A1g
      write(unit=14,fmt=530)  ns_mc(i)
      write(unit=15,fmt=530)  npx_mc(i)
      write(unit=16,fmt=530)  nLs_mc(i)
      write(unit=17,fmt=530)  Sp_mc(i)/Xi_A1g
      write(unit=18,fmt=530)  Bp_mc(i)/Xi_A1g
      write(unit=19,fmt=530)  Sp_mc(i)
      write(unit=20,fmt=530)  Bp_mc(i)
    endif

   enddo
 enddo   !end loop ix, iy

 deallocate(work)

 return
 end subroutine do_measurements


end module measurements
