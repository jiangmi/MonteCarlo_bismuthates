module measurements
 use parameters, only: nbin
 integer cnt  ! count the measurements in each bin

 ! accumulated quantites
 double precision aenergy, ax, ax2
 double precision antot, anbis, anox, anoy
 double precision anpa1g  !Number of oxygen holes with A1g symmetry
 double precision asp_site_avg
 double precision abp_site_avg
 double precision, dimension(:), allocatable  :: aspolaron
 double precision, dimension(:), allocatable  :: abpolaron
 double precision, dimension(:,:), allocatable :: aspolaron_ij

 !Below are for binned quantites
 double precision, dimension(:), allocatable  :: benergy 
 double precision, dimension(:), allocatable  :: bx
 double precision, dimension(:), allocatable  :: bx2
 double precision, dimension(:), allocatable  :: bnbis
 double precision, dimension(:), allocatable  :: bnox
 double precision, dimension(:), allocatable  :: bnoy
 double precision, dimension(:), allocatable  :: bntot
 double precision, dimension(:), allocatable  :: bnpa1g
 double precision, dimension(:), allocatable  :: bsp_site_avg
 double precision, dimension(:), allocatable  :: bbp_site_avg
 double precision, dimension(:,:), allocatable :: bspolaron
 double precision, dimension(:,:), allocatable :: bbpolaron
 double precision, dimension(:,:,:), allocatable :: bspolaron_ij
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
 allocate(aspolaron(0:Nbi-1))
 allocate(abpolaron(0:Nbi-1))
 allocate(aspolaron_ij(0:Nbi-1,0:Nbi-1))
 allocate(benergy(1:nbin))
 allocate(bx(1:nbin))
 allocate(bx2(1:nbin))
 allocate(bnbis(1:nbin))
 allocate(bnox(1:nbin))
 allocate(bnoy(1:nbin))
 allocate(bntot(1:nbin))
 allocate(bnpa1g(1:nbin))
 allocate(bsp_site_avg(1:nbin))
 allocate(bbp_site_avg(1:nbin))
 allocate(bspolaron(1:nbin, 0:Nbi-1))
 allocate(bbpolaron(1:nbin, 0:Nbi-1))
 allocate(bspolaron_ij(1:nbin, 0:Nbi-1, 0:Nbi-1))
 end subroutine allocate_quantities
 !=============================================================================
 subroutine deallocate_quantities()
 implicit none
 deallocate(aspolaron)
 deallocate(abpolaron)
 deallocate(aspolaron_ij)
 deallocate(bsp_site_avg)
 deallocate(bbp_site_avg)
 deallocate(bspolaron)
 deallocate(bbpolaron)
 deallocate(bspolaron_ij)
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
 asp_site_avg = 0.d0
 abp_site_avg = 0.d0
 aspolaron = 0.0d0
 abpolaron = 0.0d0 
 aspolaron_ij = 0.0d0
 end subroutine zero_accumulators
 !=============================================================================
 subroutine populate_bins(bin)
 implicit none
 integer bin
 print*, 'Populating bin ', bin
 benergy(bin) = aenergy/dfloat(cnt)
 bX(bin) = aX/dfloat(cnt)
 bX2(bin) = aX2/dfloat(cnt)
 bntot(bin) = antot/dfloat(cnt)
 bnbis(bin) = anbis/dfloat(cnt)
 bnox(bin) = anox/dfloat(cnt)
 bnoy(bin) = anoy/dfloat(cnt)
 bnpa1g(bin) = anpa1g/dfloat(cnt)
 bsp_site_avg(bin) = asp_site_avg/dfloat(cnt)
 bbp_site_avg(bin) = abp_site_avg/dfloat(cnt)
 bspolaron(bin,:) = aspolaron/dfloat(cnt)
 bbpolaron(bin,:) = abpolaron/dfloat(cnt)
 bspolaron_ij(bin,:,:) = aspolaron_ij/dfloat(cnt)
 return
 end subroutine populate_bins
 
!============================================================================
 ! do_measurements below
 !============================================================================
 subroutine do_measurements(X)
 use parameters
 use cluster, only: return_index_for_coordinates
 use monte_carlo, only: compute_total_E, get_H
 implicit none
 integer ix, iy, jx, jy, tau, nn, nnp, mm, nnn
 integer ixp, iyp, ixm, iym
 integer jxp, jyp, jxm, jym
 integer i, j, info
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
   Xi_A1g = 0.5d0*(X(ixp) - X(ixm) + X(iyp) - X(iym))

   !sum over eigenstates
   do nn = 0,N-1       
     fermi = 1.0d0/(exp(beta*(Ek(nn)-mu))+1.0d0)
     fac = 2.0d0*fermi/Nbi
     tmp1 = fac*U(i,nn)*U(i,nn)
     tmp2 = fac*U(ixp,nn)*U(ixp,nn)
     tmp3 = fac*U(iyp,nn)*U(iyp,nn)
     anbis = anbis + tmp1
     anox  = anox  + tmp2
     anoy  = anoy  + tmp3
     antot = antot + tmp1+tmp2+tmp3

     a_inn = 0.5d0*(U(ixp,nn) - U(ixm,nn) + U(iyp,nn) - U(iym,nn))
     s_inn = U(i,nn)
     a_inn2 = a_inn*a_inn
     s_inn2 = s_inn*s_inn

     anpa1g = anpa1g + fac*a_inn2

     ! aspolaron does not need Nbi because of (i) index
     tmp1 = 2.0d0*fermi*(s_inn2 + a_inn2)
    ! aspolaron(i) = aspolaron(i) + Xi_A1g*tmp1
     aspolaron(i) = aspolaron(i) + tmp1
     asp_site_avg = asp_site_avg + tmp1/Nbi

     ! second sum over eigenstates for bipolaron etc.
     do nnp = 0,N-1     
       fermi1 = 1.0d0/(exp(beta*(Ek(nnp)-mu))+1.0d0)
       fac1 = 2.0d0*fermi1/Nbi

       a_innp = 0.5d0*(U(ixp,nnp) - U(ixm,nnp) + U(iyp,nnp) - U(iym,nnp))
       s_innp = U(i,nnp)
       a_innp2 = a_innp*a_innp
       s_innp2 = s_innp*s_innp
       
       tmp1 = 2.0d0*fermi**fermi1* &!Xi_A1g*   &
              (a_inn2*s_innp2 + s_inn2*s_innp2 + a_inn2*a_innp2)
       abpolaron(i) = abpolaron(i) + tmp1
       abp_site_avg = abp_site_avg + tmp1/Nbi

       !!!!!!!!!!!!!!!!!!!!!!!!!!
       !!       <L_i*L_j>      !!
       !!!!!!!!!!!!!!!!!!!!!!!!!!
       do jx = 0,Nx-1
         do jy = 0,Ny-1
           j = return_index_for_coordinates(jx,jy,0)
           jxp = return_index_for_coordinates(jx  ,jy  ,1)
           jxm = return_index_for_coordinates(jx-1,jy  ,1)
           jyp = return_index_for_coordinates(jx  ,jy  ,1)
           jym = return_index_for_coordinates(jx  ,jy-1,1)

           ! TODO: classify the vector j-i 

           Xj_A1g = 0.5d0*(X(jxp) - X(jxm) + X(jyp) - X(jym))
           a_jnnp = 0.5d0*(U(jxp,nnp) - U(jxm,nnp) + U(jyp,nnp) - U(jym,nnp))
           a_jnn  = 0.5d0*(U(jxp,nn) - U(jxm,nn) + U(jyp,nn) - U(jym,nn))
           a_jnn2  = a_jnn *a_jnn
           a_jnnp2 = a_jnnp*a_jnnp

           factor = fac*fac1*Xi_A1g*Xj_A1g
           aspolaron_ij(i,j) = aspolaron_ij(i,j) + factor*(          &
                                          s_inn2*s_jnnp2             &
                                        - s_inn*s_jnn*s_innp*s_jnnp  &
                                        + s_inn2*a_jnn2              &
                                        - s_inn*a_jnn*s_innp*a_jnnp  &
                                        + a_inn2*s_jnn2              &
                                        - a_inn*s_jnn*a_innp*s_jnnp  &
                                        + a_inn2*a_jnnp2             &
                                        - a_inn*a_jnn*a_innp*a_jnnp )
         enddo
       enddo  !end loop jx, jy
     enddo !end loop nnp
    enddo  !end loop nn
   enddo
 enddo   !end loop ix, iy

 deallocate(work)

 return
 end subroutine do_measurements


end module measurements
