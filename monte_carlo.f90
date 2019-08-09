module monte_carlo
 integer accept, reject
contains
 !=============================================================================
 subroutine single_site_sweep(X,accept,reject)
 use parameters
 use random
 implicit none
 integer accept, reject
 integer site
 double precision Eold, Enew
 double precision deltaX
 double precision r, ratio
 double precision, dimension(0:N-1) :: X, Xproposed

 !sweep every site
 call compute_total_E(Eold,X)
 do site = NBi,N-1
  if (if_X_displace==1) then
    !deltaX = (ran2(iran)-0.5d0)*2.d0*dXamp
    deltaX = 0.d0
  else
    deltaX = 0.d0
  endif

  Xproposed = X
  Xproposed(site) = Xproposed(site) + deltaX

  call compute_total_E(Enew,XProposed)
  ratio = minval([exp(-beta*(Enew-Eold)),1.0d0])
  r = ran2(iran)
  !print *, "Eold = ", Eold
  !print *, "Enew = ", Enew
  !print *, "ratio = ", ratio

  if(r.le.ratio)then
   accept = accept + 1
   Eold = Enew
   X = Xproposed
   !print*, "accepted new displacement X="
   !print*, X(NBi:N-1)
  else
   !print*, "X update rejected"
   reject = reject + 1
  endif
 enddo

 return
 end subroutine single_site_sweep
 !=============================================================================
 subroutine get_H(H,X)
 use parameters, only: Norb, Nx, Ny, N, Nbi, tsp, tpp, es, ep, d0, if_X_displace
 use cluster, only: return_index_for_coordinates
 implicit none
 integer ix, iy, del, idx, i, ixp, ixm, iyp, iym
 integer ixpy, iymy, ixpymy, iypx, ixmx, ixmypx
 double precision, dimension(0:N-1) :: X
 double precision, dimension(0:N-1,0:N-1) :: H

 H = 0.0d0
 do ix = 0,Nx-1
  do iy = 0,Ny-1
   !i is a Bi index
   i = return_index_for_coordinates(ix,iy,0)
   !ixp is the index for the oxygen to the right of the Bi atom
   ixp = return_index_for_coordinates(ix,iy,1)
   iyp = return_index_for_coordinates(ix,iy,2)
   ixm = return_index_for_coordinates(ix-1,iy,1)
   iym = return_index_for_coordinates(ix,iy-1,2)
   !set the onsite energies
   H(i,i) = es
   H(ixp,ixp) = ep
   H(iyp,iyp) = ep
   H(i,ixp) =-tsp/((1.0d0 + if_X_displace* X(ixp)/d0)**(2.0d0))
   H(i,iyp) =-tsp/((1.0d0 + if_X_displace* X(iyp)/d0)**(2.0d0))
   H(i,ixm) = tsp/((1.0d0 - if_X_displace* X(ixm)/d0)**(2.0d0))
   H(i,iym) = tsp/((1.0d0 - if_X_displace* X(iym)/d0)**(2.0d0))
   H(ixp,i) = H(i,ixp)
   H(iyp,i) = H(i,iyp)
   H(ixm,i) = H(i,ixm)
   H(iym,i) = H(i,iym)

   ! set tpp for px orbital:
   ixpy = return_index_for_coordinates(ix+1,iy,2)
   iymy = return_index_for_coordinates(ix,iy-1,2)
   ixpymy = return_index_for_coordinates(ix+1,iy-1,2)
   H(ixp,ixpy)   =-tpp
   H(ixp,iyp)    = tpp
   H(ixp,iymy)   =-tpp
   H(ixp,ixpymy) = tpp
   H(ixpy,ixp)   =-tpp
   H(iyp, ixp)   = tpp
   H(iymy,ixp)   =-tpp
   H(ixpymy,ixp) = tpp
   ! set tpp for py orbital:
   iypx = return_index_for_coordinates(ix,iy+1,1)
   ixmx = return_index_for_coordinates(ix-1,iy,1)
   ixmypx = return_index_for_coordinates(ix-1,iy+1,1)
   H(iyp,iypx)   =-tpp
   H(iyp,ixmypx) = tpp
   H(iyp,ixmx)   =-tpp
   H(iyp,ixp)    = tpp
   H(iypx,  iyp) =-tpp
   H(ixmypx,iyp) = tpp
   H(ixmx,  iyp) =-tpp
   H(ixp,   iyp) = tpp
  enddo
 enddo

 do ix = 0,N-1
  do iy = ix,N-1
   if(H(ix,iy).ne.H(iy,ix))then
    print*, 'Error, H is not hermitian.'
   endif
  enddo
 enddo
 return
 end subroutine get_H
 !=============================================================================

 !=============================================================================
 double precision function elastic_energy(X)
 use parameters, only: Nbi, N, spring_Const, alpha
 implicit none
 double precision, dimension(0:N-1) :: X
 elastic_energy = 0.5d0*Spring_Const*sum(X*X) + 0.25d0*alpha*sum(X*X*X*X) 
 return
 end function elastic_energy
 !=============================================================================

 !=============================================================================
 subroutine compute_total_E(energy,X)
 use parameters, only: N, beta, mu
 implicit none
 integer info
 integer lwork
 integer i,j
 double precision energy 
 double precision fermi
 double precision, dimension(0:N-1) :: X
 double precision, dimension(0:N-1,0:N-1) :: H0
 double precision, dimension(0:N-1) :: Ek
 double precision, allocatable, dimension(:) :: Work

 lwork = 3*N-1
 allocate(work(1:lwork))

 !set the total energy to the elastic energy
 energy = elastic_energy(X)
 !print*, '-----------------------------------'
 !print *, 'elastic energy=',energy
 !Get the Hamiltonian for the current displacement pattern
 info = 0

 call get_H(H0,X)
 !get the eigenvalues for H0

 call dsyev('N','U',N,H0,N,Ek,work,lwork,info)
 !print *, "X=",X
 !print *, "sum(X)=", sum(X)/N
 !print *, "min and max eigen E=", minval(Ek), maxval(Ek)
 if (info.ne.0) then
    print *,"dsyev error", N, lwork, info
    stop
 endif

 !Get the total energy
 do i = 0,N-1
   fermi = 1.0d0/(exp(beta*(Ek(i)-mu))+1.0d0)
   energy = energy + 2.0d0*(Ek(i)-mu)*fermi
 enddo
 !print *, 'total E=',energy
 return
 end subroutine compute_total_E

end module monte_carlo
