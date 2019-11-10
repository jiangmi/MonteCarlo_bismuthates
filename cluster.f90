! The order of N orbitals is:
! first 0:NBi-1 is bismuth and then px orbital and finally py orbital 
! N = Nx*Ny*Norb; NBi = Nx*Ny
module cluster 
 use parameters, only: N, pi, Nx, Ny, Nbi, nclass
 real*8 qx, qy
 integer, dimension(:,:), allocatable :: dclass  
 integer, dimension(:),   allocatable :: dclass_F
 real*8,  dimension(:,:,:), allocatable :: expqr
 real*8,  dimension(:,:),   allocatable :: phase
contains
 
 !=============================================================================
 ! (i,j) label which unit cell, del denotes site/orbital index 
 integer function return_index_for_coordinates(i,j,d)
 implicit none
 integer ix, iy, del, i, j, d
 ix = i
 iy = j
 del = d
 if(Del.ge.3)then
  print*, 'Invalid basis site index passed to return_index_for_coordinates.'
  stop
 endif
 if(ix.lt.0) ix = ix + Nx
 if(ix.ge.Nx) ix = ix - Nx
 if(iy.lt.0) iy = iy + Ny
 if(iy.ge.Ny) iy = iy - Ny
 return_index_for_coordinates = ix + iy*Nx + del*Nx*Ny
 return
 end function return_index_for_coordinates

 !=============================================================================
 subroutine return_coordinates_for_index(idx,ix,iy,del)
 implicit none
 integer itmp, idx, ix, iy, del

 del = idx/(Nbi)
 itmp = idx - del*Nbi
 iy = itmp/Nx
 ix = itmp - iy*Nx
 return
 end subroutine return_coordinates_for_index

 !=============================================================================
 logical function is_bismuth(i)
 implicit none
 integer i
 if(i.le.Nbi-1)then
   is_bismuth = .true.
 else
   is_bismuth = .false.
 endif
 return
 end function is_bismuth

 !=============================================================================
 subroutine get_distance_class()
 ! Obtain the distance_class for two unit cell origins (Bi atoms)
 ! Also get table of exp(q*(r_i+r_j)) for computing charge chi_ch(q)
 ! In expqr(q, r_i, r_j), the number of possible q is same as distance class
 ! but q's value is (dx*2*pi/Nx, dy*2*pi/Ny)
 implicit none
 integer dx, dy, ix,iy,jx,jy,i,j,k,cnt
 allocate(dclass(0:Nbi-1,0:Nbi-1))
 allocate(dclass_F(0:nclass-1))
 allocate(expqr(0:nclass-1, 0:Nbi-1,0:Nbi-1))
 allocate(phase(0:Nbi-1,0:Nbi-1))
 expqr = 0.d0
 dclass = 0
 dclass_F = 0
 cnt = 0

 do ix = 0,Nx-1
  do iy = 0,Ny-1
   i = return_index_for_coordinates(ix,iy,0)
   do jx = 0,Nx-1
    do jy = 0,Ny-1
     j = return_index_for_coordinates(jx,jy,0)

     dx = abs(ix-jx)
     if (dx>Nx/2) then
       dx = Nx-dx
     endif
     dy = abs(iy-jy)
     if (dy>Nx/2) then
       dy = Nx-dy
     endif

     ! set additional symmetry (dx,dy)=(dy,dx)
     if (dx<dy) then
       k = dx
       dx = dy
       dy = k
     endif

     k = get_index(dx,dy)
     dclass(i,j) = k
     !print*, k, 'd=',dx, dy, 'site',i,j

     ! Note for k = get_index(dx,dy) with dx<dy, dclass_F(k)=0
     dclass_F(k) = dclass_F(k)+1 

     ! Compute exp(i*q*(r_i+r_j))     
     qx = dx*2.0d0*pi/Nx
     qy = dy*2.0d0*pi/Ny
     expqr(k,i,j) = cos(qx*(ix+jx) + qy*(iy+jy))
   !  print*, 'q=',qx,qy,'r=',ix+jx,iy+jy,'expqr=',expqr(k,i,j)

     ! staggered phase for computing staggered correlation in measurement.f90
     if (mod(dx+dy,2)==0) then
       phase(i,j) = 1.0
     else
       phase(i,j) = -1.0
     endif
    enddo
   enddo
  enddo
 enddo
 return
 end subroutine get_distance_class

 !=============================================================================                       
 ! get a unique integer to identify a distance class      
 ! (i,j) labels distance = (dx,dy) between sites                                 
 integer function get_index(i,j)                                                 
 implicit none 
 integer i, j
 if(i>Nx/2 .or. j>Ny/2)then                                                                                     
  print*, 'Distance vector between two sites exceeds half lattice size!'                          
  stop                                                                                                
 endif     
 get_index = i + j*(Nx/2+1)                                             
 return    
 end function get_index

end module cluster 
