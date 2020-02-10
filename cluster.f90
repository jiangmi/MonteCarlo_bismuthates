! The order of N orbitals is:
! first 0:NBi-1 is bismuth and then px orbital and finally py orbital 
! N = Nx*Ny*Norb; NBi = Nx*Ny
module cluster 
 use parameters, only: N, Nc, pi, Nx, Ny, Ncx, Ncy, Nbi, nclass
 real*8 qx, qy
 integer, dimension(:,:), allocatable :: dclass  
 integer, dimension(:),   allocatable :: dclass_F
 real*8,  dimension(:,:,:), allocatable :: expqr
 real*8,  dimension(:,:),   allocatable :: phase
 integer, dimension(:,:),   allocatable :: dx
 integer, dimension(:,:),   allocatable :: dy
 integer, dimension(:,:),   allocatable :: cluster_boundary
 integer, dimension(:,:),   allocatable :: cluster_sites
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
 ! similat to return_index_for_coordinates above but replace Nx by Ncx
 integer function return_index_for_coordinates_cluster(i,j,d)
 implicit none
 integer ix, iy, del, i, j, d
 ix = i
 iy = j
 del = d
 if(Del.ge.3)then
  print*, 'Invalid basis site index passed to return_index_for_coordinates.'
  stop
 endif
 if(ix.lt.0) ix = ix + Ncx
 if(ix.ge.Ncx) ix = ix - Ncx
 if(iy.lt.0) iy = iy + Ncy
 if(iy.ge.Ncy) iy = iy - Ncy
 return_index_for_coordinates_cluster = ix + iy*Ncx + del*Ncx*Ncy
 return
 end function return_index_for_coordinates_cluster

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
 ! correct index when out of lattice boundary using periodic boundary
 integer function pos(i,Ns)
 implicit none
 integer i, j, Ns
 j = i
 if(i.lt.0)  j = i + Ns
 if(i.ge.Ns) j = i - Ns
 pos = j
 return
 end function pos

 !=============================================================================
 !sort index array in ascending order
 subroutine sort(indices)
 implicit none
 integer, dimension(1:Nc) :: indices  
 integer i, k, buf
 do i= 1,Nc-1
    k = minloc(indices(i:Nc), dim=1) + i-1
    buf = indices(i)
    indices(i) = indices(k)
    indices(k) = buf
 enddo

 return
 end subroutine sort

 !=============================================================================
 subroutine find_cluster()
 !find unit cell range of the cluster around all px/py site i despite 
 !declared cluster_boundary as (0:N-1,:) but all Bi sites are set to be zero
 !px/py orbital is in the center of the cluster
 !cluster_boundary(i,4): 4 denotes the cluster located within x:[1,2]; y:[3,4]
 implicit none
 integer ix, iy, jx, jy, kx, ky, i, j, cnt, iorb, jorb
 integer, dimension(1:Nc) :: tmp  
 allocate(cluster_boundary(0:N-1, 4))
 allocate(cluster_sites(0:N-1,0:Nc-1))

 cluster_boundary = 0
 cluster_sites = 0

 do ix = 0,Nx-1
  do iy = 0,Ny-1
   do iorb = 1,2
     i = return_index_for_coordinates(ix,iy,iorb)
     cluster_boundary(i, 1) = pos(ix-(Ncx/2-1), Nx)
     cluster_boundary(i, 2) = pos(ix+Ncx/2, Nx)
     cluster_boundary(i, 3) = pos(iy-(Ncy/2-1), Ny)
     cluster_boundary(i, 4) = pos(iy+Ncy/2, Ny)
     !print*, i,'bound:', cluster_boundary(i, 1:4)

     cnt = 0
     tmp = 0

     !Follow the ordering rule of original lattice (see top):
     !first 0:NBi-1 is bismuth and then px and finally py
     do jorb = 0,2
       do jy = iy-(Ncy/2-1), iy+Ncy/2
         ky = pos(jy,Ny)

         do jx = ix-(Ncx/2-1), ix+Ncx/2
           kx = pos(jx,Nx)

           j = return_index_for_coordinates(kx,ky,jorb)
           tmp(cnt) = j
           cnt = cnt+1
         enddo
       enddo
     enddo
     
     if (cnt/=Nc) then
       print*, 'Error: # of cluster sites incorrect!'
     endif

     !print*, 'site',i, 'cluster:', tmp
     !call sort(tmp)
     cluster_sites(i,:) = tmp
     !print*, 'site',i, 'sorted cluster:', cluster_sites(i,:)
   enddo
  enddo
 enddo

 return
 end subroutine find_cluster

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
 integer ix,iy,jx,jy,i,j,k,cnt
 allocate(dclass(0:Nbi-1,0:Nbi-1))
 allocate(dclass_F(0:nclass-1))
 allocate(expqr(0:nclass-1, 0:Nbi-1,0:Nbi-1))
 allocate(phase(0:Nbi-1,0:Nbi-1))

 !vector from site i to j's (x,y) components
 !note that x and y can be negative
 allocate(dx(0:Nbi-1,0:Nbi-1))
 allocate(dy(0:Nbi-1,0:Nbi-1))

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

     dx(i,j) = jx-ix
     dy(i,j) = jy-iy
 
     ! consider periodic boundary
     if (abs(dx(i,j))>Nx/2) then
       if (dx(i,j)>0) then
         dx(i,j) = dx(i,j)-Nx
       else
         dx(i,j) = dx(i,j)+Nx
       endif
     endif

     if (abs(dy(i,j))>Ny/2) then
       if (dy(i,j)>0) then
         dy(i,j) = dy(i,j)-Ny
       else
         dy(i,j) = dy(i,j)+Ny
       endif
     endif

     ! set additional symmetry (dx,dy)=(dy,dx)
     !if (dx<dy) then
     !  k = dx
     !  dx = dy
     !  dy = k
     !endif

     k = get_index(dx(i,j),dy(i,j))
     dclass(i,j) = k
     print*, 'class=',k, 'd=',dx(i,j), dy(i,j), '  site',i,j

     ! Note for k = get_index(dx,dy) with dx<dy, dclass_F(k)=0
     dclass_F(k) = dclass_F(k)+1 

     ! Compute exp(i*q*(r_i+r_j))     
     qx = dx(i,j)*2.0d0*pi/Nx
     qy = dy(i,j)*2.0d0*pi/Ny
     expqr(k,i,j) = cos(qx*(ix+jx) + qy*(iy+jy))
   !  print*, 'q=',qx,qy,'r=',ix+jx,iy+jy,'expqr=',expqr(k,i,j)

     ! staggered phase for computing staggered correlation in measurement.f90
     if (mod(dx(i,j)+dy(i,j),2)==0) then
       phase(i,j) = 1.0
     else
       phase(i,j) = -1.0
     endif
    enddo
   enddo
  enddo
 enddo

 !check dclass_F
 do k=0, (Nx+1)*(Ny+1)
   i = mod(k, Nx+1)-Nx/2
   j = k/(Nx+1)-Ny/2
   print*, 'class ',k,'  F=',dclass_F(k), ' d=', i, j
 enddo
 return
 end subroutine get_distance_class

 !=============================================================================                       
 ! get a unique integer to identify a distance class      
 ! (i,j) labels distance = (dx,dy) between sites       
 ! note that dx, dy can be negative                          
 integer function get_index(i,j)                                                 
 implicit none 
 integer i, j
 if(abs(i)>Nx/2 .or. abs(j)>Ny/2)then                                                                                     
  print*, 'Distance vector between two sites exceeds half lattice size!'                          
  stop                                                                                                
 endif     
 get_index = i+Nx/2 + (j+Ny/2)*(Nx+1)                                             
 return    
 end function get_index

end module cluster 
