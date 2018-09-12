! The order of N orbitals is:
! first 0:NBi-1 is bismuth and then px orbital and finally py orbital 
! N = Nx*Ny*Norb
! NBi = Nx*Ny
module cluster 
 use parameters, only: N, Nx, Ny, Nbi
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
end module cluster 
