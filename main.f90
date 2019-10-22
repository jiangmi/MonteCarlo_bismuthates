include 'random.f90'
include 'parameters.f90'
include 'cluster.f90'
include 'monte_carlo.f90'
include 'measurements.f90'
!==============================================================================
!==============================================================================
program main
use random
use parameters
use cluster
use monte_carlo
use measurements
implicit none
integer ii,jj,i,j,k,ix,iy, bin, ibeta
double precision mean,std,Emean,Estd,Nmean,Nstd, &
                 spmean,spstd,bpmean,bpstd,      &
                 swave_s_mean, swave_s_std,      &
                 swave_px_mean, swave_px_std,    &
                 swave_py_mean, swave_py_std
double precision means(12)  ! store sublat dependent quantities
double precision t1, t2, Xval, Etot, Xavg
double precision, dimension(:),   allocatable :: X
double precision, dimension(:,:), allocatable :: H0

call init_parameters()     ! classify distance in cluster.f90
call allocate_quantities() ! set physical quantites
500 format(a20,i7,a3,i7,a20,f8.5)
600 format(a30,f8.5,a8,f8.5)
call get_distance_class()

call cpu_time(t1)

print*, 'Enter a random seed for this run:'
!read(5,*) iran
iran = -312434;
!allocate space for the non-interacting Hamiltonian and displacements
allocate(X(0:N-1))
allocate(H0(0:N-1,0:N-1))

!Initialize the displacement to zero and then a small value
X = 0.0d0
!if (if_X_displace==1) then
  !random initial X
  do i = 0,N-1
   if(.not.is_bismuth(i))then
    X(i) = (ran2(iran)-0.50d0)*2.d0*dXamp
   endif
  enddo

  !initial bond disproportionated lattice distortion
  !Xval = 0.15
  !do ix = 0,Nx-1
  !  do iy = 0,Ny-1
  !    !px orbital's displacement
  !    i = return_index_for_coordinates(ix,iy,1)
  !    if (mod(ix+iy,2)==1) then
  !      X(i) = Xval
  !    else
  !      X(i) = -Xval
  !    endif
  !    !py orbital's displacement
  !    i = return_index_for_coordinates(ix,iy,2)                                    
  !    if (mod(ix+iy,2)==1) then                                                    
  !      X(i) = Xval                                                                
  !    else
  !      X(i) = -Xval                                                               
  !    endif
  !  enddo
  !enddo
!endif

! Loop over temperature starting from highest T
do ibeta = 1,nbeta
 beta = betas(ibeta)
 print*, '  '
 print*, '==========================================================='
 print 600, 'Start carrying out MC for beta = ', beta, 'T = ', 1./beta
 
 ! warmup begins
 accept = 0
 reject = 0
 do i = 1,nwarms
  if(mod(i,nwarms/ninv).eq.0)then
   print 500, 'Completed Warmup', i, 'of', nwarms, &
           'Accept ratio = ', dfloat(accept)/(dfloat(accept+reject))
   accept = 0
   reject = 0
  endif
  !print*, '------------------------'
  !print*, 'warmup sweep No.',i
  call single_site_sweep(X,accept,reject)
 enddo

 print *, 'warmup finished, X='
 do ii=Nbi,N-1
   print *, X(ii)
 enddo
 print *, '<|X|> =', sum(abs(X(Nbi:N-1)))/(N*2/3)

 call compute_total_E(Etot, X)
 print *, 'warmup finished, total E=', Etot
 
 ! measurements begins
 accept = 0
 reject = 0
 call zero_accumulators()
 bin = 0
 do i = 1,nmeas
  if(mod(i,nmeas/nbin).eq.0)then
   print 500, 'Completed Nmeas', i, 'of', nmeas,  &
           'Accept ratio = ', dfloat(accept)/(dfloat(accept+reject))
   accept = 0
   reject = 0
   bin = bin + 1
   call populate_bins(bin)
   call zero_accumulators()
  endif
  !print*, '------------------------'
  !print*, 'measurement sweep No.',i
  call single_site_sweep(X,accept,reject)
  call do_measurements(X)
 enddo

 print *, 'meas finished, X='
 do ii=Nbi,N-1
   print *, X(ii)
 enddo
 Xavg = sum(abs(X(Nbi:N-1)))/(N*2/3)
 print *, '<|X|> =', Xavg

 call compute_total_E(Etot, X)
 print *, 'meas finished, total E=', Etot
 print *, '===================================='
 
 include 'output_results.f90' 
enddo

call deallocate_quantities() 

call cpu_time(t2)
write (*,*) 'Elapsed CPU time = ', t2-t1

stop
end program main
