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
integer nbeta, ii,jj,i,j, bin
double precision mean, std
double precision E, t1, t2
double precision, dimension(:),   allocatable :: X
double precision, dimension(:,:), allocatable :: H0
call init_parameters()     ! set all the relevant parameters
call allocate_quantities() ! set physical quantites
500 format(a20,i7,a3,i7,a20,f8.5)
600 format(a30,f8.5,a8,f8.5)

call cpu_time(t1)

print*, 'Enter a random seed for this run:'
!read(5,*) iran
iran = -312434;
!allocate space for the non-interacting Hamiltonian and displacements
allocate(X(0:N-1))
allocate(H0(0:N-1,0:N-1))

!Initialize the displacement to zero and then a small value
X = 0.0d0
do i = 0,N-1
 if(.not.is_bismuth(i))then
  X(i) = (ran2(iran)-0.50d0)*0.001d0
 endif
enddo

! Loop over temperature starting from highest T
do nbeta = 0,num_beta_steps
 beta = beta_min + dfloat(nbeta)*(beta_max-beta_min)/dfloat(num_beta_steps)
 print*, '  '
 print 600, 'Carrying out MC for beta = ', beta, 'T = ', 1./beta
 
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
  call single_site_sweep(X,accept,reject)
 enddo
 
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
  call single_site_sweep(X,accept,reject)
  call do_measurements(X)
 enddo
 
 include 'output_results.f90' 
enddo

call deallocate_quantities() 

call cpu_time(t2)
write (*,*) 'Elapsed CPU time = ', t2-t1

stop
end program main
