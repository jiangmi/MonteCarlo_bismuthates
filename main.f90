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
integer ii,jj,i,j,k,ix,iy, bin, ibeta, idstart
double precision mean,std,Emean,Estd,Nmean,Nstd, &
                 spmean,spstd,bpmean,bpstd,      &
                 swave_s_mean, swave_s_std,      &
                 swave_px_mean, swave_px_std,    &
                 swave_py_mean, swave_py_std
character(len=6) str1
double precision means(12)  ! store sublat dependent quantities
double precision t1, t2, t3, t4, t5, t6, t7, tsweep, teigenvector,tmeas
double precision Xval, Etot, Xavg
character(len=:), allocatable :: str2
double precision, dimension(:),   allocatable :: X
double precision, dimension(:),   allocatable :: X_input
double precision, dimension(:,:), allocatable :: H0

call init_parameters()     ! classify distance in cluster.f90
call allocate_quantities() ! set physical quantites
500 format(a20,i7,a3,i7,a20,f8.5)
600 format(a30,f8.5,a8,f8.5)

if (travel_cluster==1) then
  print*, 'This run uses traveling cluster approximation'
  call find_cluster()
endif
call get_distance_class()

print*, 'Enter a random seed for this run:'
!read(5,*) iran
iran = -312434;
!allocate space for the non-interacting Hamiltonian and displacements
allocate(X(0:N-1))
if (if_use_equilibrium_X==1) then
 allocate(X_input(2*Nbi))
endif
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

 call cpu_time(t1)

 ! Decide if using equilibrium X at the corresponding T in previous run
 ! so that the measurements and warmup can be separately performed
 ! Otherwise, do nothing means using X from the last T
 if (if_use_equilibrium_X==1) then
   write(str1,'(f6.1)') beta
   if (abs(mu)<1.e-4) then
     allocate(character(len=4) :: str2)
     write(str2,'(f4.2)') mu
   elseif (mu>0.0) then
     allocate(character(len=4) :: str2)
     write(str2,'(f4.2)') mu
   else
     allocate(character(len=5) :: str2)
     write(str2,'(f5.2)') mu
   endif

   open(35, file='X_readin_mu'//adjustl(trim(str2))//'_be'//adjustl(trim(str1)), status='old')
   do i = 1,2*Nbi
     read(35,*) X_input(i)
   enddo

 !  idstart = 7*2*Nbi+1
   X = 0.0d0
   print *, 'readin equilibrium X of previous run done, X='
   do i=Nbi,N-1
     X(i) = X_input(i-Nbi+1)
     print *, X(i)
   enddo
   print *, '<|X|> =', sum(abs(X(Nbi:N-1)))/(N*2/3)

   close(35)
 endif

 if (if_warmup/=1) then
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
    if (travel_cluster==0) then
      call single_site_sweep(X,accept,reject)
    else
      call single_site_sweep_cluster(X,accept,reject)
    endif
   enddo

   print *, 'warmup finished, X='
   do ii=Nbi,N-1
     print *, X(ii)
   enddo
   print *, '<|X|> =', sum(abs(X(Nbi:N-1)))/(N*2/3)

   call compute_total_E(Etot, X)
   print *, 'warmup finished, total E=', Etot

   call cpu_time(t2)
 endif ! end if_warmup

 ! measurements begins
 tsweep = 0.0
 teigenvector = 0.0
 tmeas = 0.0
 t5 = 0.0
 t6 = 0.0

 accept = 0
 reject = 0
 call zero_accumulators()
 bin = 0
 do i = 1,nmeas
  !print*, '------------------------'
  !print*, 'measurement sweep No.',i

  call cpu_time(t3)
  if (travel_cluster==0) then
    call single_site_sweep(X,accept,reject)
  else
    call single_site_sweep_cluster(X,accept,reject)
  endif
  call cpu_time(t4)

  tsweep = tsweep + t4-t3

  !Do measurement per measinv MC sweeps
  if (mod(i,measinv).eq.0) then
    call do_measurements(X,t5,t6,i)
    teigenvector = teigenvector + t5
    tmeas = tmeas + t6
  endif
 
  if(mod(i,nmeas/nbin).eq.0)then
   print 500, 'Completed Nmeas', i, 'of', nmeas,  &
           'Accept ratio = ', dfloat(accept)/(dfloat(accept+reject))
   accept = 0
   reject = 0
   bin = bin + 1
   call populate_bins(bin)
   call zero_accumulators()
  endif
 enddo ! end of measurement

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

 call cpu_time(t7)
 write (*,*) 'Warmup CPU time =               ', t2-t1
 write (*,*) 'Sweep time during meas =        ', tsweep
 write (*,*) 'Diag (eigenpairs) during meas = ', teigenvector
 write (*,*) 'Meas quantities time =          ', tmeas
 write (*,*) 'Total CPU time =                ', t7-t1
enddo ! end of loop over beta

call deallocate_quantities() 

stop
end program main
