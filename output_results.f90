! Define any formats needed for output
100 format(a40,' ',f13.5,' +- ',f13.5)
110 format(a40,' ',f9.5)
200 format(a6, i3, a30, ' ',f8.5,' +- ',f8.5)
300 format(i3, ' ', i3, ' ',f13.8,' ',f13.8, ' ',f13.8,' ',f13.8)
400 format(a6, f4.1, a1, f4.1, a1)
410 format(a18, f9.5, f9.5)
420 format(a5, i3, i3, a11, f9.5, f9.5) 
800 format(f5.1,' ',f10.5,' ',f10.5,' ',f10.5,' ',f10.5,' ', &
           f10.5,' ',f10.5,' ',f10.5,' ',f10.5, ' ', f10.5, ' ', &
           f10.5,' ',f10.5,' ',f10.5,' ',f10.5,' ',f10.5,' ',f10.5)
810 format(f5.1,' ',f10.5,' ',f10.5,' ',f10.5,' ',f10.5,' ', &
           f10.5,' ',f10.5,' ',f10.5,' ',f10.5, ' ', f10.5, ' ', &
           f10.5,' ',f10.5,' ',f10.5,' ',f10.5)
820 format(f13.5, f13.5)

!print results into out.txt
write(unit=6,fmt=110) '  es = ', es
write(unit=6,fmt=110) '  ep = ', ep
write(unit=6,fmt=110) '  tsp = ', tsp
write(unit=6,fmt=110) '  tpp = ', tpp
write(unit=6,fmt=110) '  mu = ', mu
write(unit=6,fmt=110) '  alpha = ', alpha
write(unit=6,fmt=110) '  spring_const = ', spring_const
write(unit=6,fmt=110) '  beta = ', beta

!print the spectrum of H
do ii = 0,N-1
  call get_err(bEk(ii,:),mean,std)
  write(unit=30,fmt=820)  mean, std
enddo

call get_err(benergy,Emean,Estd)
write(unit=6,fmt=100) '  E_avg = ', Emean, Estd
call get_err(bX,mean,std)
write(unit=6,fmt=100) '  <X> = ', mean, std
call get_err(bX2,mean,std)
write(unit=6,fmt=100) '  <X2> = ', mean, std
call get_err(bntot,Nmean,Nstd)
write(unit=6,fmt=100) '  total filling = ', Nmean, Nstd
call get_err(bnbis_avg,mean,std)
write(unit=6,fmt=100) '  Bi filling = ', mean, std
call get_err(bnox,mean,std)
write(unit=6,fmt=100) '  O(x) filling = ', mean, std
call get_err(bnoy,mean,std)
write(unit=6,fmt=100) '  O(y) filling = ', mean, std
call get_err(bnpLs_avg,mean,std)
write(unit=6,fmt=100) '  Ls Oxygen holes = ', mean, std
call get_err(bnpLd_avg,mean,std)
write(unit=6,fmt=100) '  Ld Oxygen holes = ', mean, std
call get_err(bnpLx_avg,mean,std)
write(unit=6,fmt=100) '  Lx Oxygen holes = ', mean, std
call get_err(bnpLy_avg,mean,std)
write(unit=6,fmt=100) '  Ly Oxygen holes = ', mean, std
call get_err(bsp_site_avg,spmean,spstd)
write(unit=6,fmt=100) '  Lattice-averaged single polaron = ', spmean, spstd

!compute sublat dependent quantities
call get_err(bnbi1,means(1),std)
write(unit=6,fmt=100) '  <n_Bi> sublat1 = ', means(1), std
call get_err(bnbi2,means(2),std)
write(unit=6,fmt=100) '  <n_Bi> sublat2 = ', means(2), std

call get_err(bnpx1,means(3),std)
write(unit=6,fmt=100) '  <n_px> sublat1 = ', means(3), std
call get_err(bnpx2,means(4),std)
write(unit=6,fmt=100) '  <n_px> sublat2 = ', means(4), std

call get_err(bnpy1,means(5),std)
write(unit=6,fmt=100) '  <n_py> sublat1 = ', means(5), std
call get_err(bnpy2,means(6),std)
write(unit=6,fmt=100) '  <n_py> sublat2 = ', means(6), std

call get_err(bnpLs1,means(7),std)
write(unit=6,fmt=100) '  <n_Ls> sublat1 = ', means(7), std
call get_err(bnpLs2,means(8),std)
write(unit=6,fmt=100) '  <n_Ls> sublat2 = ', means(8), std

call get_err(bsp1,means(9),std)
write(unit=6,fmt=100) '  <Sp> sublat1 = ', means(9), std
call get_err(bsp2,means(10),std)
write(unit=6,fmt=100) '  <Sp> sublat2 = ', means(10), std

call get_err(bbp1,means(11),std)
write(unit=6,fmt=100) '  <Bp> sublat1 = ', means(11), std
call get_err(bbp2,means(12),std)
write(unit=6,fmt=100) '  <Bp> sublat2 = ', means(12), std

!quantities for each site
do ii = 0,Nbi-1                                                                                       
  call get_err(bnbis(ii,:),mean,std)                                                              
  write(unit=6,fmt=200) 'cell', ii, '  <n_Bi> = ', mean, std                                  
enddo

do ii = 0,Nbi-1
  call get_err(bnpLs(ii,:),mean,std)    
  write(unit=6,fmt=200) 'cell', ii, '  <n_p_Ls> = ', mean, std        
enddo
do ii = 0,Nbi-1
  call get_err(bnpLd(ii,:),mean,std)
  write(unit=6,fmt=200) 'cell', ii, '  <n_p_Ld> = ', mean, std
enddo
do ii = 0,Nbi-1
  call get_err(bnpLx(ii,:),mean,std)
  write(unit=6,fmt=200) 'cell', ii, '  <n_p_Lx> = ', mean, std
enddo
do ii = 0,Nbi-1
  call get_err(bnpLy(ii,:),mean,std)
  write(unit=6,fmt=200) 'cell', ii, '  <n_p_Ly> = ', mean, std
enddo

do ii = 0,Nbi-1
  call get_err(bspolaron(ii,:),mean,std)
  write(unit=6,fmt=200) 'cell', ii, '  single polaron = ', mean, std
enddo

call get_err(bbp_site_avg,bpmean,bpstd)
write(unit=6,fmt=100) '  Lattice-averaged bipolaron = ', bpmean, bpstd


do ii = 0,Nbi-1
  call get_err(bbpolaron(ii,:),mean,std)
  write(unit=6,fmt=200) 'cell', ii, '  bipolaron = ', mean, std
enddo

!print s-wave susceptibility chi_sc(orb)                               
call get_err(bchi_sc(0,:), swave_s_mean, swave_s_std)
write(unit=6,fmt=410) 'orb  s,  chi_sc = ', swave_s_mean, swave_s_std
call get_err(bchi_sc(1,:), swave_px_mean, swave_px_std)
write(unit=6,fmt=410) 'orb px,  chi_sc = ', swave_px_mean, swave_px_std
call get_err(bchi_sc(2,:), swave_py_mean, swave_py_std)
write(unit=6,fmt=410) 'orb py,  chi_sc = ', swave_py_mean, swave_py_std

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!print results into data.txt
write(unit=11,fmt=800) beta, Nmean, Nstd, Emean, Estd, Xavg, &
                       spmean, spstd, bpmean, bpstd,   &
                       swave_s_mean, swave_s_std,      &
                       swave_px_mean, swave_px_std,    &
                       swave_py_mean, swave_py_std
write(unit=12,fmt=810) beta, Nmean, means(1:12)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!print (bi)polaron correlation function into separate file
do k = 0,nclass-1
  if (dclass_F(k)/=0) then
    call get_err(bspolaron_ij(k,:),spmean,spstd)
    call get_err(bbpolaron_ij(k,:),bpmean,bpstd)
    i = mod(k, Nx+1)-Nx/2
    j = k/(Nx+1)-Ny/2
    write(unit=7,fmt=300) i, j, spmean, spstd, bpmean, bpstd
  endif
enddo

!print charge susceptibility chi_ch(q)
do k = 0,nclass-1
  if (dclass_F(k)/=0) then
    i = mod(k,(Nx/2+1))
    j = k/(Nx/2+1)
    qx = i*2.0d0/Nx
    qy = j*2.0d0/Ny
    write(unit=6,fmt=400) 'q=pi*(', qx, ',', qy, ')'

    do o1 = 0,Norb-1
      do o2 = o1,Norb-1
         call get_err(bchi_ch(k,o1,o2,:),mean,std)
         write(unit=6,fmt=420) 'orb', o1, o2, '  chi_ch = ', mean, std
      enddo
    enddo
  endif
enddo

