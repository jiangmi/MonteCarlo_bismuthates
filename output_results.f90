! Define any formats needed for output
100 format(a40,' ',f13.5,' +- ',f13.5)
110 format(a40,' ',f9.5)
200 format(a6, i3, a30, ' ',f8.5,' +- ',f8.5)
300 format(a6, i3, ' ', i3, a30, ' ',f8.5,' +- ',f8.5)
400 format(a6, f4.1, a1, f4.1, a1)
410 format(a5, i3, a11, f9.5, f9.5)
420 format(a5, i3, i3, a11, f9.5, f9.5) 
800 format(f4.1,' ',f10.5,' ',f10.5,' ',f10.5,' ',f10.5,' ',f10.5,' ',f10.5,' ',f10.5,' ',f10.5)

!print results into out.txt
write(unit=6,fmt=110) '  es = ', es
write(unit=6,fmt=110) '  ep = ', ep
write(unit=6,fmt=110) '  tsp = ', tsp
write(unit=6,fmt=110) '  tpp = ', tpp
write(unit=6,fmt=110) '  mu = ', mu
write(unit=6,fmt=110) '  alpha = ', alpha
write(unit=6,fmt=110) '  spring_const = ', spring_const
write(unit=6,fmt=110) '  beta = ', beta

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
call get_err(bnpa1g_avg,mean,std)
write(unit=6,fmt=100) '  A1g Oxygen holes = ', mean, std
call get_err(bsp_site_avg,spmean,spstd)
write(unit=6,fmt=100) '  Lattice-averaged single polaron = ', spmean, spstd

do ii = 0,Nbi-1                                                                                       
  call get_err(bnbis(ii,:),mean,std)                                                              
  write(unit=6,fmt=200) 'cell', ii, '  <n_Bi> = ', mean, std                                  
enddo

do ii = 0,Nbi-1
  call get_err(bnpa1g(ii,:),mean,std)    
  write(unit=6,fmt=200) 'cell', ii, '  <n_p_A1g> = ', mean, std        
enddo

do ii = 0,Nbi-1
  call get_err(bspolaron(ii,:),mean,std)
  write(unit=6,fmt=200) 'cell', ii, '  single polaron = ', mean, std
enddo

call get_err(bbp_site_avg,bpmean,bpstd)
write(unit=6,fmt=100) '  Lattice-averaged bipolaron = ', bpmean, bpstd

!print results into data.txt
write(unit=11,fmt=800) beta, Nmean, Nstd, Emean, Estd, spmean, spstd, bpmean, bpstd

do ii = 0,Nbi-1
  call get_err(bbpolaron(ii,:),mean,std)
  write(unit=6,fmt=200) 'cell', ii, '  bipolaron = ', mean, std
enddo

!print single polaron correlation function
do k = 0,nclass-1
  if (dclass_F(k)/=0) then
    call get_err(bspolaron_ij(k,:),mean,std)
    i = mod(k,(Nx/2+1))
    j = k/(Nx/2+1)
    write(unit=6,fmt=300) 'dist', i, j, '  <L_i*L_j> = ', mean, std
  endif
enddo

!print s-wave susceptibility chi_sc(orb)                               
do o1 = 0,Norb-1
  call get_err(bchi_sc(o1,:),mean,std)
  write(unit=6,fmt=410) 'orb', o1, '  chi_sc = ', mean, std
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

