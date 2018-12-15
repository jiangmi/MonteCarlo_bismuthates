! Define any formats needed for output
100 format(a40,' ',f9.5,' +- ',f9.5)
200 format(a6, i3, a30, ' ',f8.5,' +- ',f8.5)
300 format(a6, i3, ' ', i3, a30, ' ',f8.5,' +- ',f8.5)
800 format(f4.1,' ',f10.5,' ',f10.5,' ',f10.5,' ',f10.5,' ',f10.5,' ',f10.5,' ',f10.5,' ',f10.5)

!print results into out.txt
call get_err(benergy,Emean,Estd)
write(unit=6,fmt=100) '  E_avg = ', Emean, Estd
call get_err(bX,mean,std)
write(unit=6,fmt=100) '  <X> = ', mean, std
call get_err(bX2,mean,std)
write(unit=6,fmt=100) '  <X2> = ', mean, std
call get_err(bntot,Nmean,Nstd)
write(unit=6,fmt=100) '  total filling = ', Nmean, Nstd
call get_err(bnbis,mean,std)
write(unit=6,fmt=100) '  Bi filling = ', mean, std
call get_err(bnox,mean,std)
write(unit=6,fmt=100) '  O(x) filling = ', mean, std
call get_err(bnoy,mean,std)
write(unit=6,fmt=100) '  O(y) filling = ', mean, std
call get_err(bnpa1g,mean,std)
write(unit=6,fmt=100) '  A1g Oxygen holes = ', mean, std
call get_err(bsp_site_avg,spmean,spstd)
write(unit=6,fmt=100) '  Site-averaged single polaron = ', spmean, spstd

do ii = 0,Nbi-1
  call get_err(bspolaron(:,ii),mean,std)
  write(unit=6,fmt=200) 'site', ii, '  single polaron = ', mean, std
enddo

call get_err(bbp_site_avg,bpmean,bpstd)
write(unit=6,fmt=100) '  Site-averaged bipolaron = ', bpmean, bpstd

!print results into data.txt
write(unit=11,fmt=800) beta, Nmean, Nstd, Emean, Estd, spmean, spstd, bpmean, bpstd

do ii = 0,Nbi-1
  call get_err(bbpolaron(:,ii),mean,std)
  write(unit=6,fmt=200) 'site', ii, '  bipolaron = ', mean, std
enddo
!
do k = 0,nclass-1
  if (dclass_F(k)/=0) then
    call get_err(bspolaron_ij(:,k),mean,std)
    i = mod(k,(Nx/2+1))
    j = k/(Nx/2+1)
    write(unit=6,fmt=300) 'site', i, j, '  <L_i*L_j> = ', mean, std
  endif
enddo
