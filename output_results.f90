! Define any formats needed for output
100 format(a40,' ',f8.5,' +- ',f8.5)
200 format(a6, i3, a40, ' ',f8.5,' +- ',f8.5)
300 format(a6, i3, ' ', i3, a30, ' ',f8.5,' +- ',f8.5)
400 format(a40,' ',f10.5,' +- ',f8.5)

!code to output the results
call get_err(benergy,mean,std)
write(unit=6,fmt=400) '  E_tot = ', mean, std
call get_err(bX,mean,std)
write(unit=6,fmt=100) '  <X> = ', mean, std
call get_err(bX2,mean,std)
write(unit=6,fmt=100) '  <X2> = ', mean, std
call get_err(bntot,mean,std)
write(unit=6,fmt=100) '  total filling = ', mean, std
call get_err(bnbis,mean,std)
write(unit=6,fmt=100) '  Bi filling = ', mean, std
call get_err(bnox,mean,std)
write(unit=6,fmt=100) '  O(x) filling = ', mean, std
call get_err(bnoy,mean,std)
write(unit=6,fmt=100) '  O(y) filling = ', mean, std

call get_err(bnpa1g,mean,std)
write(unit=6,fmt=100) '  A1g Oxygen holes = ', mean, std

do ii = 0,Nbi-1
  call get_err(bspolaron(:,ii),mean,std)
  write(unit=6,fmt=200) 'site', ii, '  single polaron = ', mean, std
enddo

do ii = 0,Nbi-1
  call get_err(bbpolaron(:,ii),mean,std)
  write(unit=6,fmt=200) 'site', ii, '  bipolaron = ', mean, std
enddo

!do ii = 0,Nbi-1
!  do jj = 0,Nbi-1
!    call get_err(bspolaron_ij(:,ii,jj),mean,std)
!    write(unit=6,fmt=300) 'site', ii, jj, '  <L_i*L_j> = ', mean, std
!  enddo
!enddo
