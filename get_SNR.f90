!calculate the SNR
program get_SNR
integer nptsmax
parameter (nptsmax=1000000)
real amaxf,amaxb
real sig(nptsmax)
real dist,baz,mag
real beg,dt,t1,trf,snr(10)
real f1(10),f2(10),twin1,twin2
integer i,j,k,id,num,im
integer nwin,npts,nerr,nb,nmax
integer tmpn,nn,ii,lenz,nwin1,nwin2
character(180)saclist,tr
character(180)files(2),para,output
logical ext
nargc=iargc()
if (nargc.ne.1)then
       write(*,*)'Usage: get_SNR parameter_file'
       write(*,*)'Parameter_file is like:'
       write(*,*)'list           [All z component sac file]'
       write(*,*)'twin1 twin2    [window length before and after t1 in seconds]'
       write(*,*)'id t_ref[t1/t2/t8] [ using method 1 or 2, reference time mark]'
       write(*,*)'output         [output file name]'
       write(*,*)'id num         [do band pass filter (1) or not (0), number of frequency band]'
       write(*,*)'f1 f2          [low frequency, high frequency]'
       write(*,*)'...'
       write(*,*)'Caution:'
       write(*,*)'1, all the files should be named as: eq.z eq.r eq.t'
       write(*,*)'2, P arrival time is marked in SAC as t1'
       stop
endif
call getarg(1,para)
open(11,file=para)
read(11,'(a180)')saclist
read(11,*)twin1,twin2
read(11,*)im,tr
read(11,'(a180)')output
read(11,*)id,num
!write(*,*)twin1,twin2
if(id.eq.1)then
     if (num.eq.0)then
          write(*,*)'Something wrong with the frequecy band, we will not do band pass filter!'
          goto 12
     endif
     do i=1,num
          read(11,*)f1(i),f2(i)
     enddo
else
     write(*,*)'We will not do band pass filter!'
endif
12 close(11)
open(100,file=output)
open(12,file=saclist)
13   read(12,*,err=14,end=14)files(1)
     lenz=len(trim(files(1)))-1   ! length of the sac file minus 1
     files(2)=files(1)(1:lenz)//"r"
     do i=1,2
            inquire(exist=ext,file=files(i))
            if (.not.ext)then
                  write(*,*)'Caution, file ',trim(files(i)),' doesnot exist!'
                  continue
                  goto 13
            endif
            call rsac1(files(i),sig,npts,beg,dt,nptsmax,nerr) ! read the sac file
            if(nerr .NE. 0) then          
                  write(*,*)'Error in reading file: ',trim(files(i))
                  continue
            endif
            !call rmean(sig,npts)                       ! remove mean value
            !call rtrend(sig,npts,dt)                   ! remove trend
            !call getfhv('gcarc',dist,nerr)             ! get the distance of the first sac file
            !call getfhv('baz',  baz, nerr)             ! get baz of the first file
            !call getfhv('mag',  mag, nerr)             ! get baz of the first file
            call getfhv(tr,   trf, nerr)                ! get the reference time, which is t1/t2/t8
            trf=trf-2.0
            if(nerr.ne.0)then
                  write(*,*)'T1 is not defined for file ',trim(files(i))
                  continue
            endif
            nwin1=int(twin1/dt)+1                      ! number of points in the time window before t1
            nwin2=int(twin2/dt)+1                      ! number of points in the time window after t1
            !write(*,*)nwin1,nwin2
            amaxf=0;amaxb=0
            tmpn=int((trf-beg)/dt)+1                   ! the id of the point at the initial reference time (t1)
            do ii=tmpn-nwin1+1,tmpn
                  if (ii.gt.1)then
                      amaxf=max(amaxf,abs(sig(ii)))        ! maximum of amplitude before the time window
                  endif
            enddo
            do ii=1,nwin1                              ! search the onset, works only if t1 is smaller than the onset
                  if(amaxf.lt.abs(sig(tmpn+ii)))goto 15
            enddo
      15    continue
            nb=tmpn+ii                                 ! the new id of the point at reference time
            if (id.eq.1)then                           ! do filter
                  do nf=1,num                          ! loop over frequency band
                       call filter(sig,npts,dble(f1(nf)),dble(f2(nf)),dt)   ! band pass filter the file
                       amaxf=0;amaxb=0
                       if (im.eq.1)then ! use ratio between maximum amplitude and root-mean-square
                             amaxf=sum(sig(nb-nwin1+1:nb)**2)/nwin1 ! standard error of the noise before the reference time
                             amaxf=sqrt(amaxf)
                             do ii=nb+1,nb+nwin2
                                    amaxb=max(amaxb,abs(sig(ii)))  ! maximun of the  amplitude in the time window
                             enddo
                             snr(nf)=amaxb/amaxf                  ! SNR of P wave
                       else  ! use ratio between standard variances
                             amaxf=sum(sig(nb-nwin1+1:nb))/nwin1
                             amaxf=sum((sig(nb-nwin1+1:nb)-amaxf)**2)
                             amaxb=sum(sig(nb:nwin1+nb-1))/nwin1
                             amaxb=sum((sig(nb:nwin1+nb-1)-amaxb)**2)
                             snr(nf)=amaxb/amaxf                  ! SNR of P wave
                       endif
                  enddo                                     ! end loop over frequency band
            else                                            ! if do not apply band pass filter
                  amaxf=0;amaxb=0
                  if(im.eq.1)then
                       amaxf=sum(sig(nb-nwin1+1:nb)**2)/nwin1      ! standard error of the noise before the reference time
                       amaxf=sqrt(amaxf)
                       do ii=nb+1,nb+nwin2
                             amaxb=max(amaxb,abs(sig(ii)))        ! maximun of the  amplitude in the time window
                       enddo
                       snr(1)=amaxb/amaxf                        ! SNR of P wave
                  else  ! use ratio between standard variances
                       amaxf=sum(sig(nb-nwin1+1:nb))/nwin1
                       amaxf=sum((sig(nb-nwin1+1:nb)-amaxf)**2)
                       amaxb=sum(sig(nb:nwin1+nb-1))/nwin1
                             amaxb=sum((sig(nb:nwin1+nb-1)-amaxb)**2)
                       snr(1)=amaxb/amaxf                  ! SNR of P wave
                  endif
                  !snr(1)=amaxb/amaxf                        ! SNR of P wave
                  num=1
            endif ! if do not apply band pass filter
            write(100,*)trim(files(i)),(snr(j),j=1,num)
            write(*,*)trim(files(i)),(snr(j),j=1,num)
     enddo        ! end loop over sac file
     goto 13
14   continue
close(12)
close(100)
end program
