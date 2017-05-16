! measurement the time shift (double difference) using cross-correlation
program cross_rm
implicit none
integer,parameter:: nptsmax=1000000,nsmax=1000000
integer nargc,nfile,i,ncc,j,ncc2,k,id
integer nwin,npts,nerr,nb,nmax,ncon,ns
integer do_filter,tmpn,nn,ii
character(180):: files(100),arg,opt,tr
character(180):: out_cc,out_reg,out_sacfile
real shift_t1
real dist,baz,mag
real devia0,devia1,dum
real x1tmp(nptsmax),out_sac(nptsmax)
real trf,fl,fh,twin,tlen,ccc,a0,b0,c0,tsht
real beg,dt0,t1,dt1,cc(nptsmax),snr0,snr1,snr2,amaxf,amaxb
real sig0(nptsmax),sig1(nptsmax),tmp0(nptsmax),tmp1(nptsmax)
double precision matrix(3,3),res(3),b(3),ytmp(nptsmax),xtmp(nptsmax)

nfile=0;tr='t1'  ! reference time
ncon=5;fl=0.8;fh=2;twin=16;tlen=20;trf=0;do_filter=0
nargc=iargc()
if(nargc.eq.0)then
   call usage
   stop
endif
i=1
do while(i.le.nargc)
   call getarg(i,opt)
   if(opt(1:2).eq.'-T')then            ! maximum of time shift allowed in seconds.
      read(opt(3:),'(bn,f20.0)')tlen  ! maximum of time shift in seconds.
   else if (opt(1:2).eq.'-W')then      ! window length 
      read(opt(3:),'(bn,f20.0)')twin  ! window length in seconds
   else if (opt(1:2).eq.'-R')then
      if (opt(3:3).ne.'t')then     ! reference time
         read(opt(3:),'(bn,f20.0)')trf
      else
         read(opt(3:),'(bn,1a)')tr
      endif
   else if (opt(1:2).eq.'-L')then
      read(opt(3:),'(bn,f20.0)')fl    ! low frequency
   else if (opt(1:2).eq.'-H')then
      read(opt(3:),'(bn,f20.0)')fh    ! high frequency
   else if (opt(1:2).eq.'-N')then
      read(opt(3:),'(bn,i10)')ncon    ! number of points to do regression
   else if (opt(1:2).eq.'-Y')then          ! do filter or not
      do_filter=1
   else if (opt(1:1).ne.'-')then
      nfile=nfile+1
      files(nfile)=opt                ! the sac files, where the first one is the reference
   else
      call usage
      stop
   endif
   i=i+1
enddo
ns=int((ncon-1)/2)                             ! half of the number of points to do regression
call rsac1(files(1),sig0,npts,beg,dt0,nptsmax,nerr) ! read the reference file
if(nerr .NE. 0) then          
   stop 'Error reading in reference file.'
endif
call getfhv('gcarc',dist,nerr)                  ! get the distance of the first sac file
call getfhv('baz',baz,nerr)                    ! get baz of the first file
call getfhv('mag',mag,nerr)                    ! get baz of the first file
!nwin=int(twin/dt0)+1                           ! number of points in the time window
!twin=(mag*(twin-8.0)*2.0/3.0+112.0/3.0-11.0*twin/3.0)
nwin=int((twin+4)/dt0)+1                           ! number of points in the time window
if(trf.eq.0)call getfhv(tr,trf,nerr)           ! get the reference time
if(nerr.ne.0)then
   write(*,*)trim(tr),' not defined for file ',trim(files(1))
   call exit(-1)
endif
if (do_filter.eq.1)call filter(sig0,npts,dble(fl),dble(fh),dt0)   ! band pass filter the file
!out_sacfile=trim(files(1))//'_shift'
!call wsac0(out_sacfile,dum,sig0,nerr)         ! write out the bandpass filtered file
!write(*,*)'hello read in parameter done!'
amaxf=0;amaxb=0
tmpn=int((trf-beg)/dt0)+1                      ! the id of the point at the initial reference time tr
do i=tmpn-nwin+1,tmpn
   amaxf=max(amaxf,abs(sig0(i)))               ! maximum of amplitude before the time window
enddo
devia0=trf-(tmpn-1)*dt0-beg                    ! the reference time residual respect to n times of dt0
!nb=int((trf-beg)/dt0)+1                       ! the id of the point at reference time tr
do i=1,80                                      ! find the onset, works only if t1 is smaller than the onset
       if(amaxf.lt.abs(sig0(tmpn+i)))goto 10
enddo
10 continue
nb=tmpn+i-int(4.0/dt0)                         ! the id of the point at reference time
shift_t1=i*dt0                                 ! time shift of the initial reference, suppose the t1 arrive earlier than onset
!write(*,*)'shift_t1=',shift_t1
out_sacfile=trim(files(1))//'_shift'           ! output shifted file
call setfhv('t1',trf+shift_t1,nerr)            ! set new t1 with respect to shift_t1
call wsac0(out_sacfile,dum,sig0,nerr)          ! write out the bandpass filtered file
! get the snr
amaxf=0;amaxb=0
tmpn=int((trf+shift_t1-beg)/dt0)               ! the id of the point at reference time
do i=tmpn-nwin+1,tmpn
       amaxf=amaxf+sig0(i)**2/nwin             ! maximum of amplitude before the time window
enddo
amaxf=sqrt(amaxf)
do i=nb+1,nb+nwin
       amaxb=max(amaxb,abs(sig0(i)))           ! maximun of the  amplitude in the time window
enddo
snr0=amaxb/amaxf                               ! snr of reference signal

ncc=int(tlen*2/dt0)+1                          ! number of points of the maximum time shift
ncc2=int(tlen/dt0)                             ! half of the number of the maximum time shift
!ncc=int(tlen*2/ts)+1                          ! number of points of the maximum time shift
!call newhdr
      !write(*,*)'ncc=',ncc
!write(*,*)'shift_t1=',shift_t1
do i=1,nwin
        tmp0(i)=sig0(nb+i-1)+devia0*(sig0(nb+i)-sig0(nb+i-1))*dt0     ! reference time series of the window
enddo
call rmean(tmp0,nwin)                          ! remove mean value
call rtrend(tmp0,nwin,dt0)                     ! remove trend
!call taper(tmp0,nwin)                          ! do taper
!write(*,*)'shift_t1=',shift_t1
open(40,file='cross_rm.dat')
!write(40,*)files(1)
do i=2,nfile                                   ! deal with the other sac file
   sig1=0;npts=0;beg=0;dt1=0;trf=0;nb=0;devia1=0
   call rsac1(files(i),sig1,npts,beg,dt1,nptsmax,nerr) ! read the comparison sac file
   if (nerr.ne.0)then
      write(*,*)'Error reading in file: ',trim(files(i))
      continue
   endif
   if (dt0.ne.dt1)then
      write(*,*)'Delta incompatible with reference file'
      continue
   endif
   if(trf.eq.0)call getfhv(tr,trf,nerr)                  ! get the reference time
   if(nerr.ne.0)then
      write(*,*)trim(tr),' is not defined for file ',trim(files(i))
      continue
   endif
   if(do_filter.eq.1)call filter(sig1,npts,dble(fl),dble(fh),dt1) ! bandpass filter the waveform
   nb=int((trf+shift_t1-4.0-beg)/dt1)+1      ! the id of the reference point
   devia1=trf+shift_t1-4.0-(nb-1)*dt1-beg    ! the deivation time
     !write(*,*)'devia1=',devia1,nb
     !write(*,*)'hello, i=',i
     !call filter(sig0,npts,dble(fl),dble(fh),dt0) ! band pass filter the file
   ccc=0
   !write(*,*)'hello, i=',i
   !open(30,file='test.cc')
   cc=0;tmp1=0
   do j=1,ncc           ! loop over the maximum time shift
      tmp1=0
      do k=1,nwin
               !       tmp1(k)=sig1(nb-ncc2+j+k-2)
         tmp1(k)=sig1(nb-ncc2+j+k-2)+devia1*(sig1(nb-ncc2+j+k-1)-sig1(nb-ncc2+j+k-2))*dt0     ! reference time series of the window
      enddo
      call rmean(tmp1,nwin)
      call rtrend(tmp1,nwin,dt1)
           !   call taper(tmp1,nwin,dt1)
      call cross(tmp0,tmp1,nwin,cc(j))
        !      write(40,*)-tlen+(j-1)*dt0,cc(j)
      if (cc(j)>ccc) then
         ccc=cc(j)
         nmax=j               ! id of the peak amplitude in the CC function
      endif
                     !write(30,*)-tlen+(j-1)*dt0,cc(j)
      !                write(*,*)-tlen+(j-1)*dt0,cc(j)
   enddo
      
     !close(30)
              !write(*,*)'nmax=',nmax,'ncon=',ncon,'ns=',ns
   matrix=0;res=0;b=0;ytmp=0;xtmp=0
   nn=3;!write(*,*)'nn=',nn
   do j=1,ncon
      id=nmax-ns+j-1
      ytmp(j)=dble(cc(id))
      xtmp(j)=dble((id-1)*dt0-tlen)
             !write(*,*)xtmp(j),ytmp(j)
      matrix(1,1)=matrix(1,1)+xtmp(j)**4
      matrix(1,2)=matrix(1,2)+xtmp(j)**3
      matrix(1,3)=matrix(1,3)+xtmp(j)**2
             !matrix(2,1)=matrix(2,1)+xtmp(j)**3
             !matrix(2,2)=matrix(2,2)+xtmp(j)**2
      matrix(2,3)=matrix(2,3)+xtmp(j)
             !matrix(3,1)=matrix(3,1)+xtmp(j)**2
             !matrix(3,2)=matrix(3,2)+xtmp(j)
             !matrix(3,3)=matrix(3,3)+1
      b(1)=b(1)+xtmp(j)**2*ytmp(j)
      b(2)=b(2)+xtmp(j)*ytmp(j)
      b(3)=b(3)+ytmp(j)
   enddo
   matrix(2,1)=matrix(1,2)
   matrix(2,2)=matrix(1,3)
   matrix(3,1)=matrix(1,3)
   matrix(3,2)=matrix(2,3)
   matrix(3,3)=ncon
!     continue
!     call gauss(matrix,b,res,nn)
     !write(*,*)matrix(1,1),matrix(1,2),matrix(1,3)
     !write(*,*)matrix(2,1),matrix(2,2),matrix(2,3)
     !write(*,*)matrix(3,1),matrix(3,2),matrix(3,3)
     !write(*,*)b(1),b(2),b(3)
   call gauss(matrix,b,nn,res)
   a0=res(1)                                  !y=a0(x+b0)^2+c0
   b0=res(2)/2.0/res(1)                       ! time shift
   c0=(4*res(1)*res(3)-res(2)**2)/4.0/res(1)  ! cross correlation coefficient
         !open(40,file=out_reg)
!     write(*,*)res(1),res(2),res(3)
   do j=1,20
      x1tmp(j)=(xtmp(ncon)-xtmp(1))*(j-1)/19.0+xtmp(1)
      ytmp(j)=a0*(x1tmp(j)+b0)**2+c0
      write(40,*)x1tmp(j),ytmp(j)
   enddo
!         write(*,*)trim(files(i)),ccc,-tlen+(nmax-1)*dt0,(4*res(1)*res(3)-res(2)**2)/4.0/res(1),-res(2)/2.0/res(1),res(1)
   tsht=-tlen+(nmax-1)*dt0
             !write(*,'(1a,1x,8f15.4)')trim(files(i)),dist,baz,snr,ccc,tsht,c0,-b0+devia1-devia0,a0
     !write(40,'(1a,1x,8f15.4)')trim(files(i)),dist,baz,ccc,tsht,c0,-b0
              !call setfhv('t1',trf+tsht,nerr)
   t1=trf-b0+shift_t1
   call setfhv('t1',t1,nerr)
              !call setfhv('t1',trf-b0+devia1-devia0,nerr)
   out_sacfile=trim(files(i))//'_shift'
   call wsac0(out_sacfile,dum,sig1,nerr)
     ! get the snr
   amaxf=0;amaxb=0
   tmpn=int((t1-beg)/dt0)+1                      ! the id of the point at the initial reference time tr
   do ii=tmpn-nwin+1,tmpn
      amaxf=amaxf+sig1(ii)**2/nwin             ! maximum of amplitude before the time window
   enddo
   do ii=tmpn+1,tmpn+nwin
      amaxb=max(amaxb,abs(sig1(ii)))           ! maximun of the  amplitude in the time window
   enddo
   amaxf=sqrt(amaxf)
   snr1=amaxb/amaxf                              ! snr of the signal
   write(*,'(1a,1x,10f10.3)') trim(files(i)),snr0,snr1,mag,twin,dist,baz,ccc,tsht,c0,-b0
enddo
close(40)
end program cross_rm
