      subroutine rmean(sig,npts)
      integer npts,i
      real sig(npts),avg/0/
      do i=1,npts
             avg=avg+sig(i)/npts
      enddo
      do i=1,npts
             sig(i)=sig(i)-avg
      enddo
      return
      end subroutine
