        subroutine rtrend(y,npts,dt)
        integer npts,i
        real y(npts),dt,x(npts)
        real sigmax/0/,sigmax2/0/,sigmay/0/,sigmaxy/0/
        real matrix(2,2),a(2),b(2)
        do i=1,npts
               x(i)=(i-1)*dt
               sigmax=sigmax+x(i)
               sigmay=sigmay+y(i)
               sigmax2=sigmax2+x(i)*x(i)
               sigmaxy=sigmaxy+x(i)*y(i)
        enddo
        matrix(1,1)=sigmax2
        matrix(1,2)=sigmax
        matrix(2,1)=sigmax
        matrix(2,2)=real(npts)
        b(1)=sigmaxy
        b(2)=sigmay
        call gauss(matrix,b,2,a)
        do i=1,npts
               y(i)=y(i)-a(1)*x(i)-a(2)
        enddo
        return
        end subroutine
