subroutine gauss(a,b,n,x)
!============================================================
! Solutions to a system of linear equations A*x=b
! Method: the basic elimination (simple Gauss elimination)
! Alex G. November 2009
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! b(n)   - vector of the right hand coefficients b
! n      - number of equations
! output ...
! x(n)   - solutions
! comments ...
! the original arrays a(n,n) and b(n) will be destroyed 
! during the calculation
!===========================================================
!implicit none 
integer n,i,j,k
double precision a(n,n), b(n), x(n),c
!step 1: forward elimination
!write(*,*)'n=',n
do k=1, n-1
        do i=k+1,n
               c=a(i,k)/a(k,k)
               a(i,k) = 0.0
               b(i)=b(i)- c*b(k)
               do j=k+1,n
                        a(i,j) = a(i,j)-c*a(k,j)
               end do
        end do
end do
!step 2: back
!substitution
x(n) = b(n)/a(n,n)
do i=n-1,1,-1
        c=0.0
        do j=i+1,n
              c= c + a(i,j)*x(j)
        end do 
        x(i) = (b(i)- c)/a(i,i)
enddo
end subroutine gauss
