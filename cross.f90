subroutine cross(x,y,n,r)
integer:: n
real(4):: x(1000000),y(1000000),r
real(8):: ex,ey,sumxx,sumyy,sumxy
ex=dble(sum(x(1:n))/n)
ey=dble(sum(y(1:n))/n)
sumxx=sum((dble(x(1:n))-ex)*(dble(x(1:n))-ex))
sumyy=sum((dble(y(1:n))-ey)*(dble(y(1:n))-ey))
sumxy=sum((dble(x(1:n))-ex)*(dble(y(1:n))-ey))
r=sngl(sumxy/dsqrt(sumxx*sumyy))
!write(*,*)sumxy,sumxx,sumyy
end subroutine
