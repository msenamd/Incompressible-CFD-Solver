subroutine ADI(M_old,R,m,n,beta,dx,ww)

use tridiag

implicit none

   real(8)                            :: max_error,beta,ww,dx
   integer                            :: m,n,i,j,iteration
   real                               :: start, finish, t
   real(8),dimension(m,n)             :: M_old,M_bet,M_new, difference ,R
   real(8),dimension(m)               :: a,b,c,d
   real(8),dimension(n)               :: aj,bj,cj,dj


d(1)=1.0d0
a(1)=0.0d0
b(1)=0.0d0
d(m)=1.0d0
a(m)=0.0d0
b(m)=0.0d0
do i=2,m-1
    a(i)=-ww
    b(i)=-ww
    d(i)=2*(1+beta**2)
end do

dj(1)=1.0d0
aj(1)=0.0d0
bj(1)=0.0d0
dj(n)=1.0d0
aj(n)=0.0d0
bj(n)=0.0d0
do j=2,n-1
    aj(j)=-ww*beta**2
    bj(j)=-ww*beta**2
    dj(j)=2*(1+beta**2)
end do


M_new=M_old
M_bet=M_old
iteration=0
max_error=1


do while (max_error .GE. 1e-6)

          do j=2,n-1
             c(1)=M_old(1,j)
             c(m)=M_old(m,j)

             do i=2,m-1
                c(i)=2*(1+beta**2)*(1-ww)*M_old(i,j) + ww*(beta**2)*(M_old(i,j+1)+M_bet(i,j-1)) - ww*R(i,j)*dx**2
             enddo

             M_bet(:,j)=Thomas(a,b,c,d)
          enddo

         do i=2,m-1
             cj(1)=M_old(i,1)
             cj(n)=M_old(i,n)

             do j=2,n-1
                cj(j)=2*(1+beta**2)*(1-ww)*M_bet(i,j) + ww*(M_bet(i+1,j)+M_new(i-1,j)) - ww*R(i,j)*dx**2
             enddo

             M_new(i,:)=Thomas(aj,bj,cj,dj)
          enddo

          difference=M_new-M_old
          max_error=maxval(abs(difference))
          M_old=M_new
          M_bet=M_new
          iteration=iteration+1


enddo
          print*, 'P residual= ' , max_error
          print*, 'P #iterations= ', iteration

end subroutine



