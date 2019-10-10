subroutine SOR(M_old,R,m,n,beta,dx,ww)

implicit none

   real(8)                            :: max_error,beta,ww,dx
   integer                            :: m,n,i,j,iteration
   real                               :: start, finish, t
   real(8),dimension(m,n)             :: M_old,M_new, abs_error, R

M_new=M_old
iteration=0
max_error=1


do while (max_error .GE. 1e-9)
          do i=2,m-1
            do j=2,n-1
                M_new(i,j)= (1-ww)*M_old(i,j)&
                + ww/(2*(1+beta**2))*(M_old(i+1,j)+M_new(i-1,j) + beta**2*(M_old(i,j+1)+M_new(i,j-1)) -dx**2*R(i,j))

                abs_error(i,j)=abs(M_new(i,j)-M_old(i,j))
            enddo
          enddo
          max_error=maxval(abs_error)
          M_old=M_new
          iteration=iteration+1
enddo

          print*, 'P residual= ' , max_error
          print*, 'P #iterations= ', iteration

end subroutine
