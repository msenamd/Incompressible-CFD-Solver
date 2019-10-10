program Incomp_NS

implicit none

!!!!!!-------variables declaration--------!!!!!
   integer                            :: m,n,i,j, iteration , counter ,anim_freq
   real(8)                            :: Lx,Ly,dx,dy,beta, v_w , Ren, cycles ,relaxation
   real(4)                            :: time, t_final, dt , CFL
   real(8),allocatable,dimension(:,:) :: P, u_f, u_b, v_f , v_b , u , v , R , psi , vor
   real(8),allocatable,dimension(:,:) :: F_f, F_b , G_f, G_b
   real(8),allocatable,dimension(:,:) :: uv_ff , uv_fb , uv_bf , uv_bb
   real(8),allocatable,dimension(:,:) :: P_old,P_new, difference
   real(8)                            :: max_error , u_center
   character(len=50)                  :: filename
   real(8)                            :: mm,nn


!!!!!-------------Inputs-------------!!!!!
cycles=100000
anim_freq=500
Lx=1.0d0
Ly=1.0d0
v_w=1.0d0
m=240
n=240
dx=Lx/(m-2)
dy=Ly/(n-2)
beta=dx/dy
time=0.0d0
Ren=10000.0d0
dt=0.0025*Ren*dx**2
t_final=cycles*dt
CFL=dt/dx
relaxation=1.3
mm=m
nn=n

!!!!!-------Creating variables--------!!!!!

allocate(P(m,n) , u(m,n) , v(m,n) , R(m,n) , psi(m,n) , vor(m,n))
allocate(u_f(m,n) , u_b(m,n) , v_f(m,n) , v_b(m,n))
allocate(uv_ff(m,n) , uv_fb(m,n) , uv_bf(m,n) , uv_bb(m,n))
allocate(F_f(m,n), F_b(m,n), G_f(m,n), G_b(m,n))
allocate(P_new(m,n),P_old(m,n),difference(m,n))

P=0.0d0
u=0.0d0
v=0.0d0
R=0.0d0
psi=0.0d0
vor=0.0d0

uv_ff=0.0d0
uv_fb=0.0d0
uv_bf=0.0d0
uv_bb=0.0d0

u_f=0.0d0
u_b=0.0d0
v_f=0.0d0
v_b=0.0d0

F_f=0.0d0
F_b=0.0d0
G_f=0.0d0
G_b=0.0d0

P_old=0.0d0
P_new=0.0d0


!!!!!-------Update Boundary conditions--------!!!!!
    do i=1,m
       !bottom
       v_f(i,1)=0.0d0
       u_f(i,1)= -u_f(i,2)
       P(i,1)=P(i,2)

       !top
       v_f(i,n-1)=0.0d0
       u_f(i,n)=2*v_w-u_f(i,n-1)
       P(i,n)=P(i,n-1)
    end do

    do j=1,n
        !left
        u_f(1,j)=0.0d0
        v_f(1,j)= -v_f(2,j)
        P(1,j)=P(2,j)

        !right
        u_f(m-1,j)=0.0d0
        v_f(m,j)= -v_f(m-1,j)
        P(m,j)=P(m-1,j)
    end do

    do i=2,m
        do j=1,n
            u_b(i,j)=u_f(i-1,j)
        end do
    end do

    do i=1,m
        do j=2,n
            v_b(i,j)=v_f(i,j-1)
        end do
    end do

    do i=1,m
        do j=1,n
            u(i,j)=0.5*(u_f(i,j)+u_b(i,j))
            v(i,j)=0.5*(v_f(i,j)+v_b(i,j))
        end do
    end do


!!!!!----------------------------Main Solver------------------------------!!!!!

counter=0

      write (filename,'(A,I6.6,A)') "p", counter,".csv"
      print*,filename

      open(unit=100,file=filename)
      do i=2,m-1
         write(100,*) (P(i,j),',', j = 2,n-2), P(i,n-1)
      enddo
      close(unit=100)

      write (filename,'(A,I6.6,A)') "u", counter,".csv"
      print*,filename
      open(unit=200,file=filename)
      do i=2,m-1
         write(200,*) (u(i,j),',', j = 2,n-2), u(i,n-1)
      enddo
      close(unit=200)

      write (filename,'(A,I6.6,A)') "v", counter,".csv"
      print*,filename
      open(unit=300,file=filename)
      do i=2,m-1
         write(300,*) (v(i,j),',', j = 2,n-2), v(i,n-1)
      enddo
      close(unit=300)

open(unit=10,file="Time_U.csv")


do counter=1,cycles

    time=time+dt
    u_center=u(nint(mm/2),nint(nn/4))
    print*, '*********************************************'
    print*,'cycle number =    ', counter
    print*, 'dt= ',dt
    print*, 'dx= ',dx
    print*,'current time =    ', time
    print*,'Courant number =    ', CFL
    print*,'monitoring velocity=  ', u_center

    do i=2,m-1
        do j=2,n-1
            uv_ff(i,j)=0.5*(u_f(i,j)+u_f(i,j+1)) * 0.5*(v_f(i,j)+v_f(i+1,j))
            uv_fb(i,j)=0.5*(u_f(i,j)+u_f(i,j-1)) * 0.5*(v_b(i,j)+v_b(i+1,j))
            uv_bf(i,j)=0.5*(u_b(i,j)+u_b(i,j+1)) * 0.5*(v_f(i,j)+v_f(i-1,j))
            uv_bb(i,j)=0.5*(u_b(i,j-1)+u_b(i,j)) * 0.5*(v_b(i,j)+v_b(i-1,j))
        end do
    end do

! Intermediate step (F,G)
    do i=2,m-1
        do j=2,n-1
            F_f(i,j)=u_f(i,j) +(dt/(Ren*dx**2))*(u_f(i+1,j)-2*u_f(i,j)+u_b(i,j)) &
                     +(dt/(Ren*dy**2))*(u_f(i,j-1)-2*u_f(i,j)+u_f(i,j+1)) &
                     -(dt/dx)*(u(i+1,j)**2 - u(i,j)**2) &
                     -(dt/dy)*(uv_ff(i,j)-uv_fb(i,j))

            G_f(i,j)=v_f(i,j) +(dt/(Ren*dy**2))*(v_f(i+1,j)-2*v_f(i,j)+v_f(i-1,j)) &
                     +(dt/(Ren*dy**2))*(v_f(i,j+1)-2*v_f(i,j)+v_b(i,j)) &
                     -(dt/dx)*(uv_ff(i,j) - uv_bf(i,j)) &
                     -(dt/dy)*(v(i,j+1)**2 - v(i,j)**2)

            F_b(i,j)=u_b(i,j) +(dt/(Ren*dx**2))*(u_f(i,j)-2*u_b(i,j)+u_b(i-1,j)) &
                     +(dt/(Ren*dy**2))*(u_b(i,j-1)-2*u_b(i,j)+u_b(i,j+1)) &
                     -(dt/dx)*(u(i,j)**2 - u(i-1,j)**2) &
                     -(dt/dy)*(uv_bf(i,j)-uv_bb(i,j))

            G_b(i,j)=v_b(i,j) +(dt/(Ren*dy**2))*(v_b(i+1,j)-2*v_b(i,j)+v_b(i-1,j)) &
                     +(dt/(Ren*dy**2))*(v_f(i,j)-2*v_b(i,j)+v_b(i,j-1)) &
                     -(dt/dx)*(uv_fb(i,j) - uv_bb(i,j)) &
                     -(dt/dy)*(v(i,j)**2 - v(i,j-1)**2)

        end do
    end do


! Pressure equation
    do i=1,m
        do j=1,n
            R(i,j)=( (F_f(i,j)-F_b(i,j))/dx + (G_f(i,j)-G_b(i,j))/dy ) / dt
        end do
    end do

    !call SOR(P,R,m,n,beta,dx,relaxation)
    call ADI(P,R,m,n,beta,dx,relaxation)

! Momentum equation
    do i=2,m-1
        do j=2,n-1
            u_f(i,j)=F_f(i,j) - dt/dx*(P(i+1,j) - P(i,j))
            v_f(i,j)=G_f(i,j) - dt/dy*(P(i,j+1) - P(i,j))
        end do
    end do


!!!!!-------Update Boundary conditions--------!!!!!
    do i=1,m
       !bottom
       v_f(i,1)=0.0d0
       u_f(i,1)= -u_f(i,2)
       P(i,1)=P(i,2)

       !top
       v_f(i,n-1)=0.0d0
       u_f(i,n)=2*v_w-u_f(i,n-1)
       P(i,n)=P(i,n-1)
    end do

    do j=1,n
        !left
        u_f(1,j)=0.0d0
        v_f(1,j)= -v_f(2,j)
        P(1,j)=P(2,j)

        !right
        u_f(m-1,j)=0.0d0
        v_f(m,j)= -v_f(m-1,j)
        P(m,j)=P(m-1,j)
    end do

    do i=2,m
        do j=1,n
            u_b(i,j)=u_f(i-1,j)
        end do
    end do

    do i=1,m
        do j=2,n
            v_b(i,j)=v_f(i,j-1)
        end do
    end do

    do i=1,m
        do j=1,n
            u(i,j)=0.5*(u_f(i,j)+u_b(i,j))
            v(i,j)=0.5*(v_f(i,j)+v_b(i,j))
        end do
    end do


!vorticity calculation

    do i=1,m-1
        do j=1,n-1
            psi(i,j)=psi(i-1,j)-dx*(v(i,j))
            vor(i,j)=(v(i+1,j)-v(i,j))/(dx) - (u(i,j+1)-u(i,j))/(dy)
        end do
    end do


!Animation output
   if(mod(counter,anim_freq)==0) then
      write (filename,'(A,I6.6,A)') "p", counter,".csv"
      print*,filename

      open(unit=100,file=filename)
      do i=2,m-1
         write(100,*) (P(i,j),',', j = 2,n-2), P(i,n-1)
      enddo
      close(unit=100)

      write (filename,'(A,I6.6,A)') "u", counter,".csv"
      print*,filename
      open(unit=200,file=filename)
      do i=2,m-1
         write(200,*) (u(i,j),',', j = 2,n-2), u(i,n-1)
      enddo
      close(unit=200)

      write (filename,'(A,I6.6,A)') "v", counter,".csv"
      print*,filename
      open(unit=300,file=filename)
      do i=2,m-1
         write(300,*) (v(i,j),',', j = 2,n-2), v(i,n-1)
      enddo
      close(unit=300)

   endif


write(10,*) time,',', u_center

end do

close(unit=10)


!!!!!-------Output results--------!!!!

open(unit=1, file="p.csv")
do i=2,m-1
write(1,*) (P(i,j),',', j = 2,n-2), P(i,n-1)
enddo
close(unit=1)

open(unit=2, file="u.csv")
do i=2,m-1
write(2,*) (u(i,j),',', j = 2,n-2), u(i,n-1)
enddo
close(unit=2)

open(unit=3, file="v.csv")
do i=2,m-1
write(3,*) (v(i,j),',', j = 2,n-2), v(i,n-1)
enddo
close(unit=3)

open(unit=4, file="stream.csv")
do i=2,m-1
write(4,*) (psi(i,j),',', j = 2,n-2), psi(i,n-1)
enddo
close(unit=4)

open(unit=5, file="vorticity.csv")
do i=2,m-1
write(5,*) (vor(i,j),',', j = 2,n-2), vor(i,n-1)
enddo
close(unit=5)

end program Incomp_NS
