module tridiag
    implicit none
    contains


function Thomas(a,b,c,d)

   real(8),dimension(:)               :: a,b,c,d
   real(8),allocatable,dimension(:)   :: aa,bb,cc,dd
   integer                            :: jj,kk,nn
   real(8),dimension(:),allocatable   ::Thomas,uu
    nn=size(a)
    allocate(Thomas(nn),uu(nn),aa(nn),bb(nn),cc(nn),dd(nn))

aa=a
bb=b
cc=c
dd=d


do jj=2,nn
    dd(jj)=dd(jj)-bb(jj)/dd(jj-1)*aa(jj-1)
    cc(jj)=cc(jj)-bb(jj)/dd(jj-1)*cc(jj-1)
enddo

uu(nn)=cc(nn)/dd(nn)
do kk=(nn-1),1,-1
    uu(kk)=(cc(kk)-aa(kk)*uu(kk+1))/dd(kk)
enddo

Thomas=uu
end function Thomas


end module tridiag

