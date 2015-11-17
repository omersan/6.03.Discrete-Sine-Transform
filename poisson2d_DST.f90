!-----------------------------------------------------------------------------!
! MAE 5093 DOCUMENTATION | Engineering Numerical Analysis
!-----------------------------------------------------------------------------!
! >>> Discrete Sine Transform for solving 2D Poisson equation
!     d2u/dx2 + d2u/dy2 = f(x,y)
!     Drichlet b.c.
!
!-----------------------------------------------------------------------------!
! References: 
! * Fundamentals of Engineering Numerical Analysis by P. Moin (2012)
! * Numerical Recipes: The Art of Scientific Computing, 2nd Edition (1992) 
!-----------------------------------------------------------------------------!
! Written by Omer San
!            CFDLab, Oklahoma State University, cfdlab.osu@gmail.com
!            www.cfdlab.org
! 
! Last updated: Nov. 17, 2015
!-----------------------------------------------------------------------------!

program poisson2d
implicit none
integer::i,j,nx,ny
real*8,dimension(:,:),allocatable ::u,f,ue,e
real*8,dimension(:),allocatable ::x,y
real*8 ::dx,dy,x0,xL,y0,yL

!Domain
x0 =-1.0d0 !left
xL = 1.0d0 !right

y0 =-1.0d0 !bottom
yL = 1.0d0 !up

!number of points
nx = 64  !number of grid points in x (i.e., should be power of 2)
ny = nx  !number of grid points in y

!grid spacing (spatial)
dx = (xL-x0)/dfloat(nx)
dy = (yL-y0)/dfloat(ny)

!spatial coordinates 
allocate(x(0:nx))
do i=0,nx
x(i) = x0 + dfloat(i)*dx
end do

allocate(y(0:ny))
do j=0,ny
y(j) = y0 + dfloat(j)*dy
end do

allocate(u(0:nx,0:ny))
allocate(f(0:nx,0:ny))
allocate(e(0:nx,0:ny))
allocate(ue(0:nx,0:ny))

!---------------------------------------------!
!Exact solution (test case from Moin's textbook):
!---------------------------------------------!
do j=0,ny
do i=0,nx
f(i,j) =-2.0d0*(2.0d0-x(i)*x(i)-y(j)*y(j))
ue(i,j)= (x(i)*x(i)-1.0d0)*(y(j)*y(j)-1.0d0)
end do
end do

!boundary conditions
do j=0,ny
u(0,j) = 0.0d0
u(nx,j) = 0.0d0
end do
do i=0,nx
u(i,0) = 0.0d0
u(i,ny) = 0.0d0
end do

!----------------------!
!Solver:
!----------------------!
call DST(nx,ny,dx,dy,f,u) !Discrete Sine Transform


!----------------------!
!Error analysis:
!----------------------!
do i=0,nx
do j=0,ny
e(i,j) = dabs(u(i,j)-ue(i,j))
end do 
end do

!maximum norm
write(*,*)"Max-norm =",maxval(e)


!Plot field
open(10,file='field.plt')
write(10,*) 'variables ="x","y","f","u","ue"'
write(10,*)'zone f=point i=',nx+1,',j=',ny+1
do j=0,ny
do i=0,nx
write(10,*) x(i),y(j),f(i,j),u(i,j),ue(i,j)
end do
end do
close(10)

end


!---------------------------------------------------------------------------!
!Discrete Sine Transform for Poisson Equation in 2D
!---------------------------------------------------------------------------!
subroutine DST(nx,ny,dx,dy,f,u) 
implicit none
integer::i,j,k,l,nx,ny
real*8,dimension(0:nx,0:ny)::u,f
real*8,dimension(:,:),allocatable:: ft,ut
real*8::dx,dy,pi,alpha

pi=4.0d0*datan(1.0d0)

allocate(ft(0:nx,0:ny))
allocate(ut(0:nx,0:ny))

!1. fast inverse discrete sine transform of source term:
do k=1,nx-1
do l=1,ny-1
  
    !sum
    ft(k,l) = 0.0d0
	do i=1,nx-1
	do j=1,ny-1
    ft(k,l) = ft(k,l) &
            + f(i,j)*dsin(pi*dfloat(k)*dfloat(i)/dfloat(nx)) &
                    *dsin(pi*dfloat(l)*dfloat(j)/dfloat(ny))
	end do
	end do
    
end do
end do
  
!normalize    
do k=1,nx-1
do l=1,ny-1
ft(k,l) = ft(k,l)*(2.0d0/dfloat(nx))*(2.0d0/dfloat(ny))
end do
end do
  
!2. Compute fourier coefficient of solution u:
do k=1,nx-1
do l=1,ny-1
alpha=2.0d0/(dx*dx)*(dcos(pi*dfloat(k)/dfloat(nx))-1.0d0) &
     +2.0d0/(dy*dy)*(dcos(pi*dfloat(l)/dfloat(ny))-1.0d0)

ut(k,l)=ft(k,l)/alpha
end do
end do

!3. fast forward fourier sine transform to find u:
do i=1,nx-1
do j=1,ny-1
  
    !sum
    u(i,j) = 0.0d0
	do k=1,nx-1
	do l=1,ny-1
    u(i,j) = u(i,j) &
            +ut(k,l)*dsin(pi*dfloat(k)*dfloat(i)/dfloat(nx)) &
                    *dsin(pi*dfloat(l)*dfloat(j)/dfloat(ny))
	end do
	end do

end do
end do

deallocate(ft)
deallocate(ut)

return
end


