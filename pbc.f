c-----------------------------------------------------------------------
c
      subroutine pbc
      include 'parameter.inc'
      include 'field.inc'

      call field_interpolate

      do j=1,ny
        do i=1,nx
          pb(i,j)=Re1*dz2*(data7(i,j,2)-data7(i,j,1)
     .            +dz*((data3(i,j,0)-data3(i-1,j,0))/dx+
     .                 (data6(i,j,0)-data6(i,j-1,0))/dy))
          pt(i,j)=Re1*dz2*(data7(i,j,nz-1)-data7(i,j,nz)
     .            -dz*((data3(i,j,nz)-data3(i-1,j,nz))/dx+
     .                 (data6(i,j,nz)-data6(i,j-1,nz))/dy))
        enddo
      enddo

      do k=1,nz
        do j=1,ny
          pw(j,k)=Re1*dx2*(data1(2,j,k)-data1(1,j,k)
     .            +dx*((data5(0,j,k)-data5(0,j-1,k))/dy+
     .                 (data8(0,j,k)-data8(0,j,k-1))/dz))
          pe(j,k)=Re1*dx2*(data1(nx-1,j,k)-data1(nx,j,k)
     .            -dx*((data5(nx,j,k)-data5(nx,j-1,k))/dy+
     .                 (data8(nx,j,k)-data8(nx,j,k-1))/dz))
        enddo
        do i=1,nx
          ps(i,k)=Re1*dy2*(data4(i,2,k)-data4(i,1,k)
     .            +dy*((data2(i,0,k)-data2(i-1,0,k))/dx+
     .                 (data9(i,0,k)-data9(i,0,k-1))/dz))
          pn(i,k)=Re1*dy2*(data4(i,ny-1,k)-data4(j,ny,k)
     .            -dy*((data2(i,ny,k)-data2(i-1,ny,k))/dx+
     .                 (data9(i,ny,k)-data9(i,ny,k-1))/dz))
        enddo
      enddo

      return
      end


c-----------------------------------------------------------------------
