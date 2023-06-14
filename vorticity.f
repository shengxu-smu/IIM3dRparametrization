c-----------------------------------------------------------------------
c
      subroutine vorticity
      include 'parameter.inc'
      include 'field.inc'
      real*8 u11,u12,u13,v21,v22,v23,w31,w32,w33,rr,qq,discriminant

      call u_interpolate
      call u_strain
      call v_interpolate
      call v_strain
      call w_interpolate
      call w_strain

      do k=1,nz
        do j=1,ny
          do i=1,nx
            u11=data1(i,j,k)
            u12=data2(i,j,k)
            u13=data3(i,j,k)
            v21=data4(i,j,k)
            v22=data5(i,j,k)
            v23=data6(i,j,k)
            w31=data7(i,j,k)
            w32=data8(i,j,k)
            w33=data9(i,j,k)
            qq=-0.5d0*(u11*u11+v22*v22+w33*w33+
     .          2.0d0*(u12*v21+u13*w31+v23*w32))
            rr=u11*v22*w33+u12*v23*w31+u13*v21*w32-
     .         u11*v23*w32-u12*v21*w33-u13*v22*w31
            discriminant=(rr/2.0d0)**2.0d0+(qq/3.0d0)**3.0d0
            if(discriminant.gt.0.0d0) then
              o(i,j,k)=qq
            else
              o(i,j,k)=0.0d0
            endif
          enddo
        enddo
      enddo

      do k=1,nz
        do j=1,ny
          do i=1,nx
            fecc(i,j,k)=data8(i,j,k)-data6(i,j,k)
            fcec(i,j,k)=data3(i,j,k)-data7(i,j,k)
            fcce(i,j,k)=data4(i,j,k)-data2(i,j,k)
          enddo
        enddo
      enddo

      call field_interpolate

      return
      end


c-----------------------------------------------------------------------



