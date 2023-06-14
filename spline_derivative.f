c-----------------------------------------------------------------------
c
      subroutine spline_derivative(falfa,dfalfa,isphr,i12,nalfa)
      include 'parameter.inc'
      include 'surface.inc'
      integer n1,n2,nd,nh,isphr,i12,nalfa
      parameter(n1=ns1,n2=ns2,nd=2*ns1+2,nh=ns2/2-1)
      real*8 falfa(0:n1,0:n2),dfalfa(0:n1,0:n2),cfalfa(0:n1,0:n2)
      real*8 aalfa(0:nd),ffalfa(0:nd,0:nh),cffalfa(0:nd,0:nh)
      real*8 f(0:n1),f1(0:n1),cf(0:n1),cf1(0:n1),df(0:n1),df1(0:n1)
      real*8 b(0:nh),b1(0:nh),cb(0:nh),cb1(0:nh),db(0:nh),db1(0:nh)
      real*8 g(0:n2),g1(0:n2),cg(0:n2),cg1(0:n2),dg(0:n2),dg1(0:n2)
      
      if(i12.eq.2) then

        call cubic_spline2(alfa20,falfa,cfalfa,n2,n2,n1+1)
        
        do m2=0,n2-1
          do m1=0,n1
            f(m1)=falfa(m1,m2)
            f1(m1)=falfa(m1,m2+1)
            cf(m1)=cfalfa(m1,m2)
            cf1(m1)=cfalfa(m1,m2+1)
          enddo
          call fdf(alfa20(m2),alfa20(m2+1),f,f1,cf,cf1,
     .             alfa20(m2),df,df1,n1+1)
          do m1=0,n1
            dfalfa(m1,m2)=df1(m1)
          enddo
        enddo

      endif

      if(i12.eq.1) then

        if(isphr.eq.1) then

          do m1=0,n1
            aalfa(m1)=alfa11(m1)
            aalfa(n1+1+m1)=2.0d0*sl1-alfa11(n1-m1)
          enddo
          aalfa(nd)=aalfa(0)+2.0d0*sl1

          do m2=0,nh
            do m1=0,n1
              ffalfa(m1,m2)=falfa(m1,m2)
            enddo            
            if(mod(nalfa,2).eq.1) then
              do m1=n1+1,nd-1
                ffalfa(m1,m2)=falfa(nd-1-m1,m2+nh+1)
              enddo
            else
              do m1=n1+1,nd-1
                ffalfa(m1,m2)=-falfa(nd-1-m1,m2+nh+1)
              enddo
            endif
            ffalfa(nd,m2)=ffalfa(0,m2)
          enddo

          call cubic_spline1(aalfa,ffalfa,cffalfa,nd,nd,nh+1)

          do m1=0,nd-1
            do m2=0,nh
              b(m2)=ffalfa(m1,m2)
              b1(m2)=ffalfa(m1+1,m2)
              cb(m2)=cffalfa(m1,m2)
              cb1(m2)=cffalfa(m1+1,m2)
            enddo
            call fdf(aalfa(m1),aalfa(m1+1),b,b1,cb,cb1,
     .               aalfa(m1),db,db1,nh+1)
            if(m1.le.n1) then
              do m2=0,nh
                dfalfa(m1,m2)=db1(m2)
              enddo
            else
              if(mod(nalfa,2).eq.1) then
                do m2=0,nh
                  dfalfa(nd-1-m1,m2+nh+1)=-db1(m2)
                enddo
              else
                do m2=0,nh
                  dfalfa(nd-1-m1,m2+nh+1)=db1(m2)
                enddo
              endif
            endif
          enddo

        endif

        if(isphr.eq.0) then

          call cubic_spline1(alfa10,falfa,cfalfa,n1,n1,n2+1)

          do m1=0,n1-1
            do m2=0,n2
              g(m2)=falfa(m1,m2)
              g1(m2)=falfa(m1+1,m2)
              cg(m2)=cfalfa(m1,m2)
              cg1(m2)=cfalfa(m1+1,m2)
            enddo
            call fdf(alfa10(m1),alfa10(m1+1),g,g1,cg,cg1,
     .               alfa10(m1),dg,dg1,n2+1)
            do m2=0,n2
              dfalfa(m1,m2)=dg1(m2)
            enddo
          enddo
          do m2=0,n2
            dfalfa(n1,m2)=dfalfa(0,m2)
          enddo

        endif
      endif

      do m1=0,n1
        dfalfa(m1,n2)=dfalfa(m1,0)
      enddo

      return
      end


c-----------------------------------------------------------------------
