c-----------------------------------------------------------------------
c
c     Determine if a point (0,0) is inside a polygon defined by np points
c     pnt(2,np). Point (0,0) on an edge is regarded outside.
c
      subroutine point_polygon(np,pnt,inside)
      implicit none
      integer i,np,inside,nrcount,nlcount
      real*8 x0,y0,x1,y1,x2,y2,xi
      real*8 pnt(2,np)

      write(*,*)'  !! note: point-inside-polygon test !!'

      inside=0
      nlcount=0
      nrcount=0
     
      x0=0.0d0
      y0=0.0d0
      x1=pnt(1,1)
      y1=pnt(2,1)

      do i=2,np+1
        if(i.eq.np+1) then
          x2=pnt(1,1)
          y2=pnt(2,1)
        else
          x2=pnt(1,i)
          y2=pnt(2,i)
        endif
        if(y0.gt.min(y1,y2)) then
          if(y0.le.max(y1,y2)) then
            if(x0.le.max(x1,x2)) then
              if(y1.ne.y2) then
                xi=x1+(x2-x1)*(y0-y1)/(y2-y1)
                if(x0.lt.xi) nrcount=nrcount+1
              endif
            endif
            if(x0.ge.min(x1,x2)) then
              if(y1.ne.y2) then
                xi=x1+(x2-x1)*(y0-y1)/(y2-y1)
                if(x0.gt.xi) nlcount=nlcount+1
              endif
            endif
          endif
        endif
        x1=x2
        y1=y2
      enddo

      if(mod(nrcount,2).eq.1.and.mod(nlcount,2).eq.1) inside=1

      return
      end


c-----------------------------------------------------------------------
