c-----------------------------------------------------------------------
c
c     Project vertices on plane defined by indices (m,n); Determine whether 
c     vertex pt3 and point (pp,qq) are at the same side of line formed by 
c     vertices pt1 and pt2 on the plane.
c     Note: The order of pt1 and pt2 matters due to roundoff errors!
c
      subroutine line_side(m,n,pt1,pt2,pt3,pp,qq,iflag)
      implicit none
      integer m,n,iflag
      real*8 pt1(3),pt2(3),pt3(3)
      real*8 pp,qq,p1,q1,p2,q2,p3,q3,p21,q21,p30,q30,pp0,qq0
 
      p1=pt1(m)
      q1=pt1(n)
      p2=pt2(m)
      q2=pt2(n)
      p3=pt3(m)
      q3=pt3(n)
      p21=p2-p1
      q21=q2-q1
      if(p21.eq.0.0d0.and.q21.eq.0.0d0) then
        write(*,*)'  !! erro: vertices projected to the same point !!'
        stop
      endif

      iflag=-1
      if(abs(p21).ge.abs(q21)) then
        q30=q1+q21*(p3-p1)/p21
        qq0=q1+q21*(pp-p1)/p21
        if(((q3-q30)*(qq-qq0)).gt.0.0d0) iflag=1
        if((qq-qq0).eq.0.0d0) iflag=0
      else
        p30=p1+p21*(q3-q1)/q21
        pp0=p1+p21*(qq-q1)/q21
        if(((p3-p30)*(pp-pp0)).gt.0.0d0) iflag=1
        if((pp-pp0).eq.0.0d0) iflag=0
      endif
        
      return
      end


c-----------------------------------------------------------------------

