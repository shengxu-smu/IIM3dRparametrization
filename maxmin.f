c-----------------------------------------------------------------------
c

      subroutine maxmin(l1,l2,l3,lu,ld)
      implicit none
      integer l1,l2,l3,lu,ld

      if(l1.lt.l2) then
        lu=l2
        ld=l1
      else
        lu=l1
        ld=l2
      endif
      if(l3.gt.lu) then
        lu=l3
      elseif(l3.lt.ld) then
        ld=l3
      endif

      return
      end


c-----------------------------------------------------------------------

