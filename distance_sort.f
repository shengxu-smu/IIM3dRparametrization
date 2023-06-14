c-----------------------------------------------------------------------
c
c     insertion sort: sort distance from the smallest to the largest
c
      subroutine distance_sort(np,nn,dist,iorder)
      implicit none
      integer np,nn,n,i,ifoo,iorder(nn)
      real*8 foo,dist(nn)
 
      do n=1,nn
        iorder(n)=n
      enddo

      do n=2,np
        i=n-1
        foo=dist(n)
        ifoo=iorder(n)
 1      continue
        if(i.gt.0.and.dist(i).gt.foo) then
          dist(i+1)=dist(i)
          iorder(i+1)=iorder(i)
          i=i-1
          goto 1
        else
          dist(i+1)=foo
          iorder(i+1)=ifoo
        endif
      enddo

      return
      end


c-----------------------------------------------------------------------

