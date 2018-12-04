
       subroutine  on off slow binaries
       include 'ARCparams.CH'
       ! common/slowi/aat(nmx),rix(nmx),rsch(nmx),OK(nmx,-1:1)
       ! common/rajahairio/tiny_pertu
        real*8 va(nmx3),vc(nmx3)
c       transformation of the ks-variables to the physical ones.
c       first transform to chain coordinates.
c        'obtain physical chain variables.'
       l=3*(n-2)
       do k=1,3
       va(k)=-wc(k)/mc(1)
       va(l+k+3)=wc(l+k)/mc(n)
       end do
        x(1)=-xc(1)
        x(2)=-xc(2)
        x(3)=-xc(3)
        do i=4,3*n
        x(i)=x(i-3)+xc(i-3)
        end do
       do i=2,n
       l=3*(i-1)
       do k=1,3
       va(l+k)=(wc(l+k-3)-wc(l+k))/mc(i)
       end do
       end do
         do j=1,3*(n-1)
         vc(j)= va(j+3)-va(j)
         end do

c        decide about slow-down binaries
c--------begin--------------------------------
          do k=1,n-1
          rix(k)=0
          do i=1,n
          if( (i.ne.k).and.(i.ne.(k+1)) )then
          lk=3*(k-1)
          lk1=lk+3
          li=3*(i-1)
        rk = sqrt((x(lk+1)-x(li+1))**2+(x(lk+2)-x(li+2))**2+
     & (x(lk+3)-x(li+3))**2)
        rk1= sqrt((x(lk1+1)-x(li+1))**2+(x(lk1+2)-x(li+2))**2+
     & (x(lk1+3)-x(li+3))**2)
          rix(k)=rix(k)+mc(i)*(1./rk**3+1./rk1**3)
          end if
          end do
          end do
c---------end-----------------------------------

         do i=1,n-1
         l=3*(i-1)
         r2=xc(l+1)**2+xc(l+2)**2+xc(l+3)**2
         w2=vc(l+1)**2+vc(l+2)**2+vc(l+3)**2
         r1=sqrt(r2)
         mi1=mc(i)+mc(i+1)
                     aat(i)=2.d0/r1-w2/mi1
                     perta=8./mi1*rix(i)/aat(i)!!!**3
                     pert0=tiny_pertu
                     yrwi=sqrt(abs(perta)/(pert0+1.d-33))/aat(i)
                     rwi=max(1.d0,1/yrwi)
         if(r1*aat(i).gt.0.3)then
          if(abs(rsch(i)/rwi-1).gt.0.1
     &  .or.(rwi.eq.1.0.and.rsch(i).ne.1.))then
          rsch(i)=rwi
          end if
         end if
         OK(I,0)=1.d0/rsch(i)
         if(I.gt.2)then
         OK(I,-1)=(1-OK(I,0))*MC(I-2)/(MC(I-1)+MC(I-2))
         else
         OK(I,-1)=0
         end if
         if(I.le.N-3)then
         OK(I,+1)=(1-OK(I,0))*MC(I+2)/(MC(I+1)+MC(I+2)) 
         else
         OK(I,+1)=0
         end if
         end do
c-------auxiliary variables
        iw=0
        do ik=1,N-1
        if(rsch(ik).gt.1.0)iw=iw+1
        end do

       if(iw.gt.0.)then
       k6=69
       !write(k6,*)' SLOW ',(rsch(k),k=1,N-1)
       !write(k6,*) 'AIKA ',iw,Taika+chtime,(iname(k),k=1,N)
        do K=1,N-1
        if(OK(K,0).ne.1.) 
     & write(k6,116)K,OK(K,-1),OK(K,0),OK(K,1),Taika+chtime,tiny_pertu
116     format(1x,I3,3f10.6,5x,f12.2,1p,g11.2)
        end do
       write(k6,*)' -------- '
        end if
       return
       end

