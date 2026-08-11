! From the book "Numerical Recipes"
! random number (uniform deviation): between 0 and 1
!  ran3 for iris: use static vars.
      function ran3(idum)
      integer idum
      integer mbig,mseed,mz
      double precision ran3,fac
!         parameter (mbig=4000000.,mseed=1618033.,mz=0.,fac=2.5e-7)
      parameter (mbig=1000000000,mseed=161803398,mz=0,fac=1.d-9)
      integer mj,mk,ma(55)
      integer i,iff,ii,inext,inextp,k
      save iff,inext,inextp,ma
      data iff /0/
1      if(idum.lt.0.or.iff.eq.0)then
        iff=1
        mj=mseed-iabs(idum)
        mj=mod(mj,mbig)
        ma(55)=mj
        mk=1
        do 11 i=1,54
          ii=mod(21*i,55)
          ma(ii)=mk
          mk=mj-mk
          if(mk.lt.mz)mk=mk+mbig
          mj=ma(ii)
11      continue
        do 13 k=1,4
          do 12 i=1,55
            ma(i)=ma(i)-ma(1+mod(i+30,55))
            if (ma(i).lt.mz) ma(i)=ma(i)+mbig
12        continue
13      continue
        inext=0
        inextp=31
        idum=1
      endif
      inext=inext+1
      if (inext.eq.56) inext=1
      inextp=inextp+1
      if (inextp.eq.56) inextp=1
      mj=ma(inext)-ma(inextp)
      if (mj.lt.mz) mj=mj+mbig
      ma(inext)=mj
      ran3=mj*fac
      if (ran3.le.0.or.ran3.ge.1) goto 1
      return
      end                              

