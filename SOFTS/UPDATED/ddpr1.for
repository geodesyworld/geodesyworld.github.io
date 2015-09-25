      program ddpr1

*** pseudo-range version  (dgm   1-mar-96)
*** nothing with L2 pr (so far)

*** version (0.1  test initiated 11-jan-95)
*** fixed nasty bug in pscan()
*** note to me: work in cutoff-angle
*** note to me: if orbits have NO sat clock corr., use .nav file

      implicit double precision(a-h,o-z)
      parameter(lorbit=(4*48)*(1+4*24))
      dimension orbits(lorbit),coord(20),icols1(0:4),icols2(0:4)
      character*4 aid1,aid2
      common/lus/lst,lpos,lrin1,lrin2,leph,lobs
      common/stugg/pi,rad,sol,f1,f2,w1,w2,we

      pi=4.d0*datan(1.d0)
      rad=180.d0/pi

      sol=299792458.d0
      f1 =1575.42d6
      f2 =1227.60d6
      w1 =sol/f1
      w2 =sol/f2
      we=7.2921151467d-5

*** note: input file names currently hard-coded

      lst=1
      open(lst,file='ddpr1.lst',form='formatted')
      write(lst,*) 'Pseudorange Program -- L1 only -- 2005sep17'

      lpos=10
      open(lpos,file='pos',form='formatted',status='old')

      lrin1=2
      open(lrin1,file='point1.rnx',form='formatted',status='old')

      lrin2=3
      open(lrin2,file='point2.rnx',form='formatted',status='old')

      leph=4
      open(leph,file='orbit.sp3',form='formatted',status='old')

*** gather critical information

      call recon(nsat,nepoch,tstrt,tstop,dstrt,dstop,delt,
     *           icols1,ndat1,aid1,icols2,ndat2,aid2,coord)

*** load ephemeris based on identified interval -- then close file

      rewind(leph)
      call lodorb(orbits,lorbit,nsat,nepoch,tstrt,tstop)
      close(leph)
      leph=6

*** get precise coordinates from the lpos file

      call lodpos(coord,aid1,aid2)

*** get antenna coordinates

      call getoff(coord,aid1,aid2)

      write(lst,*) 'coordinate dump'
      write(lst,'(3f14.4))') coord( 1),coord( 2),coord( 3)
      write(lst,'(3f14.4))') coord( 4),coord( 5),coord( 6)
      write(lst,'(3f14.4))') coord( 7),coord( 8),coord( 9)
      write(lst,'(3f14.4))') coord(10),coord(11),coord(12)
      write(lst,'(4f14.4))') coord(13),coord(14),coord(15),coord(16)
      write(lst,'(4f14.4))') coord(17),coord(18),coord(19),coord(20)
      basel=dsqrt( (coord( 7)-coord( 1)) **2 +
     *             (coord( 8)-coord( 2)) **2 +
     *             (coord( 9)-coord( 3)) **2   )
      write(lst,'(a,f14.4,a,f6.2,a))') 'baseline = ',basel,
     *          ' m.    interval = ',(dstop-dstrt)/3600.d0,' hr.'

*** default troposphere correction

      call inimet(.false.)
      write(lst,*) 'no troposphere model used'

***** call inimet(.true.)
***** write(lst,*) 'default troposphere model used'

*** build reference satellite table

      call makrst(dstrt,dstop,delt,coord,orbits,nsat)

*** adjust the data

      call adjst(orbits,coord,nsat,dstrt,dstop,delt,
     *           icols1,ndat1,icols2,ndat2)

*** end of processing

      write(lst,*) 'End of Processing'

      stop
      end
      subroutine adjst(orbits,coord,nsat,tstrt,tstop,delt,
     *                 icols1,ndat1,icols2,ndat2)

*** iterative adjustment routine

      implicit double precision(a-h,o-z)
      dimension orbits(*),coord(20),icols1(0:4),icols2(0:4)
      common/lus/lst,lpos,lrin1,lrin2,leph,lobs

*** read data, form observation equations, solve for first time
*** apply corrections one time only in formob()

      lobs=8
      open(lobs ,form='unformatted',status='scratch')
      call formob(orbits,coord,nsat,tstrt,tstop,delt,
     *            icols1,ndat1,icols2,ndat2,rms)

*** loop over the adjustment  (do auto-rejection)

      siglvl=2.5d0
      maxitr=4
      iter=1

10000 if(iter.le.maxitr) then
        write(*,*) ' iteration -- ',iter
        rejlvl=siglvl*rms
        call comrms(orbits,coord,rejlvl,rms,nsat)
        iter=iter+1
        go to 10000
      endif

*** provide results  (bogus past here)

      rejlvl=siglvl*rms
      lres=10
      open(lres,file='ddpr1.res',form='formatted')
      call post(orbits,coord,nsat,rejlvl,lres)
      close(lobs)

      return
      end
      subroutine post(orbits,coord,nsat,rejlvl,lres)

*** post adjustment display

      implicit double precision(a-h,o-z)
      logical getrs,xyzgeo
      parameter(maxsat=40)
      dimension nsv(maxsat),rmss(maxsat),nsv2(maxsat),rms2(maxsat)
      dimension orbits(*),coord(20),cdum(3)
      common/stugg/pi,rad,sol,f1,f2,w1,w2,we
      common/lus/lst,lpos,lrin1,lrin2,leph,lobs
      save  /xyz2o/
      common/xyz2o/x2o,y2o,z2o

*** storage scheme for coord()
***
***  1  2  3       X Y Z  point 1 marker
***  4  5  6       X Y Z  point 1 antenna, L1
***  7  8  9       X Y Z  point 2 marker       (may be approx)
*** 10 11 12       X Y Z  point 2 antenna, L1  (may be approx)

*** 13 14 15 16    dl1,dl2,dn,de  point 1 offsets
*** 17 18 19 20    dl1,dl2,dn,de  point 2 offsets

      write(lst,*) 'adjusted point 2 marker and l1 phase center'
      write(lst,'(3f14.4))') coord( 7),coord( 8),coord( 9)
      write(lst,'(3f14.4))') coord(10),coord(11),coord(12)

      if(.not.xyzgeo(coord(10),coord(11),coord(12),gla,glo,eht)) stop 25389

      dx=coord(7)-x2o
      dy=coord(8)-y2o
      dz=coord(9)-z2o
      call rge(gla,glo,dn,de,du,dx,dy,dz)
      write(lst,*) 'delta x,y,z and n,e,u at point 2 marker'
      write(lst,'(3f14.4))') dx,dy,dz
      write(lst,'(3f14.4))') dn,de,du

*** create residual file
*** loop over obs. equations -- accumulate statistics

      do 1 i=1,maxsat
      nsv(i) =0
      nsv2(i)=0
      rmss(i)=0.d0
    1 rms2(i)=0.d0

      rewind(lobs)
  500 read(lobs,end=773) iisat,trec1,trec2,r110,r210,r11,r21,rdd1
      if(iisat.le.0.or.iisat.gt.nsat) then
        continue

*** form double differences --   sr,s geometric range (m)

      else
        if(.not.getrs(iirfsv,trec1)) stop 95397
        call comdtr(iirfsv,coord,orbits,trec1,trec2,r110,r210,dtr1,dtr2)
        call comrng(iirfsv,coord,orbits,trec1,trec2,r110,r210,dtr1,dtr2,
     *              sr1,sr2,cdum)
        call comrg2(iisat ,coord,orbits,trec1,trec2,r11,r21,dtr1,dtr2,
     *              s1,s2,gla,glo,vad)
        sdd1=((s2-s1)-(sr2-sr1))
        cmo=sdd1-rdd1
        isat=invprn(iisat)
        nsv2(isat)=nsv2(isat)+1
        rms2(isat)=rms2(isat)+cmo*cmo
        write(lres,'(f9.1,f9.3,f6.2,i3)') trec1,cmo,vad,isat
        if(dabs(cmo).le.rejlvl) then
          nsv(isat)=nsv(isat)+1
          rmss(isat)=rmss(isat)+cmo*cmo
        endif
      endif
      go to 500

*** end of file encountered -- display rms after trim

  773 rewind(lobs)
      write(lst,*) 'rms per satellite    (with and without trim)'
      iobs=0
      rms=0.d0
      do 10 i=1,maxsat
        if(nsv(i).gt.0) then
          iobs=iobs+nsv(i)
          rms =rms +rmss(i)
          rmss(i)=dsqrt(rmss(i)/dble(nsv(i)))
          rms2(i)=dsqrt(rms2(i)/dble(nsv2(i)))
          write(lst,'(i2,i8,2f10.3)') i,nsv(i),rmss(i),rms2(i)

*** deleted satellite ??

        elseif(nsv2(i).gt.0) then
          rmss(i)=0.d0
          rms2(i)=dsqrt(rms2(i)/dble(nsv2(i)))
          write(lst,'(i2,i8,2f10.4,a)') i,nsv(i),rmss(i),rms2(i),' ****'
        endif
   10 continue
      rms=dsqrt(rms/dble(iobs))
      write(lst,'(a,f9.4,a,i6)') 'post-fit rms =',rms,' m.    n=',iobs

      return
      end
      subroutine comrms(orbits,coord,rejlvl,rms,nsat)

*** routine to evaluate obs. equations, form and solve normals

      implicit double precision(a-h,o-z)
      logical getrs
      dimension orbits(*),coord(20),c(3),cr(3)
      dimension ic(3),an(3,3),u(3),x(3)
      common/stugg/pi,rad,sol,f1,f2,w1,w2,we
      common/lus/lst,lpos,lrin1,lrin2,leph,lobs

*** initialize normal equation stuff

      nunk=3
      call nitil2(an,u,eltpel,nunk)
      nc=3
      ic(1)=1
      ic(2)=2
      ic(3)=3
      p=1.d0
      iobs=0
      sum =0.d0
      sum2=0.d0

*** loop over obs. equations

      rewind(lobs)
  100 read(lobs,end=777) iisat,trec1,trec2,r110,r210,r11,r21,rdd1
      if(iisat.le.0.or.iisat.gt.nsat) then
        continue

*** evaluate with latest coordinates --              sr,s geometric range
*** goal: scale direction cosines by w1 for l1 phase data

      else
*** debug stuff
        if(.not.getrs(iirfsv,trec1)) then
              write(lst,*) iisat,trec1,trec2
          stop 45397
        endif
        call comdtr(iirfsv,coord,orbits,trec1,trec2,r110,r210,dtr1,dtr2)
        call comrng(iirfsv,coord,orbits,trec1,trec2,r110,r210,dtr1,dtr2,
     *              sr1,sr2,cr)
        call comrng(iisat ,coord,orbits,trec1,trec2,r11,r21,dtr1,dtr2,
     *              s1,s2,c)

        sdd1=((s2-s1)-(sr2-sr1))
***     sdd1=((s2-s1)-(sr2-sr1))/w1
        do 10 i=1,3
   10   c(i)=(c(i)-cr(i))
***10   c(i)=(c(i)-cr(i))/w1

        cmo=sdd1-rdd1
***     write(lst,'(i3,f9.1,f14.3)') invprn(iisat),trec1,cmo
        if(dabs(cmo).le.rejlvl) then
          iobs=iobs+1
          sum =sum +cmo
          sum2=sum2+cmo*cmo
          call nrmal2(an,u,c,cmo,p,eltpel,ic,nc,nunk)
        endif
      endif
      go to 100

*** end of file encountered -- invert and update coordinates

  777 rewind(lobs)
      call fill(an,nunk)
      call invert(an,nunk)
      call ab(an,u,x,nunk,nunk,1)
      xmax=dmax1(dabs(x(1)),dabs(x(2)),dabs(x(3)))

*** compute rms from "basement window"  (v'pv=l'pl-x'x)  (ok for p=1.d0)

      call innerp(x,nunk,val)
      rms=dsqrt((eltpel-val)/dble(iobs))
      write(lst,'(a,f14.6,a,i7,a,f9.4,a)') 
     *    'max shift = ',xmax,' m.   n =',iobs,'    rms =',rms,' m.'

*** update coordinates for point no 2

      coord( 7)=coord( 7)-x(1)
      coord( 8)=coord( 8)-x(2)
      coord( 9)=coord( 9)-x(3)

      coord(10)=coord(10)-x(1)
      coord(11)=coord(11)-x(2)
      coord(12)=coord(12)-x(3)

      return
      end
************************************************************************
*** form observation equation units ************************************
************************************************************************
      subroutine formob(orbits,coord,nsat,tstrt,tstop,delt,
     *                  icols1,ndat1,icols2,ndat2,rms)

*** form observation equations for the first time (apply corrn.)

      implicit double precision(a-h,o-z)
      character*80 card
      character*4  scrnam
      external invprn,scrnam
      logical oget,dget,getrs,xyzgeo
      parameter(dnull=1.d31)
      parameter(maxprn=32)
      dimension ibuf1(maxprn)  ,ibuf2(maxprn)
      dimension dbuf1(maxprn*4),dbuf2(maxprn*4)
      dimension cs2(3),cs20(3),csdd(3),dummy(11)
      dimension orbits(*),coord(20),icols1(0:4),icols2(0:4)
      dimension ic(3),an(3,3),u(3),x(3)

      common/stugg/pi,rad,sol,f1,f2,w1,w2,we
      common/lus/lst,lpos,lrin1,lrin2,leph,lobs

      logical l2lin1,l2lin2,li1,li2
      save  /flags/
      common/flags/l2lin1,l2lin2

      if(.not.xyzgeo(coord( 4),coord( 5),coord( 6),gla1,glo1,eht1))
     *        stop 25551
      if(.not.xyzgeo(coord(10),coord(11),coord(12),gla2,glo2,eht2))
     *        stop 25552

*** use scatter/gather technique to order data into time series
*** open a file for each satellite, use logical units (11 -> 10+nsat)

      do 10 lio=11,10+nsat
   10 open(lio,file=scrnam(lio-10),form='unformatted',status='unknown')
      write(lst,*) 'opened', nsat,' scatter/gather files'

*** scan past headers of the rinex obs files

      rewind(lrin1)
   20 read(lrin1,'(a80)') card
      if(card(61:73).ne.'END OF HEADER') go to 20

      rewind(lrin2)
   30 read(lrin2,'(a80)') card
      if(card(61:73).ne.'END OF HEADER') go to 30

*** begin scatter of data to satellite-based files
*** repeat/until loop over rinex obs files
*** synchronization (+/- 0.2 sec) obtained using largest obs. interval

      write(*,*) ' begin scatter & formobs'

      tiny=0.2d0
      iflag1=0
      iflag2=0
      tmas=tstrt
      lu1=lrin1
      lu2=lrin2
      li1=l2lin1
      li2=l2lin2

   90 continue
      if(oget(lu1,li1,nsat,icols1,ndat1,tmas,trec1,ibuf1,dbuf1,iflag1)
     *   .and.
     *   oget(lu2,li2,nsat,icols2,ndat2,tmas,trec2,ibuf2,dbuf2,iflag2))
     *   then

*** only use epoch if reference satellite available

*** identification order {l1, p1/c1, l2, p2}  (2nd index)
***   p11  1st point, L1 phase
***   r11  1st point, L1 range
***   p12  1st point, L2 phase
***   p21  2nd point, L1 phase
***   r21  2nd point, L1 range
***   p22  2nd point, L2 phase

        if(getrs(iirfsv,tmas)) then
*** reference sat
          if(dget(iirfsv,ibuf1,ibuf2,dbuf1,dbuf2,nsat,
     *            p110,r110,p120,p210,r210,p220)) then
            call comdtr(iirfsv,coord,orbits,trec1,trec2,r110,r210,
     *                  dtr1,dtr2)
            call comrg0(iirfsv,coord,orbits,trec1,trec2,r110,r210,
     *                  dtr1,dtr2,s10,s20,cs20,gla1,glo1,eht1,
     *                  gla2,glo2,eht2,tdly10,tdly20)
*** correct observed range (measured too long due to tropo)
            r110=r110-tdly10
            r210=r210-tdly20
*** single difference (between stations) observed range
            rsd10=r210-r110
            ssd0=s20-s10

*** if valid data in both buffers, scatter based on internal sat id
***  bypass the reference satellite
*** goal: scale direction cosines by w1 for l1 phase data

            do 91 iisat=1,nsat
              if(iisat.ne.iirfsv) then
                if(dget(iisat,ibuf1,ibuf2,dbuf1,dbuf2,nsat,
     *                  p11,r11,p12,p21,r21,p22)) then
                  call comrg0(iisat,coord,orbits,trec1,trec2,r11,r21,
     *                        dtr1,dtr2,s1,s2,cs2,gla1,glo1,eht1,
     *                        gla2,glo2,eht2,tdly1,tdly2)
*** correct observed range (measured too long due to tropo)
                  r11 =r11 -tdly1
                  r21 =r21 -tdly2
*** double difference (now, between sats) observed range
                  rdd1=(r21-r11)-rsd10
                  sdd1=((s2-s1)-ssd0)

                  do 92 i=1,3
   92             csdd(i)=(cs2(i)-cs20(i))
                  write(10+iisat) trec1,trec2,r110,r210,r11,r21,
     *                            rdd1,sdd1,csdd
                endif
              endif
   91       continue
          endif
        endif

*** iflag=2 -> end of file,  iflag=1 -> data is future

      elseif(iflag1.eq.2.or.iflag2.eq.2) then
        write(lst,*) 'eof =',iflag1,iflag2
        go to 97
      endif
      tmas=tmas+delt

***
*** repeat/until loop until master time passes stop time
***

      if((tmas-tstop).le.tiny) go to 90

*** end of file raised -- proceed -- rewind scatter files

   97 continue
      do 98 lio=11,10+nsat
      end file(lio)
   98 rewind  (lio)

*** scatter completed -- gather data into observation equations
*** write header, data for a given satellite, and trailer
*** code should be robust if a single reference satellite

*** initialize normal equation stuff

      nunk=3
      call nitil2(an,u,eltpel,nunk)
      nc=3
      ic(1)=1
      ic(2)=2
      ic(3)=3
      p=1.d0
      iobs=0
      sum =0.d0
      sum2=0.d0

*** now write obs. equations and accumulate normal equations

      write(*,*) ' begin gather'
      do 100 iisat=1,nsat
        ihead=0
        write(lobs) ihead,dummy

  101   read (10+iisat,end=102) trec1,trec2,r110,r210,r11,r21,
     *                          rdd1,sdd1,csdd
        write(lobs    )   iisat,trec1,trec2,r110,r210,r11,r21,
     *                    rdd1
        cmo=sdd1-rdd1
        iobs=iobs+1
        sum =sum +cmo
        sum2=sum2+cmo*cmo
        call nrmal2(an,u,csdd,cmo,p,eltpel,ic,nc,nunk)
        go to 101

  102   continue
        ihead=99
        write(lobs) ihead,dummy
  100 continue

*** gather completed -- close the scatter/gather files

      do 200 lio=11,10+nsat
  200 close(lio,status='delete')

*** invert normal equations

      if(iobs.le.0) then
        write(lst,*) 'no observation equations -- fatal termination'
        write( * ,*) 'no observation equations -- fatal termination'
        stop 87854
      endif
      call fill(an,nunk)    
      call invert(an,nunk)     
      call ab(an,u,x,nunk,nunk,1)
      xmax=dmax1(dabs(x(1)),dabs(x(2)),dabs(x(3)))

*** compute rms from "basement window"  (v'pv=l'pl-x'x)  (ok for p=1.d0)

      call innerp(x,nunk,val)
      rms=dsqrt((eltpel-val)/dble(iobs))
      write(lst,'(a,f14.6,a,i7,a,f9.4,a)') 
     *    'max shift = ',xmax,' m.   n =',iobs,'    rms =',rms,' m.'

*** update coordinates for point no 2

      coord( 7)=coord( 7)-x(1)
      coord( 8)=coord( 8)-x(2)
      coord( 9)=coord( 9)-x(3)

      coord(10)=coord(10)-x(1)
      coord(11)=coord(11)-x(2)
      coord(12)=coord(12)-x(3)

      return
      end
      character*4 function scrnam(i)

*** name for scratch file  (MS fortran powerstation has dumb template)

      implicit double precision(a-h,o-z)
      character*4 scrf

      if(i.le.0.or.i.ge.36) stop 98786

*** A to Z, then 0 to 9

      if(i.le.26) then
        scrf = 'ZZ-'//char(i+64)
      else
        scrf = 'ZZ-'//char((i-27)+48)
      endif
      scrnam=scrf

      return
      end
************************************************************************
*** forward model stuff   **********************************************
************************************************************************
      subroutine comdtr(iiprn,coord,orbits,t1,t2,r1,r2,dtr1,dtr2)

*** compute receiver clock offsets at points 1 and 2
*** input times -- biased rcvr. clocks at receipt time

      implicit double precision(a-h,o-z)
      logical oterp
      external invprn
      dimension orbits(*),coord(20)
      common/lus/lst,lpos,lrin1,lrin2,leph,lobs
      common/stugg/pi,rad,sol,f1,f2,w1,w2,we

*** derive clock offsets from pseudorange/slant range discrepencies
*** r1, r2  -- input pseudoranges (also signal propagation delay, spd)

*** storage scheme for coord()
***
***  1  2  3       X Y Z  point 1 marker
***  4  5  6       X Y Z  point 1 antenna, L1
***  7  8  9       X Y Z  point 2 marker       (may be approx)
*** 10 11 12       X Y Z  point 2 antenna, L1  (may be approx)

*** 13 14 15 16    dl1,dl2,dn,de  point 1 offsets
*** 17 18 19 20    dl1,dl2,dn,de  point 2 offsets

*** ranges with respect to l1 phase center
*** translate biased receiver time to gps time at transmit
*** signal prop. delay since we interpolate orbit at transmit time

      dtr1=0.d0
      s1  =r1
      do 100 i=1,2
        spd=s1/sol
        t=(t1-spd)-dtr1
        if(.not.oterp(invprn(iiprn),t,orbits,x,y,z,dts,ierr)) then
          write(lst,*) 'oterp error in comdtr--',ierr,invprn(iiprn)
          stop 34553
        endif

*** maybe add sagnac one day -- need a good dt
***     sag=we*dt
***     dx=x+sag*y-coord( 4)
***     dy=y-sag*x-coord( 5)
        dx=x-coord( 4)
        dy=y-coord( 5)
        dz=z-coord( 6)
        s1=dsqrt(dx*dx+dy*dy+dz*dz)
        dtr1=(r1-s1)/sol + dts
  100 continue

*** second receiver

      dtr2=0.d0
      s2  =r2
      do 200 i=1,2
        spd=s2/sol
        t=(t2-spd)-dtr2
        if(.not.oterp(invprn(iiprn),t,orbits,x,y,z,dts,ierr)) then
          write(lst,*) 'oterp error in comdtr--',ierr,invprn(iiprn)
          stop 34554
        endif

*** maybe add sagnac one day
        dx=x-coord(10)
        dy=y-coord(11)
        dz=z-coord(12)
        s2=dsqrt(dx*dx+dy*dy+dz*dz)
        dtr2=(r2-s2)/sol + dts
  200 continue

      return
      end
      subroutine comrng(iiprn,coord,orbits,t1,t2,r1,r2,dtr1,dtr2,
     *                  s1,s2,cs2)

*** compute ranges to satellite from points 1 and 2  (s1,s2)
*** ITERATION VERSION
*** also compute direction cosines for point 2 to sat  (pt. 1 fixed)
*** input times -- biased rcvr. clocks at receipt time

      implicit double precision(a-h,o-z)
      logical oterp
      external invprn
      dimension orbits(*),coord(20),cs2(3)
      common/lus/lst,lpos,lrin1,lrin2,leph,lobs
      common/stugg/pi,rad,sol,f1,f2,w1,w2,we

*** r1, r2  -- input pseudoranges for signal propagation delay, spd
*** dtr1,2  -- input receiver time delays  (gps = rcvr - dtr)
*** note to me: can initialize loop using pseudoranges !!
      dummy=r1
      dummy=r2

*** storage scheme for coord()
***
***  1  2  3       X Y Z  point 1 marker
***  4  5  6       X Y Z  point 1 antenna, L1
***  7  8  9       X Y Z  point 2 marker       (may be approx)
*** 10 11 12       X Y Z  point 2 antenna, L1  (may be approx)

*** 13 14 15 16    dl1,dl2,dn,de  point 1 offsets
*** 17 18 19 20    dl1,dl2,dn,de  point 2 offsets

*** ranges with respect to l1 phase center
*** direction cosines computed using geometric (slant) range

*** iterate the light loop a fixed number of times

        dt=0.075d0
        do 1 i=1,3

*** get satellite position -- correct for receiver bias

          tt=(t1-dtr1)-dt
          if(.not.oterp(invprn(iiprn),tt,orbits,x,y,z,dts,ierr)) then
            write(lst,*) 'oterp error in comrng--',ierr,invprn(iiprn)
            stop 34551
          endif
          dx=x-coord( 4)
          dy=y-coord( 5)
          dz=z-coord( 6)
          s1=dsqrt(dx*dx+dy*dy+dz*dz)
          dt=s1/sol
    1   continue

*** now apply sagnac effect for the range

        sag=we*dt
        dx=x+sag*y-coord( 4)
        dy=y-sag*x-coord( 5)
        dz=z      -coord( 6)
        s1=dsqrt(dx*dx+dy*dy+dz*dz)

*** now do point 2  -- including direction cosines
*** iterate the light loop a fixed number of times

        dt=0.075d0
        do 2 i=1,3

*** get satellite position -- correct for receiver bias

          tt=(t2-dtr2)-dt
          if(.not.oterp(invprn(iiprn),tt,orbits,x,y,z,dts,ierr)) then
            write(lst,*) 'oterp error in comrng--',ierr,invprn(iiprn)
            stop 34552
          endif
          dx=x-coord(10)
          dy=y-coord(11)
          dz=z-coord(12)
          s2=dsqrt(dx*dx+dy*dy+dz*dz)
          dt=s2/sol
    2   continue

*** now apply sagnac effect for the range

        sag=we*dt
        dx=x+sag*y-coord(10)
        dy=y-sag*x-coord(11)
        dz=z      -coord(12)
        s2=dsqrt(dx*dx+dy*dy+dz*dz)

*** direction cosines for station 2   (only)

      cs2(1)=-dx/s2
      cs2(2)=-dy/s2
      cs2(3)=-dz/s2

      return
      end
      subroutine comrg0(iiprn,coord,orbits,t1,t2,r1,r2,dtr1,dtr2,
     *                  s1,s2,cs2,gla1,glo1,eht1,gla2,glo2,eht2,
     *                  tdly1,tdly2)

*** compute ranges to satellite from points 1 and 2  (s1,s2)
*** this version computes tropo corrections, tdly1/2
*** also compute direction cosines for point 2 to sat  (pt. 1 fixed)
*** input times -- biased rcvr. clocks at receipt time

      implicit double precision(a-h,o-z)
      logical oterp
      external invprn
      dimension orbits(*),coord(20),cs2(3)
      common/lus/lst,lpos,lrin1,lrin2,leph,lobs
      common/stugg/pi,rad,sol,f1,f2,w1,w2,we

*** r1, r2  -- input pseudoranges for signal propagation delay, spd
*** dtr1,2  -- input receiver time delays  (gps = rcvr - dtr)
*** note to me: can initialize loop using pseudoranges !!
      dummy=r1
      dummy=r2

*** storage scheme for coord()
***
***  1  2  3       X Y Z  point 1 marker
***  4  5  6       X Y Z  point 1 antenna, L1
***  7  8  9       X Y Z  point 2 marker       (may be approx)
*** 10 11 12       X Y Z  point 2 antenna, L1  (may be approx)

*** 13 14 15 16    dl1,dl2,dn,de  point 1 offsets
*** 17 18 19 20    dl1,dl2,dn,de  point 2 offsets

*** ranges with respect to l1 phase center
*** direction cosines computed using geometric (slant) range

*** iterate the light loop a fixed number of times

        dt=0.075d0
        do 1 i=1,3

*** get satellite position -- correct for receiver bias

          tt=(t1-dtr1)-dt
          if(.not.oterp(invprn(iiprn),tt,orbits,x,y,z,dts,ierr)) then
            write(lst,*) 'oterp error in comrng--',ierr,invprn(iiprn)
            stop 34551
          endif
          dx=x-coord( 4)
          dy=y-coord( 5)
          dz=z-coord( 6)
          s1=dsqrt(dx*dx+dy*dy+dz*dz)
          dt=s1/sol
    1   continue

*** now apply sagnac effect for the range

        sag=we*dt
        dx=x+sag*y-coord( 4)
        dy=y-sag*x-coord( 5)
        dz=z      -coord( 6)
        s1=dsqrt(dx*dx+dy*dy+dz*dz)

*** compute vertical angle from point 1 antenna center to satellite
*** debug -- something wrong here

      call rge(gla1,glo1,dn,de,du,dx,dy,dz)
      ds=dsqrt(dn*dn+de*de)
      va1=datan(du/ds)

*** compute tropo delay -- treat ortho ht. as ellipsoid ht.

      oht=eht1
***  call gettdl(gla1,oht,va1,tdly1)
*** debug -- something wrong here
      tdly1=0.d0

*** now do point 2  -- including direction cosines
*** iterate the light loop a fixed number of times

        dt=0.075d0
        do 2 i=1,3

*** get satellite position -- correct for receiver bias

          tt=(t2-dtr2)-dt
          if(.not.oterp(invprn(iiprn),tt,orbits,x,y,z,dts,ierr)) then
            write(lst,*) 'oterp error in comrng--',ierr,invprn(iiprn)
            stop 34552
          endif
          dx=x-coord(10)
          dy=y-coord(11)
          dz=z-coord(12)
          s2=dsqrt(dx*dx+dy*dy+dz*dz)
          dt=s2/sol
    2   continue

*** now apply sagnac effect for the range

        sag=we*dt
        dx=x+sag*y-coord(10)
        dy=y-sag*x-coord(11)
        dz=z      -coord(12)
        s2=dsqrt(dx*dx+dy*dy+dz*dz)

*** direction cosines for station 2   (only)

      cs2(1)=-dx/s2
      cs2(2)=-dy/s2
      cs2(3)=-dz/s2

*** compute vertical angle from point 2 antenna center to satellite
*** debug -- something wrong here

      call rge(gla2,glo2,dn,de,du,dx,dy,dz)
      ds=dsqrt(dn*dn+de*de)
      va2=datan(du/ds)

*** compute tropo delay -- treat ortho ht. as ellipsoid ht.

      oht=eht2
***   call gettdl(gla2,oht,va2,tdly2)
*** debug -- something wrong here
      tdly2=0.d0

      return
      end
      subroutine comrg2(iiprn,coord,orbits,t1,t2,r1,r2,dtr1,dtr2,
     *                  s1,s2,gla,glo,vad)

*** compute ranges to satellite from points 1 and 2  (s1,s2)
*** ITERATION VERSION
*** also compute direction cosines for point 2 to sat  (pt. 1 fixed)
*** input times -- biased rcvr. clocks at receipt time

      implicit double precision(a-h,o-z)
      logical oterp
      external invprn
      dimension orbits(*),coord(20)
      common/lus/lst,lpos,lrin1,lrin2,leph,lobs
      common/stugg/pi,rad,sol,f1,f2,w1,w2,we

*** r1, r2  -- input pseudoranges for signal propagation delay, spd
*** dtr1,2  -- input receiver time delays  (gps = rcvr - dtr)
*** note to me: can initialize loop using pseudoranges !!
      dummy=r1
      dummy=r2

*** storage scheme for coord()
***
***  1  2  3       X Y Z  point 1 marker
***  4  5  6       X Y Z  point 1 antenna, L1
***  7  8  9       X Y Z  point 2 marker       (may be approx)
*** 10 11 12       X Y Z  point 2 antenna, L1  (may be approx)

*** 13 14 15 16    dl1,dl2,dn,de  point 1 offsets
*** 17 18 19 20    dl1,dl2,dn,de  point 2 offsets

*** ranges with respect to l1 phase center
*** direction cosines computed using geometric (slant) range

*** iterate the light loop a fixed number of times

        dt=0.075d0
        do 1 i=1,3

*** get satellite position -- correct for receiver bias

          tt=(t1-dtr1)-dt
          if(.not.oterp(invprn(iiprn),tt,orbits,x,y,z,dts,ierr)) then
            write(lst,*) 'oterp error in comrng--',ierr,invprn(iiprn)
            stop 34551
          endif
          dx=x-coord( 4)
          dy=y-coord( 5)
          dz=z-coord( 6)
          s1=dsqrt(dx*dx+dy*dy+dz*dz)
          dt=s1/sol
    1   continue

*** now apply sagnac effect for the range

        sag=we*dt
        dx=x+sag*y-coord( 4)
        dy=y-sag*x-coord( 5)
        dz=z      -coord( 6)
        s1=dsqrt(dx*dx+dy*dy+dz*dz)

*** now do point 2  -- including direction cosines
*** iterate the light loop a fixed number of times

        dt=0.075d0
        do 2 i=1,3

*** get satellite position -- correct for receiver bias

          tt=(t2-dtr2)-dt
          if(.not.oterp(invprn(iiprn),tt,orbits,x,y,z,dts,ierr)) then
            write(lst,*) 'oterp error in comrng--',ierr,invprn(iiprn)
            stop 34552
          endif
          dx=x-coord(10)
          dy=y-coord(11)
          dz=z-coord(12)
          s2=dsqrt(dx*dx+dy*dy+dz*dz)
          dt=s2/sol
    2   continue

*** now apply sagnac effect for the range

        sag=we*dt
        dx=x+sag*y-coord(10)
        dy=y-sag*x-coord(11)
        dz=z      -coord(12)
        s2=dsqrt(dx*dx+dy*dy+dz*dz)

*** compute vertical angle from point 2 antenna center to (non-ref) satellite

      call rge(gla,glo,dn,de,du,dx,dy,dz)
      ds=dsqrt(dn*dn+de*de)
      va=datan(du/ds)
      vad=va*rad

      return
      end
************************************************************************
*** utility work routines **********************************************
************************************************************************
      logical function oget(lrin,l2line,nsat,icols,ndat,tmas,trec,ibuf,
     *                      dbuf,iflag)

*** load epochs of data from a rinex obs file until match master time

*** iflag=0   successful buffer load, synchronous with tmas
*** iflag=1   not yet loaded, non-synchronous, (future time)
*** iflag=2   unsuccessful, end of file encountered

      implicit double precision(a-h,o-z)
      character*80 card,card2
      logical getprn,l2line
      parameter(dnull=1.d31)
      parameter(mxsat=12)
      dimension iprns(mxsat),buf(5),icols(0:4)
      dimension ibuf(nsat),dbuf(nsat,4)
      common/lus/lst,lpos,lrin1,lrin2,leph,lobs

***********************************************************************
*** contols how much synchroniztion is needed  ************************
***********************************************************************

****  tiny=5.d-2
      tiny=3.d-1

*** if buffers already loaded, and time is in future of master time,
*** its because of an earlier (multi-epoch) dropout
*** action: retain iflag=1 and kick out (wait for master to catch up)

      if(iflag.eq.1.and.(trec-tmas).ge.tiny) then
        oget=.false.
        return
      endif

*** if time matches on input, its because of prior invocation
***    (an iflag=1 situation)  buffers already loaded

      if(dabs(trec-tmas).le.tiny) then
        if(iflag.ne.1) stop 34591
        iflag=0
        oget=.true.
        return
      endif

*** read data until time match or future data bundle loaded
*** note: code assumes no more than 12 satellites per epoch
*** iprns() is the buffer to hold one record of prn numbers
*** assume no more than 5 columns of data per record, in buf()
*** lli (loss of lock indicator) not used
*** iss (signal strength value ) not used

  100 read(lrin,1,end=777) iyr,imo,idy,ihr,imn,sec,
     *                     iepflg,nprn,iprns,roff
    1 format(5i3,f11.7,i3,i3,12(1x,i2),f12.9)
      if(iepflg.gt.1) then
        do 108 i=1,nprn
  108   read(lrin,'(a80)') card
        go to 100
      endif
      if(iyr.lt.80) then
        iyr=iyr+2000
      elseif(iyr.le.99) then
        iyr=iyr+1900
      endif
      if(nprn.le.0.or.nprn.gt.mxsat) stop 24303
      call civjts(iyr,imo,idy,ihr,imn,sec,trec)

*** header too old -- bypass data and read another epoch header

      if(trec.lt.tmas-tiny) then
        do 110 i=1,nprn
          read(lrin,'(a80)') card
          if(l2line) read(lrin,'(a80)') card2
  110   continue

*** data in synch or header in the future -- process data records

      elseif(dabs(trec-tmas).le.tiny.or.
     *          trec.gt.tmas+tiny       ) then
      if(trec.gt.tmas+tiny) then
      endif

*** initialize buffer flags and load the buffer
*** check empty and 0.0 fields and store null (dnull) instead
*** warning: code assumes data in first card only !!!!!!!!!!!!!!!!!!!!!!!!

        ndatx=ndat
        if(ndatx.gt.5) ndatx=5

        do 120 i=1,nsat
  120   ibuf(i)=0

        do 130 i=1,nprn
          if(iprns(i).gt.0) then
            read(lrin,'(a80)') card
            if(l2line) read(lrin,'(a80)') card2
            read(card,2) (buf(ic),lli,iss,ic=1,ndatx)
    2       format(5(f14.3,i1,i1))
            do 131 ic=1,ndatx
              if(dabs(buf(ic)).le.0.001d0) buf(ic)=dnull
  131       continue

*** if prn lookup fails -- assume due to missing satellite orbit

*** debug -- review this !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            if(.not.getprn(iprns(i),iisat)) go to 130

*** anal-retentive assertion check -- delete in beta version !!!!!!!!!
***         if(iisat.le.0.or.iisat.gt.nsat) then
***           write(lst,*) 'prn index range check fail in oget'
***           stop 32645
***         endif

*** transfer to buffer -- order l1, p1/c1, l2, p2  (j index)
***                    -- insert null if data not gathered

            ibuf(iisat)  =1
            do 132 j=1,4
              if(icols(j).gt.0) then
                dbuf(iisat,j)=buf(icols(j))

*** exception logic for p1/c1

                if(j.eq.2.and.dbuf(iisat,2).ge.dnull) then
                  if(icols(0).gt.0) then
                    dbuf(iisat,2)=buf(icols(0))
                  else
                    dbuf(iisat,2)=dnull
                  endif
                endif

*** no p1/c1 at all

              else
                dbuf(iisat,j)=dnull
              endif
  132       continue

          else
            write(lst,*) 'zero prn number was not expected'
            stop 32646
          endif
  130   continue

*** was that header in synch?

        if(dabs(trec-tmas).le.tiny) then
          iflag=0
          oget=.true.
          return

*** no, the header was in the future  (due to if test in this clause)

        else
          iflag=1
          oget=.false.
          return
        endif

*** should never get here --- some kinda logic bug

      else
        write(lst,*) 'time synch logic bug in oget()'
        stop 36649
      endif

*** read another epoch header

      go to 100

*** end of file encountered

  777 iflag=2
      oget=.false.

      return
      end
      logical function dget(iisat,ibuf1,ibuf2,dbuf1,dbuf2,nsat,
     *                      p11,r11,p12,p21,r21,p22)

*** see if valid phase data in both buffers for iisat-th satellite

      implicit double precision(a-h,o-z)
      parameter(dnull=1.d31)
      dimension ibuf1(nsat),dbuf1(nsat,4)
      dimension ibuf2(nsat),dbuf2(nsat,4)
      common/lus/lst,lpos,lrin1,lrin2,leph,lobs

*** also return l1 pseudoranges for transmission time

      if(iisat.le.0.or.iisat.gt.nsat)  stop 87473

*** assume the worst

      dget=.false.

      if(ibuf1(iisat).eq.0) return
      if(ibuf2(iisat).eq.0) return

*** order {l1, p1/c1, l2, p2}  (2nd index)
***   p11  1st point, l1 phase
***   r11  1st point, l1 range
***   p12  1st point, l2 phase
***   p21  2nd point, l1 phase
***   r21  2nd point, l1 range
***   p22  2nd point, l2 phase

      p11=dbuf1(iisat,1)
      r11=dbuf1(iisat,2)
      p12=dbuf1(iisat,3)
      p21=dbuf2(iisat,1)
      r21=dbuf2(iisat,2)
      p22=dbuf2(iisat,3)

      if(r11.ge.dnull) return
      if(r21.ge.dnull) return

*** things are all right after all

      dget=.true.

      return
      end
************************************************************************
*** initial scan routines **********************************************
************************************************************************
      subroutine recon(nsat,nepoch,tstrt,tstop,start,stop,delt,
     *                 icols1,ndat1,aid1,icols2,ndat2,aid2,coord)

*** preliminary scan of files for critical information

      implicit double precision(a-h,o-z)
      character*14 osource
      character*4 aid1,aid2
      dimension icols1(0:4),icols2(0:4),coord(20)
      common/lus/lst,lpos,lrin1,lrin2,leph,lobs

      logical l2lin1,l2lin2
      save /flags/
      common/flags/l2lin1,l2lin2

*** initialize coordinates and offsets

      do 3 i=1,20
    3 coord(i)=0.d0

*** scan orbit for start/stop times and epoch interval (sp3 format)

      read(leph,1) iyr,imo,idy,ihr,imn,sec,nepoch,osource
    1 format(3x,i4,4i3,f12.8,i8,7x,a14)
      read(leph,2) deltat
    2 format(23x,f15.8)
      write(lst,*) 'orbit ref sys/source -- ',osource

*** set initial epoch for all subsequent time computations

      call setjd0(iyr,imo,idy)
      call civjts(iyr,imo,idy,ihr,imn,sec,ostrt)
      ostop=ostrt+deltat*(nepoch-1)
      write(lst,'(a,2f10.2,f7.2,i6)') 'orbit strt/stop,delt,neph ',
     *                                 ostrt,ostop,deltat,nepoch
***** call jtsciv(ostop,iyr,imo,idy,ihr,imn,sec)

*** scan obs files, get start/stop times, obs interval, sat. list, pos.

      call newprn
      lrin=lrin1
      call oscan(lrin,aid1,tstrt1,tstop1,delt1,icols1,ndat1,x1,y1,z1,
     *           l2lin1)
      write(lst,'(a,3f10.2)') 'file1  start/stop, intrvl',
     *                                tstrt1,tstop1,delt1
      coord(1)=x1
      coord(2)=y1
      coord(3)=z1
      rewind(lrin1)

      lrin=lrin2
      call oscan(lrin,aid2,tstrt2,tstop2,delt2,icols2,ndat2,x2,y2,z2,
     *           l2lin2)
      write(lst,'(a,3f10.2)') 'file2  start/stop, intrvl',
     *                                tstrt2,tstop2,delt2
      coord(7)=x2
      coord(8)=y2
      coord(9)=z2
      rewind(lrin2)

*** find common data interval and overlap with orbit interval
*** find common obs. interval (max of all obs. intervals)
*** insure data intervals are compatible with the obs interval

      start=dmax1(tstrt1,tstrt2)
      stop =dmin1(tstop1,tstop2)
      delt =dmax1(delt1,delt2)

      delt=dble(idnint(delt))
      start=delt*idnint(start/delt)
      stop =delt*idnint(stop /delt)

      write(lst,'(a,3f10.0)') 'common start/stop, intrvl',
     *                         start,stop,delt
      if(stop.le.start) then
        write(lst,*) 'receiver data do not overlap -- fatal'
        stop 43587
      endif

      tiny=0.2d0
      if(start.lt.ostrt-tiny.or.stop.gt.ostop+tiny) then
        write(lst,*) 'orbit interval too short -- fatal'
        write(lst,'(a,2f10.0)') 'orbit start/stop ',ostrt,ostop
        stop 43588
      endif

*** identify interval and epoch number for ephemeris load
*** the ephemeris interval is reset based on received data

      tstrt=(idnint(start/deltat)-9)*deltat
      tstop=(idnint(stop /deltat)+9)*deltat

      if(tstrt.lt.ostrt) tstrt=ostrt
      if(tstop.gt.ostop) tstop=ostop
*** force these times to orbit interval granularity
      tstrt=(idnint((tstrt-ostrt)/deltat))*deltat+ostrt
      tstop=(idnint((tstop-ostrt)/deltat))*deltat+ostrt
      write(lst,'(a,2f10.2)') 'lodorb start/stop',tstrt,tstop
      nepoch=idnint((tstop-tstrt)/deltat)+1

      call numprn(nsat)
      write(lst,*) 'satellites and epochs',nsat,nepoch

      return
      end
      subroutine oscan(lrin,aid,tstrt,tstop,delt,icols,ndat,x,y,z,
     *                 l2line)

*** scan rinex obs file, accumulate start/stop times, obs. interval,
***    and satellite list

*** 4-mar-96: force multiple of 5 seconds if interval greater than 5
      implicit double precision(a-h,o-z)
      logical putprn,l2line
      character*80 card
      character*4 aid
      character*2 as(9)
      dimension icols(0:4)
      parameter(mxsat=12)
      dimension iprns(mxsat)
      common/lus/lst,lpos,lrin1,lrin2,leph,lobs

*** initialize approximate marker position

      x=0.d0
      y=0.d0
      z=0.d0

*** scan header for observation type description record

   10 read(lrin,1) card
    1 format(a80)
      
*** extract marker name

      if(card(61:71).eq.'MARKER NAME') then
        aid=card(1:4)
        write(lst,*) 'station id -- ',aid

*** get approx xyz

      elseif(card(61:79).eq.'APPROX POSITION XYZ') then
        read(card,'(3f14.4)') x,y,z

*** get antenna dh,de,dn  (this is useless due to definition problems

      elseif(card(61:80).eq.'ANTENNA: DELTA H/E/N') then
        go to 10

*** determine which data are in which columns
*** look in order l1, p1/c1, l2, p2
*** (design restriction, do not allow more than 5 types)

      elseif(card(61:79).eq.'# / TYPES OF OBSERV') then
        read(card,'(i6)') ndat
        if(ndat.le.0) stop 34989
        if(ndat.gt.5) then
          l2line=.true.
*****     write(lst,*) '2 line rinex format detected'
        else
          l2line=.false.
*****     write(lst,*) '1 line rinex format detected'
        endif
        read(card,3) (as(i),i=1,ndat)
    3   format(6x,9(4x,a2))

        call alook('L1',as,ndat,i)
        icols(1)=i

*** dual logic, some sats show P1, others show C1

        call alook('P1',as,ndat,i)
        if(i.eq.0) then
          call alook('C1',as,ndat,i)
          icols(2)=i
          icols(0)=0

*** found P1 (still check for C1)

        else
          icols(2)=i
          call alook('C1',as,ndat,i)
          icols(0)=i
        endif

        call alook('L2',as,ndat,i)
        icols(3)=i

        call alook('P2',as,ndat,i)
        icols(4)=i

*** end of header

      elseif(card(61:73).eq.'END OF HEADER') then
        go to 900
      endif
      go to 10

*** end of header encountered -- proceed with processing

  900 continue

*** read the first epoch (and get start time)
*** note: code assumes no more than mxsat(12) satellites per epoch
*** iprns() is the buffer to hold one record of prn numbers
 
  100 read(lrin,2) iyr,imo,idy,ihr,imn,sec,iepflg,nprn,iprns,roff
    2 format(5i3,f11.7,i3,i3,12(1x,i2),f12.9)
      if(iepflg.gt.1) then
        do 108 i=1,nprn
  108   read(lrin,1) card
        go to 100
      endif
      if(iyr.lt.80) then
        iyr=iyr+2000
      elseif(iyr.le.99) then
        iyr=iyr+1900
      endif
      if(nprn.le.0.or.nprn.gt.mxsat) stop 24307
      call civjts(iyr,imo,idy,ihr,imn,sec,tstrt)
      tsecx=tstrt

      do 104 i=1,nprn
        if(iprns(i).gt.0) then
          if(.not.putprn(iprns(i),idup)) then
            if(idup.eq.0) then
              write(lst,*) 'prn table overflow'
              stop 32898
            endif
          endif
        endif
  104 continue

*** bypass the data records for now

      do 105 i=1,nprn
        read(lrin,1) card
        if(l2line) read(lrin,1) card
  105 continue

*** read subsequent epoch, update the end time, accumulate delt's
*** note: code assumes no more than mxsat(12) satellites per epoch
*** iprns() is the buffer to hold one record of prn numbers

      n=0
      delt=0.d0
  101 read(lrin,2,end=777,err=666) iyr,imo,idy,ihr,imn,sec,
     *                             iepflg,nprn,iprns,roff
      if(iepflg.gt.1) then
        do 109 i=1,nprn
  109   read(lrin,1) card
        go to 101
      endif
      if(iyr.lt.80) then
        iyr=iyr+2000
      elseif(iyr.le.99) then
        iyr=iyr+1900
      endif
      if(nprn.le.0.or.nprn.gt.mxsat) stop 24308
      call civjts(iyr,imo,idy,ihr,imn,sec,tstop)
      n=n+1
      delt=delt+(tstop-tsecx)
      tsecx=tstop

      do 102 i=1,nprn
        if(iprns(i).gt.0) then
          if(.not.putprn(iprns(i),idup)) then
            if(idup.eq.0) then
              write(lst,*) 'prn table overflow'
              stop 32899
            endif
          endif
        endif
  102 continue

*** bypass the data records for now

      do 103 i=1,nprn
        read(lrin,1) card
        if(l2line) read(lrin,1) card
  103 continue
      go to 101

*** end of file encountered -- get average obs interval
*** note:  force multiple of 5 seconds if interval greater than 5

  777 delt=delt/dble(n)
      if(delt.gt.5.1d0) then
        idelt=idnint(delt/5.d0)
        delt=5.d0*idelt
      endif
      return

***
*** error handler **********************************************
***

  666 write(*,*) 'read error encountered -- see .lst'
      write(lst,*) 'last valid data header was:'
      write(lst,'(5i3,f9.3,a)') iyr,imo,idy,ihr,imn,sec,' was OK'
      stop 66663
      return
      end
      subroutine alook(a,as,ndat,ival)

*** search array of character codes for observation types

      character*2 a,as(ndat)

      do 1 i=1,ndat
        if(a.eq.as(i)) then
          ival=i
          return
        endif
    1 continue

*** fell thru loop -- code not present

      ival=0

      return
      end
************************************************************************
*** compute satellite reference table **********************************
************************************************************************
      subroutine makrst(dstrt,dstop,delt,coord,orbits,nsat)

*** compute satellite reference table
*** build table based on satellites with highest vertical angle
*** this version uses original data spacing (can be changed)
 
      implicit double precision(a-h,o-z)
      parameter(tiny=1.d-1)
      logical oterp,putrs,xyzgeo
      external invprn
      dimension coord(20)
      parameter(mxsat=12)
      dimension orbits(*)
      common/lus/lst,lpos,lrin1,lrin2,leph,lobs

*** initialize the reference satellite table

      call openrs

*** use l1 phase center for point 1  (rx -- receiver)
*** lat/lon used to rotate ant/sat. vector into l.h.s. for vert.ang.

      xrx=coord(4)
      yrx=coord(5)
      zrx=coord(6)
      if(.not.xyzgeo(xrx,yrx,zrx,gla,glo,eht)) stop 93278

*** master loop over all data epochs -- inner loop over prn table
*** use coarse increment for this table, target delt = targdt

      targdt=5.d0*60.d0
      ifac=idnint(targdt/delt)
      if(ifac.le.1) then
        delt2=delt
      else
        delt2=ifac*delt
      endif

      tsec=dstrt
  100 if(tsec.le.dstop+tiny) then
        vadmx=0.d0
        iiprn=0
        do 10 i=1,nsat
          if(.not.oterp(invprn(i),tsec,orbits,x,y,z,dt,ierr)) then
            write(lst,*) 'oterp error in makrst--',ierr,invprn(i),tsec
            stop 34958
          endif
          call comvad(gla,glo,xrx,yrx,zrx,x,y,z,vad)
          if(vad.gt.vadmx) then
            vadmx=vad
            iiprn=i
          endif
   10   continue
        if(iiprn.eq.0) stop 32875

*** store internal prn into table

        if(.not.putrs(iiprn,tsec,nsvrt)) then
          write(lst,*) 'ref sat table overflow at ',tsec
          stop 87435
        endif

        tsec=tsec+delt2
        go to 100
      endif

*** store last internal prn to force stop time (delt vs. delt2)

      if(.not.putrs(iiprn,dstop,nsvrt)) stop 87436
      write(lst,*) nsvrt,' entries in reference satellite table'
      call dumprs

*** condense the table  (not implemented)

      return
      end
      subroutine comvad(gla,glo,x1,y1,z1,x2,y2,z2,vad)

*** compute vertical angle from point 1 to point 2
*** input in meters (lat/lon in radians), output in degrees

      implicit double precision(a-h,o-z)
      common/stugg/pi,rad,sol,f1,f2,w1,w2,we

      dx=x2-x1
      dy=y2-y1
      dz=z2-z1
      call rge(gla,glo,dn,de,du,dx,dy,dz)

      ds=dsqrt(dn*dn+de*de)
      va=datan(du/ds)
      vad=va*rad

      return
      end
************************************************************************
*** troposphere handlers ***********************************************
************************************************************************
      subroutine inimet(lyesno)

*** initialize met stuff

      implicit double precision(a-h,o-z)
      logical lyesno

      save  /metstf/
      logical ltropo
      common/metstf/bpxx,rhxx,tcxx,ltropo

*** default met -- barom. press.  1013 mbar
***             -- rel. humidity  0.5
***             -- surf. temp(c)  20 deg. celsius

      bpxx=1013.d0
      rhxx=0.5d0
      tcxx=20.d0

*** do we compute met? (yes/no)

      ltropo=lyesno

      return
      end
      subroutine gettdl(gla,oht,v,tdelay)

*** compute tropo delay in meters

*** note: just use default met values

      implicit double precision(a-h,o-z)

      save  /metstf/
      logical ltropo
      common/metstf/bpxx,rhxx,tcxx,ltropo

      if(.not.ltropo) then
        tdelay=0.d0
        return
      endif

*** use default met values -- read rinex someday

      bp=bpxx
      rh=rhxx
      tc=tcxx

*** do the work and get tdelay in meters

      call tropo(bp,rh,tc,gla,oht,v,tdelay)

      return
      end
      subroutine tropo(bp,rh,ts,gla,oht,v,tdelay)

*** compute troposphere delay (saastamoinen with herring mapping)

      implicit double precision(a-h,o-z)
      double precision lhz,lwz,mh,mw

*** v      -- vertical angle (radians)
*** gla    -- geodetic latitude (radians)
*** oht    -- orthometric height of rcvr. (meters)
*** ts     -- surface temperature         (centigrade)

*** sanity checks

      if(bp .lt.   900.d0 .or. bp .gt. 1100.d0) stop 80001
      if(rh .lt.     0.d0 .or. rh .gt.    1.d0) stop 80002
      if(ts .lt.   -50.d0 .or. ts .gt.   50.d0) stop 80003
      if(gla.lt. -1.572d0 .or. gla.gt. 1.572d0) stop 80004
      if(oht.lt.  -150.d0 .or. oht.gt. 5000.d0) stop 80005
      if(v  .lt.   -0.5d0 .or. v  .gt. 1.572d0)then
         write(*,*) 'v=',v
         stop 80006
      endif

***
*** zenith delay section -- saastamoinen
***

      call saast(bp,rh,ts,gla,oht,lhz,lwz)

***
*** mapping function section
*** herring: 1992, refraction of transatmospheric signals in geodesy
*** note: mapping fcns derived: lat 27N-65N, and ht. 0-1.6 km
***

*** hs -- orthometric height of rcvr. (kilometers)

      hs=oht/1000.d0

      sine=dsin(v)
      cgla=dcos(gla)
      ts10=ts-10.d0

*** hydrostatic mapping fcn  (herring eq. 8)  (RMS ~ 3-10 mm)

      ah=(1.2330d0+0.0139d0*cgla-0.0209d0*hs+0.00215d0*ts10)*1.d-3
      bh=(3.1612d0-0.1600d0*cgla-0.0331d0*hs+0.00206d0*ts10)*1.d-3
      ch=(71.244d0-4.2930d0*cgla-0.1490d0*hs-0.00210d0*ts10)*1.d-3

      mh=(1.d0+ah/(1.d0+bh/(1.d0+ch)))/
     *   (sine+ah/(sine+bh/(sine+ch)))


*** wet mapping fcn  (herring eq. 9)

      aw=(0.583d0-0.011d0*cgla-0.052d0*hs+0.0014d0*ts10)*1.d-3
      bw=(1.402d0-0.102d0*cgla-0.101d0*hs+0.0020d0*ts10)*1.d-3
      cw=(45.85d0-1.910d0*cgla-1.290d0*hs+0.0150d0*ts10)*1.d-3

      mw=(1.d0+aw/(1.d0+bw/(1.d0+cw)))/
     *   (sine+aw/(sine+bw/(sine+cw)))

***
*** total atmospheric delay  (herring eq 3)
***

      tdelay=lhz*mh+lwz*mw

      return
      end
      subroutine saast(bp,rh,tc,gla,oht,dryzen,wetzen)

*** compute wet and dry zenith delay using the saastamoinen formula

************************************************************************
* input
* bp    barometric pressure [millibars]
* rh    relative humidity   [0.0 --> 1.0]
* tc    temperature         [centigrade]
* gla   site latitude       [radians]
* oht   orthometric height  [m]
*
* output
* dryzen  dry atmosphere zenith path delay [m]
* wetzen  wet zenith delay [m]
*
*  87-02-15, written by:   j. l. davis
*  91-01-09, mss, standard header.  remove vlbi specific code
*  96-04-11, dgm, combined wet and dry, and deleted derivative stuff
************************************************************************

      implicit double precision (a-h,o-z)

*** saastamoinen's function for gravity

      fgrav=1.d0-0.00266d0*dcos(2.d0*gla)-0.00028d-03*oht

*** water vapor partial pressure in mbars   (proportional to rh)

      tk=tc+273.15d0
      pp=rh*6.11d0*(tk/273.15d0)**(-5.3d0)*dexp(25.2d0*tc/tk)

*** zenith delays

      wetzen=0.0022768d0*(1255.d0/tk+0.05d0)*pp/fgrav
      dryzen=0.0022768d0*                    bp/fgrav

      return
      end
************************************************************************
*** marker and antenna position ****************************************
************************************************************************
      subroutine lodpos(coord,aid1,aid2)

*** get precise coordinates from the lpos file
*** precise coordinates may not be available for point2

      implicit double precision(a-h,o-z)
      character*4 aid1,aid2
      logical pscan
      dimension coord(20)
      common/lus/lst,lpos,lrin1,lrin2,leph,lobs
      save  /xyz2o/
      common/xyz2o/x2o,y2o,z2o

*** storage scheme for coord()
***
***  1  2  3       X Y Z  point 1 marker
***  4  5  6       X Y Z  point 1 antenna, L1
***  7  8  9       X Y Z  point 2 marker       (may be approx)
*** 10 11 12       X Y Z  point 2 antenna, L1  (may be approx)

*** 13 14 15 16    dl1,dl2,dn,de  point 1 offsets
*** 17 18 19 20    dl1,dl2,dn,de  point 2 offsets

      rewind(lpos)
      if(.not.pscan(aid1,x,y,z)) then
        write(lst,*) 'warning --',aid1,' -- id not found in pos file'
        write(lst,*) 'this is not good -- point 2 unreliable'
      else
        coord(1)=x
        coord(2)=y
        coord(3)=z
      endif

*** point2

      rewind(lpos)
      if(.not.pscan(aid2,x,y,z)) then
        write(lst,*) 'info    --',aid2,' -- id not found in pos file'
      else
        coord(7)=x
        coord(8)=y
        coord(9)=z
      endif

*** point2 -- save marker in common /xyz2o/ for pos. shifts

      x2o=coord(7)
      y2o=coord(8)
      z2o=coord(9)

      return
      end
      logical function pscan(aid,x,y,z)

*** scan the pos file for 4 char id and get x,y,z

      implicit double precision(a-h,o-z)
      character*80 card1,card2
      character*4 aid
      character*7 asla,aslo,aeht
      character*1 adla,adlo
      common/lus/lst,lpos,lrin1,lrin2,leph,lobs
      common/stugg/pi,rad,sol,f1,f2,w1,w2,we

*** loop until end of file -- cards must be pairwise

  100 read(lpos,1,end=777) card1
      read(lpos,1)         card2
    1 format(a80)
      if(card1(7:10).ne.'*80*') stop 43588
      if(card2(7:10).ne.'*86*') stop 43589

      if(aid.eq.card1(1:4)) then
        if(card1(45:55).eq.'           ' .or.
     *     card1(57:68).eq.'            '.or.
     *     card2(46:52).eq.'       '     ) then
          x=0.d0
          y=0.d0
          z=0.d0
        else
          read(card1,51) idla,imla,asla,adla,idlo,imlo,aslo,adlo
   51     format(44x,2i2,a7,a1,i3,i2,a7,a1)
          read(card2,52) aeht
   52     format(45x,a7)

          call nblank(asla,5,iblk)
          call nblank(aslo,5,iblk)
          read(asla,54) isla
          read(aslo,54) islo
   54     format(i7)

          call nblank(aeht,3,iblk)
          read(aeht,55) eht
   55     format(f7.3)

*** latitudes are positive north in decimal degrees
*** longitudes are positive east in decimal degrees [0 -> 360]
 
          dla=dble(idla)+imla/60.d0+isla/3600.d5
          if(adla.eq.'s'.or.adla.eq.'S') dla=     -dla
*****     if(adla.eq.'s'.or.adla.eq.'S') dla=90.d0-dla  !*** awful bug
          dlo=dble(idlo)+imlo/60.d0+islo/3600.d5
          if(adlo.eq.'w'.or.adlo.eq.'W') dlo=360.d0-dlo
 
*** convert to radians and get x,y,z

          gla=dla/rad
          glo=dlo/rad
          call geoxyz(gla,glo,eht,x,y,z)
          write(lst,*) 'precise position found for ',aid
        endif
        pscan=.true.
        return
      endif
      go to 100

*** end of file encountered -- marker id not found

  777 pscan=.false.

      return
      end
      subroutine getoff(coord,aid1,aid2)

*** get antenna offsets, establish antenna coordinates
*** offset file may contain 2 records in either order
*** warning -- this version defaults to no l1-l2 if offset file absent

      implicit double precision(a-h,o-z)
      character*4 aid1,aid2,aid
      logical lexist
      dimension coord(20)
      common/lus/lst,lpos,lrin1,lrin2,leph,lobs

*** storage scheme for coord()
***
***  1  2  3       X Y Z  point 1 marker
***  4  5  6       X Y Z  point 1 antenna, L1
***  7  8  9       X Y Z  point 2 marker       (may be approx)
*** 10 11 12       X Y Z  point 2 antenna, L1  (may be approx)

*** 13 14 15 16    dl1,dl2,dn,de  point 1 offsets
*** 17 18 19 20    dl1,dl2,dn,de  point 2 offsets

*** default values -- marker is l1 phase center -- no l1-l2 data

      do 10 i=13,20
   10 coord(i)=0.d0

*** see if offset file exists -- if so, borrow the lpos unit

      inquire(file='off',exist=lexist)
      if(lexist) then
        close(lpos)
        open(lpos,file='off',form='formatted',status='old')

*** scan for first id -- if found, override defaults

  100   read(lpos,1,end=199) aid,dl1,dl2,dn,de
    1   format(a4,4f10.4)
        if(aid.eq.aid1) then
          coord(13)=dl1
          coord(14)=dl2
          coord(15)=dn
          coord(16)=de
          go to 199
        endif
        go to 100
  199   continue

*** now scan for second id

        rewind(lpos)
  200   read(lpos,1,end=299) aid,dl1,dl2,dn,de
        if(aid.eq.aid2) then
          coord(17)=dl1
          coord(18)=dl2
          coord(19)=dn
          coord(20)=de
          go to 299
        endif
        go to 200
  299   continue

*** return the logical unit back to coordinate file

        close(lpos)
        open(lpos,file='pos',form='formatted',status='old')
      endif

*** sanity check on marker coordinates -- just being extra sure

      x=coord(1)
      y=coord(2)
      z=coord(3)
      r2=x*x+y*y+z*z
      if(r2.lt.3.6d13) then
        write(lst,*) 'invalid point 1 xyz',x,y,z
        stop 23674
      endif

      x=coord(7)
      y=coord(8)
      z=coord(9)
      r2=x*x+y*y+z*z
      if(r2.lt.3.6d13) then
        write(lst,*) 'invalid point 2 xyz -- set to point 1'
        coord(7)=coord(1)
        coord(8)=coord(2)
        coord(9)=coord(3)
      endif

*** compute the l1 centers 

      call mktol1(coord)

      return
      end
      subroutine mktol1(coord)

*** compute l1 phase centers using marker coordinates

      implicit double precision(a-h,o-z)
      logical xyzgeo
      dimension coord(20)

*** storage scheme for coord()
***
***  1  2  3       X Y Z  point 1 marker
***  4  5  6       X Y Z  point 1 antenna, L1
***  7  8  9       X Y Z  point 2 marker       (may be approx)
*** 10 11 12       X Y Z  point 2 antenna, L1  (may be approx)

*** 13 14 15 16    dl1,dl2,dn,de  point 1 offsets
*** 17 18 19 20    dl1,dl2,dn,de  point 2 offsets

*** point 1 -- rotate offsets into ECF and add to marker coordinates.
 
      x=coord(1)
      y=coord(2)
      z=coord(3)
      dn=coord(15)
      de=coord(16)
      du=coord(13)

      if(.not.xyzgeo(x,y,z,gla,glo,eht)) stop 43588
      call reg(gla,glo,dx,dy,dz,dn,de,du)

      coord( 4)=coord( 1)+dx
      coord( 5)=coord( 2)+dy
      coord( 6)=coord( 3)+dz

*** point 2 -- rotate offsets into ECF and add to marker coordinates.
 
      x=coord(7)
      y=coord(8)
      z=coord(9)
      dn=coord(19)
      de=coord(20)
      du=coord(17)

      if(.not.xyzgeo(x,y,z,gla,glo,eht)) stop 43589
      call reg(gla,glo,dx,dy,dz,dn,de,du)

      coord(10)=coord( 7)+dx
      coord(11)=coord( 8)+dy
      coord(12)=coord( 9)+dz

      return
      end
      subroutine l1tomk(coord)

*** compute point 2 marker coordinate from l1 phase center
*** note: inverse function of mktol1(), but only for point 2

      implicit double precision(a-h,o-z)
      logical xyzgeo
      dimension coord(20)

*** storage scheme for coord()
***
***  1  2  3       X Y Z  point 1 marker
***  4  5  6       X Y Z  point 1 antenna, L1
***  7  8  9       X Y Z  point 2 marker       (should be good)
*** 10 11 12       X Y Z  point 2 antenna, L1  (should be good)

*** 13 14 15 16    dl1,dl2,dn,de  point 1 offsets
*** 17 18 19 20    dl1,dl2,dn,de  point 2 offsets

*** point 2 -- rotate offsets into ECF and subtract from l1 center
 
      x=coord(10)
      y=coord(11)
      z=coord(12)
      dn=coord(19)
      de=coord(20)
      du=coord(17)

      if(.not.xyzgeo(x,y,z,gla,glo,eht)) stop 43587
      call reg(gla,glo,dx,dy,dz,dn,de,du)

      coord( 7)=coord(10)-dx
      coord( 8)=coord(11)-dy
      coord( 9)=coord(12)-dz

      return
      end
***********************************************************************t
*** core geodetic routines *********************************************
************************************************************************
      subroutine rge(gla,glo,u,v,w,x,y,z)
 
*** given a rectangular cartesian system (x,y,z)
*** compute a geodetic h cartesian sys   (u,v,w)
 
      implicit double precision(a-h,o-z)
 
      sb=dsin(gla)
      cb=dcos(gla)
      sl=dsin(glo)
      cl=dcos(glo)

      u=-sb*cl*x-sb*sl*y+cb*z
      v=-   sl*x+   cl*y
      w= cb*cl*x+cb*sl*y+sb*z

      return
      end
      subroutine reg(gla,glo,x,y,z,u,v,w)

*** given geodetic  cartesian coordinates  (u,v,w)  (lhs)
*** compute rectangular cartesian coordinates (x,y,z)  (rhs)
 
      implicit double precision(a-h,o-z)
 
      sb=dsin(gla)
      cb=dcos(gla)
      sl=dsin(glo)
      cl=dcos(glo)
 
      x=-sb*cl*u-sl*v+cb*cl*w
      y=-sb*sl*u+cl*v+cb*sl*w
      z= cb   *u      +sb   *w
 
      return
      end
      subroutine geoxyz(gla,glo,eht,x,y,z)
 
*** convert geodetic lat, long, ellip ht. to x,y,z
*** input units of radians, output in meters
******NOTE: e2 is from wgs84  ***********CHECK THIS
 
      implicit double precision(a-h,o-z)
      parameter(a=6378137.d0,e2=0.00669437999013d0)
 
      sla=dsin(gla)
      cla=dcos(gla)
      w2=1.d0-e2*sla*sla
      w=dsqrt(w2)
      en=a/w
 
      x=(en+eht)*cla*dcos(glo)
      y=(en+eht)*cla*dsin(glo)
      z=(en*(1.d0-e2)+eht)*sla
 
      return
      end
      logical function xyzgeo(x,y,z,glat,glon,eht)
 
*** convert x,y,z into geodetic lat, lon, and ellip. ht
*** input units of meters, output in radians
******NOTE: e2 is from wgs84  ***********CHECK THIS
*** ref: eq a.4b, p. 132, appendix a, osu #370
*** ref: geom geod notes gs 658, rapp
 
      implicit double precision(a-h,o-z)
      parameter(maxint=10,tol=5.d-15)
      parameter(a=6378137.d0,e2=0.00669437999013d0)
 
      ae2=a*e2
 
*** compute initial estimate of reduced latitude  (eht=0)
 
      p=dsqrt(x*x+y*y)
      icount=0
      tgla=z/p/(1.d0-e2)
 
*** iterate to convergence, or to max # iterations
 
    1 if(icount.le.maxint) then
        tglax=tgla
        tgla=z/(p-(ae2/dsqrt(1.d0+(1.d0-e2)*tgla*tgla)))
        icount=icount+1
        if(dabs((tgla-tglax)/tgla).gt.tol) go to 1
 
*** convergence achieved
 
        xyzgeo=.true.
        glat=datan(tgla)
        slat=dsin(glat)
        clat=dcos(glat)
        glon=datan2(y,x)
        w=dsqrt(1.d0-e2*slat*slat)
        en=a/w
        if(dabs(glat).le.0.7854d0) then
          eht=p/clat-en
        else
          eht=z/slat-en+e2*en
        endif
        glon=datan2(y,x)
 
*** too many iterations
 
      else
        xyzgeo=.false.
        glat=0.d0
        glon=0.d0
        eht=0.d0
      endif
 
      return
      end
      subroutine nblank(a,inum,iblk)
 
*** return precision of number and zero fill
 
      character*(*) a
      logical count
 
      leng=len(a) 
      l1=leng-inum+1
      iblk=0
      count=.true.
 
      do 1 i=leng,l1,-1
        if(a(i:i).eq.' ') then
          if(count) iblk=iblk+1
          a(i:i)='0'
        else
          count=.false.
        endif    
    1 continue
 
      return
      end  
************************************************************************
*** orbit load and interpolation ***************************************
************************************************************************
      subroutine lodorb(orbits,lorbit,nsat,nepoch,tstrt,tstop)

*** load ephmeris in sp3 ascii format
*** load only interval and satellites required

      implicit double precision(a-h,o-z)
      parameter(tiny=1.d-6)
      character*60 card
      logical getprn,delprn,scanbf
      parameter(ltab=32)
      dimension itlate(ltab),iprbuf(34)
      dimension orbits(lorbit)
      save /satstf/
      common/satstf/tstr,tstp,dlt,t0mx,imx,npt,n2,nst,nep,it,ix,iy,iz,id
      common/lus/lst,lpos,lrin1,lrin2,leph,lobs

*** skip over the sp3 format header
*** nloop is the number of satellites for each epoch, prns saved in ibuf

      read(leph,'(a60)') card
      read(leph,'(a60)') card
      read(leph,4) nloop,iprbuf
    4 format(4x,i2,3x,17i3/9x,17i3)
      do 1 i=5,22
    1 read(leph,'(a60)') card

*** if prn in prn table but not in orbit iprbuf(), delete prn from table
***   (situations where satellite not in orbit due to maneuver)
***   ltab is the maximum value for a prn  (32 possible sv prns)

      do 6 i=1,ltab
        if(getprn(i,iiprn).and..not.scanbf(i,iprbuf,nloop)) then
          if(.not.delprn(i,nsat)) stop 78439
          write(lst,5) nsat
    5     format('lodorb reset number of sats ==> ',i3)
        endif
    6 continue

*** now safe to build translation table for fast lookups

      do 2 i=1,ltab
        if(.not.getprn(i,iiprn)) then
          itlate(i)=0
        else
          itlate(i)=iiprn
        endif
    2 continue

*** compute indicies for the orbit data block and load common

      it    =1
      ix    =it+     nepoch
      iy    =ix+nsat*nepoch
      iz    =iy+nsat*nepoch
      id    =iz+nsat*nepoch
      itotal=id+nsat*nepoch-1
      nst=nsat
      nep=nepoch

*** initialize the orbit data block to 0   (for safety)

      do 3 i=1,itotal
    3 orbits(i)=0.d0

*** master loop over ephemeris file
*** read the epoch headers, then find position records
*** exit the loop after loading 'nepoch' epochs   (counted by nt)

      nt=0
      tstr=+1.d38
      tstp=-1.d38
  100 read(leph,'(a60)',end=777) card
      if(card(1:1).eq.'*') then
        read(card,103) iyr,imo,idy,ihr,imn,sec
  103   format(3x,i4,4i3,f12.8)
        call civjts(iyr,imo,idy,ihr,imn,sec,tsec)

*** check if within obs data common interval   (tiny for safety)
*** gps time stored in initial part of orbit data block (locations 1 --> nepoch)
*** module assumes orbits in monotonically ascending time order

        if(tsec.lt.(tstrt-tiny).or.tsec.gt.(tstop+tiny)) then
          continue
        else
          nt=nt+1
          orbits(nt)=tsec
          if(tsec.lt.tstr) tstr=tsec
          if(tsec.gt.tstp) tstp=tsec

*** bypass satellites not observed  (indicated by 0 entry in table)
*** immediately make ephemeris data into meters and seconds

          do 10 i=1,nloop
  101       read(leph,'(a60)') card
            if(card(1:1).eq.'*') stop 89756
            if(card(1:1).ne.'P'.and.card(1:1).ne.'p') go to 101
            read(card,102) is,x,y,z,dt
  102       format(1x,i3,4f14.6)
            x=x*1000.d0
            y=y*1000.d0
            z=z*1000.d0
            dt=dt/1.d6

            iis=itlate(is)
            if(iis.lt.0.or.iis.gt.ltab) stop 89757
            if(iis.ne.0) then
              call lodor2(nt,iis,x,y,z,dt,
     *                    orbits(ix),orbits(iy),orbits(iz),orbits(id))
            endif
   10     continue
          if(nt.eq.nepoch) go to 777
        endif
      endif
      go to 100
 
*** end of file encountered
 
  777 continue

*** set remaining constants for index work

      dlt=dble( idnint((tstp-tstr)/(nep-1)) )
      npt=9
      n2=(npt-1)/2
      imx=nep-npt+1
      t0mx=(imx-1)*dlt+tstr

*** initialize interpolation coefficients  (hard coded for 9)

      call itrini

      return
      end
      subroutine lodor2(nt,iis,x,y,z,dt,xs,ys,zs,ds)

*** place data into storage locations in orbit data block

      implicit double precision(a-h,o-z)
      dimension xs(nep,nst),ys(nep,nst),zs(nep,nst),ds(nep,nst)
      save /satstf/
      common/satstf/tstr,tstp,dlt,t0mx,imx,npt,n2,nst,nep,it,ix,iy,iz,id

*** sanity checks on indicies in lodorb()

      xs(nt,iis)=x
      ys(nt,iis)=y
      zs(nt,iis)=z
      ds(nt,iis)=dt

      return
      end
      logical function scanbf(iprn,iprbuf,n)

*** scan buffer of orbit header prn numbers for iprn

      implicit double precision(a-h,o-z)
      dimension iprbuf(n)

      do 1 i=1,n
        if(iprn.eq.iprbuf(i)) then
          scanbf=.true.
          return
        endif
    1 continue

*** fell thru loop -- iprn not found

      scanbf=.false.

      return
      end
      logical function oterp(iprn,tsec,orbits,x,y,z,dt,ierr)

*** interpolation of orbit
*** iprn is prn identifier for a satellite
*** input time, tsec, is gps time (maintained by usaf) w.r.t. mjd0
*** output (xyz) units of m in orbit reference system
*** output (dt)  units of sec
*** data in /satstf/ loaded by lodsat()

      implicit double precision(a-h,o-z)
      logical getprn,oterp2
      dimension orbits(*)
      save /satstf/
      common/satstf/tstr,tstp,dlt,t0mx,imx,npt,n2,nst,nep,it,ix,iy,iz,id

*** get internal sequence number of a satellite prn

      if(.not.getprn(iprn,iiprn)) then
        ierr=1
        oterp=.false.
        return
      endif

*** continue with interpolation

      if(.not.oterp2(iiprn,tsec,
     *           orbits(it),orbits(ix),orbits(iy),orbits(iz),orbits(id),
     *           x,y,z,dt,ierr)) then
        oterp=.false.
      else
        oterp=.true.
      endif

      return
      end
      logical function oterp2(is,tsec,ts,xs,ys,zs,ds,x,y,z,dt,ierr)

*** orbit interpolation subroutine  (called by oterp)

      implicit double precision(a-h,o-z)
      parameter(tiny=1.d-2)
      dimension ts(nep),xs(nep,nst),ys(nep,nst),zs(nep,nst),ds(nep,nst)
      save /satstf/
      common/satstf/tstr,tstp,dlt,t0mx,imx,npt,n2,nst,nep,it,ix,iy,iz,id
      common/lus/lst,lpos,lrin1,lrin2,leph,lobs

*** is  -- internal sequence number for a satellite (1 <= is <= nst)
*** nst -- number of satellite prns that were loaded
*** nep -- number of orbit epochs that were loaded

*** this line if bogus ******************

      ts(1)=ts(1)

*** satellite prn number not loaded from the sp3 file (1st epoch)

      r=xs(1,is)*xs(1,is)+ys(1,is)*ys(1,is)+zs(1,is)*zs(1,is)
      if(r.le.tiny) then
        ierr=2
        oterp2=.false.
        return
      endif

*** prevent extrapolation

      if(tsec.lt.tstr.or.tsec.gt.tstp) then
        ierr=3
        oterp2=.false.
        return
      endif

*** routine for finding index into table
*** produces normalized time [0,1,...,8] --> should be 3.5 < tnorm < 4.5

      i=(idnint((tsec-tstr)/dlt)+1)-n2
      t0=(i-1)*dlt+tstr

*** boundary conditions

      if(i.lt.1) then
        i=1
        t0=tstr
      elseif(i.gt.imx) then
        i=imx
        t0=t0mx
      endif
      t=(tsec-t0)/dlt

      call iterp9(t,xs(i,is),x)
      call iterp9(t,ys(i,is),y)
      call iterp9(t,zs(i,is),z)
      call iterp9(t,ds(i,is),dt)
      ierr=0
      oterp2=.true.

      return
      end
***********************************************************************
      subroutine iterp9(t,xs,x)

*** 9-th order interpolator                  -- equally-spaced points
*** time normalized 0,1,...,8                -- should be 3.5 < t < 4.5
*** coefficents, c(), initialized by itrini  -- see Ben Remondi dissertation

*** patched 19-jan-96 for t exactly equal to 4.0

      implicit double precision(a-h,o-z)
      parameter(eps=1.d-13)
      dimension xs(9)
      save /coeff9/
      common/coeff9/c(9)

*** return stored value if t is exact integer time (avoid division by 0)

      intt=idnint(t)
      tint=dble(intt)
      if(dabs(t-tint).le.eps) then
        if(intt.ge.0.and.intt.le.8) then
          x=xs(intt+1)
          return
        else
          stop 88998
        endif
      endif

*** normal interpolation

      powr=1.d0
      x=0.d0

      do 100 i=1,9
        dt=(t-dble(i-1))
        powr=dt*powr
        x=x+xs(i)/(c(i)*dt)
  100 continue
      x=x*powr/24.d0

      return
      end
***********************************************************************
*     subroutine iterp9(t,xs,x)

*** 9-th order interpolator                  -- equally-spaced points
*** time normalized 0,1,...,8                -- should be 3.5 < t < 4.5
*** coefficents, c(), initialized by itrini  -- see Ben Remondi dissertation

*     implicit double precision(a-h,o-z)
*     dimension xs(9)
*     save /coeff9/
*     common/coeff9/c(9)

*     powr=1.d0
*     x=0.d0

*     do 100 i=1,9
*       dt=(t-dble(i-1))
*       powr=dt*powr
*       x=x+xs(i)/(c(i)*dt)
* 100 continue
*     x=x*powr/24.d0

*     return
*     end
***********************************************************************
      subroutine itrini

*** initialize coefficients for 9-th order "lagrangean" interpolator
*** follows theory of Ben Remondi  (c.f. his dissertation)

      implicit double precision(a-h,o-z)
      save /coeff9/
      common/coeff9/c(9)

      c(1)=40320.d0 /24.d0
      c(2)=-5040.d0 /24.d0
      c(3)= 1440.d0 /24.d0
      c(4)= -720.d0 /24.d0
      c(5)=  576.d0 /24.d0
      c(6)= -720.d0 /24.d0
      c(7)= 1440.d0 /24.d0
      c(8)=-5040.d0 /24.d0
      c(9)=40320.d0 /24.d0

      return
      end
************************************************************************
*** time conversion ****************************************************
************************************************************************
      subroutine setjd0(iyr,imo,idy)

*** set the integer part of a modified julian date as epoch, mjd0
*** the modified julian day is derived from civil time as in civmjd()
*** allows single number expression of time in seconds w.r.t. mjd0

      implicit double precision(a-h,o-z)
      integer y
      save /mjdoff/
      common/mjdoff/mjd0

      if(iyr.lt.1900) stop 34587

      if(imo.le.2) then
        y=iyr-1
        m=imo+12
      else
        y=iyr
        m=imo
      endif

      it1=365.25d0*y
      it2=30.6001d0*(m+1)
      mjd=it1+it2+idy-679019

*** now set the epoch for future time computations

      mjd0=mjd

      return
      end
      subroutine civjts(iyr,imo,idy,ihr,imn,sec,tsec)

*** convert civil date to time in seconds past mjd epoch, mjd0
*** requires initialization of mjd0 by setjd0()

*** imo in range 1-12, idy in range 1-31
*** only valid in range mar-1900 thru feb-2100     (leap year protocols)
*** ref: hofmann-wellenhof, 2nd ed., pg 34-35
*** adapted from civmjd()

      implicit double precision(a-h,o-z)
      integer y
      save /mjdoff/
      common/mjdoff/mjd0

      if(iyr.lt.1900) stop 34589

      if(imo.le.2) then
        y=iyr-1
        m=imo+12
      else
        y=iyr
        m=imo
      endif

      it1=365.25d0*y
      it2=30.6001d0*(m+1)
      mjd=it1+it2+idy-679019

      tsec=(mjd-mjd0)*86400.d0+3600*ihr+60*imn+sec

      return
      end
      subroutine jtsciv(tsec,iyr,imo,idy,ihr,imn,sec)

*** convert time in seconds past mjd0 epoch into civil date
*** requires initialization of mjd0 by setjd0()

*** imo in range 1-12, idy in range 1-31
*** only valid in range mar-1900 thru feb-2100
*** ref: hofmann-wellenhof, 2nd ed., pg 34-35
*** adapted from mjdciv()

      implicit double precision(a-h,o-z)
      save /mjdoff/
      common/mjdoff/mjd0

      mjd=mjd0+tsec/86400.d0
*** the following equation preserves significant digits
      fmjd=dmod(tsec,86400.d0)/86400.d0

      rjd=mjd+fmjd+2400000.5d0
      ia=(rjd+0.5d0)
      ib=ia+1537
      ic=(ib-122.1d0)/365.25d0
      id=365.25d0*ic
      ie=(ib-id)/30.6001d0

*** the fractional part of a julian day is (fractional mjd + 0.5)
*** therefore, fractional part of julian day + 0.5 is (fractional mjd)

      it1=ie*30.6001d0
      idy=ib-id-it1+fmjd
      it2=ie/14.d0
      imo=ie-1-12*it2
      it3=(7+imo)/10.d0
      iyr=ic-4715-it3

      tmp=fmjd*24.d0
      ihr=tmp
      tmp=(tmp-ihr)*60.d0
      imn=tmp
      sec=(tmp-imn)*60.d0

      return
      end
********************************************************************8
      subroutine civmjd(iyr,imo,idy,ihr,imn,sec,mjd,fmjd)

*** convert civil date to modified julian date

*** imo in range 1-12, idy in range 1-31
*** only valid in range mar-1900 thru feb-210      (leap year protocols)
*** ref: hofmann-wellenhof, 2nd ed., pg 34-35
*** operation confirmed against table 3.3 values on pg.34

      implicit double precision(a-h,o-z)
      integer y

      if(iyr.lt.1900) stop 34588

      if(imo.le.2) then
        y=iyr-1
        m=imo+12
      else
        y=iyr
        m=imo
      endif

      it1=365.25d0*y
      it2=30.6001d0*(m+1)
      mjd=it1+it2+idy-679019

      fmjd=(3600*ihr+60*imn+sec)/86400.d0

      return
      end
      subroutine mjdciv(mjd,fmjd,iyr,imo,idy,ihr,imn,sec)

*** convert modified julian date to civil date

*** imo in range 1-12, idy in range 1-31
*** only valid in range mar-1900 thru feb-2100
*** ref: hofmann-wellenhof, 2nd ed., pg 34-35
*** operation confirmed for leap years (incl. year 2000)

      implicit double precision(a-h,o-z)

      rjd=mjd+fmjd+2400000.5d0
      ia=(rjd+0.5d0)
      ib=ia+1537
      ic=(ib-122.1d0)/365.25d0
      id=365.25d0*ic
      ie=(ib-id)/30.6001d0

*** the fractional part of a julian day is fractional mjd + 0.5
*** therefore, fractional part of julian day + 0.5 is fractional mjd

      it1=ie*30.6001d0
      idy=ib-id-it1+fmjd
      it2=ie/14.d0
      imo=ie-1-12*it2
      it3=(7+imo)/10.d0
      iyr=ic-4715-it3

      tmp=fmjd*24.d0
      ihr=tmp
      tmp=(tmp-ihr)*60.d0
      imn=tmp
      sec=(tmp-imn)*60.d0

      return
      end
************************************************************************
*** prn table handlers *************************************************
************************************************************************
      subroutine newprn 
  
*** initialize the table to zero length 
  
      parameter(maxprn=32)
      save /prntab/
      integer prns
      common/prntab/nprns,prns(maxprn) 
  
      nprns=0
  
      return
      end 
      subroutine numprn(isize)
  
*** return current size of table
  
      parameter(maxprn=32)
      save /prntab/
      integer prns
      common/prntab/nprns,prns(maxprn) 
  
      isize=nprns
  
      return
      end 
      logical function getprn(prn,iprn) 
  
*** table lookup
  
      parameter(maxprn=32)
      integer prns,prn 
      save /prntab/
      common/prntab/nprns,prns(maxprn) 
  
      iprn=1 
  100 if(iprn.le.nprns) then
        if(prn.eq.prns(iprn)) then 
          getprn=.true. 
          return
        else
          iprn=iprn+1 
        endif 
        goto 100
      endif 
  
*** fell thru table 
  
      getprn=.false.
  
      return
      end 
      logical function putprn(prn,idup)
  
*** add entry to table
*** idup is location of duplicate entry in list 
*** idup=0 indicates exceeded maximum length of list
  
      parameter(maxprn=32)
      integer prns,prn 
      logical getprn
      save /prntab/
      common/prntab/nprns,prns(maxprn) 
  
      if(getprn(prn,iprn)) then 
        idup=iprn
        putprn=.false.
        return
      else
        idup=0
        if(iprn.gt.maxprn) then
          putprn=.false.
          return
        else
          nprns=iprn
          prns(nprns)=prn
          putprn=.true. 
        endif 
      endif 
  
      return
      end 
      logical function delprn(prn,isize)
  
*** delete entry in table, return new size of table
  
      parameter(maxprn=32)
      integer prns,prn 
      logical getprn
      save /prntab/
      common/prntab/nprns,prns(maxprn) 
  
      if(getprn(prn,iprn)) then 
        do 1 i=iprn+1,nprns
    1   prns(i-1)=prns(i)
        nprns=nprns-1

        delprn=.true.
        isize=nprns
        return

*** entry not in table -- no action

      else
        delprn=.false.
        isize=nprns
      endif 
  
      return
      end 
      integer function invprn(iprn) 
  
*** inverse table lookup
  
      parameter(maxprn=32)
      integer prns
      save /prntab/
      common/prntab/nprns,prns(maxprn) 

*** given the index, return the prn number

      if(iprn.lt.1.or.iprn.gt.nprns) stop 74358
  
      invprn=prns(iprn)
  
      return
      end 
************************************************************************
*** reference satellite table handlers *********************************
************************************************************************
      subroutine openrs

*** initialize reference satellite table

      implicit double precision(a-h,o-z)
      parameter(maxrst=100)
      save /rstab/
      common/rstab/nrs,irs,tsvrfs(maxrst),isvrfs(maxrst)

      nrs=0
      irs=0

      return
      end
      subroutine sizers(n)

*** report size of reference satellite table

      implicit double precision(a-h,o-z)
      parameter(maxrst=100)
      save /rstab/
      common/rstab/nrs,irs,tsvrfs(maxrst),isvrfs(maxrst)

      n=nrs

      return
      end
      subroutine dumprs

*** dump reference satellite table to output

      implicit double precision(a-h,o-z)
      parameter(maxrst=100)
      save /rstab/
      common/rstab/nrs,irs,tsvrfs(maxrst),isvrfs(maxrst)
      common/lus/lst,lpos,lrin1,lrin2,leph,lobs

***** THIS IS DEBUG STUFF -- debug
***** write(lst,*) 'WARNING!! REF. SAT. RESET!!! -- 1'
***** isvrfs(1)=1
***** THIS IS DEBUG STUFF -- debug

*** report external sat numbers -- not internal

      write(lst,*) 'dump ref. sat. table  nrs =',nrs
      write(lst,'(f9.2,i6)') tsvrfs(1),invprn(isvrfs(1))
      do 1 i=2,nrs
    1 write(lst,'(f9.2,i6,f12.2)') tsvrfs(i),invprn(isvrfs(i)),
     *                             tsvrfs(i)-tsvrfs(i-1)

      return
      end
      logical function putrs(isv,tsv,n)

*** add reference satellite to end of list

      implicit double precision(a-h,o-z)
      parameter(maxrst=100)
      save /rstab/
      common/rstab/nrs,irs,tsvrfs(maxrst),isvrfs(maxrst)

*** first entry

      if(nrs.le.0) then
        nrs=1
        irs=1
        tsvrfs(nrs)=tsv
        isvrfs(nrs)=isv

*** subsequent entry -- test if new sv

      elseif(isv.ne.isvrfs(nrs)) then
        if(nrs.ge.maxrst) then
          putrs=.false.
        endif
        nrs=nrs+1
        isvrfs(nrs)=isv
        tsvrfs(nrs)=tsv

*** same sv, so just update the time

      else
        tsvrfs(nrs)=tsv
      endif

*** return current size of table

      n=nrs
      putrs=.true.

      return
      end
      logical function getrs(isv,tsv)

*** get reference satellite based on time
*** 'irs' is current index -- used to optimize sequential search
*** note: receivers should be synchronized to within "tiny"

      implicit double precision(a-h,o-z)
      parameter(tiny=0.2d0)
      parameter(maxrst=100)
      save /rstab/
      common/rstab/nrs,irs,tsvrfs(maxrst),isvrfs(maxrst)

      if(nrs.le.0) stop 13598
      if(nrs.eq.1) then
        isv=isvrfs(nrs)
        getrs=.true.
        return
      endif

*** check for index reset -- look at prior time

      if(irs.ge.2) then
        if(tsv.lt.tsvrfs(irs-1)) irs=1
      endif

*** search current location

  100 if(tsv.le.tsvrfs(irs)+tiny) then
        isv=isvrfs(irs)
        getrs=.true.
        return

*** increment the index --- retry

      elseif(irs.lt.nrs) then
        irs=irs+1
        go to 100
      endif

*** table exhausted

      getrs=.false.

      return
      end
************************************************************************
*** matrix mashers *****************************************************
************************************************************************
      subroutine writev(a,m)

*** routine to write a vector 

      implicit double precision(a-h,o-z)
      dimension a(m)
      common/lus/lst,lpos,lrin1,lrin2,leph,lobs

      do 10 i=1,m 
   10 write(lst,1) a(i) 
    1 format(' ',1pd30.20)
      return
      end 
      subroutine writem(a,m,n)
 
*** routine to write a matrix 

      implicit double precision(a-h,o-z)
      dimension a(m,n)
      common/lus/lst,lpos,lrin1,lrin2,leph,lobs

      do 10 i=1,m 
   10 write(lst,1) (a(i,j),j=1,n) 
    1 format(' ',1p4d18.10) 
      return
      end 
      subroutine innerp(x,n,val)

*** compute inner product of vector

      implicit double precision(a-h,o-z)
      dimension x(n)

      val=0.d0
      do 1 i=1,n
    1 val=val+x(n)*x(n)

      return
      end
      subroutine ab(a,b,r,l,m,n)
 
*** form the matrix product r=ab 
*** the matrices a and b are returned unchanged

      implicit double precision(a-h,o-z)
      dimension a(l,m),b(m,n),r(l,n)

      do 5 i=1,l
      do 5 j=1,n
      r(i,j)=0.d0
      do 5 k=1,m
    5 r(i,j)=r(i,j)+a(i,k)*b(k,j) 
      return
      end 
      subroutine invert(a,n)

*** routine to invert a matrix in place 

      implicit double precision(a-h,o-z)
      dimension a(n,n)

      do 2 i=1,n
      aii=1.d0/a(i,i)
      a(i,i)=aii
      do 2 j=1,n
      if(i.eq.j) go to 2
      aji=a(j,i)*aii
      a(j,i)=aji
      do 1 k=1,n
      if(i.eq.k) go to 1
      a(j,k)=a(j,k)-aji*a(i,k)
      if(j.ne.n) go to 1
      a(i,k)=-aii*a(i,k)
    1 continue
    2 continue
      ann=-a(n,n) 
      k=n-1 
      do 3 j=1,k
    3 a(n,j)=ann*a(n,j) 
      return
      end 
      subroutine nrmal2(an,u,d,f,p,bw,lc,nd,n) 

***  form normals for observation equations one observation at a time 
***            *****(nx + u = 0)***** 
***        n=a'pa         (n matrix stored in an) 
***        u=a'pl         (l stored in f) 
***  calling sequence 
***    1)  call nitial  (one time to clear normals) 
***    2)  call normal  (once for each observation) 
***    3)  call fill    (one time to fill in lower triangular)
*** 
***  d = non-zero elements of a row of the a matrix 
***  lc= column numbers of non-zero elements in d.  elements
***      of lc must be in ascending order and correspond to 
***      d, i.e., d(j) is the non-zero element in column lc(j)
***  nd= number of elements in d and lc 
***  n = order of the system (number of unknowns) 
***  bw= l'pl accumulated (basement window)

      implicit double precision(a-h,o-z) 
      dimension an(n,n),u(n),d(nd),lc(nd) 

      pf=p*f
      bw=bw+f*pf
      do 5 i=1,nd 
      li=lc(i)
      u(li)=u(li)+d(i)*pf 
      do 5 j=i,nd 
      lj=lc(j)
    5 an(li,lj)=an(li,lj) + d(i)*d(j)*p 
      return

*** clear normals

      entry nitil2(an,u,bw,n)
      do 10 i=1,n 
      u(i)=0.d0
      do 10 j=i,n 
   10 an(i,j)=0.d0 
      bw=0.d0 
      return

*** fill in the lower triangular portion of an 

      entry fill(an,n)
      do 15 i=1,n 
      do 15 j=i,n 
   15 an(j,i)=an(i,j) 
      return
      end 
