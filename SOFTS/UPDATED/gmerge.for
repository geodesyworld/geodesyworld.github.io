      program gmerge

*** merge and reformat data

*** fixed nasty bug in pscan()

      implicit double precision(a-h,o-z)
      dimension coord(2,17),icols1(0:4),icols2(0:4)
      dimension tena1(2,19),tena2(2,19)
      character*16 ant1,ant2
      character*4 aid1,aid2
      logical xyzgeo
      common/lus/lst,lpos,lrin1,lrin2
      common/stugg/pi,rad,sol,f1,f2,w1,w2,we

      pi=4.d0*datan(1.d0)
      rad=180.d0/pi

      sol=299792458.d0
      f1 =1575.42d6
      f2 =1227.60d6
      w1 =sol/f1
      w2 =sol/f2
      we=7.2921151467d-5

      write(*,*) 'gmerge -- 2005sep17'

*** debug: input file names currently hard-coded

      lst=1
      open(lst,file='gmerge.out',form='formatted')

      lpos=10
      open(lpos,file='pos',form='formatted',status='old')

      lrin1=2
      open(lrin1,file='point1.rnx',form='formatted',status='old')

      lrin2=3
      open(lrin2,file='point2.rnx',form='formatted',status='old')

*** gather critical information

      call recon(nsat,dstrt,dstop,delt,
     *           icols1,ndat1,aid1,ant1,icols2,ndat2,aid2,ant2,coord)

*** get precise coordinates from the lpos file

      call lodpos(coord,aid1,aid2)

*** get arp and phase coordinates

      call getoff(coord,ant1,tena1,ant2,tena2)

      if(.not.xyzgeo(coord(1, 1),coord(1, 2),coord(1, 3),
     *   glat,glon,eht)) stop 77689
      if(glon.lt.pi) glon=glon+pi+pi
      write(lst,3) aid1,       coord(1, 1),coord(1, 2),coord(1, 3),
     *           glat*rad,glon*rad,eht
    3 format(a4,12x,3f14.4,2f14.9,f10.4)

      if(.not.xyzgeo(coord(1, 4),coord(1, 5),coord(1, 6),
     *   glat,glon,eht)) stop 77688
      if(glon.lt.pi) glon=glon+pi+pi
      write(lst,1) ant1,       coord(1, 4),coord(1, 5),coord(1, 6),
     *           glat*rad,glon*rad,eht
    1 format(a16,3f14.4,2f14.9,f10.4)

      if(.not.xyzgeo(coord(1, 7),coord(1, 8),coord(1, 9),
     *   glat,glon,eht)) stop 77687
      if(glon.lt.pi) glon=glon+pi+pi
      write(lst,2) coord(1,16),coord(1, 7),coord(1, 8),coord(1, 9),
     *           glat*rad,glon*rad,eht

      if(.not.xyzgeo(coord(1,10),coord(1,11),coord(1,12),
     *   glat,glon,eht)) stop 77686
      if(glon.lt.pi) glon=glon+pi+pi
      write(lst,2) coord(1,17),coord(1,10),coord(1,11),coord(1,12),
     *           glat*rad,glon*rad,eht
    2 format(f7.4,9x,3f14.4,2f14.9,f10.4)

      write(lst,'(19i5)') (idnint(10000.d0*tena1(1,j)),j=1,19)
      write(lst,'(19i5)') (idnint(10000.d0*tena1(2,j)),j=1,19)

      if(.not.xyzgeo(coord(2, 1),coord(2, 2),coord(2, 3),
     *   glat,glon,eht)) stop 77685
      if(glon.lt.pi) glon=glon+pi+pi
      write(lst,3) aid2,       coord(2, 1),coord(2, 2),coord(2, 3),
     *           glat*rad,glon*rad,eht

      if(.not.xyzgeo(coord(2, 4),coord(2, 5),coord(2, 6),
     *   glat,glon,eht)) stop 77684
      if(glon.lt.pi) glon=glon+pi+pi
      write(lst,1) ant2,       coord(2, 4),coord(2, 5),coord(2, 6),
     *           glat*rad,glon*rad,eht

      if(.not.xyzgeo(coord(2, 7),coord(2, 8),coord(2, 9),
     *   glat,glon,eht)) stop 77683
      if(glon.lt.pi) glon=glon+pi+pi
      write(lst,2) coord(2,16),coord(2, 7),coord(2, 8),coord(2, 9),
     *           glat*rad,glon*rad,eht

      if(.not.xyzgeo(coord(2,10),coord(2,11),coord(2,12),
     *   glat,glon,eht)) stop 77682
      if(glon.lt.pi) glon=glon+pi+pi
      write(lst,2) coord(2,17),coord(2,10),coord(2,11),coord(2,12),
     *           glat*rad,glon*rad,eht

      write(lst,'(19i5)') (idnint(10000.d0*tena2(1,j)),j=1,19)
      write(lst,'(19i5)') (idnint(10000.d0*tena2(2,j)),j=1,19)


*** reformat data

      write(lst,'(f9.3,f10.3,f6.1)') dstrt,dstop,delt
      call reform(nsat,dstrt,dstop,delt,icols1,ndat1,icols2,ndat2)

      stop
      end
************************************************************************
*** reformat data ******************************************************
************************************************************************
      subroutine reform(nsat,tstrt,tstop,delt,icols1,ndat1,icols2,ndat2)

*** merge and reformat data

      implicit double precision(a-h,o-z)
      character*80 card
      external invprn
      logical oget,dget,dget2
      parameter(dnull=1.d31)
      parameter(maxprn=32)
      dimension ibuf1(maxprn)  ,ibuf2(maxprn)
      dimension dbuf1(maxprn*4),dbuf2(maxprn*4)
      dimension icols1(0:4),icols2(0:4)

      common/stugg/pi,rad,sol,f1,f2,w1,w2,we
      common/lus/lst,lpos,lrin1,lrin2

      logical l2lin1,l2lin2,li1,li2
      save  /flags/
      common/flags/l2lin1,l2lin2

*** scan past headers of the rinex obs files

      rewind(lrin1)
   20 read(lrin1,'(a80)') card
      if(card(61:73).ne.'END OF HEADER') go to 20

      rewind(lrin2)
   30 read(lrin2,'(a80)') card
      if(card(61:73).ne.'END OF HEADER') go to 30

*** repeat/until loop over rinex obs files
*** synchronization (+/- 0.2 sec) obtained using largest obs. interval

      tiny=0.2d0
      iflag1=0
      iflag2=0
      tmas=tstrt
      lu1=lrin1
      lu2=lrin2
      li1=l2lin1
      li2=l2lin2

*** back to original

   90 continue
      if(oget(lu1,li1,nsat,icols1,ndat1,tmas,trec1,ibuf1,dbuf1,iflag1)
     *   .and.
     *   oget(lu2,li2,nsat,icols2,ndat2,tmas,trec2,ibuf2,dbuf2,iflag2))
     *   then

*** only use epoch if reference satellite available (count epochs)
*** will settle for just L1 pseudorange

        nboth=0
        do 92 iisat=1,nsat
          if(ibuf1(iisat).eq.0) go to 92
          if(ibuf2(iisat).eq.0) go to 92
          if(dget2(iisat,dbuf1,dbuf2,nsat)) then
            nboth=nboth+1            
          endif
   92   continue

*** dget() order
***   r11  1st point, L1 range
***   r12  1st point, L2 range
***   r21  2nd point, L1 range
***   r22  2nd point, L2 range

***   p11  1st point, L1 phase
***   p12  1st point, L2 phase
***   p21  2nd point, L1 phase
***   p22  2nd point, L2 phase

*** debug realign data ???????????????????????????????????????????

        write(lst,'(f9.3,f10.3,i3)') trec1,trec2,nboth
        do 91 iisat=1,nsat
          if(dget(iisat,ibuf1,ibuf2,dbuf1,dbuf2,nsat,
     *            r11,r12,r21,r22,p11,p12,p21,p22)) then

      write(lst,76543) invprn(iisat),r11,r21,p11,p21,r12,r22,p12,p22
76543 format(i2,2(2f13.3,2f14.3))           

          endif
   91   continue

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

*** end of file raised -- proceed

   97 continue

      return
      end
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
      common/lus/lst,lpos,lrin1,lrin2

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
     *                      r11,r12,r21,r22,p11,p12,p21,p22)

*** see if valid phase data in both buffers for iisat-th satellite
*** note: update this in synch with dget2()

      implicit double precision(a-h,o-z)
      parameter(dnull=1.d31,dnull2=99999999.999d0)
      dimension ibuf1(nsat),dbuf1(nsat,4)
      dimension ibuf2(nsat),dbuf2(nsat,4)
      common/lus/lst,lpos,lrin1,lrin2

      if(iisat.le.0.or.iisat.gt.nsat)  stop 87473

*** assume the worst

      dget=.false.

      if(ibuf1(iisat).eq.0) return
      if(ibuf2(iisat).eq.0) return

*** order {l1, p1/c1, l2, p2}  (2nd index)
***   p11  1st point, l1 phase
***   p12  1st point, l2 phase
***   r11  1st point, l1 range
***   r12  1st point, l2 range

***   p21  2nd point, l1 phase
***   p22  2nd point, l2 phase
***   r21  2nd point, l1 range
***   r22  2nd point, l2 range

      r11=dbuf1(iisat,2)
      r12=dbuf1(iisat,4)
      r21=dbuf2(iisat,2)
      r22=dbuf2(iisat,4)

      p11=dbuf1(iisat,1)
      p12=dbuf1(iisat,3)
      p21=dbuf2(iisat,1)
      p22=dbuf2(iisat,3)

*** will settle for just L1 pseudorange

      if(r11.ge.dnull) return
      if(r21.ge.dnull) return

***   if(r12.ge.dnull) return
***   if(r22.ge.dnull) return

***   if(p11.ge.dnull) return
***   if(p21.ge.dnull) return
***   if(p12.ge.dnull) return
***   if(p22.ge.dnull) return

      if(r12.ge.dnull) r12=dnull2
      if(r22.ge.dnull) r22=dnull2

      if(p11.ge.dnull) p11=dnull2
      if(p21.ge.dnull) p21=dnull2
      if(p12.ge.dnull) p12=dnull2
      if(p22.ge.dnull) p22=dnull2

*** things are all right after all

      dget=.true.

      return
      end
      logical function dget2(iisat,dbuf1,dbuf2,nsat)

*** see if valid data in both data buffers for iisat-th satellite
*** note: update this in synch with dget()

      implicit double precision(a-h,o-z)
      parameter(dnull=1.d31)
      dimension dbuf1(nsat,4),dbuf2(nsat,4)

*** assume the worst

      dget2=.false.

*** order {l1, p1/c1, l2, p2}  (2nd index)
*** will settle for just L1 pseudorange

      if(dbuf1(iisat,2).ge.dnull) return
      if(dbuf2(iisat,2).ge.dnull) return

*** things are all right after all

      dget2=.true.

      return
      end
************************************************************************
*** initial scan routines **********************************************
************************************************************************
      subroutine recon(nsat,start,stop,delt,icols1,ndat1,aid1,ant1,
     *                                     icols2,ndat2,aid2,ant2,coord)

*** preliminary scan of files for critical information

      implicit double precision(a-h,o-z)
      character*80 card
      character*16 ant1,ant2
      character*4  aid1,aid2
      dimension coord(2,17),icols1(0:4),icols2(0:4)
      common/lus/lst,lpos,lrin1,lrin2

      logical l2lin1,l2lin2
      save /flags/
      common/flags/l2lin1,l2lin2

*** initialize coordinates and offsets

      do 3 i=1,2
      do 3 j=1,17
    3 coord(i,j)=0.d0

*** set initial epoch for all subsequent time computations

    1 read(lrin1,'(a80)') card
      if(card(61:73).ne.'END OF HEADER') go to 1
 
      read(lrin1,2) iyr,imo,idy
    2 format(3i3)
      if(iyr.lt.80) then
        iyr=iyr+2000
      elseif(iyr.le.99) then
        iyr=iyr+1900
      endif
      rewind(lrin1)

      call setjd0(iyr,imo,idy)
      write(lst,'(i4,2i3)') iyr,imo,idy

*** scan obs files, get start/stop times, obs interval, sat. list, pos.

      call newprn
      lrin=lrin1
      call oscan(lrin,aid1,ant1,tstrt1,tstop1,delt1,icols1,ndat1,
     *           x1,y1,z1,dh1,de1,dn1,l2lin1)
      write(*,'(a,3f10.2)') ' file1  start/stop, intrvl',
     *                                tstrt1,tstop1,delt1
      coord(1, 1)=x1
      coord(1, 2)=y1
      coord(1, 3)=z1
      coord(1,13)=dh1
      coord(1,14)=de1
      coord(1,15)=dn1
      rewind(lrin1)

      lrin=lrin2
      call oscan(lrin,aid2,ant2,tstrt2,tstop2,delt2,icols2,ndat2,
     *           x2,y2,z2,dh2,de2,dn2,l2lin2)
      write(*,'(a,3f10.2)') ' file2  start/stop, intrvl',
     *                                tstrt2,tstop2,delt2
      coord(2,1)=x2
      coord(2,2)=y2
      coord(2,3)=z2
      coord(2,13)=dh2
      coord(2,14)=de2
      coord(2,15)=dn2
      rewind(lrin2)

*** find common data interval
*** find common obs. interval (max of all obs. intervals)
*** insure data intervals are compatible with the obs interval

      start=dmax1(tstrt1,tstrt2)
      stop =dmin1(tstop1,tstop2)
      delt =dmax1(delt1,delt2)

      delt=dble(idnint(delt))
      start=delt*idnint(start/delt)
      stop =delt*idnint(stop /delt)

      write(*,'(a,3f10.0)') ' common start/stop, intrvl',
     *                         start,stop,delt
      if(stop.le.start) then
        write(lst,*) 'receiver data do not overlap -- fatal'
        stop 43587
      endif

      call numprn(nsat)

      return
      end
      subroutine oscan(lrin,aid,ant,tstrt,tstop,delt,icols,ndat,x,y,z,
     *                 dh,de,dn,l2line)

*** scan rinex obs file, accumulate start/stop times, obs. interval,
***    and satellite list

*** 4-mar-96: force multiple of 5 seconds if interval greater than 5

      implicit double precision(a-h,o-z)
      logical putprn,l2line
      character*80 card
      character*16 ant
      character*4 aid
      character*2 as(9)
      dimension icols(0:4)
      parameter(mxsat=12)
      dimension iprns(mxsat)
      common/lus/lst,lpos,lrin1,lrin2

*** initialize approximate marker position and offsets

      x=0.d0
      y=0.d0
      z=0.d0

      dh=0.d0
      de=0.d0
      dn=0.d0

*** scan header for observation type description record

   10 read(lrin,1) card
    1 format(a80)
      
*** extract marker name

      if(card(61:71).eq.'MARKER NAME') then
        aid=card(1:4)

*** get antenna name

      elseif(card(61:72).eq.'ANT # / TYPE') then
        ant=card(21:40)

*** get approx xyz

      elseif(card(61:79).eq.'APPROX POSITION XYZ') then
        read(card,'(3f14.4)') x,y,z

*** get antenna dh,de,dn  (this is from mark to arp)

      elseif(card(61:80).eq.'ANTENNA: DELTA H/E/N') then
        read(card,'(3f14.4)') dh,de,dn

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
      if(nprn.le.0.or.nprn.gt.mxsat) then
        write(lst,*) '********** fatal error **********'
        write(lst,*) 'iyr,imo,idy,ihr,imn,sec,iepflg,nprn'
        write(lst,662) iyr,imo,idy,ihr,imn,sec,iepflg,nprn
  662   format(1x,i5,4i3,f11.7,i3,i3)
        stop 24308
      endif
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
      write(*,'(a)') card
      write(lst,*) 'lrin=',lrin,'  last valid data header was:'
      write(lst,'(i5,4i4,f9.3,a)') iyr,imo,idy,ihr,imn,sec,' was OK'
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
*** marker and antenna position ****************************************
************************************************************************
      subroutine lodpos(coord,aid1,aid2)

*** get precise coordinates from the lpos file
*** precise coordinates may not be available for point2

      implicit double precision(a-h,o-z)
      character*4 aid1,aid2
      logical pscan
      dimension coord(2,17)
      common/lus/lst,lpos,lrin1,lrin2

*** storage scheme for coord()
***
***  1  2  3       X Y Z  marker (may be arp if offset is 0)
***  4  5  6       X Y Z  arp
***  7  8  9       X Y Z  L1 phase center
*** 10 11 12       X Y Z  L2 phase center
*** 13 14 15       dh,de,dn offset  (arp - mark)
***       16       L1 offset  (L1 - arp)
***       17       L2 offset  (L2 - arp)

      rewind(lpos)
      if(.not.pscan(aid1,x,y,z)) then
        write(*,*) 'warning --',aid1,' -- id not found in pos file'
        write(*,*) 'this is not good -- point 2 unreliable'
      else
        coord(1,1)=x
        coord(1,2)=y
        coord(1,3)=z
      endif

*** point2

      rewind(lpos)
      if(.not.pscan(aid2,x,y,z)) then
        write(*,*) 'info    --',aid2,' -- id not found in pos file'
      else
        coord(2,1)=x
        coord(2,2)=y
        coord(2,3)=z
      endif

      return
      end
      logical function pscan(aid,x,y,z)

*** scan the pos file for 4 char id and get x,y,z

      implicit double precision(a-h,o-z)
      character*80 card1,card2
      character*4 aid
      character*7 asla,aslo,aeht
      character*1 adla,adlo
      common/lus/lst,lpos,lrin1,lrin2
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
          write(*,*) 'precise position found for ',aid
        endif
        pscan=.true.
        return
      endif
      go to 100

*** end of file encountered -- marker id not found

  777 pscan=.false.

      return
      end
      subroutine getoff(coord,ant1,tena1,ant2,tena2)

*** get offsets to complete coordinate sets

      implicit double precision(a-h,o-z)
      dimension coord(2,17),tena1(2,19),tena2(2,19)
      character*16 ant1,ant2
      common/lus/lst,lpos,lrin1,lrin2

*** sanity check on marker coordinates -- just being extra sure

      x=coord(1,1)
      y=coord(1,2)
      z=coord(1,3)
      r2=x*x+y*y+z*z
      if(r2.lt.3.6d13) then
        write(*,*) 'invalid point 1 xyz',x,y,z
        stop 23674
      endif

      x=coord(2,1)
      y=coord(2,2)
      z=coord(2,3)
      r2=x*x+y*y+z*z
      if(r2.lt.3.6d13) then
        write(*,*) 'invalid point 2 xyz -- set to point 1'
        coord(2,1)=coord(1,1)
        coord(2,2)=coord(1,2)
        coord(2,3)=coord(1,3)
      endif

*** reuse lpos for antenna file

      close(lpos)
      open(lpos,file='ant',form='formatted',status='old')
      call ascan(ant1,dn11,de11,du11,dn21,de21,du21,tena1)
      coord(1,16)=du11
      coord(1,17)=du21
      rewind(lpos)
      call ascan(ant2,dn12,de12,du12,dn22,de22,du22,tena2)
      coord(2,16)=du12
      coord(2,17)=du22
      close(lpos)

*** storage scheme for coord()
***
***  1  2  3       X Y Z  marker (may be arp if offset is 0)
***  4  5  6       X Y Z  arp
***  7  8  9       X Y Z  L1 phase center
*** 10 11 12       X Y Z  L1 phase center
*** 13 14 15       dh,de,dn offset  (arp - mark)
***       16       L1 offset  (L1 - arp)
***       17       L2 offset  (L2 - arp)

*** compute relative offsets

      call aoffst(coord(1, 1), coord(1, 2), coord(1, 3),
     *            coord(1, 4), coord(1, 5), coord(1, 6),
     *            coord(1,15), coord(1,14), coord(1,13))

      call aoffst(coord(1, 4), coord(1, 5), coord(1, 6),
     *            coord(1, 7), coord(1, 8), coord(1, 9),
     *                   dn11,        de11, coord(1,16))

      call aoffst(coord(1, 4), coord(1, 5), coord(1, 6),
     *            coord(1,10), coord(1,11), coord(1,12),
     *                   dn21,        de21, coord(1,17))

      call aoffst(coord(2, 1), coord(2, 2), coord(2, 3),
     *            coord(2, 4), coord(2, 5), coord(2, 6),
     *            coord(2,15), coord(2,14), coord(2,13))

      call aoffst(coord(2, 4), coord(2, 5), coord(2, 6),
     *            coord(2, 7), coord(2, 8), coord(2, 9),
     *                   dn12,        de12, coord(2,16))

      call aoffst(coord(2, 4), coord(2, 5), coord(2, 6),
     *            coord(2,10), coord(2,11), coord(2,12),
     *                   dn22,        de22, coord(2,17))

      return
      end
      subroutine ascan(ant,dn1,de1,du1,dn2,de2,du2,tena)

*** scan the antenna file for 20 char id and get data

      implicit double precision(a-h,o-z)
      character*16 ant
      dimension tena(2,19)
      character*80 card
      common/lus/lst,lpos,lrin1,lrin2

*** bypass header

      do 10 i=1,11
   10 read(lpos,1) card
    1 format(a80)

*** loop until antenna data found

  100 read(lpos,1,end=777) card
      if(card(1:16).eq.ant) then
        read(lpos,'(3f10.1)') dn1,de1,du1
        read(lpos,'(10f6.1)') (tena(1,j),j= 1,10)
        read(lpos,'( 9f6.1)') (tena(1,j),j=11,19)
        read(lpos,'(3f10.1)') dn2,de2,du2
        read(lpos,'(10f6.1)') (tena(2,j),j= 1,10)
        read(lpos,'( 9f6.1)') (tena(2,j),j=11,19)

*** rescale mm to meters

        dn1=dn1/1000.d0
        de1=de1/1000.d0
        du1=du1/1000.d0
        dn2=dn2/1000.d0
        de2=de2/1000.d0
        du2=du2/1000.d0
        do 20 i=1,2
        do 20 j=1,19
   20   tena(i,j)=tena(i,j)/1000.d0
        return
      else
        do 2 i=1,6
    2   read(lpos,1,end=777) card
      endif
      go to 100

*** end of file encountered -- marker id not found

  777 dn1=0.d0
      de1=0.d0
      du1=0.d0
      dn2=0.d0
      de2=0.d0
      du2=0.d0
      do 30 i=1,2
      do 30 j=1,19
   30 tena(i,j)=0.d0

      return
      end
      subroutine aoffst(x1,y1,z1,x2,y2,z2,dn,de,du)

*** rotate offsets into ECF and add to marker coordinates.

      implicit double precision(a-h,o-z)
      logical xyzgeo

      if(.not.xyzgeo(x1,y1,z1,gla,glo,eht)) stop 43588
      call reg(gla,glo,dx,dy,dz,dn,de,du)

      x2=x1+dx
      y2=y1+dy
      z2=z1+dz

      return
      end
************************************************************************
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
      subroutine todmsp(val,id,im,s,isign)
 
*** convert position radians to deg,min,sec
*** range is [-pi to +pi]
 
      implicit double precision(a-h,o-z)
      common/stugg/pi,rad,sol,f1,f2,w1,w2,we
 
    1 if(val.gt.pi) then
        val=val-pi-pi
        go to 1
      endif
 
    2 if(val.lt.-pi) then
        val=val+pi+pi
        go to 2
      endif
 
      if(val.lt.0.d0) then
        isign=-1
      else
        isign=+1
      endif
 
      s=dabs(val*rad)
      id=idint(s)
      s=(s-id)*60.d0
      im=idint(s)
      s=(s-im)*60.d0
 
*** account for rounding error
 
      is=idnint(s*1.d5)
      if(is.ge.6000000) then
        s=0.d0
        im=im+1
      endif
      if(im.ge.60) then
        im=0
        id=id+1
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
