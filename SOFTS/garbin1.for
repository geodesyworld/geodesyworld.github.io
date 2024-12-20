      program garbin1

*** test version of garbin for craig bailey
*** this version retains the garmin protocol wrappers and stuffing

      USE DFLIB
      USE DFPORT

      implicit double precision(a-h,o-z)
      integer*2      lenbuf
      character*7    abuf
      CHARACTER*1024 sbuf,rbuf
*** lines added to hold message for asynk2()
      integer*1 b1,b2
      parameter(b1=z'FF',b2=z'FF')
*** lines added to hold message for asynk2()

      lout=1
      open(lout,file='garbin1.bin',form='formatted')

*** argument for log time

      narg=NARGS()-1
      if(narg.eq.1) then
        call GETARG(1,abuf,lenbuf)
        if(lenbuf.lt.1) stop 20002
        read(abuf,'(i7)') logsec
      else
        write(*,*) 'not single argument'
        logsec=5
      endif
      seclog=dble(logsec)

*** set up the port parameters after connecting to the port and 
*** leave them unchanged until after the port has been released

*** port is com1, 19200 baud, 8-N-1 (no parity = 0, 1 stop = 0)

      iport=1
      ibaud=19200
      ipar=0
      idat=8
      istop=0

      if(SPORT_CONNECT(iport,0)                      .ne.0) stop 10001
      if(SPORT_CANCEL_IO(iport)                      .ne.0) stop 10002
      if(SPORT_SET_STATE(iport,ibaud,ipar,idat,istop).ne.0) stop 10003
      if(SPORT_SET_TIMEOUTS(iport,100,0,1000).ne.0) stop 10005
***********************************************************************
*** the following is the message to begin asynchronous comms        ***
***********************************************************************
      call asynk2(b1,b2,sbuf,lbuf)
***********************************************************************
      if(SPORT_WRITE_DATA(iport,sbuf,lbuf).ne.0) stop 10006

*** start the timer -- log data until time elapses

      elapsed=TIMEF( )

      elapsed=0.d0
  100 if(elapsed.lt.seclog) then
        if(SPORT_READ_DATA(iport,rbuf,lbuf).ne.0) stop 10007
        do 90 i=1,lbuf
   90   write(lout,'(a$)') rbuf(i:i)
        elapsed=TIMEF( )
        write(*,'(a,f9.1,a$)') 't=',elapsed,char(13)
        go to 100
      endif

      if(SPORT_RELEASE(iport).ne.0) stop 10008
      close(lout)

      end
      subroutine asynk2(b1,b2,sbuf,lbuf)

*** build 2 byte message (b1,b2) in type x1C packet
*** note: does not implement DLE stuffing (c.f. checksum)

      integer*1 b1,b2,cksum
      character*(*) sbuf

      integer*1 EOD,DLE,ETX,ACK,NAK
      parameter(DLE = z'10')
      parameter(ETX = z'03')
      parameter(ACK = z'06')
      parameter(NAK = z'15')
      parameter(EOD = z'0C')       !*** End of Data (in request commands)

      lbuf=8                       !*** data length = 2 bytes

      sbuf(1:1)=char(DLE)
      sbuf(2:2)=char(z'1C')        !*** message id -- asynchronous comms
      sbuf(3:3)=char(z'02')        !*** data length = 2 bytes

      sbuf(4:4)=char(b1)
      sbuf(5:5)=char(b2)

      sbuf(6:6)=char(cksum(sbuf,lbuf))
      sbuf(7:7)=char(DLE)
      sbuf(8:8)=char(ETX)

      return
      end
      integer*1 function cksum(sbuf,lbuf)

*** checksum on Garmin protocol buffer

      character*(*) sbuf
      integer*1 chksum

      chksum = 0
      do 100 i=2,lbuf-3
  100 chksum=chksum+ichar(sbuf(i:i))
      cksum=not(chksum)+1

      return
      end
