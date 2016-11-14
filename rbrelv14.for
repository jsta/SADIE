C     -------------------------------------------------------------------
C
C                                RBREL14
C
C              Spatial Analysis by Distance IndicEs (for counts)
C
C                             Version 1.4
C
C                         by Kelvin F. Conrad
C
C                          created 23 Oct 2007
C
C
C                      Original Version (v 1.3) by
C                           by Joe N. Perry
C
C                             2 June 1999
C
C
C                 Additions/Changes logged at end of file
C
C     -------------------------------------------------------------------

C     This code and the associated software program is the copyright (C) 2001
C     of J.N.Perry & his employers IACR Rothamsted Experimental Station AL5 2JQ UK

C     This program is free software; you can redistribute it and/or modify
C     it under the terms of the GNU General Public License as published by
C     the Free Software Foundation; either version 2 of the License, or
C     (at your option) any later version.

C     This program is distributed in the hope that it will be useful,
C     but WITHOUT ANY WARRANTY; without even the implied warranty of
C     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C     GNU General Public License for more details.

C     You should have received a copy of the GNU General Public License
C     along with this program reproduced below); if not, write to the Free Software
C     Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

C     A copy of this license appears at the end of this program code

C     You may contact Joe Perry, the author, as follows:
C     Professor J.N. Perry  DSc
C     Plant & Invertebrate Ecology Division
C     Rothamsted Experimental Station
C     Harpenden
C     Herts AL5 2JQ
C     UK
C     email:  joe.perry@bbsrc.ac.uk
C     Fax: +44 1582 760981
C     Phone:  +44 1582 763133 (extension 2375)

C     -------------------------------------------------------------------

C     Input/output details

C     Channels used:
C          5 - input for data in x,y pairs
C          6 - detailed output from program
C          7 - minimal output from program
C          8 - input for further data relating to
C                 random number seed, no.
C                 of simulations, etc.
C          9 - output for graphics only
C          10 - abbreviated output of indices and probabilities
C          11 - file for reading into Surfer, containing
C                 ordered cluster data

C     Schedule for program:

C                                                            IFLAG2 = 0
C     Read in data
C
C     Calculate distance to regularity of actual data
C                                                            IFLAG2 = 1
C     Calculate distance to regularity for permuted counts
C        and indices of clustering for observed data
C        using red-blue techniques.
C                                                           IFLAG  = 11
C     Go back for another NSIMS loops around to get test
C        statistics for clustering indices
C                                                     IFLAG2 = 12,13,14
C     Use ICSORT with minor and temporary changes to
C         IFLAG2,to sort test results for clustering
C         indices.
C                                                            IFLAG2 = 2
C     Calculate distance to crowding for actual data

C     The program differs from both nsadiec.for and xnsadiec.for in that it
C     calculates average distances of flow for each sample unit, for use in
C     the red-blue plots and the associated EDFs of those averages.  Clustering
C     is defined by a doubly standardized index, allowing for observed count
C     and unit, thus allowing for numerical effects and edge-effects.
C     Unlike earlier versions of the red-blue software, this version does
C     not calculate the upper and lower confidence envelopes for the EDF
C     plots of the average flow distances - these are now considered redundant.
C     -------------------------------------------------------------------

      PROGRAM RB14KFC

C     for command-line arguments
      USE PORTLIB  !comment out for g77 compiler
 
      IMPLICIT NONE

      integer*4 NOP,NSIMS,ICOUNT(2000),IFLAG1,IFLAG2,KOST,I,J,II,JJ,K
      integer*4 IRANKD(2000),KN,KZ,KP,IMNODE,INCAT(2000),WC(2000)
      integer*4 INFROM(2600000),INTO(2600000),INCOST(2600000)
      integer*4 ITOTSY,ISUPLY(2000),NUMNOD,IDUMMYU,INMARC
      integer*4 ITOT,RN(2000),NP
      integer*4 JRANKD(2000),KIAFLN,KIAFLP,IOUTER,IINNER,KIA,KIAFLC
      integer*4 IDUMMYM,IDUMMYL,IX,IY,IZ,ISEED,K5PSIM,NN,D2,D1,IDUMMY
      integer*4 RC(2000),TEMPWC,RUNMIN,RUNMAX,count,status
      integer*4 errorflag ! transfres error codes among sbr's

      INTEGER(2) cmdno	!for command-line arguments

      double precision X(2000),Y(2000),DIST(2600000),DTREGS(6000)
      double precision DTCROS(6000),DTOCRO,WX(2000),WY(2000)
      double precision EDFO(2000),EDFS(2000),TEDFS(2000)
      double precision ASFC(2000),DASFC(2000),ASFP(2000),DASFP(2000)
      double precision TEMPAP(2000),SORTKIA(78000),IAFSP(6000)
      double precision IAFSN(6000),IAFSA(6000),IAFDUM(2000),INK5PU(153)
      double precision IPK5PU(153),IAFMAXN(6000),IAFMAXP(6000)
      double precision INK5PM(153),IPK5PL(153),IAF(2000)
      double precision IPK5PM(153),AIPK5PM,INK5PL(153)
      double precision DUMMYX(2000),DUMMYY(2000),SORTNOP(2000)
      double precision DISSCL,MINDIS,MAXDIS,DTOREG,DTREGO,DTCROO
      double precision POINCX,POINCY,COUNCX,COUNCY,DPCCC,MAXD
      double precision DASCP,ASCP,IAFOPS,IAFONS,AINK5PM,IAFOAS
      double precision AINK5PU,AIPK5PU,IAFMXP5U,IAFMXN5U
      double precision DUM,DUMD,RATIO1,DUMMY1,DUMMY2,AINK5PL,AIPK5PL
      double precision TEMPX,TEMPY

      Logical RFLAG,intflag

      CHARACTER(80) buf       ! for command-line arguments
      CHARACTER(80) messerror ! for error messages
      character(256) infile5  ! filename for input (usu rbni5.dat)
	character(256) rbni8    ! filename for control parameters
      character(256) projpath ! project path (usu current folder)
	character(256) rbno6    ! output file pathname for rbno6
	character(256) rbno7    ! output file pathname for rbno7
	character(256) rbno9    ! output file pathname for rbno9
	character(256) rbno10   ! output file pathname for rbno10
	character(256) cluster  ! output file pathname for cluster.dat
      
      COMMON /ORIG/ X,Y,WX,WY,ICOUNT,WC,IRANKD,NOP
      COMMON /LINK/IMNODE,INCAT,ITOTSY,ISUPLY,INFROM,INTO,INCOST,INMARC
      common /errs/ messerror, errorflag ! transfres error codes among sbr's

	projpath = ' '
      intflag = .false.      ! don't use interactive mode
      infile5 = 'rbni5.dat'  ! assume local rbni5 for input
	i = getcwd(projpath)   ! projpath defaults to current directory


C     Process command line input
C     __________________________
 
      count=IARGC()
      DO cmdno = 1, count
        CALL GETARG(cmdno, buf, status)
        IF (buf(2:status) .EQ.'r') THEN
          RFLAG=.TRUE.
        else if (buf(2:status).EQ.'h') THEN
            call displayhelp()
            stop ' '
        elseif (buf(2:).eq.'i') then
          intflag=.true.
        ELSE
          write(*,*)  ' '
          WRITE(*,*)  '    Invalid command-line switch'
          write(*,*)  ' '
          call displayhelp()
          stop ' '
        END IF
      END DO

C     g77 version -- comment out unneed version
C      __________________________

C      count=IARGC()
C      DO cmdno = 1, count
C        CALL GETARG(cmdno, buf)
C        IF (buf(2:) .EQ.'r') THEN
C          RFLAG=.TRUE.
C        else if (buf(2:).EQ.'h') THEN
C          call displayhelp()
C          stop ' '
C        elseif (buf(2:).eq.'i') then
C          intflag=.true.
C        ELSE
C          write(*,*)  ' '
C          WRITE(*,*)  '    Invalid command-line switch'
C          write(*,*)  ' '
C          call displayhelp()
C          stop
C        END IF
C      END DO
C     -------- end command-line code -----------

C     Set errorflag (0 = no errors)
C     _______________________________
      errorflag = 0
      
C     If interactive mode, collect user input
C     _______________________________________________
      if (intflag) then
	  Call Interactor(projpath,infile5,rflag)
      endif

C     Open unit 12 to log data read and errors
C     _______________________________________________
      rbno6 = projpath(:LNBLNK(projpath)) // '\' // 'rbno6.dat'
      OPEN (UNIT=12,FILE=rbno6, err=25, STATUS='REPLACE')

C     Open the input file (usu. rbni5.dat)
C     _______________________________

      OPEN (UNIT=15, FILE=infile5, err=26)

      goto 27            !skip error-handling
      
   25 errorflag = 1      !rbno6 I/O error
      messerror = '  Error opening/creating RBNO6.DAT'
      goto 6000          !end the program
      
   26 errorflag = 2      !input file I/O error
      messerror = '  Error opening input files'
      goto 6000          !end the program
      
   27 continue           !no input file I/O errors

C     Read in actual data
C     ___________________

C     used for debugging of npar method
C     RFLAG=.TRUE.

      IFLAG2=0          !Mode = Read in data

C     Read in data and store number of values
C         in NOP; x,y pairs - one per line
      DO I = 1, 2001
         READ (15,*,END=55, err=50) X(I), Y(I), ICOUNT(I)
         IF (ICOUNT(I).LT.0) THEN
             messerror = '  Negative numbers not allowed in data'
             errorflag = 3   ! error: negative counts in the data
             goto 6000       !end the program
         ENDIF
      end do

C     Too much data for this version (NOP>2000) so report and stop
      messerror = '  More than 2000 data values '//
     +            '- Too many sampling points'
      errorflag = 4      ! Too many sampling points
      goto 6000          ! end the program
      
  50  errorflag = 5      ! error reading unit 5
      messerror =  '  Error reading input file!'
      goto 6000          !end the program

  55  continue           ! no inpput errors have occurred

      NOP = I - 1        ! record the number of rows entered

C     Check to ensure there are some data values read in
      IF (NOP .LE. 0) THEN
         errorflag = 6      !  No numbers in data
         messerror = '  No numbers in data'
         goto 6000          ! end the program
      ENDIF

      DO I = 1,NOP
C        create working copies of the input data
         WC(I) = ICOUNT(I)
         WX(I) = X(I)
         WY(I) = Y(I)
C        initialize some values concerning the red-blue EDF plots
         EDFO(I) = 0.0D0
         EDFS(I) = 0.0D0
      end do

C     Write out actual data
C     _______________________________

      write (12,*)
	write (12,109) infile5
 109  format (1x,'Data read from: ', A)
      write (12,*)

      WRITE (12,*) 'Original data'
      WRITE (12,*)
      WRITE (12,*) 'There are ',NOP,' counts, which are:'
      WRITE (12,*) '           X                   Y         ',
     +            '              Count      Ref. No.'
      DO I = 1,NOP
         WRITE (12,110) X(I), Y(I), ICOUNT(I), I
      end do

110   FORMAT (1X ,F14.5,6X,F14.5,20X,I7,5X,I7)

      WRITE (12,*)

C     Check flag to see if Non-parametric version is required (RFLAG=TRUE)
      If (RFLAG) then
C     RFLAG=TRUE is so begin:
C     Reorder original data in ascending order according to algorithm of King
C     Rank counts, WC() and order WX(), WY() according to this ranking.
        DO K = 2, NOP
           I = K
   44       IF (I .GT. 1) THEN        
                IF (WC(I-1) .GT. WC(I)) THEN
                   TEMPWC = WC(I-1)
                   TEMPX = WX(I-1)
                   TEMPY = WY(I-1)
                   WC(I-1) = WC(I)
                   WX(I-1) = WX(I)
                   WY(I-1) = WY(I)
                   WC(I)  = TEMPWC
                   WX(I)  = TEMPX
                   WY(I)  = TEMPY
                   I = I - 1
                ELSE
                  I = 1
                ENDIF
                GOTO 44
             ENDIF 
        end do

C       Run through list finding runs of same numbers.  Assign (rankingsx2),
C       instead of counts, from 1 (lowest) to NOP (highest)
C       First some initialisation
        RUNMIN=1
        RUNMAX=2
        
C       Test for a run
   66   IF (WC(RUNMIN).EQ.WC(RUNMAX)) THEN
          RUNMAX=RUNMAX+1
          IF (RUNMAX.LE.NOP) THEN
            GOTO 66
          ENDIF
        ENDIF

        RUNMAX=RUNMAX-1

        DO J = RUNMIN,RUNMAX
          RC(J)=RUNMAX+RUNMIN
        end do

        RUNMIN=RUNMAX+1
        RUNMAX=RUNMIN+1

        IF (RUNMAX.LE.NOP) THEN
          GOTO 66
        ENDIF

        IF (RUNMIN.EQ.NOP) THEN
          RC(J)=2*NOP
        ENDIF

C       Write out new non-parametric and reordered original data, side by side
        WRITE (12,*)
        WRITE (12,*)
        WRITE (12,*) 'New data'
        WRITE (12,*)
        WRITE (12,*) 'There are ',NOP,' counts, which are:'
        WRITE (12,*) '         X               Y       ',
     +            ' Ref. No.   Orig. Count     New ranked count'
11000   FORMAT (1X ,F14.5,2X,F14.5,3X,I7,3X,I7,11X,I7)

        DO I = 1,NOP
          WRITE (12,11000) WX(I), WY(I), I, WC(I), RC(I)
        end do

        WRITE (12,*)
        WRITE (12,*)
C       Now overwrite original data completely
        DO I = 1,NOP
          WC(I) = RC(I)
          ICOUNT(I) = RC(I)
          X(I) = WX(I)
          Y(I) = WY(I)
        end do
       endif   !if RFLAG

C     calculate the units centroid, POINTC, and the centroid of these
C     original counts, COUNTC, and the distance between them.
      POINCX=0.0D0
      POINCY=0.0D0
      COUNCX=0.0D0
      COUNCY=0.0D0
      DUM = 0.0D0

      DO I = 1,NOP
         POINCX = POINCX + X(I)
         POINCY = POINCY + Y(I)
         COUNCX = COUNCX + (ICOUNT(I)*X(I))
         COUNCY = COUNCY + (ICOUNT(I)*Y(I))
         DUM = DUM + DFLOAT(ICOUNT(I))
      end do

      POINCX = POINCX/DFLOAT(NOP)
      POINCY = POINCY/DFLOAT(NOP)
      COUNCX = COUNCX/DUM
      COUNCY = COUNCY/DUM

C     Now calculate the distance between the units centroid and the counts centroid,
C     denoted dpccc
      DPCCC = ((POINCX-COUNCX)*(POINCX-COUNCX))+
     +        ((POINCY-COUNCY)*(POINCY-COUNCY))
      DPCCC = SQRT(DPCCC)
      MAXD = 0.0D0
      DO I = 1,NOP
        DO J = 1,NOP
           DUMD = ((X(I)-X(J))*(X(I)-X(J)))+((Y(I)-Y(J))*(Y(I)-Y(J)))
           DUMD = DSQRT(DUMD)
           IF (DUMD.GT.MAXD) THEN
              MAXD = DUMD
           ENDIF
         end do
      end do

      if (intflag) then
        print *, ' '
   90   print *, 'Enter a seed for the random number generator '//
     +           '(1 < seed < 30000)'
	    WRITE ( *, 91 )
   91   FORMAT (' --> ', $ )
        read (*,*) iseed
        IF (ISEED.GT.30000) THEN
          print *, 'Value of ISEED too large - maximum is 30000'
          goto 90
        endif
        print *, ' '
   95   print *, 'Enter a number for k5psim (max=153) '
	    WRITE ( *, 91 )
        read (*,*) k5psim
        IF (k5psim.GT.153) THEN
          print *, 'Too many simulations - max. k5psim is 153'
          goto 95
        endif
        NSIMS = 39*K5PSIM
        print *, ' '
        print *, 'Program running ...'
        print *, ' '
      
      else     ! not in interactive mode
      
C       Open rbni8.dat for input or kp5sims and iseed
        rbni8 = projpath(:LNBLNK(projpath)) // '\' // 'rbni8.dat'
        OPEN (UNIT=8, FILE=rbni8, err=96)

C       Read a seed for random number generation from a different channel
C        from that in which data was read; integer. If no data, assume 100.
C        Value can be between 1 and 30000
        READ (8,*,END=5,err=100) ISEED
        IF (ISEED.GT.30000) THEN
          messerror = '  Value of ISEED too large - maximum is 30000'
          errorflag = 8    !  Value of ISEED too large
          goto 6000        ! end program
        ENDIF

C       Read the value of K5PSIM, from which shall immediately be calculated NSIMS,
C       the number of simulations of a random pattern to be run; integer.
C       Note that NSIMS is exactly 39 times K5PSIM.
C       If no data, assume K5PSIM=25,NSIMS=975.  Array DTREGS() in main program assumes maximum
C       value of 6000 for NSIMS; this equates to a maximum value of 153 for K5PSIM.
        READ (8,*,END=6, err=100) K5PSIM
        NSIMS = 39*K5PSIM
        IF (K5PSIM.GT.153) THEN
          messerror = '  Too many simulations for array bounds'//
     +                ' - maximum is 6000'
          errorflag = 9     ! Too many simulations
          goto 6000
        ENDIF

        goto 101            ! Skip error handling

   96   errorflag = 2      !input file I/O error
        messerror = '  Error opening input files'
        goto 6000          !end the program

  100   errorflag = 10      ! Error reading unit 8
        write(12,*)
        write(12,*) '  Error reading unit 8'
        goto 6000           ! end the program

  101   continue            ! no read errors

        GOTO 105            ! Skip iseed/kp5sim errors
    5   ISEED = 100
    6   K5PSIM = 25
        NSIMS = 39*K5PSIM
  105   continue
      endif

      WRITE(12,*) 'ISEED = ', ISEED
      WRITE(12,*) 'K5PSIM = ', K5PSIM
      WRITE(12,*) 'NSIMS = ', NSIMS

C     Write descriptive statistics to log file
      WRITE(12,*)
      WRITE(12,*)

 9999 FORMAT(1X ,'x co-ordinate, centroid of sample units (location P)',
     +': ', F12.6)
 9998 FORMAT(1X ,'y co-ordinate, centroid of sample units (location P)',
     +': ',F12.6)
 9997 FORMAT(1X ,'x co-ordinate, centroid of counts       (location C)',
     +': ',F12.6)
 9996 FORMAT(1X ,'y co-ordinate, centroid of counts       (location C)',
     +': ',F12.6)
 9995 FORMAT(1X ,'Distance between centroids of units and counts',
     +', (delta): ', F12.6)
      WRITE(12,9999) POINCX
      WRITE(12,9998) POINCY
      WRITE (12,*)
      WRITE(12,9997) COUNCX
      WRITE(12,9996) COUNCY
      WRITE (12,*)
      WRITE(12,9995) DPCCC
      WRITE(12,*)
      WRITE(12,*) 'The maximum distance between sample units is: ',MAXD
      WRITE(12,*)
      WRITE(12,*)

C     Open the output files
C     _______________________________

      rbno7 = projpath(:LNBLNK(projpath)) // '\' // 'rbno7.dat'
      rbno9 = projpath(:LNBLNK(projpath)) // '\' // 'rbno9.dat'
      rbno10 = projpath(:LNBLNK(projpath)) // '\' // 'rbno10.dat'
      cluster = projpath(:LNBLNK(projpath)) // '\' // 'cluster.dat'


      OPEN (UNIT=7,FILE=rbno7,STATUS='NEW',err=75)
      OPEN (UNIT=9,FILE=rbno9,STATUS='NEW',err=75)
      OPEN (UNIT=10,FILE=rbno10,STATUS='NEW',err=75)
      OPEN (UNIT=11,FILE=cluster, STATUS='NEW',err=75)

      goto 76            !skip output file error-handling

   75 errorflag = 7      ! error creating output files
      messerror = '  Error creating output files'
      goto 6000          ! end the program

   76 continue           ! no errors creating output files

C     Now calculate distance to regularity for actual data
C     ____________________________________________________

C     Scale counts by inflation of integers in working data WC(), to ensure that
C       mean of working counts is an integer;
C       only if IFLAG1 is returned from CSCALR as zero was no scaling done.
      IFLAG1=0
      
      CALL CSCALR(IFLAG1,ITOT,IFLAG2)
      
      WRITE(12,*)
      WRITE(12,*)
      WRITE(12,*)'Total counts = ',ITOT
      WRITE(7,*)
      WRITE(7,*)
      WRITE (9,*) 'I.A.F. plot information for graphics'
      WRITE (9,*) 'There follows in column order'
      WRITE (9,*) '(but in no particular unit order):'
      WRITE (9,*)'1: x-coordinate of FROM unit'
      WRITE (9,*)'2: y-coordinate of FROM unit'
      WRITE (9,*)'3: x-coordinate of TO unit'
      WRITE (9,*)'4: y-coordinate of TO unit'
      WRITE (9,*)'5: flow'
      WRITE (9,*)
      
C     Rank the counts in WC(); rank WX() and WY() in parallel;
C     count the number postive, zero and negative;
C     identify from and to arrays, get INCAT(I);
C     create arrays IRANKD() and JRANKD()
      KN = 0
      KZ = 0
      KP = 0
      MINDIS = 0.0D0
      MAXDIS = 0.0D0
      NUMNOD = 0

      CALL RANKCO(KN,KZ,KP,DIST,MINDIS,MAXDIS,JRANKD,IFLAG2)

C     if an error has occurred in RANKO, exit the program
      if ((errorflag .ge. 11) .and. (errorflag .le. 15)) goto 6000

C     Get and order distances in working data, then scale
      CALL DSCALR(DIST,DISSCL,MINDIS,MAXDIS,NOP)

C     If IFLAG2 is zero then some output will be printed by TRANSP -
C     otherwise (if using working on simulated not actual data) not.
      CALL TRANSP(IFLAG2,KOST,IRANKD,X,Y,NUMNOD)
      
C     if an error has occurred in TRANSP, exit the program
      if ((errorflag .ge. 16) .and. (errorflag .le. 23)) goto 6000

      CALL AVFLPP(IFLAG2,EDFO,NOP,X,Y,RN,ICOUNT,DASFC,NUMNOD)
      
C     On return from AVFLPP, EDFO() contains	ranked values of average
C     flow distances for original data; DASFC() is an unused dummy variable;
C     RN() has original reference numbers, reordered according to EDFO().
C     Calculate NP, the number of (positive) outflows and NN the number of (negative) inflows
      NP=0
      NN=0

      DO J = 1,NOP
        if (EDFO(J).GE.0.0D0) then
          NP = NP + 1
        else
          NP = NP + 0  ! not strictly necessary
        endif
        if (EDFO(J).LT.0.0D0) then
          NN = NN + 1
        else
          NN = NN + 0  ! not strictly necessary
        endif
      end do

C     Calculate KIAFLN as 39*NN, KIAFLP as 39*NP, and KIA as 39*NOP
      KIA = 39*NOP
      KIAFLN = 39*NN
      KIAFLP = 39*NP

      CALL UNSCAL(KOST,ITOTSY,MAXDIS,DTREGO,IFLAG1,NOP,IFLAG2)

      WRITE(9,*)
      WRITE(9,*)
      WRITE(9,*) 'The distance to regularity of the actual data is :'
      WRITE(9,1100) DTREGO
1100  format (1X , F21.10)
      WRITE(9,*)
      WRITE(9,*)

      IFLAG2=1

      WRITE(9,*) 'Permuted data'
      WRITE(9,*)
      WRITE(9,*) 'Distances to regularity of the permuted data are :'

C     Now start simulation of the NSIMS sets of permuted counts for
C      distance to regularity
C     ______________________________________________________________

      WRITE(7,1004)
      WRITE(7,1005)
 1004 FORMAT('Results for regularity - permutations ',
     +'conditioning on observed counts')
 1005 FORMAT('**************************************',
     +'*******************************')

C     Initialize values for random number generation
      IZ = ISEED
      IX = 30001 - IZ
      IY = 15000
      
C     Initialize values for EDF calculations for red-blue plots;

C     TEDFS() is the total of the values of ordered average flow distances over all simulations
C     ASFC() is the unranked total of the average flow distances (for that observed count)
C            over all simulations
C     ASFC() is the unranked total of the average flow distances (for that observed count)
C            over all simulations
C     ASFP() is the unranked total of the average flow distances (for that observed unit)
C            over all simulations
C     DASCP will be the average value of the average flow distance over all observations and
C           over all simulations for all flows (similarly ASCP)
C           (cf TEDFS(), but TEDFS() is ranked and ASF values are all unranked)

      DO I = 1,NOP
        TEDFS(I) = 0.0D0
        ASFC(I)   = 0.0D0
        ASFP(I)   = 0.0D0
        ASCP   = 0.0D0
      end do

C     Loop round to get NSIMS new values of DTOREG
      DO I = 1,NSIMS
         DASCP = 0.0D0
         
         CALL PERMUT(IX,IY,IZ)
         CALL CSCALR(IFLAG1,ITOT,IFLAG2)
         CALL RANKCO(KN,KZ,KP,DIST,MINDIS,MAXDIS,JRANKD,IFLAG2)

C        if an error has occurred in RANKO, exit the program
         if ((errorflag .ge. 11) .and. (errorflag .le. 15)) goto 6000
      
         CALL DSCALR(DIST,DISSCL,MINDIS,MAXDIS,NOP)
         CALL TRANSP(IFLAG2,KOST,IRANKD,WX,WY,NUMNOD)

C        if an error has occurred in TRANSP, exit the program
         if ((errorflag .ge. 16) .and. (errorflag .le. 23)) goto 6000

         CALL AVFLPP(IFLAG2,EDFS,NOP,WX,WY,RN,ICOUNT,DASFC,NUMNOD)
         
C        On return from AVFLPP, EDFS() contains	ranked values of average
C        flow distances for simulated data; DASFC() contains unranked values
C        of average flow distances for this set of simulated data, with the tie to
C        the relevant value of WC() retained;
C        RN() has not been altered by the subroutine, since this is not original data.

C       Now, the values of DASFP() are found, which are unranked average flow distances
C        for this set of simulated data, with a tie to the observed unit.
C        Also calculated are: DASCP, which is the overall mean average flow distance
C        for all flows for this set of simulated data.

         DO II = 1,NOP
           DO JJ = 1,NOP
            IF ((WX(JJ).EQ.X(II)).AND.(WY(JJ).EQ.Y(II))) THEN
              DASFP(II) = DABS(DASFC(JJ))
              DASCP = DASCP + (DABS(DASFC(JJ))/DFLOAT(NOP))
            ENDIF
   	       end do
         end do
         
C        End of calculation of DASFP() and DASCP

         CALL UNSCAL(KOST,ITOTSY,MAXDIS,DTOREG,IFLAG1,NOP,IFLAG2)
         
         WRITE(9,1111) DTOREG
 1111    format (1X , F21.10)

C        Store the value of distance to regularity for this Ith simulation in DTREGS(I).
         DTREGS(I) = DTOREG
C        Loop over all sample units
         DO J = 1,NOP
C        Update total of EDFS() for this simulation - later the total will be
C        replaced by the mean
            TEDFS(J)=TEDFS(J)+EDFS(J)
C           Update totals of ASFC, ASFP for this simulation - later these totals will be
C           replaced by the mean, then used in subroutine OUTPU4 to construct
C           IAF(), the index of clustering for each sample unit
            ASFC(J) = ASFC(J) + DASFC(J)
C           and assign the correct sign to ASFP()
            ASFP(J) = ASFP(J) + DASFP(J)
         end do
C        and update the total of ASCP for this simulation -
         ASCP = ASCP + DASCP
C        Later these totals will be replaced by means, then used to
C        construct IAF(), the index of clustering for each sample unit
      end do

      CALL DTSORT(IFLAG2,DTREGS,DTREGO,NSIMS)
      
C     Now calculate various quantities concerning the distribution of average distances
C     of flow, for the simulated and then observed data.
C     First, replace TEDFS(), a total, by the mean, EDFS()
C     Note that, regardless of its previous meaning, EDFS() is now the mean value of
C     ordered average flow distances over all simulations.
C     Similarly replace ASFC(), ASFP() and ASCP totals by means.
      DO J = 1,NOP
         EDFS(J)=TEDFS(J)/DFLOAT(NSIMS)
         ASFC(J)=ASFC(J)/DFLOAT(NSIMS)
         ASFP(J)=ASFP(J)/DFLOAT(NSIMS)
      end do
      
      ASCP=ASCP/DFLOAT(NSIMS)

C     Output IAF(), scaled unitless index of clustering, formed by division of
C     (a.f.d. in EDFO()  x  a.f.d. all simulations) by
C     (a.f.d. simulations with that count  x  a.f.d. simulations at that unit)
C     Input are:
C       EDFO()  (ordered average flow distances for original data)
C       EDFS()  (mean value of ordered average flow distances over all sim./cons.)
C       AFC(I)  (mean value of the average flow distances over all sim./cons.,
C                tied to the value of the Ith count in WC(), i.e. ICOUNT after ranking
C       AFP(I)  (mean value of the average flow distances over all sim./cons.,
C                tied to the value of the Ith unit, (X(I),Y(I))
C       AFCP    (mean value of the average flow distances over all sim./cons.)
C      The reason that AFC() is tied to the count at the observed unit is that although
C      PERMUT changes only WC(), RANKCO ranks WC(), so the WC() values are always
C      the same values output in the same (ranked) order.  Since RANKCO also orders
C      WX() and WY() in parallel to WC(), it is WX() and WY() that change in this
C      subroutine, and not the counts.  ICOUNT is not actually used, aside from
C      printing, but the fact that AFC() is tied to counts is of paramount importance
C      regarding scaling.
C      JRANKD() (array relating to rank of each of the original counts)
C      Recall that if JRANKD(1)=4 this means that the count in reference number 1 was
C      ranked number 4 (i.e. was the fourth smallest), and if JRANKD(K)=NOP then the
C      largest count was in reference number K.
C      Outside this subroutine, back in the main program, another variable, AFP()
C      was derived, tied to the observed unit rather than the observed count.
C      Also, for output is some limited information concerning the EDFs of simulated data,
C      but since these are now thought redundant, the 95% CIs are no longer calculated
C      Other variables - unchanged from previous definitions
C      NP and NN absolutely MUST be initialised

      NP=0
      NN=0

      WRITE(7,*)
      WRITE(7,*)
 1000 FORMAT('Results for clustering')
 1001 FORMAT('**********************')
      WRITE(7,1000)
      WRITE(7,1001)
      WRITE(7,*)
      WRITE (7,*) 'There follows in column order, '
      WRITE (7,*)'with minus signs denoting direction of flow ',
     +'where appropriate,'
      WRITE(7,*)'where average flow distance is abbreviated to a.f.d.,'
      WRITE(7,*)'upper-case greek gamma subscript i is denoted as: Y_i,'
      WRITE(7,*)'and lower-case greek nu subscript i is denoted as: v_i'
      WRITE(7,*)
      WRITE(7,*)'1: Observed a.f.d. (after ranking):         Y_i or Y_j'
      WRITE(7,*)'2: Overall average of permuted a.f.d.s:     o_Y'
      WRITE(7,*)'3: Average permuted a.f.d. for this count:  c_Y'
      WRITE(7,*)'4: Average permuted a.f.d. for this unit:   i_Y or j_Y'
      WRITE(7,*)'5: Standardized index of clustering:        v_i or v_j'
      WRITE (7,*)'6: x-coordinate of unit'
      WRITE (7,*)'7: y-coordinate of unit'
      WRITE (7,*)'8: Rank of ordered distance in col. 1'
      WRITE (7,*)'9: Count of unit'
      WRITE (7,*)'10: Rank of unit'
      WRITE (7,*)'11: Orig. unit no. of data as read'
      WRITE(7,*)
      WRITE(7,*)'The most important columns are 5 (clustering index),',
     +' 6 (x) and 7 (y)'
      WRITE(7,*)
      WRITE(7,*)'                                            Cluster',
     +'ing'
      WRITE(7,*)'                                               Index',
     +'      x          y'
8888  FORMAT(1X ,4(E10.4,1X),2X,F7.3,1X,F10.4,1X,F10.4,1X,
     +           I4,1X,I4,1X,I4,1X,I4)
      WRITE (9,*)
      WRITE (9,*)
      WRITE (9,*) 'There follows in column order:'
      WRITE (9,*)'1: x-coordinate of unit'
      WRITE (9,*)'2: y-coordinate of unit'
      WRITE (9,*)'3: Index of clustering'
8886  FORMAT(1X ,F18.8,1X,F18.8,1X,F10.6)
      WRITE (9,*)

C     There follows what is at the heart of the red-blue program, the construction of the
C     clustering index, IAF().  The idea is that IAF() be scaled to allow both for the
C     observed count and the observed unit, because large counts tend to have much greater
C     average flow distances, and edge-units, or units relatively far from the other units
C     of a set also have relatively large a.f.d.s.  Hence, if the average over all
C     simulations is denoted AFCP; if the mean for that count is denoted AFC(); and for that
C     unit is denoted AFP(); then the index is constructed as:

C                IAF() = EDFO()   x     AFCP      x         AFCP
C                       ------         -------             -----
C                        AFCP           AFC()               AFP()
c
c
C                           = EDFO() x  AFCP
C                           ------    -------
C                            AFC()     AFP()

C     Also, since EDFO(), ordered according to X(),Y(), has been reordered according to RN(),
C     then so must AFP(), currently ordered w.r.t. X(),Y() be reordered by RN().  Also, AFC()
C     must be ordered not only by RN(), but first to account for the ordering
C     done on ICOUNT() to get WC(); this reordering is achieved by JRANKD().

      DO I = 1,NOP
        D1 = RN(I)
        D2 = JRANKD(D1)
C       IAF is the ordered scaled average flow distance index of clustering
        IF (EDFO(I).NE.0.0D0) THEN
          IAF(I)=(DSIGN(1.0D0,EDFO(I)))*EDFO(I)*ASCP/(ASFC(D2)*ASFP(D1))
        ENDIF
        IF (EDFO(I).EQ.0.0D0) THEN
           IAF(I)=0.0D0
        ENDIF
        DUMMY1 = DSIGN(1.0D0,EDFO(I))*ASCP
        DUMMY2 = DSIGN(1.0D0,EDFO(I))*ASFP(D1)
        WRITE (7,8888) EDFO(I),DUMMY1,ASFC(D2),DUMMY2,IAF(I),
     +                 X(D1),Y(D1),I,ICOUNT(D1),D2,D1
        WRITE (9,8886) X(D1),Y(D1),IAF(I)
        DUMMYX(I) = X(D1)
        DUMMYY(I) = Y(D1)
      end do

C     IAFOP is the mean of the positive outflow values of the clustering index, IAF()
C     IAFON is the mean of the negative inflow values of the clustering index, IAF()
C     IAFOA is the mean over all units of the absolute values of the clustering index, IAF()
C     Calculate number of outflows and inflows and then mean IAF() over them separately

      IAFOPS = 0.0D0
      IAFONS = 0.0D0
      IAFOAS = 0.0D0

      DO I = 1,NOP
        if (IAF(I).GE.0.0D0) then
          NP = NP + 1
        else
          NP = NP + 0  ! not strictly necessary
        endif
        if (IAF(I).LT.0.0D0) then
          NN = NN + 1
        else
          NN = NN + 0
        endif
        if (IAF(I).GE.0.0D0) then
          IAFOPS = IAFOPS + IAF(I) 
C         IAFOPS = IAFOPS + (IAF(I)*(IAF(I).GE.0.0D0))
        endif
        if (IAF(I).LT.0.0D0) then
          IAFONS = IAFONS + IAF(I)
C         IAFONS = IAFONS + (IAF(I)*(IAF(I).LT.0.0D0))
        endif
        if (IAF(I).LT.0.0D0) then
          IAFOAS=IAFOAS-IAF(I)
C         IAFOAS=IAFOAS-(IAF(I)*(IAF(I).LT.0.0D0))+
C    1       (IAF(I)*(IAF(I).GE.0.0D0))
        endif
        if (IAF(I).GE.0.0D0) then
          IAFOAS=IAFOAS+IAF(I)
C         IAFOAS=IAFOAS-(IAF(I)*(IAF(I).LT.0.0D0))+
C    1       (IAF(I)*(IAF(I).GE.0.0D0))
        endif
      end do

C     Now order and write out clustering indices with additional
C     information, as explained in subroutine ORDERCL

      CALL ORDERCL(DUMMYX,DUMMYY,IAF,NN,NP,NOP)

      IAFOPS = IAFOPS/DFLOAT(NP)
      IAFONS = IAFONS/DFLOAT(NN)
      IAFOAS = IAFOAS/DFLOAT(NOP)

      WRITE(7,*)
      WRITE(7,*)
      WRITE(10,2000) IAFONS
      WRITE(10,2000) IAFOPS
 2000 FORMAT(1X,F7.3)
 
      WRITE(7,8887) IAFONS
 8887 FORMAT(1X,'Mean of clustering index, v_j, over inflows:  ',F7.3)
      WRITE(7,8885) IAFOPS
 8885 FORMAT(1X,'Mean of clustering index, v_i, over outflows: ',F7.3)
      WRITE(7,8884) IAFOAS
 8884 FORMAT(1X,'Mean of v_i and abs. |v_j|, over all flows: ',F7.3)

      WRITE(9,*)
      WRITE(9,*)
      WRITE(9,*)'Limited information for EDF plots, '
      WRITE(9,*)'these are now calculated without the CIs'
      WRITE (9,*) 'There follows in column order:'
      WRITE (9,*)'1: Ordered unscaled average flow distances - observed'
      WRITE (9,*)'2: Ordered scaled by average of EDF values'
      WRITE (9,*)'3: Rank of ordered distance in col. 1'
      WRITE (9,*)'4: Average of simulated ordered distances ',
     +            'with this rank (unscaled)'
      WRITE (9,*)'5: Orig. unit no. of ordered distance in col. 1'
      WRITE(9,*)

      DO I = 1,NOP
        D1 = RN(I)
        IF (EDFS(I).NE.0.0D0) THEN
          RATIO1 = DSIGN(1.0D0,EDFO(I))*EDFO(I)/EDFS(I)
        ENDIF
        IF (EDFS(I).EQ.0.0D0) THEN
          RATIO1 = 0.0D0
        ENDIF
        WRITE (9,7778) EDFO(I),RATIO1,I,EDFS(I),D1
      end do
 7778	FORMAT(1X,2(F15.7,1X),I4,1X,F15.7,1X,I4)

      WRITE(9,*)
      WRITE(9,*)
      WRITE(9,8887) IAFONS
      WRITE(9,8885) IAFOPS
      WRITE(9,8884) IAFOAS
      WRITE(9,*)
      WRITE(9,*)
      WRITE(9,*) 'Mean clustering indices of permuted data are: '
      WRITE(9,*) '(-ve inflows, +ve outflows, abs(flows))'
      write(9,*) ' '
      
C     At this point, IAFOPS is the mean of the positive outflow values
C     of the clustering index, IAF(), scaled by the standard permutations.
C     IAFONS is the similar mean of the negative inflow values.
C     IAFOAS is the similar mean of the absolute values of all flows.
C     This is the start of a new phase of the program, and IFLAG2 is changed so
C     that within AVFLPP, RN() can be altered, as it was for the original data
C     but not for the first batch of simulations.
C     The program now takes another NSIMS randomizations to get a distribution of
C     values of the clustering index under randomizations
C     Loop round to get NSIMS new values of IAFSP(), IAFSN() and IAFSA()
C     Initialise variables AINK5PU, AIPK5PU, AINK5PM, AIPK5PM, AINK5PL, AIPK5L

      IFLAG2=11
      AINK5PU = 0.0D0
      AIPK5PU = 0.0D0
      AINK5PM = 0.0D0
      AIPK5PM = 0.0D0
      AINK5PL = 0.0D0
      AIPK5PL = 0.0D0

      DO IOUTER = 1,K5PSIM
C     Initialise variable KIAFLC
        KIAFLC = 0
        DO IINNER = 1,39
          I = IINNER + ((IOUTER-1)*39)
          IAFSN(I) = 0.0D0
          IAFSP(I) = 0.0D0
          
          CALL PERMUT(IX,IY,IZ)
          CALL CSCALR(IFLAG1,ITOT,IFLAG2)
          CALL RANKCO(KN,KZ,KP,DIST,MINDIS,MAXDIS,JRANKD,IFLAG2)

C        if an error has occurred in RANKO, exit the program
         if ((errorflag .ge. 11) .and. (errorflag .le. 15)) goto 6000

          CALL DSCALR(DIST,DISSCL,MINDIS,MAXDIS,NOP)
          CALL TRANSP(IFLAG2,KOST,IRANKD,WX,WY,NUMNOD)

C         if an error has occurred in TRANSP, exit the program
          if ((errorflag .ge. 16) .and. (errorflag .le. 23)) goto 6000

          CALL AVFLPP(IFLAG2,EDFS,NOP,WX,WY,RN,ICOUNT,DASFC,NUMNOD)
          
C     Calculate  values of the clustering index for each unit IAFDUM() using similar formula to
C       that in subroutine OUTPU4, and noting that:	(i)	Each value of EDFS() will have its
C       own reordering of WX() and WY() from PERMUT to allow for; (ii) the ordering of AFC()
C       and ASCP remains tied to WC(), so no reordering of those two arrays is required;
C       (iii) these EDFS() values, for IFLAG=11, unlike those above for IFLAG2=1 where the
C       observed clustering indices were derived, have not been ordered in subroutine AVFLPP,
C       thus simplifying the operations required somewhat.  First, find the ordering link
C       between (X(),Y()) and (WX(),WY()) and store values in temporary TEMPAP().

          DO II = 1,NOP
            DO JJ = 1,NOP
              IF ((WX(JJ).EQ.X(II)).AND.(WY(JJ).EQ.Y(II))) THEN
                TEMPAP(JJ) = ASFP(II)
              ENDIF
   	        end do
          end do

          DO J = 1,NOP
            IF (EDFS(J).NE.0.0D0) THEN
              IAFDUM(J)=(DSIGN(1.0D0,EDFS(J)))*EDFS(J)*ASCP
     +                   /(ASFC(J)*TEMPAP(J))
              KIAFLC = KIAFLC + 1
              SORTKIA(KIAFLC) = IAFDUM(J)
              SORTNOP(J) = IAFDUM(J)
            ELSEIF (EDFS(J).EQ.0.0D0) THEN
              IAFDUM(J)=0.0D0
              KIAFLC = KIAFLC + 1
              SORTKIA(KIAFLC) = 0.0D0
              SORTNOP(J) = 0.0D0
            ENDIF
   	      end do

C     Split IAFDUM() into negative and positive values
          DO J = 1,NOP
	        if (IAFDUM(J).LT.0.0D0) then
	          IAFSN(I)=IAFSN(I)+(IAFDUM(J)/DFLOAT(NN))
C             IAFSN(I)=IAFSN(I)+(IAFDUM(J)*(IAFDUM(J).LT.0.0D0)/DFLOAT(NN))
            endif
            if (IAFDUM(J).GE.0.0D0) then
              IAFSP(I)=IAFSP(I)+(IAFDUM(J)/DFLOAT(NP))
C             IAFSP(I)=IAFSP(I)+(IAFDUM(J)*(IAFDUM(J).GE.0.0D0)/DFLOAT(NP))
            endif
   	      end do

          IAFSA(I)=((DFLOAT(NP)*IAFSP(I))+(DFLOAT(NN)*DABS(IAFSN(I))))
     +            /DFLOAT(NOP)

 9994     FORMAT(1X,F8.4,5X,F8.4,5X,F8.4,1X)
          WRITE(9,9994) IAFSN(I),IAFSP(I),IAFSA(I)


C     There follows a block of code that calculates two values for each of
C     the blue inflows and the red outflows, based on the null
C     hypothesis, i.e. on the randomised distribution.  These are:
C     (i) the 95th centile of all randomized indices, and
C     (ii) the 95th centile of the maximum value for each set
C     Details of the variables used are given below.  This is followed by
C     a rough breakdown of the steps involved in their calculation.

C     KIAFLN is the number of negative inflow indices, equal to 39*NN.
C     KIAFLP is the number of positive outflow indices, equal to 39*NP.
C     KIA is the integer that is equal to 39*NOP
C     SORTKIA() is an array, of length KIA, that holds sorted (KIAFLN
C       inflow and KIAFLP outflow) indices, used to store indices over a block of
C       39 simulations.
C     KIAFLC is a counter running from 1 to KIA that indexes the array SORTKIA().
C     INK5PU() is an array that holds the 95th centile of the negative inflow
C       indices, i.e. the int[(KIAFLN+1)/20]th point for each batch of KIA
C       values; there are K5PSIM such batches
C     INK5PM() is an array that holds the 90th centile of the negative inflow
C       indices, i.e. the int[(KIAFLN+1)/10]th point for each batch of KIA
C       values; there are K5PSIM such batches
C     INK5PL() is an array that holds the 75th centile of the negative inflow
C       indices, i.e. the int[(KIAFLN+1)/4]th point for each batch of KIA
C       values; there are K5PSIM such batches
C     AINK5PU,M,L are the (geometric mean) averages of the K5PSIM INK5PU,M,L() values
C     IPK5PU() is an array that holds the 95th centile, (i.e. the
C       {NOP + 1 - int[(KIAFLP+1)/20]}th point) of the positive outflow indices
C       in array SORTKIA(), for each batch of KIA values; there are K5PSIM such batches
C     IPK5PM() is an array that holds the 90th centile, (i.e. the
C       {NOP + 1 - int[(KIAFLP+1)/10]}th point) of the positive outflow indices
C       in array SORTKIA(), for each batch of KIA values; there are K5PSIM such batches
C     IPK5PL() is an array that holds the 75th centile, (i.e. the
C       {NOP + 1 - int[(KIAFLP+1)/4]}th point) of the positive outflow indices
C       in array SORTKIA(), for each batch of KIA values; there are K5PSIM such batches
C     AIPK5PU,M,L are the (geometric mean) averages of the K5PSIM IPK5PU,M,L() values
C	  IOUTER is an index that goes from 1 to K5PSIM
C     IINNER is an index that goes from 1 to 39
C     I is the index used in previous versions, where I = IINNER + ((IOUTER-1)*39),
C       goes from 1 to NSIMS
C     SORTNOP() is an array, of length NOP, that holds sorted (NN inflow and
C       NP outflow) indices, used temporarily for each of the NSIMS simulations.
C     IAFMAXN() is an array, indexed by I, that holds each of the maximum values of
C       the negative inflow indices over the NOP points for a single simulation.
C     IAFMAXP() is an array, indexed by I, that holds each of the maximum values of
C       the positive outflow indices over the NOP points for a single simulation.
C     IAFMXP5U is the 95th centile (the NSIMS + 1 - int[(NSIMS+1)/20]th point) of the ordered
C       array IAFMAXN()
C     IAFMXN5U is the 95th centile (the int[(NSIMS+1)/20]th point) of the ordered
C       array IAFMAXP()

C     Schedule for code:

C      Outside both loops:
C         calculate KIAFLN as 39*NN, KIAFLP as 39*NP, and KIA as 39*NOP
C         initialise variables AINK5PU,M,L and AIPK5PU,M,L

C      Outside inner loop, but inside outer loop:
C         initialise variable KIAFLC

C      Inside the inner loop:
C         assign indices as they are produced to SORTNOP() and SORTKIA().
C         sort array SORTNOP() and assign maximum and minimum values
C              to IAFMAXN() and IAFMAXP(), respectively

C      Outside inner loop, but inside outer loop:
C         sort the array SORTKIA()
C         find the positive and negative 95th centiles of sorted array SORTKIA()
C         assign these values to further arrays INK5PU,M,L() and IPK5PU,M,L()
C         increment values of AINK5PU,M,L and AIPK5PU,M,L by addition (later, to be
C             replaced by geometric mean)

C      Outside outer loop
C         divide current values of AINK5PU,M,L and AIPK5PU,M,L by K5PSIM to get means
C         sort arrays IAFMAXN() and IAFMAXP()
C         find 95th centile of sorted arrays IAFMAXN() and IAFMAXP()

C     End of schedule for code:

          CALL SORTX(SORTNOP,NOP)
          IAFMAXN(I) = SORTNOP(1)
          IAFMAXP(I) = SORTNOP(NOP)
        end do

        CALL SORTX(SORTKIA,KIA)

        IDUMMYU = (KIAFLN+1)/20
        IDUMMYM = (KIAFLN+1)/10
        IDUMMYL = (KIAFLN+1)/4
        INK5PU(IOUTER) = SORTKIA(IDUMMYU)
        INK5PM(IOUTER) = SORTKIA(IDUMMYM)
        INK5PL(IOUTER) = SORTKIA(IDUMMYL)
        AINK5PU = AINK5PU + LOG(-INK5PU(IOUTER))
        AINK5PM = AINK5PM + LOG(-INK5PM(IOUTER))
        AINK5PL = AINK5PL + LOG(-INK5PL(IOUTER))

        IDUMMYU = KIA + 1 - ((KIAFLP+1)/20)
        IDUMMYM = KIA + 1 - ((KIAFLP+1)/10)
        IDUMMYL = KIA + 1 - ((KIAFLP+1)/4)
        IPK5PU(IOUTER) = SORTKIA(IDUMMYU)
        IPK5PM(IOUTER) = SORTKIA(IDUMMYM)
        IPK5PL(IOUTER) = SORTKIA(IDUMMYL)
        AIPK5PU = AIPK5PU + LOG(IPK5PU(IOUTER))
        AIPK5PM = AIPK5PM + LOG(IPK5PM(IOUTER))
        AIPK5PL = AIPK5PL + LOG(IPK5PL(IOUTER))

      end do

      AINK5PU = AINK5PU/DFLOAT(K5PSIM)
      AINK5PU = -EXP(AINK5PU)
      AIPK5PU = AIPK5PU/DFLOAT(K5PSIM)
      AIPK5PU = EXP(AIPK5PU)

      AINK5PM = AINK5PM/DFLOAT(K5PSIM)
      AINK5PM = -EXP(AINK5PM)
      AIPK5PM = AIPK5PM/DFLOAT(K5PSIM)
      AIPK5PM = EXP(AIPK5PM)

      AINK5PL = AINK5PL/DFLOAT(K5PSIM)
      AINK5PL = -EXP(AINK5PL)
      AIPK5PL = AIPK5PL/DFLOAT(K5PSIM)
      AIPK5PL = EXP(AIPK5PL)

      CALL SORTX(IAFMAXN,NSIMS)
      CALL SORTX(IAFMAXP,NSIMS)
      
      IDUMMY = NSIMS + 1 - ((NSIMS+1)/20)
      IAFMXP5U = IAFMAXP(IDUMMY)
      IDUMMY = (NSIMS+1)/20
      IAFMXN5U = IAFMAXN(IDUMMY)

C     The value of IFLAG2 is now changed temporarily to help distinguish between
C     the three index averages, negative, positive and absolute flows,
C     when passing values to subroutine ICSORT.
C       IFLAG2=12 denotes negative inflows; usual randomizations
C       IFLAG2=13 denotes positive outflows; usual randomizations
C       IFLAG2=14 denotes all flows; usual randomizations

      WRITE(7,*)
      WRITE(7,*)
      IFLAG2=12
      WRITE(7,*) 'Negative inflows:'
      CALL ICSORT(IFLAG2,IAFSN,IAFONS,NSIMS)
      IFLAG2=13
      WRITE(7,*)
      WRITE(7,*) 'Positive outflows:'
      CALL ICSORT(IFLAG2,IAFSP,IAFOPS,NSIMS)
      IFLAG2=14
      WRITE(7,*)
      WRITE(7,*) 'Absolute values of flows:'
      CALL ICSORT(IFLAG2,IAFSA,IAFOAS,NSIMS)
      IFLAG2 = 2
      WRITE(7,*)
      WRITE(7,*) 'Statistics relating to null hypothesis distribution'
      WRITE(7,*) 'of randomized counts'
      WRITE(7,*)
      WRITE(7,*)'Centiles:               75th         90th         95th'
 1234 FORMAT (1X ,'Negative inflows     ',F8.4,5X,F8.4,5X,F8.4)
      WRITE(7,1234) AINK5PL,AINK5PM,AINK5PU
 1235 FORMAT (1X ,'Positive outflows     ',F7.4,6X,F7.4,6X,F7.4)
      WRITE(7,1235) AIPK5PL,AIPK5PM,AIPK5PU
      WRITE(7,*)
      WRITE(7,*)
 1236 FORMAT(1X ,'95th %ile of the maxima of negative inflows = ',F8.4)
      WRITE(7,1236) IAFMXN5U
 1237 FORMAT (1X ,'95th %ile of the maxima of positive outflows = '
     +          ,F8.4)
      WRITE(7,1237) IAFMXP5U

C     Now calculate distance to crowding for actual data
C     __________________________________________________


 1003 FORMAT('-------------------------------------------------',
     +'---------------------------------------------------------------')
      WRITE (7,*)
      WRITE (7,*)
      WRITE(7,1003)
      WRITE(7,*)
      WRITE(7,*)
      WRITE(7,1008)
      WRITE(7,1009)
 1008 FORMAT('Results for crowding - actual data')
 1009 FORMAT('**********************************')
      WRITE(7,*)
      WRITE(7,*)'WARNING: CROWDING INDICES ONLY GIVE INFORMATIVE VALUES'
      WRITE(7,*)'WHEN THERE IS JUST A SINGLE CLUSTER IN THE DATA - '
      WRITE(7,*)'THEY REQUIRE CAUTION IN INTERPRETATION WHEN THE'
      WRITE(7,*)'DATA COMPRISE SEVERAL CLUSTERS'
      WRITE(7,*)

      DO I = 1,NOP
         WX(I) = X(I)
         WY(I) = Y(I)
         WC(I) = ICOUNT(I)
      end do

      CALL DISCRO(DTCROO,IFLAG2)

      WRITE(9,*)
      WRITE(9,*)
      WRITE(9,*) 'The distance to crowding of the actual data is :'
      WRITE(9,1025) DTCROO
 1025 format (1X ,F21.10)
      WRITE(9,*)

      IFLAG2 = 3

C     Now start simulation of the NSIMS sets of permuted counts for distance
C     to crowding
C     ______________________________________________________________________


      WRITE(7,1010)
      WRITE(7,1011)
 1010 FORMAT('Results for crowding - permuted counts')
 1011 FORMAT('**************************************')
 
C     Loop round to get NSIMS new values of DTOCRO
C     Initialize values for random number generation so that identical
C     permutations are generated as for distance to regularity, above.
      IZ = ISEED
      IX = 30001 - IZ
      IY = 15000
      
      WRITE(9,*) 'Distances to crowding of the permuted data are :'

      DO I = 1,NSIMS
         CALL PERMUT(IX,IY,IZ)
         CALL DISCRO(DTOCRO,IFLAG2)
         DTCROS(I) = DTOCRO
         WRITE(9,1110) DTOCRO
      end do
 1110 format (1X, F21.10)

      CALL DTSORT(IFLAG2,DTCROS,DTCROO,NSIMS)

C     Report on any errors, if they occurred
 6000 continue
 
      if (errorflag .eq. 0) then     ! no errors found
        WRITE (*,*) ' '
	    WRITE (*,*) 'Program executed successfully!'
	    write (12,*)
	    write (12,*) '*** Program executed successfully! ***'
	  elseif (errorflag .gt. 0) then ! errors detected and rbno6
        write (*,*) ' '
        write (*,6100) messerror
 6100   format (2X,A)
        write (*,*) ' '
        write (*,*) 'Program execution failed!'
        if (errorflag .gt. 1) then     ! rbno6 is availble for logging
          write (12,*)
          write (12,6100) messerror
          write (12,*)
          write (12,*) '*** Program execution failed! ***'
        endif
      endif

      END

C     The end of the program.  If you have been, thank you for following
C     __________________________________________________________________




C     The start of the subroutines
C     __________________________________________________________________



C     ------------------------------------------------------------------
C                             CSCALR
C     ------------------------------------------------------------------

C     A subroutine to scale the counts, if required, via working version WC().

      SUBROUTINE CSCALR(IFLAG1,ITOT,IFLAG2)

      implicit none

      integer*4 NOP,I,ICOUNT(2000),WC(2000),IFLAG1
      integer*4 ITOT,IRANKD(2000),IFLAG2
      double precision X(2000),Y(2000),MEAN,WX(2000),WY(2000),VAR,IOD
      COMMON  /ORIG/ X,Y,WX,WY,ICOUNT,WC,IRANKD,NOP

	SAVE

      MEAN = 0.0D0
      ITOT=0

      DO I = 1,NOP
         ITOT = ITOT + ICOUNT(I)
         MEAN = MEAN + DFLOAT(ICOUNT(I))
      end do

      MEAN = MEAN/DFLOAT(NOP)

      IF (IFLAG2.EQ.0) THEN
C     Calculate sample variance and output with sample mean
         VAR = 0.0D0

         DO I = 1,NOP
            VAR=VAR+((DFLOAT(ICOUNT(I))-MEAN)*(DFLOAT(ICOUNT(I))-MEAN))
         end do

         IOD = VAR / MEAN
         VAR = VAR / (DFLOAT(NOP-1))

 1000    FORMAT (1X ,'Mean of original data = ',F12.5)
 1001    FORMAT (1X ,'Variance of original data = ',F12.5)
 1002    FORMAT (1X ,'Index of dispersion = ',F12.5, ' with ',I5,
     +                  ' degrees of freedom')
         WRITE(7,1000) MEAN
         WRITE(12,1000) MEAN
         WRITE(12,1001) VAR
         WRITE(12,1002) IOD,(NOP-1)
      ENDIF
      IF ((MEAN-DFLOAT(INT(MEAN))).GT.1.0D-12) THEN
C     Mean of counts not an integer, so scaling is required
         IFLAG1=1
         DO I = 1, NOP
            WC(I) = (NOP*WC(I)) - ITOT
         end do
      ELSE
         DO I = 1, NOP
            WC(I) = WC(I) - INT(MEAN)
         end do
      ENDIF
      RETURN
      END

C     ------------------------------------------------------------------
C                                RANKCO
C     ------------------------------------------------------------------

      SUBROUTINE RANKCO(KN,KZ,KP,DIST,MINDIS,MAXDIS,JRANKD,IFLAG2)

      implicit none

      integer*4 NOP,I,ICOUNT(2000),ISUPLY(2000),IFLAG2,WC(2000)
      integer*4 TEMPWC,TEMPI,IRANKD(2000),KN,KZ,KP,INMARC,JRANKD(2000)
      integer*4 IPOS,INEGS,IZERO,IARC,INFROM(2600000),INTO(2600000)
      integer*4 INCAT(2000),IMNODE,ITOTN,ITOTSY,INCOST(2600000),K,J
      integer*4 errorflag ! transfres error codes among sbr's

      double precision X(2000),Y(2000),DIST(2600000),MINDIS,MAXDIS
      double precision TEMPX,TEMPY,WX(2000),WY(2000)

      CHARACTER(80) messerror !for error messages

      COMMON  /ORIG/ X,Y,WX,WY,ICOUNT,WC,IRANKD,NOP
      COMMON/LINK/IMNODE,INCAT,ITOTSY,ISUPLY,INFROM,INTO,INCOST,INMARC
      common /errs/ messerror, errorflag ! transfres error codes among sbr's

	SAVE

C     Store the number of nodes, identical to the number of units
      IMNODE  = NOP
      
C     Rank counts, WC() and order WX(), WY() and original reference number
C     according to this ranking.
C     Store orig. ordered ref. nos. in IRANKD for later use.	If IRANKD(1)=4
C     this means that the 1st ranked count (i.e. the smallest) was reference number 4
C     of the list in ICOUNT().  Similarly, if IRANKD(NOP)=K, then the largest count was
C     reference number K of the original list.
C     By contrast, if JRANKD(1)=4 this means that the count in reference number 1 was
C     ranked number 4 (i.e. was the fourth smallest), and if JRANKD(K)=NOP then the largest
C     count was in reference number K.  JRANKD() will be used in the two crucial statements
C     to calculate the clustering index IAF() and IAFDUM(), in subroutine OUTPU4 and in the
C     main program respectively.

      DO K = 1, NOP
         IRANKD(K) = K
      end do

      DO K = 2, NOP
         I = K
    4       IF (I .GT. 1) THEN
               IF (WC(I-1) .GT. WC(I)) THEN
                  TEMPWC = WC(I-1)
                  TEMPX = WX(I-1)
                  TEMPY = WY(I-1)
                  TEMPI = IRANKD(I-1)
                  WC(I-1) = WC(I)
                  WX(I-1) = WX(I)
                  WY(I-1) = WY(I)
                  IRANKD(I-1) = IRANKD(I)
                  WC(I)  = TEMPWC
                  WX(I)  = TEMPX
                  WY(I)  = TEMPY
                  IRANKD(I)  = TEMPI
                  I = I - 1
               ELSE
                  I = 1
               ENDIF
               GOTO 4
            ENDIF
      end do

C     Now calculate JRANKD() from IRANKD() for original data only
      IF (IFLAG2.EQ.0) THEN
        DO I=1,NOP
          DO J=1,NOP
            IF (J.EQ.IRANKD(I)) THEN
              JRANKD(J) = I
            ENDIF
            end do
        end do
      ENDIF

C     First set up the array ISUPLY() for use by the transportation algorithm
C     Then count the number of negative counts (KN), zeroes (KZ) and positives (KP)
C     Also, get the total supply and check that it equals the total demand

      KN = 0
      KZ = 0
      KP = 0
      ITOTN = 0
      ITOTSY = 0

      DO I = 1, NOP
         ISUPLY(I) = WC(I)
         IF (WC(I).LT.0) THEN
            KN = KN + 1
            ITOTN = ITOTN + WC(I)
         ENDIF
         IF (WC(I).EQ.0) KZ = KZ + 1
         IF (WC(I).GT.0) THEN
             KP = KP + 1
             ITOTSY = ITOTSY + WC(I)
         ENDIF
      end do

      IF ((-ITOTN).NE.ITOTSY) THEN
         messerror =  '  Total supply does not equal total demand '
         errorflag = 11  !Total supply <> total demand
         return
      ENDIF
      
C     Now overwrite ITOTSY with the maximum of the supplies and the
C     absolute demands; this value will then be passed to the
C     subroutine TRANSP as an upper bound for capacity and also
C     used for scaling in DSCALR.
      ITOTSY = ABS(WC(1))
      IF (ITOTSY.LT.WC(NOP)) ITOTSY = WC(NOP)

C     Now check that number of negative, positive and zero nodes
C     sum to the total number of units.
      IF ((KN+KZ+KP).NE.NOP) THEN
         messerror =  '  Number of -, 0 and + nodes do not equal
     + total number of units'
         errorflag = 12  ! sum of neg, zero and pos nodes <> total no. units
         return
      ENDIF

      IF (KZ.EQ.NOP) THEN
        messerror = '  All data values identical - no analysis possible'
         errorflag = 13     ! All data values identical
         return
      ENDIF

C     Now set up the array INCAT() for use by the transportation algorithm
      DO INEGS = 1,KN
         INCAT(INEGS) = KP
      end do

      IF (KZ.NE.0) THEN
         DO IZERO = KN+1, KN+KZ
            INCAT(IZERO) = 0
         end do
      ENDIF

      DO IPOS = KN+KZ+1, NOP
         INCAT(IPOS) = KN
      end do

C     Now set up the arrays INFROM and INTO for use by the transportation
C     algorithm, and get the distances (costs) between the from (supply)
C     and to (demand) nodes.  Find the max and min of these, for scaling.
      IARC = 1
      DO INEGS = 1, KN
         DO IPOS = KN+KZ+1, NOP
           INFROM(IARC) = IPOS
           INTO(IARC) = INEGS
           DIST(IARC) = (WX(INEGS)-WX(IPOS))*(WX(INEGS)-WX(IPOS))
           DIST(IARC) = DIST(IARC) +
     +                  (WY(INEGS)-WY(IPOS))*(WY(INEGS)-WY(IPOS))
           DIST(IARC) = SQRT(DIST(IARC))
           IF (IARC.EQ.1) MINDIS = DIST(1)
           IF (MINDIS.GT.DIST(IARC)) MINDIS = DIST(IARC)
           IF (IARC.EQ.1) MAXDIS = DIST(1)
           IF (MAXDIS.LT.DIST(IARC)) MAXDIS = DIST(IARC)
           IARC = IARC + 1
         end do
      end do

C     Check to ensure that two units have not been entered with identical
C     coordinates
      IF (MINDIS.EQ.0.0D0) THEN
         messerror = '  Two units have same coordinates - check data'
         errorflag = 14     !  Two units have same coordinates
         return
      ENDIF

      IARC = IARC - 1

C     Check to ensure that IARC, the number of arcs, has the value KP*KN
      IF (IARC.NE.(KP*KN)) THEN
         messerror = '  No. of arcs <> no. supply nodes X'//
     +   ' no. of demand nodes'
         errorflag = 15   !  Invalid number of arcs
         return
      ENDIF

C     Store the number of arcs for later use by the transportation algorithm
      INMARC = IARC

      RETURN
      END

C     ------------------------------------------------------------------
C                                DSCALR
C     ------------------------------------------------------------------

      SUBROUTINE DSCALR(DIST,DISSCL,MINDIS,MAXDIS,NOP)

      implicit none

      integer*4 IMNODE,INCAT(2000),INFROM(2600000),INTO(2600000)
      integer*4 ITOTSY,ISUPLY(2000),INMARC,NOP,INCOST(2600000),I
      double precision DIST(2600000),MAXDIS,MINDIS,DISSCL
      COMMON/LINK/IMNODE,INCAT,ITOTSY,ISUPLY,INFROM,INTO,INCOST,INMARC

	SAVE

C     Check to see whether the ratio of max and min distances is sufficient
C     to cause a problem after scaling and translation to integer values
C     and, if so, issue a warning.

      IF ((MAXDIS/MINDIS).GT.2147483647.0D0) THEN
         WRITE(12,*) 'WARNING   ***   WARNING   ***   WARNING'
         WRITE(12,*) 'WARNING   ***   WARNING   ***   WARNING'
         WRITE(12,*)
         WRITE(12,*) 'The distances between the sample units in'
         WRITE(12,*) 'your data are so disparate - some large, '
         WRITE(12,*) 'some small - that the integer programming'
         WRITE(12,*) 'algorithm will only give approximate '
         WRITE(12,*) 'results; the minimum distance between'
         WRITE(12,*) 'units (possibly in more than one case)'
         WRITE(12,*) 'has been taken to be zero, while the'
         WRITE(12,*) '(scaled) maximum is 2,147,483,647'
         WRITE(12,*)
         WRITE(12,*) 'WARNING   ***   WARNING   ***   WARNING'
         WRITE(12,*) 'WARNING   ***   WARNING   ***   WARNING'
      ENDIF
      
C     Calculate the scaling factor, DISSCL, store and use to scale distances
C     prior to forming an integer variable INCOST() to replace DIST().

      DISSCL = 2147483646.0D0/(MAXDIS*(DFLOAT(NOP*ITOTSY)))
      
      DO I = 1, INMARC
           INCOST(I) = INT((DIST(I)*DISSCL)+0.5D0)
      end do

      RETURN
      END

C     ------------------------------------------------------------------
C                                UNSCAL
C     ------------------------------------------------------------------

      SUBROUTINE UNSCAL(KOST,ITOTSY,MAXDIS,DTOREG,IFLAG1,NOP,IFLAG2)

      implicit none

      integer*4 ITOTSY,KOST,NOP,IFLAG1,IFLAG2
      double precision DTOREG,MAXDIS

	SAVE

      DTOREG = (DFLOAT(KOST) - 0.5D0)/2147483646.0D0
      DTOREG = DTOREG * MAXDIS * DFLOAT(ITOTSY*NOP)

C     If counts have been scaled (indicated by IFLAG1), then descale.

      IF (IFLAG1.EQ.1) THEN
         DTOREG = DTOREG / DFLOAT(NOP)
      ENDIF

      IF (IFLAG2.EQ.0) THEN
 1003    FORMAT('-------------------------------------------------',
     +'---------------------------------------------------------------')
         WRITE(7,1003)
         WRITE(7,*)
         WRITE(7,*)
         WRITE(7,1000)
         WRITE(7,1001)
 1000    FORMAT('Results for regularity - actual data')
 1001    FORMAT('************************************')
         WRITE(7,1002) DTOREG
 1002    FORMAT('Distance to regularity of original data = ',F12.5)
         WRITE(7,*)
      ENDIF

      RETURN
      END

C     ------------------------------------------------------------------
C                               RANDOM
C     ------------------------------------------------------------------

      DOUBLE PRECISION FUNCTION RANDOM(IX,IY,IZ)

      implicit none

      integer*4 IX,IY,IZ
      
      save

C     AS 183 Applied Statistics (1982), 31, 188.  Wichmann & Hill's algorithm
     
      IX = 171 * MOD(IX,177) - 2 * (IX/177)
      IY = 172 * MOD(IY,176) - 35 * (IY/176)
      IZ = 170 * MOD(IZ,178) - 63 * (IZ/178)
      IF (IX.LT.0) IX = IX + 30269
      IF (IY.LT.0) IY = IY + 30307
      IF (IZ.LT.0) IZ = IZ + 30323
      RANDOM = AMOD(FLOAT(IX) / 30269.0 + FLOAT(IY) / 30307.0 +
     +              FLOAT(IZ) / 30323.0, 1.0)
     
      RETURN
      END

C     ------------------------------------------------------------------
C                                PERMUT
C     ------------------------------------------------------------------

      SUBROUTINE PERMUT(IX,IY,IZ)

      implicit none

      integer*4 ICOUNT(2000),IRANKD(2000),TEMPWC(2000)
      integer*4 I,NOP,K,ID,IE,WC(2000)
      integer*4 IX,IY,IZ,J

      double precision X(2000),Y(2000),WX(2000),WY(2000),RANDOM,DUM

      LOGICAL FULLWC(2000)
      COMMON  /ORIG/ X,Y,WX,WY,ICOUNT,WC,IRANKD,NOP

	SAVE

C     Permute the counts in ICOUNT() so that they are distributed to
C     different units from those they occupied originally.
C      L = 1

      DO I = 1,NOP
         FULLWC(I) = .FALSE.
      end do

      DO I = 1,NOP
         DUM = RANDOM(IX,IY,IZ)
         DUM = DUM*FLOAT(NOP-I+1)
         K = INT(DUM) + 1
         IE = 0
         ID = 0
           DO J = 1,NOP
              IF (.NOT.FULLWC(J)) THEN
                 ID = ID + 1
                 IF (K.EQ.ID) THEN
                    K = K + IE
                    TEMPWC(K) =  ICOUNT(I)
                    FULLWC(K) = .TRUE.
                    GOTO 2
                 ELSE
                    GOTO 3
                 ENDIF
              ELSE
                 IE = IE + 1
              ENDIF
    3      end do
    2 end do

      DO I = 1,NOP
         WC(I) = TEMPWC(I)
      end do

      RETURN
      END

C     ------------------------------------------------------------------
C                               ORDERCL
C     ------------------------------------------------------------------

      SUBROUTINE ORDERCL(DX,DY,DIAF,NN,NP,NOP)

      implicit none

      integer*4 I,K,NOP,NP,NN,J,INCR

      double precision DIAF(2000),ITEMP,XTEMP,YTEMP
      double precision DX(2000),DY(2000),ACI(2000)

	SAVE

C     This subroutine sorts the array of clustering indices and sorts the
C     X() and Y() arrays in parallel to it.  The clustering array ends
C     in ascending order, using ALGORITHM AS 304.8 APPL.STATIST. (1996),
C     VOL.45, NO.3, i.e. the first value becomes the most
C     extreme inflow index value and the NOPth becomes
C     the most extreme outflow index value.  These reordered arrays are
C     then written to a file on channel 11, but, additionally, the
C     the average cumulative value ACI() is first found, for outflows and
C     inflows separately, but placed into the same array.  This array is
C     written to column 4.  These average cumulative value, are used in
C     SURFER, or similar programs, to select an appropriate contour level.
C
      INCR = 1
C
C        Loop : calculate the increment
C
   10 INCR = 3 * INCR + 1
      IF (INCR .LE. NOP) GOTO 10

C
C        Loop : Shell-Metzner sort
C
   20 INCR = INCR / 3
      I = INCR + 1
   30 IF (I .GT. NOP) GOTO 60
        ITEMP = DIAF(I)
        XTEMP = DX(I)
        YTEMP = DY(I)
      J = I
      
   40 IF (DIAF(J - INCR) .LT. ITEMP) GOTO 50
        DIAF(J) = DIAF(J - INCR)
        DX(J) = DX(J - INCR)
        DY(J) = DY(J - INCR)
      J = J - INCR
      IF (J .GT. INCR) GOTO 40
   50   DIAF(J) = ITEMP
        DX(J) = XTEMP
        DY(J) = YTEMP
        I = I + 1
      GOTO 30
   60 IF (INCR .GT. 1) GOTO 20

C     Clustering indices now ordered
C     Find average cumulative value for negative inflows, after initialisation

      ACI(1) = DIAF(1)

      DO I = 2,NN
        ACI(I) = ((ACI(I-1)*DFLOAT(I-1))+DIAF(I))/DFLOAT(I)
      end do

C     Find average cumulative value for positive outflows, after initialisation
      ACI(NOP) = DIAF(NOP)

      DO I = 2,NP
        K = NOP - I + 1
        ACI(K) = ((ACI(K+1)*DFLOAT(I-1))+DIAF(K))/DFLOAT(I)
      end do

 1000 FORMAT(1X ,F21.10,2X,F21.10,2X,F10.5,2X,F10.5)

      DO I = 1,NOP
         WRITE(11,1000) DX(I),DY(I),DIAF(I),ACI(I)
      end do

      RETURN
      END

C     ------------------------------------------------------------------
C                               SORTX
C     ------------------------------------------------------------------

C     This subroutine just sorts a real array X() of length N of
C     positive values into ascending order
C     i.e. S(1) becomes the minimum value and S(L) the maximum
C     This replaces SORT78K, SORT6K and SORT2K in the original code -- KFC

      SUBROUTINE SORTX (X, N)

      implicit none
C
C        ALGORITHM AS 304.8 APPL.STATIST. (1996), VOL.45, NO.3
C
      integer*4 N
      double precision X(N)
C
      integer*4 I, J, INCR
      double precision TEMP

	SAVE

      INCR = 1
C
C        Loop : calculate the increment
C
   10 INCR = 3 * INCR + 1
      IF (INCR .LE. N) GOTO 10

C        Loop : Shell-Metzner sort
C
   20 INCR = INCR / 3
      I = INCR + 1
   30 IF (I .GT. N) GOTO 60
      TEMP = X(I)
      J = I
   40 IF (X(J - INCR) .LT. TEMP) GOTO 50
      X(J) = X(J - INCR)
      J = J - INCR
      IF (J .GT. INCR) GOTO 40
   50 X(J) = TEMP
      I = I + 1
      GOTO 30
   60 IF (INCR .GT. 1) GOTO 20
C
      RETURN
      END


C     ------------------------------------------------------------------
C                               DTSORT
C     ------------------------------------------------------------------

      SUBROUTINE DTSORT(IFLAG2,DTS,DTO,NSIMS)

      implicit none

      integer*4 NSIMS,I,K,JM1,NMJP1,IFLAG2,INCR,J

      double precision DTO,DTS(6000),TEMP,PA,IOA


C     Get average value of DTS() to calculate index of aggregation

      IOA = DTS(1)

      DO K = 2, NSIMS
         IOA = IOA + DTS(K)
      end do

C        Sort DTS() into ascending order using
C        ALGORITHM AS 304.8 APPL.STATIST. (1996), VOL.45, NO.3
C
C        Sorts the N values stored in array X in ascending order

      INCR = 1
C
C        Loop : calculate the increment
C
   10 INCR = 3 * INCR + 1
      IF (INCR .LE. NSIMS) GOTO 10

C        Loop : Shell-Metzner sort
C
   20 INCR = INCR / 3
      I = INCR + 1
   30 IF (I .GT. NSIMS) GOTO 60
        TEMP = DTS(I)
        J = I
   40 IF (DTS(J - INCR) .LT. TEMP) GOTO 50
        DTS(J) = DTS(J - INCR)
        J = J - INCR
      IF (J .GT. INCR) GOTO 40
   50   DTS(J) = TEMP
        I = I + 1
      GOTO 30
   60 IF (INCR .GT. 1) GOTO 20

C     Now find position w.r.t. this ordering of DTO

      IF (DTS(NSIMS) .LT. DTO) THEN
C     All simulated values are less than the observed, so the observed
C      pattern is as aggregated as it is possible to be
        JM1 = NSIMS
        NMJP1 = 0
C        PR = 1.0D0
        PA = 0.0D0
      ELSE
C     Not all simulated values are less than the observed
        DO I = 1, NSIMS
            IF (DTS(I) .LT. DTO) THEN
               GOTO 2
            ELSE
               JM1 = I - 1
               NMJP1 = NSIMS - I + 1
C               PR = DBLE(JM1)/DBLE(NSIMS)
               PA = DBLE(NMJP1)/DBLE(NSIMS)
               GOTO 3
            ENDIF
    2   end do
    3 ENDIF

C     Now calculate index of aggregation
      IOA = IOA/DBLE(NSIMS)

C     Now output results
      IF (IFLAG2.EQ.3) THEN
         WRITE(7,*)  'Results for distance to crowding'
         WRITE (7,*) 'Permutations conditioning on observed counts'
         WRITE(7,*)
      ENDIF

      WRITE (7,999) NSIMS
      WRITE (7,997) JM1
      WRITE (7,998) NMJP1
  997 format (7x, I5,'  had total distances less than',
     +' the observed, actual data, while ')
  998 format (7x, I5,'  had total distances greater',
     +' than the observed.')
  999 FORMAT (1X ,'Of the',I5,' permutations: ')
      write (7,*)

      IF (IFLAG2.EQ.1) THEN
       write (10,*) 'P sub a, I sub a, mean v_j, mean v_i,',
     +' P(mean v_j), P(mean v_i)'
 1000 FORMAT(1X ,'The probability that the observed data is no more',/,
     +' aggregated than expected from a random permutation of ',/,
     +' the counts, i.e. P sub a, is: ', F7.4)
 1001  FORMAT(1X ,'The probability that the observed data is no more',/,
     +' aggregated than expected from a random permutation of ',/,
     +'the counts, i.e. P sub a, is less than: ', F7.4)
 1002  FORMAT(1X ,'The probability that the observed data is no more',/,
     +' aggregated than expected from a random permutation of ',/,
     +' the counts, i.e. P sub a, is more than: ', F7.4)
 2000	 FORMAT(1X ,F7.4)

      IF (PA.NE.0.0D0.AND.PA.NE.1.0D0) THEN
          WRITE (7,1000)	PA
          WRITE (10,2000)	PA
       ENDIF
       IF (PA.EQ.0.0D0) THEN
          PA = 1.0D0/DBLE(NSIMS)
          WRITE (7,1001) PA
          WRITE (10,2000) PA
       ENDIF
       IF (PA.EQ.1.0D0) THEN
          PA = 1.0D0 - (1.0D0/DBLE(NSIMS))
          WRITE (7,1002) PA
          WRITE (10,2000) PA
       ENDIF
      ENDIF

      IF (IFLAG2.EQ.3) THEN
       IF (PA.NE.0.0D0.AND.PA.NE.1.0D0) THEN
       WRITE (7,*) 'The probability that the observed data is no more',
     +   ' aggregated than expected from a random permutation of ',
     +   'the counts, i.e. Q sub a, is: ', 1.0D0-PA
       ENDIF

       IF (PA.EQ.0.0D0) THEN
          PA = 1.0D0/DBLE(NSIMS)
       WRITE (7,*) 'The probability that the observed data is no more',
     +   ' aggregated than expected from a random permutation of ',
     +   'the counts, i.e. Q sub a, is more than: ', 1.0D0-PA
       ENDIF

       IF (PA.EQ.1.0D0) THEN
          PA = 1.0D0 - (1.0D0/DBLE(NSIMS))
       WRITE (7,*) 'The probability that the observed data is no more',
     +   ' aggregated than expected from a random permutation of ',
     +   ' the counts, i.e. Q sub a, is less than: ', 1.0D0-PA
       ENDIF
      ENDIF

      WRITE (7,*)
      IF (IFLAG2.EQ.1) THEN
         WRITE (7,1003) IOA
 1003 FORMAT(1X ,'Average distance to regularity of simulated data = ',
     +F12.5)
        IOA = DTO/IOA
        WRITE (7,*)
        WRITE (7,1004)  IOA
        WRITE (10,2001)  IOA
 1004   FORMAT(1X ,'The index of aggregation, I sub a, is: ', F6.3)
 2001	FORMAT(1X ,F6.3)
        WRITE (7,*)
        WRITE (7,*)
        WRITE(7,1005)
 1005 FORMAT('-------------------------------------------------',
     +'---------------------------------------------------------------')
      ENDIF

      IF (IFLAG2.EQ.3) THEN
         write (7,*) ' The average distance to crowding of the',
     +             ' simulated data was: '
         WRITE (7,1303) IOA
 1303    FORMAT(1X, F21.10)
         IOA = IOA/DTO
         WRITE (7,*)
         WRITE (7,1304) IOA
 1304    FORMAT(1X ,'The index of aggregation, J sub a, is: ',F21.10)
      ENDIF

      RETURN
      END

C     ------------------------------------------------------------------
C                               ICSORT
C     ------------------------------------------------------------------

      SUBROUTINE ICSORT(IFLAG2,ICS,ICO,NSIMS)

      implicit none

      integer*4 NSIMS,I,K,JM1,NMJP1,IFLAG2,INCR,J

      double precision ICO,ICS(6000),TEMP,PA,PR,IOC

	SAVE

C     Sort ICS() into ascending order
C     Get average value of ICS() to compare with ICO from
C     original data.  The value of IFLAG2 is used to distinguish between
C     the three index averages, negative, positive and absolute flows:
C       IFLAG2=12 denotes negative inflows
C       IFLAG2=13 denotes positive outflows
C       IFLAG2=14 denotes all flows
C
      IOC = ICS(1)

      DO K = 2, NSIMS
         IOC = IOC + ICS(K)
      end do

C        Sort ICS() into ascending order using
C        ALGORITHM AS 304.8 APPL.STATIST. (1996), VOL.45, NO.3
C
      INCR = 1
C
C        Loop : calculate the increment
C
   10 INCR = 3 * INCR + 1
      IF (INCR .LE. NSIMS) GOTO 10

C        Loop : Shell-Metzner sort
C
   20 INCR = INCR / 3
      I = INCR + 1
   30 IF (I .GT. NSIMS) GOTO 60
        TEMP = ICS(I)
        J = I
   40 IF (ICS(J - INCR) .LT. TEMP) GOTO 50
        ICS(J) = ICS(J - INCR)
        J = J - INCR
      IF (J .GT. INCR) GOTO 40
   50   ICS(J) = TEMP
        I = I + 1
      GOTO 30
   60 IF (INCR .GT. 1) GOTO 20

C     Now find position w.r.t. this ordering of ICO
      IF (ICS(NSIMS) .LT. ICO) THEN
C       All simulated values are less than the observed ICO
        JM1 = NSIMS
        NMJP1 = 0
        PR = 1.0D0
        PA = 0.0D0
      ELSE
C       Not all simulated values are less than the observed ICO
        DO I = 1, NSIMS
            IF (ICS(I) .LT. ICO) THEN
               GOTO 2
            ELSE
               JM1 = I - 1
               NMJP1 = NSIMS - I + 1
               PR = DBLE(JM1)/DBLE(NSIMS)
               PA = DBLE(NMJP1)/DBLE(NSIMS)
               GOTO 3
            ENDIF
    2     end do
    3 ENDIF

C     Now calculate average degree of clustering
      IOC = IOC/DBLE(NSIMS)

C     Now output results

      IF (IFLAG2.EQ.12) THEN
  999   FORMAT (1X ,'Of the',I5,' permutations, ',I5, ' had average ')
        WRITE (7,999) NSIMS,JM1
        WRITE (7,*) 'clustering negatively greater or equal to the',
     + ' observed data,'
        WRITE (7,101) PR
 2000   FORMAT(1X ,F6.4)
        WRITE (10,2000) PR
 101    FORMAT (1X ,'a proportion of ',F6.4,' of the permutations.')
      ENDIF
c
      IF (IFLAG2.EQ.13) THEN
        WRITE (7,999) NSIMS,NMJP1
        WRITE (7,*) 'clustering greater or equal to the observed data,'
        WRITE (7,101) PA
        WRITE (10,2000) PA
      ENDIF

      IF (IFLAG2.EQ.14) THEN
        WRITE (7,999) NSIMS,NMJP1
        WRITE (7,*) 'clustering greater or equal to the observed data,'
        WRITE (7,101) PA
      ENDIF

      RETURN
      END

C     ------------------------------------------------------------------
C                                DISCRO
C     ------------------------------------------------------------------

      SUBROUTINE DISCRO(DTOCRO,IFLAG2)

      implicit none

      integer*4 I,NOP,ICOUNT(2000),IRANKD(2000),IFLAG2,IDUM,J,WC(2000)

      double precision DTOCRO,DUM1,DUM2
      double precision X(2000),Y(2000),WY(2000),WX(2000)

      COMMON  /ORIG/ X,Y,WX,WY,ICOUNT,WC,IRANKD,NOP


C     Select each unit in turn, calculate the distance to crowding for each,
C     and store the minimum in DTOCRO
      DO I = 1,NOP
         DUM1 = 0.0D0
         DO J = 1,NOP
           DUM2 = (WX(J)-WX(I))*(WX(J)-WX(I))
           DUM2 = DUM2 + ((WY(J)-WY(I))*(WY(J)-WY(I)))
           DUM2 = SQRT(DUM2)
           DUM1 = DUM1 + (DUM2*WC(J))
         end do
         IF (I.EQ.1) THEN
            DTOCRO = DUM1
            IDUM = 1
         ELSE
            IF (DUM1.LT.DTOCRO) THEN
               DTOCRO = DUM1
               IDUM = I
            ENDIF
         ENDIF
      end do

      IF (IFLAG2.EQ.2) THEN
C     Output results for actual data
         WRITE(7,1000) DTOCRO
1000     format ('Distance to crowding for actual data = ',F21.10)
         WRITE(7,*) 'with unit ',IDUM,' of the original data as the ',
     +              'focus of the crowding'
         WRITE(7,*)
         WRITE(7,*)
      ENDIF
      RETURN
      END

C     ------------------------------------------------------------------
C                                AVFLPP
C     ------------------------------------------------------------------

      SUBROUTINE AVFLPP(IFLAG2,EDF,NOP,X,Y,RN,ICOUNT,DAF,NUMNOD)

      implicit none

      integer*4 I,J,K,NOP,NUMNOD,STFROM(2000),STTO(2000),STFLOW(2000)
      integer*4 TNF(4000),DUMFRO(4000),DUMTO(4000),DUMFLO(4000)
      integer*4 D1,D2,DCOUNT,RN(2000),IFLAG2,TNN,TF(4000)
      integer*4 ICOUNT(2000),DUM1,DUM2,ITEMP

      double precision EDF(2000),TFTD(2000),DD,X(2000),Y(2000)
      double precision DINCL(4000),DUMDIS(4000),DAF(2000)
      double precision DINCX(4000),DINCY(4000),TEMP

      common  /FLOWS/STFROM,STTO,STFLOW


C     NOTE: TEMP and ITEMP definitions added 30/10/07,kfc -- prior bug?

C     Subroutine to calculate the average flow per sample unit, and determine
C     whether positive (outflow) or negative (inflow)
c
C     Number of nodes calculated in TRANSP and read into NUMNOD.
C     Calculate TNN as twice NUMNOD
      TNN = 2 * NUMNOD
C     Initialize internal dummy counter for printed output
      DCOUNT = 0

C     Loop around the list of NOP units
      DO I = 1,NOP
C     Then initialize integer total number of flows, total flow x dist, and total flow
C     for each unit
         TNF(I) = 0
         TFTD(I) = 0.0D0
         TF(I) = 0
C     and the real DINCX,the scaled incremental vector flow for X
         DINCX(I) = 0.0D0
C     and the real DINCY,the scaled incremental vector flow for Y
         DINCY(I) = 0.0D0
     	end do

C     Loop around the list of NOP units
      DO I = 1,NOP
C     and within each unit, loop around the list of NUMNOD nodes
        DO K = 1,NUMNOD
C     First find if the next flow in the list of NUMNOD flows is concerned with
C     unit I and, if so, whether it is a 'from' (positive) or 'to' (negative) flow.
            IF (STFROM(K).EQ.I) THEN
              DCOUNT = DCOUNT + 1
              J = STTO(K)
C     I is the index of the original unit, and the flow is from I to J (outflow)
C     Increment total number of flows for unit I; total (flow x dist) for unit I; and
C     total flow for unit I; similarly for unit J.
              TNF(I) = TNF(I) + 1
              TNF(J) = TNF(J) + 1
C     Note that TF() is defined with positive values for outflows and negative for inflows.
              TF(I) = TF(I) + STFLOW(K)
              TF(J) = TF(J) - STFLOW(K)
C     now calculate DD the distance FROM to TO (I to J)
              DD = 0.0D0
              DD = DD + ((X(I)-X(J))*(X(I)-X(J)))
              DD = DD + ((Y(I)-Y(J))*(Y(I)-Y(J)))
              DD = DSQRT(DD)
C     Note that TFTD() is always positive for inflows and outflows
              TFTD(I) = TFTD(I) + (DD*DFLOAT(STFLOW(K)))
              TFTD(J) = TFTD(J) + (DD*DFLOAT(STFLOW(K)))
C     now store some variables temporarily for printing
              DUMFRO(DCOUNT) = I
              DUMTO(DCOUNT) = J
              DUMFLO(DCOUNT) = STFLOW(K)
              DUMDIS(DCOUNT) = DD
              GOTO 3
            ENDIF
            IF (STTO(K).EQ.I) THEN
               DCOUNT = DCOUNT + 1
               J = STFROM(K)
C     I is the index of the original unit, and the flow is from J TO I (inflow)
C     Increment total number of flows for unit I; total (flow x dist) for unit I; and
C     total flow for unit I; similarly for unit J.
               TNF(I) = TNF(I) + 1
               TNF(J) = TNF(J) + 1
C     Note that TF() is defined with positive values for outflows and negative for inflows.
                TF(I) = TF(I) - STFLOW(K)
                TF(J) = TF(J) + STFLOW(K)
C     Now calculate DD, the distance of the vector of FROM to TO (J to I)
                DD = 0.0D0
                DD = DD + ((X(I)-X(J))*(X(I)-X(J)))
                DD = DD + ((Y(I)-Y(J))*(Y(I)-Y(J)))
                DD = DSQRT(DD)
C     Note that again TFTD() is always positive for inflows and outflows
                TFTD(I) = TFTD(I) + (DD*DFLOAT(STFLOW(K)))
                TFTD(J) = TFTD(J) + (DD*DFLOAT(STFLOW(K)))
C     now store some variables temporarily for printing
                DUMFRO(DCOUNT) = J
                DUMTO(DCOUNT) = I
                DUMFLO(DCOUNT) = STFLOW(K)
                DUMDIS(DCOUNT) = DD
            ENDIF
C     End of calculations for this particular nodal flow
    3    end do
C     End of loop around units
      end do

C     Determine the total number of flows and the total flow for each sample unit,
C     allowing for the fact that each flow has been counted twice in the double loop above
      DO I=1,NOP
         TNF(I) = TNF(I)/2
         TF(I)  =  TF(I)/2
      end do

C     Now, if original data, print some output ordered by ref. no. of units
      IF (IFLAG2.EQ.0) THEN
        WRITE(9,*)
        WRITE(9,*)
        WRITE(9,*) 'All flows, unit by unit'
        WRITE (9,*) 'There follows in column order:'
        WRITE (9,*)'1: Ref. No. of unit'
        WRITE (9,*)'2: Ref. No. of FROM unit'
        WRITE (9,*)'3: x-coordinate of FROM unit'
        WRITE (9,*)'4: y-coordinate of FROM unit'
        WRITE (9,*)'5: Ref. No. of TO unit'
        WRITE (9,*)'6: x-coordinate of TO unit'
        WRITE (9,*)'7: y-coordinate of TO unit'
        WRITE (9,*)'8: flow'
        WRITE (9,*)'9: distance of flow'
        WRITE(9,*)
 9999   FORMAT(1X ,I4,1X,2(I4,1X,2(F16.7,1X)),I10,1X,F12.6)
C     The index variable DUM1 refers to the sample unit concerned
        DUM1=1
C       and the integer variable DUM2 refers to the incremental number of flows
        DUM2=1
        DO K = 1,TNN
          D1 = DUMFRO(K)
          D2 = DUMTO(K)
          IF (TNF(DUM1).NE.0) THEN
       DINCX(DUM1) = DINCX(DUM1) + ((X(D2) - X(D1))*DUMFLO(K)/TNF(DUM1))
       DINCY(DUM1) = DINCY(DUM1) + ((Y(D2) - Y(D1))*DUMFLO(K)/TNF(DUM1))
          ENDIF
      WRITE(9,9999) DUM1,DUMFRO(K),X(D1),Y(D1),DUMTO(K),X(D2),Y(D2),
     +                  DUMFLO(K),DUMDIS(K)
             IF (DUM2.EQ.TNF(DUM1)) THEN
                WRITE(9,*)
C      now reduce the x- and y-vector flows to unit length
                DINCL(DUM1) = SQRT((DINCX(DUM1)*DINCX(DUM1))+
     +                       (DINCY(DUM1)*DINCY(DUM1)))
                IF (DINCL(DUM1).EQ.0.0D0) THEN
                  DINCX(DUM1) = 0.0D0
                  DINCY(DUM1) = 0.0D0
                ELSE
                DINCX(DUM1) = DINCX(DUM1)/DINCL(DUM1)
                DINCY(DUM1) = DINCY(DUM1)/DINCL(DUM1)
              ENDIF
              DUM1 = DUM1 + 1
              DUM2 = 1
            ELSE
              DUM2 = DUM2 + 1
            ENDIF
        end do
      ENDIF

C     And, again if original data, print some output, but now ordered by average flow distance
C     First some printed headings:
      IF (IFLAG2.EQ.0) THEN
        WRITE(9,*)
        WRITE(9,*)
        WRITE(9,*)  'Summary of Flows - unit by unit'
        WRITE (9,*) 'There follows in column order:'
        WRITE (9,*)'1: Ref. No. of unit'
        WRITE (9,*)'2: x-coordinate of unit'
        WRITE (9,*)'3: y-coordinate of unit'
        WRITE (9,*)'4: count'
        WRITE (9,*)'5: number of flows'
        WRITE (9,*)'6: total flow'
        WRITE (9,*)'7: average distance of flow'
        WRITE (9,*)'8: vector flow in x-direction with unit length'
        WRITE (9,*)'9: vector flow in y-direction with unit length'
        WRITE(9,*)
      ENDIF

 9998 FORMAT(1X ,I4,1X,F16.7,1X,F16.7,1X,I6,1X,I3,1X,I10,1X,
     +       F16.7,2(1X,F9.6))
C         and, if the original data, then print out.
      DO I = 1,NOP
C       Bear in mind that flows have been counted twice when stored in TFTD() and now correct
C       before storing average flow distance in EDF.
C       Also that while TFTD() is always positive, TF() differs according to inflow or outflow
        IF (TF(I).NE.0) THEN
          EDF(I) = TFTD(I)/(2.0D0*DFLOAT(TF(I)))
        ENDIF
        IF (TF(I).EQ.0) THEN
          EDF(I) = 0.0D0
        ENDIF
C       Hence EDF(I) is positive for outflow and negative for inflow.
C       At this stage EDF(I) contains the average flow distance for the Ith sample unit
C       For the original data only, calculate RN() for use throughout the next simulation phase
C       of the program. This provides unique reference numbers for later use.
        IF (IFLAG2.EQ.0) THEN
           RN(I) = I
        ENDIF
C       Print out flow variables for each unit, if original data.
        IF (IFLAG2.EQ.0) THEN
            WRITE(9,9998) I,X(I),Y(I),ICOUNT(I),
     +                  TNF(I),TF(I),EDF(I),DINCX(I),DINCY(I)
C       Explanation of these variables is given above underneath the section write
C        statement that prints: Summary of Flows - unit by unit
        ENDIF
C       Store values of EDF() in DAF().  DAF() is passed back to main program
C       as DASF() if simulated data
C       If this is original data, then DAF() is also passed back as DASF()
C        but then the variable is unused.
C       DAF(I) is, like current EDF(I), the average flow distance for the Ith sample unit
        DAF(I) = EDF(I)
      end do

C     Now, unless IFLAG2=11,  the array EDF(), which holds the average distance of flow
C     for each unit will be ordered.  If IFLAG2=11, there is no point in ordering EDFS,
C     since there is no detailed output.
C     Note that:
C      (i)  DAF remains unordered.  Its value relates to a constant ICOUNT but different
C          (x,y) coordinates on each permutation, through the WX() and WY() values.
C      (ii) RN() is a reference number array holding the unit numbers of the
C              ordered average distances; RN() is only relevant to the original data,
C              so it is not altered unless this is original data.  If that IS the case
C              then it too is sorted in parallel to EDF().
c
C        Sort EDF() into ascending order using
C        ALGORITHM AS 304.8 APPL.STATIST. (1996), VOL.45, NO.3
C        Note that this splits all outflows from all inflows automatically.
c
c
c      IF (IFLAG2.NE.11) THEN
c
c        INCR = 1
c
cC       Loop : calculate the increment
c
c   15   INCR = 3 * INCR + 1
c        IF (INCR .LE. NOP) GOTO 15
c
cC         Loop : Shell-Metzner sort
c
c   20     INCR = INCR / 3
c          I = INCR + 1
c   30     IF (I .GT. NOP) GOTO 60
c            TEMP = EDF(I)
c            ITEMP = RN(I)
c            J = I
c   40       IF (EDF(J - INCR) .LT. TEMP) GOTO 50
c              EDF(J) = EDF(J - INCR)
c              IF (IFLAG2.EQ.0) THEN
c                RN(J) = RN(J - INCR)
c              ENDIF
c              J = J - INCR
c            IF (J .GT. INCR) GOTO 40
c   50       EDF(J) = TEMP
c            IF (IFLAG2.EQ.0) THEN
c             RN(J) = ITEMP
c            ENDIF
c            I = I + 1
c            GOTO 30
c   60     IF (INCR .GT. 1) GOTO 20
c
c        ENDIF

c         Restored the bubble sort after finding that the shell sort appeared
c            to increase discrepancy with older version because of TEMP
c            rounding error  -- KFC 22/07/08

      IF (IFLAG2.NE.11) THEN
         DO 12 K = 2, NOP
	      I = K
   13       IF (I .GT. 1) THEN        
               IF (EDF(I-1) .GT. EDF(I)) THEN
                  TEMP = EDF(I-1)
                  IF (IFLAG2.EQ.0) THEN
                     ITEMP = RN(I-1)
                  ENDIF
                  EDF(I-1) = EDF(I)
                  IF (IFLAG2.EQ.0) THEN
				   RN(I-1) = RN(I)
                  ENDIF
                  EDF(I)  = TEMP
	            IF (IFLAG2.EQ.0) THEN
                     RN(I) = ITEMP
                  ENDIF
                  I = I - 1
               ELSE
                  I = 1
               ENDIF
               GOTO 13
            ENDIF 
   12    CONTINUE
      ENDIF


      RETURN
      END


C     ------------------------------------------------------------------
C                          displayhelp
C     ------------------------------------------------------------------

      SUBROUTINE displayhelp ()

      write(*,*)  ' '
      write(*,*)  '    Program RBRel14'
      write(*,*)  ' '
      write(*,*)  '    Usage: rbrel14 [-r] [-h] [-i]'
      write(*,*)  ' '
      write(*,*)  '      Performs a Red-Blue SADIE analysis'
      write(*,*)  ' '
      write(*,*)  '      Optionally uses input files rbni5.dat',
     +             ' rbni8.dat'
      write(*,*)  ' '
      write(*,*)  '      Optional command-line parameters:'
      write(*,*)  ' '
      write(*,*)  '       -r Analyses ranks of counts rather',
     +                   ' than the counts'
      write(*,*)  '          (The "non-parametric method")'
      write(*,*)  ' '
      write(*,*)  '       -i Enters interactive mode with',
     +                   ' prompts for input'
      write(*,*)  ' '
      write(*,*)  '       -h Displays this text'
      write(*,*)  ' '
      write(*,*)  '      Only one option may be used at a time'
      write(*,*)  ' '
      write(*,*)  ' '

      return
      end


C     ------------------------------------------------------------------
C                                Interactor
C     ------------------------------------------------------------------
C     Collects 'interactive mode' input from the user

      subroutine Interactor(projpath,infile5,rflag)

      USE PORTLIB            !comment out for g77 compiler

      integer*4 istat        ! for finding working directory

      Logical RFLAG

      character(256) infile5  ! filename for input (usu rbni5.dat)
      character(1) yesno      ! input y or n
      character(256) projpath ! project path (usu current folder)
      CHARACTER(256) currdir  ! placeholder for working directory

	SAVE

	currdir = ' '
	istat = getcwd(currdir)
      print *, ' '
      print *, '*** Program Red-Blue version 1.4 *** '
      print *, ' '
      print *, 'Enter project path for output files'
	print *, 'Use . for current directory'
	write (*,20)
  20  FORMAT (' --> ', $)
      READ (*,*) projpath
	istat = LNBLNK(projpath)
      if (projpath(:istat) .eq. '.') then
        projpath = currdir(:LNBLNK(currdir))
      endif
 	print *, ' '
      print *, 'Enter the file name for input'
	print *, '(full pathnames permitted)'
	WRITE (*, 20)
      READ (*,*) infile5
	istat = scan(infile5,':\')
	if (istat .eq. 0) then
	  infile5 = projpath(:LNBLNK(projpath)) // '\'
     +            // infile5(:LNBLNK(infile5))
	endif
      print *, ' '
      print *, 'Use non-parametric method? (y/n)'
	WRITE (*, 20)
      READ (*,*) yesno
      if (yesno .eq. 'y') then
        rflag=.true.
      else
        rflag=.false.
      endif

	return
	end


C     ------------------------------------------------------------------
C                                TRANSP
C     ------------------------------------------------------------------
C
C     The code for this subroutine and those which it calls
C     was adapted from that very kindly supplied by Les Proll (address below)
C
C     ------------------------------------------------------------------
C     PROGRAM TO SOLVE MINIMUM COST NETWORK FLOW PROBLEMS
C     ADAPTED FROM THE CODE NETFLO GIVEN IN
C     KENNINGTON & HELGASON - ALGORITHMS FOR NETWORK PROGRAMMING
C     (WILEY,1980), P244-277.
C
C     LIBRARIES REQUIRED - LGPLIB
C
C     L.G.PROLL
C     DEPT. OF COMPUTER STUDIES
C     UNIVERSITY OF LEEDS
C     VERSION 2.3
C     APRIL 1983
C     ------------------------------------------------------------------
C
C     NODE FORMAT
C
C     1 - PREDECESSOR OR DOWN POINTER
C     2 - THREAD OR NEXT POINTER
C     3 - LEVEL NUMBER
C     4 - ASSOCIATED ARC IDENTIFIER
C     ( + = ORIENTATION OPPOSITE TO DOWN POINTER )
C     ( - = ORIENTATION SAME AS DOWN POINTER )
C     5 - FLOW ON ARC
C     6 - DUAL VARIABLE VALUE
C
C     ARC FORMAT
C
C     1 - FROM NODE
C     2 - COST
C     3 - CAPACITY(-,IF AT UB)
C     4 - LOWER BOUND
C     5 - NAME
C
C     NSTOP - EXIT CONDITION
C     MNODE - NUMBER OF NODES
C     NET - TOTAL BALANCE IN NETWORK
C     MSORC - NUMBER OF SOURCES
C     MSINK - NUMBER OF SINKS
C     MARC - NUMBER OF ARCS
C     MTREE - NUMBER OF BRANCHES ON TREE (EXCLUDING ROOT)
C     THD - POINTER MOVING ALONG THREAD
C     TRY - VARIABLE ENCOUNTERED DURING SETUP OR PRICING
C     PRICE - REDUCED COST FOR TRY
C     NEWARC - BEST VARIABLE TO ENTER
C     NEWPR - PRICE FOR NEWARC
C     DW - DOWN PTRS FOR RATIO TEST (FROM STEM, TO STEM)
C     CH - PATH CONDITIONS FOR RATIO TEST (FROM STEM, TO STEM)
C     DWN - POINTER MOVING ALONG DOWN PATH
C     CHG - PATH CONDITION
C     THETA - MINIMUM RATIO IN RATIO TEST
C     JTHETA - UB (-) OR LB (+) CONDITION FOR MIN THETA
C     KTHETA - MIN THETA OCCURS ON FROM STEM (1) OR TO STEM (2)
C     POSS - CANDIDATE FOR MIN THETA
C     DWE - ROOT OF CYCLE
C     ------------------------------------------------------------------

      SUBROUTINE TRANSP(IFLAG2,KOST,IRANKD,X,Y,NUMNOD)

      implicit none

      LOGICAL OPTIM

      integer*4 ZERO,BIG,IFLAG2,IRANKD(2000),NUMNOD,LARCP1
      integer*4 MTREE,MSLK,MREG,KOST0,MNODP1,MNODP2,LNODP1
      integer*4 MSORC,ITER
      integer*4 LNODE,MNODE
      integer*4 LARC,MARC,KOST
      integer*4 errorflag ! transfers error codes among sbr's

      double precision X(2000),Y(2000)
      
      CHARACTER(80) messerror !for error messages

      COMMON/PARAMI/LNODE,LARCP1,BIG
      common /errs/ messerror, errorflag ! transfres error codes among sbr's

	SAVE

      KOST = 0
      LNODE=2000
      LARCP1=2600000
      BIG=1000000000
      MTREE=0

      CALL INDAT(MNODE,MARC,MSORC,LARC,LNODP1,MNODP1,MNODP2,KOST0)
      
      if ((errorflag .ge. 18) .and. (errorflag .le. 23)) then
        return
      endif

      CALL TRANSF(MARC,MSORC,LARC,LNODP1,MNODP1,MNODP2,MREG,MSLK,ZERO)

      if (errorflag .eq. 16) then
        return
      endif
      
      CALL LWBND(MNODE,MARC,MNODP1,MTREE)
      
      if (errorflag .eq. 17) then
        return
      endif
      
      CALL START(MNODE,MARC,LNODP1,MNODP1,MTREE,ZERO)
      CALL KHFLOW(LNODP1,MNODP2,ITER,OPTIM)
      CALL OUT(MNODE,ZERO,LNODP1,MNODP1,MSLK,KOST0,OPTIM,IFLAG2,
     +        KOST,IRANKD,X,Y,NUMNOD)

      RETURN
      END


C     ------------------------------------------------------------------
C                                   XOR
C     ------------------------------------------------------------------

      integer*4 FUNCTION TXOR(I,J)

      implicit none

      integer*4 I,J

      TXOR = I*J
      IF (TXOR.EQ.0) TXOR = I+J
      RETURN
      END


C     -------------------------------------------------------------------
C                                 TRANSF
C     -------------------------------------------------------------------

C     SUBROUTINE TRANSF TRANSFORMS THE NETWORK BY REMOVING REDUNDANT DUMMY
C     ARCS, ADDING SLACKS AND INTRODUCING AN ARTIFICIAL ARC TO PROVIDE AN
C     INITIAL BASIS

      SUBROUTINE TRANSF(MARC,MSORC,LARC,LNODP1,MNODP1,MNODP2,MREG,MSLK
     +                  ,ZERO)

      implicit none

      integer*4 ZERO,BIG,TXOR,THD
      integer*4 LNODE,DOWN(2000),NEXT(2000),LEVEL(2000)
      integer*4 ARCID(2000),FLOW(2000),DUAL(2000),CAT(2000)
      integer*4 FROM(2600000),COST(2600000),CAPAC(2600000)
      integer*4 I,J,K,L,J400,MREG,MARC,LARC,MNODP1,NXT,MSLK,MSORC
      integer*4 MNODP2,LARCP1,LNODP1,TFLOOR(2600000)
      integer*4 errorflag ! transfres error codes among sbr's

      CHARACTER(8) SLACK,ARTIF,DUMMY,XCESS,NAME
      CHARACTER(80) messerror !for error messages

      COMMON/NODINF/DOWN,NEXT,LEVEL,ARCID,FLOW,DUAL,CAT
      COMMON/ARCINI/FROM,COST,CAPAC,TFLOOR
      COMMON/PARAMI/LNODE,LARCP1,BIG
      COMMON/ARCINR/NAME(2600000)
      common /errs/ messerror, errorflag ! transfres error codes among sbr's

	SAVE

      EXTERNAL TXOR

C     ELIMINATE NON-ESSENTIAL DUMMY ARCS
C     CAT(N) WILL POINT TO LOCATION OF FIRST ARC TERMINAL AT NODE N

      I = LNODP1
      K = 0
      L = 0
      SLACK='SLACK'
      ARTIF='ARTIFIC'
      DUMMY='DUMMY'
      XCESS='EXCESS'

      MARC = MARC-1

      DO J400 = 1,MARC
        J = FROM(J400)
        IF (TXOR(I,J).LE.0) then
          I = -I
          L = L+1
          CAT(L) = K+1
          K = K+1
          IF (K.NE.J400) then
            FROM(K) = FROM(J400)
            COST(K) = COST(J400)
            CAPAC(K) = CAPAC(J400)
            TFLOOR(K) = TFLOOR(J400)
            NAME(K) = NAME(J400)
          endif
        else IF (IABS(J).NE.L) then
          K = K+1
          IF (K.NE.J400) then
            FROM(K) = FROM(J400)
            COST(K) = COST(J400)
            CAPAC(K) = CAPAC(J400)
            TFLOOR(K) = TFLOOR(J400)
            NAME(K) = NAME(J400)
          endif
        endif
      end do

      MARC = K
      MREG = K

C      NSTOP = 15
C      IF (MARC+MAX0(1,MSORC)+1.GT.LARC) CALL ERRMES(NSTOP)
      IF (MARC+MAX0(1,MSORC)+1.GT.LARC) then
        errorflag = 16     !Too many arcs for this version
        messerror = '  Too many arcs for this version'
        return
      endif

C     ADD REGULAR SLACKS

      I = -FROM(MARC)
      THD = NEXT(MNODP1)
      NEXT(MNODP1) = MNODP1

      IF (THD.EQ.MNODP1) then
C       NO REGULAR SLACKS, ADD DUMMY
        MARC = MARC+1
        FROM(MARC) = ISIGN(MNODP1,I)
        COST(MARC) = 0
        CAPAC(MARC) = -BIG
        TFLOOR(MARC) = 0
        NAME(MARC) = DUMMY
      endif

C     FOLLOW LIST
      do while (THD.NE.MNODP1)
        MARC = MARC+1
        NAME(MARC) = SLACK
        FROM(MARC) = ISIGN(THD,I)
        COST(MARC) = 0
        CAPAC(MARC) = LEVEL(THD)
        LEVEL(THD) = 0
        TFLOOR(MARC) = 0
        NXT = NEXT(THD)
        NEXT(THD) = 0
        THD = NXT
      end do

      MSLK = MARC

C     ADD EXCESS ARC AT END OF REGULAR SLACKS
      MARC = MARC+1
      FROM(MARC) = ISIGN(MNODP2,-I)
      COST(MARC) = BIG
      CAPAC(MARC) = 0
      TFLOOR(MARC) = 0
      NAME(MARC) = XCESS

C     INITIALIZE ARTIFICIAL
      ZERO = MARC+1
      FROM(ZERO) = MNODP1
      COST(ZERO) = BIG
      CAPAC(ZERO) = 0
      TFLOOR(ZERO) = 0
      NAME(ZERO) = ARTIF

      RETURN
      END
      
      
C     ------------------------------------------------------------------
C                             LWBND
C     ------------------------------------------------------------------

C     SUBROUTINE LWBND ADJUSTS THE NETWORK FLOW PROBLEM FOR LOWER BOUNDS
C     ON THE ARC FLOWS

      SUBROUTINE LWBND(MNODE,MARC,MNODP1,MTREE)

      implicit none

      integer*4 BIG
      integer*4 LNODE,MNODE,DOWN(2000),NEXT(2000),LEVEL(2000)
      integer*4 ARCID(2000),FLOW(2000),DUAL(2000),CAT(2000)
      integer*4 NET,I500,DWN,MTREE,THD,MNODP1,J,NXT,MARC,LARCP1
      integer*4 errorflag ! transfres error codes among sbr's

      CHARACTER(80) messerror !for error messages

      COMMON/NODINF/DOWN,NEXT,LEVEL,ARCID,FLOW,DUAL,CAT
      COMMON/PARAMI/LNODE,LARCP1,BIG
      common /errs/ messerror, errorflag ! transfres error codes among sbr's

	SAVE

C     LOCATE SOURCES AND SINKS FOR PROBLEM
C     ADJUSTED FOR LOWER BOUNDS

      NET=0
      THD = MNODP1
      DO I500 = 1,MNODE
        J = FLOW(I500)
        NET = NET+J

        if (J .NE. 0) then
          if (J < 0) then
C         SINK
            FLOW(I500) = -J
C           LINK DEMANDS IN DECREASING SIZE ORDER
            DWN = MNODP1
  200       CONTINUE
            NXT = DOWN(DWN)
            IF (FLOW(NXT)+J.GT.0) then
              DWN = NXT
              GOTO 200
            endif
            DOWN(DWN) = I500
            DOWN(I500) = NXT
            LEVEL(I500) = -1
          else if (J > 0) then
C           SOURCE
            MTREE = MTREE+1
            ARCID(I500) = -MARC
            FLOW(I500) = J
            NEXT(THD) = I500
            DOWN(I500) = MNODP1
            NEXT(I500) = MNODP1
            LEVEL(I500) = 1
            DUAL(I500) = BIG
            THD = I500
          endif
        endif
      end do

C     CHECK FOR FEASIBILITY
C      NSTOP = 16
C      IF (NET.LT.0) CALL ERRMES(NSTOP)
      IF (NET.LT.0) then
        errorflag = 17    ! infeasible - lower bounds on arc flow cannot be met
        messerror = '  Infeasible - lower bounds on arc flow
     + cannot be met'
        return
      endif
      
      RETURN
      END


C     ------------------------------------------------------------------
C                              START
C     ------------------------------------------------------------------
      SUBROUTINE START(MNODE,MARC,LNODP1,MNODP1,MTREE,ZERO)

      implicit none

C     SUBROUTINE START PROVIDES AN INITIAL BASIS USING ALG F.1 OF
C     KENNINGTON AND HELGASON P.245-248

      integer*4 TO,TOO,THD,TXOR,TRY,PRICE,NEWARC,NEWPR,ZERO,BIG
     +       ,FRM,LVJ,FM,LST,FLW
      integer*4 LNODE,MNODE,DOWN(2000),NEXT(2000)
     +       ,LEVEL(2000),ARCID(2000),FLOW(2000),DUAL(2000),CAT(2000)
      integer*4 FROM(2600000),COST(2600000),CAPAC(2600000)
      integer*4 I,J,K,L,I2500,MNODP1,MTREE,M
      integer*4 MARC,LNODP1,LARCP1,TFLOOR(2600000)

      COMMON/NODINF/DOWN,NEXT,LEVEL,ARCID,FLOW,DUAL,CAT
      COMMON/ARCINI/FROM,COST,CAPAC,TFLOOR
      COMMON/PARAMI/LNODE,LARCP1,BIG

	SAVE

      EXTERNAL TXOR

C     ADVANCED START
C
C     SELECT HIGHEST RANK DEMAND ON LIST

  100 CONTINUE

      TO = DOWN(MNODP1)

C     IS LIST EXHAUSTED?
      IF (TO.ne.MNODP1) then

  200   CONTINUE

C       SET TO LINK TO ARTIFICIAL
        NEWARC = ZERO
        NEWPR = BIG

C       ANY DEMAND LEFT
        IF (FLOW(TO).EQ.0) GOTO 1300

C       LOOK FOR SOURCES FIRST

        TRY = CAT(TO)
        FRM = FROM(TRY)
        LST = ISIGN(LNODP1,FRM)

  300   CONTINUE

C		IS IT UNAVAILABLE?
        IF (CAPAC(TRY).GT.0) then
          FM = IABS(FRM)
C           IS IT FROM A NON-SOURCE?
          IF ((LEVEL(FM).eq.1) .or. (ARCID(FM).ne.0)) then
            PRICE = COST(TRY)
C             IS COST WORSE?
            IF (PRICE.lt.NEWPR) then
C               DOES CAPACITY EXCEED DEMAND?
              IF (CAPAC(TRY).gt.FLOW(TO)) then
C               IS THERE NOT ENOUGH SUPPLY FOR DEMAND?
                IF (FLOW(FM).ge.FLOW(TO)) then
                  NEWARC = TRY
                  NEWPR = PRICE
                endif
C               IS THERE NOT ENOUGH SUPPLY FOR CAPACITY?
              elseif (FLOW(FM).ge.CAPAC(TRY)) then
                NEWARC = -TRY
                NEWPR = PRICE
              endif
            endif
          endif
        endif

        IF (NEWPR.ne.0) then
          TRY = TRY+1
          FRM = FROM(TRY)
          IF (TXOR(FRM,LST).GT.0) GOTO 300
          IF (NEWARC.EQ.ZERO) GOTO 800
        endif

C       ARE WE SENDING DEMAND?
        IF (NEWARC.le.0) then
C           SEND CAPACITY
          NEWARC = -NEWARC
          FM = IABS(FROM(NEWARC))
C         GET CAPACITY
          FLW = CAPAC(NEWARC)
C         MARK UNAVAILABLE
          CAPAC(NEWARC) = -FLW
C         ADJUST FLOWS
          FLOW(FM) = FLOW(FM)-FLW
          FLOW(TO) = FLOW(TO)-FLW
          GOTO 200
C         SEND DEMAND
        end if

C       MARK UNAVAILABLE
        CAPAC(NEWARC) = -CAPAC(NEWARC)

C       ADJUST FLOWS
        FM = IABS(FROM(NEWARC))
        FLOW(FM) = FLOW(FM)-FLOW(TO)
        K = BIG

        GOTO 1400

C       LOOK FOR TRANSHIPMENT POINTS

  800   CONTINUE

        TRY = CAT(TO)
        FRM = FROM(TRY)

  900   CONTINUE !might be suitable for do...while loop

C       IS IT UNAVAILABLE?
        IF (CAPAC(TRY).GT.0) then
          FM = IABS(FRM)
C         IS IT ALREADY LINKED
          IF (LEVEL(FM).EQ.0) then
            PRICE = COST(TRY)
C           IS COST WORSE
            IF (PRICE.LT.NEWPR) then
              NEWARC = TRY
              NEWPR = PRICE
              IF (NEWPR.EQ.0) GOTO 1100
            endif
          endif
        endif

        TRY = TRY + 1
        FRM = FROM(TRY)
        IF (TXOR(FRM,LST).GT.0) GOTO 900
        IF (NEWARC.EQ.ZERO) GOTO 1300

 1100   CONTINUE

        FM = IABS(FROM(NEWARC))

C       DOES CAPACITY EXCEED DEMAND?
        IF (CAPAC(NEWARC).LE.FLOW(TO)) then
C         GET CAPACITY
          FLW = CAPAC(NEWARC)
C         MARK UNAVAILABLE
          CAPAC(NEWARC) = -FLW
C         ADJUST FLOWS
          FLOW(FM) = FLW
          FLOW(TO) = FLOW(TO)-FLW
C         LINK IN TO DEMAND LIST
          DOWN(FM) = TO
          DOWN(MNODP1) = FM
C         START NEW CHAIN
          LEVEL(FM) = -1
          GOTO 100
        endif

C       MARK UNAVAILABLE
        CAPAC(NEWARC) = -CAPAC(NEWARC)
C       PASS ALONG FLOW DEMAND
        FLOW(FM) = FLOW(TO)
C       LINK IN NEW DEMAND NODE
        DOWN(FM) = DOWN(TO)
        DOWN(TO) = FM
        DOWN(MNODP1) = FM
        NEXT(FM) = TO
        ARCID(TO) = NEWARC
        LEVEL(FM) = LEVEL(TO)-1
        DUAL(TO) = NEWPR
        MTREE = MTREE+1
        GOTO 100

C       ADD CHAIN TO TREE
 1300   CONTINUE

        K = 0

 1400   CONTINUE

C       REMOVE FROM DEMAND LIST
        DOWN(MNODP1) = DOWN(TO)
C       LINK IN AS LEFTMOST BRANCH
        FM = IABS(FROM(NEWARC))
        ARCID(TO) = NEWARC
        DUAL(TO) = NEWPR
        DOWN(TO) = FM
        I = NEXT(FM)
        NEXT(FM) = TO
        J = LEVEL(FM)-LEVEL(TO)+1
        THD = FM

 1500   CONTINUE

C       MOVE ALONG CHAIN
        THD = NEXT(THD)
C       ADJUST LEVEL AND DUAL VARIABLES
        L = LEVEL(THD)
        LEVEL(THD) = L+J
        K = K-DUAL(THD)
        DUAL(THD) = K

        IF (L.NE.-1) GOTO 1500

        NEXT(THD) = I
        MTREE = MTREE+1

        GOTO 100
      endif	!1600

C     SET UP TO EXPAND TREE
      TO = 1
      TRY = 1
      FRM = FROM(TRY)

 1700 CONTINUE

C     DO WE NEED TO EXPAND TREE TO REACH ALL NODES?
      IF (MTREE.ne.MNODE) then

        TOO = TO
        NEWPR = BIG
C		SEARCH FOR LEAST COST CONNECTABLE ARC IN CURRENT TO-GROUP

 1800   CONTINUE

        LVJ = LEVEL(TO)
        LST = ISIGN(LNODP1,FRM)

 1900 CONTINUE

      IF (CAPAC(TRY).gt.0) then
        M = COST(TRY)
        IF (NEWPR.ge.M) then
          FM = IABS(FRM)
          IF (LEVEL(FM).ne.0) then
            IF (LVJ.eq.0) then
C             TO END IS CONNECTABLE
              I = FM
              J = TO
              K = M
              L = TRY
              NEWPR = M
            endif
          else
            IF (LVJ.ne.0) then
C             FROM END IS CONNECTABLE
              I = TO
              J = FM
              K = -M
              L = -TRY
            endif
          endif
        endif
      endif

      TRY = TRY+1
      FRM = FROM(TRY)

C     IS TO-GROUP EXHAUSTED?
      IF (TXOR(FRM,LST).GT.0) GOTO 1900

C     PREPARE FOR NEXT TO-GROUP
      TO = TO+1

      IF (TO.eq.MNODP1) then
        TO = 1
        TRY = 1
        FRM = FROM(TRY)
      endif

      IF (NEWPR.eq.BIG) then
        IF (TO.NE.TOO) GOTO 1800
C         NOT ALL NODES CONNECTABLE - CHECK FOR ISOLATED POINTS
          DO I2500 = 1,MNODE
            IF (LEVEL(I2500).eq.0) then
C             TEST FOR ARCS RUNNING FROM IT
C             IF (ARCID(I2500).ne.-1) then
C               CHECK FOR DUMMY - NO ARCS RUNNING TO IT
C               J2500 = CAT(I2500)
C             endif
C             ADD ARTIFICIAL
              ARCID(I2500) = 0
              FLOW(I2500) = 0
              NEXT(I2500) = NEXT(MNODP1)
              NEXT(MNODP1) = I2500
              DOWN(I2500) = MNODP1
              LEVEL(I2500) = 1
              DUAL(I2500) = -BIG
            endif
          end do
        else
          ARCID(J) = L
          DOWN(J) = I
          NEXT(J) = NEXT(I)
          NEXT(I) = J
          LEVEL(J) = LEVEL(I)+1
          DUAL(J) = DUAL(I)-K
          NEWARC = IABS(L)
          CAPAC(NEWARC) = -CAPAC(NEWARC)
          MTREE = MTREE+1
          GOTO 1700
        endif
      endif !2700

C     CLEAR UPPER BOUND FLAGS ON BASIC ARCS
      DO I = 1,MNODE
        J = IABS(ARCID(I))
        CAPAC(J) = -CAPAC(J)
      end do

C     CLEAR OUT UPPER BOUND FLAGS ON DUMMY ARCS
      DO I = 1,MARC
        IF (CAPAC(I)+BIG.EQ.0) then
          CAPAC(I) = 0
        endif
      end do

C     SET UPPER BOUND FOR ARTIFICIAL AND EXCESS
      CAPAC(ZERO) = BIG
      CAPAC(MARC) = BIG

      RETURN
      END


C     ------------------------------------------------------------------
C                                    OUT
C     ------------------------------------------------------------------

C     SUBROUTINE OUT READJUSTS THE SOLUTION FOR LOWER BOUNDS ON ARC FLOWS
C     AND OUTPUTS RESULTS

      SUBROUTINE OUT(MNODE,ZERO,LNODP1,MNODP1,MSLK,KOST0,
     +OPTIM,IFLAG2,KOST,IRANKD,X,Y,NUMNOD)

      implicit none

      LOGICAL INFEAS,OPTIM,PRES

      integer*4 STFROM(2000),STTO(2000),STFLOW(2000),FROM(2600000)
      integer*4 COST(2600000),CAPAC(2600000),TFLOOR(2600000)
      integer*4 TXOR,TO,TRY,BIG,FRM,FM,LST,FLW,ZERO,IFLAG2
      integer*4 SUMFLW,KOST,IRANKD(2000),IDUM,NUMNOD
      integer*4 UPPER(2600000),IAMADF,IAMADT,KOST0
      integer*4 LNODE,MNODE,DOWN(2000),NEXT(2000)
      integer*4 LEVEL(2000),ARCID(2000),FLOW(2000),DUAL(2000),CAT(2000)
      integer*4 I,J,MSLK,LARCP1,LNODP1
      integer*4 I13,J15,J16,J1700,L1700,MNODP1,I1700

      double precision X(2000),Y(2000)

      COMMON/PARAMI/LNODE,LARCP1,BIG
      COMMON/NODINF/DOWN,NEXT,LEVEL,ARCID,FLOW,DUAL,CAT
      COMMON/ARCINI/FROM,COST,CAPAC,TFLOOR
      common/FLOWS/STFROM,STTO,STFLOW

	SAVE

      EXTERNAL TXOR

C     Initialise STTO(), STFROM() & STFLOW() to ensure overwriting of previous
C     values, even when flow is zero.  Use maximum array size for this.
      DO I = 1,2000
           STFROM(I) = 0
           STTO(I) = 0
           STFLOW(I) = 0
      end do

C     Initialise IDUM
      IDUM = 0

C     Now this is the start of the subroutine proper

      INFEAS = .FALSE.
      PRES = .FALSE.
      KOST = KOST0

      DO I = 1,MNODE
        J = IABS(ARCID(I))
        IF ((FLOW(I).NE.0).AND.(COST(J).EQ.BIG)) then
          INFEAS = .TRUE.
        endif
        KOST = KOST+COST(J)*FLOW(I)
      end do

      DO I = 1,MSLK
        IF (CAPAC(I) < 0) then
          J = -CAPAC(I)
          KOST = KOST+COST(I)*J
        endif
      end do

      if (iflag2 .eq. 0) then
         IF (OPTIM.AND.INFEAS) WRITE(12,5100)
      end if

      DO I = 1,LARCP1
        UPPER(I) = TFLOOR(I)+IABS(CAPAC(I))
      end do

      PRES = .TRUE.

C     MOVE BASIC FLOWS TO CAPACITY CELLS FOR BASIC ARCS
      IF (PRES) then
        DO I = 1,MNODE
          J = IABS(ARCID(I))
          CAPAC(J) = -FLOW(I)
        end do
      endif

      TO = 1
      TRY = 1
      FRM = FROM(TRY)

  600 CONTINUE

      LST = ISIGN(LNODP1,FRM)

  700 CONTINUE

      FLW = MAX0(0,-CAPAC(TRY))+TFLOOR(TRY)
      FM = IABS(FRM)
C     CST = FLW*COST(TRY)
      IF (FM.GT.MNODE) GOTO 750
      IF (TO.GT.MNODE) GOTO 730
C     next two statements could be commented out b/c they are redundant?
      IF (UPPER(TRY).LE.0) GOTO 800
      if (iflag2.NE.0) goto 800

      GOTO 800

  730 IF (FLW.GT.0) then
         if (iflag2.NE.0) then
            goto 800
         endif
      endif

      GOTO 800

  750 IF (TO.GT.MNODE) GOTO 770
      IF (FLW.GT.0) then
         if (iflag2.NE.0) then
             goto 800
         endif
      endif

      GOTO 800

  770 IF (FLW.GT.0) then
         if (iflag2.NE.0) then
            goto 800
         endif
      endif

  800 CONTINUE

      TRY = TRY+1
      IF (TRY.GT.ZERO) GOTO 950
      FRM = FROM(TRY)
      if (frm.ge.100.and.lst.ge.100) goto 700
      IF (TXOR(FRM,LST).GT.0) GOTO 700
      TO = TO+1

      GOTO 600

  950 CONTINUE

C     SUM FLOWS IN EACH TO-GROUP
      TO = 1
      TRY = 1
      FRM = FROM(TRY)

  900 CONTINUE

      LST = ISIGN(LNODP1,FRM)

C     CLEAR OUT ALL FLOW AND DOWN CELLS
C     FOR SUMS OF FLOWS AND COSTS RESPECTIVELY
      DO I = 1,MNODE
        FLOW(I) = 0
        DOWN(I) = 0
      end do

      SUMFLW = 0

 1100 CONTINUE

      FLW = MAX0(0,-CAPAC(TRY))+TFLOOR(TRY)

      IF (FLW.NE.0) then
        FM = IABS(FRM)
        FLOW(FM) = FLOW(FM)+FLW
        DOWN(FM) = DOWN(FM)+FLW*COST(TRY)
        SUMFLW = SUMFLW+FLW
      endif

      TRY = TRY+1
      FRM = FROM(TRY)

      if (frm.ge.100.and.lst.ge.100) goto 1100

      IF (TXOR(FRM,LST).GT.0) GOTO 1100

C     OUTPUT NON-ZERO VALUES
      IF (SUMFLW.NE.0) then
        DO I13 = 1,MNODE
          IF (FLOW(I13).NE.0) then
            IDUM = IDUM + 1
            IF (IFLAG2.EQ.0) THEN
              IAMADF = IRANKD(I13)
              IAMADT = IRANKD(TO)
            ELSE
              IAMADF = I13
              IAMADT = TO
            ENDIF
            STFROM(IDUM) = IAMADF
            STTO(IDUM) = IAMADT
            STFLOW(IDUM) = FLOW(I13)
            if (iflag2.eq.0) then
              WRITE(9,9999)X(IAMADF),Y(IAMADF)
     +             ,X(IAMADT),Y(IAMADT),FLOW(I13)
            endif
          endif
        end do
      endif

 9999 FORMAT(1X ,4(F16.7,2X),I12)

      TO = TO+1

      IF (TO.NE.MNODP1) GOTO 900

C     NORMALISE NODE POTENTIALS SO THAT MINIMUM ABSOLUTE VALUE = 0
      IF (.NOT.(OPTIM.AND.INFEAS)) then
        I = 1

        DO J15 = 2,MNODE
          IF (IABS(DUAL(J15)).LT.IABS(DUAL(I))) I = J15
        end do

        DO J16 = 1,MNODE
          IF (J16.NE.I) DUAL(J16) = DUAL(J16)-DUAL(I)
        end do

        DUAL(I) = 0
        J1700 = 1+(MNODE-1)/8
        L1700 = 0

        DO I1700 = 1,J1700
C         K1700 = L1700+1
          L1700 = MIN0(MNODE,L1700+8)
        end do

        NUMNOD = IDUM

C        WRITE(9,*)

        IF (IFLAG2.EQ.0) THEN
          write(9,*)
          WRITE(9,5095) IDUM
        ENDIF
      endif
 5095 format (1X , 'The number of nodal flows output above = ' , I8)
 
 5100 FORMAT(//' ******** INFEASIBLE ********'///
     +       ' ******** FLOW PATTERN OBTAINED AT TERMINATION *******'/)

      Return
      END

C     ------------------------------------------------------------------
C                               KHFLOW
C     ------------------------------------------------------------------
C     ZERO LOWER BOUNDS ON ARC FLOWS, USING THE PREDECESSOR/THREAD/DISTANCE
C     DATA STRUCTURE DESCRIBED IN KENNINGTON AND HELGASON, APPENDIX B,
C     P. 207-221.

      SUBROUTINE KHFLOW(LNODP1,MNODP2,ITER,OPTIM)

      implicit none

      LOGICAL OPTIM

      integer*4 TXOR,TO,PRICE0,TOO,I1300,LSAVE,L,BIG
      integer*4 TRY,PRICE,NEWARC,NEWPR,NEWFRM,NEWTO,LNODP1
      integer*4 THD,DW(2),CH(2),DWN,CHG,THETA,JTHETA,KTHETA,POSS,JPOSS
      integer*4 FRM,FM,LST,DWE,FLW,AID,Q1,Q2,DIR,REF,LARCP1
      integer*4 LNODE,DOWN(2000),NEXT(2000),ITER,MNODP2
      integer*4 LEVEL(2000),ARCID(2000),FLOW(2000),DUAL(2000),CAT(2000)
      integer*4 LDIFF,I,I800,K800,L1100,ITHETA,J,K,NXT,LSTAR,N
      integer*4 FROM(2600000),COST(2600000),CAPAC(2600000)
      integer*4 TFLOOR(2600000)

      COMMON/PARAMI/LNODE,LARCP1,BIG
      COMMON/NODINF/DOWN,NEXT,LEVEL,ARCID,FLOW,DUAL,CAT
      COMMON/ARCINI/FROM,COST,CAPAC,TFLOOR

	SAVE

      EXTERNAL TXOR

C     INITIALIZE PRICING
      TO = 1
      TRY = 1
      FRM = FROM(TRY)
      ITER = 0
      OPTIM = .FALSE.

C     NEW ITERATION
  100 CONTINUE

      ITER = ITER+1
C     PRICING
C     NOTE THAT WE ARE PRICING OUT BASIC ARCS
      TOO = TO
      NEWPR = 0

  200 CONTINUE

      PRICE0 = -DUAL(TO)
      LST = ISIGN(LNODP1,FRM)

  250 CONTINUE

      FM = IABS(FRM)
      PRICE = DUAL(FM)+PRICE0-COST(TRY)

      if (CAPAC(TRY) < 0) then
        PRICE = -PRICE
        IF (PRICE > NEWPR) then
          NEWARC = TRY
          NEWPR = PRICE
          NEWTO = TO
        endif
      elseif (CAPAC(TRY) > 0) then
        IF (PRICE > NEWPR) then
          NEWARC = TRY
          NEWPR = PRICE
          NEWTO = TO
        endif
      endif

      TRY = TRY+1
      FRM = FROM(TRY)

      IF (TXOR(FRM,LST).GT.0) GOTO 250

      TO = TO+1

      IF (TO .eq. MNODP2) then
        TO = 1
        TRY = 1
        FRM = FROM(TRY)
      endif

      IF (NEWPR.EQ.0) then
        IF (TO.NE.TOO) GOTO 200
C       OPTIMAL INDICATION
          OPTIM = .TRUE.
        RETURN
      endif

      NEWFRM = IABS(FROM(NEWARC))
C     RATIO TEST
      THETA = IABS(CAPAC(NEWARC))
      JTHETA = 0
C     SET FOR CYCLE SEARCH
      CH(2) = ISIGN(LARCP1,CAPAC(NEWARC))
      CH(1) = -CH(2)
      DW(1) = NEWFRM
      DW(2) = NEWTO
      LDIFF = LEVEL(NEWFRM)-LEVEL(NEWTO)
      KTHETA = 1

      if (LDIFF .NE. 0) then
        if (LDIFF < 0) then
          KTHETA = 2
        endif
        DWN = DW(KTHETA)
        CHG = CH(KTHETA)
        K800 = IABS(LDIFF)

        DO I800 = 1,K800
          if (chg.lt.-500.and.arcid(dwn).gt.500) goto 1234
          if (chg.gt.500.and.arcid(dwn).lt.-500) goto 1234
          if (chg.gt.500.and.arcid(dwn).gt.500) goto 650
          if (chg.lt.-500.and.arcid(dwn).lt.-500) goto 650
          IF (TXOR(CHG,ARCID(DWN)).GT.0) GOTO 650

 1234     continue

C         INCREASING FLOW
          I = IABS(ARCID(DWN))
          POSS = CAPAC(I)-FLOW(DWN)
          JPOSS = -DWN

          GOTO 700

C         DECREASING FLOW
  650     POSS = FLOW(DWN)
          JPOSS = DWN
C         FIND MIN

  700     CONTINUE

          IF (THETA.GT.POSS) then
            THETA = POSS
            JTHETA = JPOSS
            IF (THETA.EQ.0) GOTO 1200
          endif

          DWN = DOWN(DWN)

        end do

        DW(KTHETA) = DWN

C       AT COMMON LEVEL
      endif	!if (LDIFF .NE. 0) then

C     SEARCH FOR CYCLE END
  900 CONTINUE

      IF (DW(1).EQ.DW(2)) GOTO 1150

      DO L1100 = 1,2
        DWN = DW(L1100)
        if (ch(l1100).lt.-500.and.arcid(dwn).gt.500) goto 1235
        if (ch(l1100).gt.500.and.arcid(dwn).lt.-500) goto 1235
        if (ch(l1100).gt.500.and.arcid(dwn).gt.500) goto 950
        if (ch(l1100).lt.-500.and.arcid(dwn).lt.-500) goto 950
        IF (TXOR(CH(L1100),ARCID(DWN)).GT.0) GOTO 950

 1235   continue

C       INCREASING FLOW
        I = IABS(ARCID(DWN))
        POSS = CAPAC(I)-FLOW(DWN)
        JPOSS = -DWN

        GOTO 1000

  950   POSS = FLOW(DWN)
        JPOSS = DWN

C		FIND MIN
 1000   CONTINUE

        IF (THETA.GT.POSS) then

          THETA = POSS
          JTHETA = JPOSS
          KTHETA = L1100

          IF (THETA.EQ.0) GOTO 1200

        endif

        DW(L1100) = DOWN(DWN)

      end do

      GOTO 900

 1150 DWE = DW(1)

C     RATIO TEST COMPLETE

 1200 CONTINUE

      IF (THETA.NE.0) then
C     UPDATE FLOWS ON INTACT PORTION OF CYCLE
        DW(1) = NEWFRM
        DW(2) = NEWTO
        IF (JTHETA.NE.0) DW(KTHETA) = IABS(JTHETA)

        DO I1300 = 1,2
          DWN = DW(I1300)
          CHG = ISIGN(THETA,CH(I1300))

 1250     CONTINUE

          IF (DWN.NE.DWE) then
            FLOW(DWN) = FLOW(DWN)-CHG*ISIGN(1,ARCID(DWN))
            DWN = DOWN(DWN)
            GOTO 1250
          endif
        end do
      endif

      IF (JTHETA.EQ.0) then
C     CHANGE OF BOUNDS ONLY
        CAPAC(NEWARC) = -CAPAC(NEWARC)
        GOTO 100
      endif

      ITHETA = IABS(JTHETA)

      IF (JTHETA.LE.0) then
        J = IABS(ARCID(ITHETA))
C       SET OLD ARC TO UPPER BOUND
        CAPAC(J) = -CAPAC(J)
      endif

C     FLOW ON NEWARC
      FLW = THETA

      IF (CAPAC(NEWARC).LE.0) then
        CAPAC(NEWARC) = -CAPAC(NEWARC)
        FLW = CAPAC(NEWARC)-FLW
        NEWPR = -NEWPR
      endif

      IF (KTHETA.NE.2) then
        Q1 = NEWFRM
        Q2 = NEWTO
        AID = -NEWARC
        NEWPR = -NEWPR
      else
        Q1 = NEWTO
        Q2 = NEWFRM
        AID = NEWARC
      endif

C     UPDATE TREE
      I = Q1
      J = DOWN(I)
      LSTAR = LEVEL(Q2)+1

      IF (THETA.EQ.0) GOTO 1900

C     FLOWS NEED TO BE UPDATED
      CHG = ISIGN(THETA,CH(KTHETA))

 1650 CONTINUE

C     UPDATE DUAL VARIABLE ON STEM
      DUAL(I) = DUAL(I)+NEWPR
C     UPDATE FLOW ON STEM
      N = FLOW(I)
      FLOW(I) = FLW
C     UPDATE ARC ID ON STEM
      DIR = ISIGN(1,ARCID(I))
      REF = IABS(ARCID(I))
      ARCID(I) = AID
C     PREPARE FOR LEVEL UPDATES
      LSAVE = LEVEL(I)
      LDIFF = LSTAR-LSAVE
      LEVEL(I) = LSTAR
      THD = I

 1700 CONTINUE

      NXT = NEXT(THD)

      IF (LEVEL(NXT).Gt.LSAVE) then
C       UPDATE LEVEL
        LEVEL(NXT) = LEVEL(NXT)+LDIFF
C       UPDATE DUAL VARIABLE
        DUAL(NXT) = DUAL(NXT)+NEWPR
        THD = NXT
        GOTO 1700
      endif

      K=J

 1800 CONTINUE

      L = NEXT(K)

      IF (L.NE.I) then
        K = L
        GOTO 1800
      endif

C     TEST FOR LEAVING ARC
      IF (I.EQ.ITHETA) GOTO 2200

C       PREPARE FOR NEXT UPDATE ON STEM
        FLW = N-CHG*DIR
        AID = -ISIGN(REF,DIR)
C       MOVE DOWN STEM
        NEXT(K) = NXT
        NEXT(THD) = J
        K = I
        I = J
        J = DOWN(J)
        DOWN(I) = K
        LSTAR = LSTAR+1

        GOTO 1650

 1900   CONTINUE

C       ONLY ON NEW ARC CHANGES
 1950   CONTINUE

C       UPDATE DUAL VARIABLE ON STEM
        DUAL(I) = DUAL(I)+NEWPR
C       UPDATE FLOW ON STEM
        N = FLOW(I)
        FLOW(I) = FLW
C       UPDATE ARC ID ON STEM
        DIR = ISIGN(1,ARCID(I))
        REF = IABS(ARCID(I))
        ARCID(I) = AID
C       PREPARE FOR LEVEL UPDATES
        LSAVE = LEVEL(I)
        LDIFF = LSTAR-LSAVE
        LEVEL(I) = LSTAR
        THD = I

 2000   CONTINUE

        NXT = NEXT(THD)

        IF (LEVEL(NXT).GT.LSAVE) then
C         UPDATE LEVEL
          LEVEL(NXT) = LEVEL(NXT)+LDIFF
C         UPDATE DUAL VARIABLE
          DUAL(NXT) = DUAL(NXT)+NEWPR
          THD = NXT
          GOTO 2000
        endif

        K = J

 2100   CONTINUE

        L = NEXT(K)

        IF (L.NE.I) then
        K = L
        GOTO 2100
      endif

C     TEST FOR LEAVING ARC
      IF (I.EQ.ITHETA) goto 2200

C     PREPARE FOR NEXT UPDATE ON STEM
      FLW = N
      AID = -ISIGN(REF,DIR)
C     MOVE DOWN STEM
      NEXT(K) = NXT
      NEXT(THD) = J
      K = I
      I = J
      J = DOWN(J)
      DOWN(I) = K
      LSTAR = LSTAR+1
      GOTO 1950

 2200 CONTINUE

      NEXT(K) = NXT
      NEXT(THD) = NEXT(Q2)
      NEXT(Q2) = Q1
      DOWN(Q1) = Q2

      GOTO 100

      END


C     ------------------------------------------------------------------
C                              INDAT
C     ------------------------------------------------------------------

C     SUBROUTINE INDAT PERFORMS INPUT AND VALIDATION OF ALL DATA
C     IT ALSO PERFORMS A SIMPLE FEASIBILITY TEST AND VARIOUS
C     INITIALISATION TASKS

      SUBROUTINE INDAT(MNODE,MARC,MSORC,LARC,LNODP1,MNODP1,MNODP2,KOST0)

      implicit none

      integer*4 TXOR,ISUPLY(2000),ITOTSY,INFROM(2600000),INTO(2600000)
      integer*4 INCOST(2600000),IMNODE,INCAT(2000),INMARC
      integer*4 LNODE,MNODE,DOWN(2000),NEXT(2000),MSORC,M,L,KOST0
      integer*4 LEVEL(2000),ARCID(2000),FLOW(2000),DUAL(2000),CAT(2000)
      integer*4 BIG,I,J,K,II,JJ,J300,K200,I200,N,K1000,M1000
      integer*4 L1000,LARCP1,LNODP1,MNODP1,MNODP2,NET
      integer*4 LARC,MARC,FROM(2600000),COST(2600000),CAPAC(2600000)
      integer*4 TFLOOR(2600000),KK,J1000,J1100
      integer*4 errorflag ! transfres error codes among sbr's

      CHARACTER(8) NAME,TAG,DUMMY
      CHARACTER(80) messerror !for error messages

      COMMON/PARAMI/LNODE,LARCP1,BIG
      COMMON/NODINF/DOWN,NEXT,LEVEL,ARCID,FLOW,DUAL,CAT
      COMMON/ARCINI/FROM,COST,CAPAC,TFLOOR
      COMMON/ARCINR/NAME(2600000)
      COMMON/LINK/IMNODE,INCAT,ITOTSY,ISUPLY,INFROM,INTO,INCOST,INMARC
      common/errs/messerror,errorflag ! transfres error codes among sbr's

	SAVE

      EXTERNAL TXOR

      TAG=' '
      DUMMY='DUMMY'

      MNODE = IMNODE

      DO I = 1,MNODE
        CAT(I) = INCAT(I)
      end do

C     ALLOW FOR ARTIFICIAL ARC
      LARC = LARCP1-1

C     UNUSED NODE NUMBER
      LNODP1 = LNODE+1

C     Now get information on number of nodes and number of arcs for
C     each node
      MNODP1 = MNODE+1
      MNODP2 = MNODE+2

C      NSTOP = 4
C      IF (MNODP1.GT.LNODE) CALL ERRMES(NSTOP)
      IF (MNODP1.GT.LNODE) then
        errorflag = 18       !  Too many nodes for this version
        messerror = '  Too many nodes for this version'
        return
      endif

C     INITIALIZE NODE ARRAYS
      DO J = 1,MNODP1
        DOWN(J) = 0
        NEXT(J) = 0
        LEVEL(J) = 0
        ARCID(J) = 0
        FLOW(J) = 0
        DUAL(J) = 0
      end do

      MARC = 0
C
C     RESERVE LOCATIONS FOR INPUT ARCS BY FILLING WITH DUMMIES
C     CAT(N) WILL POINT TO THE NEXT OPEN LOCATION FOR STORING
C     ARCS WHOSE TERMINAL NODE IS N.

C      NSTOP = 10

      I = 1
      J = 0

      DO J300 = 1,MNODE
        I = -I
        K200 = MAX0(1,CAT(J300))

C       IF (J+K200.GT.LARC) CALL ERRMES(NSTOP)
        IF (J+K200.GT.LARC) then
          errorflag = 19         ! Too many arcs for this version
          messerror = '  Too many arcs for this version'
          return
        endif

        CAT(J300) = ISIGN(J+1,I)

        DO I200 = 1,K200
          J = J+1
C         NOTE THAT DUMMY ARCS HAVE SAME INITIAL AND TERMINAL NODE
          FROM(J) = ISIGN(J300,I)
          COST(J) = 0
          CAPAC(J) = -BIG
          TFLOOR(J) = 0
          NAME(J) = DUMMY
        end do
      end do

      MARC = J+1

C     IF (MARC.GT.LARC) CALL ERRMES(NSTOP)
      IF (MARC.GT.LARC) then
        errorflag = 19         !
        messerror = '  Too many arcs for this version'
        return
      endif

      FROM(MARC) = ISIGN(MNODP1,-I)
C
C     DURING SETUP WE USE THE ABSENCE OF THE UPPER BOUND FLAG
C     ON AN ARC TO INDICATE THAT THE ARC IS ELIGIBLE TO BECOME
C     PART OF THE STARTING BASIS, OTHERWISE IT IS NOT
C
      NET = 0
      MSORC = 0
C     INITIALIZE SUPPLY AND DEMAND LISTS WITH SELF-POINTER
      NEXT(MNODP1) = MNODP1
      DOWN(MNODP1) = MNODP1

      DO I = 1,MNODE
        J=ISUPLY(I)

C        NSTOP = 7
C        IF (FLOW(I).NE.0) CALL ERRMES(NSTOP)
        IF (FLOW(I).NE.0) then
          errorflag = 20     ! Supply/demand already specified for node
          messerror = '  Supply/demand already specified for node'
          return
        endif

        FLOW(I) = J
        NET = NET+J
        IF (J.GE.0) then
C         ALLOW FOR SOURCES/SINKS WITH ZERO BALANCE
C         I have changed this next statement from cat(j) to cat(i), since
C         that is what I think Les Proll intended
          IF (.not.(J.EQ.0.AND.CAT(I).GT.0)) then
            MSORC = MSORC+1
C           SAVE ORIGINAL SUPPLY IN LEVEL
            LEVEL(I) = J
            NEXT(I) = NEXT(MNODP1)
            NEXT(MNODP1) = I
          endif
        endif
      end do

C     TEST FOR FEASIBILITY
C      NSTOP = 8
C      IF (NET.LT.0) CALL ERRMES(NSTOP)
      IF (NET.LT.0) then
        errorflag = 21      !
        messerror = '  Problem infeasible - demand exceeds supply'
        return
      endif

C     Now get information on lower and upper bounds for each arc
      M = 0
      L = ITOTSY
      KOST0 = 0
C     Now get information on cost for each arc, and define arc in terms
C     of to and from nodes.
      DO N = 1,INMARC
        I = INFROM(N)
        J = INTO(N)
        K = INCOST(N)

C        NSTOP = 13
C        IF (L.GE.BIG) CALL ERRMES(NSTOP)
        IF (L.GE.BIG) then
          errorflag = 22      ! Improper bounds on arc flow
          messerror = '  Improper bounds on arc flow'
          return
        endif

C       FOLLOWING STATEMENT IN K&H VERSION BUT CANT SEE WHY
C       IF (L.EQ.0) L = BIG
        IF (L.LT.0) L = 0
        II = CAT(J)
        JJ = IABS(II)
C       TEST TO SEE IF CATEGORY IS FULL
        KK = ISIGN(LNODP1,II)

        IF (TXOR(KK,FROM(JJ)).LE.0) then

C         MOVE REST OF ARCS DOWN TO ACCOMODATE

C          NSTOP = 14
C          IF (MARC.EQ.LARC) CALL ERRMES(NSTOP)
          IF (MARC.EQ.LARC) then
            errorflag = 23      !
            messerror = '  Too many arcs for this version'
            return
          endif

          MARC = MARC+1
          K1000 = MARC-JJ
          M1000 = MARC
          DO J1000 = 1,K1000
            L1000 = M1000-1
            FROM(M1000) = FROM(L1000)
            COST(M1000) = COST(L1000)
            CAPAC(M1000) = CAPAC(L1000)
            TFLOOR(M1000) = TFLOOR(L1000)
            NAME(M1000) = NAME(L1000)
            M1000 = L1000
          end do

          DO J1100 = J,MNODE
            CAT(J1100) = CAT(J1100)+ISIGN(1,CAT(J1100))
          end do
        endif

C       INSERT NEW ARC
        FROM(JJ) = ISIGN(I,II)
        COST(JJ) = K
        KOST0 = KOST0+K*M
        CAPAC(JJ) = L-M
        TFLOOR(JJ) = M
        FLOW(I) = FLOW(I)-M
        FLOW(J) = FLOW(J)+M
        NAME(JJ) = TAG
        CAT(J) = ISIGN(JJ+1,II)
C       MARK FROM NODE WITH ARC RUNNING FROM IT FLAG
        ARCID(I) = -1
      end do

      RETURN
      END

C########################## End of Program Code ################################

C ******************************************************************************
C ******************************************************************************
C **                                                                          **
C **                                                                          **
C **          History of additions and modifications by KFC                   **
C **                                                                          **
C **                                                                          **
C ******************************************************************************
C **                                                                          **
C **                                                                          **
C **    2 June 1999 -- Original version by J. N. Perry                        **
C **                                                                          **
C **    29 OCT 2007 -- Recoded to avoid use of EQUIVALENCE of NODE and ARC    **
C **                   in the TRANSP subroutine                               **
C **                -- Removed BLOCK DATA and DATA statments and used         **
C **                   simple assignment to initialize variables              **
C **                -- Replaced the varible FLOOR with TFLOOR b/c FLOOR       **
C **                   is a reserved word                                     **
C **                -- Replaced the varible XOR with TXOR b/c XOR is a        **
C **                   reserved word                                          **
C **                -- Replaced the bubble sorts in SORT2K, SORT6K, SORT78K   **
C **                   and DTSORT  with slightly faster Shell sorts (all      **
C **                   could be replaced eventually with one SBR)             **
C **                                                                          **
C **    30 OCT 2007 -- Replaced the bubble sort in ICSORT with a slightly     **
C **                   faster Shell sort                                      **
C **                -- Replaced the bubble sort in AVFLPP with a slightly     **
C **                   faster Shell sort                                      **
C **                -- Declared TEMP and ITEMP in AVFLPP, which may havd      **
C **                   caused rounding errors in EDF in previous versions     **
C **                                                                          **
C **    28 NOV 2007 -- Changed Holerith variables to characters and NAME      **
C **                   to a character array                                   **
C **                                                                          **
C **    29 NOV 2007 -- Reduced sorting sbr (sort2k,sort6k,sort78k) to use     **
C **                   one universal one                                      **
C **                                                                          **
C **    30 NOV 2007 -- imposed IMPLICIT NONE on all sbr and explicitly        **
C **                   declared all variables                                 **
C **                                                                          **
C **    02 DEC 2007 -- added initialization of CAT(J100) = 0 in sbr INDAT     **
C **                   Moved STFROM(2000),STTO(2000),STFLOW(2000) to a        **
C **                   common block                                           **
C **                   Increased the dimension of DUMFRO,DUMTO,TNF,DINCX,     **
C **                   DUMFLO,DINCY,DUMDIS,DINCL in AVFLPP from 2000 to       **
C **                   4000 to account for the fact that they could           **
C **                   require indices up to 2*NOP (=TNN)                     **
C **                                                                          **
C **                   NOTE: program now tested successufully with 1300       **
C **                         sampling points                                  **
C **                                                                          **
C **    30 DEC 2007 -- Revised code to provide both parametric and non-       **
C **                   parametric versions from one code base using           **
C **                   RFLAG (=TRUE when non-parametric analysis needed       **
C **                                                                          **
C **                   Added the command-line parameter -r to switch to       **
C **                   to the "non-parametric" mode                           **
C **                                                                          **
C **    30 JAN 2008 -- Revised ERRMES so that errors from the trans-          **
C **                   portation algorithm are written to rbno6 as well       **
C **                   as the console window                                  **
C **                                                                          **
C **                   Increased the allowable number of arcs from 444444     **
C **                   to 2600000. A better estimate of the max arcs is       **
C **                   needed.                                                **
C **                                                                          **
C **    01 FEB 2008 -- Reduced CONTINUES in the Transportation algorithm      **
C **                   by restructuring DO and IF...GOTO statements           **
C **                                                                          **
C **                   Removed Arithmetic IF statements                       **
C **                                                                          **
C **                                                                          **
C **    18 FEB 2008 -- Recoded sections of the Transportation Algorithm for   **
C **                   g77 and gnufortran compliance                          **
C **                                                                          **
C **                   Reduced rudundant common blocks                        **
C **                   Reduced rudundant variable declarations                **
C **                                                                          **
C **                   converted GOTOs in  the Transportation Algorithm to    **
C **                   if...then structures                                   **
C **                                                                          **
C **                   Tarted up the code                                     **
C **                                                                          **
C **                                                                          **
C **    12 MAR 2008 -- Added SAVE sbrs                                        **
C **                                                                          **
C **                   Tried to enforce stricter formatting of output         **
C **                                                                          **
C **                                                                          **
C **    21 MAR 2008 -- Changed rbno6 from I/O unit 6 to unit 12 to avoid      **
C **                   conflicts with console output                          **
C **                                                                          **
C **                   Changed rbno6 from STATUS=NEW to STATUS=REPLACE so     **
C **                   it can be used for logging errors                      **
C **                                                                          **
C **                   Consolidated and simplified (??) error-handling by     **
C **                   moving error handling to the end of the main prog.     **
C **                   Also added the errs common block to share error        **
C **                   messages among modules                                 **
C **                                                                          **
C **                                                                          **
C **    23 MAR 2008 -- Changed rbni5 from I/O unit 5 to unit 15 to avoid      **
C **                   conflicts with console output                          **
C **                                                                          **
C **                   Added the -i command-line option to prompt the user    **
C **                   for the input file and iseed and kp5sim                **
C **                                                                          **
C **                                                                          **
C **    26 MAR 2008 -- Implemented rudimentary support for project folders    **
C **                   Note that no folders can be created and there is       **
C **                   no checking validity of project folders                **
C **                                                                          **
C **                   Added the filename read to rbno6                       **
C **                                                                          **
C **                   Note that implementing folder handling means that      **
C **                   code is no longer portable between MS Fortran and      **
C **                   g77 (and linux) b/c the system and CHAR functions      **
C **                   are non-standard Fortran                               **
C **                                                                          **
C **                                                                          **
C **    23 JLY 2008 -- Re-defined INTEGER as INTEGER*4                        **
C **                   Removed SAVE from main routine DTSORT, DISCRO          **
C **                   Switched from Shell to Bubble sort in AVFLPP           **
C **                   these changes removed discrepancies with distance to   **
C **                   crowding compared with RBRELV13                        **
C **                                                                          **
C **                   Verified rounding error discrepancies caused by        **
C **                    previously inplicitly defined TEMP in AVFLPP          **
C **                    TEMP now declared DOUBLE PRECISION                    **
C **                                                                          **
C **                                                                          **
C ******************************************************************************
C ******************************************************************************

C ==============================================================================
C ===                                                                        ===
C ===  To Do:                                                                ===
C ===                                                                        ===
C ===                                                                        ===
C ===                                                                        ===
C ===                                                                        ===
C ==============================================================================



c                    GNU GENERAL PUBLIC LICENSE
c                       Version 3, 29 June 2007
c
c Copyright (C) 2007 Free Software Foundation, Inc. <http://fsf.org/>
c Everyone is permitted to copy and distribute verbatim copies
c of this license document, but changing it is not allowed.
c
c                            Preamble
c
c  The GNU General Public License is a free, copyleft license for
c software and other kinds of works.
c
c  The licenses for most software and other practical works are designed
c to take away your freedom to share and change the works.  By contrast,
c the GNU General Public License is intended to guarantee your freedom to
c share and change all versions of a program--to make sure it remains free
c software for all its users.  We, the Free Software Foundation, use the
c GNU General Public License for most of our software; it applies also to
c any other work released this way by its authors.  You can apply it to
c your programs, too.
c 
c   When we speak of free software, we are referring to freedom, not
c price.  Our General Public Licenses are designed to make sure that you
c have the freedom to distribute copies of free software (and charge for
c them if you wish), that you receive source code or can get it if you
c want it, that you can change the software or use pieces of it in new
c free programs, and that you know you can do these things.
c 
c   To protect your rights, we need to prevent others from denying you
c these rights or asking you to surrender the rights.  Therefore, you have
c certain responsibilities if you distribute copies of the software, or if
c you modify it: responsibilities to respect the freedom of others.
c 
c   For example, if you distribute copies of such a program, whether
c gratis or for a fee, you must pass on to the recipients the same
c freedoms that you received.  You must make sure that they, too, receive
c or can get the source code.  And you must show them these terms so they
c know their rights.
c 
c   Developers that use the GNU GPL protect your rights with two steps:
c (1) assert copyright on the software, and (2) offer you this License
c giving you legal permission to copy, distribute and/or modify it.
c 
c   For the developers' and authors' protection, the GPL clearly explains
c that there is no warranty for this free software.  For both users' and
c authors' sake, the GPL requires that modified versions be marked as
c changed, so that their problems will not be attributed erroneously to
c authors of previous versions.
c 
c   Some devices are designed to deny users access to install or run
c modified versions of the software inside them, although the manufacturer
c can do so.  This is fundamentally incompatible with the aim of
c protecting users' freedom to change the software.  The systematic
c pattern of such abuse occurs in the area of products for individuals to
c use, which is precisely where it is most unacceptable.  Therefore, we
c have designed this version of the GPL to prohibit the practice for those
c products.  If such problems arise substantially in other domains, we
c stand ready to extend this provision to those domains in future versions
c of the GPL, as needed to protect the freedom of users.
c 
c   Finally, every program is threatened constantly by software patents.
c States should not allow patents to restrict development and use of
c software on general-purpose computers, but in those that do, we wish to
c avoid the special danger that patents applied to a free program could
c make it effectively proprietary.  To prevent this, the GPL assures that
c patents cannot be used to render the program non-free.
c 
c   The precise terms and conditions for copying, distribution and
c modification follow.
c 
c                        TERMS AND CONDITIONS
c 
c   0. Definitions.
c 
c   "This License" refers to version 3 of the GNU General Public License.
c 
c   "Copyright" also means copyright-like laws that apply to other kinds of
c works, such as semiconductor masks.
c 
c   "The Program" refers to any copyrightable work licensed under this
c License.  Each licensee is addressed as "you".  "Licensees" and
c "recipients" may be individuals or organizations.
c  
c   To "modify" a work means to copy from or adapt all or part of the work
c in a fashion requiring copyright permission, other than the making of an
c exact copy.  The resulting work is called a "modified version" of the
c earlier work or a work "based on" the earlier work.
c 
c   A "covered work" means either the unmodified Program or a work based
c on the Program.
c 
c   To "propagate" a work means to do anything with it that, without
c permission, would make you directly or secondarily liable for
c infringement under applicable copyright law, except executing it on a
c computer or modifying a private copy.  Propagation includes copying,
c distribution (with or without modification), making available to the
c public, and in some countries other activities as well.
c 
c   To "convey" a work means any kind of propagation that enables other
c parties to make or receive copies.  Mere interaction with a user through
c a computer network, with no transfer of a copy, is not conveying.
c 
c   An interactive user interface displays "Appropriate Legal Notices"
c to the extent that it includes a convenient and prominently visible
c feature that (1) displays an appropriate copyright notice, and (2)
c tells the user that there is no warranty for the work (except to the
c extent that warranties are provided), that licensees may convey the
c work under this License, and how to view a copy of this License.  If
c the interface presents a list of user commands or options, such as a
c menu, a prominent item in the list meets this criterion.
c 
c   1. Source Code.
c 
c   The "source code" for a work means the preferred form of the work
c for making modifications to it.  "Object code" means any non-source
c form of a work.
c 
c   A "Standard Interface" means an interface that either is an official
c standard defined by a recognized standards body, or, in the case of
c interfaces specified for a particular programming language, one that
c is widely used among developers working in that language.
c 
c   The "System Libraries" of an executable work include anything, other
c than the work as a whole, that (a) is included in the normal form of
c packaging a Major Component, but which is not part of that Major
c Component, and (b) serves only to enable use of the work with that
c Major Component, or to implement a Standard Interface for which an
c implementation is available to the public in source code form.  A
c "Major Component", in this context, means a major essential component
c (kernel, window system, and so on) of the specific operating system
c (if any) on which the executable work runs, or a compiler used to
c produce the work, or an object code interpreter used to run it.
c 
c   The "Corresponding Source" for a work in object code form means all
c the source code needed to generate, install, and (for an executable
c work) run the object code and to modify the work, including scripts to
c control those activities.  However, it does not include the work's
c System Libraries, or general-purpose tools or generally available free
c programs which are used unmodified in performing those activities but
c which are not part of the work.  For example, Corresponding Source
c includes interface definition files associated with source files for
c the work, and the source code for shared libraries and dynamically
c linked subprograms that the work is specifically designed to require,
c such as by intimate data communication or control flow between those
c subprograms and other parts of the work.
c 
c   The Corresponding Source need not include anything that users
c can regenerate automatically from other parts of the Corresponding
c Source.
c 
c   The Corresponding Source for a work in source code form is that
c same work.
c 
c   2. Basic Permissions.
c 
c   All rights granted under this License are granted for the term of
c copyright on the Program, and are irrevocable provided the stated
c conditions are met.  This License explicitly affirms your unlimited
c permission to run the unmodified Program.  The output from running a
c covered work is covered by this License only if the output, given its
c content, constitutes a covered work.  This License acknowledges your
c rights of fair use or other equivalent, as provided by copyright law.
c 
c   You may make, run and propagate covered works that you do not
c convey, without conditions so long as your license otherwise remains
c in force.  You may convey covered works to others for the sole purpose
c of having them make modifications exclusively for you, or provide you
c with facilities for running those works, provided that you comply with
c the terms of this License in conveying all material for which you do
c not control copyright.  Those thus making or running the covered works
c for you must do so exclusively on your behalf, under your direction
c and control, on terms that prohibit them from making any copies of
c your copyrighted material outside their relationship with you.
c 
c   Conveying under any other circumstances is permitted solely under
c the conditions stated below.  Sublicensing is not allowed; section 10
c makes it unnecessary.
c 
c   3. Protecting Users' Legal Rights From Anti-Circumvention Law.
c 
c   No covered work shall be deemed part of an effective technological
c measure under any applicable law fulfilling obligations under article
c 11 of the WIPO copyright treaty adopted on 20 December 1996, or
c similar laws prohibiting or restricting circumvention of such
c measures.
c 
c   When you convey a covered work, you waive any legal power to forbid
c circumvention of technological measures to the extent such circumvention
c is effected by exercising rights under this License with respect to
c the covered work, and you disclaim any intention to limit operation or
c modification of the work as a means of enforcing, against the work's
c users, your or third parties' legal rights to forbid circumvention of
c technological measures.
c 
c   4. Conveying Verbatim Copies.
c 
c   You may convey verbatim copies of the Program's source code as you
c receive it, in any medium, provided that you conspicuously and
c appropriately publish on each copy an appropriate copyright notice;
c keep intact all notices stating that this License and any
c non-permissive terms added in accord with section 7 apply to the code;
c keep intact all notices of the absence of any warranty; and give all
c recipients a copy of this License along with the Program.
c 
c   You may charge any price or no price for each copy that you convey,
c and you may offer support or warranty protection for a fee.
c 
c   5. Conveying Modified Source Versions.
c 
c   You may convey a work based on the Program, or the modifications to
c produce it from the Program, in the form of source code under the
c terms of section 4, provided that you also meet all of these conditions:
c 
c     a) The work must carry prominent notices stating that you modified
c     it, and giving a relevant date.
c 
c     b) The work must carry prominent notices stating that it is
c     released under this License and any conditions added under section
c     7.  This requirement modifies the requirement in section 4 to
c     "keep intact all notices".
c 
c     c) You must license the entire work, as a whole, under this
c     License to anyone who comes into possession of a copy.  This
c     License will therefore apply, along with any applicable section 7
c     additional terms, to the whole of the work, and all its parts,
c     regardless of how they are packaged.  This License gives no
c     permission to license the work in any other way, but it does not
c     invalidate such permission if you have separately received it.
c 
c     d) If the work has interactive user interfaces, each must display
c     Appropriate Legal Notices; however, if the Program has interactive
c     interfaces that do not display Appropriate Legal Notices, your
c     work need not make them do so.
c 
c   A compilation of a covered work with other separate and independent
c works, which are not by their nature extensions of the covered work,
c and which are not combined with it such as to form a larger program,
c in or on a volume of a storage or distribution medium, is called an
c "aggregate" if the compilation and its resulting copyright are not
c used to limit the access or legal rights of the compilation's users
c beyond what the individual works permit.  Inclusion of a covered work
c in an aggregate does not cause this License to apply to the other
c parts of the aggregate.
c 
c   6. Conveying Non-Source Forms.
c 
c   You may convey a covered work in object code form under the terms
c of sections 4 and 5, provided that you also convey the
c machine-readable Corresponding Source under the terms of this License,
c in one of these ways:
c 
c     a) Convey the object code in, or embodied in, a physical product
c     (including a physical distribution medium), accompanied by the
c     Corresponding Source fixed on a durable physical medium
c     customarily used for software interchange.
c 
c     b) Convey the object code in, or embodied in, a physical product
c     (including a physical distribution medium), accompanied by a
c     written offer, valid for at least three years and valid for as
c     long as you offer spare parts or customer support for that product
c     model, to give anyone who possesses the object code either (1) a
c     copy of the Corresponding Source for all the software in the
c     product that is covered by this License, on a durable physical
c     medium customarily used for software interchange, for a price no
c     more than your reasonable cost of physically performing this
c     conveying of source, or (2) access to copy the
c     Corresponding Source from a network server at no charge.
c 
c     c) Convey individual copies of the object code with a copy of the
c     written offer to provide the Corresponding Source.  This
c     alternative is allowed only occasionally and noncommercially, and
c     only if you received the object code with such an offer, in accord
c     with subsection 6b.
c 
c     d) Convey the object code by offering access from a designated
c     place (gratis or for a charge), and offer equivalent access to the
c     Corresponding Source in the same way through the same place at no
c     further charge.  You need not require recipients to copy the
c     Corresponding Source along with the object code.  If the place to
c     copy the object code is a network server, the Corresponding Source
c     may be on a different server (operated by you or a third party)
c     that supports equivalent copying facilities, provided you maintain
c     clear directions next to the object code saying where to find the
c     Corresponding Source.  Regardless of what server hosts the
c     Corresponding Source, you remain obligated to ensure that it is
c     available for as long as needed to satisfy these requirements.
c 
c     e) Convey the object code using peer-to-peer transmission, provided
c     you inform other peers where the object code and Corresponding
c     Source of the work are being offered to the general public at no
c     charge under subsection 6d.
c 
c   A separable portion of the object code, whose source code is excluded
c from the Corresponding Source as a System Library, need not be
c included in conveying the object code work.
c 
c   A "User Product" is either (1) a "consumer product", which means any
c tangible personal property which is normally used for personal, family,
c or household purposes, or (2) anything designed or sold for incorporation
c into a dwelling.  In determining whether a product is a consumer product,
c doubtful cases shall be resolved in favor of coverage.  For a particular
c product received by a particular user, "normally used" refers to a
c typical or common use of that class of product, regardless of the status
c of the particular user or of the way in which the particular user
c actually uses, or expects or is expected to use, the product.  A product
c is a consumer product regardless of whether the product has substantial
c commercial, industrial or non-consumer uses, unless such uses represent
c the only significant mode of use of the product.
c 
c   "Installation Information" for a User Product means any methods,
c procedures, authorization keys, or other information required to install
c and execute modified versions of a covered work in that User Product from
c a modified version of its Corresponding Source.  The information must
c suffice to ensure that the continued functioning of the modified object
c code is in no case prevented or interfered with solely because
c modification has been made.
c 
c   If you convey an object code work under this section in, or with, or
c specifically for use in, a User Product, and the conveying occurs as
c part of a transaction in which the right of possession and use of the
c User Product is transferred to the recipient in perpetuity or for a
c fixed term (regardless of how the transaction is characterized), the
c Corresponding Source conveyed under this section must be accompanied
c by the Installation Information.  But this requirement does not apply
c if neither you nor any third party retains the ability to install
c modified object code on the User Product (for example, the work has
c been installed in ROM).
c 
c   The requirement to provide Installation Information does not include a
c requirement to continue to provide support service, warranty, or updates
c for a work that has been modified or installed by the recipient, or for
c the User Product in which it has been modified or installed.  Access to a
c network may be denied when the modification itself materially and
c adversely affects the operation of the network or violates the rules and
c protocols for communication across the network.
c 
c   Corresponding Source conveyed, and Installation Information provided,
c in accord with this section must be in a format that is publicly
c documented (and with an implementation available to the public in
c source code form), and must require no special password or key for
c unpacking, reading or copying.
c 
c   7. Additional Terms.
c 
c   "Additional permissions" are terms that supplement the terms of this
c License by making exceptions from one or more of its conditions.
c Additional permissions that are applicable to the entire Program shall
c be treated as though they were included in this License, to the extent
c that they are valid under applicable law.  If additional permissions
c apply only to part of the Program, that part may be used separately
c under those permissions, but the entire Program remains governed by
c this License without regard to the additional permissions.
c 
c   When you convey a copy of a covered work, you may at your option
c remove any additional permissions from that copy, or from any part of
c it.  (Additional permissions may be written to require their own
c removal in certain cases when you modify the work.)  You may place
c additional permissions on material, added by you to a covered work,
c for which you have or can give appropriate copyright permission.
c 
c   Notwithstanding any other provision of this License, for material you
c add to a covered work, you may (if authorized by the copyright holders of
c that material) supplement the terms of this License with terms:
c 
c     a) Disclaiming warranty or limiting liability differently from the
c     terms of sections 15 and 16 of this License; or
c 
c     b) Requiring preservation of specified reasonable legal notices or
c     author attributions in that material or in the Appropriate Legal
c     Notices displayed by works containing it; or
c 
c     c) Prohibiting misrepresentation of the origin of that material, or
c     requiring that modified versions of such material be marked in
c     reasonable ways as different from the original version; or
c 
c     d) Limiting the use for publicity purposes of names of licensors or
c     authors of the material; or
c 
c     e) Declining to grant rights under trademark law for use of some
c     trade names, trademarks, or service marks; or
c 
c     f) Requiring indemnification of licensors and authors of that
c     material by anyone who conveys the material (or modified versions of
c     it) with contractual assumptions of liability to the recipient, for
c     any liability that these contractual assumptions directly impose on
c     those licensors and authors.
c 
c   All other non-permissive additional terms are considered "further
c restrictions" within the meaning of section 10.  If the Program as you
c received it, or any part of it, contains a notice stating that it is
c governed by this License along with a term that is a further
c restriction, you may remove that term.  If a license document contains
c a further restriction but permits relicensing or conveying under this
c License, you may add to a covered work material governed by the terms
c of that license document, provided that the further restriction does
c not survive such relicensing or conveying.
c 
c  If you add terms to a covered work in accord with this section, you
c must place, in the relevant source files, a statement of the
c additional terms that apply to those files, or a notice indicating
c where to find the applicable terms.
c 
c   Additional terms, permissive or non-permissive, may be stated in the
c form of a separately written license, or stated as exceptions;
c the above requirements apply either way.
c 
c   8. Termination.
c 
c   You may not propagate or modify a covered work except as expressly
c provided under this License.  Any attempt otherwise to propagate or
c modify it is void, and will automatically terminate your rights under
c this License (including any patent licenses granted under the third
c paragraph of section 11).
c 
c   However, if you cease all violation of this License, then your
c license from a particular copyright holder is reinstated (a)
c provisionally, unless and until the copyright holder explicitly and
c finally terminates your license, and (b) permanently, if the copyright
c holder fails to notify you of the violation by some reasonable means
c prior to 60 days after the cessation.
c 
c   Moreover, your license from a particular copyright holder is
c reinstated permanently if the copyright holder notifies you of the
c violation by some reasonable means, this is the first time you have
c received notice of violation of this License (for any work) from that
c copyright holder, and you cure the violation prior to 30 days after
c your receipt of the notice.
c 
c   Termination of your rights under this section does not terminate the
c licenses of parties who have received copies or rights from you under
c this License.  If your rights have been terminated and not permanently
c reinstated, you do not qualify to receive new licenses for the same
c material under section 10.
c 
c   9. Acceptance Not Required for Having Copies.
c 
c   You are not required to accept this License in order to receive or
c run a copy of the Program.  Ancillary propagation of a covered work
c occurring solely as a consequence of using peer-to-peer transmission
c to receive a copy likewise does not require acceptance.  However,
c nothing other than this License grants you permission to propagate or
c modify any covered work.  These actions infringe copyright if you do
c not accept this License.  Therefore, by modifying or propagating a
c covered work, you indicate your acceptance of this License to do so.
c 
c   10. Automatic Licensing of Downstream Recipients.
c 
c   Each time you convey a covered work, the recipient automatically
c receives a license from the original licensors, to run, modify and
c propagate that work, subject to this License.  You are not responsible
c for enforcing compliance by third parties with this License.
c 
c   An "entity transaction" is a transaction transferring control of an
c organization, or substantially all assets of one, or subdividing an
c organization, or merging organizations.  If propagation of a covered
c work results from an entity transaction, each party to that
c transaction who receives a copy of the work also receives whatever
c licenses to the work the party's predecessor in interest had or could
c give under the previous paragraph, plus a right to possession of the
c Corresponding Source of the work from the predecessor in interest, if
c the predecessor has it or can get it with reasonable efforts.
c 
c   You may not impose any further restrictions on the exercise of the
c rights granted or affirmed under this License.  For example, you may
c not impose a license fee, royalty, or other charge for exercise of
c rights granted under this License, and you may not initiate litigation
c (including a cross-claim or counterclaim in a lawsuit) alleging that
c any patent claim is infringed by making, using, selling, offering for
c sale, or importing the Program or any portion of it.
c 
c   11. Patents.
c 
c   A "contributor" is a copyright holder who authorizes use under this
c License of the Program or a work on which the Program is based.  The
c work thus licensed is called the contributor's "contributor version".
c 
c   A contributor's "essential patent claims" are all patent claims
c owned or controlled by the contributor, whether already acquired or
c hereafter acquired, that would be infringed by some manner, permitted
c by this License, of making, using, or selling its contributor version,
c but do not include claims that would be infringed only as a
c consequence of further modification of the contributor version.  For
c purposes of this definition, "control" includes the right to grant
c patent sublicenses in a manner consistent with the requirements of
c this License.
c 
c   Each contributor grants you a non-exclusive, worldwide, royalty-free
c patent license under the contributor's essential patent claims, to
c make, use, sell, offer for sale, import and otherwise run, modify and
c propagate the contents of its contributor version.
c 
c   In the following three paragraphs, a "patent license" is any express
c agreement or commitment, however denominated, not to enforce a patent
c (such as an express permission to practice a patent or covenant not to
c sue for patent infringement).  To "grant" such a patent license to a
c party means to make such an agreement or commitment not to enforce a
c patent against the party.
c 
c   If you convey a covered work, knowingly relying on a patent license,
c and the Corresponding Source of the work is not available for anyone
c to copy, free of charge and under the terms of this License, through a
c publicly available network server or other readily accessible means,
c then you must either (1) cause the Corresponding Source to be so
c available, or (2) arrange to deprive yourself of the benefit of the
c patent license for this particular work, or (3) arrange, in a manner
c consistent with the requirements of this License, to extend the patent
c license to downstream recipients.  "Knowingly relying" means you have
c actual knowledge that, but for the patent license, your conveying the
c covered work in a country, or your recipient's use of the covered work
c in a country, would infringe one or more identifiable patents in that
c country that you have reason to believe are valid.
c 
c   If, pursuant to or in connection with a single transaction or
c arrangement, you convey, or propagate by procuring conveyance of, a
c covered work, and grant a patent license to some of the parties
c receiving the covered work authorizing them to use, propagate, modify
c or convey a specific copy of the covered work, then the patent license
c you grant is automatically extended to all recipients of the covered
c work and works based on it.
c 
c   A patent license is "discriminatory" if it does not include within
c the scope of its coverage, prohibits the exercise of, or is
c conditioned on the non-exercise of one or more of the rights that are
c specifically granted under this License.  You may not convey a covered
c work if you are a party to an arrangement with a third party that is
c in the business of distributing software, under which you make payment
c to the third party based on the extent of your activity of conveying
c the work, and under which the third party grants, to any of the
c parties who would receive the covered work from you, a discriminatory
c patent license (a) in connection with copies of the covered work
c conveyed by you (or copies made from those copies), or (b) primarily
c for and in connection with specific products or compilations that
c contain the covered work, unless you entered into that arrangement,
c or that patent license was granted, prior to 28 March 2007.
c 
c   Nothing in this License shall be construed as excluding or limiting
c any implied license or other defenses to infringement that may
c otherwise be available to you under applicable patent law.
c 
c   12. No Surrender of Others' Freedom.
c 
c   If conditions are imposed on you (whether by court order, agreement or
c otherwise) that contradict the conditions of this License, they do not
c excuse you from the conditions of this License.  If you cannot convey a
c covered work so as to satisfy simultaneously your obligations under this
c License and any other pertinent obligations, then as a consequence you may
c not convey it at all.  For example, if you agree to terms that obligate you
c to collect a royalty for further conveying from those to whom you convey
c the Program, the only way you could satisfy both those terms and this
c License would be to refrain entirely from conveying the Program.
c 
c   13. Use with the GNU Affero General Public License.
c 
c   Notwithstanding any other provision of this License, you have
c permission to link or combine any covered work with a work licensed
c under version 3 of the GNU Affero General Public License into a single
c combined work, and to convey the resulting work.  The terms of this
c License will continue to apply to the part which is the covered work,
c but the special requirements of the GNU Affero General Public License,
c section 13, concerning interaction through a network will apply to the
c combination as such.
c 
c   14. Revised Versions of this License.
c 
c   The Free Software Foundation may publish revised and/or new versions of
c the GNU General Public License from time to time.  Such new versions will
c be similar in spirit to the present version, but may differ in detail to
c address new problems or concerns.
c 
c   Each version is given a distinguishing version number.  If the
c Program specifies that a certain numbered version of the GNU General
c Public License "or any later version" applies to it, you have the
c option of following the terms and conditions either of that numbered
c version or of any later version published by the Free Software
c Foundation.  If the Program does not specify a version number of the
c GNU General Public License, you may choose any version ever published
c by the Free Software Foundation.
c 
c   If the Program specifies that a proxy can decide which future
c versions of the GNU General Public License can be used, that proxy's
c public statement of acceptance of a version permanently authorizes you
c to choose that version for the Program.
c 
c   Later license versions may give you additional or different
c permissions.  However, no additional obligations are imposed on any
c author or copyright holder as a result of your choosing to follow a
c later version.
c 
c   15. Disclaimer of Warranty.
c 
c   THERE IS NO WARRANTY FOR THE PROGRAM, TO THE EXTENT PERMITTED BY
c APPLICABLE LAW.  EXCEPT WHEN OTHERWISE STATED IN WRITING THE COPYRIGHT
c HOLDERS AND/OR OTHER PARTIES PROVIDE THE PROGRAM "AS IS" WITHOUT WARRANTY
c OF ANY KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO,
c THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
c PURPOSE.  THE ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE PROGRAM
c IS WITH YOU.  SHOULD THE PROGRAM PROVE DEFECTIVE, YOU ASSUME THE COST OF
c ALL NECESSARY SERVICING, REPAIR OR CORRECTION.
c 
c   16. Limitation of Liability.
c 
c   IN NO EVENT UNLESS REQUIRED BY APPLICABLE LAW OR AGREED TO IN WRITING
c WILL ANY COPYRIGHT HOLDER, OR ANY OTHER PARTY WHO MODIFIES AND/OR CONVEYS
c THE PROGRAM AS PERMITTED ABOVE, BE LIABLE TO YOU FOR DAMAGES, INCLUDING ANY
c GENERAL, SPECIAL, INCIDENTAL OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE
c USE OR INABILITY TO USE THE PROGRAM (INCLUDING BUT NOT LIMITED TO LOSS OF
c DATA OR DATA BEING RENDERED INACCURATE OR LOSSES SUSTAINED BY YOU OR THIRD
c PARTIES OR A FAILURE OF THE PROGRAM TO OPERATE WITH ANY OTHER PROGRAMS),
c EVEN IF SUCH HOLDER OR OTHER PARTY HAS BEEN ADVISED OF THE POSSIBILITY OF
c SUCH DAMAGES.
c 
c   17. Interpretation of Sections 15 and 16.
c 
c   If the disclaimer of warranty and limitation of liability provided
c above cannot be given local legal effect according to their terms,
c reviewing courts shall apply local law that most closely approximates
c an absolute waiver of all civil liability in connection with the
c Program, unless a warranty or assumption of liability accompanies a
c copy of the Program in return for a fee.
c 
c                      END OF TERMS AND CONDITIONS