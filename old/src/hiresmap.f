      SUBROUTINE HIRESMAP(IUNIT)                                                 HIRESMAP.1
      DIMENSION XLATMAP(10000),XLONMAP(10000)                                    HIRESMAP.2
C                                                                                HIRESMAP.3
c  This routine was adapted from Graph (hence all the ugly upper case).
c  It assumes the datafile has been opened and all the line width, color,
c  etc. has already been set.
c  Jim Bresch 5/6/99.
c
C     ... THIS ROUTINE CALLS ROUTINE TO READ THE HIGH RESOLUTION                 HIRESMAP.18
C         DATA ONCE FOR EACH BLOCK OF VALUES                                     HIRESMAP.19
C     ... THE NCARG ROUTINES MAPIT AND MAPIQ PLOT OUT THE LINES                  HIRESMAP.20
C         AND FLUSH THE BUFFERS AFTER EACH LINE SEGMENT IS DRAWN                 HIRESMAP.21
C     ... THE ROUTINE RETURNS CONTROL WHEN ALL OF THE DATA HAS                   HIRESMAP.22
C         BEEN PROCESSED                                                         HIRESMAP.23
                                                                                 HIRESMAP.25
      IBLOCK=1                                                                   HIRESMAP.29
10    CONTINUE                                                                   HIRESMAP.30
      CALL READMAP(NUMPOINTS,                                                    HIRESMAP.31
     *                   XLATMIN,XLATMAX,XLONMIN,XLONMAX,                        HIRESMAP.32
     *                   XLATMAP,XLONMAP,                                        HIRESMAP.33
     *                   IUNIT,IEND)                                             HIRESMAP.34
                                                                                 HIRESMAP.35
      IF(IEND.EQ.1) GOTO 1000                                                    HIRESMAP.36
C     PRINT *,'FINISHED BLOCK ',IBLOCK                                           HIRESMAP.37
      CALL MAPIT(XLATMAP(1),XLONMAP(1),0)                                        HIRESMAP.38
      DO 200 IPOINT=2,NUMPOINTS/2                                                HIRESMAP.39
         CALL MAPIT(XLATMAP(IPOINT),XLONMAP(IPOINT),2)                           HIRESMAP.40
200   CONTINUE                                                                   HIRESMAP.41
      CALL MAPIQ                                                                 HIRESMAP.42
      IBLOCK=IBLOCK+1                                                            HIRESMAP.43
      GOTO 10                                                                    HIRESMAP.44
1000  CONTINUE                                                                   HIRESMAP.45
c     PRINT *,'MAPDRV - END OF DATA FOR HIPONE.ASCII'                            HIRESMAP.46
      RETURN                                                                     HIRESMAP.48
      END                                                                        HIRESMAP.49
                                                                                 HIRESMAP.50
                                                                                 HIRESMAP.51
      SUBROUTINE READMAP(NUMPOINTS,                                              HIRESMAP.52
     *                   XLATMIN,XLATMAX,XLONMIN,XLONMAX,                        HIRESMAP.53
     *                   XLAT,XLON,                                              HIRESMAP.54
     *                   IUNIT,IEND)                                             HIRESMAP.55
C                                                                                HIRESMAP.56
C     ... THIS PROGRAM READS THE ASCII DATA AND STORES IT IN                     HIRESMAP.57
C         ARRAYS OF LAT AND LON, AND THE NUMBER OF POINTS READ                   HIRESMAP.58
C     ... THE BOX SURROUNDING THE AREA ARE READ IN                               HIRESMAP.59
C                                                                                HIRESMAP.60
      DIMENSION XLAT(NUMPOINTS/2),XLON(NUMPOINTS/2),VAL(8)                       HIRESMAP.61
C                                                                                HIRESMAP.62
      READ(IUNIT,100,END=1000) NUMPOINTS,                                        HIRESMAP.63
     *                         XLATMIN,XLATMAX,XLONMIN,XLONMAX,                  HIRESMAP.64
     *                         XLATSTART,XLONSTART                               HIRESMAP.65
      XLAT(1)=XLATSTART                                                          HIRESMAP.66
      XLON(1)=XLONSTART                                                          HIRESMAP.67
      NUMLINESFULL=(NUMPOINTS-2)/8                                               HIRESMAP.68
C     PRINT *,'NUMBER OF FULL LINES=',NUMLINESFULL                               HIRESMAP.69
      NUMLEFT=NUMPOINTS-2-8*NUMLINESFULL                                         HIRESMAP.70
C     PRINT *,'NUMBER OF LEFT OVERS=',NUMLEFT                                    HIRESMAP.71
      ISTART=1                                                                   HIRESMAP.72
      DO 50 LINE=1,NUMLINESFULL                                                  HIRESMAP.73
         READ(IUNIT,108) VAL                                                     HIRESMAP.74
         DO 40 LOC=1,4                                                           HIRESMAP.75
            XLAT((LINE-1)*4+LOC+ISTART)=VAL(LOC*2-1)                             HIRESMAP.76
            XLON((LINE-1)*4+LOC+ISTART)=VAL(LOC*2  )                             HIRESMAP.77
40       CONTINUE                                                                HIRESMAP.78
50    CONTINUE                                                                   HIRESMAP.79
      IF     (NUMLEFT.EQ.6) THEN                                                 HIRESMAP.80
         READ(IUNIT,106) (VAL(II),II=1,NUMLEFT)                                  HIRESMAP.81
      ELSE IF(NUMLEFT.EQ.4) THEN                                                 HIRESMAP.82
         READ(IUNIT,104) (VAL(II),II=1,NUMLEFT)                                  HIRESMAP.83
      ELSE IF(NUMLEFT.EQ.2) THEN                                                 HIRESMAP.84
         READ(IUNIT,102) (VAL(II),II=1,NUMLEFT)                                  HIRESMAP.85
      ENDIF                                                                      HIRESMAP.86
      DO 60 LOC=1,NUMLEFT/2                                                      HIRESMAP.87
         XLAT(NUMLINESFULL*4+LOC+ISTART)=VAL(LOC*2-1)                            HIRESMAP.88
         XLON(NUMLINESFULL*4+LOC+ISTART)=VAL(LOC*2  )                            HIRESMAP.89
60    CONTINUE                                                                   HIRESMAP.90
                                                                                 HIRESMAP.91
C     PRINT *,'READ ALL DATA FROM THIS BLOCK'                                    HIRESMAP.92
      IEND=0                                                                     HIRESMAP.93
      RETURN                                                                     HIRESMAP.94
                                                                                 HIRESMAP.95
1000  CONTINUE                                                                   HIRESMAP.96
C     PRINT *,'END OF DATA'                                                      HIRESMAP.97
      IEND=1                                                                     HIRESMAP.98
      RETURN                                                                     HIRESMAP.99
                                                                                 HIRESMAP.100
100   FORMAT(I4,14X,6F9.3)                                                       HIRESMAP.101
102   FORMAT(2F9.3)                                                              HIRESMAP.102
104   FORMAT(4F9.3)                                                              HIRESMAP.103
106   FORMAT(6F9.3)                                                              HIRESMAP.104
108   FORMAT(8F9.3)                                                              HIRESMAP.105
                                                                                 HIRESMAP.106
      END                                                                        HIRESMAP.107
