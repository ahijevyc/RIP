subroutine find_maxima(v, ix, jx, map_dx, map_dy, &
 maxes, sq_maxes, big_maxes, sm_maxes, scale, &
 debug_print)

!  USE GRIDINFO
!  USE STORAGE_MODULE
  IMPLICIT NONE

  LOGICAL :: debug_print 
  integer :: ix, jx
  integer :: i, j, i0, j0, dist, c, went_up, not_round, too_big, too_small
  integer :: is_negative, positive_Laplacian, number_of_MCVs
  INTEGER, PARAMETER        :: max_bin = 10
  INTEGER, PARAMETER        :: max_radius = 300, min_radius = 50
  INTEGER, PARAMETER        :: mcv_info_file_unit   = 50 
  REAL, PARAMETER           :: max_std_dev_percent = 30.
  REAL, PARAMETER           :: bin_average_value_thresh = .1
  REAL                      :: AVERAGE_DISTANCE_IN_BIN

  real, dimension(max_bin)    :: sumv, sumv2, sigv, bin_mean, total_dist_in_bin
  real, dimension(max_bin)    :: avg_dist_in_bin
  integer, dimension(max_bin) :: n

  real :: outside4, inner4, Laplacian, bin_spacing, variance
  real :: map_dx, map_dy
  real :: grid_spacing, tmp
  real, dimension(ix,jx) :: v, maxes, sq_maxes, big_maxes, sm_maxes, scale

  LOGICAL           :: opened
  LOGICAL           :: found_one_mcv = .FALSE.
  INTEGER           :: iostat, k
  CHARACTER(LEN=8)  :: now_date
  CHARACTER(LEN=10) :: now_time
  CHARACTER(LEN=5)  :: now_zone
  CHARACTER(LEN=500):: format
  CHARACTER(LEN=256):: mcv_info_file = 'mcv_info_file.txt'
  CHARACTER(LEN=19) :: hdate
  CHARACTER(LEN=19) :: hstart
  INTEGER           :: idts



  ! Go to the end of the mcv information file so we can append to it. 

  INQUIRE(mcv_info_file_unit, OPENED=opened)

  if(debug_print)write(*,*)'opened(T/F):',opened
  if(debug_print)write(*,*)'mcv_info_file:'//mcv_info_file
  IF (.NOT.opened) OPEN(mcv_info_file_unit,FILE=mcv_info_file)

  iostat = 0
  DO WHILE (iostat .EQ. 0)
     READ(mcv_info_file_unit,*,IOSTAT=iostat)
  END DO

  CALL DATE_AND_TIME(DATE = now_date, TIME = now_time, ZONE = now_zone)
  format = "(A2,':',A2,1X,A5,3X,A4,'/',A2,'/',A2)"
  WRITE(mcv_info_file_unit,format,ADVANCE='NO') now_time(1:2), now_time(3:4), &
        now_zone,  now_date(1:4), now_date(5:6), now_date(7:8)
  WRITE(mcv_info_file_unit,*)" entering subroutine find_maxima"
  WRITE(mcv_info_file_unit,*)
  




  bin_spacing = AMAX1(map_dx,map_dy)    ! AMAX1 input real output real
  grid_spacing = bin_spacing

  if(debug_print)WRITE(*,*) 'map_dx = ', map_dx, ' map_dy = ', map_dy

  ! In the RUC2 grib data, the map_dx, map_dy maximum is 40.6, not the
  ! expected 40. -- Dave

  !max_bin = nint(max_radius/bin_spacing)  ! AINT truncates fraction part

  WRITE(*,*) 'max_bin = ', max_bin
  WRITE(*,*) 'bin spacing = ', bin_spacing
  WRITE(mcv_info_file_unit,*) 'number of bins = ', max_bin
  WRITE(mcv_info_file_unit,*) 'bin spacing = ', bin_spacing
  format = '(1X,A16,I4," (actually ",F6.1,")")'
  WRITE(mcv_info_file_unit,format) 'minimum radius = ', min_radius, &
     AVERAGE_DISTANCE_IN_BIN(IFIX(min_radius/bin_spacing)+1,bin_spacing,grid_spacing)
  WRITE(mcv_info_file_unit,format) 'maximum radius = ', max_radius, &
     AVERAGE_DISTANCE_IN_BIN(IFIX(max_radius/bin_spacing),bin_spacing,grid_spacing)
  WRITE(mcv_info_file_unit,*) 'maximum standard deviation of values &
                              &in a bin as percentage of central value = ', &
                               max_std_dev_percent
  WRITE(mcv_info_file_unit,*) 'average value in points in bin must drop to this &
                              &fraction of central value = ', &
                               bin_average_value_thresh



  maxes     = 0      ! Initialize maxima field to zeroes.
  big_maxes = 0      ! Initialize maxima field to zeroes.
  sq_maxes  = 0      ! Initialize maxima field to zeroes.
  scale     = 0      ! Initialize maxima field to zeroes.
  sm_maxes  = 0      ! Initialize maxima field to zeroes.
  went_up   = 0
  not_round = 0
  too_big   = 0
  too_small = 0
  is_negative        = 0
  positive_Laplacian = 0
  number_of_MCVs     = 0


  ILOOP : do i = 1+1, ix-1
     JLOOP : do j = 1+1, jx-1

        if (.false.)write(*,'("v(",i3,",",i3,")=",E12.3)')i,j,v(i,j)
        IF (mod(i,25)==0 .AND. mod(j,25)==0) then

           ! Create the SCALE field
           ! Every 25 grid points, make scale(i0,j0)=1 between the min/max radiiof MCVs.

           do i0 = 1,ix
              do j0 = 1,jx
                 dist = SQRT( (map_dx*(i-i0))**2 + (map_dy*(j-j0))**2 )
                 IF (dist <= max_radius .AND. &
                    dist >  min_radius ) scale(i0,j0)=1
              ENDDO
           ENDDO
        ENDIF
                     
        IF (v(i,j) < 0) THEN
           is_negative = is_negative + 1 
           CYCLE JLOOP   !  If vorticity is negative, skip the point.
        ENDIF

        ! A quick check to see if the Laplacian at the point in 
        ! question is negative (denoting a local maximum).
        ! If it is positive (denoting a local minimum),
        ! then skip to the next point.

        ! Compute the Laplacian
        

        outside4 = v(i,j+1) + v(i,j-1) + v(i-1,j) + v(i+1,j)
        inner4   = 4 * v(i,j)
        Laplacian = outside4 - inner4
        

        IF (debug_print.and. i.eq.15.AND.j.eq.15) then
           WRITE(*,*)
           WRITE(*,'(12x,E12.3/)')          v(i,j+1)
           WRITE(*,'(3E12.3/)')    v(i-1,j), v(i,j), v(i+1,j)
           WRITE(*,'(12x,E12.3/)')          v(i,j-1)
           print*
           print*, ' outside4 = ', outside4
           print*, ' inner4   = ', inner4
           print*
           WRITE(*,'("Laplacian = ",E12.3)') Laplacian
        endif
        
        IF (Laplacian > 0) THEN
           positive_Laplacian = positive_Laplacian + 1
           CYCLE JLOOP ! If point is local min, skip it.
        ENDIF

        !  Loop thru the horizontal grid domain.
        !  If the point belongs in one of the range bins, which one?
        !  And keep running totals of v, v squared, the number of v
        !  in each range bin, and the sum of the distances from the
        !  center point for each bin.

        n     = 0
        sumv  = 0.
        sumv2 = 0.
        total_dist_in_bin = 0.
        sigv  = 999
        bin_mean   = -999

        INNER_ILOOP : do i0 = 1, ix

           ! If one component of the distance (zonal or meridional) exceeds the max_radius,
           ! then the total distance must exceed max_radius.  Therefore, check the compenents
           ! and cut our losses if one exceeds max_radius.

           ! Move on if the zonal distance exceeds max_radius.
           ! Note, iabs() is the absolute value operation for integers
           if ( map_dx * iabs(i0-i) > max_radius) CYCLE INNER_ILOOP


           INNER_JLOOP : do j0 = 1, jx

              ! Move on if meridional distance exceeds max_radius.
              if ( map_dy * iabs(j0-j) > max_radius) CYCLE INNER_JLOOP

              dist = (map_dx*(i0-i))**2 + (map_dy*(j0-j))**2
              if (dist .LT. 0.)then
                  write(*,*)'i0,j0,i,j=',i0,j0,i,j
                  write(*,*)'dist squared=',dist
                  STOP
              endif
              dist = SQRT( (map_dx*(i0-i))**2 + (map_dy*(j0-j))**2 )

              ! Move on if total distance exceeds max_radius.
              IF (dist > max_radius) CYCLE INNER_JLOOP

              ! Assign the point to a bin and adjust the bin statistics.
              BIN_ASSIGN_LOOP: do k = 1, max_bin

                 IF ( dist > (k-1)*bin_spacing .AND. dist <= k*bin_spacing ) then


                    n(k)     = n(k)     + 1
                    sumv(k)  = sumv(k)  + v(i0,j0)
                    sumv2(k) = sumv2(k) + v(i0,j0)**2

                    total_dist_in_bin(k) = total_dist_in_bin(k) + dist

                    if(.FALSE.)then
                       write(*,*)"BIN_ASSIGN_LOOP: i,j,i0,j0=",i,j,i0,j0, " bin ", k
                       print '(" pt ",i3," in bin ",i2)', n(k),k
                       print '("sumv  ",E10.3)', sumv(k)
                       print '("sumv2 ",E10.3)', sumv2(k)
                    endif

                    CYCLE INNER_JLOOP

                 ENDIF


              END DO BIN_ASSIGN_LOOP


           END DO INNER_JLOOP
        END DO INNER_ILOOP



        !  Loop thru the range bins, from smallest to largest.
        !  Calculate the standard deviation of v within each range bin.
        !  
        !  Abandon the point in question if anything in the
        !  v pattern or sig(v) pattern suggests that it is
        !  not a local maxima associated with an MCV.

        BIN_STATS_LOOP: do k = 1, max_bin

           if(.FALSE.)write(*,*)"BIN_STATS_LOOP: k=", k

           if(n(k).le.0) CYCLE BIN_STATS_LOOP

           bin_mean(k) = sumv(k)/n(k)


           avg_dist_in_bin(k) = total_dist_in_bin(k)/n(k)

           !----------------------------------------------------------------------
           ! If average value of v in the current range bin is greater than the
           ! average value in the adjacent, interior range bin, skip the point.
           ! Average v should decrease monotonically with distance from center
           ! point.
           ! I could allow for 5% deviations.         

           IF (k.GT.1 .AND. bin_mean(k) .GT. bin_mean(k-1) ) THEN
              IF(.FALSE.)print*, "not centered on maximum!"
              went_up = went_up + 1    
              EXIT BIN_STATS_LOOP
           ENDIF

           !----------------------------------------------------------------------



         

           !----------------------------------------------------------------------
           ! If the standard deviation of v in the current range bin is
           ! too high, the point in question is probably not the
           ! center of a quasi-circular maximum in the v field.
           ! Therefore, skip the point. 

           ! Compare sigma to the value of the point in question, 
           ! not the average v within the current bin.  
           ! The average value could easily drop below zero.

           IF (n(k).le.1) CYCLE BIN_STATS_LOOP ! Can't divide by zero and
                                         ! need more than 1 value to
                                         ! calculate sigma
           
           sigv(k) = SQRT( (sumv2(k) - sumv(k)**2/n(k)) / (n(k)-1) )
           IF (sigv(k) > max_std_dev_percent/100. *v(i,j)) then
               sq_maxes(i,j)=1
               !print*, "v: ",v(i,j)," bin ",k," sigma: ",sigv(k)," avg: ", bin_mean(k)
               not_round = not_round + 1
               EXIT BIN_STATS_LOOP
           endif

           !-----------------------------------------------------------------------

           IF (avg_dist_in_bin(k) > max_radius) then
           !IF (k * bin_spacing > max_radius) then
              format = '("(",i3,",",i3,")")'
              if(debug_print)then
                 WRITE(*,format,ADVANCE='NO') i, j
                 WRITE(*,*) " too big"
              endif
              too_big = too_big + 1
              big_maxes(i,j)=1 ! If we're at the furthest range and 
                               ! still haven't dropped 90%, remember (i,j)
                               ! as a "big max".
              EXIT BIN_STATS_LOOP
           ENDIF







           ! The point in question qualifies as a local maxima if
           ! the average v drops to a fraction of its value within 
           ! a distance that is not too small (we don't want every
           ! little jiggle in the data to be labelled as a max)
           ! and not too large.
 
           ! If the point is qualified,
           ! assign a value of 1 to the maxes field at (i,j).

           IF (bin_mean(k) < bin_average_value_thresh*v(i,j) ) then

              IF ( avg_dist_in_bin(k) > min_radius) then
              !IF ( k * bin_spacing > min_radius) then
                 maxes(i,j) = 1
                 number_of_MCVs = number_of_MCVs +1
                 IF(.NOT. found_one_mcv) THEN
                    format = '(/"grid coord center value     size      bin #")'
                    WRITE(mcv_info_file_unit,format)
                    found_one_mcv = .TRUE.
                 ENDIF
                 format = '("(",i3,",",i3,")",E12.3)'
                 WRITE(mcv_info_file_unit,*)
                 WRITE(mcv_info_file_unit,format,ADVANCE='NO') i, j, v(i,j)
                 WRITE(*,format,ADVANCE='NO') i, j, v(i,j)
                 WRITE(*,*) " MCV"  
                 format = '(F12.1)' 
                 WRITE(mcv_info_file_unit,format,ADVANCE='NO') avg_dist_in_bin(k)
                 WRITE(mcv_info_file_unit,'(99(1X,I9))') (i0, i0=1,k)
                 WRITE(mcv_info_file_unit,'(20X,"num of pts = ")',ADVANCE='NO')
                 WRITE(mcv_info_file_unit,'(99(1X,I9))')    (n(i0), i0=1,k)
                 WRITE(mcv_info_file_unit,'(20X,"   average = ")',ADVANCE='NO')
                 WRITE(mcv_info_file_unit,'(99(1X,E9.2))')  (bin_mean(i0), i0=1,k)
                 WRITE(mcv_info_file_unit,'(20X,"std deviat = ")',ADVANCE='NO')
                 WRITE(mcv_info_file_unit,'(99(1X,E9.2))')  (sigv(i0), i0=1,k)
                 !print '("sumv    ",10E9.2)', sumv
                 !print '("sumv2   ",10E9.2)', sumv2
              else
                 !print*, "mean dropped off too quickly with distance"
                 sm_maxes(i,j)=k
                 too_small = too_small + 1
                 format = '("(",i3,",",i3,")",E12.3)'
                 WRITE(*,format,ADVANCE='NO') i, j, v(i,j)
                 WRITE(*,*) " too small ", avg_dist_in_bin(k)  
              endif

              EXIT BIN_STATS_LOOP ! If v(i,j) dropped 90%, we're done with it.
           ENDIF



        END DO BIN_STATS_LOOP


        IF(i.EQ.96 .AND. (j.EQ.51 .OR. j.EQ.54)) THEN
           WRITE(*,*)
           WRITE(*,'("(",i3,",",i3,")",E12.3)') i,j,v(i,j)
           WRITE(*,*)
           WRITE(*,'("         bin ",99(1X,I9)  )') (k, k=1,max_bin)
           WRITE(*,'("avg bin dist ",99(1X,F9.1))') (avg_dist_in_bin(k), k=1,max_bin)
           WRITE(*,'("num of pts   ",99(1X,I9)  )') (n(k), k=1,max_bin)
           WRITE(*,'("    mean     ",99(1X,E9.2))') (bin_mean(k), k=1,max_bin)
           WRITE(*,'(" std deviat  ",99(1X,E9.2))') (sigv(k), k=1,max_bin)
           WRITE(*,*)
        ENDIF


     END DO JLOOP
  END DO ILOOP



  WRITE(mcv_info_file_unit,*)
  WRITE(mcv_info_file_unit,'(a16,2x,i5,2x,F4.1,"%")') "less than zero",  is_negative, &
     100.*is_negative/ix/jx
  WRITE(mcv_info_file_unit,'(a16,2x,i5,2x,F4.1,"%")') "Laplacian > 0", positive_Laplacian,&
     100.*positive_Laplacian/ix/jx
  WRITE(mcv_info_file_unit,'(a16,2x,i5,2x,F4.1,"%")') "incr with dist",   went_up, &
     100.*went_up/ix/jx
  WRITE(mcv_info_file_unit,'(a16,2x,i5,2x,F4.1,"%")') "not round", not_round, &
     100.*not_round/ix/jx
  WRITE(mcv_info_file_unit,'(a16,2x,i5,2x,F4.1,"%")') "too big",   too_big, &
     100.*too_big/ix/jx
  WRITE(mcv_info_file_unit,'(a16,2x,i5,2x,F4.1,"%")') "too small", too_small, &
     100.*too_small/ix/jx
  WRITE(mcv_info_file_unit,'(a16,2x,i5,2x,F4.1,"%")') "MCVs", number_of_MCVs, &
     100.*number_of_MCVs/ix/jx
  WRITE(mcv_info_file_unit,*)


!  call put_storage (201300, 'MAXES', maxes, ix, jx)
!  call put_storage (201300, 'BIGMAX', big_maxes, ix, jx)
!  call put_storage (201300, 'SMLMAX', sm_maxes, ix, jx)
!  call put_storage (201300, 'SQMAX', sq_maxes, ix, jx)
!  call put_storage (201300, 'SCALE', scale, ix, jx)



  INQUIRE(mcv_info_file_unit, OPENED=opened)

  CALL DATE_AND_TIME(DATE = now_date, TIME = now_time, ZONE = now_zone)
  format = "(A2,':',A2,1X,A5,3X,A4,'/',A2,'/',A2)"
  WRITE(mcv_info_file_unit,format,ADVANCE='NO') now_time(1:2), now_time(3:4), &
              now_zone,  now_date(1:4), now_date(5:6), now_date(7:8)

  WRITE(mcv_info_file_unit,*)" exiting subroutine find_maxima"
  WRITE(mcv_info_file_unit,*)


  IF (opened) CLOSE(mcv_info_file_unit)

  write(*,*) 'exiting subroutine find_maxima'



end subroutine find_maxima





REAL FUNCTION AVERAGE_DISTANCE_IN_BIN(bin,bin_spacing,grid_spacing)


    IMPLICIT NONE

    REAL                    :: bin_spacing
    REAL                    :: grid_spacing
    REAL                    :: dist, total_dist_in_bin

     
    INTEGER                 :: bin, i, j, err, x, y
    INTEGER                 :: num_in_bin
    
    x=50
    y=50
    num_in_bin=0
    total_dist_in_bin=0.
    DO i = 1,100
       DO j = 1,100

          dist = SQRT( (grid_spacing*(i-x))**2 + (grid_spacing*(j-y))**2 )
           
          IF ( (bin-1)*bin_spacing .LE. dist .AND. &
               ( bin )*bin_spacing .GT. dist   ) THEN

             num_in_bin          = num_in_bin         + 1
             total_dist_in_bin   = total_dist_in_bin  + dist


          ENDIF
       ENDDO
    ENDDO

    AVERAGE_DISTANCE_IN_BIN = total_dist_in_bin/num_in_bin

    RETURN

END FUNCTION AVERAGE_DISTANCE_IN_BIN

