PROGRAM gbm
	IMPLICIT NONE
	
	INTEGER, PARAMETER :: period=100000, m=2000
	REAL :: drift, vola, start_price(1:period), new_price(1:period), change_price(1:period)
	REAL :: mean, stdev, median,var, mean2, total(1:period), f_kelly, vary_vola(1:period)
	REAL :: start_money(1:period), new_money(1:period), money(1:period), vary_drift(1:period)
	REAL :: gmean(1:m), f(1:m), x, please, kelly(1:period)
	INTEGER :: i, iseed, first, t
	REAL, EXTERNAL :: ran_gauss
	DOUBLE PRECISION :: share_price(period)

    !setting parameters
	vola = 0.095
	drift = 0.003
	mean = 0.0
	stdev = 1.0
	f_kelly = drift/vola**2
    start_price = 300.0
	first=0
	iseed=99455

    !loop for stochastic vola
	DO i = 1,period
	    vary_vola(i) = 0.075 * (ran_gauss(iseed,first,mean,stdev))
	    vary_drift(i) = 0.003 * (ran_gauss(iseed,first,mean,stdev))
    END DO

	!loop for random share prices over a certain period
	DO i = 1, period
		change_price(i) = drift * start_price(i) + vola * start_price(i) * (ran_gauss(iseed,first,mean,stdev))
		new_price(i) = start_price(i) + change_price(i)
		share_price(i) = new_price(i)
		start_price(i) = new_price(i)
	END DO

	OPEN(20,file='gbm.dat')

	DO i = 1, period
		WRITE(20,*) i, share_price(i)
	END DO

	!loop for the cumulative distribution function
	CALL Bubble_sort(share_price,period)
    median = share_price(period/2)

	OPEN(19,file='cdf.dat')

	DO i = 1,period
		WRITE(19,*) share_price(i),dble(period-i)/dble(period)
	END DO

	!Kelly Criterion
	OPEN(17,file='kelly2.dat')

    DO t= 1,m

		f(t) = real(t)/real(m)
		x = f(t)
		start_money = 100000.0
		start_price =  100000.0 * x

            !loop for the change in funds
            DO i = 1,period

                change_price(i) = drift * start_price(i) + vola * start_price(i) * (ran_gauss(iseed,first,mean,stdev))
                new_price(i) = start_price(i) + change_price(i)
                share_price(i) = new_price(i)
                start_price(i) = new_price(i)
                new_money(i) = start_money(i) + x*change_price(i)
                money(i) = new_money(i)
                start_money(i) = new_money(i)
                kelly(i) = money(i)/100000.0

			END DO

			please = kelly(1)

            !loop for geometric mean
			DO i = 1,(period)
			    please = kelly(i)*please
			END DO
			gmean(t) = please**(1.0/(period)) - 1.0
	END DO

	DO t = 1,m
		WRITE(17,*) f(t), gmean(t)
	END DO
	
END PROGRAM


!random gaussian number
real function ran_gauss(iseed,first,mean,stdev)

  implicit none
  integer :: iseed,first
  real :: mean,stdev
  real :: v1,v2,rsq,fac
  real, external :: rand
  1 continue
  v1=2.0*rand(iseed,first)-1.0
  v2=2.0*rand(iseed,first)-1.0
  rsq=v1**2+v2**2
  if(rsq.ge.1.0.or.rsq.eq.0.0)goto 1
  fac=sqrt(-2.0*log(rsq)/rsq)
  ran_gauss=v2*fac * stdev + mean
end function ran_gauss

real function rand(iseed,first)

  implicit none
  integer, parameter :: MPLIER=16807
  integer, parameter :: MODLUS=2147483647
  integer, parameter :: MOBYMP=127773
  integer, parameter :: MOMDMP=2836
  integer hvlue,lvlue,testv,nextn,first,iseed
  save nextn

  if(first == 0) THEN
    nextn=iseed
    first=1
  endif

  hvlue=nextn/mobymp
  lvlue=mod(nextn,mobymp)
  testv=mplier*lvlue-momdmp*hvlue
  if(testv > 0)then
    nextn=testv
  else
    nextn=testv+modlus
  endif
  rand = real(nextn)/real(modlus)

end function rand

!bubble sort algorithm
subroutine Bubble_Sort(a,nit)
  
  double precision :: a(nit)
  double precision :: temp
  integer :: i, j
  logical :: swapped
  do j = nit-1, 1, -1
    swapped = .FALSE.
    do i = 1, j
      if (a(i) > a(i+1)) then
        temp = a(i)
        a(i) = a(i+1)
        a(i+1) = temp
        swapped = .TRUE.
      endif
    enddo
    if (.not. swapped) exit
  enddo
end subroutine Bubble_Sort

































