	!****************************************************************************
    !     Utility subroutines
    !****************************************************************************
	! -----------------------------------------------
    SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,RPL,DDSDDT, &
        DRPLDE,DRPLDT,STRAN,DSTRAN,TIME,DTIME,TEMP, &
        DTEMP,PREDEF,DPRED,MATERL,NDI,NSHR,NTENS,NSTATV,PROPS, &
        NPROPS,COORDS,DROT,PNEWDT,CELENT, &
        DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)

    include 'aba_param.inc'

    dimension stress(ntens),statev(nstatv),ddsdde(ntens,ntens,2), &
        ddsddt(ntens),drplde(ntens),stran(ntens),dstran(ntens), &
        predef(*),dpred(*),props(nprops),coords(3),drot(3,3), &
        time(2),dfgrd0(3,3),dfgrd1(3,3)
    character*80 materl
    parameter (zero=0.d0,one=1.d0,two=2.d0,half=0.5d0)
    

    DOUBLE PRECISION FREQ, FAC1, FAC2, TWOKAPPA_RE_K, TWOKAPPA_IM_K, LMBDA_RE_K, LMBDA_IM_K
	DOUBLE PRECISION N_RE_K, N_IM_K, TWOMU_RE_K, TWOMU_IM_K, TWOGAMMA_RE_K, TWOGAMMA_IM_K

    INTEGER N_FREQ,K,KK,IND1,IND2

    DOUBLE PRECISION FREQM(250),TWOKAPPA_RE(250),LMBDA_RE(250),N_RE(250),TWOMU_RE(250),TWOGAMMA_RE(250)
    DOUBLE PRECISION TWOKAPPA_IM(250),LMBDA_IM(250),N_IM(250),TWOMU_IM(250),TWOGAMMA_IM(250)
	DOUBLE PRECISION C_RE(6,6), C_IM(6,6)


	offset = 8
	n_freq = int(props(1))

C     In practice we're reading data each increment 
	do i=1,n_freq
		freqm(i)=props(1+offset+(i-1)*11)
		twokappa_Re(i)=props(2 + offset+(i-1)*11)
		twokappa_Im(i)=props(3 + offset+(i-1)*11)
		lmbda_Re(i)=props(4 + offset+(i-1)*11)
		lmbda_Im(i)=props(5 + offset+(i-1)*11)
		n_Re(i)=props(6 + offset +(i-1)*11)
		n_Im(i)=props(7 + offset +(i-1)*11)
		twomu_Re(i)=props(8 + offset +(i-1)*11)
		twomu_Im(i)=props(9 + offset +(i-1)*11)
		twogamma_Re(i)=props(10 + offset +(i-1)*11)
		twogamma_Im(i)=props(11 + offset +(i-1)*11)
	enddo
	

C     Frequence imposee
    freq = time(1)
	if ( freq .le. freqm(1) ) then 
          ind1 = 1
          ind2 = 1
          fac1 = 1.
          fac2 = 0.
    else if ( freq .ge. freqm(n_freq) ) then 
          ind1 = n_freq
          ind2 = n_freq
          fac1 = 1.
          fac2 = 0. 
    else 
          do nbr=1,n_freq-1
          
              if ( freq .ge. freqm(nbr) ) then 
                  if ( freq .lt. freqm(nbr+1) ) then 
                      ind1 = nbr
                      ind2 = nbr + 1
                      fac2 = ( freq - freqm(nbr))  &
                      / ( freqm(nbr+1) - freqm(nbr))
                      fac1 = 1. - fac2
                  end if
              end if
          end do
    end if

C     Generating the fourth order complex moduli tensor for the given frquency
	
	C_Re = 0.0d0
	C_Im = 0.0d0
	
	C_Re(1,1) = 0.5d0*(twokappa_Re(ind1) + twomu_Re(ind1))*fac1 + 0.5d0*(twokappa_Re(ind2) + twomu_Re(ind2))*fac2
	C_Re(2,2) = C_Re(1,1)
	C_Re(1,2) = 0.5d0*(twokappa_Re(ind1) - twomu_Re(ind1))*fac1 + 0.5d0*(twokappa_Re(ind2) - twomu_Re(ind2))*fac2
	C_Re(2,1) = C_Re(1,2)
	C_Re(1,3) = lmbda_Re(ind1)*fac1 + lmbda_Re(ind2)*fac2
	C_Re(2,3) = C_Re(1,3)
	C_Re(3,1) = C_Re(1,3)
	C_Re(3,2) = C_Re(1,3)
	C_Re(3,3) = n_Re(ind1)*fac1 + n_Re(ind2)*fac2
	C_Re(4,4) = 0.5d0*twomu_Re(ind1)*fac1 + 0.5d0*twomu_Re(ind2)*fac2
	C_Re(5,5) = 0.5d0*twogamma_Re(ind1)*fac1 + 0.5d0*twogamma_Re(ind2)*fac2
	C_Re(6,6) = C_Re(5,5)	
	C_Im(1,1) = 0.5d0*(twokappa_Im(ind1) + twomu_Im(ind1))*fac1 + 0.5d0*(twokappa_Im(ind2) + twomu_Im(ind2))*fac2
	C_Im(2,2) = C_Im(1,1)
	C_Im(1,2) = 0.5d0*(twokappa_Im(ind1) - twomu_Im(ind1))*fac1 + 0.5d0*(twokappa_Im(ind2) - twomu_Im(ind2))*fac2
	C_Im(2,1) = C_Im(1,2)
	C_Im(1,3) = lmbda_Im(ind1)*fac1 + lmbda_Im(ind2)*fac2
	C_Im(2,3) = C_Im(1,3)
	C_Im(3,1) = C_Im(1,3)
	C_Im(3,2) = C_Im(1,3)
	C_Im(3,3) = n_Im(ind1)*fac1 + n_Im(ind2)*fac2
	C_Im(4,4) = 0.5d0*twomu_Im(ind1)*fac1 + 0.5d0*twomu_Im(ind2)*fac2
	C_Im(5,5) = 0.5d0*twogamma_Im(ind1)*fac1 + 0.5d0*twogamma_Im(ind2)*fac2
	C_Im(6,6) = C_Im(5,5)
	

	!-------------------------------------------------------------------------------
	!-------------------------------------------------------------------------------
	!-------------------------------------------------------------------------------
C     Udating the fourth complex moduli tensor 
    ddsdde(:, :, 1)=C_Re
    ddsdde(:, :, 2)=C_Im
	
    END SUBROUTINE