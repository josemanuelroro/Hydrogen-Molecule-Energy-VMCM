module mt19937
  implicit none
  integer,private :: N, N1, M, MATA, UMASK, LMASK, TMASKB, TMASKC
  parameter( &
             & N = 624, &
             & N1 = 625, &
             & M = 397, &
             & MATA = -1727483681, &
             & UMASK = -2147483, &
             & LMASK = 2147483647, &
             & TMASKB = -1658038656, &
             & TMASKC = -272236544 &
             & )
  integer,private :: mti = N1, mt(0:N-1), mag01(0:1) = (/0, MATA/)

  contains

  subroutine sgrnd(seed)
    integer,intent(in) :: seed


    mt(0) = iand(seed, -1)
    do mti = 1, N - 1
      mt(mti) = iand(69069 * mt(mti - 1), -1)
    end do
  end subroutine sgrnd

  real(8) function grnd()
    integer :: y, kk

    if(mti >= N) then

      if(mti == N + 1) then

        call sgrnd(4357)

      endif

      do kk = 0, N - M - 1
        y = ior(iand(mt(kk), UMASK), iand(mt(kk + 1), LMASK))
        mt(kk) = ieor(ieor(mt(kk + M), ishft(y, -1)), mag01(iand(y, 1)))
      end do

      do kk = N - M, N - 2
        y = ior(iand(mt(kk), UMASK), iand(mt(kk + 1), LMASK))
        mt(kk) = ieor(ieor(mt(kk + (M - N)), ishft(y, -1)), mag01(iand(y, 1)))
      end do

      y = ior(iand(mt(N - 1), UMASK), iand(mt(0), LMASK))
      mt(N - 1) = ieor(ieor(mt(M - 1), ishft(y, -1)), mag01(iand(y, 1)))
      mti = 0
    endif

    y = mt(mti)
    mti = mti + 1
    y = ieor(y, ishft(y, -11))
    y = ieor(y, iand(ishft(y, 7), TMASKB))
    y = ieor(y, iand(ishft(y, 15), TMASKC))
    y = ieor(y, ishft(y, -18))

    if(y < 0) then
      grnd = (dble(y) + 2.0d0 ** 32) / (2.0d0 ** 32)
    else
      grnd = dble(y) / (2.0d0 ** 32)
    endif
  end function grnd

end module mt19937