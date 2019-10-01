! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! IMSL name:  dbsnak (double precision version)
!
! purpose:    compute the `not-a-knot' spline knot sequence.
!
! usage:      call dbsnak(ndata, xdata, korder, xknot)
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      use bspline

      implicit none

!
!        specifications for parameters
!

      integer    kmax, kmin, ndata
      parameter  (kmax=8, kmin=3, ndata=20)

      integer    i, k, korder

      double precision  bscoef(ndata), dif, difmax,                     &
     &           fdata(ndata), ft, st, t, x, xdata(ndata),              &
     &           xknot(kmax+ndata), xt, f

!
!        define function and tau function
!

      f(x) = sin(10.0d0*x*x*x)
      t(x) = 1.0d0 - x*x

!
!        set up data
!

      do i = 1, ndata
         xt               = dble(i-1)/dble(ndata-1)
         xdata(ndata-i+1) = t(xt)
      end do

      xdata(1) = 0.0d0

      do i = 1, ndata
         fdata(i) = f(xdata(i))
      end do

!
!        write heading
!

      write (6,99999)

!
!        loop over different orders
!

      do k = kmin, kmax
         korder = k

!
!           generate knots
!

         call dbsnak(ndata, xdata, korder, xknot)

!
!           interpolate
!

         call dbsint(ndata, xdata, fdata, korder, xknot, bscoef)

         difmax = 0.0d0

         do i = 1, 100
            xt     = dble(i-1)/99.1d0

!
!                                  evaluate spline
!

            st     = dbsval(xt,korder,xknot,ndata,bscoef)
            ft     = f(xt)
            dif    = abs(ft-st)

!
!              compute maximum difference
!

            difmax = max(dif,difmax)

         end do

!
!           print maximum difference
!

         write (6,99998) korder, difmax

      end do

99998 format (' ', i3, 5x, f9.4)
99999 format (' korder', 5x, 'maximum difference', /)

      end

!  korder     maximum difference
!
!    3        0.0081
!    4        0.0026
!    5        0.0004
!    6        0.0008
!    7        0.0010
!    8        0.0004
