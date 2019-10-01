c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c IMSL name:  dbsnak (double precision version)
c
c purpose:    compute the `not-a-knot' spline knot sequence.
c
c usage:      call dbsnak(ndata, xdata, korder, xknot)
c
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


      implicit none

c
c        specifications for parameters
c

      integer    kmax, kmin, ndata
      parameter  (kmax=8, kmin=3, ndata=20)

      integer    i, k, korder

      double precision  bscoef(ndata), dbsval, dif, difmax,
     &           fdata(ndata), ft, st, t, x, xdata(ndata),
     &           xknot(kmax+ndata), xt, f

      external   dbsint, dbsnak, dbsval

c
c        define function and tau function
c

      f(x) = sin(10.0d0*x*x*x)
      t(x) = 1.0d0 - x*x

c
c        set up data
c

      do i = 1, ndata
         xt               = dble(i-1)/dble(ndata-1)
         xdata(ndata-i+1) = t(xt)
      end do

      xdata(1) = 0.0d0

      do i = 1, ndata
         fdata(i) = f(xdata(i))
      end do

c
c        write heading
c

      write (6,99999)

c
c        loop over different orders
c

      do k = kmin, kmax
         korder = k

c
c           generate knots
c

         call dbsnak(ndata, xdata, korder, xknot)

c
c           interpolate
c

         call dbsint(ndata, xdata, fdata, korder, xknot, bscoef)

         difmax = 0.0d0

         do i = 1, 100
            xt     = dble(i-1)/99.1d0

c
c                                  evaluate spline
c

            st     = dbsval(xt,korder,xknot,ndata,bscoef)
            ft     = f(xt)
            dif    = abs(ft-st)

c
c              compute maximum difference
c

            difmax = max(dif,difmax)

         end do

c
c           print maximum difference
c

         write (6,99998) korder, difmax

      end do

99998 format (' ', i3, 5x, f9.4)
99999 format (' korder', 5x, 'maximum difference', /)

      end

c  korder     maximum difference
c
c    3        0.0081
c    4        0.0026
c    5        0.0004
c    6        0.0008
c    7        0.0010
c    8        0.0004
