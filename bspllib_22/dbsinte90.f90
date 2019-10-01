! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! IMSL name:  dbsint (double precision version)
!
! purpose:    compute the spline interpolant, returning the B-spline
!             coefficients.
!
! usage:      call dbsint(ndata, xdata, fdata, korder, xknot, bscoef)
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      use bspline

      implicit none

!
!        specifications for parameters
!

      integer    korder, ndata, nknot
      parameter  (korder=3, ndata=5, nknot=ndata+korder)

      integer    i, ncoef

      double precision bscoef(ndata), bt, f, fdata(ndata),              &
     &           x, xdata(ndata), xknot(nknot), xt

!
!        define function
!

      f(x) = sqrt(x)

!
!        set up interpolation points
!

      do i = 1, ndata
         xdata(i) = dble(i-1)/dble(ndata-1)
         fdata(i) = f(xdata(i))
      end do

!
!        generate knot sequence
!

      call dbsnak(ndata, xdata, korder, xknot)

!
!       interpolate
!

      call dbsint (ndata, xdata, fdata, korder, xknot, bscoef)

!
!        write heading
!

      write (6,99999)

!
!        print on a finer grid
!

      ncoef = ndata
      xt    = xdata(1)

!
!        evaluate spline
!

      bt = dbsval(xt,korder,xknot,ncoef,bscoef)

      write (6,99998) xt, bt, f(xt) - bt

      do i = 2, ndata
         xt = (xdata(i-1)+xdata(i))/2.0

!
!           evaluate spline
!

         bt = dbsval(xt,korder,xknot,ncoef,bscoef)

         write (6,99998) xt, bt, f(xt) - bt

         xt = xdata(i)

!
!           evaluate spline
!

         bt = dbsval(xt,korder,xknot,ncoef,bscoef)
         write (6,99998) xt, bt, f(xt) - bt

      end do

99998 format (' ', f6.4, 15x, f8.4, 12x, f11.6)
99999 format (/, 6x, 'x', 19x, 's(x)', 18x, 'error', /)

      end

!
!       x                   s(x)                  error
!
!  0.0000                 0.0000               0.000000
!  0.1250                 0.2918               0.061781
!  0.2500                 0.5000               0.000000
!  0.3750                 0.6247              -0.012311
!  0.5000                 0.7071               0.000000
!  0.6250                 0.7886               0.002013
!  0.7500                 0.8660               0.000000
!  0.8750                 0.9365              -0.001092
!  1.0000                 1.0000               0.000000
