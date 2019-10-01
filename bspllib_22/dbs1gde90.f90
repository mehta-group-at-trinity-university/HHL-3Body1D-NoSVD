! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! IMSL name:  dbs1gd (double precision version)
!
! purpose:    evaluate the derivative of a spline on a grid, given
!             its B-spline representation.
!
! usage:      call dbs1gd(ideriv, n, xvec, korder, xknot, ncoef,
!                         bscoef, value)
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


      use bspline

      implicit none

!
!        specifications for parameters
!

      integer    korder, ndata, nknot, nfgrid
      parameter  (korder=3, ndata=5, nknot=ndata+korder, nfgrid = 9)

!
!        specifications for local variables
!

      integer    i, ncoef
      double precision ans0(nfgrid), ans1(nfgrid), bscoef(ndata),       &
     &           fdata(ndata),                                          &
     &           x, xdata(ndata), xknot(nknot), xvec(nfgrid)

!
!        specifications for subroutines
!

      double precision df, f


      f(x)  = sqrt(x)
      df(x) = 0.5/sqrt(x)

!
!        set up interpolation points
!

      do i = 1, ndata
         xdata(i) = dble(i)/dble(ndata)
         fdata(i) = f(xdata(i))
      end do

      call dbsnak (ndata, xdata, korder, xknot)

!
!        interpolate
!
      call dbsint (ndata, xdata, fdata, korder, xknot, bscoef)
      write (6,99999)

!
!        print on a finer grid
!

      ncoef   = ndata
      xvec(1) = xdata(1)

      do i = 2, 2*ndata - 2, 2
         xvec(i)   = (xdata(i/2+1)+xdata(i/2))/2.0
         xvec(i+1) = xdata(i/2+1)
      end do

      call dbs1gd (0, 2*ndata-1, xvec, korder, xknot, ncoef, bscoef,    &
     &            ans0)
      call dbs1gd (1, 2*ndata-1, xvec, korder, xknot, ncoef, bscoef,    &
     &            ans1)

      do i = 1, 2*ndata - 1
         write (6,99998) xvec(i), ans0(i), f(xvec(i)) - ans0(i),        &
     &        ans1(i), df(xvec(i)) - ans1(i)
      end do

99998 format (' ', f6.4, 5x, f7.4, 5x, f8.4, 5x, f8.4, 5x, f8.4)
99999 format (6x, 'x', 8x, 's(x)', 7x, 'error', 8x, 's''(x)', 8x,       &
     &       'error', /)

      end

!       x        s(x)       error        s'(x)        error
!
!  0.2000      0.4472       0.0000       1.0423       0.0757
!  0.3000      0.5456       0.0021       0.9262      -0.0133
!  0.4000      0.6325       0.0000       0.8101      -0.0196
!  0.5000      0.7077      -0.0006       0.6940       0.0131
!  0.6000      0.7746       0.0000       0.6446       0.0009
!  0.7000      0.8366       0.0001       0.5952       0.0024
!  0.8000      0.8944       0.0000       0.5615      -0.0025
!  0.9000      0.9489      -0.0002       0.5279      -0.0008
!  1.0000      1.0000       0.0000       0.4942       0.0058
