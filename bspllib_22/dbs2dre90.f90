! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! IMSL name:  dbs2dr (double precision version)
!
! purpose:    evaluate the derivative of a two-dimensional
!             tensor-product spline, given its tensor-product
!             B-spline representation.
!
! usage:      dbs2dr(ixder, iyder, x, y, kxord, kyord, xknot, yknot,
!                    nxcoef, nycoef, bscoef)
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      use bspline

      implicit none

!
!        specifications for parameters
!

      integer    kxord, kyord, ldf, nxdata, nxknot, nydata, nyknot
      parameter  (kxord=5, kyord=3, nxdata=21, nydata=6, ldf=nxdata,    &
     &           nxknot=nxdata+kxord, nyknot=nydata+kyord)

      integer    i, j, nxcoef, nycoef

      double precision bscoef(nxdata,nydata), f, f21,                   &
     &     fdata(ldf,nydata), s21, x, xdata(nxdata),                    &
     &     xknot(nxknot), y, ydata(nydata), yknot(nyknot)

!
!        define function and (2,1) derivative
!

      f(x,y)   = x*x*x*x + x*x*x*y*y
      f21(x,y) = 12.0*x*y

!
!        set up interpolation points
!

      do i = 1, nxdata
         xdata(i) = dble(i-11)/10.0
      end do

!
!        generate knot sequence
!

      call dbsnak (nxdata, xdata, kxord, xknot)

!
!        set up interpolation points
!

      do i = 1, nydata
         ydata(i) = dble(i-1)/5.0
      end do

!
!        generate knot sequence
!

      call dbsnak (nydata, ydata, kyord, yknot)

!
!        generate fdata
!

      do i = 1, nydata
         do j = 1, nxdata
            fdata(j,i) = f(xdata(j),ydata(i))
         end do
      end do

!
!        interpolate
!

      call dbs2in (nxdata, xdata, nydata, ydata, fdata, ldf, kxord,     &
     &            kyord, xknot, yknot, bscoef)

      nxcoef = nxdata
      nycoef = nydata

!
!        write heading
!

      write (6,99999)

!
!        print (2,1) derivative over a grid of [0.0,1.0] x [0.0,1.0]
!        at 16 points.
!

      do i = 1, 4
         do j = 1, 4
            x   = dble(i-1)/3.0
            y   = dble(j-1)/3.0

!
!              evaluate spline
!

            s21 = dbs2dr(2,1,x,y,kxord,kyord,xknot,yknot,nxcoef,nycoef, &
     &           bscoef)
            write (6,'(3f15.4, f15.6)') x, y, s21, f21(x,y) - s21
         end do
      end do

99999 format (39x, '(2,1)', /, 13x, 'x', 14x, 'y', 10x, 's    (x,y)',   &
     &        5x, 'error')

      end

!                                        (2,1)
!              x              y          s    (x,y)     error
!          0.0000         0.0000         0.0000       0.000000
!          0.0000         0.3333         0.0000       0.000000
!          0.0000         0.6667         0.0000       0.000000
!          0.0000         1.0000         0.0000       0.000001
!          0.3333         0.0000         0.0000       0.000000
!          0.3333         0.3333         1.3333       0.000002
!          0.3333         0.6667         2.6667      -0.000002
!          0.3333         1.0000         4.0000       0.000008
!          0.6667         0.0000         0.0000       0.000006
!          0.6667         0.3333         2.6667      -0.000011
!          0.6667         0.6667         5.3333       0.000028
!          0.6667         1.0000         8.0001      -0.000134
!          1.0000         0.0000        -0.0004       0.000439
!          1.0000         0.3333         4.0003      -0.000319
!          1.0000         0.6667         7.9996       0.000363
!          1.0000         1.0000        12.0005      -0.000458
