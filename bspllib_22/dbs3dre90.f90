! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! IMSL name:  dbs3dr (double precision version)
!
! purpose:    evaluate the derivative of a three-dimensional
!             tensor-product spline, given its tensor-product
!             B-spline representation.
!
! usage:      dbs3dr(ixder, iyder, izder, x, y, z, kxord, kyord,
!                    kzord, xknot, yknot, zknot, nxcoef, nycoef, nzcoef,
!                    bscoef)
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      use bspline

      implicit none

!
!        specifications for parameters
!

      integer    kxord, kyord, kzord, ldf, mdf, nxdata, nxknot,         &
     &           nydata, nyknot, nzdata, nzknot
      parameter  (kxord=5, kyord=2, kzord=3, nxdata=21, nydata=6,       &
     &           nzdata=8, ldf=nxdata, mdf=nydata,                      &
     &           nxknot=nxdata+kxord, nyknot=nydata+kyord,              &
     &           nzknot=nzdata+kzord)

      integer    i, j, k, l, nxcoef, nycoef, nzcoef
      double precision bscoef(nxdata,nydata,nzdata), f, f201,           &
     &           fdata(ldf,mdf,nzdata), s201, x, xdata(nxdata),         &
     &           xknot(nxknot), y, ydata(nydata), yknot(nyknot), z,     &
     &           zdata(nzdata), zknot(nzknot)

!
!        define function and (2,0,1) derivative
!

      f(x,y,z)    = x*x*x*x + x*x*x*y*z*z*z
      f201(x,y,z) = 18.0*x*y*z

!
!        set up x-interpolation points
!

      do i = 1, nxdata
         xdata(i) = dble(i-11)/10.0
      end do

!
!        set up y-interpolation points
!

      do i = 1, nydata
         ydata(i) = dble(i-1)/dble(nydata-1)
      end do

!
!        set up z-interpolation points
!

      do i = 1, nzdata
         zdata(i) = dble(i-1)/dble(nzdata-1)
      end do

!
!        generate knots
!

      call dbsnak (nxdata, xdata, kxord, xknot)
      call dbsnak (nydata, ydata, kyord, yknot)
      call dbsnak (nzdata, zdata, kzord, zknot)

!
!        generate fdata
!

      do k = 1, nzdata
         do i = 1, nydata
            do j = 1, nxdata
               fdata(j,i,k) = f(xdata(j),ydata(i),zdata(k))
            end do
         end do
      end do

!
!        interpolate
!

      call dbs3in (nxdata, xdata, nydata, ydata, nzdata, zdata, fdata,  &
     &            ldf, mdf, kxord, kyord, kzord, xknot, yknot, zknot,   &
     &            bscoef)

      nxcoef = nxdata
      nycoef = nydata
      nzcoef = nzdata

!
!        write heading
!

      write (6,99999)
!
!        print over a grid of [-1.0,1.0] x [0.0,1.0] x [0.0,1.0]
!        at 32 points.
!
      do i = 1, 4
         do j = 1, 4
            do l = 1, 2
               x    = 2.0*(dble(i-1)/3.0) - 1.0
               y    = dble(j-1)/3.0
               z    = dble(l-1)

!
!                 evaluate spline
!

               s201 = dbs3dr(2,0,1,x,y,z,kxord,kyord,kzord,xknot,yknot, &
     &                zknot,nxcoef,nycoef,nzcoef,bscoef)
               write (6,'(3f12.4,2f12.6)') x, y, z, s201,               &
     &                f201(x,y,z) - s201
            end do
         end do
      end do

99999 format (38x, '(2,0,1)', /, 9x, 'x', 11x,                          &
     &       'y', 11x, 'z', 4x, 's     (x,y,z)    error')

      end

!                                       (2,0,1)
!          x           y           z    s     (x,y,z)    error
!      -1.0000      0.0000      0.0000   -0.000107    0.000107
!      -1.0000      0.0000      1.0000    0.000053   -0.000053
!      -1.0000      0.3333      0.0000    0.064051   -0.064051
!      -1.0000      0.3333      1.0000   -5.935941   -0.064059
!      -1.0000      0.6667      0.0000    0.127542   -0.127542
!      -1.0000      0.6667      1.0000  -11.873034   -0.126966
!      -1.0000      1.0000      0.0000    0.191166   -0.191166
!      -1.0000      1.0000      1.0000  -17.808527   -0.191473
!      -0.3333      0.0000      0.0000   -0.000002    0.000002
!      -0.3333      0.0000      1.0000    0.000000    0.000000
!      -0.3333      0.3333      0.0000    0.021228   -0.021228
!      -0.3333      0.3333      1.0000   -1.978768   -0.021232
!      -0.3333      0.6667      0.0000    0.042464   -0.042464
!      -0.3333      0.6667      1.0000   -3.957536   -0.042464
!      -0.3333      1.0000      0.0000    0.063700   -0.063700
!      -0.3333      1.0000      1.0000   -5.936305   -0.063694
!       0.3333      0.0000      0.0000   -0.000003    0.000003
!       0.3333      0.0000      1.0000    0.000000    0.000000
!       0.3333      0.3333      0.0000   -0.021229    0.021229
!       0.3333      0.3333      1.0000    1.978763    0.021238
!       0.3333      0.6667      0.0000   -0.042465    0.042465
!       0.3333      0.6667      1.0000    3.957539    0.042462
!       0.3333      1.0000      0.0000   -0.063700    0.063700
!       0.3333      1.0000      1.0000    5.936304    0.063697
!       1.0000      0.0000      0.0000   -0.000098    0.000098
!       1.0000      0.0000      1.0000    0.000053   -0.000053
!       1.0000      0.3333      0.0000   -0.063855    0.063855
!       1.0000      0.3333      1.0000    5.936146    0.063854
!       1.0000      0.6667      0.0000   -0.127631    0.127631
!       1.0000      0.6667      1.0000   11.873067    0.126933
!       1.0000      1.0000      0.0000   -0.191442    0.191442
!       1.0000      1.0000      1.0000   17.807940    0.192060
