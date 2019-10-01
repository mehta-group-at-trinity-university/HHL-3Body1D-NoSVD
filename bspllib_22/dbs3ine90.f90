! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! IMSL name:  dbs3in (double precision version)
!
! purpose:    compute a three-dimensional tensor-product spline
!             interpolant, returning the tensor-product B-spline
!             coefficients.
!
! usage:      call dbs3in(nxdata, xdata, nydata, ydata, nzdata,
!                         zdata, fdata, ldf, mdf, kxord, kyord,
!                         kzord, xknot, yknot, zknot, bscoef)
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      use bspline

      implicit none

!
!        specifications for parameters
!

      integer    kxord, kyord, kzord, ldf, mdf, nxdata, nxknot, nxvec,  &
     &           nydata, nyknot, nyvec, nzdata, nzknot, nzvec
      parameter  (kxord=5, kyord=2, kzord=3, nxdata=21, nxvec=4,        &
     &           nydata=6, nyvec=4, nzdata=8, nzvec=2, ldf=nxdata,      &
     &           mdf=nydata, nxknot=nxdata+kxord, nyknot=nydata+kyord,  &
     &           nzknot=nzdata+kzord)

      integer    i, j, k, nxcoef, nycoef, nzcoef
      double precision  bscoef(nxdata,nydata,nzdata), f,                &
     &           fdata(ldf,mdf,nzdata), value(nxvec,nyvec,nzvec)        &
     &           , x, xdata(nxdata), xknot(nxknot), xvec(nxvec), y,     &
     &           ydata(nydata), yknot(nyknot), yvec(nyvec), z,          &
     &           zdata(nzdata), zknot(nzknot), zvec(nzvec)

!
!        define function.
!

      f(x,y,z) = x*x*x + x*y*z

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

!        print over a grid of [-1.0,1.0] x [0.0,1.0] x [0.0,1.0]
!        at 32 points.

      do i = 1, nxvec
         xvec(i) = 2.0*(dble(i-1)/3.0) - 1.0
      end do

      do i = 1, nyvec
         yvec(i) = dble(i-1)/3.0
      end do

      do i = 1, nzvec
         zvec(i) = dble(i-1)
      end do

!
!        call the evaluation routine.
!

      call dbs3gd (0, 0, 0, nxvec, xvec, nyvec, yvec, nzvec, zvec,      &
     &            kxord, kyord, kzord, xknot, yknot, zknot, nxcoef,     &
     &            nycoef, nzcoef, bscoef, value, nxvec, nyvec)

      do i = 1, nxvec
         do j = 1, nyvec
            do k = 1, nzvec
               write (6,'(4f13.4, f13.6)') xvec(i), yvec(k),            &
     &                                       zvec(k), value(i,j,k),     &
     &                                       f(xvec(i),yvec(j),zvec(k)) &
     &                                        - value(i,j,k)
            end do
         end do
      end do

99999 format (10x, 'x', 11x, 'y', 10x, 'z', 10x, 's(x,y,z)', 7x,        &
     &       'error')

      end

!           x           y          z          s(x,y,z)       error
!       -1.0000       0.0000       0.0000      -1.0000     0.000000
!       -1.0000       0.3333       1.0000      -1.0000     0.000000
!       -1.0000       0.0000       0.0000      -1.0000     0.000000
!       -1.0000       0.3333       1.0000      -1.3333     0.000000
!       -1.0000       0.0000       0.0000      -1.0000     0.000000
!       -1.0000       0.3333       1.0000      -1.6667     0.000000
!       -1.0000       0.0000       0.0000      -1.0000     0.000000
!       -1.0000       0.3333       1.0000      -2.0000     0.000000
!       -0.3333       0.0000       0.0000      -0.0370     0.000000
!       -0.3333       0.3333       1.0000      -0.0370     0.000000
!       -0.3333       0.0000       0.0000      -0.0370     0.000000
!       -0.3333       0.3333       1.0000      -0.1481     0.000000
!       -0.3333       0.0000       0.0000      -0.0370     0.000000
!       -0.3333       0.3333       1.0000      -0.2593     0.000000
!       -0.3333       0.0000       0.0000      -0.0370     0.000000
!       -0.3333       0.3333       1.0000      -0.3704     0.000000
!        0.3333       0.0000       0.0000       0.0370     0.000000
!        0.3333       0.3333       1.0000       0.0370     0.000000
!        0.3333       0.0000       0.0000       0.0370     0.000000
!        0.3333       0.3333       1.0000       0.1481     0.000000
!        0.3333       0.0000       0.0000       0.0370     0.000000
!        0.3333       0.3333       1.0000       0.2593     0.000000
!        0.3333       0.0000       0.0000       0.0370     0.000000
!        0.3333       0.3333       1.0000       0.3704     0.000000
!        1.0000       0.0000       0.0000       1.0000     0.000000
!        1.0000       0.3333       1.0000       1.0000     0.000000
!        1.0000       0.0000       0.0000       1.0000     0.000000
!        1.0000       0.3333       1.0000       1.3333     0.000000
!        1.0000       0.0000       0.0000       1.0000     0.000000
!        1.0000       0.3333       1.0000       1.6667     0.000000
!        1.0000       0.0000       0.0000       1.0000     0.000000
!        1.0000       0.3333       1.0000       2.0000     0.000000
