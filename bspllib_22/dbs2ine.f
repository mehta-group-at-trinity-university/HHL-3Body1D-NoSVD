c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c IMSL name:  dbs2in (double precision version)
c
c purpose:    compute a two-dimensional tensor-product spline
c             interpolant, returning the tensor-product B-spline
c             coefficients.
c
c usage:      call dbs2in(nxdata, xdata, nydata, ydata, fdata, ldf,
c                         kxord, kyord, xknot, yknot, bscoef)
c
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none

c
c        specifications for parameters
c

      integer    kxord, kyord, ldf, nxdata, nxknot, nxvec, nydata,
     &           nyknot, nyvec
      parameter  (kxord=5, kyord=2, nxdata=21, nxvec=4, nydata=6,
     &           nyvec=4, ldf=nxdata, nxknot=nxdata+kxord,
     &           nyknot=nydata+kyord)

      integer    i, j, nxcoef, nycoef

      double precision bscoef(nxdata,nydata), f, fdata(ldf,nydata),
     &           value(nxvec,nyvec), x, xdata(nxdata), xknot(nxknot),
     &           xvec(nxvec), y, ydata(nydata), yknot(nyknot),
     &           yvec(nyvec)

c
c        define function
c

      f(x,y) = x*x*x + x*y

c
c        set up interpolation points
c

      do i = 1, nxdata
         xdata(i) = dble(i-11)/10.0
      end do

c
c        generate knot sequence
c

      call dbsnak (nxdata, xdata, kxord, xknot)

c
c        set up interpolation points
c

      do i = 1, nydata
         ydata(i) = dble(i-1)/5.0
      end do

c
c        generate knot sequence
c

      call dbsnak (nydata, ydata, kyord, yknot)

c
c        generate fdata
c

      do i = 1, nydata
         do j = 1, nxdata
            fdata(j,i) = f(xdata(j),ydata(i))
         end do
      end do

c
c        interpolate
c

      call dbs2in (nxdata, xdata, nydata, ydata, fdata, ldf, kxord,
     &            kyord, xknot, yknot, bscoef)

      nxcoef = nxdata
      nycoef = nydata

c
c                                  write heading
c

      write (6,99999)

c
c        print over a grid of
c        [0.0,1.0] x [0.0,1.0] at 16 points.
c

      do i = 1, nxvec
         xvec(i) = dble(i-1)/3.0
      end do

      do i = 1, nyvec
         yvec(i) = dble(i-1)/3.0
      end do

c
c        evaluate spline
c

      call dbs2gd (0, 0, nxvec, xvec, nyvec, yvec, kxord, kyord, xknot,
     &     yknot, nxcoef, nycoef, bscoef, value, nxvec)

      do i = 1, nxvec
         do j = 1, nyvec
            write (6,'(3f15.4,f15.6)') xvec(i), yvec(j),
     &                                   value(i,j),
     &                                   (f(xvec(i),yvec(j))-
     &                                   value(i,j))
         end do
      end do

99999 format (13x, 'x', 14x, 'y', 10x, 's(x,y)', 9x, 'error')
      end

c              x              y          s(x,y)         error
c          0.0000         0.0000         0.0000       0.000000
c          0.0000         0.3333         0.0000       0.000000
c          0.0000         0.6667         0.0000       0.000000
c          0.0000         1.0000         0.0000       0.000000
c          0.3333         0.0000         0.0370       0.000000
c          0.3333         0.3333         0.1481       0.000000
c          0.3333         0.6667         0.2593       0.000000
c          0.3333         1.0000         0.3704       0.000000
c          0.6667         0.0000         0.2963       0.000000
c          0.6667         0.3333         0.5185       0.000000
c          0.6667         0.6667         0.7407       0.000000
c          0.6667         1.0000         0.9630       0.000000
c          1.0000         0.0000         1.0000       0.000000
c          1.0000         0.3333         1.3333       0.000000
c          1.0000         0.6667         1.6667       0.000000
c          1.0000         1.0000         2.0000       0.000000
