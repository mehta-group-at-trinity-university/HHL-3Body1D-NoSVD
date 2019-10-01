c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c IMSL name:  dbs2dr (double precision version)
c
c purpose:    evaluate the derivative of a two-dimensional
c             tensor-product spline, given its tensor-product
c             B-spline representation.
c
c usage:      dbs2dr(ixder, iyder, x, y, kxord, kyord, xknot, yknot,
c                    nxcoef, nycoef, bscoef)
c
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none

c
c        specifications for parameters
c


      integer    kxord, kyord, ldf, nxdata, nxknot, nydata, nyknot
      parameter  (kxord=5, kyord=3, nxdata=21, nydata=6, ldf=nxdata,
     &           nxknot=nxdata+kxord, nyknot=nydata+kyord)

      integer    i, j, nxcoef, nycoef

      double precision bscoef(nxdata,nydata), f, f21,
     &     fdata(ldf,nydata), s21, x, xdata(nxdata), dbs2dr,
     &     xknot(nxknot), y, ydata(nydata), yknot(nyknot)

c
c        define function and (2,1) derivative
c

      f(x,y)   = x*x*x*x + x*x*x*y*y
      f21(x,y) = 12.0*x*y

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
c        write heading
c

      write (6,99999)

c
c        print (2,1) derivative over a  grid of [0.0,1.0] x [0.0,1.0]
c        at 16 points.
c

      do i = 1, 4
         do j = 1, 4
            x   = dble(i-1)/3.0
            y   = dble(j-1)/3.0

c
c              evaluate spline
c

            s21 = dbs2dr(2,1,x,y,kxord,kyord,xknot,yknot,nxcoef,nycoef,
     &           bscoef)
            write (6,'(3f15.4, f15.6)') x, y, s21, f21(x,y) - s21
         end do
      end do

99999 format (39x, '(2,1)', /, 13x, 'x', 14x, 'y', 10x, 's    (x,y)',
     &        5x, 'error')
      end

c                                        (2,1)
c              x              y          s    (x,y)     error
c          0.0000         0.0000         0.0000       0.000000
c          0.0000         0.3333         0.0000       0.000000
c          0.0000         0.6667         0.0000       0.000000
c          0.0000         1.0000         0.0000       0.000001
c          0.3333         0.0000         0.0000       0.000000
c          0.3333         0.3333         1.3333       0.000002
c          0.3333         0.6667         2.6667      -0.000002
c          0.3333         1.0000         4.0000       0.000008
c          0.6667         0.0000         0.0000       0.000006
c          0.6667         0.3333         2.6667      -0.000011
c          0.6667         0.6667         5.3333       0.000028
c          0.6667         1.0000         8.0001      -0.000134
c          1.0000         0.0000        -0.0004       0.000439
c          1.0000         0.3333         4.0003      -0.000319
c          1.0000         0.6667         7.9996       0.000363
c          1.0000         1.0000        12.0005      -0.000458
