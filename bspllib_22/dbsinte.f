c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c IMSL name:  dbsint (double precision version)
c
c purpose:    compute the spline interpolant, returning the B-spline
c             coefficients.
c
c usage:      call dbsint(ndata, xdata, fdata, korder, xknot, bscoef)
c
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


      implicit none

c
c        specifications for parameters
c

      integer    korder, ndata, nknot
      parameter  (korder=3, ndata=5, nknot=ndata+korder)

      integer    i, ncoef

      double precision bscoef(ndata), bt, f, fdata(ndata),
     &           x, xdata(ndata), xknot(nknot), xt, dbsval

      external dbsint, dbsnak, dbsval

c
c        define function
c

      f(x) = sqrt(x)

c
c        set up interpolation points
c

      do i = 1, ndata
         xdata(i) = dble(i-1)/dble(ndata-1)
         fdata(i) = f(xdata(i))
      end do

c
c        generate knot sequence
c

      call dbsnak(ndata, xdata, korder, xknot)

c
c       interpolate
c

      call dbsint (ndata, xdata, fdata, korder, xknot, bscoef)

c
c        write heading
c

      write (6,99999)

c
c        print on a finer grid
c

      ncoef = ndata
      xt    = xdata(1)

c
c        evaluate spline
c

      bt = dbsval(xt,korder,xknot,ncoef,bscoef)

      write (6,99998) xt, bt, f(xt) - bt

      do i = 2, ndata
         xt = (xdata(i-1)+xdata(i))/2.0

c
c           evaluate spline
c

         bt = dbsval(xt,korder,xknot,ncoef,bscoef)

         write (6,99998) xt, bt, f(xt) - bt

         xt = xdata(i)

c
c           evaluate spline
c

         bt = dbsval(xt,korder,xknot,ncoef,bscoef)
         write (6,99998) xt, bt, f(xt) - bt

      end do

99998 format (' ', f6.4, 15x, f8.4, 12x, f11.6)
99999 format (/, 6x, 'x', 19x, 's(x)', 18x, 'error', /)

      end

c
c       x                   s(x)                  error
c
c  0.0000                 0.0000               0.000000
c  0.1250                 0.2918               0.061781
c  0.2500                 0.5000               0.000000
c  0.3750                 0.6247              -0.012311
c  0.5000                 0.7071               0.000000
c  0.6250                 0.7886               0.002013
c  0.7500                 0.8660               0.000000
c  0.8750                 0.9365              -0.001092
c  1.0000                 1.0000               0.000000
