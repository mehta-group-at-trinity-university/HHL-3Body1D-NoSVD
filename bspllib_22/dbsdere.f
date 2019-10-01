c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c IMSL name:  dbsder (double precision version)
c
c purpose:    evaluate the derivative of a spline, given its B-spline
c             representation.
c
c usage:      dbsder(ideriv, x, korder, xknot, ncoef, bscoef)
c
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none

c
c        specifications for parameters
c

      integer    korder, ndata, nknot
      parameter  (korder=3, ndata=5, nknot=ndata+korder)

      integer    i, ncoef
      double precision  bscoef(ndata), dbsder, bt0, bt1, df, f,
     &           fdata(ndata),
     &           x, xdata(ndata), xknot(nknot), xt

      external   dbsder, dbsint, dbsnak

c
c       define function and derivative
c

      f(x)  = sqrt(x)
      df(x) = 0.5/sqrt(x)

c
c        set up interpolation points
c

      do i = 1, ndata
         xdata(i) = dble(i)/dble(ndata)
         fdata(i) = f(xdata(i))
      end do

c
c        generate knot sequence
c

      call dbsnak (ndata, xdata, korder, xknot)

c
c        interpolate
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

      bt0   = dbsder(0,xt,korder,xknot,ncoef,bscoef)
      bt1   = dbsder(1,xt,korder,xknot,ncoef,bscoef)

      write (6,99998) xt, bt0, f(xt) - bt0, bt1, df(xt) - bt1

      do i = 2, ndata
         xt  = (xdata(i-1)+xdata(i))/2.0

c
c           evaluate spline
c

         bt0 = dbsder(0,xt,korder,xknot,ncoef,bscoef)
         bt1 = dbsder(1,xt,korder,xknot,ncoef,bscoef)
         write (6,99998) xt, bt0, f(xt) - bt0, bt1, df(xt) - bt1
         xt  = xdata(i)

c
c           evaluate spline
c

         bt0 = dbsder(0,xt,korder,xknot,ncoef,bscoef)
         bt1 = dbsder(1,xt,korder,xknot,ncoef,bscoef)
         write (6,99998) xt, bt0, f(xt) - bt0, bt1, df(xt) - bt1
      end do

99998 format (' ', f6.4, 5x, f7.4, 3x, f10.6, 5x, f8.4, 3x, f10.6)
99999 format (6x, 'x', 8x, 's(x)', 7x, 'error', 8x, 's''(x)', 8x,
     &       'error', /)

      end

c       x        s(x)       error        s'(x)        error
c
c  0.2000      0.4472     0.000000       1.0423     0.075738
c  0.3000      0.5456     0.002084       0.9262    -0.013339
c  0.4000      0.6325     0.000000       0.8101    -0.019553
c  0.5000      0.7077    -0.000557       0.6940     0.013071
c  0.6000      0.7746     0.000000       0.6446     0.000869
c  0.7000      0.8366     0.000071       0.5952     0.002394
c  0.8000      0.8944     0.000000       0.5615    -0.002525
c  0.9000      0.9489    -0.000214       0.5279    -0.000818
c  1.0000      1.0000     0.000000       0.4942     0.005814
