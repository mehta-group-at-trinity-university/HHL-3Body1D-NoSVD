c     234567890
      program HHL1DHyperspherical
      implicit none
      integer LegPoints,xNumPoints
      integer NumStates,PsiFlag,Order,Left,Right
      integer RSteps,CouplingFlag,CalcNewBasisFunc
      double precision alpha,mass,Shift,Shift2,NumStateInc,m1,m2,m3,phi23,phi13,phi12,mgamma
      double precision RLeft,RRight,RDerivDelt,DD,L
      DOUBLE PRECISION RFirst,RLast,XFirst,XLast,StepX
      double precision xMin,xMax
      double precision, allocatable :: R(:)
      double precision, allocatable :: xPoints(:)

      logical, allocatable :: Select(:)

      integer iparam(11),ncv,info
      integer i,j,k,iR,NumFirst,NumBound
      integer LeadDim,MatrixDim,HalfBandWidth
      integer xDim
      integer, allocatable :: iwork(:)
      integer, allocatable :: xBounds(:)
      double precision Tol,RChange
      double precision TotalMemory
      double precision mu, mu12, mu123, r0diatom, dDiatom, etaOVERpi, Pi

      double precision, allocatable :: LUFac(:,:),workl(:)
      double precision, allocatable :: workd(:),Residuals(:)
      double precision, allocatable :: xLeg(:),wLeg(:)
      double precision, allocatable :: u(:,:,:),uxx(:,:,:)
      double precision, allocatable :: S(:,:),H(:,:)
      double precision, allocatable :: lPsi(:,:),mPsi(:,:),rPsi(:,:),
     >     Energies(:,:)
      double precision, allocatable :: P(:,:),Q(:,:),dP(:,:)
      double precision ur(1:50000),acoef,bcoef,diff
      double precision sec,time,Rinitial,secp,timep,Rvalue
      character*64 LegendreFile
      common /Rvalue/ Rvalue      

c     read in number of energies and states to print
      read(5,*)
      read(5,*) NumStates,PsiFlag,CouplingFlag
      write(6,*) NumStates,PsiFlag,CouplingFlag
      
c     read in Gauss-Legendre info
      read(5,*)
      read(5,*)
      read(5,1002) LegendreFile
      write(6,1002) LegendreFile
      read(5,*)
      read(5,*)
      read(5,*) LegPoints
      write(6,*) LegPoints,' LegPoints'
      
c     read in boundary conditions
      read(5,*)
      read(5,*)
      read(5,*) Shift,Shift2,Order,Left,Right
      print*, 'Shift,Shift2, Order, Left, Right'
      print*, Shift,Shift2,Order,Left,Right

c     read in potential parameters
      read(5,*)
      read(5,*)
      read(5,*) alpha,m1,m2,m3,DD,L
      write(6,*) alpha,m1,m2,m3,DD,L

      mu12=m1*m2/(m1+m2)
      mu123=(m1+m2)*m3/(m1+m2+m3)
      mu=dsqrt(mu12*mu123)
      Pi=dacos(-1.d0)
      write(6,*) 'Pi=',Pi, 'mu12 = ', mu12, 'mu123 = ',mu123, 'mu = ', mu
      mgamma = mu/m1
      phi12=Pi/2
      phi23=datan(mgamma)
c     read in grid information
      read(5,*)
      read(5,*)
      read(5,*) xNumPoints,xMin,xMax
      write(6,*) xNumPoints,xMin,xMax

      read(5,*)
      read(5,*)
      read(5,*) RSteps,RDerivDelt,RFirst,RLast
      write(6,*) RSteps,RDerivDelt,RFirst,RLast

c     c	XFirst=dsqrt(RFirst)
c     c	XLast=dsqrt(RLast)
c     c	XFirst=RFirst**(1.d0/3.d0)
c     c	XLast=RLast**(1.d0/3.d0)

c     SET UP A LOG-GRID IN THE HYPERRADIUS
      XFirst = dlog10(RFirst)
      XLast = dlog10(RLast)
      StepX=(XLast-XFirst)/(RSteps-1.d0)
      
      allocate(R(RSteps))
      do i = 1,RSteps
c     read(5,*) R(i)
c     R(i)= (XFirst+(i-1)*StepX)**3
         R(i)= 10.d0**(XFirst+(i-1)*StepX)
      enddo

c      if (mod(xNumPoints,2) .ne. 0) then
c         write(6,*) 'xNumPoints not divisible by 2'
c         xNumPoints = (xNumPoints/2)*2
c         write(6,*) '   truncated to ',xNumPoints
c      endif

      allocate(xLeg(LegPoints),wLeg(LegPoints))
      call GetGaussFactors(LegendreFile,LegPoints,xLeg,wLeg)

      xDim = xNumPoints+Order-3
      if (Left .eq. 2) xDim = xDim + 1
      if (Right .eq. 2) xDim = xDim + 1

      MatrixDim = xDim
      HalfBandWidth = Order
      LeadDim = 3*HalfBandWidth+1

      TotalMemory = 2.0d0*(HalfBandWidth+1)*MatrixDim ! S, H
      TotalMemory = TotalMemory + 2.0d0*LegPoints*(Order+2)*xDim ! x splines
      TotalMemory = TotalMemory + 2.0d0*NumStates*NumStates ! P and Q matrices
      TotalMemory = TotalMemory + LeadDim*MatrixDim ! LUFac
      TotalMemory = TotalMemory + 4.0d0*NumStates*MatrixDim ! channel functions
      TotalMemory = TotalMemory + 4*xDim*xDim ! (CalcHamiltonian)
      TotalMemory = TotalMemory + LegPoints**2*xNumPoints ! (CalcHamiltonian)
      TotalMemory = 8.0d0*TotalMemory/(1024.0d0*1024.0d0)

      write(6,*)
      write(6,*) 'MatrixDim ',MatrixDim
      write(6,*) 'HalfBandWidth ',HalfBandWidth
      write(6,*) 'Approximate peak memory usage (in Mb) ',TotalMemory
      write(6,*)

      allocate(xPoints(xNumPoints))
      allocate(xBounds(xNumPoints+2*Order))
      allocate(u(LegPoints,xNumPoints,xDim),uxx(LegPoints,xNumPoints,
     >     xDim))
      allocate(S(HalfBandWidth+1,MatrixDim),H(HalfBandWidth+1,
     >     MatrixDim))
      allocate(P(NumStates,NumStates),Q(NumStates,NumStates),
     >     dP(NumStates,NumStates))

      ncv = 2*NumStates
      LeadDim = 3*HalfBandWidth+1
      allocate(iwork(MatrixDim))
      allocate(Select(ncv))
      allocate(LUFac(LeadDim,MatrixDim))
      allocate(workl(ncv*ncv+8*ncv))
      allocate(workd(3*MatrixDim))
      allocate(lPsi(MatrixDim,ncv),mPsi(MatrixDim,ncv),
     >     rPsi(MatrixDim,ncv))
      allocate(Residuals(MatrixDim))
      allocate(Energies(ncv,2))
      info = 0
      iR=1
      CalcNewBasisFunc=1
      Tol=1e-20

      NumBound=0

      RChange=100.d0
      do iR = 1,RSteps
         NumFirst=NumStates
         if (R(iR).gt.RChange) then
            NumFirst=NumBound
         endif
         NumStateInc=NumStates-NumFirst
c----------------------------------------------------------------------------------------
c     must move this block inside the loop over iR if the grid is adaptive
         
         print*, 'calling GridMaker'
c     call GridMaker(mu,R(iR),2.0d0, xNumPoints,xMin,xMax,xPoints,CalcNewBasisFunc)
         call GridMakerHHL(mu,mu12,mu123,phi23,R(iR),2.0d0, xNumPoints,xMin,xMax,xPoints,CalcNewBasisFunc)
         if(CalcNewBasisFunc.eq.1) then
            print*, 'done... Calculating Basis functions'
            call CalcBasisFuncs(Left,Right,Order,xPoints,LegPoints,xLeg,
     >           xDim,xBounds,xNumPoints,0,u)
            call CalcBasisFuncs(Left,Right,Order,xPoints,LegPoints,xLeg,
     >           xDim,xBounds,xNumPoints,2,uxx)
         endif
         print*, 'done... Calculating overlap matrix'
c     must move this block inside the loop if the grid is adaptive
c----------------------------------------------------------------------------------------
         call CalcOverlap(Order,xPoints,LegPoints,xLeg,wLeg,xDim,
     >        xNumPoints,u,xBounds,HalfBandWidth,S)

         write(6,*) 'Calculating Hamiltonian with R = ', R(iR)
         call CalcHamiltonian(alpha,R(iR),mu,mgamma,DD,L,Order,xPoints,
     >        LegPoints,xLeg,wLeg,xDim,
     >        xNumPoints,u,uxx,xBounds,
     >        HalfBandWidth,H)
         
         write(6,*) 'Calling MyLargeDsband'
         call MyLargeDsband(NumFirst,Shift2,NumStateInc,Energies,mPsi,
     >        Shift,MatrixDim,H,S,LUFac,LeadDim,HalfBandWidth,NumStates)
         write(6,*) 'done...'

         if (CouplingFlag .ne. 0) then
 
            if(iR.gt.1) then
               write(6,*) 'Calling FixPhase'
               call FixPhase(NumStates,HalfBandWidth,MatrixDim,S,NumStates,rPsi,mPsi)
            endif
         endif
         
c     call MyDsband(Select,Energies,mPsi,MatrixDim,Shift,MatrixDim,
c     >        H,S,HalfBandWidth+1,LUFac,LeadDim,HalfBandWidth,NumStates,
c     >        Tol,Residuals,ncv,mPsi,MatrixDim,iparam,workd,workl,
c     >        ncv*ncv+8*ncv,iwork,info)
c     if (CouplingFlag .ne. 0) call FixPhase(NumStates,HalfBandWidth,
c     >        MatrixDim,S,ncv,rPsi,mPsi)
         
c     call CalcEigenErrors(info,iparam,MatrixDim,H,HalfBandWidth+1,S,
c     >        HalfBandWidth,NumStates,mPsi,Energies,ncv)
         write(6,*) 'writing the energies'
         write(200,20) R(iR),(Energies(i,1), i = 1,min(NumStates,
     >        iparam(5)))
         
         write(6,*)
         write(6,*) 'RMid = ', R(iR)
         do i = 1,min(NumStates,iparam(5))
            write(6,*) 'Energy(',i,') = ',Energies(i,1),'  Error = ', Energies(i,2)
         enddo

         if (CouplingFlag .ne. 0) then
            
            RLeft = 0.99d0*R(iR) !R(iR)-RDerivDelt
            write(6,*) 'Calculating Hamiltonian'
            call CalcHamiltonian(alpha,RLeft,mu,mgamma,DD,L,Order,xPoints,
     >           LegPoints,xLeg,wLeg,xDim,xNumPoints,u,
     >           uxx,xBounds,HalfBandWidth,H)
            write(6,*) 'Calling MyLargeDsband'
            call MyLargeDsband(NumFirst,Shift2,NumStateInc,Energies,lPsi,
     >           Shift,MatrixDim,H,S,LUFac,LeadDim,HalfBandWidth,NumStates)
c     call MyDsband(Select,Energies,lPsi,MatrixDim,Shift,MatrixDim,
c     >           H,S,HalfBandWidth+1,LUFac,LeadDim,HalfBandWidth,
c     >           NumStates,Tol,Residuals,ncv,lPsi,MatrixDim,iparam,workd,workl,
c     >           ncv*ncv+8*ncv,iwork,info)
            if (iR .gt. 1) call FixPhase(NumStates,HalfBandWidth,
     >           MatrixDim,S,ncv,mPsi,lPsi)
c     call CalcEigenErrors(info,iparam,MatrixDim,H,
c     >           HalfBandWidth+1,
c     >           S,HalfBandWidth,NumStates,lPsi,Energies,ncv)
c     IF(R(iR).GT. 20.d0) Shift = 1.05d0*Energies(1,1)
c            IF(iR.GT.1) Shift = 1.05d0*Energies(1,1)
c     do i = 1,min(NumStates,iparam(5))
c     Energies(i,1) = Energies(i,1) + 1.875d0/(mu*RLeft*RLeft)
c     enddo
c     write(100,10) RLeft,(Energies(i,1), i = 1,NumStates)
            write(6,*)
            write(6,*) 'RLeft = ', RLeft
            do i = 1,min(NumStates,iparam(5))
               write(6,*) 'Energy(',i,') = ',Energies(i,1),'  Error = ', Energies(i,2)
            enddo

            RRight = 1.01d0*R(iR) !R(iR)+RDerivDelt
            call CalcHamiltonian(alpha,RRight,mu,mgamma,DD,L,Order,xPoints,
     >           LegPoints,xLeg,wLeg,xDim,
     >           xNumPoints,u,uxx,xBounds,
     >           HalfBandWidth,H)
            call MyLargeDsband(NumFirst,Shift2,NumStateInc,Energies,rPsi,
     >           Shift,MatrixDim,H,S,LUFac,LeadDim,HalfBandWidth,NumStates)
            call FixPhase(NumStates,HalfBandWidth,MatrixDim,S,NumStates,lPsi,rPsi)

c            call MyDsband(Select,Energies,rPsi,MatrixDim,Shift,
c     >           MatrixDim,
c     >           H,S,HalfBandWidth+1,LUFac,LeadDim,HalfBandWidth,
c     >           NumStates,Tol,
c     >           Residuals,ncv,rPsi,MatrixDim,iparam,workd,workl,
c     >           ncv*ncv+8*ncv,iwork,info)
c            call FixPhase(NumStates,HalfBandWidth,MatrixDim,S,ncv,
c     >           lPsi,rPsi)
c            call CalcEigenErrors(info,iparam,MatrixDim,H,
c     >           HalfBandWidth+1,
c     >           S,HalfBandWidth,NumStates,rPsi,Energies,ncv)
c            IF(R(iR).GT.2.2d0) Shift = 0.93d0*Energies(1,1)
c     do i = 1,min(NumStates,iparam(5))
c     enddo
c     write(100,10) RRight,(Energies(i,1), i = 1,NumStates)
            write(6,*)
            write(6,*) 'RRight = ', RRight
            do i = 1,min(NumStates,iparam(5))
               write(6,*) 'Energy(',i,') = ',Energies(i,1),'  Error = ', Energies(i,2)
            enddo
            
         endif

         
c     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c     Adjusting Shift
         ur(iR) = dreal(Energies(1,1))
         if (iR.ge.2) then
            if ((R(iR+1).lt.Rchange).and.((ur(iR)-ur(iR-1)).gt.0.d0)) then
               Shift = dreal(Energies(1,1))
               if (dreal(Energies(1,1)).gt.0.d0) Shift = 0.1d0*Shift
               if (dreal(Energies(1,1)).lt.0.d0) Shift = 10.d0*Shift
            endif
            if (R(iR+1).gt.Rchange) then
c     First group of energies
               Shift = dreal(Energies(1,1))
               if (dreal(Energies(1,1)).gt.0.d0) then 
                  if (dreal(Energies(1,1)).lt.1.d-10) then
                     Shift = -1.d-9
                  else
                     Shift = -10.d0*Shift  
c     if (Shift.lt.1.d-10) Shift = -1.d-10
                  endif
               endif   
               if (dreal(Energies(1,1)).lt.0.d0) Shift = 10.d0*Shift

c     if ((dreal(Energies(1,1)).lt.0.d0)) then
c     if ((R(iR).ge.300.d0).and.(ncount.le.10)) then
c     Shift = -1.d-9
c     ncount = ncount+1
c     else
c     Shift = 10.d0*Shift
c     endif
c     endif
               
c     Second group of energies
               Shift2 = dreal(Energies(NumBound+1,1))
               if (dreal(Energies(NumBound+1,1)).gt.0.d0) Shift2 = 0.01d0*Shift2
               if (dreal(Energies(NumBound+1,1)).lt.0.d0) Shift2 = 100.0d0*Shift2
               if (dreal(Shift2).le.dreal(Energies(NumBound,1))) Shift2 = dreal(Energies(NumBound,1))
            endif
         endif   
c     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         
         call cpu_time(timep)
         secp = timep 
         
         
         if (CouplingFlag .ne. 0) then
            call CalcCoupling(NumStates,HalfBandWidth,MatrixDim,
     >           RDerivDelt,lPsi,mPsi,rPsi,S,P,Q,dP)

            write(101,*) R(iR)
            write(102,*) R(iR)
            write(103,*) R(iR)
            do i = 1,min(NumStates,iparam(5))
               write(101,20) (P(i,j), j = 1,min(NumStates,iparam(5)))
               write(102,20) (Q(i,j), j = 1,min(NumStates,iparam(5)))
               write(103,20) (dP(i,j), j = 1,min(NumStates,iparam(5)))
            enddo
         endif
c         write(400,20) R(iR),R(iR)**3.0d0*(Energies(2,1)-Q(2,2))
         if (PsiFlag .ne. 0) then
            do i = 1,xNumPoints
               write(97,*) xPoints(i)
            enddo
            do i = 1,MatrixDim
               write(999+iR,20) (mPsi(i,j), j = 1,NumStates)
            enddo
            close(unit=999+iR)
         endif

      enddo

      deallocate(S,H)
      deallocate(Energies)
      deallocate(iwork)
      deallocate(Select)
      deallocate(LUFac)
      deallocate(workl)
      deallocate(workd)
      deallocate(lPsi,mPsi,rPsi)
      deallocate(Residuals)
      deallocate(P,Q,dP)
      deallocate(xPoints)
      deallocate(xLeg,wLeg)
      deallocate(xBounds)
      deallocate(u,uxx)

      deallocate(R)

 10   format(1P,100e25.15)
 20   format(1P,100e16.8)
 1002 format(a64)

      stop
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine CalcOverlap(Order,xPoints,LegPoints,xLeg,wLeg,xDim,
     >     xNumPoints,u,xBounds,HalfBandWidth,S)
      implicit none
      integer Order,LegPoints,xDim,xNumPoints,xBounds(xNumPoints+2*Order),HalfBandWidth
      double precision xPoints(*),xLeg(*),wLeg(*)
      double precision S(HalfBandWidth+1,xDim)
      double precision u(LegPoints,xNumPoints,xDim)

      integer ix,ixp,kx,lx
      integer i1,i1p
      integer Row,NewRow,Col
      integer, allocatable :: kxMin(:,:),kxMax(:,:)
      double precision a,b,m
      double precision xTempS
      double precision ax,bx
      double precision, allocatable :: xIntScale(:),xS(:,:)
   
      allocate(xIntScale(xNumPoints),xS(xDim,xDim))
      allocate(kxMin(xDim,xDim),kxMax(xDim,xDim))

      S = 0.0d0

      do kx = 1,xNumPoints-1
         ax = xPoints(kx)
         bx = xPoints(kx+1)
         xIntScale(kx) = 0.5d0*(bx-ax)
      enddo

c      do ix=1,xNumPoints+2*Order
c         print*, ix, xBounds(ix)
c      enddo

      do ix = 1,xDim
         do ixp = 1,xDim
            kxMin(ixp,ix) = max(xBounds(ix),xBounds(ixp))
            kxMax(ixp,ix) = min(xBounds(ix+Order+1),xBounds(ixp+Order+1))-1
         enddo
      enddo

      do ix = 1,xDim
         do ixp = max(1,ix-Order),min(xDim,ix+Order)
            xS(ixp,ix) = 0.0d0
            do kx = kxMin(ixp,ix),kxMax(ixp,ix)
               xTempS = 0.0d0
               do lx = 1,LegPoints
                  a = wLeg(lx)*xIntScale(kx)*u(lx,kx,ix)
                  b = a*u(lx,kx,ixp)
                  xTempS = xTempS + b
               enddo
               xS(ixp,ix) = xS(ixp,ix) + xTempS
            enddo
         enddo
      enddo

      do ix = 1,xDim
         Row=ix
         do ixp = max(1,ix-Order),min(xDim,ix+Order)
            Col = ixp
            if (Col .ge. Row) then
               NewRow = HalfBandWidth+1+Row-Col
               S(NewRow,Col) = xS(ixp,ix)
c               write(26,*) ix,ixp,S(NewRow,Col)
            endif
         enddo
      enddo

      deallocate(xIntScale,xS)
      deallocate(kxMin,kxMax)

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine CalcHamiltonian(alpha,R,mu,mgamma,DD,L,
     >     Order,xPoints,LegPoints,xLeg,wLeg,xDim,
     >     xNumPoints,u,uxx,xBounds,HalfBandWidth,H)
      implicit none
      double precision, external :: VSech
      integer Order,LegPoints,xDim,xNumPoints,xBounds(*),HalfBandWidth
      double precision alpha,R,mu,mgamma,phi23,phi13,phi12,DD,L
      double precision xPoints(*),xLeg(*),wLeg(*)
      double precision H(HalfBandWidth+1,xDim)
      double precision u(LegPoints,xNumPoints,xDim),uxx(LegPoints,xNumPoints,xDim)

      integer ix,ixp,kx,lx
      integer i1,i1p
      integer Row,NewRow,Col
      integer, allocatable :: kxMin(:,:),kxMax(:,:)
      double precision a,b,m,Pi
      double precision Rall,r12,r12a,r23,r23a,r23b,r23c,r13,r13a,r13b,r13c,r14,r24,r34
      double precision u1,sys_ss_pot,V12,V23,V31
      double precision VInt,VTempInt,potvalue, xTempV
c     double precision TempPot,VInt,VTempInt
      double precision x,ax,bx,xScaledZero,xTempT,xTempS,xInt
      double precision, allocatable :: Pot(:,:)
      double precision, allocatable :: xIntScale(:),xT(:,:),xV(:,:)
      double precision, allocatable :: cosx0(:,:),cosxp(:,:),cosxm(:,:),cosx(:,:),sinx(:,:)
      double precision, allocatable :: sin12(:,:),sin23(:,:),sin13(:,:)

      double precision mu12,r0diatom,dDiatom



      allocate(xIntScale(xNumPoints),xT(xDim,xDim),xV(xDim,xDim))
      allocate(cosx0(LegPoints,xNumPoints),cosxp(LegPoints,xNumPoints),cosxm(LegPoints,xNumPoints))
      allocate(cosx(LegPoints,xNumPoints))
      allocate(sin12(LegPoints,xNumPoints))
      allocate(sin23(LegPoints,xNumPoints))
      allocate(sin13(LegPoints,xNumPoints))
      allocate(kxMin(xDim,xDim),kxMax(xDim,xDim))
      allocate(Pot(LegPoints,xNumPoints))

      Pi = 3.1415926535897932385d0

      m = -1.0d0/(2.0d0*mu*R*R)
      phi23 = datan(mgamma)
      phi12 = 0.5d0*Pi
      phi13 = -phi23

      do kx = 1,xNumPoints-1
         ax = xPoints(kx)
         bx = xPoints(kx+1)
         xIntScale(kx) = 0.5d0*(bx-ax)
         xScaledZero = 0.5d0*(bx+ax)
         do lx = 1,LegPoints
            x = xIntScale(kx)*xLeg(lx)+xScaledZero
c            cosx(lx,kx) = dcos(x)
c            cosxp(lx,kx) = dcos(x+Pi/3.0d0)
c            cosxm(lx,kx) = dcos(x-Pi/3.0d0)
            sin12(lx,kx) = dsin(x-phi12)
            sin13(lx,kx) = dsin(x-phi13)
            sin23(lx,kx) = dsin(x-phi23)

         enddo
      enddo

      do ix = 1,xDim
         do ixp = 1,xDim
            kxMin(ixp,ix) = max(xBounds(ix),xBounds(ixp))
            kxMax(ixp,ix) = min(xBounds(ix+Order+1),xBounds(ixp+Order+1))-1
         enddo
      enddo

      do kx = 1,xNumPoints-1
         do lx = 1,LegPoints
c
c            r13 = dsqrt(2.d0/dsqrt(3.0d0))*R*dabs(cosxm(lx,kx))
c            r23 = dsqrt(2.d0/dsqrt(3.0d0))*R*dabs(cosxp(lx,kx))
            r12 = dsqrt(2.d0*mgamma)*R*dabs(sin12(lx,kx))
            r13 = dsqrt( (1.0d0 + mgamma**2.d0)/(2.d0*mgamma) )*R*dabs(sin13(lx,kx))
            r23 = dsqrt( (1.0d0 + mgamma**2.d0)/(2.d0*mgamma) )*R*dabs(sin23(lx,kx))
c            call  sumpairwisepot(r12, r13, r23, potvalue)
            potvalue = VSech(r23,DD,L) + VSech(r12,DD,L) + VSech(r13,DD,L) 
            Pot(lx,kx) = alpha*potvalue
c            write(24,*) kx, lx, Pot(lx,kx)
         enddo
      enddo
      
      do ix = 1,xDim
         do ixp = max(1,ix-Order),min(xDim,ix+Order)
            xT(ix,ixp) = 0.0d0
            xV(ix,ixp) = 0.0d0
            do kx = kxMin(ixp,ix),kxMax(ixp,ix)
               xTempT = 0.0d0
               xTempV = 0.0d0
               do lx = 1,LegPoints
                  a = wLeg(lx)*xIntScale(kx)*u(lx,kx,ix)
                  xTempT = xTempT + a*uxx(lx,kx,ixp)
                  xTempV = xTempV + a*(Pot(lx,kx))*u(lx,kx,ixp)
               enddo
               xT(ix,ixp) = xT(ix,ixp) + xTempT
               xV(ix,ixp) = xV(ix,ixp) + xTempV
            enddo
         enddo
      enddo

      H = 0.0d0      
      do ix = 1,xDim
         Row=ix
         do ixp = max(1,ix-Order),min(xDim,ix+Order)
            Col = ixp
            if (Col .ge. Row) then
               NewRow = HalfBandWidth+1+Row-Col
               H(NewRow,Col) = (m*xT(ix,ixp)+xV(ix,ixp))
c     write(25,*) ix,ixp,H(NewRow,Col)
            endif
         enddo
      enddo

      deallocate(Pot)
      deallocate(xIntScale,xT,xV)
      deallocate(cosx0,cosxp,cosxm,cosx,sin12,sin23,sin13)
      deallocate(kxMin,kxMax)


      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine CalcPMatrix(NumStates,HalfBandWidth,MatrixDim,RDelt,lPsi,mPsi,rPsi,S,P)
      implicit none
      integer NumStates,HalfBandWidth,MatrixDim
      double precision RDelt
      double precision lPsi(MatrixDim,NumStates),mPsi(MatrixDim,NumStates),rPsi(MatrixDim,NumStates)
      double precision S(HalfBandWidth+1,MatrixDim)
      double precision P(NumStates,NumStates)

      integer i,j,k
      double precision a,ddot
      double precision, allocatable :: TempPsi1(:),TempPsi2(:)

      allocate(TempPsi1(MatrixDim),TempPsi2(MatrixDim))

      a = 0.5d0/RDelt

      do j = 1,NumStates
         do k = 1,MatrixDim
            TempPsi1(k) = rPsi(k,j)-lPsi(k,j)
         enddo
         call dsbmv('U',MatrixDim,HalfBandWidth,1.0d0,S,HalfBandWidth+1,TempPsi1,1,0.0d0,TempPsi2,1)
         do i = 1,NumStates
            P(i,j) = a*ddot(MatrixDim,TempPsi2,1,mPsi(1,i),1)
         enddo
      enddo

      deallocate(TempPsi1,TempPsi2)

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine CalcQMatrix(NumStates,HalfBandWidth,MatrixDim,RDelt,lPsi,mPsi,rPsi,S,Q)
      implicit none
      integer NumStates,HalfBandWidth,MatrixDim
      double precision RDelt
      double precision lPsi(MatrixDim,NumStates),mPsi(MatrixDim,NumStates),rPsi(MatrixDim,NumStates)
      double precision S(HalfBandWidth+1,MatrixDim)
      double precision Q(NumStates,NumStates)
      
      integer i,j,k
      double precision a,ddot
      double precision, allocatable :: TempPsi1(:),TempPsi2(:)
      
      allocate(TempPsi1(MatrixDim),TempPsi2(MatrixDim))
      
      a = 1.0d0/(RDelt**2)
      
      do j = 1,NumStates
         do k = 1,MatrixDim
            TempPsi1(k) = lPsi(k,j)+rPsi(k,j)-2.0d0*mPsi(k,j)
         enddo
         call dsbmv('U',MatrixDim,HalfBandWidth,1.0d0,S,HalfBandWidth+1,TempPsi1,1,0.0d0,TempPsi2,1)
         do i = 1,NumStates
            Q(i,j) = a*ddot(MatrixDim,TempPsi2,1,mPsi(1,i),1)
         enddo
      enddo
      
      deallocate(TempPsi1,TempPsi2)
      
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
      subroutine FixPhase(NumStates,HalfBandWidth,MatrixDim,S,ncv,mPsi,rPsi)
      implicit none
      integer NumStates,HalfBandWidth,MatrixDim,ncv
      double precision S(HalfBandWidth+1,MatrixDim),Psi(MatrixDim,ncv)
      double precision mPsi(MatrixDim,ncv),rPsi(MatrixDim,ncv)

      integer i,j
      double precision Phase,ddot
      double precision, allocatable :: TempPsi(:)

      allocate(TempPsi(MatrixDim))

      do i = 1,NumStates
         call dsbmv('U',MatrixDim,HalfBandWidth,1.0d0,S,HalfBandWidth+1,rPsi(1,i),1,0.0d0,TempPsi,1)
         Phase = ddot(MatrixDim,mPsi(1,i),1,TempPsi,1)
         if (Phase .lt. 0.0d0) then
            do j = 1,MatrixDim
               rPsi(j,i) = -rPsi(j,i)
            enddo
         endif
      enddo

      deallocate(TempPsi)

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine GridMakerHHL(mu,mu12,mu123,phi23,R,r0,xNumPoints,xMin,xMax,xPoints,CalcNewBasisFunc)
      implicit none
      integer xNumPoints,CalcNewBasisFunc
      double precision mu,R,r0,xMin,xMax,xPoints(xNumPoints)
      double precision mu12,mu123,phi23
      integer i,j,k,OPGRID
      double precision Pi
      double precision r0New
      double precision xRswitch
      double precision xDelt,x0,x1,x2,x3,x4,deltax


      Pi = 3.1415926535897932385d0
      r0New=r0*2.0d0
      deltax = r0/R

      xRswitch = 20.0d0*r0New/Pi
      OPGRID=1
      write (6,*) 'xMin = ', xMin, 'xMax = ', xMax
      if((OPGRID.eq.1).and.(R.gt.xRswitch)) then
         print*, 'R>xRswitch!! using modified grid!!'
         x0 = xMin
         x1 = phi23 - deltax  
         x2 = phi23 + deltax
         x3 = xMax-deltax
         x4 = xMax
         k = 1
         xDelt = (x1-x0)/dfloat(xNumPoints/4)
         do i = 1,xNumPoints/4
            xPoints(k) = (i-1)*xDelt + x0
c     write(6,*) k, xPoints(k)
            k = k + 1
         enddo
         xDelt = (x2-x1)/dfloat(xNumPoints/4)
         do i = 1,xNumPoints/4
            xPoints(k) = (i-1)*xDelt + x1
c     write(6,*) k, xPoints(k)
            k = k + 1
         enddo
         xDelt = (x3-x2)/dfloat(xNumPoints/4)
         do i = 1,xNumPoints/4
            xPoints(k) = (i-1)*xDelt + x2
c     write(6,*) k, xPoints(k)
            k = k + 1
         enddo
         xDelt = (x4-x3)/dfloat(xNumPoints/4-1)
         do i = 1, xNumPoints/4
            xPoints(k) = (i-1)*xDelt + x3
c     write(6,*) k, xPoints(k)
            k = k + 1
         enddo
         
c     FOR SMALL R, USE A LINEAR GRID
      else
         k = 1
         xDelt = (xMax-xMin)/dfloat(xNumPoints-1)
         do i = 1,xNumPoints
            xPoints(k) = (i-1)*xDelt + x0
            k = k + 1
         enddo
      endif
      
c     Smooth Grid 
      
      write(20,*) 1, xPoints(1)
      do i = 2, xNumPoints-1
         xPoints(i)=(xPoints(i-1)+2.d0*xPoints(i)+xPoints(i+1))/4.d0
      write(20,*) i, xPoints(i)
      enddo
      write(20,*) xNumPoints, xPoints(xNumPoints)
      write(20,*) ' ' 
c      write(96,15) (xPoints(k),k=1,xNumPoints)
      
 15   format(6(1x,1pd12.5))
      


      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine GridMaker(mu,R,r0,xNumPoints,xMin,xMax,xPoints,CalcNewBasisFunc)
      implicit none
      integer xNumPoints,CalcNewBasisFunc
      double precision mu,R,r0,xMin,xMax,xPoints(xNumPoints)

      integer i,j,k,OPGRID
      double precision Pi
      double precision r0New
      double precision xRswitch
      double precision xDelt,x0,x1,x2


      Pi = 3.1415926535897932385d0

      x0 = xMin
      x1 = xMax
c     write(96,*) 'x0,x1=',x0,x1
      r0New=10.0d0*r0           !/2.0d0
c     r0New=Pi/12.0*R
      xRswitch = 12.0d0*r0New/Pi
      OPGRID=1
      
      if((OPGRID.eq.1).and.(R.gt.xRswitch)) then
         print*, 'R>xRswitch!! using modified grid!!'
         x0 = xMin
         x1 = xMax - r0New/R  
         x2 = xMax
         k = 1
         xDelt = (x1-x0)/dfloat(xNumPoints/2)
         do i = 1,xNumPoints/2
            xPoints(k) = (i-1)*xDelt + x0
c            print*, k, xPoints(k), xDelt
            k = k + 1
         enddo
         xDelt = (x2-x1)/dfloat(xNumPoints/2-1)
         do i = 1,xNumPoints/2
            xPoints(k) = (i-1)*xDelt + x1
c            print*, k, xPoints(k), xDelt
            k = k + 1
         enddo
      else
         x0 = xMin
         x1 = xMax
         k = 1
         xDelt = (x1-x0)/dfloat(xNumPoints-1)
         do i = 1,xNumPoints
            xPoints(k) = (i-1)*xDelt + x0
c            print*, k, xPoints(k), xDelt
            k = k + 1
         enddo
      endif
      
c      write(96,15) (xPoints(k),k=1,xNumPoints)
 15   format(6(1x,1pd12.5))
      


      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine CalcCoupling(NumStates,HalfBandWidth,MatrixDim,RDelt,lPsi,mPsi,rPsi,S,P,Q,dP)
      implicit none
      integer NumStates,HalfBandWidth,MatrixDim
      double precision RDelt
      double precision lPsi(MatrixDim,NumStates),mPsi(MatrixDim,NumStates),rPsi(MatrixDim,NumStates)
      double precision S(HalfBandWidth+1,MatrixDim),testorth
      double precision P(NumStates,NumStates),Q(NumStates,NumStates),dP(NumStates,NumStates)

      integer i,j,k
      double precision aP,aQ,ddot
      double precision, allocatable :: lDiffPsi(:),rDiffPsi(:),TempPsi(:),TempPsiB(:),rSumPsi(:)
      double precision, allocatable :: TempmPsi(:)

      allocate(lDiffPsi(MatrixDim),rDiffPsi(MatrixDim),TempPsi(MatrixDim),
     >     TempPsiB(MatrixDim),rSumPsi(MatrixDim))
      allocate(TempmPsi(MatrixDim))

      aP = 0.5d0/RDelt
      aQ = aP*aP

      do j = 1,NumStates
         do k = 1,MatrixDim
            rDiffPsi(k) = rPsi(k,j)-lPsi(k,j)
            rSumPsi(k)  = lPsi(k,j)+mPsi(k,j)+rPsi(k,j)
c            rSumPsi(k)  = lPsi(k,j)-2.0d0*mPsi(k,j)+rPsi(k,j)
c            rSumPsi(k)  = lPsi(k,j)+rPsi(k,j)
         enddo
         call dsbmv('U',MatrixDim,HalfBandWidth,1.0d0,S,HalfBandWidth+1,rDiffPsi,1,0.0d0,TempPsi,1)   ! Calculate the vector S*rDiffPsi
         call dsbmv('U',MatrixDim,HalfBandWidth,1.0d0,S,HalfBandWidth+1,rSumPsi,1,0.0d0,TempPsiB,1)   ! Calculate the vector S*rSumPsi
         call dsbmv('U',MatrixDim,HalfBandWidth,1.0d0,S,HalfBandWidth+1,mPsi(1,j),1,0.0d0,TempmPsi,1) ! Calculate the vector S*mPsi(1,j)

         do i = 1,NumStates

c            testorth=ddot(MatrixDim,mPsi(1,i),1,TempmPsi,1)
c            write(309,*) i,j, '   testorth=',testorth

            P(i,j) = aP*ddot(MatrixDim,mPsi(1,i),1,TempPsi,1)
            dP(i,j)= ddot(MatrixDim,mPsi(1,i),1,TempPsiB,1)

            do k = 1,MatrixDim
               lDiffPsi(k) = rPsi(k,i)-lPsi(k,i)
            enddo
            Q(i,j) = -aQ*ddot(MatrixDim,lDiffPsi,1,TempPsi,1)
         enddo
      enddo

      do j=1,NumStates
	 do i=j,NumStates
            dP(i,j)=2.d0*aQ*(dP(i,j)-dP(j,i))
            dP(j,i)=-dP(i,j)
	 enddo
      enddo

      deallocate(lDiffPsi,rDiffPsi,TempPsi,rSumPsi,TempPsiB,TempmPsi)

      return
      end
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      double precision function VSech(rij,DD,L)
      
      double precision rij,DD,L
      VSech = -DD/dcosh(rij/L)**2.d0
      end 

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      double precision function phirecon(R,beta,evec,left,right,RDim,MatrixDim,RNumPoints,RPoints,order)
      implicit none
      double precision, external :: BasisPhi
      integer MatrixDim,RDim,nch,beta,i,RNumPoints,left,right,order
      double precision R,evec(MatrixDim,MatrixDim),RPoints(RNumPoints)
      phirecon = 0.0d0
      do i = 1,RDim
      phirecon = phirecon + evec(i,beta)*BasisPhi(R,left,right,order,RDim,RPoints,
     >        RNumPoints,0,i)
      enddo
      return
      end function phirecon
