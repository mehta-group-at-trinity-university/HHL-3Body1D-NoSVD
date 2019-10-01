c     234567890
      program OneDimChannels
      implicit none
      integer LegPoints,xNumPoints
      integer NumStates,PsiFlag,Order,Left,Right
      integer RSteps,CouplingFlag,CalcNewBasisFunc
      double precision alpha,mass,Shift
      double precision RLeft,RRight,RDerivDelt
      DOUBLE PRECISION RFirst,RLast,XFirst,XLast,StepX
      double precision xMin,xMax
      double precision, allocatable :: R(:)
      double precision, allocatable :: xPoints(:)

      logical, allocatable :: Select(:)

      integer iparam(11),ncv,info
      integer i,j,k,iR
      integer LeadDim,MatrixDim,HalfBandWidth
      integer xDim
      integer, allocatable :: iwork(:)
      integer, allocatable :: xBounds(:)
      double precision Tol
      double precision TotalMemory
      double precision mu, mu12, r0diatom, dDiatom, etaOVERpi, Pi

      double precision, allocatable :: LUFac(:,:),workl(:)
      double precision, allocatable :: workd(:),Residuals(:)
      double precision, allocatable :: xLeg(:),wLeg(:)
      double precision, allocatable :: u(:,:,:),uxx(:,:,:)
      double precision, allocatable :: S(:,:),H(:,:)
      double precision, allocatable :: lPsi(:,:),mPsi(:,:),rPsi(:,:),
     >     Energies(:,:)
      double precision, allocatable :: P(:,:),Q(:,:),dP(:,:)

      common/MassInfo/mu12,r0diatom,dDiatom
      
      character*64 LegendreFile
      
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
      read(5,*) Shift,Order,Left,Right
      print*, 'Shift, Order, Left, Right'
      print*, Shift,Order,Left,Right

c     read in potential parameters
      read(5,*)
      read(5,*)
      read(5,*) alpha,mass
      write(6,*) alpha,mass

c     mu12=m/2.d0
c     mu12=40.d0*m/2.d0
      Pi=dacos(-1.d0)
      write(6,*) 'Pi=',Pi

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

      RLeft = 0.0d0


c     mu=0.5
      mu=mass/dsqrt(3.0d0)
c      mu=2.d0*mu12
c      mu=mass
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

      do iR = 1,RSteps
c     must move this block inside the loop if the grid is adaptive
         print*, 'calling GridMaker'
         call GridMaker(mu,R(iR),2.0d0, xNumPoints,xMin,xMax,xPoints,CalcNewBasisFunc)
         if(CalcNewBasisFunc.eq.1) then
            print*, 'done... Calculating Basis functions'
            call CalcBasisFuncs(Left,Right,Order,xPoints,LegPoints,xLeg,
     >           xDim,xBounds,xNumPoints,0,u)
            call CalcBasisFuncs(Left,Right,Order,xPoints,LegPoints,xLeg,
     >           xDim,xBounds,xNumPoints,2,uxx)
         endif
         print*, 'done... Calculating overlap matrix'
c     must move this block inside the loop if the grid is adaptive
         
         call CalcOverlap(Order,xPoints,LegPoints,xLeg,wLeg,xDim,
     >        xNumPoints,u,xBounds,HalfBandWidth,S)

         if (CouplingFlag .ne. 0) then

            RLeft = R(iR)-RDerivDelt
            call CalcHamiltonian(alpha,RLeft,mu,Order,xPoints,
     >           LegPoints,xLeg,wLeg,xDim,xNumPoints,u,
     >           uxx,xBounds,HalfBandWidth,H)
            call MyDsband(Select,Energies,lPsi,MatrixDim,Shift,MatrixDim,
     >           H,S,HalfBandWidth+1,LUFac,LeadDim,HalfBandWidth,
     >           NumStates,Tol,Residuals,ncv,lPsi,MatrixDim,iparam,workd,workl,
     >           ncv*ncv+8*ncv,iwork,info)

            if (iR .gt. 1) call FixPhase(NumStates,HalfBandWidth,
     >           MatrixDim,S,ncv,mPsi,lPsi)
            call CalcEigenErrors(info,iparam,MatrixDim,H,
     >           HalfBandWidth+1,
     >           S,HalfBandWidth,NumStates,lPsi,Energies,ncv)
c            IF(R(iR).GT. 20.d0) Shift = 1.05d0*Energies(1,1)
            IF(iR.GT.1) Shift = 1.05d0*Energies(1,1)
            do i = 1,min(NumStates,iparam(5))
c     Energies(i,1) = Energies(i,1) + 1.875d0/(mu*RLeft*RLeft)
            enddo
c            write(100,10) RLeft,(Energies(i,1), i = 1,NumStates)
            write(6,*)
            write(6,*) RLeft
            do i = 1,min(NumStates,iparam(5))
               write(6,*) i,Energies(i,1),Energies(i,2)
            enddo

            RRight = R(iR)+RDerivDelt
            call CalcHamiltonian(alpha,RRight,mu,Order,xPoints,
     >           LegPoints,xLeg,wLeg,xDim,
     >           xNumPoints,u,uxx,xBounds,
     >           HalfBandWidth,H)
            call MyDsband(Select,Energies,rPsi,MatrixDim,Shift,
     >           MatrixDim,
     >           H,S,HalfBandWidth+1,LUFac,LeadDim,HalfBandWidth,
     >           NumStates,Tol,
     >           Residuals,ncv,rPsi,MatrixDim,iparam,workd,workl,
     >           ncv*ncv+8*ncv,iwork,info)
            call FixPhase(NumStates,HalfBandWidth,MatrixDim,S,ncv,
     >           lPsi,rPsi)
            call CalcEigenErrors(info,iparam,MatrixDim,H,
     >           HalfBandWidth+1,
     >           S,HalfBandWidth,NumStates,rPsi,Energies,ncv)
            IF(R(iR).GT. 2.2d0) Shift = 0.93d0*Energies(1,1)
c            do i = 1,min(NumStates,iparam(5))
c            enddo
c     write(100,10) RRight,(Energies(i,1), i = 1,NumStates)
            write(6,*)
            write(6,*) RRight
            do i = 1,min(NumStates,iparam(5))
               write(6,*) i,Energies(i,1),Energies(i,2)
            enddo

         endif
         
         call CalcHamiltonian(alpha,R(iR),mu,Order,xPoints,
     >        LegPoints,xLeg,wLeg,xDim,
     >        xNumPoints,u,uxx,xBounds,
     >        HalfBandWidth,H)
         call MyDsband(Select,Energies,mPsi,MatrixDim,Shift,MatrixDim,
     >        H,S,HalfBandWidth+1,LUFac,LeadDim,HalfBandWidth,NumStates,
     >        Tol,Residuals,ncv,mPsi,MatrixDim,iparam,workd,workl,
     >        ncv*ncv+8*ncv,iwork,info)
         if (CouplingFlag .ne. 0) call FixPhase(NumStates,HalfBandWidth,
     >        MatrixDim,S,ncv,rPsi,mPsi)
         call CalcEigenErrors(info,iparam,MatrixDim,H,HalfBandWidth+1,S,
     >        HalfBandWidth,NumStates,mPsi,Energies,ncv)
         IF(R(iR).GT. 2.2d0) Shift = 0.93d0*Energies(1,1)
c         do i = 1,min(NumStates,iparam(5))
c     Energies(i,1) = Energies(i,1) + 1.875d0/(mu*R(iR)*R(iR))
c         enddo
         write(200,20) R(iR),(Energies(i,1), i = 1,min(NumStates,
     >        iparam(5)))

         write(6,*)
         write(6,*) R(iR)
         do i = 1,min(NumStates,iparam(5))
            write(6,*) i,Energies(i,1),Energies(i,2)
         enddo

         if (CouplingFlag .ne. 0) then
c     call CalcPMatrix(min(NumStates,iparam(5)),HalfBandWidth,MatrixDim,RDerivDelt,lPsi,mPsi,rPsi,S,P)
c     call CalcQMatrix(min(NumStates,iparam(5)),HalfBandWidth,MatrixDim,RDerivDelt,lPsi,mPsi,rPsi,S,Q)
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
      subroutine CalcHamiltonian(alpha,R,mu,
     >     Order,xPoints,LegPoints,xLeg,wLeg,xDim,
     >     xNumPoints,u,uxx,xBounds,HalfBandWidth,H)
      implicit none
      integer Order,LegPoints,xDim,xNumPoints,xBounds(*),HalfBandWidth
      double precision alpha,R,mu
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

      double precision mu12,r0diatom,dDiatom
c      common/MassInfo/mu12,r0diatom,dDiatom

      allocate(xIntScale(xNumPoints),xT(xDim,xDim),xV(xDim,xDim))
      allocate(cosx0(LegPoints,xNumPoints),cosxp(LegPoints,xNumPoints),cosxm(LegPoints,xNumPoints))
      allocate(cosx(LegPoints,xNumPoints))
      allocate(kxMin(xDim,xDim),kxMax(xDim,xDim))
      allocate(Pot(LegPoints,xNumPoints))

      Pi = 3.1415926535897932385d0

      m = -1.0d0/(2.0d0*mu*R*R)

      do kx = 1,xNumPoints-1
         ax = xPoints(kx)
         bx = xPoints(kx+1)
         xIntScale(kx) = 0.5d0*(bx-ax)
         xScaledZero = 0.5d0*(bx+ax)
         do lx = 1,LegPoints
            x = xIntScale(kx)*xLeg(lx)+xScaledZero
            cosx(lx,kx) = dcos(x)
            cosxp(lx,kx) = dcos(x+Pi/3.0d0)
            cosxm(lx,kx) = dcos(x-Pi/3.0d0)
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
            r12 = dsqrt(2.d0/dsqrt(3.0d0))*R*dabs(cosx(lx,kx))
            r13 = dsqrt(2.d0/dsqrt(3.0d0))*R*dabs(cosxm(lx,kx))
            r23 = dsqrt(2.d0/dsqrt(3.0d0))*R*dabs(cosxp(lx,kx))
            call  sumpairwisepot(r12, r13, r23, potvalue)
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
      deallocate(cosx0,cosxp,cosxm,cosx)
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
      subroutine GridMaker222(m,mu,R,r0,xNumPoints,xMin,xMax,yNumPoints,yMin,yMax,xPoints,yPoints)
      implicit none
      integer xNumPoints,yNumPoints
      double precision m,mu,R,r0,xMin,xMax,yMin,yMax,xPoints(xNumPoints),yPoints(yNumPoints)

      integer i,j,k
      double precision Pi
      double precision r0New
      double precision xRswitch,yRswitch
      double precision xDelt,x0,x1,x2
      double precision yDelt,y0,y1,y2

      Pi = 3.1415926535897932385d0

      r0New = 3.0d0*r0

      xRswitch = dsqrt(dsqrt(3.0d0)*r0New**2/(1.0d0+dcos(2.0d0*Pi*(1/12.0d0+1/3.0d0))))

      if (R .gt. xRswitch) then
         x0 = xMin
         x1 = 0.5d0*dacos(dsqrt(3.0d0)*r0New**2/R**2-1.0d0) - Pi/3.0d0
         x2 = xMax
         k = 1
         xDelt = (x1-x0)/dfloat(xNumPoints/2)
         do i = 1,xNumPoints/2
            xPoints(k) = (i-1)*xDelt + x0
            k = k + 1
         enddo
         xDelt = (x2-x1)/dfloat(xNumPoints/2-1)
         do i = 1,xNumPoints/2
            xPoints(k) = (i-1)*xDelt + x1
            k = k + 1
         enddo
      else
         x0 = xMin
         x1 = xMax
         k = 1
         xDelt = (x1-x0)/dfloat(xNumPoints-1)
         do i = 1,xNumPoints
            xPoints(k) = (i-1)*xDelt + x0
            k = k + 1
         enddo
      endif

      yRswitch = dsqrt(dsqrt(3.0d0)*r0New**2/(1.0d0-dcos(Pi/4.0d0)))

      if (R .gt. yRswitch) then
         y0 = yMin
         y1 = 0.5d0*dacos(1.0d0-dsqrt(3.0d0)*r0New**2/R**2)
         y2 = yMax
         k = 1
         yDelt = (y1-y0)/dfloat(yNumPoints/2)
         do i = 1,yNumPoints/2
            yPoints(k) = (i-1)*yDelt + y0
            k = k + 1
         enddo
         yDelt = (y2-y1)/dfloat(yNumPoints/2-1)
         do i = 1,yNumPoints/2
            yPoints(k) = (i-1)*yDelt + y1
            k = k + 1
         enddo
      else
         y0 = yMin
         y1 = yMax
         k = 1
         yDelt = (y1-y0)/dfloat(yNumPoints-1)
         do i = 1,yNumPoints
            yPoints(k) = (i-1)*yDelt + y0
            k = k + 1
         enddo
      endif

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
      subroutine GridMaker111(m,mu,R,r0,xNumPoints,xMin,xMax,yNumPoints,yMin,yMax,xPoints,yPoints)
      implicit none
      integer xNumPoints,yNumPoints
      double precision m,mu,R,r0,xMin,xMax,yMin,yMax,xPoints(xNumPoints),yPoints(yNumPoints)

      
      integer i,j,k
      double precision Pi
      double precision r0New
      double precision xRswitch,yRswitch
      double precision xDelt,x0,x1,x2
      double precision yDelt,y0,y1,y2
      
      Pi = 3.1415926535897932385d0
      

      r0New = 3.0d0*r0
      
      xRswitch = dsqrt(dsqrt(3.0d0)*r0New**2/(1.0d0+dcos(2.0d0*Pi*(1/12.0d0+1/3.0d0)) ))
      
      
      if (R .gt. xRswitch) then
         x0 = xMin
         x1 = 0.5d0*dacos(dsqrt(3.0d0)*r0New**2/R**2-1.0d0) - Pi/3.0d0
         x2 = xMax
         k = 1
         xDelt = (x1-x0)/dfloat(xNumPoints/2)
         do i = 1,xNumPoints/2
            xPoints(k) = (i-1)*xDelt + x0
            k = k + 1
         enddo
         xDelt = (x2-x1)/dfloat(xNumPoints/2-1)
         do i = 1,xNumPoints/2
            xPoints(k) = (i-1)*xDelt + x1
            k = k + 1
         enddo
      else
         x0 = xMin
         x1 = xMax
         k = 1
         xDelt = (x1-x0)/dfloat(xNumPoints-1)
         do i = 1,xNumPoints
            xPoints(k) = (i-1)*xDelt + x0
            k = k + 1
         enddo
      endif
      
      yRswitch = dsqrt(dsqrt(3.0d0)*r0New**2/(1.0d0-dcos(Pi/4.0d0)))
      
      if (R .gt. yRswitch) then
         y0 = yMin
         y1 = 0.5d0*dacos(1.0d0-dsqrt(3.0d0)*r0New**2/R**2)
         y2 = yMax
         k = 1
         yDelt = (y1-y0)/dfloat(yNumPoints/2)
         do i = 1,yNumPoints/2
            yPoints(k) = (i-1)*yDelt + y0
            k = k + 1
         enddo
         yDelt = (y2-y1)/dfloat(yNumPoints/2-1)
         do i = 1,yNumPoints/2
            yPoints(k) = (i-1)*yDelt + y1
            k = k + 1
         enddo
      else
         y0 = yMin
         y1 = yMax
         k = 1
         yDelt = (y1-y0)/dfloat(yNumPoints-1)
         do i = 1,yNumPoints
            yPoints(k) = (i-1)*yDelt + y0
            k = k + 1
         enddo
      endif
      
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
         call dsbmv('U',MatrixDim,HalfBandWidth,1.0d0,S,HalfBandWidth+1,rDiffPsi,1,0.0d0,TempPsi,1)
         call dsbmv('U',MatrixDim,HalfBandWidth,1.0d0,S,HalfBandWidth+1,rSumPsi,1,0.0d0,TempPsiB,1)
c         call dsbmv('U',MatrixDim,HalfBandWidth,1.0d0,S,HalfBandWidth+1,mPsi(1,j),1,0.0d0,TempmPsi,1)

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



