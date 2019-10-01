      program plot
      implicit none
      double precision, allocatable :: R(:), coupling(:,:,:)
      double complex, allocatable :: SS(:,:,:)
      integer NumChan, NumPoints, file
      integer i,j,k,l
      NumChan=4
      NumPoints=400
      file=103

      allocate(SS(NumChan,NumChan,NumPoints))
      allocate(R(NumPoints))
      allocate(coupling(NumChan,NumChan,NumPoints))

      open(file)
      do i=1,NumPoints
         read(file,*) R(i)
c         read(file,*)
         R(i)=R(i)
         do k=1,NumChan
            read(file,*) (coupling(k,j,i), j = 1,NumChan)
            do l=1, NumChan
c               if (coupling(k,l,i).le.1e-16) coupling(k,l,i)=0.0d0
            enddo
         enddo
c         read(file,*)
      enddo
      close(file)


c      do i=1,NumPoints
c         do j=1,NumChan
c            coupling(1,j,i) = dsqrt(real(SS(1,j,i)*conjg(SS(1,j,i))))
c            coupling(2,2,i) = dsqrt(real(SS(2,2,i)*conjg(SS(2,2,i))))
c         enddo
c      enddo




      do i=1,NumPoints
c         write(6,20) R(i), (SS(1,j,i), j = 1,NumChan), SS(2,2,i)!, SS(3,3), SS(4,4), SS(5,5) ! 
         write(6,20), R(i), (coupling(1,j,i), j = 1,2)!, coupling(2,2,i) !, SS(3,3), SS(4,4), SS(5,5) ! 
c         write(6,20), R(i), (coupling(j,j,i), j = 1,2)!, coupling(2,2,i) !, SS(3,3), SS(4,4), SS(5,5) ! 
c         write(6,20) R(i), (SS(1,j)*conjg(SS(1,j)), j = 1,NumChan)!, sqrt(SS(2,2)*conjg(S(2,2))), ! 
!     >        sqrt(SS(3,3)*conjg(S(3,3))), sqrt(SS(4,4)*conjg(S(4,4))),sqrt(SS(5,5)*conjg(S(5,5)))! 
      enddo

      deallocate(R)
      deallocate(SS,coupling)
      
 20   format(1P,100e16.8)
 21   format(1P,100e16.7)
      end
