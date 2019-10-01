
      
      subroutine sumpairwisepot(r12, r13, r23, potvalue) 

      
c     returns the value of the sum of pair-wise interactions 
      
      implicit none
      double precision c0, cutoff, potvalue, r12, r13, r14, r23, r24, r34
      double precision mv, ms, Vv, Vs, L, dd,mu 
c      common /vindex/ index(67,3)
c      mu=1.0d0
c      dimension vcof(67)

c      potvalue=0.d0
c      potvalue=Vs*exp(-ms*abs(r12))+Vv*exp(-mv*abs(r12))
c      potvalue=potvalue+Vs*exp(-ms*abs(r13))+Vv*exp(-mv*abs(r13))
c      potvalue=potvalue+Vs*exp(-ms*abs(r23))+Vv*exp(-mv*abs(r23))

c      c0=-8.694192e+00 ! a2=100
c      c0=-8.605163e+00 ! a2=-100
c      c0=-7.321700e+00 ! a2=-2
c      c0=-8.887028e+00 ! a2=20
c      cutoff=1.d0
c      potvalue=c0*cutoff*(
c     >     dexp(-(r12*cutoff)**2.d0) +
c     >     dexp(-(r13*cutoff)**2.d0) + 
c     >     dexp(-(r23*cutoff)**2.d0))  

c     a2even=20, with one deep (even) E=-4.218858455372532 one shallow  E=-0.00291449 bound state L=1 
c      L=1.0d0
c      dd=6.272844447382345d0 ! 
      


c     a2=Infinity
c      L=4.0d0
c      dd=0.375d0

c      a2=Infinity
c      L=2.0d0
c      dd=1.5d0

c     a2=Infinity 
c      L=1.0d0
c      dd=6.0d0 ! a2even=Infinity with one deep bound state L=1

c      a2even=Infinity with one deep bound state L=1/2
c      L=0.5d0
c      dd=24.0d0                 ! 
      
c      a2even=Infinity with one deep bound state L=1/4
c      L=0.25d0
c      dd=96.0d0
      
c      a2even=2 no deep bound state
c      L=0.125d0
c      dd=4.2632742048833

c      a2even=2 nodeep bound state
c      L=0.25d0
c      dd=2.277328123623222
      
c      a2even=2 nodeep bound state
c      L=0.5d0
c      dd=1.3068628111568357

c      a2even=2 nodeep bound state
      L=1.0d0
      dd=0.853138729379683

c      a2even=2 nodeep bound state
c      L=2.0d0
c      dd=0.5d0

c      a2even=2 nodeep bound state
c      L=4.0d0 
c      dd=0.17569845216574295

c      a2even=2 nodeep bound state
c      L=8.0d0
c      dd=0.04939720050161369

c      a2even=2 nodeep bound state
c      L=16.0d0
c      dd=0.01295818465997589d0

c     a2=20 no deep bound state  B2=0.0025096021704961394
c      L=1.0d0
c      dd=0.05260553185042266d0
      
c     a2=10 no deep bound state B2=0.010144578990854353
c      L=1.0d0
c      dd=0.11086487977899903d0

c      a2=5 no deep bound state B2=0.0420734544977747
c      L=1.0d0
c      dd=0.24719160215364333d0

      potvalue=-dd*(dcosh(r12/L)**(-2.0d0) + 
     >     dcosh(r13/L)**(-2.0d0) + 
     >     dcosh(r23/L)**(-2.0d0)) ! 
c      potvalue=-dd*(dcosh(2.0d0*r12/L)**(-2.0d0) + dcosh(2.0d0*r13/L)**(-2.0d0) + dcosh(2.0d0*r23/L)**(-2.0d0)) ! 

      end 

      
