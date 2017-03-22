      subroutine color(numberOfMassesWithColors,ntrk,massOfWD,
     &           luminosity,color_U,color_B,color_R,color_V,color_I)
C=======================================================================
C     This subroutine reads the colors by Rene and interpolates 
C     according to the vector of reference time
C-----------------------------------------------------------------------
C     Input parameters
C       ntrk: number of elements vector of reference time
C
C-----------------------------------------------------------------------
C     Output parameters:
C       numberOfMassesWithColors: number of masses for which colors are 
C                                 calculated
C       massOfWD: mass of the WD. [M0]
C       luminosity: luminosity [log(L/L_0)]
C       color_U: M_U
C       color_B: M_B
C       color_V: M_V
C       color_R: M_R
C       color_I: M_I
C
C=======================================================================
      implicit double precision (a-h,m,o-z)

C     ---   Parameters   ---
      integer nrow
      parameter (nrow=650)
      integer irdr
      parameter (irdr=60)
      double precision vtanmin
      parameter (vtanmin=30)

C     ---   Declaration of variables   ---
      integer i,k,numberOfMassesWithColors,
     &        ntrk(numberOfMassesWithColors)
      double precision a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14
      double precision a15,a16,a17,a18,a19,a20,a21,a22,a23,a24,a25,a26
      double precision a27,a28,massOfWD(numberOfMassesWithColors)
      double precision tprewdda(numberOfMassesWithColors),
     &                 ttrk(numberOfMassesWithColors,nrow)
      double precision luminosity(numberOfMassesWithColors,nrow),
     &                 color_U(numberOfMassesWithColors,nrow),
     &                 color_B(numberOfMassesWithColors,nrow)
      double precision color_V(numberOfMassesWithColors,nrow),
     &                 color_R(numberOfMassesWithColors,nrow),
     &                 color_I(numberOfMassesWithColors,nrow)
      double precision Hg(nrow),ggii(nrow),Hbj(nrow),bjr(nrow)
      double precision g,gr,bj,term

C     ---  Commons  ---
C     Curve to delimit RPMD, for Mwd=0.61Mo
      common /Hgcurve/ Hg,ggii
      common /Hbjcurve/ Hbj,bjr

C     read masses
      massOfWD(1)=0.524
      massOfWD(2)=0.570
      massOfWD(3)=0.593
      massOfWD(4)=0.610
      massOfWD(5)=0.632
      massOfWD(6)=0.659
      massOfWD(7)=0.705
      massOfWD(8)=0.767
      massOfWD(9)=0.837
      massOfWD(10)=0.877
      
C     --- Tprew_WD LPCODE ---
      tprewdda(1)=11.117
      tprewdda(2)=2.7004
      tprewdda(3)=1.699
      tprewdda(4)=1.2114
      tprewdda(5)=0.9892
      tprewdda(6)=0.7422
      tprewdda(7)=0.4431
      tprewdda(8)=0.2686
      tprewdda(9)=0.200
      tprewdda(10)=0.114

C     read values from files
      do k=1,numberOfMassesWithColors
        do i=1,nrow
          read(irdr+k,*,end=15) a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,
     &                          a13,a14,a15,a16,a17,a18,a19,a20,a21,a22,
     &                          a23,a24,a25,a26,a27,a28
          luminosity(k,i)=a3
          ttrk(k,i)=a4-tprewdda(k)
          color_U(k,i)=a24
          color_B(k,i)=a25
          color_V(k,i)=a26
          color_R(k,i)=a27
          color_I(k,i)=a28
        enddo
15      ntrk(k)=i-1
      end do
      
      term=5*log10(vtanmin)-3.379
      do i=1,ntrk(4)
C       SDSS
        g=color_V(4,i)+0.630*(color_B(4,i)-color_V(4,i))-0.124
C       Hg=Mg+5*log10(v_tan)-3.379
        Hg(i)=g+term
        ggii(i)=1.646*(color_V(4,i)-color_R(4,i))+1.007*(color_R(4,i)-
     &          color_I(4,i))-0.375
C       SuperCosmos
        gr=1.646*(color_V(4,i)-color_R(4,i))-0.139
        bj=0.15+0.13*gr+g
C       Hbj=Mbj+5*log10(v_tan)-3.379
        Hbj(i)=bj+term
        bjr(i)=0.28+1.13*gr
      end do
      
      return
      end
C***********************************************************************










      subroutine color2(table)
C=======================================================================
C     This subroutine reads the colors by Rene and interpolates 
C     according to the vector of reference time
C-----------------------------------------------------------------------
C     I/O parameters:
C       table: instance of FileGroupInfo
C=======================================================================
      use external_types
      implicit double precision (a-h,m,o-z)

      double precision vtanmin
      parameter (vtanmin=30)

      TYPE(FileGroupInfo) :: table
      integer i,k
      double precision a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14
      double precision a15,a16,a17,a18,a19,a20,a21,a22,a23,a24,a25,a26
      double precision a27,a28
C     TODO: size should be table%nrow      
      double precision Hg(650),ggii(650),Hbj(650),
     &                 bjr(650)
      double precision g,gr,bj,term

C     Curve to delimit RPMD, for Mwd=0.61Mo
      common /Hgcurve/ Hg,ggii
      common /Hbjcurve/ Hbj,bjr

C     Read masses
      table%mass(1)=0.524
      table%mass(2)=0.570
      table%mass(3)=0.593
      table%mass(4)=0.610
      table%mass(5)=0.632
      table%mass(6)=0.659
      table%mass(7)=0.705
      table%mass(8)=0.767
      table%mass(9)=0.837
      table%mass(10)=0.877
      
C     --- Tprew_WD LPCODE ---
      table%prevTime(1)=11.117
      table%prevTime(2)=2.7004
      table%prevTime(3)=1.699
      table%prevTime(4)=1.2114
      table%prevTime(5)=0.9892
      table%prevTime(6)=0.7422
      table%prevTime(7)=0.4431
      table%prevTime(8)=0.2686
      table%prevTime(9)=0.200
      table%prevTime(10)=0.114

C     read values from files
      do k=1,table%ncol
        do i=1,table%nrow
          read(table%initLink+k,*,end=15) a1,a2,a3,a4,a5,a6,a7,a8,a9,
     &        a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20,a21,a22,a23,
     &        a24,a25,a26,a27,a28
          table%luminosity(k,i)=a3
          table%coolingTime(k,i)=a4-table%prevTime(k)
          table%color_U(k,i)=a24
          table%color_B(k,i)=a25
          table%color_V(k,i)=a26
          table%color_R(k,i)=a27
          table%color_I(k,i)=a28
        enddo
15      table%ntrk(k)=i-1
      end do
      
      term=5*log10(vtanmin)-3.379
      do i=1,table%ntrk(4)
C       SDSS
        g = table%color_V(4,i)+0.630*(table%color_B(4,i)-
     &      table%color_V(4,i))-0.124
C       Hg=Mg+5*log10(v_tan)-3.379
        Hg(i)=g+term
        ggii(i)=1.646*(table%color_V(4,i)-table%color_R(4,i))+
     &          1.007*(table%color_R(4,i)-table%color_I(4,i))-0.375
C       SuperCosmos
        gr=1.646*(table%color_V(4,i)-table%color_R(4,i))-0.139
        bj=0.15+0.13*gr+g
C       Hbj=Mbj+5*log10(v_tan)-3.379
        Hbj(i)=bj+term
        bjr(i)=0.28+1.13*gr
      end do
      
      return
      end
C***********************************************************************