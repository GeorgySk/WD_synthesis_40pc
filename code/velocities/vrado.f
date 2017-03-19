C***********************************************************************
C     TODO:rewrite      
      subroutine vrado
C***********************************************************************
C     This subroutine calculates the heliocentric velocities, starting
C     from the proper motions in galactic coordinates, making zero
C     the component of radial velocity
C***********************************************************************
      implicit double precision (a-h,m,o-z)
       
      integer numberOfStars,i,ntwd
      double precision k,xcb,xsb,xcl,xsl,r
      double precision a1,a2,b1,b2,b3,c1,c2,c3
C     ---   Parameters  ---
      parameter (numberOfStars=6000000)
      parameter (k=4.74d0)
C     ---   Dimensiones   --- 
      double precision mpl(numberOfStars),mpb(numberOfStars),
     &                 vr(numberOfStars)
      double precision rgac(numberOfStars),lgac(numberOfStars),
     &                 bgac(numberOfStars)
      double precision iwd(numberOfStars)

C    ---   Commons  ---
      common /lb/ lgac,bgac
      common /paral/ rgac
      common /mopro/ mpb,mpl,vr
      common /index/ iwd,ntwd

C     ---  Calculating the heliocentric velocity  
C          making zero the radial velocity ---
      do 1 i=1,ntwd
        xcb=cos(bgac(i))
        xsb=sin(bgac(i))
        xcl=cos(lgac(i))
        xsl=sin(lgac(i)) 
        r=rgac(i)*1000.0
        a1=-k*xcb*xsl
        b1=-k*xsb*xcl
        c1=0.0
        a2=k*xcb*xcl
        b2=-k*xsb*xsl
        c2=0.0
        b3=k*xcb
        c3=0.0
1     continue
         
      return
      end
C***********************************************************************
