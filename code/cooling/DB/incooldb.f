C***********************************************************************
C     TODO: rewrite      
      subroutine incooldb(param,irdr,ncol,ntrk,mtrk,ttrk,tprewd,
     &           luminosity,tefftrk,gtrk)
C=======================================================================
C
C     This subroutine reads the cooling tables by:
C     QUESTION: personal what?
C     * Althaus for Z=0.001 (personal mail/post 01/2013)
C     * Althaus PG-DB, Z=0.01
C     * Althaus for Z=0.06
C
C     Created by S. Torres
C     Modifitions 04/2013 (ER Cojocaru)
C-----------------------------------------------------------------------
C     Input parameters:
C       param: parameter deciding what group of sequences to read
C       irdr: first file of the respected group
C       ncol: number of sequences in group (number of columns)
C-----------------------------------------------------------------------
C     Output parameters
C       ntrk: vector with number of points of every sequence of the 
C             group 
C       mtrk: vector of masses
C       ttrk: matrix with the vectors tcool (Gyr) as columns
C       tprewddb: vector with previous times that must be substracted
C       luminosity: matrix with vectors log(L/Lo) as columns
C       tefftrk: matrix with vectors Teff as columns
C       gtrk: matrix with vectors log(g) as columns
C=======================================================================
      implicit double precision (a-h,m,o-z)

C     ---   Parameters   ---
      integer nrow
      parameter (nrow=400)

C     ---   Declaration of variables   ---
      integer i,j,ncol,irdr,param,ntrk(ncol)
      double precision a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,
     &                 a15,a16,a17,a18,a19,a20,a21,a22,a23,a24,a25,a26
      double precision mtrk(ncol),tprewd(ncol)
      double precision ttrk(ncol,nrow),luminosity(ncol,nrow)
      double precision tefftrk(ncol,nrow),gtrk(ncol,nrow)

      if(param.eq.1) then
C       --- Mass (Z=0.001) new by Leandro ---
        mtrk(1)=0.5047
        mtrk(2)=0.5527
        mtrk(3)=0.59328
        mtrk(4)=0.62738
        mtrk(5)=0.6602
        mtrk(6)=0.69289
        mtrk(7)=0.8637
      else if(param.eq.2) then
C       --- Mass (Z=0.01) PG-DB ---
        mtrk(1)=0.514
        mtrk(2)=0.53
        mtrk(3)=0.542
        mtrk(4)=0.565
        mtrk(5)=0.584
        mtrk(6)=0.61
        mtrk(7)=0.664
        mtrk(8)=0.741
        mtrk(9)=0.869
      else
C       --- Mass (Z=0.06) LPCODE ---
        mtrk(1)=0.524
        mtrk(2)=0.570
        mtrk(3)=0.593
        mtrk(4)=0.61
        mtrk(5)=0.632
        mtrk(6)=0.659
        mtrk(7)=0.70
        mtrk(8)=0.76
        mtrk(9)=0.87
      end if

      if(param.eq.1) then
C       --- Tprew_WD (Z=0.001) new ---
        tprewd(1)=0.0
        tprewd(2)=0.0
        tprewd(3)=0.0
        tprewd(4)=0.0
        tprewd(5)=0.0
        tprewd(6)=0.0
        tprewd(7)=0.0     
      else
C       --- Tprew_WD (Z=0.06) LPCODE ---
        tprewd(1)=11.117
        tprewd(2)=2.7004
        tprewd(3)=1.699
        tprewd(4)=1.2114
        tprewd(5)=0.9892
        tprewd(6)=0.7422
        tprewd(7)=0.4431
        tprewd(8)=0.287
        tprewd(9)=0.114
      end if
 
      if(param.eq.1.OR.param.eq.3) then
C       ---   Reading the tables   ---
        do i=1,ncol
C         ---   Reading unit   ---
          irdr=irdr+1
C         ---   Reading the cooling curves   --
          do j=1,nrow
            read(irdr,*,end=2) a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,
     &                         a13,a14,a15,a16,a17,a18,a19,a20,a21,a22,
     &                         a23,a24,a25,a26
            ttrk(i,j)=10.0**a9/1000.0
C           --- Restricting the lifetime of pre-WD's ---
            ttrk(i,j)=ttrk(i,j)-tprewd(i)
            tefftrk(i,j)=10.0**a2
            gtrk(i,j)=a23
            luminosity(i,j)=a1      
          end do
2         ntrk(i)=j-1       
        end do
      else
C       ---  Reading the tables   ---
        do i=1,ncol
C         ---   Reading unit   ---
          irdr=irdr+1 
C         ---   Reading the cooling curves   --
          do j=1,nrow
            read (irdr,*,end=20) a1,a2,a3,a4      
            ttrk(i,j)=a3
            tefftrk(i,j)=10.0**a2
            gtrk(i,j)=a4
            luminosity(i,j)=a1
          end do
20        ntrk(i)=j-1
        end do
      end if
 
      return
      end
C***********************************************************************