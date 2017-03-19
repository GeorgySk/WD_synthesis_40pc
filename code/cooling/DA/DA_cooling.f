C***********************************************************************
C     TODO: give better names to variables and the subr.itself  
C     NOTE: example of terrible coding
      subroutine incoolda(param,irdr,ntrk,ncol,mtrk,ttrk,tprewdda,
     &           luminosity,tefftrk,gtrk)
C=======================================================================
C
C     This subroutine reads the cooling tables for the sequences of
C     WD DA by Althaus et al. (2009) and Renedo et al. (2010) and
C     interpolates/extrapolates according to the vector of time 
C     reference
C     This subroutine can be used to read the tables for one 
C     metallicity in particular, according to the input index.
C
C     Modification at 07/2012 (ER Cojocaru)
C
C-----------------------------------------------------------------------
C
C     Input parameters:
C
C       param: parameter that decides type of reading according to the 
C              file
C       irdr: initial index for the series of files that it's going to 
C             read
C       ntrk: number of files according to vector of time reference
C
C-----------------------------------------------------------------------
C
C     Output parameters
C
C       ncol: number of masses for which the cooling sequence has been 
C             calculated
C       mtrk: mass of the WD. [M0]
C       luminosity [log(L/L_0)]
C       tefftrk: effective temeprature
C       gtrk: g
C
C=======================================================================
      implicit double precision (a-h,m,o-z)

C     ---   Parameters   ---
      integer nrow
      parameter (nrow=650)

C     ---   Declaration of the variables   ---

      integer j,k,param,ncol,irdr,ntrk(ncol)
      double precision a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14
      double precision a15,a16,a17,a18,a19,a20,a21,a22,a23,a24,a25,a26
      double precision mtrk(ncol),tprewdda(ncol)
      double precision ttrk(ncol,nrow),luminosity(ncol,nrow)
      double precision tefftrk(ncol,nrow),gtrk(ncol,nrow)

C     reading the masses, depending on the format of the input files,
C     because they are different

      if(param.eq.1) then
C       --- Masses  (Z=0.001) ---
        mtrk(1)=0.505
        mtrk(2)=0.553
        mtrk(3)=0.593
        mtrk(4)=0.627
        mtrk(5)=0.660
        mtrk(6)=0.692
        mtrk(7)=0.863
C       --- Tprew_WD LPCODE ---
        tprewdda(1)=0.0
        tprewdda(2)=0.0
        tprewdda(3)=0.0
        tprewdda(4)=0.0
        tprewdda(5)=0.0
        tprewdda(6)=0.0
        tprewdda(7)=0.0
      end if

      if(param.eq.2) then
C       --- Masses (Z=0.01) ---
        mtrk(1)=0.524
        mtrk(2)=0.570
        mtrk(3)=0.593
        mtrk(4)=0.609
        mtrk(5)=0.632
        mtrk(6)=0.659
        mtrk(7)=0.705
        mtrk(8)=0.767
        mtrk(9)=0.837
        mtrk(10)=0.877
C       --- Tprew_WD LPCODE ---
        tprewdda(1)=0.0
        tprewdda(2)=0.0
        tprewdda(3)=0.0
        tprewdda(4)=0.0
        tprewdda(5)=0.0
        tprewdda(6)=0.0
        tprewdda(7)=0.0
        tprewdda(8)=0.0
        tprewdda(9)=0.0
        tprewdda(10)=0.0
      end if

      if(param.eq.3) then
C       --- Masses  (Z=0.03/0.06) ---
        mtrk(1)=0.524
        mtrk(2)=0.570
        mtrk(3)=0.593
        mtrk(4)=0.609
        mtrk(5)=0.632
        mtrk(6)=0.659
        mtrk(7)=0.705
        mtrk(8)=1.000
C       --- Tprew_WD LPCODE ---
        tprewdda(1)=11.117
        tprewdda(2)=2.7004
        tprewdda(3)=1.699
        tprewdda(4)=1.2114
        tprewdda(5)=0.9892
        tprewdda(6)=0.7422
        tprewdda(7)=0.4431
        tprewdda(8)=0.0
      end if

C     ---   Choosing the format of the input file   ---
      if(param.eq.1.OR.param.eq.2) then
C       ---   Reading the tables   ---
        do k=1,ncol
C         ---  Reading the cooling curves  ---
          do j=1,nrow
            read(irdr+k,*,end=5) a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,
     &                           a13
            ttrk(k,j)=a5/1000.0
            tefftrk(k,j)=10**a2
            gtrk(k,j)=a12
            luminosity(k,j)=a1
          end do
5         ntrk(k)=j-1
        end do
      else
C       ---   Reading the tables   ---
        do k=1,ncol
C         ---   Reading the cooling curves   ---
          do j=1,nrow
            read(irdr+k,*,end=10) a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,
     &                            a12,a13,a14,a15,a16,a17,a18,a19,a20,
     &                            a21,a22,a23,a24,a25,a26
            ttrk(k,j)=10.0**a9/1000.0
C           --- Restricting the lifetime of pre-WD ---
            ttrk(k,j)=ttrk(k,j)-tprewdda(k)
            tefftrk(k,j)=10**a2
            gtrk(k,j)=a23
            luminosity(k,j)=a1   
          end do
10        ntrk(k)=j-1
        end do
      end if

      return
      end
C***********************************************************************