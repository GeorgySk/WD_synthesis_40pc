      subroutine colordb(ncolb,numberOfPointsInSequence,
     &           massSequence,luminosityDB,colorDB_U,colorDB_B,
     &           colorDB_V,colorDB_R,colorDB_I)
C=======================================================================
C     This subroutine reads the color tables by Bergeron for the  
C     sequences of WD DB
C-----------------------------------------------------------------------
C     Input parameters:
C       ncolb: number of sequences
C-----------------------------------------------------------------------
C     Output parameters:
C       massSequence: masses seq.
C       numberOfPointsInSequence
C       luminosityDB: luminosity
C       colorDB_U,colorDB_B,colorDB_V,colorDB_R,colorDB_I: colors DBs 
C=======================================================================
      implicit double precision (a-h,m,o-z)

C     ---   Declaration of variables   ---
      integer i,j,ncolb,nrowb,irdr
      double precision a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,
     &                 a15,a16,a17,a18,a19,a20,a21,a22,a23,a24,a25,a26,
     &                 a27

C     ---   Parameters   ---
      parameter (nrowb=60)

C     ---   Dimensiones   ---
      integer numberOfPointsInSequence(ncolb)
      double precision massSequence(ncolb),luminosityDB(ncolb,nrowb)
      double precision colorDB_U(ncolb,nrowb),colorDB_B(ncolb,nrowb),
     &                 colorDB_V(ncolb,nrowb),colorDB_R(ncolb,nrowb),
     &                 colorDB_I(ncolb,nrowb)
      
      massSequence(1)=0.5
      massSequence(2)=0.6
      massSequence(3)=0.7
      massSequence(4)=0.8
      massSequence(5)=0.9
      massSequence(6)=1.0
      massSequence(7)=1.2

C     ---   Initialization   ---
      irdr=131
 
C     ---   reading the tables   ---
      do 3 i=1,ncolb
C       ---  Reading unit   ---
        irdr=irdr+1
C       ---  Reading the cooling curves   --
        do 1 j=1,nrowb
          read(irdr,*,end=2) a1,a2,a3,a4,a5,a6,a7,a8,
     &                       a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,
     &                       a20,a21,a22,a23,a24,a25,a26,a27
          luminosityDB(i,j)= -(a3-4.72)/2.5
          colorDB_U(i,j)=a5
          colorDB_B(i,j)=a6
          colorDB_V(i,j)=a7
          colorDB_R(i,j)=a8
          colorDB_I(i,j)=a9
1       continue
2       numberOfPointsInSequence(i)=j-1    
3     continue
 
      return
      end
C***********************************************************************










      subroutine colordb2(table)
C=======================================================================
C     This subroutine reads the color tables by Bergeron for the  
C     sequences of WD DB
C-----------------------------------------------------------------------
C     I/O parameters:
C       table: instance of FileGroupInfo
C=======================================================================
      use external_types
      implicit double precision (a-h,m,o-z)

      integer i,j
      double precision a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,
     &                 a15,a16,a17,a18,a19,a20,a21,a22,a23,a24,a25,a26,
     &                 a27
      TYPE(FileGroupInfo) :: table

      table%mass(1)=0.5
      table%mass(2)=0.6
      table%mass(3)=0.7
      table%mass(4)=0.8
      table%mass(5)=0.9
      table%mass(6)=1.0
      table%mass(7)=1.2
 
C     ---   reading the tables   ---
      do 3 i=1,table%ncol
C       ---  Reading unit   ---
        table%initLink=table%initLink+1
C       ---  Reading the cooling curves   --
        do 1 j=1,table%nrow
          read(table%initLink,*,end=2) a1,a2,a3,a4,a5,a6,a7,a8,
     &                       a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,
     &                       a20,a21,a22,a23,a24,a25,a26,a27
          table%luminosity(i,j)= -(a3-4.72)/2.5
          table%color_U(i,j)=a5
          table%color_B(i,j)=a6
          table%color_V(i,j)=a7
          table%color_R(i,j)=a8
          table%color_I(i,j)=a9
1       continue
2       table%ntrk(i)=j-1    
3     continue
 
      return
      end