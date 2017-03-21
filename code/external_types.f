      MODULE external_types

      TYPE FileInfo
        SEQUENCE
        CHARACTER(len=7) :: tableType
        CHARACTER(len=3) :: WDtype
        INTEGER :: flag
        REAL :: Z
        INTEGER :: initLink
        INTEGER :: link
        INTEGER :: ncol
        INTEGER :: nrow
C       should be ntrk(ncol) but compiler gets angry. 10 is max of ncol
        INTEGER :: ntrk(10)
C       should be mass(ncol)..        
        DOUBLE PRECISION :: mass(10)
C       should be coolingTime(ncol,nrow)..   
        DOUBLE PRECISION :: coolingTime(10,650)
C       should be prevTime(ncol)..        
        DOUBLE PRECISION :: prevTime(10)
C       should be luminosity(ncol,nrow)..   
        DOUBLE PRECISION :: luminosity(10,650)
      ENDTYPE FileInfo

      END MODULE external_types