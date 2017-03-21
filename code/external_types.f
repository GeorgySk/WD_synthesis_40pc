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
      ENDTYPE FileInfo

      END MODULE external_types