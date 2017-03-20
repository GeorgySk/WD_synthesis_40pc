      TYPE FileInfo
        CHARACTER(len=7),DIMENSION(1) :: tableType
        CHARACTER(len=3),DIMENSION(1) :: WDtype
        INTEGER,DIMENSION(1) :: flag
        REAL,DIMENSION(1) :: Z
        INTEGER,DIMENSION(1) :: initLink
        INTEGER,DIMENSION(1) :: link
        INTEGER,DIMENSION(1) :: ncol
        INTEGER,DIMENSION(1) :: nrow
      ENDTYPE FileInfo