      integer ii

      TYPE(FileInfo),DIMENSION(86) :: table

      do 1 ii=1,7
        table(ii)%link=10+ii
        table(ii)%tableType="cooling"
        table(ii)%WDtype="DA"
        table(ii)%Z=0.001
        table(ii)%ncol=7
        table(ii)%nrow=650
1     continue

      do 2 ii=8,17
        table(ii)%link=13+ii
        table(ii)%tableType="cooling"
        table(ii)%WDtype="DA"
        table(ii)%Z=0.01
        table(ii)%ncol=10
        table(ii)%nrow=650
2     continue

      do 3 ii=18,25
        table(ii)%link=13+ii
        table(ii)%tableType="cooling"
        table(ii)%WDtype="DA"
        table(ii)%Z=0.03
        table(ii)%ncol=8
        table(ii)%nrow=650
3     continue

      do 4 ii=26,33
        table(ii)%link=15+ii
        table(ii)%tableType="cooling"
        table(ii)%WDtype="DA"
        table(ii)%Z=0.06
        table(ii)%ncol=8
        table(ii)%nrow=650
4     continue

      do 5 ii=34,40
        table(ii)%link=57+ii
        table(ii)%tableType="cooling"
        table(ii)%WDtype="DB"
        table(ii)%Z=0.001
        table(ii)%ncol=7
        table(ii)%nrow=400
5     continue

      do 6 ii=41,49
        table(ii)%link=60+ii
        table(ii)%tableType="cooling"
        table(ii)%WDtype="DB"
        table(ii)%Z=0.01
        table(ii)%ncol=9
        table(ii)%nrow=400
6     continue

      do 77 ii=50,58
        table(ii)%link=61+ii
        table(ii)%tableType="cooling"
        table(ii)%WDtype="DB"
        table(ii)%Z=0.06
        table(ii)%ncol=9
        table(ii)%nrow=400
77     continue

      do 8 ii=59,64
        table(ii)%link=62+ii
        table(ii)%tableType="color"
        table(ii)%WDtype="ONe"
8     continue

      do 9 ii=65 ,69
        table(ii)%link=62+ii
        table(ii)%tableType="cooling"
        table(ii)%WDtype="ONe"
9     continue

      do 10 ii=70,76
        table(ii)%link=62+ii
        table(ii)%tableType="colors"
        table(ii)%WDtype="DB"
        table(ii)%ncol=7
        table(ii)%nrow=60
10     continue

      do 11 ii=77,86
        table(ii)%link=ii-16
        table(ii)%tableType="colors"
        table(ii)%WDtype="DA"
        table(ii)%ncol=10
        table(ii)%nrow=650
11     continue