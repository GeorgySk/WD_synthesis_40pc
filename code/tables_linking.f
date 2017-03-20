      TYPE(FileInfo),DIMENSION(86) :: table

      do 1 i=1,7
        table(i)%initLink=10
        table(i)%link=10+i
        table(i)%tableType="cooling"
        table(i)%WDtype="DA"
        table(i)%flag=1
        table(i)%Z=0.001
        table(i)%ncol=7
        table(i)%nrow=650
1     continue

      do 2 i=8,17
        table(i)%initLink=20
        table(i)%link=13+i
        table(i)%tableType="cooling"
        table(i)%WDtype="DA"
        table(i)%flag=2
        table(i)%Z=0.01
        table(i)%ncol=10
        table(i)%nrow=650
2     continue

      do 3 i=18,25
        table(i)%initLink=30
        table(i)%link=13+i
        table(i)%tableType="cooling"
        table(i)%WDtype="DA"
        table(i)%flag=3
        table(i)%Z=0.03
        table(i)%ncol=8
        table(i)%nrow=650
3     continue

      do 4 i=26,33
        table(i)%initLink=40
        table(i)%link=15+i
        table(i)%tableType="cooling"
        table(i)%WDtype="DA"
        table(i)%flag=3
        table(i)%Z=0.06
        table(i)%ncol=8
        table(i)%nrow=650
4     continue

      do 5 i=34,40
        table(i)%initLink=90
        table(i)%link=57+i
        table(i)%tableType="cooling"
        table(i)%WDtype="DB"
        table(i)%flag=1
        table(i)%Z=0.001
        table(i)%ncol=7
        table(i)%nrow=400
5     continue

      do 6 i=41,49
        table(i)%initLink=100
        table(i)%link=60+i
        table(i)%tableType="cooling"
        table(i)%WDtype="DB"
        table(i)%flag=2
        table(i)%Z=0.01
        table(i)%ncol=9
        table(i)%nrow=400
6     continue

      do 77 i=50,58
        table(i)%initLink=110
        table(i)%link=61+i
        table(i)%tableType="cooling"
        table(i)%WDtype="DB"
        table(i)%flag=3
        table(i)%Z=0.06
        table(i)%ncol=9
        table(i)%nrow=400
77     continue

      do 8 i=59,64
        table(i)%initLink=120
        table(i)%link=62+i
        table(i)%tableType="color"
        table(i)%WDtype="ONe"
8     continue

      do 9 i=65,69
        table(i)%initLink=127
        table(i)%link=62+i
        table(i)%tableType="cooling"
        table(i)%WDtype="ONe"
9     continue

      do 10 i=70,76
        table(i)%initLink=131
        table(i)%link=62+i
        table(i)%tableType="colors"
        table(i)%WDtype="DB"
        table(i)%ncol=7
        table(i)%nrow=60
10     continue

      do 11 i=77,86
        table(i)%initLink=60
        table(i)%link=i-16
        table(i)%tableType="colors"
        table(i)%WDtype="DA"
        table(i)%ncol=10
        table(i)%nrow=650
11     continue