C----------------------------------------------------------------------
C The program monte.f simulates the population of white dwarfs in the 
C solar environment. It distributes randomly n points following a 
C uniform distribution in the galactic plane and an exponential 
C distribution in the perpendicular plane. It adds velocity distribution 
C at each point. From the cooling tables it interpolates brightness and 
C cooling times. Version B: In this version, the SFR function, within 
C each interval in which the time of the galaxy is divided, determines 
C the mass of stars to share. Each mass that is generated, following 
C the IMF, is assigned with a birth time according to the SFR. Finally  
C the maximum volume method is used to calculate the density function of 
C white dwarfs.
C----------------------------------------------------------------------
C     adding FileInfo type which carries all the info about
C     files of cooling and color tables: fort.xx links, numbers of rows
C     and columns, metallicities
      include 'code/external_types.f'

      program monte
      use external_types
       
      implicit double precision (a-h,m,o-z)
   
C     'external' statement specifies that 'ran' function is no longer
C     intrinsic and must be defined in program
      external ran
      real ran

C     minimumSectorRadius - min radius of the sector of considered stars
C     maximumSectorRadius - max radius of the sector of considered stars
C     angleCoveringSector - angle in degrees, which covers the sector
C     radiusOfSector: radius (kpc) of the sector centered at Sun
C     TODO give better names or create class (if possible)    
C     zDistribution_zo(zo),heightPattern(h) - parameters of distribution 
C            of z; zo*exp(-z/h)
C     galacticDiskAge (Gyr)
C     parameterOfSFR (taus): parameter of the SFR; Y=exp(-t/taus)
C     solarGalactocentricDistance: distance from Sun to Galaxy center;
C----------------------------------------------------------------------
      integer numberOfStars
      double precision galacticDiskAge,parameterOfSFR,
     &                 solarGalactocentricDistance,minimumSectorRadius,
     &                 maximumSectorRadius,angleCoveringSector,
     &                 parameterIMF
      double precision radiusOfSector,scaleLength,areaOfSector,pi
      
      parameter (numberOfStars=6000000)
      parameter (solarGalactocentricDistance=8.5)
      parameter (minimumSectorRadius=8.45)
      parameter (maximumSectorRadius=8.55)
      parameter (angleCoveringSector=0.674)
      parameter (radiusOfSector=0.050)
      parameter (parameterOfSFR=25.0)
      parameter (scaleLength=3.5)

C     nrow - number of rows in DA color and cooling tables
C     nrowb - number of rows in DB cooling tables
C     nrowb2 - number of rows in DB color tables
      integer nrow,nrowb,nrowb2
      parameter (nrow=650)
      parameter (nrowb=400)
      parameter (nrowb2=60)

C     TODO: move all this to hash-table/map/dictionary/array or smth
C     flags - determine what metallicity
      integer flag_1,flag_2,flag_3
C     these are links to files of DA cooling tables
      integer initialCoolSeqIndex_1,
     &        initialCoolSeqIndex_2,initialCoolSeqIndex_3,
     &        initialCoolSeqIndex_4
C     numbers of columns in DA cooling tables
      integer numberOfMassesWithCoolSeq_1,
     &        numberOfMassesWithCoolSeq_2,numberOfMassesWithCoolSeq_3,
     &        numberOfMassesWithCoolSeq_4
C     number of columns in DA color table
      integer numberOfMassesWithColors
C     these are links to files of DB cooling tables
      integer firstFileOfTheGroup_1,firstFileOfTheGroup_2,
     &        firstFileOfTheGroup_3
C     number of columns in DB cooling tables
      integer numberOfSequencesInGroup_1,
     &        numberOfSequencesInGroup_2,numberOfSequencesInGroup_3
C     number of columns in DB color table
      integer numberOfSequences

C     QUESTION: what are these variables?
C     ntrkda - color DA; numberOfRows - isn't number of rows -DA cooling
      integer ntrkda(10),numberOfRows_1(7),numberOfRows_2(10),
     &        numberOfRows_3(8),numberOfRows_4(8)
C     same as ntrkda and numberOfRows but for DB
      integer vectorOfPointsNumberOfSeq_1(7),
     &        vectorOfPointsNumberOfSeq_2(9),
     &        vectorOfPointsNumberOfSeq_3(9),
     &        numberOfPointsInSequence(7)

      integer i,ISEED1,ISEED2,iseed,numberOfStarsInSample
      double precision randomNumber,fractionOfDB

C     TODO: create classes and take them away from commons
      double precision coolingTimes_1(7,nrow),coolingTimes_2(10,nrow),
     &                 coolingTimes_3(8,nrow)
      double precision coolingTimes_4(8,nrow)
      double precision massOfWD(10),massOfWD_1(7),massOfWD_2(10),
     &                 massOfWD_3(8),massOfWD_4(8)
      double precision luminosity_1(7,nrow),
     &                 effectiveTemperature_1(7,nrow),
     &                 gravitationalAcceleration_1(7,nrow)
      double precision luminosity_2(10,nrow),
     &                 effectiveTemperature_2(10,nrow),
     &                 gravitationalAcceleration_2(10,nrow)
      double precision luminosity_3(8,nrow),
     &                 effectiveTemperature_3(8,nrow),
     &                 gravitationalAcceleration_3(8,nrow)
      double precision luminosity_4(8,nrow),
     &                 effectiveTemperature_4(8,nrow),
     &                 gravitationalAcceleration_4(8,nrow)
      double precision tprewdda1(7),tprewdda2(10),tprewdda3(8)
      double precision tprewdda4(8)
      double precision luminosity(10,nrow),color_U(10,nrow),
     &                 color_B(10,nrow),color_V(10,nrow)
      double precision color_R(10,nrow),color_I(10,nrow)
      double precision vectorOfMasses_1(7),vectorOfMasses_2(9),
     &                 vectorOfMasses_3(9),massSequence(7)
      double precision matrixOfCoolingTimes_1(7,nrowb),
     &                 matrixOfLuminosities_1(7,nrowb)
      double precision matrixOfEffectiveTemperatures_1(7,nrowb),
     &                 matrixOfLog_g_1(7,nrowb)
      double precision matrixOfCoolingTimes_2(9,nrowb),
     &                 matrixOfLuminosities_2(9,nrowb)
      double precision matrixOfEffectiveTemperatures_2(9,nrowb),
     &                 matrixOfLog_g_2(9,nrowb)
      double precision matrixOfCoolingTimes_3(9,nrowb),
     &                 matrixOfLuminosities_3(9,nrowb)
      double precision matrixOfEffectiveTemperatures_3(9,nrowb),
     &                 matrixOfLog_g_3(9,nrowb)
      double precision vectorOfPreviousTimes_1(7),
     &                 vectorOfPreviousTimes_2(9),
     &                 vectorOfPreviousTimes_3(9)
      double precision luminosityDB(7,nrowb2),colorDB_U(7,nrowb2),
     &                 colorDB_B(7,nrowb2)
      double precision colorDB_V(7,nrowb2),colorDB_R(7,nrowb2),
     &                 colorDB_I(7,nrowb2)

C     QUESTION: what is gamma?    
      double precision parameterIFMR,variationOfGravConst,gamma

      TYPE(FileInfo),DIMENSION(86) :: table

C     --- Dimensions of variables of S.Torres  ---
C--------------------------------------------------------------------
C     QUESTION: is z_cyl=z_cart?
C     (coordinate_R,coordinate_Theta,coordinate_Zcylindr): cylindrical 
C         coordinates 
C     (coordinate_X,coordinate_Y,z): coordenadas cartesianas
C     (heliocentricVelocity_U,heliocentricVelocity_V,
C         heliocentricVelocity_W): heliocentric velocities
C     heliocentricVelocities_sigma(3): dispersion of the velocities
C     NOTE: there are a lot of variables with this name
C     TODO: give better names to these variables
C     m: masa de la estrella
C     QUESTION: what is altura patron?
C     heightPattern: altura patron
C--------------------------------------------------------------------      
    
C     ---  Commons S.Torres   ---
       
      common /RSEED/ ISEED1,ISEED2
      
C     ---  Commons Rux ---
C     QUESTION: what is the maximum limit of symbols for groups names?
C     TODO: give better names to common groups
      common /nums/ numberOfMassesWithColors,
     &              numberOfMassesWithCoolSeq_1,
     &              numberOfMassesWithCoolSeq_2,
     &              numberOfMassesWithCoolSeq_3,
     &              numberOfMassesWithCoolSeq_4
      common /nums2/ ntrkda,numberOfRows_1,numberOfRows_2,
     &               numberOfRows_3,numberOfRows_4
      common /masses/ massOfWD,massOfWD_1,massOfWD_2,massOfWD_3,
     &                massOfWD_4
      common /datrks1/ coolingTimes_1,luminosity_1,
     &                 effectiveTemperature_1,
     &                 gravitationalAcceleration_1
      common /datrks2/ coolingTimes_2,luminosity_2,
     &                 effectiveTemperature_2,
     &                 gravitationalAcceleration_2
      common /datrks3/ coolingTimes_3,luminosity_3,
     &                 effectiveTemperature_3,
     &                 gravitationalAcceleration_3
      common /datrks4/ coolingTimes_4,luminosity_4,
     &                 effectiveTemperature_4,
     &                 gravitationalAcceleration_4
      common /datprewd/ tprewdda1,tprewdda2,tprewdda3,tprewdda4
      common /dbnums/ numberOfSequences,numberOfSequencesInGroup_1,
     &                numberOfSequencesInGroup_2,
     &                numberOfSequencesInGroup_3
      common /dbnums2/ vectorOfPointsNumberOfSeq_1,
     &                 vectorOfPointsNumberOfSeq_2,
     &                 vectorOfPointsNumberOfSeq_3,
     &                 numberOfPointsInSequence
      common /massesdb/ massSequence,vectorOfMasses_1,vectorOfMasses_2,
     &                  vectorOfMasses_3
      common /dbtrks1/ matrixOfLuminosities_1,
     &                 matrixOfEffectiveTemperatures_1,matrixOfLog_g_1,
     &                 matrixOfCoolingTimes_1
      common /dbtrks2/ matrixOfLuminosities_2,
     &                 matrixOfEffectiveTemperatures_2,matrixOfLog_g_2,
     &                 matrixOfCoolingTimes_2
      common /dbtrks3/ matrixOfLuminosities_3,
     &                 matrixOfEffectiveTemperatures_3,matrixOfLog_g_3,
     &                 matrixOfCoolingTimes_3
      common /dbtprewd/ vectorOfPreviousTimes_1,vectorOfPreviousTimes_2,
     &                  vectorOfPreviousTimes_3
      common /colors/ luminosity,color_U,color_B,color_V,color_R,color_I
      common /dbcolors/ luminosityDB,colorDB_U,colorDB_B,colorDB_V,
     &                  colorDB_R,colorDB_I
      common /vargra/ variationOfGravConst,gamma
      common /param/ fractionOfDB,galacticDiskAge,parameterIMF,
     &               parameterIFMR,timeOfBurst

      common /tables/ table
      
      include 'code/tables_linking.f'

C     ---  Initialization of parameters of Rux ---
      flag_1=1
      flag_2=2
      flag_3=3
      initialCoolSeqIndex_1=10
      initialCoolSeqIndex_2=20
      initialCoolSeqIndex_3=30
      initialCoolSeqIndex_4=40
      firstFileOfTheGroup_1=90
      firstFileOfTheGroup_2=100
      firstFileOfTheGroup_3=110
      numberOfMassesWithColors=10
      numberOfMassesWithCoolSeq_1=7
      numberOfMassesWithCoolSeq_2=10
      numberOfMassesWithCoolSeq_3=8
      numberOfMassesWithCoolSeq_4=8
      numberOfSequencesInGroup_1=7
      numberOfSequencesInGroup_2=9
      numberOfSequencesInGroup_3=9
      numberOfSequences=7

C     TODO: if input is from file then these parameters will be zeroes
C     ---  Screen output of used parameters  ---
C     NOTE: wrong logic, some of these variables are overwritten later
      write(6,*) '=========================================='
      write(6,*) ' '
      write(6,*) '            Programa monte.f'
      write(6,*) '          by S.Torres, 14.02.11 '
      write(6,*) ' '
      write(6,*) '            Used parameters:'
      write(6,*) 'numberOfStars=',numberOfStars
      write(6,*) 'SFR: parameterOfSFR=',parameterOfSFR,'Gyr'
      write(6,*) 'galacticDiskAge=',galacticDiskAge,'Gyr'
      write(6,*) 'minimumSectorRadius=',minimumSectorRadius,
     &           'kpc; maximumSectorRadius=',maximumSectorRadius,'kpc'
      write(6,*) 'radiusOfSector=',radiusOfSector,'kpc'
      write(6,*) ' '
      write(6,*) '=========================================='
      write(6,*) ' '
      write(6,*) '          Start of calculations:'
      
C     --- Initializing random number generator and reading the seeds ---
      
      iseed=-9
      read(72,100) iseed1,iseed2
      write(6,*) 'iseed1=',iseed1
      write(6,*) 'iseed2=',iseed2

C     QUESTION: why do we need this part?      
      do 8123 i=1,10
      randomNumber=ran(iseed)
      write (6,*) i,randomNumber
8123  continue      

C     ---  Reading free parameters ---
      
      
C ======================================================================
C  
C      fractionOfDB: fraction of DB's
C      galacticDiskAge (Gyr)
C      parameterIMF (alpha): M^{alpha}
C      Initial-to-Final Mass Relation : 
C         mfinal_new=parameterIFMR*mfinal_old
        
C     Reading parameters from the file parameter.dat:
      read (10,*) fractionOfDB,galacticDiskAge,parameterIMF,
     &            parameterIFMR,timeOfBurst

C     Fiducial values (trusted):
C        fractionOfDB=0.20
C        galacticDiskAge=8.9
C        parameterIMF=-2.35
C        parameterIFMR=1.0
C        timeOfBurst=0.6
      fractionOfDB=0.20 
      galacticDiskAge=8.9
      parameterIMF=-2.35
      parameterIFMR=1.0
      timeOfBurst=0.6     

      write (157,157) fractionOfDB,galacticDiskAge,parameterIMF,
     &                parameterIFMR,timeOfBurst
 157  format(5(f6.3,2x))
C                      
C                                                                      |
C ======================================================================
C                                                                      |
C                        Parameters of G   

      variationOfGravConst=-1.0d-11
      gamma=3.6

C                                                                      |
C ======================================================================
C     ---    Calculating the area of the sector  ---

C     QUESTION: isn't there a better way to get Pi?
      pi=4.0*atan(1.0)
C     QUESTION: what about square function?
      areaOfSector=pi*radiusOfSector*radiusOfSector

      

C     ---   Reading the cooling tables  ---

      write(6,*) '1. Reading the cooling tables (1/10)'

      write(6,*) '   1.1 Tracks of CO DA WD Z=0.001;0.01;0.03;0.06'

C     TODO: rename the function 'incoolda'     
C     NOTE: ALl this can be place in one more compact block
C     Calling the function 'incoolda' for 4 metalicities that we have
      call incoolda(table(1)%flag,table(1)%initLink,table(1)%ntrk,
     &     table(1)%ncol,massOfWD_1,
     &     coolingTimes_1,tprewdda1,luminosity_1,effectiveTemperature_1,
     &     gravitationalAcceleration_1)
      call incoolda(table(8)%flag,table(8)%initLink,table(8)%ntrk,
     &     table(8)%ncol,massOfWD_2,coolingTimes_2,
     &     tprewdda2,luminosity_2,effectiveTemperature_2,
     &     gravitationalAcceleration_2)
      call incoolda(table(18)%flag,table(18)%initLink,table(18)%ntrk,
     &     table(18)%ncol,massOfWD_3,coolingTimes_3,
     &     tprewdda3,luminosity_3,effectiveTemperature_3,
     &     gravitationalAcceleration_3)
      call incoolda(table(26)%flag,table(26)%initLink,table(26)%ntrk,
     &     table(26)%ncol,massOfWD_4,coolingTimes_4,
     &     tprewdda4,luminosity_4,effectiveTemperature_4,
     &     gravitationalAcceleration_4)
      
      write(6,*) '   1.2 Tracks of CO non-DA (DB) WD'

C     TODO: rename the function 'incooldb'
      call incooldb(table(34)%flag,table(34)%initLink,
     &     table(34)%ncol,vectorOfPointsNumberOfSeq_1,
     &     vectorOfMasses_1,matrixOfCoolingTimes_1,
     &     vectorOfPreviousTimes_1,matrixOfLuminosities_1,
     &     matrixOfEffectiveTemperatures_1,matrixOfLog_g_1)
      call incooldb(table(41)%flag,table(41)%initLink,
     &     table(41)%ncol,vectorOfPointsNumberOfSeq_2,
     &     vectorOfMasses_2,matrixOfCoolingTimes_2,
     &     vectorOfPreviousTimes_2,matrixOfLuminosities_2,
     &     matrixOfEffectiveTemperatures_2,matrixOfLog_g_2)
      call incooldb(table(50)%flag,table(50)%initLink,
     &     table(50)%ncol,vectorOfPointsNumberOfSeq_3,
     &     vectorOfMasses_3,matrixOfCoolingTimes_3,
     &     vectorOfPreviousTimes_3,matrixOfLuminosities_3,
     &     matrixOfEffectiveTemperatures_3,matrixOfLog_g_3)

      write(6,*) '   1.3 Tracks of ONe DA WD'

C     TODO: rename the function 'incoolone'
      call incoolone
      
      write(6,*) '   1.4 Reading the colors table of Rene(DAs) and Berge
     &ron(DBs)'
C     TODO: rename these functions      
      call color(table(77)%ncol,ntrkda,massOfWD,luminosity,
     &     color_U,color_B,color_R,color_V,color_I)      
      call colordb(table(70)%ncol,numberOfPointsInSequence,
     &     massSequence,luminosityDB,colorDB_U,colorDB_B,colorDB_V,
     &     colorDB_R,colorDB_I)

      write (6,*) '   1.5 Reading the tables of CO DA with G variable'
C     TODO: rename this function      
      call incoolea

C     QUESTION: what is t_b?
C     ---  Calling subroutine calculating IMF, t_b, heightPattern and z 

      write(6,*) '2. Calling the IMF and SFR (2/10)'
C     TODO: rename it
      call gen(iseed,parameterOfSFR,areaOfSector,numberOfStarsInSample,
     &     galacticDiskAge,timeOfBurst)
      
      write(6,*) "numberOfStarsInSample=", numberOfStarsInSample

C     ---   Calling subroutine calculating luminosities ---

      write(6,*) '3. Calculating luminosities (3/10)'
C     TODO: rename it
      call lumx(iseed,numberOfStarsInSample)      

C     ---  Calling subroutine calculating polar coordinates  ---

      write(6,*) '4. Calculating polar coordinates (4/10)'
C     TODO: rename it
      call polar(iseed,minimumSectorRadius,maximumSectorRadius,
     &     angleCoveringSector,radiusOfSector,
     &     solarGalactocentricDistance,scaleLength)

C     ---   Calling subroutine calculating heliocentric velocities ---

      write(6,*) '5. Generating heliocentric velocities (5/10)'
C     TODO: rename it
      call velh(iseed,numberOfStarsInSample)

C     QUESTION: why are we missing the next step?
      goto 7
C     QUESTION: what does this mean?
C     ---   Calculating the trajectories according to/along z-coordinate ---
      
      write(6,*) '6. Integrating trajectories (6/10)'
C     TODO: give a better name
      call traject(galacticDiskAge)

C     ---   Calling subroutine calculating coordinates ---

7     write(6,*) '7. Calculating coordinates (7/10)'
C     TODO: rename this
      call coor(solarGalactocentricDistance)

C     ---   Calling subroutine calculating visual magnitude ---

      write(6,*) '8. Determinating visual magnitudes (8/10)'
C     TODO: rename it
      call magi(fractionOfDB) 

C     ---   Writing the data referring to the WD's   ---

      write(6,*) '9. Generating the Luminosity Function (9/10)'
C     TODO: rename it
      call volum_40pc
  

C     TODO: give a better name to vrad      
C     ---   Testing vrad=0   ---
      
      write(6, *) '10. Making vrad to be null (10/10)'

      call vrado

      write (6,*) 'The end'

C     ---   Formatss   ---

100   format(I6,2x,I6)           

      stop

      end



C***********************************************************************
C***********************************************************************
C***********************************************************************
C***********************************************************************
      include 'code/star_generation/generator.f'
     
      include 'code/cooling/DA/DA_cooling.f'

      include 'code/cooling/unknown/incoolea.f'

      include 'code/colors/DA/byRenedo.f'

      include 'code/luminosities/luminosities.f'

      include 'code/coordinates/polar.f'

      include 'code/velocities/velocities.f'      

      include 'code/math/random_number_generators.f'

      include 'code/trajectories/trajectories.f'
      
      include 'code/coordinates/coor.f'

      include 'code/magnitudes/magi.f'      

      include 'code/magnitudes/DA/interlumda.f'

      include 'code/DA_DB_fraction/dbd_fid.f'

      include 'code/magnitudes/DA/intermag.f'

      include 'code/magnitudes/DB/interlumdb.f'

      include 'code/colors/DB/byBergeron.f'

      include 'code/cooling/ONe/incoolone.f'

      include 'code/cooling/DB/incooldb.f'

      include 'code/magnitudes/ONe/interlumONe.f'

      include 'code/magnitudes/interp.f'

      include 'code/colors/chanco.f'

      include 'code/samples/volum_40pc.f'
      
      include 'code/velocities/vrado.f'

      include 'code/math/toSort.f'