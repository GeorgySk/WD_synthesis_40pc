C----------------------------------------------------------------------
C The program monte.f simulates the population of white dwarfs in the 
C solar environment. It distributes randomly n points following a 
C uniform distribution in the galactic plane and an exponential 
C distribution in the perpendicular plane. Add a velocity distribution 
C at each point. From the cooling tables interpolates brightness and 
C cooling times. Version B: In this version, the SFR function, within 
C each interval in which the time of the galaxy is divided, determines 
C the mass of stars to share. Each mass that is generated, following 
C the IMF, is assigned with a birth time according to the SFR. Finally  
C the maximum volume method is used to calculate the density function of 
C white dwarfs.
C----------------------------------------------------------------------
      program monte
C
C=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
C   5         15        25        35        45        55        65
C---*->--+----*----+----*----+----*----+----*----+----*----+----*----+->
C     7  10        20        30        40        50        60        70
C=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
       
      implicit double precision (a-h,m,o-z)
   
C     --- Declaration of the variables  ---
C     QUESTION: what does the next 2 lines mean?      
      external ran
      real ran

C     NOTE: variables cant be personal
C     --- Parameters by S.Torres  ---
C----------------------------------------------------------------------

C     numberOfStars: number of stars
C     minimumSectorRadius - min radius of the sector of considered stars
C     maximumSectorRadius - max radius of the sector of considered stars
C     angleCoveringSector - angle in degrees, which covers the sector
C     TODO give better names or create class (if possible)    
C     zDistribution_zo(zo),heightPattern(h) - parameters of distribution 
C            of z; zo*exp(-z/h)
C     galacticDiskAge: galactic disk age; (Gyr)
C     NOTE zi,zf,tauh are never used + give better names or create class    
C     hDistribution_zi(zi),hDistribution_zf(zf),hDistribution_t(tauh): 
C           parameters of distribution of heightPattern(h)
C           h=zi*exp(-t/tauh)+zf
C         parameterOfSFR (taus): parameter of the SFR; Y=exp(-t/taus)
C     solarGalactocentricDistance: distance from Sun to Galaxy center;
C     radiusOfSector: radius (kpc) of the sector centered at Sun
C----------------------------------------------------------------------
      integer numberOfStars
      double precision galacticDiskAge,hDistribution_zi,hDistribution_t,
     &                 hDistribution_zf,parameterOfSFR,
     &                 solarGalactocentricDistance,minimumSectorRadius,
     &                 maximumSectorRadius,angleCoveringSector,
     &                 parameterIMF
      double precision radiusOfSector,scaleLength,areaOfSector,pi
      
      parameter (numberOfStars=6000000)
      parameter (hDistribution_zi=242.5,hDistribution_t=0.7,
     &          hDistribution_zf=0.500) 
      parameter (parameterOfSFR=25.0)
      parameter (solarGalactocentricDistance=8.5)
      parameter (minimumSectorRadius=8.45)
      parameter (maximumSectorRadius=8.55)
      parameter (angleCoveringSector=0.674)
      parameter (radiusOfSector=0.050)
      parameter (scaleLength=3.5)

C     ---  Parameters of Rux ---
C     QUESTION: what are nrow, nrowb, nrowb2?
C     TODO: rename these variables
      integer nrow,nrowb,nrowb2
      parameter (nrow=650)
      parameter (nrowb=400)
      parameter (nrowb2=60)

C     ---  Declaration of variables of Rux ---
C     QUESTION: is it necessary to keep the flags?
C     NOTE: there is no need to keep these indexes (they are used to 
C           open files)
      integer flag_1,flag_2,flag_3,initialCoolSeqIndex_1,
     &        initialCoolSeqIndex_2,initialCoolSeqIndex_3,
     &        initialCoolSeqIndex_4
C     TODO: give better names to these variables + understand meaning 
      integer numberOfMassesWithColors,numberOfMassesWithCoolSeq_1,
     &        numberOfMassesWithCoolSeq_2,numberOfMassesWithCoolSeq_3,
     &        numberOfMassesWithCoolSeq_4,numberOfSequencesInGroup_1,
     &        numberOfSequencesInGroup_2,numberOfSequencesInGroup_3,
     &        numberOfSequences
C     TODO: get rid of these variables
      integer firstFileOfTheGroup_1,firstFileOfTheGroup_2,
     &        firstFileOfTheGroup_3
C     QUESTION: what is ntrkda?
C     numberOfRows_1/2/3/4 - number of rows according to reference 
C                            time vector
C     QUESTION: what is reference time vector?
      integer ntrkda(10),numberOfRows_1(7),numberOfRows_2(10),
     &        numberOfRows_3(8),numberOfRows_4(8)
C     TODO: rename all these + figure out indexes meanings
      integer vectorOfPointsNumberOfSeq_1(7),
     &        vectorOfPointsNumberOfSeq_2(9),
     &        vectorOfPointsNumberOfSeq_3(9),
     &        numberOfPointsInSequence(7)
      integer i,ISEED1,ISEED2,iseed,numberOfStarsInSample
      double precision randomNumber,fractionOfDB
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
C     QUESTION: what is the meaning of the next values? 
C     TODO: rename them 
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
       
      write(6,*) '=========================================='
      write(6,*) ' '
      write(6,*) '            Programa monte.f'
      write(6,*) '          by S.Torres, 14.02.11 '
      write(6,*) ' '
      write(6,*) '            Used parameters:'
      write(6,*) 'numberOfStars=',numberOfStars
      write(6,*) 'SFR: parameterOfSFR=',parameterOfSFR,'Gyr'
      write(6,*) 'IMF: hDistribution_zi=',hDistribution_zi,
     &           'Kpc; hDistribution_t=',hDistribution_t,
     &           'Gyr; hDistribution_zf=',hDistribution_zf,'kpc'
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
      call incoolda(flag_1,initialCoolSeqIndex_1,numberOfRows_1,
     &     numberOfMassesWithCoolSeq_1,massOfWD_1,
     &     coolingTimes_1,tprewdda1,luminosity_1,effectiveTemperature_1,
     &     gravitationalAcceleration_1)
      call incoolda(flag_2,initialCoolSeqIndex_2,numberOfRows_2,
     &     numberOfMassesWithCoolSeq_2,massOfWD_2,coolingTimes_2,
     &     tprewdda2,luminosity_2,effectiveTemperature_2,
     &     gravitationalAcceleration_2)
      call incoolda(flag_3,initialCoolSeqIndex_3,numberOfRows_3,
     &     numberOfMassesWithCoolSeq_3,massOfWD_3,coolingTimes_3,
     &     tprewdda3,luminosity_3,effectiveTemperature_3,
     &     gravitationalAcceleration_3)
      call incoolda(flag_3,initialCoolSeqIndex_4,numberOfRows_4,
     &     numberOfMassesWithCoolSeq_4,massOfWD_4,coolingTimes_4,
     &     tprewdda4,luminosity_4,effectiveTemperature_4,
     &     gravitationalAcceleration_4)
      
      write(6,*) '   1.2 Tracks of CO non-DA (DB) WD'

C     TODO: rename the function 'incooldb'
      call incooldb(flag_1,firstFileOfTheGroup_1,
     &     numberOfSequencesInGroup_1,vectorOfPointsNumberOfSeq_1,
     &     vectorOfMasses_1,matrixOfCoolingTimes_1,
     &     vectorOfPreviousTimes_1,matrixOfLuminosities_1,
     &     matrixOfEffectiveTemperatures_1,matrixOfLog_g_1)
      call incooldb(flag_2,firstFileOfTheGroup_2,
     &     numberOfSequencesInGroup_2,vectorOfPointsNumberOfSeq_2,
     &     vectorOfMasses_2,matrixOfCoolingTimes_2,
     &     vectorOfPreviousTimes_2,matrixOfLuminosities_2,
     &     matrixOfEffectiveTemperatures_2,matrixOfLog_g_2)
      call incooldb(flag_3,firstFileOfTheGroup_3,
     &     numberOfSequencesInGroup_3,vectorOfPointsNumberOfSeq_3,
     &     vectorOfMasses_3,matrixOfCoolingTimes_3,
     &     vectorOfPreviousTimes_3,matrixOfLuminosities_3,
     &     matrixOfEffectiveTemperatures_3,matrixOfLog_g_3)

      write(6,*) '   1.3 Tracks of ONe DA WD'

C     TODO: rename the function 'incoolone'
      call incoolone
      
      write(6,*) '   1.4 Reading the colors table of Rene(DAs) and Berge
     &ron(DBs)'
C     TODO: rename these functions      
      call color(numberOfMassesWithColors,ntrkda,massOfWD,luminosity,
     &     color_U,color_B,color_R,color_V,color_I)      
      call colordb(numberOfSequences,numberOfPointsInSequence,
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





C     TODO: give better names
      subroutine incoolea
C=======================================================================
C
C     This subroutine reads the cooling tables by Leandro Althaus
C
C-----------------------------------------------------------------------
C     Input parameters
C
C       none
C
C-----------------------------------------------------------------------
C     Output parameters
C
C       luminosity:   Luminosuty of the WD. [log(L/L0)]
C       mtrk:  Mass of the WD. [M0]
C       ttrk:  Cooling time. [Gyr]
C
C=======================================================================
      implicit double precision (a-h,o-z)

C     ---   Declaration of variables   ---
      double precision mtrk,luminosity

C     ---   Parameters   ---
      parameter (ncol=3)
      parameter (nbank=3)
      parameter (nrow=900)

C     ---   Dimensions   ---
      dimension mtrk(ncol),ntrk(ncol,nbank)
      dimension ttrk(ncol,nrow,nbank),luminosity(ncol,nrow,nbank)
      dimension ginic(nbank)

C     ---   Commons   ---
      common /tracks/ ginic,luminosity,mtrk,ttrk,ntrk
      
C     --- Values of the mass and G for A=-1.1d-11 ---
      mtrk(1)=0.52
      mtrk(2)=0.6
      mtrk(3)=1.0

      ginic(1)=1.020
      ginic(2)=1.050
      ginic(3)=1.100

C     ---   Initialization   ---
      irdr=299
 
C     ---   Reading the tables  Z=0.05 ---
      do 4 k=1,nbank 
        do 3 i=1,ncol
C         ---   Reading unit   ---
          irdr=irdr+1
C         ---   Reading the cooling curves   ---
          do 1 j=1,nrow
C           QUESTION: what does end=2 mean?            
            read(irdr,*,end=2) luminosity(i,j,k),a2,a3,a4,ttrk(i,j,k),
     &                         a6,a7,a8,a9,goverg
            ttrk(i,j,k)=10.0**(ttrk(i,j,k)-3.0)
    1     continue
    2     ntrk(i,k)=j-1
    3   continue 
    4 continue

      return
      end
C***********************************************************************






C***********************************************************************
C     TODO: understand and reformat
      subroutine color(numberOfMassesWithColors,ntrk,massOfWD,
     &           luminosity,color_U,color_B,color_R,color_V,color_I)
C=======================================================================
C
C     This subroutine reads the colors by Rene and interpolates 
C     according to the vector of reference time
C     
C     Modification at 07/2012 (ER Cojocaru)
C
C-----------------------------------------------------------------------
C
C     Inpu parameters
C       QUESTION: what is it supposed to mean?
C       ntrk: number of elements vector of reference time
C
C-----------------------------------------------------------------------
C
C     Output parameters:
C
C       numberOfMassesWithColors: number of masses for which colors are 
C                                 calculated
C       massOfWD: mass of the WD. [M0]
C       luminosity: luminosity [log(L/L_0)]
C       QUESTION: what does M_ mean?
C       color_U: M_U
C       color_B: M_B
C       color_V: M_V
C       color_R: M_R
C       color_I: M_I
C
C=======================================================================

      implicit double precision (a-h,m,o-z)

C     ---   Parameters   ---
      integer nrow
      parameter (nrow=650)
      integer irdr
      parameter (irdr=60)
      double precision vtanmin
      parameter (vtanmin=30)

C     ---   Declaration of variables   ---
      integer i,k,numberOfMassesWithColors,
     &        ntrk(numberOfMassesWithColors)
      double precision a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14
      double precision a15,a16,a17,a18,a19,a20,a21,a22,a23,a24,a25,a26
      double precision a27,a28,massOfWD(numberOfMassesWithColors)
      double precision tprewdda(numberOfMassesWithColors),
     &                 ttrk(numberOfMassesWithColors,nrow)
      double precision luminosity(numberOfMassesWithColors,nrow),
     &                 color_U(numberOfMassesWithColors,nrow),
     &                 color_B(numberOfMassesWithColors,nrow)
      double precision color_V(numberOfMassesWithColors,nrow),
     &                 color_R(numberOfMassesWithColors,nrow),
     &                 color_I(numberOfMassesWithColors,nrow)
      double precision Hg(nrow),ggii(nrow),Hbj(nrow),bjr(nrow)
      double precision g,gr,bj,term

C     ---  Commons  ---
C     Curve to delimit RPMD, for Mwd=0.61Mo
      common /Hgcurve/ Hg,ggii
      common /Hbjcurve/ Hbj,bjr

C     read masses
      massOfWD(1)=0.524
      massOfWD(2)=0.570
      massOfWD(3)=0.593
      massOfWD(4)=0.610
      massOfWD(5)=0.632
      massOfWD(6)=0.659
      massOfWD(7)=0.705
      massOfWD(8)=0.767
      massOfWD(9)=0.837
      massOfWD(10)=0.877
      
C     --- Tprew_WD LPCODE ---
      tprewdda(1)=11.117
      tprewdda(2)=2.7004
      tprewdda(3)=1.699
      tprewdda(4)=1.2114
      tprewdda(5)=0.9892
      tprewdda(6)=0.7422
      tprewdda(7)=0.4431
      tprewdda(8)=0.2686
      tprewdda(9)=0.200
      tprewdda(10)=0.114

C     read values from files
      do k=1,numberOfMassesWithColors
        do i=1,nrow
          read(irdr+k,*,end=15) a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,
     &                          a13,a14,a15,a16,a17,a18,a19,a20,a21,a22,
     &                          a23,a24,a25,a26,a27,a28
          luminosity(k,i)=a3
          ttrk(k,i)=a4-tprewdda(k)
          color_U(k,i)=a24
          color_B(k,i)=a25
          color_V(k,i)=a26
          color_R(k,i)=a27
          color_I(k,i)=a28
        enddo
15      ntrk(k)=i-1
      end do
      
      term=5*log10(vtanmin)-3.379
      do i=1,ntrk(4)
C       SDSS
        g=color_V(4,i)+0.630*(color_B(4,i)-color_V(4,i))-0.124
C       Hg=Mg+5*log10(v_tan)-3.379
        Hg(i)=g+term
        ggii(i)=1.646*(color_V(4,i)-color_R(4,i))+1.007*(color_R(4,i)-
     &          color_I(4,i))-0.375
C       SuperCosmos
        gr=1.646*(color_V(4,i)-color_R(4,i))-0.139
        bj=0.15+0.13*gr+g
C       Hbj=Mbj+5*log10(v_tan)-3.379
        Hbj(i)=bj+term
        bjr(i)=0.28+1.13*gr
      end do
      
      return
      end
C***********************************************************************


C***********************************************************************
C     TODO: delete iseed from dummies + rewrite all this
      subroutine lumx(iseed,numberOfStarsInSample)
C=======================================================================
C     This subroutine determines what stars are WDs and calculates
C     the cooling time and the luminosity of it.
C
C     Revised at 26.09.07 by S. Torres
C     Introduced metalicity 08/2012 (ER Cojocaru)
C     Adapted for G_Var 14.05.15 by S.Torres
C-----------------------------------------------------------------------
C     input parameters:
C       galacticDiskAge
C       numberOfStars
C       numberOfStarsInSample
C-----------------------------------------------------------------------
C     output parameters
C       leb: luminosities of the WDs
C       ntwd: total number of the WDs
C=======================================================================
      implicit double precision (a-h,m,o-z)

C     ---   Definition of variables  ---
      integer iseed,numberOfStarsInSample,ntwd,numberOfStars,ntwdone,
     &        igorda,i,k,param
      double precision galacticDiskAge,mebmin,mebmax,xntwd,mone,
     &                 parameterIFMR

C     ---   Parameters   ---
      parameter (numberOfStars=6000000)
      parameter (mone=1.14)

C     ---   Declaration of variables   ---
C-----------------------------------------------------------------------
C       starBirthTime: time of birth of the star
C       tms: lifetime
C       tcool: time of cooling
C       leb: luminosity of the WD
C       meb: mass of the WD
C       iwd: 0 - it's not WD, 1 - it's a WD
C       m: mass in the main sequence
C       ntwd: total number of WDs
C-----------------------------------------------------------------------
      double precision starBirthTime(numberOfStars),tms(numberOfStars),
     &                 tcool(numberOfStars)
      double precision leb(numberOfStars),meb(numberOfStars),
     &                 zeb(numberOfStars),teb(numberOfStars)
      double precision iwd(numberOfStars)
      double precision m(numberOfStars)
      double precision ztcool(numberOfStars),zmeb(numberOfStars),
     &                 zzeb(numberOfStars),ztms(numberOfStars)
      double precision ztborn(numberOfStars)

C     ---  Commons  ---
      common /tm/ starBirthTime,m
      common /enanas/ leb,meb,zeb,teb
      common /index/ iwd,ntwd
      common /cool/ tcool
      common /tms/ tms
      common /param/ fractionOfDB,galacticDiskAge,parameterIMF,
     &               parameterIFMR,timeOfBurst

C-----------------------------------------------------------------------
C     ---  Deciding what stars are WDs  ---
C-----------------------------------------------------------------------
      mebmin=0.2
      mebmax=1.2
      ntwd=0
      ntwdone=0
      igorda=0      
      param=1

C     ---  Deciding if the star is a WD  ---
      do 1 i=1,numberOfStarsInSample
        iwd(i)=0
C       ---  Deciding if the star is a WD  ---
C       Progenitor star that generates a ONe WD: 8.5<M_MS<10.5
C       WD of CO: m_WD<1.14; of ONe: m_wd>1.14
        if (m(i).le.10.5) then
C         ---  Attributing a solar metallicity to all the stars ---
          zeb(i)=0.01
C         ---  Calculating the lifetime in the main sequence ---
          call tsp(m(i),zeb(i),tms(i))
C         ---  Calculating of the cooling time  ---
          tcool(i)=galacticDiskAge-starBirthTime(i)-tms(i)
          if (tcool(i).gt.0.0) then
C           ---- IFMR -----
            call mmswd(m(i),meb(i))
C           Using Z solar z=0.01 
            write (667,*) m(i),meb(i)
            meb(i)=parameterIFMR*meb(i)
            if(meb(i).le.1.4) then 
              iwd(i)=1
              ntwd=ntwd+1
              if(meb(i).gt.mone) then
                ntwdone=ntwdone+1
              endif
            else
              iwd(i)=0
            endif
          else
            iwd(i)=0
          endif
        else
          iwd(i)=0
        endif
 1    continue

      write (6,*) '******** Data   ***********'
      write (6,*) ' Number of WDs: ', ntwd
      xntwd=dfloat(ntwd)
      write (6,*) ' Number of ONe: ', ntwdone
      write (6,*) ' ONe percentage: ',100.0*dfloat(ntwdone)/xntwd, '%'

C     ---   Taking the stars that are WDs ---
      do 2 i=1,numberOfStarsInSample
        ztcool(i)=tcool(i)
        zmeb(i)=meb(i)
        zzeb(i)=zeb(i)
        ztborn(i)=starBirthTime(i)
        ztms(i)=tms(i)
 2    continue

      k=0
C     ---   Making the transfer   ---
      do 3 i=1,numberOfStarsInSample
        if (iwd(i).eq.1) then
          k=k+1
          tcool(k)=ztcool(i)
          meb(k)=zmeb(i)
          zeb(k)=zzeb(i)
          starBirthTime(k)=ztborn(i)
          tms(k)=ztms(i)
          write (81,81) meb(k),zeb(k),starBirthTime(k),tms(k),tcool(k)
        endif
 3    continue

 81   format (5(f7.4,2x))

C     ---   Writing data   ---
      write(6,*) '      Total number of WDs ntwd=',ntwd
      write(6,*) '      WDs of ONe=',ntwdone

      return
      end
C***********************************************************************



C     TODO: rewrite 
      subroutine tsp(m,z,t)
C===================================================================
C
C     This subroutine calculates the lifetime in the main sequence
C     for a given metallicity Z€[0.01,0.001]
C     for standart helium content
C     According to model by Leandro(private com.) & Renedo et al.(2010)
C     Data in solar masses and Gyr
C
C-------------------------------------------------------------------
C     Input parameters:
C       m: mass of the star
C       Z: metallicity
C-------------------------------------------------------------------
C     Output parameters:
C       t: lifetime in the SP.
C
C====================================================================
      implicit double precision (a-h,m,o-z)

C     ---   Dimensions  ---
      dimension mms(10),tms(10)
      dimension mms2(7),tms2(7)

C-------------------------------------------------------------------
C     ---  Table of values Z solar --
C-------------------------------------------------------------------
      mms(1)=1.00
      mms(2)=1.50    
      mms(3)=1.75   
      mms(4)=2.00   
      mms(5)=2.25 
      mms(6)=2.50 
      mms(7)=3.00 
      mms(8)=3.50  
      mms(9)=4.00 
      mms(10)=5.00   

C     Althaus priv. comm X=0.725, Y=0.265
      tms(1)=8.614
      tms(2)=1.968
      tms(3)=1.249  
      tms(4)=0.865  
      tms(5)=0.632
      tms(6)=0.480
      tms(7)=0.302
      tms(8)=0.226 
      tms(9)=0.149
      tms(10)=0.088   

C     ---  Interpolating ---
      if (m.lt.mms(1)) then
C       --- Mass less than the first, linear extrapolation, 2 last 
C           points 
        pen=(tms(2)-tms(1))/(mms(2)-mms(1))
        tsol=pen*m+(tms(1)-pen*mms(1))
      else
        if (m.gt.mms(10)) then
C       --- Mass greater than the last, taking the fraction last point
          tsol=(mms(10)/m)*tms(10)
        else
C         QUESTION:--- Interpolation properly/itself? ---    
          k=1
  1       k=k+1
          if(m.lt.mms(k)) then
            pen=(tms(k)-tms(k-1))/(mms(k)-mms(k-1))
            tsol=pen*m+(tms(k)-pen*mms(k))
          else
            goto 1
          endif 
        endif
      endif
C-------------------------------------------------------------------
C     ---  Tabla of values Z Sub-Solar --
C-------------------------------------------------------------------
      mms2(1)=0.85
      mms2(2)=1.00    
      mms2(3)=1.25   
      mms2(4)=1.50   
      mms2(5)=1.75 
      mms2(6)=2.00 
      mms2(7)=3.00 

C Althaus priv. comm X=0.752, Y=0.247
      tms2(1)=10.34
      tms2(2)=5.756
      tms2(3)=2.623  
      tms2(4)=1.412  
      tms2(5)=0.905
      tms2(6)=0.639
      tms2(7)=0.245

C     ---  Interpolating ---
      if (m.lt.mms2(1)) then
C       --- Mass less than the first, linear extrapolation, 2 last 
C           points 
        pen=(tms2(2)-tms2(1))/(mms2(2)-mms2(1))
        tsub=pen*m+(tms2(1)-pen*mms2(1))
      else
        if (m.gt.mms2(7)) then
C         --- Masa greater than the last, extrapolating 2 last points
          tsub=(mms(7)/m)*tms(7)
        else
C         QUESTION:--- Interpolation properly/itself? ---   
          k=1
  2       k=k+1
          if(m.lt.mms2(k)) then
            pen=(tms2(k)-tms2(k-1))/(mms2(k)-mms2(k-1))
            tsub=pen*m+(tms2(k)-pen*mms2(k))
          else
            goto 2
          endif 
        endif
      endif

C-------------------------------------------------------------------
C     ---  Interpolating for the value of Z --
C     z solar 10
C-------------------------------------------------------------------
      t=tsub+((tsol-tsub)/(0.01-0.001))*(z-0.001)

      return
      end
C********************************************************************
      




C     TODO: rewrite
      subroutine mmswd(mass,meb)
C ======================================================================
C     IFMR according to model by Catalán et al.2008
C     combination with the model by Iben for Mi>6Mo
C=======================================================================
      implicit double precision (a-h,m,o-z)

      double precision mass,meb
      
      if(mass.lt.2.7)then
        meb=0.096d0*mass+0.429d0
      elseif((mass.ge.2.7).and.(mass.le.6.0))then
C       Small correction for continuity
        meb=0.137d0*mass+0.3183d0
      else
C       Slope of Iben + continuity in Mi=6Mo
        meb=0.1057d0*mass+0.5061d0
      end if

      return
      end
C***********************************************************************




      
C***********************************************************************
C     TODO: rewrite      
      subroutine polar(iseed,minimumSectorRadius,maximumSectorRadius,
     &           angleCoveringSector,radiusOfSector,
     &           solarGalactocentricDistance,scaleLength)
C=======================================================================
C
C     This subrutina generates the positiones, in polar coordinates, of
C     the stars.
C   
C     Revised in 22.09.07 by S. Torres
C
C-----------------------------------------------------------------------
C     Input parameters:
C     minimumSectorRadius: minimum radius of the sector; in Kpc from the 
C                          Galactic Center
C     maximumSectorRadius: maximum radius
C     angleCoveringSector: angle covering the sector in degrees;
C                          from the GC
C     QUESTION: did i fail to name this variable correctly?
C     radiusOfSector: radial distance to the Sun
C     solarGalactocentricDistance: galactocentric distance of the Sun
C     scaleLength: scale length
C     ntwd: total number of WDs
C=======================================================================
      implicit double precision (a-h,m,o-z)

C     --- Declaration of variables  ---
      external ran       
      real ran
      integer ntwd,iseed,j
      double precision dospi,minimumSectorRadius,maximumSectorRadius,
     &                 angleCoveringSector,radiusOfSector,
     &                 solarGalactocentricDistance,scaleLength,asr
      double precision drsun2,dist,pi,xmin,xmax,zzz,zzr,zzy,zz,xx,xc,yc
      double precision xcte,xinc
       
C     ---  Parameters  ---
      integer numberOfStars
      parameter (numberOfStars=6000000)
           
C     ---  Dimensions  ---
      double precision coordinate_R(numberOfStars),
     &                 coordinate_Theta(numberOfStars),
     &                 coordinate_Zcylindr(numberOfStars)
      double precision x(numberOfStars),y(numberOfStars)
      double precision iwd(numberOfStars)

C     ---  Commons  ---
      common /coorcil/ coordinate_R,coordinate_Theta,coordinate_Zcylindr
      common /plano/ x,y
      common /index/ iwd,ntwd

C     ---  Inicialization of pi and sigma ---
      pi=4.0*atan(1.0)
      dospi=2.0*pi
      drsun2=radiusOfSector*radiusOfSector
      
C     --- Calculating the angle in the sector
C         -angleCoveringSector/2 y +angleCoveringSector/2 degrees
C      and radius between minimumSectorRadius and maximumSectorRadius ---
      asr=(angleCoveringSector*pi)/180.0
      xmax=(maximumSectorRadius*maximumSectorRadius)
      xmin=(minimumSectorRadius*minimumSectorRadius)
      xcte=xmax-xmin
      xinc=maximumSectorRadius-minimumSectorRadius
                
      do 2 j=1,ntwd
3       coordinate_Theta(j)=asr*ran(iseed)-(asr/2)
        if (coordinate_Theta(j).lt.0.0) then
          coordinate_Theta(j)=dospi+coordinate_Theta(j)
        endif
31      zzz=minimumSectorRadius+xinc*ran(iseed)
        zzy=0.16*ran(iseed)
        zzr=dexp(-zzz/scaleLength)
        if (zzy.gt.zzr) then
          goto 31
        else
        endif
        zz=(zzz-minimumSectorRadius)/xinc
        xx=xmin+(xcte*zz)
        coordinate_R(j)=dsqrt(xx)  
        xc=coordinate_R(j)*dcos(coordinate_Theta(j))
        yc=coordinate_R(j)*dsin(coordinate_Theta(j))
        dist=((xc-solarGalactocentricDistance)*
     &       (xc-solarGalactocentricDistance)+yc*yc)
C       QUESTION: what does it mean?       
C       --- Sol no hay más que uno ---
        if (dist.gt.drsun2.or.dist.lt.0.0000015) then 
          goto 3
        endif
        x(j)=xc
        y(j)=yc
 2    continue

      return
      end
C***********************************************************************




C***********************************************************************
C     TODO: rewrite      
      subroutine velh(iseed,numberOfStarsInSample)
C=======================================================================
C
C     This subroutine calculates the heliocentrical velocities. 
C     From the height pattern it calculates the dispersions.
C
C     Revised in 22.09.07 by S. Torres
C
C--------------------------------------------------------------------
C     Input parameters:
C       NOTE: I should give iseed following descriptive name
C       iseed: random number generator parameter
C       numberOfStarsInSample
C=======================================================================
      implicit double precision (a-h,m,o-z)
      
      integer numberOfStars,iseed,numberOfStarsInSample,i,k,ntwd
      double precision uo,vo,wo,a,b,solarGalactocentricDistance,uop,vop,
     &                 yy,gasdev,uom,vom
          
C     ---   Parameters   ---
C--------------------------------------------------------------------
C     a,b: Oort constants (values by Kerr and Lynden-Bell 1986)
C          a=14.4 Km/sKpc
C          b=-12.8 Km/sKpc
C     solarGalactocentricDistance=8.5 Kpc
C     uo,vo,wo: peculiar velocities of the Sun
C--------------------------------------------------------------------
      parameter (numberOfStars=6000000)
      parameter (a=14.4,b=-12.8,solarGalactocentricDistance=8.5d0)
      
C     ---   Dimensiones   ---
C--------------------------------------------------------------------
C     (uu,vv,ww): heliocentric velocities, B3 system
C     sigma(3): dispersion of velocities
C--------------------------------------------------------------------
      double precision uu(numberOfStars),vv(numberOfStars),
     &                 ww(numberOfStars)
      double precision sigma(3)
      double precision coordinate_R(numberOfStars),
     &                 coordinate_Theta(numberOfStars),
     &                 coordinate_Zcylindr(numberOfStars)
      double precision heightPattern(numberOfStars)
      double precision iwd(numberOfStars)
      double precision zz(numberOfStars),zh(numberOfStars)
      
C     ---   Commons  ---
      common /vel/ uu,vv,ww
      common /patron/ heightPattern
      common /coorcil/ coordinate_R,coordinate_Theta,coordinate_Zcylindr
      common /index/ iwd,ntwd
      
      uo=-10.0
      vo=-5.2
      wo=-7.2

C     ---  Making the transfer of heightPattern, of z  ---
      do 2 i=1,numberOfStarsInSample
        zh(i)=heightPattern(i)
        zz(i)=coordinate_Zcylindr(i)
2     continue

      k=0      
      do 3 i=1,numberOfStarsInSample
        if (iwd(i).eq.1) then
          k=k+1
          heightPattern(k)=zh(i)
          coordinate_Zcylindr(k)=zz(i)
        else
      endif
3     continue

      do 1 i=1,ntwd
C-------------------------------------------------------------------
C       ---  Calculating the dispersions  ---
C-------------------------------------------------------------------      
C       --  THIN DISK  --
        sigma(3)=dsqrt(heightPattern(i)/(6.25d-4))          
        sigma(1)=(dsqrt(2.0d0))*sigma(3)
        sigma(2)=(dsqrt(0.32+(1.67d-5)*sigma(1)*sigma(1)))*sigma(1)
C       ---   Calling to the function of gaussian distribution  ---
        yy=gasdev(iseed)
        uop=uom(coordinate_R(i),coordinate_Theta(i),a,b,
     &      solarGalactocentricDistance,uo)
        uu(i)=sigma(1)*yy+uop
        yy=gasdev(iseed)
        vop=vom(coordinate_R(i),coordinate_Theta(i),a,b,
     &      solarGalactocentricDistance,vo)
        vv(i)=sigma(2)*yy+vop-sigma(1)*sigma(1)/120.0
        yy=gasdev(iseed)
        ww(i)=sigma(3)*yy+wo
        write(158,*) i,uu(i), vv(i), ww(i)
        write(400,*) uu(i)
        write(401,*) vv(i)
        write(402,*) ww(i)
1     continue

      return
      end
C***********************************************************************




C***********************************************************************
C     TODO: rewrite
      function gasdev(iseed)
C=======================================================================
C
C     Returns a normally distributed deviate with zero mean and unit
C     variance.
C
C=======================================================================
      implicit double precision (a-h,m,o-z)
      external ran
      real ran
C     QUESTION: what is save?
      save
      
      integer iseed,iset
      double precision v1,v2,r,fac,gset,gasdev
      
      data iset/0/
      if (iset.eq.0) then
  1     v1=2.0d0*(ran(iseed))-1.0d0
        v2=2.0d0*(ran(iseed))-1.0d0
        r=v1*v1+v2*v2
        if (r.ge.1.0d0.or.r.eq.0.0d0) goto 1
        fac=dsqrt(-2.0d0*dlog(r)/r)
        gset=v1*fac
        gasdev=v2*fac
        iset=1
      else
        gasdev=gset
        iset=0
      endif

      return
      end
C***********************************************************************      
      


      
C***********************************************************************
      function uom(r,th,a,b,solarGalactocentricDistance,uo)
C=======================================================================
C
C     Calculating uo taking into account the effect of galactic rotation
C
C=======================================================================
      implicit double precision (a-h,m,o-z)
      
      double precision r,th,a,b,solarGalactocentricDistance,uo,uom
      
      uom=uo+((3.0-(2.0*r)/solarGalactocentricDistance)*a-b)*r*sin(th)
      
      return
      end
C***********************************************************************



C***********************************************************************
      function vom(r,th,a,b,solarGalactocentricDistance,vo)
C=======================================================================
C
C     Calculating vo taking into account the effect of galactic rotation
C
C=======================================================================
      implicit double precision (a-h,m,o-z)
      
      double precision r,th,a,b,solarGalactocentricDistance,vo,vom

      vom=vo+((3.0-(2.0*r)/solarGalactocentricDistance)*a-b)*r*cos(th)-
     &    (a-b)*solarGalactocentricDistance

      return
      end
C***********************************************************************


C     NOTE: the following subroutine is actually missed 
C***********************************************************************
C     TODO: rewrite
      subroutine traject(galacticDiskAge)
C=======================================================================
C     This subroutine calculates the trajectories of the WDs according 
C     to the z-axis. Using a 4th order Runge-Kuttta.
C
C     Revised in 22.09.07. by S. Torres
C
C-------------------------------------------------------------------
C     Input parameters
C       galacticDiskAge
C       numberOfStarsInSample
C=======================================================================
      implicit double precision (a-h,m,o-z)
      
      integer numberOfStars,ntwd,njumps,n,i,NOK,NBAD
      double precision galacticDiskAge,wosun,hmin,eps,fcgys,tf,xcar,ycar
      double precision wo,zDistribution_zo,ecini,tb,tinc,htry,ti,ecinf
      double precision epotf,epoti,f

C     ---   Parameters  ---
      parameter (numberOfStars=6000000)
      parameter (wosun=-8.0d0)
      
C     ---   Dimensions   ---
      double precision uu(numberOfStars),vv(numberOfStars),
     &                 ww(numberOfStars)
      double precision coordinate_R(numberOfStars),
     &                 coordinate_Theta(numberOfStars),
     &                 coordinate_Zcylindr(numberOfStars)
      double precision starBirthTime(numberOfStars),m(numberOfStars)
      double precision iwd(numberOfStars)
      double precision yscal(2),y(2),dydx(2)         
      double precision xpla(numberOfStars),ypla(numberOfStars)
        
C     ---   Commons   ---
      common /vel/ uu,vv,ww
      common /coorcil/ coordinate_R,coordinate_Theta,coordinate_Zcylindr
      common /tm/ starBirthTime,m
      common /index/ iwd,ntwd         
      common /plano/ xpla,ypla
      common /carte/ xcar,ycar

C     ---   Externals  ---
C     QUESTION: what is this?
      EXTERNAL DERIVS
      EXTERNAL RKQC
            
C     ---   Test   ---
      njumps=100
      hmin=0.0
      eps=1.0d-4
      n=2
      fcgys=(1.0d+9)*365.25*24.0*3600.0
      tf=galacticDiskAge*fcgys 

C     ---   Integrating trajectories   ---
      do 1 i=1,ntwd
        xcar=xpla(i)
        ycar=ypla(i)
        wo=ww(i)+8.0
        zDistribution_zo=coordinate_Zcylindr(i)*(3.086d+16)
        ecini=0.5*wo*wo
        call epot(zDistribution_zo,epoti)
        tb=starBirthTime(i)
        tinc=(galacticDiskAge-tb)/dfloat(njumps)
C       ---  The time in seconds  ---
        htry=tinc*fcgys
C       ---  Initial conditions  ---
        y(1)=zDistribution_zo
        y(2)=wo
        dydx(1)=wo
        call fuerza(zDistribution_zo,f)
        dydx(2)=f
        ti=tb*fcgys 
C       ---  Calling to the Runge-Kutta integrator ---

        call ODEINT(y,n,ti,tf,eps,htry,hmin,NOK,NBAD,DERIVS,RKQC,yscal,
     &       y,dydx)      
        ecinf=0.5*y(2)*y(2)
        call epot(y(1),epotf)     
        coordinate_Zcylindr(i)=y(1)/(3.086d+16)
        ww(i)=y(2)+wosun
 1    continue

      return
      end     
C***********************************************************************
     


C***********************************************************************      
C     TODO: rewrite      
      subroutine epot(z,e)
C=======================================================================
C
C     This function calculates the force along z-coordinate. 
C--------------------------------------------------------------------
C     Input parameters:
C       z: coordinate z (km)
C--------------------------------------------------------------------
C     Output parameters:
C       e: potential energy (km2/s2)
C=======================================================================
      implicit double precision (a-h,m,o-z)    
      
      double precision ro,vh,rc1,mc1,rc2,mc2,b,md1,md2,md3,a1,a2,a3,g
      double precision xcar,ycar,xpla,ypla,rpla,zsig,z,vh2,ro2
      double precision rc12,rc22,b2,rpla2,r2,poth,xa,xb,potc,xx
      double precision xd1,xd2,xd3,potd,potd1,potd2,potd3,pot,e

C     ---   Parameters   ---
      parameter(ro=8.5,vh=220.0)
      parameter(rc1=2.7,mc1=3.0d+09)
      parameter(rc2=0.42,mc2=1.6d+10)
      parameter(b=0.3)
      parameter(md1=6.6d+10,a1=5.81)
      parameter(md2=-2.9d+10,a2=17.43)
      parameter(md3=3.3d+09,a3=34.86)
      parameter(g=4.30026d-6)
      
C     ---   Common  
      common /carte/ xcar,ycar
                  
C     ---  Some calculations of interest   ---
      xpla=xcar
      ypla=ycar
      rpla=dsqrt(xpla*xpla+ypla*ypla)
      zsig=z
      z=dabs(z/(3.086d+16))                  
      vh2=0.5*vh*vh
      ro2=ro*ro
      rc12=rc1*rc1
      rc22=rc2*rc2
      b2=b*b
      rpla2=rpla*rpla
      r2=rpla2+z*z

C     ---   Calculating the potentials  ---           
C     ---   Dark halo  ---
      poth=vh2*dlog(r2+ro2)
C     ---   Central component  ---      
      xa=dsqrt(r2+rc12)
      xb=dsqrt(r2+rc22)
      potc=((-g*mc1)/xa)+((-g*mc2)/xb)
C     ---   Disk  ---
      xx=dsqrt(z*z+b2)
      xd1=rpla2+((a1+xx)*(a1+xx))
      xd2=rpla2+((a2+xx)*(a2+xx))
      xd3=rpla2+((a3+xx)*(a3+xx))
      potd1=(g*md1)/(dsqrt(xd1))
      potd2=(g*md2)/(dsqrt(xd2))
      potd3=(g*md3)/(dsqrt(xd3))
      potd=-potd1-potd2-potd3
C     ---   Total potential  ---
      pot=poth+potc+potd       
      e=pot
      z=zsig
      
      return
      end
C***********************************************************************



C     TODO: rewrite
C***********************************************************************
      subroutine fuerza(z,f)
C=======================================================================
C
C     This function calculates the force along z-coordinate.
C-------------------------------------------------------------------
C     Input parameters: z(km) xpla,ypla (kpc)
C-------------------------------------------------------------------
C     Output parameters: f
C=======================================================================
      implicit double precision (a-h,m,o-z)
      
      double precision ro,vh,rc1,mc1,rc2,mc2,b,md1,md2,md3,a1,a2,a3,g
      double precision xcar,ycar,xpla,ypla,rpla,zsig,z,vh2,ro2
      double precision rc12,rc22,b2,rpla2,r2,foh,foc1,foc2,foc,bzr
      double precision fod1,fod2,fod3,fod,ftot,fcv,f
      
C     ---   Parameters   ---
      parameter(ro=8.5,vh=220.0)
      parameter(rc1=2.7,mc1=3.0d+09)
      parameter(rc2=0.42,mc2=1.6d+10)
      parameter(b=0.3)
      parameter(md1=6.6d+10,a1=5.81)
      parameter(md2=-2.9d+10,a2=17.43)
      parameter(md3=3.3d+09,a3=34.86)
      parameter(g=4.30026d-6)
      
C     ---  Common  --- 
      common /carte/ xcar,ycar
                 
C     ---  Calculating some useful variables ---
      xpla=xcar
      ypla=ycar
      zsig=z
      z=dabs(z)/(3.086d+16)
      rpla2=xpla*xpla+ypla*ypla
      rpla=dsqrt(rpla2)
      vh2=vh*vh
      ro2=ro*ro
      rc12=rc1*rc1
      rc22=rc2*rc2
      b2=b*b
      r2=rpla*rpla+z*z
                   
C     ---   Calculating the forces  ---
C     ---   Dark halo  ---
      foh=vh2*z/(ro2+r2)
C     ---   Central component  ---      
      foc1=g*mc1*z/((rc12+r2)**(1.5d0))   
      foc2=g*mc2*z/((rc22+r2)**(1.5d0))      
      foc=foc1+foc2         
C     ---   Disk  ---
      bzr=dsqrt(b2+z*z)
      fod1=g*md1*z*(a1+bzr)/(bzr*(rpla2+(a1+bzr)*(a1+bzr))**(1.5d0))    
      fod2=g*md2*z*(a2+bzr)/(bzr*(rpla2+(a2+bzr)*(a2+bzr))**(1.5d0))     
      fod3=g*md3*z*(a3+bzr)/(bzr*(rpla2+(a3+bzr)*(a3+bzr))**(1.5d0))    
      fod=fod1+fod2+fod3
C     ---  Total force  ---
      ftot=foh+foc+fod 
C     ---  If we want the result in km/s²  ---
      fcv=1.0d0/(3.086d+16)
      ftot=fcv*dabs(ftot)
C     ---  The sign of z will be  ---
      f=-dsign(ftot,zsig)
      z=zsig
      
      return
      end
C***********************************************************************
      


      
C***********************************************************************
C     TODO: rewrite      
      subroutine coor(solarGalactocentricDistance)
C=======================================================================
C     This subroutine calculates the coordinates and proper motions
C     in the Galactic coordinate system.
C     Also it calculates the right ascension and the declination.
C     galactic coordinates calculated in radians
C     Distances in Kpc
C     Velocities in Km/s
C     Proper motions in  arcsec/yr
C
C     Revised in 22.09.07 by S. Torres
C-------------------------------------------------------------------
C     Input parameters:
C       solarGalactocentricDistance
C       ntwd: number of WDs
C=======================================================================
      implicit double precision (a-h,m,o-z)
      
C     ---   Definition of variables  ---
      integer numberOfStars,i,ntwd
      double precision solarGalactocentricDistance,pi,fi,zsdg,zcdg,ros,
     &                 zz,zzx
      double precision k,alfag,deltag,theta
      double precision zsl,zcl,zsb,zcb,zkr,zkri
      double precision sinpsi,cospsi,xc,xs
      
C     ---   Parameters   ---
      parameter (numberOfStars=6000000)
      parameter (k=4.74d+3)
      parameter (alfag=3.35,deltag=0.478,theta=2.147)
      
C     ---   Dimensions   ---
      double precision coordinate_R(numberOfStars),
     &                 coordinate_Theta(numberOfStars),
     &                 coordinate_Zcylindr(numberOfStars)
      double precision lgac(numberOfStars),bgac(numberOfStars),
     &                 rgac(numberOfStars)
      double precision mpb(numberOfStars),mpl(numberOfStars),
     &                 vr(numberOfStars)
      double precision uu(numberOfStars),vv(numberOfStars),
     &                 ww(numberOfStars)
      double precision mu(numberOfStars)
      double precision arec(numberOfStars),dec(numberOfStars)
      double precision iwd(numberOfStars)
       
C     ---   Commons   ---
      common /coorcil/ coordinate_R,coordinate_Theta,coordinate_Zcylindr
      common /vel/ uu,vv,ww
      common /paral/ rgac
      common /mad/ mu,arec,dec
      common /index/ iwd,ntwd
      common /mopro/ mpb,mpl,vr
      common /lb/ lgac,bgac
      
C     ---   Calculating some parameters  ---      
      pi=4.0*atan(1.0d0)
      fi=180.0/pi
      zsdg=dsin(deltag)
      zcdg=dcos(deltag)

C--------------------------------------------------------------------
C     ---   Calculating galactocentric coordinates (r,l,b)  ---
C--------------------------------------------------------------------
      i=1
      do 1 i=1,ntwd
C        ---   Galactic coordinate r (Kpc) ---
        ros=((solarGalactocentricDistance*solarGalactocentricDistance)+
     &      (coordinate_R(i)*coordinate_R(i))-(2*coordinate_R(i)*
     &      solarGalactocentricDistance*dcos(coordinate_Theta(i))))
        rgac(i)=dsqrt(ros+(coordinate_Zcylindr(i)*
     &          coordinate_Zcylindr(i)))
C     ---   Galactic coordinate lgac ---
        ros=dsqrt(ros)
        zz=(coordinate_R(i)/ros)*dsin(coordinate_Theta(i))
        lgac(i)=dasin(zz)
        if (coordinate_R(i)*dcos(coordinate_Theta(i)).gt.
     &     solarGalactocentricDistance) then 
          lgac(i)=pi-lgac(i)
        else
          if (dsin(coordinate_Theta(i)).lt.0.0) then
            lgac(i)=2.0*pi+lgac(i)
          endif
          continue
        endif
C        ---   galactic coordinate bgac ---
        zzx=dabs(coordinate_Zcylindr(i)/ros)
        bgac(i)=datan(zzx)
        if (coordinate_Zcylindr(i).lt.0.0) then
          bgac(i)=-bgac(i)
        else
         continue
        endif
C--------------------------------------------------------------------        
C       ---   Calculating the proper motions in galactic coordinates ---
C--------------------------------------------------------------------      
C       ---   Calculating some values that are going to be used later---
        zsl=dsin(lgac(i))
        zcl=dcos(lgac(i))
        zsb=dsin(bgac(i))
        zcb=dcos(bgac(i))
        zkr=k*rgac(i)
        zkri=(1.0/zkr)
C       ---   Calculating the components of the proper motion ---
        mpl(i)=(-zkri*(zsl/zcb)*uu(i))+(zkri*(zcl/zcb)*vv(i))
        mpb(i)=(-zkri*zcl*zsb*uu(i))+(-zkri*zsb*zsl*vv(i))+
     &         (zkri*zcb*ww(i))
        vr(i)=(zcb*zcl*uu(i))+(zcb*zsl*vv(i))+(zsb*ww(i)) 
        mu(i)=dsqrt(mpl(i)*mpl(i)+mpb(i)*mpb(i))
C-------------------------------------------------------------------
C       ---   Calculating right ascension and the declination  ---
C-------------------------------------------------------------------      
C       ---   Calculating the declination   ---
        zz=zsdg*zsb+zcdg*zcb*dcos(theta-lgac(i))
        dec(i)=dasin(zz)
C       ---   Calculating the right ascension ---
        xs= ((zcb*dsin(theta-lgac(i)))/dcos(dec(i)))
        xc= ((zcdg*zsb-zsdg*zcb*dcos(theta-lgac(i)))/dcos(dec(i)))
C       --Looking at the sign that corresponds to the right ascension---
        if (xs.ge.0.0) then
          if (xc.ge.0.0) then 
            arec(i)=dasin(xs)+alfag
          else
            arec(i)=dacos(xc)+alfag
          endif    
        else
          if (xc.lt.0.0) then
            arec(i)=pi-dasin(xs)+alfag
          else
            arec(i)=2*pi+dasin(xs)+alfag
          endif
        endif            
        if (arec(i)*fi.gt.360.0) then
          arec(i)=arec(i)-2*pi
        endif
C--------------------------------------------------------------------
C       ---  Calculating the proper motion in ecuatorial coordinates ---
C--------------------------------------------------------------------
        sinpsi=dsin(theta)
        cospsi=dcos(theta)    
C-------------------------------------------------------------------
C       ---   Calculating the proper motion  (arc sec/yr)   ---
C-------------------------------------------------------------------
  1   continue
           
      return
      end
C***********************************************************************


C***********************************************************************
C     TODO: rewrite      
      subroutine magi(fractionOfDB)
C=======================================================================
C
C     This subroutine calculates ltc,cbv,cvi,cvr,cuv visual absolute 
C     and apparent magnitude of the WDs.
C     
C     Created by S.Torres
C     Introduced metallicity 08/2012 (ER Cojocaru)
C
C-----------------------------------------------------------------------
C     Input parameters:
C       numberOfStarsInSample
C-----------------------------------------------------------------------
C     Output parameters
C       none
C=======================================================================
      implicit double precision (a-h,m,o-z)

      integer numberOfStars,iseed,i,ntwd,in
      
C     ---   Variables  ---
      double precision lum,teff,xlog,c1,c2,c3,c4,c5,n1,n2,n3,n4,n5
      double precision UB,BV,VR,RI,xg,xug,xgr,xri,xiz,xgi,fractionOfDB
      double precision mone

C     ---   Parameters  ---
      parameter (numberOfStars=6000000)
      parameter (mone=1.14)

C     ---   Dimensions  ---
      double precision leb(numberOfStars),meb(numberOfStars),
     &                 zeb(numberOfStars),teb(numberOfStars)
      double precision iwd(numberOfStars)
      double precision rgac(numberOfStars)
      double precision g(numberOfStars),go(numberOfStars),
     &                 gr(numberOfStars),v(numberOfStars)
      double precision gi(numberOfStars),ur(numberOfStars),
     &                 rz(numberOfStars)
      double precision tcool(numberOfStars)
      double precision idb(numberOfStars)

C     ---   Commons   ---
      common /enanas/ leb,meb,zeb,teb
      common /index/ iwd,ntwd
      common /paral/ rgac
      common /photo/ go,gr,gi,ur,rz
      common /johnson/ v
      common /cool/ tcool
      common /indexdb/ idb

      n1=0
      n2=0
      n3=0
      n4=0
      n5=0

C     ---  Interpolating Mv, luminosity, colors and other variables
C          from tcool and the mwd  ---
C     ---  Start DO on all the stars
      do i=1,ntwd
C       ---  ATENTION! choosing only if .lt.1.1!!!  ---
C       ---  Start IF mass <1.4  ----
        if(meb(i).le.1.4) then
C         ---  IF CO core ---
          if(meb(i).lt.mone) then  
C           --- Atention We put the "old" ones to 0.6Msol ---
C           --- Distribucion DA/DB ---
            call dbd_fid(iseed,fractionOfDB,in)
C           --- End of distribution ---
C           ---  IF DA ---
            if(in.eq.0) then
              idb(i)=0
              n1=n1+1
              call interlumda(tcool(i),meb(i),zeb(i),lum,teff,xlog,c1,
     &             c2,c3,c4,c5)
C           ---  ELSE DB  ---
            else
              n3=n3+1
              idb(i)=1    
              call interlumdb(tcool(i),meb(i),zeb(i),lum,c1,c2,c3,c4,c5,
     &             teff,xlog)
              if(teff.lt.6000) n5=n5+1
            end if
C           ---  END IF DB/NON-DB
C         ---  ELSE ONe ---
          else
            n2=n2+1
            idb(i)=2
            call interlumone(tcool(i),meb(i),lum,c1,c2,c3,c4,c5,teff,
     &           xlog)
          end if
C         ---  END IF CO/ONe ---
          if(teff.lt.6000) n4=n4+1
          leb(i)=-lum
          teb(i)=teff            
          V(i)=c3
          UB=c1-c2
          BV=c2-c3
          VR=c3-c4
          RI=c4-c5
          call chanco(V(i),UB,BV,VR,RI,xg,xug,xgr,xri,xiz,xgi)
          g(i)=xg
          gr(i)=xgr
          gi(i)=xgi
          ur(i)=xug+xgr
          rz(i)=xri+xiz
C         ---  Making g and V apparent magnitude ---
          go(i)=g(i)-5.0d0+5.0d0*(dlog10(rgac(i))+3.0d0)
          V(i)=V(i)-5.0d0+5.0d0*(dlog10(rgac(i))+3.0d0)
C       ---  ELSE mass >= 1.4  --- EXPLOTA, exceeding Chandrasekar limit
        else
          idb(i)=5
        end if
C       ---  END IF about WD mass ---
      end do
C     ---  END DO about all the stars ---

      write(*,*) "DA CO ",n1
      write(*,*) "DA ONe ",n2
      write(*,*) "DB",n3
      write(*,*) "<6000 DA", n4
      write(*,*) "<6000 DB", n5

      return
      end
C***********************************************************************




C***********************************************************************
C     TODO: rewrite      
      subroutine interlumda(tcool,mass,Z,lum,teff,logg,c1,c2,c3,c4,c5)
C=======================================================================
C
C     This subroutine interpolates the luminosity of DA WD according
C     to its mass, cooling time and metallicity.
C     It uses the cooling tables by Althaus et al. (2009) and
C     Renedo et al. (2010)
C
C     Modifications in 07/2012 (ER Cojocaru)
C
C-------------------------------------------------------------------
C     Input parameters:
C       tcool: cooling time
C       mass: mass of WD
C-------------------------------------------------------------------
C     Output parameters
C       lum: luminosity
C       teff: effective temperature [K]
C       logg: logarithm of the superficial gravity
C       c1,c2,c3,c4,c5: UBVRI colors 
C=======================================================================
      implicit double precision (a-h,m,o-z)
      
C     ---   Parameters   ---
      integer nrow,model,modlog
      parameter (nrow=650)
      double precision zet1,zet2,zet3,zet4
      parameter (zet1=0.001)
      parameter (zet2=0.01)
      parameter (zet3=0.03)
      parameter (zet4=0.06)

C     ---   Declaration of variables   ---
      integer numberOfMassesWithColors,numberOfMassesWithCoolSeq_1,
     &        numberOfMassesWithCoolSeq_2,numberOfMassesWithCoolSeq_3,
     &        numberOfMassesWithCoolSeq_4
      integer ntrkda(10),numberOfRows_1(7),numberOfRows_2(10),
     &        numberOfRows_3(8),numberOfRows_4(8)
      double precision tcool,mass,Z,lum,teff,logg,c1,c2,c3,c4,c5
      double precision Z1,lum1,teff1,logg1
      double precision Z2,lum2,teff2,logg2
      double precision massOfWD(10),massOfWD_1(7),massOfWD_2(10),
     &                 massOfWD_3(8),massOfWD_4(8)
      double precision tprewdda1(7),tprewdda2(10),tprewdda3(8)
      double precision tprewdda4(8),coolingTimes_1(7,nrow),
     &                 coolingTimes_2(10,nrow)
      double precision coolingTimes_3(8,nrow),coolingTimes_4(8,nrow)
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
      double precision luminosity(10,nrow),color_U(10,nrow),
     &                 color_B(10,nrow),color_V(10,nrow)
      double precision color_R(10,nrow),color_I(10,nrow)

C     ---   Commons   ---
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
      common /colors/ luminosity,color_U,color_B,color_V,color_R,color_I
      
      model=0
      modlog=0
C     QUESTION: what does it mean?
C     miro entre que valores de metalicidad se encuentra el Z de entrada
C     calculo lum,teff,logg para esas valores Z1 y Z2 y después hago la 
C     Interpolación para Z
      if(Z.ge.zet1.AND.Z.lt.zet2) then
        Z1=zet1
        call interp(model,modlog,tcool,mass,numberOfMassesWithCoolSeq_1,
     &       numberOfRows_1,coolingTimes_1,tprewdda1,massOfWD_1,
     &       luminosity_1,lum1)
        modlog=1
        call interp(model,modlog,tcool,mass,numberOfMassesWithCoolSeq_1,
     &       numberOfRows_1,coolingTimes_1,tprewdda1,massOfWD_1,
     &       effectiveTemperature_1,teff1)
        modlog=0
        call interp(model,modlog,tcool,mass,numberOfMassesWithCoolSeq_1,
     &       numberOfRows_1,coolingTimes_1,tprewdda1,massOfWD_1,
     &       gravitationalAcceleration_1,logg1)
        Z2=zet2
        call interp(model,modlog,tcool,mass,
     &       numberOfMassesWithCoolSeq_2,numberOfRows_2,coolingTimes_2,
     &       tprewdda2,massOfWD_2,luminosity_2,lum2)
        modlog=1
        call interp(model,modlog,tcool,mass,
     &       numberOfMassesWithCoolSeq_2,numberOfRows_2,coolingTimes_2,
     &       tprewdda2,massOfWD_2,effectiveTemperature_2,teff2)
        modlog=0
        call interp(model,modlog,tcool,mass,
     &       numberOfMassesWithCoolSeq_2,numberOfRows_2,coolingTimes_2,
     &       tprewdda2,massOfWD_2,gravitationalAcceleration_2,logg2)
      end if
      
      if(Z.ge.zet2.AND.Z.lt.zet3) then
        Z1=zet2
        call interp(model,modlog,tcool,mass,
     &       numberOfMassesWithCoolSeq_2,numberOfRows_2,coolingTimes_2,
     &       tprewdda2,massOfWD_2,luminosity_2,lum1)
        modlog=1
        call interp(model,modlog,tcool,mass,
     &       numberOfMassesWithCoolSeq_2,numberOfRows_2,coolingTimes_2,
     &       tprewdda2,massOfWD_2,effectiveTemperature_2,teff1)
        modlog=0
        call interp(model,modlog,tcool,mass,
     &       numberOfMassesWithCoolSeq_2,numberOfRows_2,coolingTimes_2,
     &       tprewdda2,massOfWD_2,gravitationalAcceleration_2,logg1)
        Z2=zet3
        call interp(model,modlog,tcool,mass,
     &       numberOfMassesWithCoolSeq_3,numberOfRows_3,coolingTimes_3,
     &       tprewdda3,massOfWD_3,luminosity_3,lum2)
        modlog=1
        call interp(model,modlog,tcool,mass,
     &       numberOfMassesWithCoolSeq_3,numberOfRows_3,coolingTimes_3,
     &       tprewdda3,massOfWD_3,effectiveTemperature_3,teff2)
        modlog=0
        call interp(model,modlog,tcool,mass,
     &       numberOfMassesWithCoolSeq_3,numberOfRows_3,coolingTimes_3,
     &       tprewdda3,massOfWD_3,gravitationalAcceleration_3,logg2)
      end if
      
      if(Z.ge.zet3.AND.Z.le.zet4) then
        Z1=zet3
        call interp(model,modlog,tcool,mass,
     &       numberOfMassesWithCoolSeq_3,numberOfRows_3,coolingTimes_3,
     &       tprewdda3,massOfWD_3,luminosity_3,lum1)
        modlog=1
        call interp(model,modlog,tcool,mass,
     &       numberOfMassesWithCoolSeq_3,numberOfRows_3,coolingTimes_3,
     &       tprewdda3,massOfWD_3,effectiveTemperature_3,teff1)
        modlog=0
        call interp(model,modlog,tcool,mass,
     &       numberOfMassesWithCoolSeq_3,numberOfRows_3,coolingTimes_3,
     &       tprewdda3,massOfWD_3,gravitationalAcceleration_3,logg1)
        Z2=zet4
        call interp(model,modlog,tcool,mass,
     &       numberOfMassesWithCoolSeq_4,numberOfRows_4,coolingTimes_4,
     &       tprewdda4,massOfWD_4,luminosity_4,lum2)
        modlog=1
        call interp(model,modlog,tcool,mass,
     &       numberOfMassesWithCoolSeq_4,numberOfRows_4,coolingTimes_4,
     &       tprewdda4,massOfWD_4,effectiveTemperature_4,teff2)
        modlog=0
        call interp(model,modlog,tcool,mass,
     &       numberOfMassesWithCoolSeq_4,numberOfRows_4,coolingTimes_4,
     &       tprewdda4,massOfWD_4,gravitationalAcceleration_4,logg2)
      end if
      
      lum=lum1+(lum2-lum1)*(Z-Z1)/(Z2-Z1)
      teff=teff1+(teff2-teff1)*(Z-Z1)/(Z2-Z1)
      logg=logg1+(logg2-logg1)*(Z-Z1)/(Z2-Z1)
      
C     interpolation of colors, here there is no metalicity
      call intermag(mass,lum,numberOfMassesWithColors,ntrkda,massOfWD,
     &     luminosity,color_U,color_B,color_V,color_R,color_I,c1,c2,c3,
     &     c4,c5)

      return
      end
C***********************************************************************





C***********************************************************************
      subroutine dbd_fid(iseed,fractionOfDB,in)
C=======================================================================
C
C     This subroutine determines if the WD is DA or non-DA (DB),
C     using the model A (20/80) by Torres (A&A, 2010) 
C    
C     in=0 DA
C     in=1 non-DA
C
C     Created by ER Cojocaru (11/2012)
C
C=======================================================================
      implicit double precision (a-h,m,o-z)
      
C     ---   Declaration of variables  ---
      external ran
      real ran
      integer iseed,in
      double precision x,fractionOfDB
      
C     --- Fiducial model: 20% DB, 80% DA
      x=ran(iseed)
      if(x.lt.fractionOfDB) then 
        in=1
      else
        in=0
      endif

      return
      end
C***********************************************************************





C***********************************************************************
C     TODO: rewrite
      subroutine intermag(mass,lumi,numberOfMassesWithColors,ntrk,mtrk,
     &           luminosity,color_U,color_B,color_V,color_R,color_I,c1,
     &           c2,c3,c4,c5)
C=======================================================================
C
C     This subroutine interpolates the luminosity of a DA WD
C     according to its mass and cooling time using the cooling sequence 
C     from input (corresponding to certain metallicity)
C
C     Modifications in 07/2012 (ER Cojocaru)
C
C-------------------------------------------------------------------
C     Input parameters:
C       mass: WD mass
C       lumi: luminosity
C       ltrx,color_U,color_B,color_V,color_R,color_I,ntrk,
C         numberOfMassesWithColors: information 
C           seq. colors (DA or DB)     
C
C-------------------------------------------------------------------
C     Output parameters:blanca
C       c1,c2,c3,c4,c5: Johnson colors
C
C=======================================================================
      implicit double precision (a-h,m,o-z)

C     ---   Declaraеion of variables   ---
      integer numberOfMassesWithColors,i,k,check1,check2,check3
      integer i1,i2,ntrk(numberOfMassesWithColors),ns1,ns2
      double precision mass,lumi,c1,c2,c3,c4,c5
      double precision c_1,c_2,a1,a2,b1,b2
      double precision mtrk(numberOfMassesWithColors)
      double precision luminosity(numberOfMassesWithColors,*),
     &                 color_U(numberOfMassesWithColors,*),
     &                 color_B(numberOfMassesWithColors,*)
      double precision color_V(numberOfMassesWithColors,*),
     &                 color_R(numberOfMassesWithColors,*),
     &                 color_I(numberOfMassesWithColors,*)

      check1=0
      check2=0
      check3=0

C     TODO: check if i mixed interpolation with extrapolation further
C     Smaller mass than known --> linear 2D extrapolation 
C     (using luminosity and mass)

      if(mass.le.mtrk(1)) then
        ns1=ntrk(1)
        ns2=ntrk(2)
C       QUESTION: what does it mean?
C       lum mayor que las conocidas --> linear 2D extrapolation
        if(lumi.gt.luminosity(1,1).OR.lumi.gt.luminosity(2,1)) then
          call extrap1(lumi,color_U(1,1),color_U(1,2),luminosity(1,1),
     &         luminosity(1,2),c_1)
          call extrap1(lumi,color_U(2,1),color_U(2,2),luminosity(2,1),
     &         luminosity(2,2),c_2)
          call extrap1(mass,c_1,c_2,mtrk(1),mtrk(2),c1)
          if(c1.lt.0.0) c1=0.0
          call extrap1(lumi,color_B(1,1),color_B(1,2),luminosity(1,1),
     &         luminosity(1,2),c_1)
          call extrap1(lumi,color_B(2,1),color_B(2,2),luminosity(2,1),
     &         luminosity(2,2),c_2)
          call extrap1(mass,c_1,c_2,mtrk(1),mtrk(2),c2)
          if(c2.lt.0.0) c2=0.0
          call extrap1(lumi,color_V(1,1),color_V(1,2),luminosity(1,1),
     &         luminosity(1,2),c_1)
          call extrap1(lumi,color_V(2,1),color_V(2,2),luminosity(2,1),
     &         luminosity(2,2),c_2)
          call extrap1(mass,c_1,c_2,mtrk(1),mtrk(2),c3)
          if(c3.lt.0.0) c3=0.0
          call extrap1(lumi,color_R(1,1),color_R(1,2),luminosity(1,1),
     &         luminosity(1,2),c_1)
          call extrap1(lumi,color_R(2,1),color_R(2,2),luminosity(2,1),
     &         luminosity(2,2),c_2)
          call extrap1(mass,c_1,c_2,mtrk(1),mtrk(2),c4)
          if(c4.lt.0.0) c4=0.0
          call extrap1(lumi,color_I(1,1),color_I(1,2),luminosity(1,1),
     &         luminosity(1,2),c_1)
          call extrap1(lumi,color_I(2,1),color_I(2,2),luminosity(2,1),
     &         luminosity(2,2),c_2)
          call extrap1(mass,c_1,c_2,mtrk(1),mtrk(2),c5)
          if(c5.lt.0.0) c5=0.0
          check3=1
          GOTO 45
C       lum menor que las conocidas --> linear 2D extrapolation
        elseif (lumi.lt.luminosity(1,ns1).OR.lumi.lt.luminosity(2,ns2)) 
     &  then
          call extrap1(lumi,color_U(1,ns1-1),color_U(1,ns1),
     &         luminosity(1,ns1-1),luminosity(1,ns1),c_1)
          call extrap1(lumi,color_U(2,ns2-1),color_U(2,ns2),
     &         luminosity(2,ns2-1),luminosity(2,ns2),c_2)
          call extrap1(mass,c_1,c_2,mtrk(1),mtrk(2),c1)
          if(c1.lt.0.0) c1=0.0
          call extrap1(lumi,color_B(1,ns1-1),color_B(1,ns1),
     &         luminosity(1,ns1-1),luminosity(1,ns1),c_1)
          call extrap1(lumi,color_B(2,ns2-1),color_B(2,ns2),
     &         luminosity(2,ns2-1),luminosity(2,ns2),c_2)
          call extrap1(mass,c_1,c_2,mtrk(1),mtrk(2),c2)
          if(c2.lt.0.0) c2=0.0
          call extrap1(lumi,color_V(1,ns1-1),color_V(1,ns1),
     &         luminosity(1,ns1-1),luminosity(1,ns1),c_1)
          call extrap1(lumi,color_V(2,ns2-1),color_V(2,ns2),
     &         luminosity(2,ns2-1),luminosity(2,ns2),c_2)
          call extrap1(mass,c_1,c_2,mtrk(1),mtrk(2),c3)
          if(c3.lt.0.0) c3=0.0
          call extrap1(lumi,color_R(1,ns1-1),color_R(1,ns1),
     &         luminosity(1,ns1-1),luminosity(1,ns1),c_1)
          call extrap1(lumi,color_R(2,ns2-1),color_R(2,ns2),
     &         luminosity(2,ns2-1),luminosity(2,ns2),c_2)
          call extrap1(mass,c_1,c_2,mtrk(1),mtrk(2),c4)
          if(c4.lt.0.0) c4=0.0      
          call extrap1(lumi,color_I(1,ns1-1),color_I(1,ns1),
     &         luminosity(1,ns1-1),luminosity(1,ns1),c_1)
          call extrap1(lumi,color_I(2,ns2-1),color_I(2,ns2),
     &         luminosity(2,ns2-1),luminosity(2,ns2),c_2)
          call extrap1(mass,c_1,c_2,mtrk(1),mtrk(2),c5)
          if(c5.lt.0.0) c5=0.0
          check3=1
          GOTO 45
C       lum contenida entre las conocidas
        else
          do i=1,ns1-1
            if(lumi.ge.luminosity(1,i+1).AND.lumi.le.luminosity(1,i)) 
     &      then
              i1=i
              check1=1
              GOTO 5
            end if
          end do
5         continue
          do i=1,ns2-1
            if(lumi.ge.luminosity(2,i+1).AND.lumi.le.luminosity(2,i)) 
     &      then
              i2=i
              check2=1
              GOTO 10
            end if
          end do
10        continue
          if(check1.eq.1.AND.check2.eq.1) then
            check3=1
            a1=lumi-luminosity(1,i1)
            a2=lumi-luminosity(2,i2)
            b1=luminosity(1,i1+1)-luminosity(1,i1)
            b2=luminosity(2,i2+1)-luminosity(2,i2)
            c_1=color_U(1,i1)+(color_U(1,i1+1)-color_U(1,i1))*a1/b1
            c_2=color_U(2,i2)+(color_U(2,i2+1)-color_U(2,i2))*a2/b2
            call extrap1(mass,c_1,c_2,mtrk(1),mtrk(2),c1)
            c_1=color_B(1,i1)+(color_B(1,i1+1)-color_B(1,i1))*a1/b1
            c_2=color_B(2,i2)+(color_B(2,i2+1)-color_B(2,i2))*a2/b2
            call extrap1(mass,c_1,c_2,mtrk(1),mtrk(2),c2)
            c_1=color_V(1,i1)+(color_V(1,i1+1)-color_V(1,i1))*a1/b1
            c_2=color_V(2,i2)+(color_V(2,i2+1)-color_V(2,i2))*a2/b2
            call extrap1(mass,c_1,c_2,mtrk(1),mtrk(2),c3)
            c_1=color_R(1,i1)+(color_R(1,i1+1)-color_R(1,i1))*a1/b1
            c_2=color_R(2,i2)+(color_R(2,i2+1)-color_R(2,i2))*a2/b2
            call extrap1(mass,c_1,c_2,mtrk(1),mtrk(2),c4)
            c_1=color_I(1,i1)+(color_I(1,i1+1)-color_I(1,i1))*a1/b1
            c_2=color_I(2,i2)+(color_I(2,i2+1)-color_I(2,i2))*a2/b2
            call extrap1(mass,c_1,c_2,mtrk(1),mtrk(2),c5)
            GOTO 45
          end if
        end if
      end if
      
C     masa mayor que las conocidas --> linear interpolation
      if(mass.gt.mtrk(numberOfMassesWithColors)) then
        ns1=ntrk(numberOfMassesWithColors-1)
        ns2=ntrk(numberOfMassesWithColors)
C       lum mayor que las conocidas --> linear 2D extrapolation
        if(lumi.gt.luminosity(numberOfMassesWithColors-1,1).OR.lumi.gt.
     &  luminosity(numberOfMassesWithColors,1)) then
          call extrap1(lumi,color_U(numberOfMassesWithColors-1,1),
     &         color_U(numberOfMassesWithColors-1,2),
     &         luminosity(numberOfMassesWithColors-1,1),
     &         luminosity(numberOfMassesWithColors-1,2),c_1)
          call extrap1(lumi,color_U(numberOfMassesWithColors,1),
     &         color_U(numberOfMassesWithColors,2),
     &         luminosity(numberOfMassesWithColors,1),
     &         luminosity(numberOfMassesWithColors,2),c_2)
          call extrap1(mass,c_1,c_2,mtrk(numberOfMassesWithColors-1),
     &         mtrk(numberOfMassesWithColors),c1)
          if(c1.lt.0.0) c1=0.0    
          call extrap1(lumi,color_B(numberOfMassesWithColors-1,1),
     &         color_B(numberOfMassesWithColors-1,2),
     &         luminosity(numberOfMassesWithColors-1,1),
     &         luminosity(numberOfMassesWithColors-1,2),c_1)
          call extrap1(lumi,color_B(numberOfMassesWithColors,1),
     &         color_B(numberOfMassesWithColors,2),
     &         luminosity(numberOfMassesWithColors,1),
     &         luminosity(numberOfMassesWithColors,2),c_2)
          call extrap1(mass,c_1,c_2,mtrk(numberOfMassesWithColors-1),
     &         mtrk(numberOfMassesWithColors),c2)
          if(c2.lt.0.0) c2=0.0
          call extrap1(lumi,color_V(numberOfMassesWithColors-1,1),
     &         color_V(numberOfMassesWithColors-1,2),
     &         luminosity(numberOfMassesWithColors-1,1),
     &         luminosity(numberOfMassesWithColors-1,2),c_1)
          call extrap1(lumi,color_V(numberOfMassesWithColors,1),
     &         color_V(numberOfMassesWithColors,2),
     &         luminosity(numberOfMassesWithColors,1),
     &         luminosity(numberOfMassesWithColors,2),c_2)
          call extrap1(mass,c_1,c_2,mtrk(numberOfMassesWithColors-1),
     &         mtrk(numberOfMassesWithColors),c3)
          if(c3.lt.0.0) c3=0.0
          call extrap1(lumi,color_R(numberOfMassesWithColors-1,1),
     &         color_R(numberOfMassesWithColors-1,2),
     &         luminosity(numberOfMassesWithColors-1,1),
     &         luminosity(numberOfMassesWithColors-1,2),c_1)
          call extrap1(lumi,color_R(numberOfMassesWithColors,1),
     &         color_R(numberOfMassesWithColors,2),
     &         luminosity(numberOfMassesWithColors,1),
     &         luminosity(numberOfMassesWithColors,2),c_2)
          call extrap1(mass,c_1,c_2,mtrk(numberOfMassesWithColors-1),
     &         mtrk(numberOfMassesWithColors),c4)
          if(c4.lt.0.0) c4=0.0
          call extrap1(lumi,color_I(numberOfMassesWithColors-1,1),
     &         color_I(numberOfMassesWithColors-1,2),
     &         luminosity(numberOfMassesWithColors-1,1),
     &         luminosity(numberOfMassesWithColors-1,2),c_1)
          call extrap1(lumi,color_I(numberOfMassesWithColors,1),
     &         color_I(numberOfMassesWithColors,2),
     &         luminosity(numberOfMassesWithColors,1),
     &         luminosity(numberOfMassesWithColors,2),c_2)
          call extrap1(mass,c_1,c_2,mtrk(numberOfMassesWithColors-1),
     &         mtrk(numberOfMassesWithColors),c5)
          if(c5.lt.0.0) c5=0.0
          check3=1
          GOTO 45
C       lum menor que las conocidas --> linear 2D extrapolation
        elseif (lumi.lt.luminosity(numberOfMassesWithColors-1,ns1).OR.
     &  lumi.lt.luminosity(numberOfMassesWithColors,ns1)) then
          call extrap1(lumi,color_U(numberOfMassesWithColors-1,ns1-1),
     &         color_U(numberOfMassesWithColors-1,ns1),
     &         luminosity(numberOfMassesWithColors-1,ns1-1),
     &         luminosity(numberOfMassesWithColors-1,ns1),c_1)
          call extrap1(lumi,color_U(numberOfMassesWithColors,ns2-1),
     &         color_U(numberOfMassesWithColors,ns2),
     &        luminosity(numberOfMassesWithColors,ns2-1),
     &        luminosity(numberOfMassesWithColors,ns2),c_2)
          call extrap1(mass,c_1,c_2,mtrk(1),mtrk(2),c1)
          call extrap1(lumi,color_B(numberOfMassesWithColors-1,ns1-1),
     &         color_B(numberOfMassesWithColors-1,ns1),
     &         luminosity(numberOfMassesWithColors-1,ns2-1),
     &         luminosity(numberOfMassesWithColors-1,ns1),c_1)
          call extrap1(lumi,color_B(numberOfMassesWithColors,ns2-1),
     &         color_B(numberOfMassesWithColors,ns2),
     &         luminosity(numberOfMassesWithColors,ns2-1),
     &         luminosity(numberOfMassesWithColors,ns2),c_2)
          call extrap1(mass,c_1,c_2,mtrk(1),mtrk(2),c2)
          call extrap1(lumi,color_V(numberOfMassesWithColors-1,ns1-1),
     &         color_V(numberOfMassesWithColors-1,ns1),
     &         luminosity(numberOfMassesWithColors-1,ns1-1),
     &         luminosity(numberOfMassesWithColors-1,ns1),c_1)
          call extrap1(lumi,color_V(numberOfMassesWithColors,ns2-1),
     &         color_V(numberOfMassesWithColors,ns2),
     &         luminosity(numberOfMassesWithColors,ns2-1),
     &         luminosity(numberOfMassesWithColors,ns2),c_2)
          call extrap1(mass,c_1,c_2,mtrk(1),mtrk(2),c3)
          call extrap1(lumi,color_R(numberOfMassesWithColors-1,ns1-1),
     &         color_R(numberOfMassesWithColors-1,ns1),
     &         luminosity(numberOfMassesWithColors-1,ns1-1),
     &         luminosity(numberOfMassesWithColors-1,ns1),c_1)
          call extrap1(lumi,color_R(numberOfMassesWithColors,ns2-1),
     &         color_R(numberOfMassesWithColors,ns2),
     &         luminosity(numberOfMassesWithColors,ns2-1),
     &         luminosity(numberOfMassesWithColors,ns2),c_2)
          call extrap1(mass,c_1,c_2,mtrk(1),mtrk(2),c4)
          call extrap1(lumi,color_I(numberOfMassesWithColors-1,ns1-1),
     &         color_I(numberOfMassesWithColors-1,ns1),
     &         luminosity(numberOfMassesWithColors-1,ns1-1),
     &         luminosity(numberOfMassesWithColors-1,ns1),c_1)
          call extrap1(lumi,color_I(numberOfMassesWithColors,ns2-1),
     &         color_I(numberOfMassesWithColors,ns2),
     &         luminosity(numberOfMassesWithColors,ns2-1),
     &         luminosity(numberOfMassesWithColors,ns2),c_2)
          call extrap1(mass,c_1,c_2,mtrk(1),mtrk(2),c5)      
          check3=1
          GOTO 45
C       lum contenida entre las conocidas
        else
          do i=1,ns1-1
            if(lumi.ge.luminosity(numberOfMassesWithColors-1,i+1).AND.
     &      lumi.le.luminosity(numberOfMassesWithColors-1,i)) then
              i1=i
              check1=1
              GOTO 15
            end if
          end do
15        continue
          do i=1,ns2-1
            if(lumi.ge.luminosity(numberOfMassesWithColors,i+1).AND.
     &      lumi.le.luminosity(numberOfMassesWithColors,i)) then
              i2=i
              check2=1
              GOTO 20
            end if
          end do
20        continue
          if(check1.eq.1.AND.check2.eq.1) then
            check3=1
            a1=lumi-luminosity(numberOfMassesWithColors-1,i1)
            a2=lumi-luminosity(numberOfMassesWithColors,i2)
            b1=luminosity(numberOfMassesWithColors-1,i1+1)-
     &         luminosity(numberOfMassesWithColors-1,i1)
            b2=luminosity(numberOfMassesWithColors,i2+1)-
     &         luminosity(numberOfMassesWithColors,i2)
            c_1=color_U(numberOfMassesWithColors-1,i1)+
     &          (color_U(numberOfMassesWithColors-1,i1+1)-
     &          color_U(numberOfMassesWithColors-1,i1))*a1/b1
            c_2=color_U(numberOfMassesWithColors,i2)+
     &          (color_U(numberOfMassesWithColors,i2+1)-
     &          color_U(numberOfMassesWithColors,i2))*a2/b2
            call extrap1(mass,c_1,c_2,mtrk(1),
     &           mtrk(numberOfMassesWithColors),c1)
            c_1=color_B(numberOfMassesWithColors-1,i1)+
     &          (color_B(numberOfMassesWithColors-1,i1+1)-
     &          color_B(numberOfMassesWithColors-1,i1))*a1/b1
            c_2=color_B(numberOfMassesWithColors,i2)+
     &          (color_B(numberOfMassesWithColors,i2+1)-
     &          color_B(numberOfMassesWithColors,i2))*a2/b2
            call extrap1(mass,c_1,c_2,mtrk(numberOfMassesWithColors-1),
     &           mtrk(numberOfMassesWithColors),c2)
            c_1=color_V(numberOfMassesWithColors-1,i1)+
     &          (color_V(numberOfMassesWithColors-1,i1+1)-
     &          color_V(numberOfMassesWithColors-1,i1))*a1/b1
            c_2=color_V(numberOfMassesWithColors,i2)+
     &          (color_V(numberOfMassesWithColors,i2+1)-
     &          color_V(numberOfMassesWithColors,i2))*a2/b2
            call extrap1(mass,c_1,c_2,mtrk(1),
     &           mtrk(numberOfMassesWithColors),c3)
            c_1=color_R(numberOfMassesWithColors-1,i1)+
     &          (color_R(numberOfMassesWithColors-1,i1+1)-
     &          color_R(numberOfMassesWithColors-1,i1))*a1/b1
            c_2=color_R(numberOfMassesWithColors,i2)+
     &          (color_R(numberOfMassesWithColors,i2+1)-
     &          color_R(numberOfMassesWithColors,i2))*a2/b2
            call extrap1(mass,c_1,c_2,mtrk(numberOfMassesWithColors-1),
     &           mtrk(numberOfMassesWithColors),c4)
            c_1=color_I(numberOfMassesWithColors-1,i1)+
     &          (color_I(numberOfMassesWithColors-1,i1+1)-
     &          color_I(numberOfMassesWithColors-1,i1))*a1/b1
            c_2=color_I(numberOfMassesWithColors,i2)+
     &          (color_I(numberOfMassesWithColors,i2+1)-
     &          color_I(numberOfMassesWithColors,i2))*a2/b2
            call extrap1(mass,c_1,c_2,mtrk(numberOfMassesWithColors-1),
     &           mtrk(numberOfMassesWithColors),c5)
            GOTO 45
          end if
        end if
      end if

C     masa contenida entre las conocidas --> linear interpolation
      do k=1,numberOfMassesWithColors-1
        if(mass.gt.mtrk(k).AND.mass.le.mtrk(k+1)) then
          ns1=ntrk(k)
          ns2=ntrk(k+1)
C         lum mayor que las conocidas --> linear 2D extrapolation
          if(lumi.gt.luminosity(k,1).OR.lumi.gt.luminosity(k+1,1)) then      
            call extrap1(lumi,color_U(k,1),color_U(k,2),luminosity(k,1),
     &           luminosity(k,2),c_1)
            call extrap1(lumi,color_U(k+1,1),color_U(k+1,2),
     &           luminosity(k+1,1),luminosity(k+1,2),c_2)
            call extrap1(mass,c_1,c_2,mtrk(k),mtrk(k+1),c1)
            if(c1.lt.0.0) c1=0.0      
            call extrap1(lumi,color_B(k,1),color_B(k,2),luminosity(k,1),
     &           luminosity(k,2),c_1)
            call extrap1(lumi,color_B(k+1,1),color_B(k+1,2),
     &           luminosity(k+1,1),luminosity(k+1,2),c_2)
            call extrap1(mass,c_1,c_2,mtrk(k),mtrk(k+1),c2)
            if(c2.lt.0.0) c2=0.0
            call extrap1(lumi,color_V(k,1),color_V(k,2),luminosity(k,1),
     &           luminosity(k,2),c_1)
            call extrap1(lumi,color_V(k+1,1),color_V(k+1,2),
     &           luminosity(k+1,1),luminosity(k+1,2),c_2)
            call extrap1(mass,c_1,c_2,mtrk(k),mtrk(k+1),c3)
            if(c3.lt.0.0) c3=0.0
            call extrap1(lumi,color_R(k,1),color_R(k,2),luminosity(k,1),
     &           luminosity(k,2),c_1)
            call extrap1(lumi,color_R(k+1,1),color_R(k+1,2),
     &           luminosity(k+1,1),luminosity(k+1,2),c_2)
            call extrap1(mass,c_1,c_2,mtrk(k),mtrk(k+1),c4)
            if(c4.lt.0.0) c4=0.0
            call extrap1(lumi,color_I(k,1),color_I(k,2),luminosity(k,1),
     &           luminosity(k,2),c_1)
            call extrap1(lumi,color_I(k+1,1),color_I(k+1,2),
     &           luminosity(k+1,1),luminosity(k+1,2),c_2)
            call extrap1(mass,c_1,c_2,mtrk(k),mtrk(k+1),c5)
            if(c5.lt.0.0) c5=0.0
            check3=1      
            GOTO 45
C         lum menor que las conocidas --> linear 2D extrapolation
          elseif (lumi.lt.luminosity(k,ns1).OR.lumi.lt.
     &    luminosity(k+1,ns2)) then
            call extrap1(lumi,color_U(k,ns1-1),color_U(k,ns1),
     &           luminosity(k,ns1-1),luminosity(k,ns1),c_1)
            call extrap1(lumi,color_U(k+1,ns2-1),color_U(k+1,ns2),
     &           luminosity(k+1,ns2-1),luminosity(k+1,ns2),c_2)
            call extrap1(mass,c_1,c_2,mtrk(k),mtrk(k+1),c1)
            call extrap1(lumi,color_B(k,ns1-1),color_B(k,ns1),
     &           luminosity(k,ns1-1),luminosity(k,ns1),c_1)
            call extrap1(lumi,color_B(k+1,ns2-1),color_B(k+1,ns2),
     &           luminosity(k+1,ns2-1),luminosity(k+1,ns2),c_2)
            call extrap1(mass,c_1,c_2,mtrk(k),mtrk(k+1),c2)
            call extrap1(lumi,color_V(k,ns1-1),color_V(k,ns1),
     &          luminosity(k,ns1-1),luminosity(k,ns1),c_1)
            call extrap1(lumi,color_V(k+1,ns2-1),color_V(k+1,ns2),
     &           luminosity(k+1,ns2-1),luminosity(k+1,ns2),c_2)
            call extrap1(mass,c_1,c_2,mtrk(k),mtrk(k+1),c3)
            call extrap1(lumi,color_R(k,ns1-1),color_R(k,ns1),
     &           luminosity(k,ns1-1),luminosity(k,ns1),c_1)
            call extrap1(lumi,color_R(k+1,ns2-1),color_R(k+1,ns2),
     &           luminosity(k+1,ns2-1),luminosity(k+1,ns2),c_2)
            call extrap1(mass,c_1,c_2,mtrk(k),mtrk(k+1),c4)      
            call extrap1(lumi,color_I(k,ns1-1),color_I(k,ns1),
     &           luminosity(k,ns1-1),luminosity(k,ns1),c_1)
            call extrap1(lumi,color_I(k+1,ns2-1),color_I(k+1,ns2),
     &           luminosity(k+1,ns2-1),luminosity(k+1,ns2),c_2)
            call extrap1(mass,c_1,c_2,mtrk(k),mtrk(k+1),c5)      
            check3=1
            GOTO 45
C         lum contenida entre las conocidas
          else
            do i=1,ns1-1
              if(lumi.ge.luminosity(k,i+1).AND.lumi.le.luminosity(k,i)) 
     &        then
                i1=i
                check1=1
                GOTO 25
              end if
            end do      
25          continue      
            do i=1,ns2-1
              if(lumi.ge.luminosity(k+1,i+1).AND.lumi.le.
     &        luminosity(k+1,i)) then
               i2=i
                check2=1
                GOTO 30
              end if
            end do
30          continue
            if(check1.eq.1.AND.check2.eq.1) then
              check3=1
              a1=lumi-luminosity(k,i1)
              a2=lumi-luminosity(k+1,i2)
              b1=luminosity(k,i1+1)-luminosity(k,i1)
              b2=luminosity(k+1,i2+1)-luminosity(k+1,i2)
              c_1=color_U(k,i1)+(color_U(k,i1+1)-color_U(k,i1))*a1/b1
              c_2=color_U(k+1,i2)+(color_U(k+1,i2+1)-color_U(k+1,i2))*
     &            a2/b2
              call extrap1(mass,c_1,c_2,mtrk(k),mtrk(k+1),c1)
              c_1=color_B(k,i1)+(color_B(k,i1+1)-color_B(k,i1))*a1/b1
              c_2=color_B(k+1,i2)+(color_B(k+1,i2+1)-color_B(k+1,i2))*
     &            a2/b2
             call extrap1(mass,c_1,c_2,mtrk(k),mtrk(k+1),c2)
              c_1=color_V(k,i1)+(color_V(k,i1+1)-color_V(k,i1))*a1/b1
              c_2=color_V(k+1,i2)+(color_V(k+1,i2+1)-color_V(k+1,i2))*
     &            a2/b2
              call extrap1(mass,c_1,c_2,mtrk(k),mtrk(k+1),c3)
              c_1=color_R(k,i1)+(color_R(k,i1+1)-color_R(k,i1))*a1/b1
              c_2=color_R(k+1,i2)+(color_R(k+1,i2+1)-color_R(k+1,i2))*
     &            a2/b2
              call extrap1(mass,c_1,c_2,mtrk(k),mtrk(k+1),c4)
              c_1=color_I(k,i1)+(color_I(k,i1+1)-color_I(k,i1))*a1/b1
              c_2=color_I(k+1,i2)+(color_I(k+1,i2+1)-color_I(k+1,i2))*
     &            a2/b2
              call extrap1(mass,c_1,c_2,mtrk(k),mtrk(k+1),c5)
              GOTO 45
            end if
          end if
        end if
      end do
45    continue

      if(check3.eq.0) write(1007,*) 'ERROR',lumi,mass

      return
      end
C***********************************************************************




C***********************************************************************
C     TODO:rewrite      
      subroutine extrap1(lumi,x1,x2,l1,l2,c)
      
      implicit double precision (a-h,m,o-z)
      double precision lumi,x1,x2,l1,l2,c,s,b
      
      s=(x2-x1)/(l2-l1)
      b=x1-s*l1
      c=s*lumi+b
         
      return
      end
C***********************************************************************


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


C***********************************************************************
C     TODO: rewrite
      subroutine colordb(ncolb,numberOfPointsInSequence,
     &           massSequence,luminosityDB,colorDB_U,colorDB_B,
     &           colorDB_V,colorDB_R,colorDB_I)
C=======================================================================
C     This subroutine reads the color tables by Bergeron for the  
C     sequences of WD DB
C
C     Created by ER Cojocaru (10/12)
C-----------------------------------------------------------------------
C     Input parameters:
C       ncolb: number of sequences
C-----------------------------------------------------------------------
C     Output parameters:
C       QUESTION: did O give right name to this variable?
C       massSequence: masses seq.
C       numberOfPointsInSequence
C       luminosityDB: luminosity
C       colorDB_U,colorDB_B,colorDB_V,colorDB_R,colorDB_I: colors DBs 
C
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


C***********************************************************************
C     TODO: rewrite 
      SUBROUTINE incoolone
C=======================================================================
C
C       Reading the cooling tables of ONe WD's
C
C       Revised in 27.09.07 by S. Torres      
C
C-----------------------------------------------------------------------
C       mtabone:  mass of the WD of ONE 
C       ltabone:  logarithm of the luminosity
C       mvtabone: visual absolute magnitude visual
C       lgtabone: logarithm of the cooling time in years
C
C=======================================================================
      implicit double precision (a-h,m,o-z)
      
      integer nrow,ncol,nrow2,ncol2
      double precision a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13 

C   ---  Parameter ---
      parameter(ncol=6)
      parameter(ncol2=5)
      parameter(nrow=300)
      parameter(nrow2=300)

C   ---  Dimension  ---
      integer i,j,ird
      integer ndatsone(ncol),ndatsone2(ncol2)
      double precision mtabone(ncol),mtabone2(ncol2)
      double precision lgtabone(ncol,nrow),ltabone(ncol,nrow)
      double precision mvtabone(ncol,nrow),lgtetabone(ncol,nrow)
      double precision bvtabone(ncol,nrow),vitabone(ncol,nrow)
      double precision vrtabone(ncol,nrow),uvtabone(ncol,nrow)
      double precision lgrtabone(ncol2,nrow2),lgt2tabone(ncol2,nrow2)

C  ---  Commons ---
      common /fredone/ lgtabone,ltabone,mvtabone,lgtetabone
      common /fredone2/ mtabone,ndatsone
      common /colorsone/ bvtabone,vitabone,vrtabone,uvtabone
      common /newone/ lgrtabone,lgt2tabone
      common /newone2/ mtabone2,ndatsone2

C   ---  ird is the unit of reading archives ---
      ird=121
 
C   ---  Masses of each archive ---
      mtabONE(1)=1.06d0
      mtabONE(2)=1.10d0
      mtabONE(3)=1.16d0
      mtabONE(4)=1.20d0
      mtabONE(5)=1.24d0
      mtabONE(6)=1.28d0
      mtabone2(1)=1.06d0
      mtabone2(2)=1.10d0
      mtabone2(3)=1.16d0
      mtabone2(4)=1.20d0
      mtabone2(5)=1.28d0

C   ---   Reading the files ---
      do  i=1,ncol
        do j=1,nrow
          read(ird,*,end=2) a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13
          ltabone(i,j)=a1
          mvtabone(i,j)=a12
          lgtabone(i,j)=a13
          lgtetabone(i,j)=a2
          bvtabone(i,j)=a3
          vitabone(i,j)=a9
          vrtabone(i,j)=a4
          uvtabone(i,j)=a10
        end do      
2       ndatsONE(i)=j-1  
        ird=ird+1
      end do

C      ---  Reading the data of log Teff and log g ---
      ird=127      
      do  i=1,ncol2
        do j=1,nrow2
          read(ird,*,end=3) a1,a2,a3,a4,a5,a6
C         ---  Reconverting radii in cm to radii in solar radius ---
          a3=10.0**a3
          a3=a3/(6.96d10)
          a3=dlog10(a3)
          lgrtabone(i,j)=a3
          lgt2tabone(i,j)=a5 
        end do      
3       continue       
        ndatsONE2(i)=j-1  
        ird=ird+1
      end do

      return
      end
C***********************************************************************



C***********************************************************************
C     TODO: rewrite       
      subroutine interlumdb(tcool,mass,Z,lum,c1,c2,c3,c4,c5,teff,logg)
C=======================================================================
C
C     This subroutine interpolates luminosity of DB WD's
C     using the cooling tables by Althaus et al. (2009)
C     for metallicities 0.001, 0.01, 0.06
C
C     Created by S. Torres
C     Modifications 10/2012 (ER Cojocaru):
C     * interpolation for subroutine interp
C     * interpolation using different metallicities
C
C-------------------------------------------------------------------
C     Input parameters
C       tcool: cooling time
C       mass: mass of the WD
C-------------------------------------------------------------------
C     Output parameters
C       lum: luminosity
C       teff: effective temperature
C       xlog: log(g)
C       c1,c2,c3,c4,c5: Johnson colors Johnson (UBVRI)
C
C=======================================================================
      implicit double precision (a-h,m,o-z)

C     ---   Declaration of variables   ---
      integer numberOfSequencesInGroup_1,numberOfSequencesInGroup_2,
     &        numberOfSequencesInGroup_3,numberOfSequences,nrowb,nrowb2,
     &        model,modlog   
      double precision tcool,mass,Z,lum,teff,logg
      double precision Z1,lum1,teff1,logg1
      double precision Z2,lum2,teff2,logg2
      double precision c1,c2,c3,c4,c5
      double precision zet1,zet2,zet3
      
C     ---   Parameters   ---
      parameter (nrowb=400)
      parameter (nrowb2=60)

C     ---   Dimensions   ---
      integer vectorOfPointsNumberOfSeq_1(7),
     &        vectorOfPointsNumberOfSeq_2(9),
     &        vectorOfPointsNumberOfSeq_3(9),
     &        numberOfPointsInSequence(7)
      double precision vectorOfMasses_1(7),vectorOfPreviousTimes_1(7)
      double precision matrixOfCoolingTimes_1(7,nrowb),
     &                 matrixOfLuminosities_1(7,nrowb)
      double precision matrixOfEffectiveTemperatures_1(7,nrowb),
     &                 matrixOfLog_g_1(7,nrowb)
      double precision vectorOfMasses_2(9),vectorOfPreviousTimes_2(9)
      double precision matrixOfCoolingTimes_2(9,nrowb),
     &                 matrixOfLuminosities_2(9,nrowb)
      double precision matrixOfEffectiveTemperatures_2(9,nrowb),
     &                 matrixOfLog_g_2(9,nrowb)
      double precision vectorOfMasses_3(9),vectorOfPreviousTimes_3(9)
      double precision matrixOfCoolingTimes_3(9,nrowb),
     &                 matrixOfLuminosities_3(9,nrowb)
      double precision matrixOfEffectiveTemperatures_3(9,nrowb),
     &                 matrixOfLog_g_3(9,nrowb)
      double precision massSequence(7),luminosityDB(7,nrowb2)
      double precision colorDB_U(7,nrowb2),colorDB_B(7,nrowb2),
     &                 colorDB_V(7,nrowb2),colorDB_R(7,nrowb2),
     &                 colorDB_I(7,nrowb2)

C     ---   Commons   ---
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
      common /dbtprewd/ vectorOfPreviousTimes_1,
     &                  vectorOfPreviousTimes_2,vectorOfPreviousTimes_3
      common /dbcolors/ luminosityDB,colorDB_U,colorDB_B,colorDB_V,
     &                  colorDB_R,colorDB_I
     
      zet1=0.001
      zet2=0.01
      zet3=0.06
     
C     ---  Determinating luminosities, magnitudes for every star ---
      model=0
      modlog=0
      
C     Checking between what values of metallicity is the Z of input
C     calculating lum,teff,logg for those values Z1 and Z2 and then 
C     interpolating for Z 
      if(Z.ge.zet1.AND.Z.lt.zet2) then
        Z1=zet1
        call interp(model,modlog,tcool,mass,
     &       numberOfSequencesInGroup_1,vectorOfPointsNumberOfSeq_1,
     &       matrixOfCoolingTimes_1,vectorOfPreviousTimes_1,
     &       vectorOfMasses_1,matrixOfLuminosities_1,lum1)
        modlog=1
        call interp(model,modlog,tcool,mass,
     &       numberOfSequencesInGroup_1,vectorOfPointsNumberOfSeq_1,
     &       matrixOfCoolingTimes_1,vectorOfPreviousTimes_1,
     &       vectorOfMasses_1,matrixOfEffectiveTemperatures_1,teff1)
        modlog=0
        call interp(model,modlog,tcool,mass,
     &       numberOfSequencesInGroup_1,vectorOfPointsNumberOfSeq_1,
     &       matrixOfCoolingTimes_1,vectorOfPreviousTimes_1,
     &       vectorOfMasses_1,matrixOfLog_g_1,logg1)

        Z2=zet2
        call interp(model,modlog,tcool,mass,
     &       numberOfSequencesInGroup_2,vectorOfPointsNumberOfSeq_2,
     &       matrixOfCoolingTimes_2,vectorOfPreviousTimes_2,
     &       vectorOfMasses_2,matrixOfLuminosities_2,lum2)
        modlog=1
        call interp(model,modlog,tcool,mass,
     &       numberOfSequencesInGroup_2,vectorOfPointsNumberOfSeq_2,
     &       matrixOfCoolingTimes_2,vectorOfPreviousTimes_2,
     &       vectorOfMasses_2,matrixOfEffectiveTemperatures_2,teff2)
        modlog=0
        call interp(model,modlog,tcool,mass,
     &       numberOfSequencesInGroup_2,vectorOfPointsNumberOfSeq_2,
     &       matrixOfCoolingTimes_2,vectorOfPreviousTimes_2,
     &       vectorOfMasses_2,matrixOfLog_g_2,logg2)
      end if
      
      if(Z.ge.zet2.AND.Z.lt.zet3) then
        Z1=zet2
        call interp(model,modlog,tcool,mass,
     &       numberOfSequencesInGroup_2,vectorOfPointsNumberOfSeq_2,
     &       matrixOfCoolingTimes_2,vectorOfPreviousTimes_2,
     &       vectorOfMasses_2,matrixOfLuminosities_2,lum1)
        modlog=1
        call interp(model,modlog,tcool,mass,
     &       numberOfSequencesInGroup_2,vectorOfPointsNumberOfSeq_2,
     &       matrixOfCoolingTimes_2,vectorOfPreviousTimes_2,
     &       vectorOfMasses_2,matrixOfEffectiveTemperatures_2,teff1)
        modlog=0
        call interp(model,modlog,tcool,mass,
     &       numberOfSequencesInGroup_2,vectorOfPointsNumberOfSeq_2,
     &       matrixOfCoolingTimes_2,vectorOfPreviousTimes_2,
     &       vectorOfMasses_2,matrixOfLog_g_2,logg1)
      Z2=zet3
        call interp(model,modlog,tcool,mass,numberOfSequencesInGroup_3,
     &       vectorOfPointsNumberOfSeq_3,matrixOfCoolingTimes_3,
     &       vectorOfPreviousTimes_3,vectorOfMasses_3,
     &       matrixOfLuminosities_3,lum2)
        modlog=1
        call interp(model,modlog,tcool,mass,numberOfSequencesInGroup_3,
     &       vectorOfPointsNumberOfSeq_3,matrixOfCoolingTimes_3,
     &       vectorOfPreviousTimes_3,vectorOfMasses_3,
     &       matrixOfEffectiveTemperatures_3,teff2)
        modlog=0
        call interp(model,modlog,tcool,mass,numberOfSequencesInGroup_3,
     &       vectorOfPointsNumberOfSeq_3,matrixOfCoolingTimes_3,
     &       vectorOfPreviousTimes_3,vectorOfMasses_3,matrixOfLog_g_3,
     &       logg2)
      end if
      
C     interpolation 3D in function of Z
      lum=lum1+(lum2-lum1)*(Z-Z1)/(Z2-Z1)
      teff=teff1+(teff2-teff1)*(Z-Z1)/(Z2-Z1)
      logg=logg1+(logg2-logg1)*(Z-Z1)/(Z2-Z1)
      
C     colors interpolation, here there is no metallicity
      call intermag(mass,lum,numberOfSequences,
     &     numberOfPointsInSequence,massSequence,luminosityDB,colorDB_U,
     &     colorDB_B,colorDB_V,colorDB_R,colorDB_I,c1,c2,c3,c4,c5)

      return
      end
C***********************************************************************




C***********************************************************************
C     TODO: rewrite
      subroutine interlumONe(tcool,mass,lum,c1,c2,c3,c4,c5,teff,xlog)
C=======================================================================
C
C     This subroutine interpolates luminosity of the DA ONe WD's
C     together with Johnson colors (U,B,V,R,I)
C
C     Created by S. Torres
C     Modifications 10/2012 (ER Cojocaru) - subtracted subroutine 
C                                             wdcoolone,
C                                    interpolation for subroutine interp
C
C-------------------------------------------------------------------
C     Input parameters:
C       tcool: cooling time
C       mass: mass of the WD
C-------------------------------------------------------------------
C     Output parameters:
C       lum: luminosity
C       c1,c2,c3,c4,c5: Johnson colors (U,B,V,R,I)
C       teff: effective temperature [K]
C       log g: logarithm of the superficial gravity
C
C=======================================================================
      implicit double precision (a-h,m,o-z)

C     ---   Declaration of variables   ---
      integer ncol,ncol2,nrow,nrow2,i,model,modlog
      double precision mass,tcool,y
      double precision mv,cbv,cvi,cvr,cuv,zlte,zlr,G
      double precision lum,c1,c2,c3,c4,c5,teff,xlog

C     ---   Parameters   ---
      parameter (ncol=6)
      parameter (ncol2=5)
      parameter (nrow=300)
      parameter (nrow2=300)

C     ---   Dimensions   ---
      integer ndatsone(ncol),ndatsone2(ncol2)
      double precision mtabone(ncol),mtabone2(ncol2)
      double precision mvtabone(ncol,nrow),ltabone(ncol,nrow)
      double precision lgtabone(ncol,nrow),lgtetabone(ncol,nrow)
      double precision bvtabone(ncol,nrow),vitabone(ncol,nrow)
      double precision vrtabone(ncol,nrow),uvtabone(ncol,nrow)
      double precision lgrtabone(ncol2,nrow2),lgt2tabone(ncol2,nrow2)
      double precision tprewd1(ncol),tprewd2(ncol2)
 
C     ---   Commons   ---
      common /fredone/ lgtabone,ltabone,mvtabone,lgtetabone
      common /fredone2/ mtabone,ndatsone
      common /colorsone/ bvtabone,vitabone,vrtabone,uvtabone
      common /newone/ lgrtabone,lgt2tabone
      common /newone2/ mtabone2,ndatsone2
      
C     ---   Interpolation  ---
      model=1
      modlog=0
      y=dlog10(tcool)+9.0

      do i=1,ncol
        tprewd1(i)=0.0
      end do

      do i=1,ncol2
        tprewd2(i)=0.0
      end do

      call interp(model,modlog,y,mass,ncol,ndatsone,lgtabone,tprewd1,
     &     mtabone,ltabone,lum)
      call interp(model,modlog,y,mass,ncol,ndatsone,lgtabone,tprewd1,
     &     mtabone,mvtabone,mv)
      call interp(model,modlog,y,mass,ncol,ndatsone,lgtabone,tprewd1,
     &     mtabone,bvtabone,cbv)
      call interp(model,modlog,y,mass,ncol,ndatsone,lgtabone,tprewd1,
     &     mtabone,vitabone,cvi)
      call interp(model,modlog,y,mass,ncol,ndatsone,lgtabone,tprewd1,
     &     mtabone,vrtabone,cvr)
      call interp(model,modlog,y,mass,ncol,ndatsone,lgtabone,tprewd1,
     &     mtabone,uvtabone,cuv)
      call interp(model,modlog,y,mass,ncol,ndatsone,lgtabone,tprewd1,
     &     mtabone,lgtetabone,zlte)
      call interp(model,modlog,y,mass,ncol2,ndatsone2,lgt2tabone,
     &     tprewd2,mtabone2,lgrtabone,zlr)

      teff=10.0**zlte

C   --- G in cm/s² kg; M in kg and R in cm ---
      G=6.67d-5
      xlog=dlog10(G)+dlog10(mass*(1.989d30))-2.0*(zlr+dlog10(6.696d10))
      c1=cuv+mv
      c2=cbv+mv
      c3=mv
      c4=mv-cvr
      c5=mv-cvi
       
      return
      end
C***********************************************************************




C***********************************************************************
C     TODO: rewrite      
      subroutine interp(model,modlog,y,m,ncol,ntrk,ttrkk,tprewd,mtrk,
     &           xtrk,xout)
C=======================================================================
C
C     This subroutine interpolates a certain related measure with the EB
C     DB or ONe, according to input data, considering all the cases of
C     interpolation and extrapolation needed.
C
C     Created by ER Cojocaru (10/2012)
C
C-------------------------------------------------------------------
C     Input parameters:
C       model: DA,DB (0) or ONe (1)
C       QUESTION: what does it mean?
C       modlog: hay que usar log de la medida (1 - teff) o no (0) (DA/DB)
C       y: cooling time
C       m: mass of the WD
C       ncol: number of sequences(columns) available
C       ntrk: number de rows for sequence
C       ttrkk: known cooling times in the sequence
C       mtrk: known masses in the sequence
C       xtrk: values of the measure of interest known in sequence
C-------------------------------------------------------------------
C     Output parameters:
C       xout: measure of interest
C=======================================================================
      implicit double precision (a-h,m,o-z)

C     ---   Declaration de variables   ---
      integer ncol,ii,jj,k,ns,j1,j2,case1,case2,model,modlog
C     TIME      
      double precision y,y1,y2,y3,y4,ym1,ym2
C     MASS
      double precision m,m1,m2
C     Measure of interest
      double precision x1,x2,x3,x4,xm1,xm2,xout
      
      double precision deltf,den,s,b,t

C     ---   Dimensions   ---
      integer ntrk(ncol)
      double precision mtrk(ncol),ttrkk(ncol,*)
      double precision tprewd(ncol),xtrk(ncol,*)
     
C     ---  Interpolating measure from tcool and mwd  ---

C     ---   Mass less than minimum   ---
      if (m .lt. mtrk(1)) then
        jj=1
        m1=mtrk(jj)
C       --- Checking if the time is less than the minimum---
        if(y.lt.ttrkk(jj,1)) then
C         making linear extrapolation in log-log
          if(model.eq.0) then
            if(modlog.eq.1) then
              deltf=log10((ttrkk(jj,2)+tprewd(jj))/(ttrkk(jj,1)+
     &              tprewd(jj)))
              s=log10(xtrk(jj,2)/xtrk(jj,1))/deltf
              b=log10(xtrk(jj,2))-s*log10(ttrkk(jj,2)+tprewd(jj))      
              xm1=10.0**(s*log10(y+tprewd(jj))+b)
            else
              deltf=log10((ttrkk(jj,2)+tprewd(jj))/(ttrkk(jj,1)+
     &              tprewd(jj)))
              s=(xtrk(jj,2)-xtrk(jj,1))/deltf
              b=xtrk(jj,2)-s*log10(ttrkk(jj,2)+tprewd(jj))      
              xm1=s*log10(y+tprewd(jj))+b
            endif
          else
            deltf=ttrkk(jj,2)-ttrkk(jj,1)
            s=(xtrk(jj,2)-xtrk(jj,1))/deltf
            b=xtrk(jj,2)-s*ttrkk(jj,2)
            xm1=s*y+b
          end if
          goto 20
        endif      
C       --- Checking if the time is greater than maximum ---
        ns=ntrk(jj)
        if(y.gt.ttrkk(jj,ns)) then
C         making linear extrapolation in log-log
          if(model.eq.0) then
            if(modlog.eq.1) then
              deltf=log10((ttrkk(jj,ns)+tprewd(jj))/(ttrkk(jj,ns-1)+
     &              tprewd(jj)))
              s=log10(xtrk(jj,ns)/xtrk(jj,ns-1))/deltf
              b=log10(xtrk(jj,ns))-s*log10(ttrkk(jj,ns)+tprewd(jj))      
              xm1=10.0**(s*log10(y+tprewd(jj))+b)
            else
              deltf=log10((ttrkk(jj,ns)+tprewd(jj))/(ttrkk(jj,ns-1)+
     &              tprewd(jj)))
              s=(xtrk(jj,ns)-xtrk(jj,ns-1))/deltf
              b=xtrk(jj,ns)-s*log10(ttrkk(jj,ns)+tprewd(jj))      
              xm1=s*log10(y+tprewd(jj))+b
            endif
          else
            deltf=ttrkk(jj,ns)-ttrkk(jj,ns-1)
            s=(xtrk(jj,ns)-xtrk(jj,ns-1))/deltf
            b=xtrk(jj,ns)-s*ttrkk(jj,ns)
            xm1=s*y+b
          end if
          goto 20
        endif
        do k=1,ns-1
          ii=k
C         ---   Interpolation of the tiempo   ---
          if (y .ge. ttrkk(jj,ii) .and. y .lt. ttrkk(jj,ii+1)) then            
            y1=ttrkk(jj,ii)
            y2=ttrkk(jj,ii+1)
            x1=xtrk(jj,ii)
            x2=xtrk(jj,ii+1)       
            deltf=(y-y1)/(y2-y1)
            xm1=x1+deltf*(x2-x1)   
            goto 20
          end if
        end do
20      CONTINUE
        jj=2
        m2=mtrk(jj)
C       --- Checking if the time is less than minimum---
        if(y.lt.ttrkk(jj,1)) then
C         hacemos un extrapolacion lineal en log-log
          if(model.eq.0) then
            if(modlog.eq.1) then
              deltf=log10((ttrkk(jj,2)+tprewd(jj))/(ttrkk(jj,1)+
     &              tprewd(jj)))
              s=log10(xtrk(jj,2)/xtrk(jj,1))/deltf
              b=log10(xtrk(jj,2))-s*log10(ttrkk(jj,2)+tprewd(jj))      
              xm2=10.0**(s*log10(y+tprewd(jj))+b)
            else
              deltf=log10((ttrkk(jj,2)+tprewd(jj))/(ttrkk(jj,1)+
     &              tprewd(jj)))
              s=(xtrk(jj,2)-xtrk(jj,1))/deltf
              b=xtrk(jj,2)-s*log10(ttrkk(jj,2)+tprewd(jj))      
              xm2=s*log10(y+tprewd(jj))+b
            endif
          else
            deltf=ttrkk(jj,2)-ttrkk(jj,1)
            s=(xtrk(jj,2)-xtrk(jj,1))/deltf
            b=xtrk(jj,2)-s*ttrkk(jj,2)
            xm2=s*y+b
          end if
          goto 30
        endif 
C       --- Checking if the time is greater than maximum---
        ns=ntrk(jj)
        if(y.gt.ttrkk(jj,ns)) then
C         making linear extrapolation in log-log
          if(model.eq.0) then
            if(modlog.eq.1) then
              deltf=log10((ttrkk(jj,ns)+tprewd(jj))/(ttrkk(jj,ns-1)+
     &              tprewd(jj)))
              s=log10(xtrk(jj,ns)/xtrk(jj,ns-1))/deltf
              b=log10(xtrk(jj,ns))-s*log10(ttrkk(jj,ns)+tprewd(jj))      
              xm2=10.0**(s*log10(y+tprewd(jj))+b)
            else
              deltf=log10((ttrkk(jj,ns)+tprewd(jj))/(ttrkk(jj,ns-1)+
     &              tprewd(jj)))
              s=(xtrk(jj,ns)-xtrk(jj,ns-1))/deltf
              b=xtrk(jj,ns)-s*log10(ttrkk(jj,ns)+tprewd(jj))      
              xm2=s*log10(y+tprewd(jj))+b
            endif
          else
            deltf=ttrkk(jj,ns)-ttrkk(jj,ns-1)
            s=(xtrk(jj,ns)-xtrk(jj,ns-1))/deltf
            b=xtrk(jj,ns)-s*ttrkk(jj,ns)
            xm2=s*y+b
          end if
          goto 30
        endif
        do k=1,ntrk(jj)-1
          ii=k
C         ---   Interpolation of the time   ---
          if (y .ge. ttrkk(jj,ii) .and. y .lt. ttrkk(jj,ii+1)) then
            y1=ttrkk(jj,ii)
            y2=ttrkk(jj,ii+1)
            x1=xtrk(jj,ii)
            x2=xtrk(jj,ii+1)  
            deltf=(y-y1)/(y2-y1)
            xm2=x1+deltf*(x2-x1) 
            goto 30
          end if
        end do
30      CONTINUE
        s=(xm2-xm1)/(m2-m1)
        t=xm2-s*m2
        xout=s*m+t
        goto 100
      end if
        
C     ---   Mass greater than maximum---
      if (m .ge. mtrk(ncol)) then
        jj=ncol-1
        m1=mtrk(jj)
C       --- Checking if time is less than minimum---
        if(y.lt.ttrkk(jj,1)) then
C         making linear extrapolation in log-log
          if(model.eq.0) then
            if(modlog.eq.1) then
              deltf=log10((ttrkk(jj,2)+tprewd(jj))/(ttrkk(jj,1)+
     &              tprewd(jj)))
              s=log10(xtrk(jj,2)/xtrk(jj,1))/deltf
              b=log10(xtrk(jj,2))-s*log10(ttrkk(jj,2)+tprewd(jj))      
              xm1=10.0**(s*log10(y+tprewd(jj))+b)
            else
              deltf=log10((ttrkk(jj,2)+tprewd(jj))/(ttrkk(jj,1)+
     &              tprewd(jj)))
              s=(xtrk(jj,2)-xtrk(jj,1))/deltf
              b=xtrk(jj,2)-s*log10(ttrkk(jj,2)+tprewd(jj))      
              xm1=s*log10(y+tprewd(jj))+b
            endif
          else
            deltf=ttrkk(jj,2)-ttrkk(jj,1)
            s=(xtrk(jj,2)-xtrk(jj,1))/deltf
            b=xtrk(jj,2)-s*ttrkk(jj,2)
            xm1=s*y+b
          end if
          goto 200
        endif      
C       --- Checking if time is greater than maximum---
        ns=ntrk(jj)
        if(y.gt.ttrkk(jj,ns)) then
C         making linear extrapolation in log-log
          if(model.eq.0) then
            if(modlog.eq.1) then
              deltf=log10((ttrkk(jj,ns)+tprewd(jj))/(ttrkk(jj,ns-1)+
     &              tprewd(jj)))
              s=log10(xtrk(jj,ns)/xtrk(jj,ns-1))/deltf
              b=log10(xtrk(jj,ns))-s*log10(ttrkk(jj,ns)+tprewd(jj))      
              xm1=10.0**(s*log10(y+tprewd(jj))+b)
            else
              deltf=log10((ttrkk(jj,ns)+tprewd(jj))/(ttrkk(jj,ns-1)+
     &              tprewd(jj)))
              s=(xtrk(jj,ns)-xtrk(jj,ns-1))/deltf
              b=xtrk(jj,ns)-s*log10(ttrkk(jj,ns)+tprewd(jj))      
              xm1=s*log10(y+tprewd(jj))+b
            endif
          else
            deltf=ttrkk(jj,ns)-ttrkk(jj,ns-1)
            s=(xtrk(jj,ns)-xtrk(jj,ns-1))/deltf
            b=xtrk(jj,ns)-s*ttrkk(jj,ns)
            xm1=s*y+b
          end if
          goto 200
        endif
        do k=1,ns-1
          ii=k
C         ---   Interpolation of time   ---
          if (y .ge. ttrkk(jj,ii) .and. y .lt. ttrkk(jj,ii+1)) then
            y1=ttrkk(jj,ii)
            y2=ttrkk(jj,ii+1)
            x1=xtrk(jj,ii)
            x2=xtrk(jj,ii+1)       
            deltf=(y-y1)/(y2-y1)
            xm1=x1+deltf*(x2-x1)  
            goto 200
          end if
        end do
200     CONTINUE
        jj=ncol
        m2=mtrk(jj)
C       --- Checking if time is less than minimum ---
        if(y.lt.ttrkk(jj,1)) then
C       making linear extrapolation in log-log
          if(model.eq.0) then
            if(modlog.eq.1) then
              deltf=log10((ttrkk(jj,2)+tprewd(jj))/(ttrkk(jj,1)+
     &              tprewd(jj)))
              s=log10(xtrk(jj,2)/xtrk(jj,1))/deltf
              b=log10(xtrk(jj,2))-s*log10(ttrkk(jj,2)+tprewd(jj))      
              xm2=10.0**(s*log10(y+tprewd(jj))+b)
            else
              deltf=log10((ttrkk(jj,2)+tprewd(jj))/(ttrkk(jj,1)+
     &              tprewd(jj)))
              s=(xtrk(jj,2)-xtrk(jj,1))/deltf
              b=xtrk(jj,2)-s*log10(ttrkk(jj,2)+tprewd(jj))      
              xm2=s*log10(y+tprewd(jj))+b
            endif
          else
            deltf=ttrkk(jj,2)-ttrkk(jj,1)
            s=(xtrk(jj,2)-xtrk(jj,1))/deltf
            b=xtrk(jj,2)-s*ttrkk(jj,2)
            xm2=s*y+b
          end if      
          goto 300
        endif 
C       --- Checking if time is greater than maximum ---
        ns=ntrk(jj)
        if(y.gt.ttrkk(jj,ns)) then
C         making linear extrapolation in log-log
          if(model.eq.0) then
            if(modlog.eq.1) then
              deltf=log10((ttrkk(jj,ns)+tprewd(jj))/(ttrkk(jj,ns-1)+
     &              tprewd(jj)))
              s=log10(xtrk(jj,ns)/xtrk(jj,ns-1))/deltf
              b=log10(xtrk(jj,ns))-s*log10(ttrkk(jj,ns)+tprewd(jj))      
              xm2=10.0**(s*log10(y+tprewd(jj))+b)
            else
              deltf=log10((ttrkk(jj,ns)+tprewd(jj))/(ttrkk(jj,ns-1)+
     &              tprewd(jj)))
              s=(xtrk(jj,ns)-xtrk(jj,ns-1))/deltf
              b=xtrk(jj,ns)-s*log10(ttrkk(jj,ns)+tprewd(jj))      
              xm2=s*log10(y+tprewd(jj))+b
            endif
          else
            deltf=ttrkk(jj,ns)-ttrkk(jj,ns-1)
            s=(xtrk(jj,ns)-xtrk(jj,ns-1))/deltf
            b=xtrk(jj,ns)-s*ttrkk(jj,ns)
            xm2=s*y+b
          end if
          goto 300
        endif
        do k=1,ntrk(jj)-1
          ii=k
C         ---  Interpolation of time   ---
          if (y .ge. ttrkk(jj,ii) .and. y .lt. ttrkk(jj,ii+1)) then      
            y1=ttrkk(jj,ii)
            y2=ttrkk(jj,ii+1)
            x1=xtrk(jj,ii)
            x2=xtrk(jj,ii+1)  
            deltf=(y-y1)/(y2-y1)
            xm2=x1+deltf*(x2-x1)  
            goto 300
          end if
        end do
300     CONTINUE
        s=(xm2-xm1)/(m2-m1)
        t=xm2-s*m2
        xout=s*m+t
        goto 100
      end if
  
C     ---   Search for masses of interpolation ---
      do 3 k=1,ncol-1
        j1=k
        if (m .ge. mtrk(k) .and. m .lt. mtrk(k+1)) then
          goto 4
        end if
3     continue
4     j2=j1+1

C     ---   Search for the times of interpolation (1)   ---
C     ---   Times less than minimum  ---
      if (y. lt. ttrkk(j1,1)) then
C       making linear extrapolation in log-log
        if(model.eq.0) then
          if(modlog.eq.1) then
            deltf=log10((ttrkk(j1,2)+tprewd(j1))/(ttrkk(j1,1)+
     &            tprewd(j1)))
            s=log10(xtrk(j1,2)/xtrk(j1,1))/deltf
            b=log10(xtrk(j1,2))-s*log10(ttrkk(j1,2)+tprewd(j1))
            x1=10.0**(s*log10(y+tprewd(j1))+b)
          else
            deltf=log10((ttrkk(j1,2)+tprewd(j1))/(ttrkk(j1,1)+
     &            tprewd(j1)))
            s=(xtrk(j1,2)-xtrk(j1,1))/deltf
            b=xtrk(j1,2)-s*log10(ttrkk(j1,2)+tprewd(j1))
            x1=s*log10(y+tprewd(j1))+b
          end if
        else
          deltf=ttrkk(j1,2)-ttrkk(j1,1)
          s=(xtrk(j1,2)-xtrk(j1,1))/deltf
          b=xtrk(j1,2)-s*ttrkk(j1,2)
          x1=s*y+b
        end if
        case1=1
        goto 6      
      end if

C     ---   Times greater than maximum ---
      if(y .ge. ttrkk(j1,ntrk(j1))) then            
        ns=ntrk(j1)
C       making linear extrapolation in log-log
        if(model.eq.0) then
          if(modlog.eq.1) then
            deltf=log10((ttrkk(j1,ns)+tprewd(j1))/(ttrkk(j1,ns-1)+
     &            tprewd(j1)))
            s=log10(xtrk(j1,ns)/xtrk(j1,ns-1))/deltf
            b=log10(xtrk(j1,ns))-s*log10(ttrkk(j1,ns)+tprewd(j1))
            x1=10.0**(s*log10(y+tprewd(j1))+b)
          else
            deltf=log10((ttrkk(j1,ns)+tprewd(j1))/(ttrkk(j1,ns-1)+
     &            tprewd(j1)))
            s=(xtrk(j1,ns)-xtrk(j1,ns-1))/deltf
            b=xtrk(j1,ns)-s*log10(ttrkk(j1,ns)+tprewd(j1))
            x1=s*log10(y+tprewd(j1))+b
          end if
        else
          deltf=ttrkk(j1,ns)-ttrkk(j1,ns-1)
          s=(xtrk(j1,ns)-xtrk(j1,ns-1))/deltf
          b=xtrk(j1,ns)-s*ttrkk(j1,ns)
          x1=s*y+b
        end if
        case1=1
        goto 6
      end if

C     ---   Times between the minimum and maximum ---
      do 5 k=1,ntrk(j1)-1
        if (y .ge. ttrkk(j1,k) .and. y .lt. ttrkk(j1,k+1)) then
          y1=ttrkk(j1,k  )
          y2=ttrkk(j1,k+1)
          x1=xtrk(j1,k  )
          x2=xtrk(j1,k+1)
          case1=0
          goto 6
        end if
5     continue
 
            
C     ---   Search for interpolation times (2)   ---
C     ---   Times less than minimum---
6     if (y. lt. ttrkk(j2,1)) then
C       making linear extrapolation in log-log
        if(model.eq.0) then
          if(modlog.eq.1) then
            deltf=log10((ttrkk(j2,2)+tprewd(j2))/(ttrkk(j2,1)+
     &            tprewd(j2)))
            s=log10(xtrk(j2,2)/xtrk(j2,1))/deltf
            b=log10(xtrk(j2,2))-s*log10(ttrkk(j2,2)+tprewd(j2))
            x3=10.0**(s*log10(y+tprewd(j2))+b)
          else
            deltf=log10((ttrkk(j2,2)+tprewd(j2))/(ttrkk(j2,1)+
     &            tprewd(j2)))
            s=(xtrk(j2,2)-xtrk(j2,1))/deltf
            b=xtrk(j2,2)-s*log10(ttrkk(j2,2)+tprewd(j2))
            x3=s*log10(y+tprewd(j2))+b
          end if
        else
          deltf=ttrkk(j2,2)-ttrkk(j2,1)
            s=(xtrk(j2,2)-xtrk(j2,1))/deltf
            b=xtrk(j2,2)-s*ttrkk(j2,2)
            x3=s*y+b
        end if
        case2=1
        goto 8
      endif

C     ---   Times greater than maximum   ---
      if(y .ge. ttrkk(j2,ntrk(j2))) then  
        ns=ntrk(j2)
C       making linear extrapolation in log-log
        if(model.eq.0) then
          if(modlog.eq.1) then
            deltf=log10((ttrkk(j2,ns)+tprewd(j2))/(ttrkk(j2,ns-1)+
     &            tprewd(j2)))
            s=log10(xtrk(j2,ns)/xtrk(j2,ns-1))/deltf
            b=log10(xtrk(j2,ns))-s*log10(ttrkk(j2,ns)+tprewd(j2))
            x3=10.0**(s*log10(y+tprewd(j2))+b)
          else
            deltf=log10((ttrkk(j2,ns)+tprewd(j2))/(ttrkk(j2,ns-1)+
     &            tprewd(j2)))
            s=(xtrk(j2,ns)-xtrk(j2,ns-1))/deltf
            b=xtrk(j2,ns)-s*log10(ttrkk(j2,ns)+tprewd(j2))
            x3=s*log10(y+tprewd(j2))+b
          endif
        else
          deltf=ttrkk(j2,ns)-ttrkk(j2,ns-1)
          s=(xtrk(j2,ns)-xtrk(j2,ns-1))/deltf
          b=xtrk(j2,ns)-s*ttrkk(j2,ns)
          x3=s*y+b
        end if
        case2=1
        goto 8
      end if
      
C     ---   Times between the minimum and maximum ---
      do 7 k=1,ntrk(j2)-1
        if (y .ge. ttrkk(j2,k) .and. y .lt. ttrkk(j2,k+1)) then
          y3=ttrkk(j2,k  )
          y4=ttrkk(j2,k+1)
          x3=xtrk(j2,k  )
          x4=xtrk(j2,k+1)
          case2=0      
          goto 8
        end if
7     continue

C     ---   Bilinear interpolation ---
8     m1=mtrk(j1)
      m2=mtrk(j2)
      den=(m-m1)/(m2-m1)
      if(case1.eq.0.AND.case2.eq.0) then
        ym1=y1+(y3-y1)*den
        ym2=y2+(y4-y2)*den
        xm1=x1+(x3-x1)*den
        xm2=x2+(x4-x2)*den
        xout=xm1+(y-ym1)/(ym2-ym1)*(xm2-xm1)
      elseif(case1.eq.0.AND.case2.eq.1) then
        xm1=x1+(x2-x1)*(y-y1)/(y2-y1)
        xout=xm1+(x3-xm1)*den
      elseif(case1.eq.1.AND.case2.eq.0) then
        xm2=x3+(x4-x3)*(y-y3)/(y4-y3)
        xout=x1+(xm2-x1)*den
      else
        xout=x1+(x3-x1)*den
      end if
100   continue

      if(isNaN(xout)) then
        write(*,*) "Error interp(tcool,mass,model,log):",y,m,model,
     &             modlog
        stop
      end if

      return
      end
C***********************************************************************




C***********************************************************************
      subroutine chanco(V,UB,BV,VR,RI,g,xug,xgr,xri,xiz,xgi)
C---------------------------------------------------------------------
C     Transformation between the Johnson-Cousins UBVRI photometry 
C     system and the SDSS ugriz system.
C     Follow the eq. 1 to 8 from Jordi,Grebel & Ammon, 2006, A&A, 460
C---------------------------------------------------------------------
      implicit double precision (a-h,m,o-z)
      
      double precision a1,a2,a3,a5,a7,b1,b2,b3,b5,b7,c5
      double precision V,UB,BV,VR,RI,g,xug,xgr,xri,xiz,xgi

C     --- Used parameters (Table 3) ----
      parameter(a1=0.630, b1=-0.124)
      parameter(a2=1.007, b2=-0.236)
      parameter(a3=1.584, b3=-0.386)
      parameter(a5=0.750, b5=0.770, c5=0.720)
      parameter(a7=1.646, b7=-0.139)      

C     --- Performing transformations ----
      g=V+a1*BV+b1
      xug=a5*UB+b5*BV+c5
      xgr=a7*VR+b7
      xri=a2*RI+b2
      xiz=(a3-a2)*RI+(b3-b2)
      xgi=xgr+xri
      
      return
      end
C***********************************************************************



      

C***********************************************************************
C     TODO: rewrite
      subroutine volum_40pc(iseed)
C=======================================================================
C
C     Determining the Luminosity Function of 40 pc, 
C     using criteria (Limoges et al. 2015)
C
C     Based on version of 22.09.07 by S. Torres
C
C     Restrictions:
C        dec>0 Nort-Hemisphere SUPERBLINK survey
C        mu>40 mas/yr-1
C        V band limit: V<19
C     Method:
C        Number density per pc^3 and half bolometric magnitude
C        No 1/Vmax correction applied
C-------------------------------------------------------------------
C     Input parameters:
C       galacticDiskAge
C       numberOfStarsInSample
C
C     ---   Parameters  ---
C     NOTE: smth wrong with the name here
C     parameterIFMR: covered sky area
C     restrictions: muo: proper motion
C                   deco: declination
C                   pio: parallax
C
C=======================================================================
      implicit double precision (a-h,m,o-z)
      external ran
      real ran

C     ---   Declaration of variables  ---
      integer numberOfStars,ntwd,i,j
      integer r1,r2,r3,r4,r5,r32
      double precision parameterIFMR,muo,mumax,deco,pio
      double precision mbolmin,mbolinc,mbolmax,vtanmin
      double precision cont,d,hrm,gz
      double precision dmonte,errinfa,errsupa,mbol
C     NOTE: vvv is always 0
      double precision fnora,fnor,pi,rg,vmedia,vo,fi,vvv,x,xx
      double precision xya,zmedia,voinv
      
      parameter (numberOfStars=6000000)
      parameter (betaV=0.5)
C     declination (Only-north Hemisphere)
      parameter (deco=0.0)
C     MINIMUM PARALLAX BELOW WHICH WE DISCARD RESULTS (0.025<=>40 pc)
      parameter (pio=0.025)
C     BINNING WDLF old values
      parameter (mbolmin=5.75,mbolmax=20.75,mbolinc=0.5)
c      new values 
C      parameter (mbolmin=6.0,mbolmax=21.0,mbolinc=0.5)
C     MIN. TANGENCIAL VEL. BELOW WHICH WE DISCARD RESULTS (30)
      parameter (vtanmin=0)
C     MINIMUM PROPER MONTION 
      parameter (muo=0.04)
C     MAXIMUM PROPER MONTION 
      parameter (mumax=1000.0)
C     Parameters histogram masses
      parameter (xmasi=0.1)
      parameter (xmasf=1.4)
      parameter (xmasinc=0.05)

C     ---   Dimension   ---
      double precision mu(numberOfStars),arec(numberOfStars),
     &                 dec(numberOfStars)
      double precision leb(numberOfStars),meb(numberOfStars),
     &                 zeb(numberOfStars),teb(numberOfStars)
      double precision iwd(numberOfStars)
      double precision rgac(numberOfStars)
      double precision tcool(numberOfStars)
      integer nbin(70),nbins
      double precision coordinate_R(numberOfStars),
     &                 coordinate_Theta(numberOfStars),
     &                 coordinate_Zcylindr(numberOfStars)
      double precision parj(numberOfStars)
      double precision vtan(numberOfStars)
      double precision error(70),errora(70),ndfa(70)
      double precision vvmax(70)
      double precision heightPattern(numberOfStars)
      double precision go(numberOfStars),gr(numberOfStars),
     &                 v(numberOfStars)
      double precision gi(numberOfStars),ur(numberOfStars),
     &                 rz(numberOfStars)
      double precision mbin(70)
      double precision idb(numberOfStars)
      double precision xfl(19),xflo(19),xflcut(3),xflocut(3)
      double precision xflhot(11),xflohot(11)
C     bins for max-region: for synthetic and observational samples
      double precision xflMaxRegion(6), xfloMaxRegion(6)
      double precision nbinmass(26)
      double precision uu(numberOfStars), vv(numberOfStars), 
     &                 ww(numberOfStars)
C     sum of WDs velocities in specific bin, _u/_v/_w - components
      double precision sumOfWDVelocitiesInBin_u(70),
     &                 sumOfWDVelocitiesInBin_v(70),
     &                 sumOfWDVelocitiesInBin_w(70)
C     average velocity for WDs in specific bin      
      double precision averageWDVelocityInBin_u(70),
     &                 averageWDVelocityInBin_v(70),
     &                 averageWDVelocityInBin_w(70)  
C     this is used to calculate sigma (SD)
      double precision sumOfSquareDifferences_u,
     &                 sumOfSquareDifferences_v,
     &                 sumOfSquareDifferences_w
C     SD for velocities in each bin
C     NOTE: this 70 is random. array should be declared as dynamic
      double precision standardDeviation_u(70),standardDeviation_v(70),
     &                 standardDeviation_w(70) 
C     2D-array of velocities (nº of bin; newly assigned to WD nº in bin)
C     needed to calculate Standart Deviation (SD) for velocities in each 
C     bin
C     NOTE: I need dynamic array, static takes too much memory
      double precision arrayOfVelocitiesForSD_u(25,50000)
      double precision arrayOfVelocitiesForSD_v(25,50000)
      double precision arrayOfVelocitiesForSD_w(25,50000)
C     2D-array of bolometric magnitudes for each WD; indexes are the 
C     same as for arrayOfVelocitiesForSD_u/v/w. (For cloud)
      double precision arrayOfMagnitudes(25,50000)


C     ---   Commons  ---
      common /enanas/ leb,meb,zeb,teb
      common /index/ iwd,ntwd      
      common /mad/ mu,arec,dec
      common /paral/ rgac
      common /coorcil/ coordinate_R,coordinate_Theta,coordinate_Zcylindr
      common /cool/ tcool
      common /veltan/ vtan
      common /patron/ heightPattern
      common /photo/ go,gr,gi,ur,rz
      common /indexdb/ idb
      common /johnson/ V
      common /param/ fractionOfDB,galacticDiskAge,parameterIMF,
     &               parameterIFMR,timeOfBurst
      common /vel/ uu,vv,ww
      
      nbins=(mbolmax-mbolmin)/mbolinc
     
C-------------------------------------------------------------------
C     ---   Initialization of variables   ---
C-------------------------------------------------------------------
C     ---  Defining pi  ---
      pi=4.0*atan(1.0d0)
      fi=180.0/pi

C     ---  Initializating ndf's  ---
      do 1 i=1,70
        error(i)=0.0d0
        ndfa(i)=0.0d0
        nbin(i)=0
        vvmax(i)=0.0
        mbin(i)=0.0
        sumOfWDVelocitiesInBin_u(i)=0.0
        sumOfWDVelocitiesInBin_v(i)=0.0
        sumOfWDVelocitiesInBin_w(i)=0.0
        averageWDVelocityInBin_u(i)=0.0
        averageWDVelocityInBin_v(i)=0.0
        averageWDVelocityInBin_w(i)=0.0
 1    continue

      do 2 i=1,26
        nbinmass(i)=0
 2    continue

C     ---  Writing data   ---
      cont=0.0
      zmedia=0.0
      vmedia=0.0
      r1=0
      r2=0
      r32=0
      r3=0
      r4=0
      r5=0
           
C-------------------------------------------------------------------
C     ---   Initiating loop  ---
C-------------------------------------------------------------------
      do 5 i=1,ntwd
        d=dec(i)*fi
        rg=rgac(i)*1000.0
        parj(i)=1.0/rg
        vtan(i)=4.74*mu(i)*rg

C       --- Full sample : we  eliminate the observational cuts ---
C        go to 93212

C     ---   Restrictions of the sample  ---

C     ---   1) Eliminate WDs with parallax <pio'' 
        if (parj(i).lt.pio) then   
          r1=r1+1
          go to 5        
        end if

C       ---   4) Eliminate declination  ---          
        if (dec(i).lt.deco) then   
          r2=r2+1
          go to 5        
        end if
        write (73,*) hrm,gz
C        goto 93212

C       ---   2) Minimum proper motion cut 
        if (mu(i).lt.muo) then   
          r3=r3+1
          go to 5  
        end if

C       --- Reduced proper motion ---
        hrm=go(i)+5.0*dlog10(mu(i))+5.0
        gz=gr(i)+rz(i)
        if(gz.lt.-0.33) then
          if(hrm.lt.14.0) then 
            r32=r32+1
            go to 5
          endif
        else 
          if(hrm.lt.(3.559*gz+15.17)) then
            if(v(i).gt.12.0) then
              xxx=ran(iseed)
            endif
            go to 5
            r32=r32+1
          endif
        endif
 453    write (160,*) 2.5*leb(i)+4.75,rgac(i)  
C       QUESTION: what does it mean? 

C       ---   3) Restriction V (de momento lo hacemos con go)---
        if(v(i).ge.19.0) then 
          r4=r4+1
          goto 5
        endif      
C       QUESTION: is it really z-coordinate or maybe metallicity?        
C       1    2   3 4    5    6   7   8   9   10   11  12   13   14 15   
C       Meb -Lum Z Mbol Gap0 g-i g-r u-r r-z arec dec rgac parj mu vtan   
C       16    17   18  19 20 21 22
C       tcool temp idb coordinate_Zcylindr  uu vv ww 
        write(156,*)  meb(i),leb(i),zeb(i),2.5*leb(i)+4.75,go(i),gi(i),
     &                gr(i),ur(i),rz(i),arec(i),dec(i),rgac(i),parj(i),
     &                mu(i),vtan(i),tcool(i),teb(i),idb(i),
     &                coordinate_Zcylindr(i),uu(i),vv(i),ww(i)
C       velocities output
93212   write(1156,*)  uu(i),vv(i),ww(i)     
      continue
C     ------------------------------------------------------------------


C     ---  Making histogram of the mass---
      K=0
  40  k=k+1 
      
      xi=xmasi+dfloat(k-1)*xmasinc
      xf=xi+xmasinc
      if(meb(i).gt.xi.and.meb(i).lt.xf) then 
        nbinmass(k)=nbinmass(k)+1
        goto 41
      else
        goto 40
      endif

C     ---   Calculating the luminosity function--- 

C     QUESTION: what does it mean?
C     ---   Calculos de turno varios  ---
 41   cont=cont+1.0
      vmedia=vmedia+vo*voinv
      zmedia=zmedia+dabs(coordinate_Zcylindr(i))
      j=0
      mbol=2.5*leb(i) + 4.75 
4     j=j+1

C     ---   Calculating luminosity function of the WD's---
      if (mbol.le.mbolmin+mbolinc*dfloat(j).and.mbol.ge.mbolmin) then
          ndfa(j)=ndfa(j)+1
          nbin(j)=nbin(j)+1
          mbin(j)=mbin(j)+meb(i)
C         calculating sum of velocities of WD in bin Nºj (only from 
C         restricted sample). We will need it for calculating average 
C         velocities of WD for each bin (only from restricted sample)
          sumOfWDVelocitiesInBin_u(j)=sumOfWDVelocitiesInBin_u(j)+uu(i)
          sumOfWDVelocitiesInBin_v(j)=sumOfWDVelocitiesInBin_v(j)+vv(i)
          sumOfWDVelocitiesInBin_w(j)=sumOfWDVelocitiesInBin_w(j)+ww(i)
C         filling arrays of velocites for calculating SD
          arrayOfVelocitiesForSD_u(j,nbin(j))=uu(i)
          arrayOfVelocitiesForSD_v(j,nbin(j))=vv(i)
          arrayOfVelocitiesForSD_w(j,nbin(j))=ww(i)
C         filling array of bolometric magnitudes for each WD in restr.s.
          arrayOfMagnitudes(j,nbin(j))=mbol
      else 
        if (j.eq.nbins) then
          goto 5
        endif
        goto 4
      endif
5     continue

C     NOTE that next loops can be put in one
C     NOTE I need to make subroutines and not to mix all this
      do 50 j=1,nbins
C       calculating average velocities of WD for each bin (only from 
C       restricted sample) 
        averageWDVelocityInBin_u(j)=sumOfWDVelocitiesInBin_u(j)/nbin(j)
        averageWDVelocityInBin_v(j)=sumOfWDVelocitiesInBin_v(j)/nbin(j)
        averageWDVelocityInBin_w(j)=sumOfWDVelocitiesInBin_w(j)/nbin(j)
C       calculating Standart Deviation for velocities in each bin
C       TODO: place all this code for SD in subroutine
        sumOfSquareDifferences_u = 0.0
        sumOfSquareDifferences_v = 0.0
        sumOfSquareDifferences_w = 0.0
        do 51 i=1,nbin(j)
          sumOfSquareDifferences_u=sumOfSquareDifferences_u+
     &                             (arrayOfVelocitiesForSD_u(j,i)-
     &                             averageWDVelocityInBin_u(j))**2
          sumOfSquareDifferences_v=sumOfSquareDifferences_v+
     &                             (arrayOfVelocitiesForSD_v(j,i)-
     &                             averageWDVelocityInBin_v(j))**2
          sumOfSquareDifferences_w=sumOfSquareDifferences_w+
     &                             (arrayOfVelocitiesForSD_w(j,i)-
     &                             averageWDVelocityInBin_w(j))**2
          write(815,444) arrayOfVelocitiesForSD_u(j,i), 
     &                     arrayOfVelocitiesForSD_v(j,i), 
     &                     arrayOfVelocitiesForSD_w(j,i), 
     &                     arrayOfMagnitudes(j,i)
51      continue
        standardDeviation_u(j)=(sumOfSquareDifferences_u/
     &                         dfloat(nbin(j)))**0.5
        standardDeviation_v(j)=(sumOfSquareDifferences_v/
     &                         dfloat(nbin(j)))**0.5
        standardDeviation_w(j)=(sumOfSquareDifferences_w/
     &                         dfloat(nbin(j)))**0.5
50    continue



C-------------------------------------------------------------------
C     ---  Write data of the LF of the WD's
C-------------------------------------------------------------------
      dmonte=0.0

      do 6 i=1,nbins
        if (ndfa(i).le.0.0) then
          ndfa(i)=10.0d-40
        else
          dmonte=dmonte+ndfa(i)*mbolinc
        endif 
6     continue

C     --- Volume of "North Hemisphere" in 40 pc ----
C          V_NH(40 pc)=134041.29
C     normalizing to the bins n=16+17+18, total 220 objects
C       n=17 is lum=-3.8 Mbol=14.25 with 72 objetos       
      fnor=(ndfa(16)+ndfa(17)+ndfa(18))/220.0
      fnora=(134041.29*fnor) 
      write (6,*) 'Factor de normalización:', fnor

C     ---   Recalculating LF   ---
C     QUESTION: what does it mean?
C     ojo factor norma
      do 71 i=1,nbins
        ndfa(i)=ndfa(i)/fnora
71    continue

      do 7 i=1,nbins
C       ---   Calculating error bars final touches  ---
C       old definition of binning
        x=mbolmin+(mbolinc)*dfloat(i)
C       new definition of binning -- when changing, go to BINNING WDLF 
C       and check the values
C        x=mbolmin+(mbolinc)*dfloat(i)-mbolinc/2.0
C       QUESTION: Why is this line here?
        xx=(x-4.75)/2.5
      
        if (nbin(i) .eq. 0) then
          xya=0.0d0
          errsupa=0.0d0
          errinfa=0.0d0
          vvv=0.0
          go to 9
        end if 
      
        xya=dlog10(ndfa(i))
        errsupa=dlog10(ndfa(i)+errora(i))-xya
        if (nbin(i).eq.1) then
          errinfa=-25.0
        else
          errinfa=dlog10(ndfa(i)-errora(i))-xya
        endif

        vvv=0.000
        mbin(i)=mbin(i)/dfloat(nbin(i))
      
C       NOTE I added output for average velocities + SD in each bin here      
9       write(155,200) vvv,xx,xya,errsupa,errinfa,nbin(i),i,
     &                 averageWDVelocityInBin_u(i),
     &                 averageWDVelocityInBin_v(i),
     &                 averageWDVelocityInBin_w(i),
     &                 standardDeviation_u(i),standardDeviation_v(i),
     &                 standardDeviation_w(i)
        write(161,*) xx,mbin(i),nbin(i)

7     continue     

C-----------------------------------------------------------------------
C     --- Writing data of histogram of masses ---
      ntotmass=0
      do 80 i=1,26
        ntotmass=ntotmass+nbinmass(i)
 80   continue

      do 8 iii=1,26
        xxb=dfloat(ntotmass)
        write(162 ,*) xmasi+xmasinc*(dfloat(iii)-0.5),nbinmass(iii)/xxb
 8    continue

C-----------------------------------------------------------------------
C     ---  Reading the LF teoretical/observational for performing chi² 
C          test ---
C     LF global bins i=1,19; k=i+3
C     LF fit cutoff, 3 last bins
C     LF hot, 11 bins hot ones

C     LF global 
      do 10 i=1,19
        k=i+3
        xfl(i)=ndfa(k)*134041.29
        read (71,*) a1,i2
        xflo(i)=dfloat(i2)
        write (6,*) i,xfl(i),xflo(i)
  10  continue
       
      xflcut(1)=xfl(17)
      xflcut(2)=xfl(18)
      xflcut(3)=xfl(19)
      xflocut(1)=xflo(17)
      xflocut(2)=xflo(18)
      xflocut(3)=xflo(19)

      do 12 i=1,11
        xflhot(i)=xfl(i)
        xflohot(i)=xflo(i)
  12  continue

C     choosing region of maximum for calculating chi^2
      do 13 i=1,6
        xflMaxRegion(i)=xfl(i+11)
        xfloMaxRegion(i)=xflo(i+11)
  13  continue

      call chstwo(xfl,xflo,19,0,df,chsq,prob)
      write (6,*) '----------- Chi2 global LF-------------'
      write(6,*) 'df=',df
      write(6,*) 'chsq=',chsq
      write(6,*) 'prob=',prob
      write (6,*) '----------------------------------'

      OPEN (UNIT=20, STATUS='OLD',ACCESS = 'APPEND',
     &     file='./output_data/chisquare_test.out')
      write (20,222) fractionOfDB,galacticDiskAge,parameterIMF,
     &               parameterIFMR,timeOfBurst,df,chsq,prob
      CLOSE (22)

      call chstwo(xflhot,xflohot,11,0,df,chsq,prob)
      write (6,*) '----------- Chi2 HOT branch-------------'
      write(6,*) 'df=',df
      write(6,*) 'chsq=',chsq
      write(6,*) 'prob=',prob
      write (6,*) '----------------------------------'

      call chstwo(xflcut,xflocut,3,0,df,chsq,prob)
      write (6,*) '----------- Chi2 cut-off----------------'
      write(6,*) 'df=',df
      write(6,*) 'chsq=',chsq
      write(6,*) 'prob=',prob
      write (6,*) '----------------------------------'

      call chstwo(xflMaxRegion,xfloMaxRegion,6,0,df,chsq,prob)
      write (6,*) '----------- Chi2 maximum-region---------'
      write(6,*) 'df=',df
      write(6,*) 'chsq=',chsq
      write(6,*) 'prob=',prob
      write (6,*) '----------------------------------'
C     output for chi² of maximum-region vs galactic disk age - test
      OPEN (UNIT=30, STATUS='OLD',ACCESS = 'APPEND',
     &     file='./output_data/chi_sq_test_maxregion.out')
      write (30,333) galacticDiskAge,chsq,prob
      CLOSE (30)

 222   format(6(f6.3,1x),f6.2,1x,f7.4)
 333   format(f4.1,2x,f8.4,2x,f8.4)
 444   format(3(f7.2,2x),f6.3)

C-----------------------------------------------------------------------
C     ----- Some results of the sample ----
             
      write(6,*) 'Initial number of WDs:               ',ntwd
      write(6,*) 'Eliminated by parallax:              ',r1
      write(6,*) '    "       "       declination:     ',r2
      write(6,*) 'Initial number northern hemisphere:  ',ntwd-r1-r2
      write(6,*) '    "       " proper motion:         ',r3
      write(6,*) '    "       " reduced proper motion: ',r32
      write(6,*) '    "       " apparent magnitude:    ',r4
      write(6,*) 'Restricted sample           :        ',ntwd-r1-r2-r3-
     &           r32-r4

200   format(f6.3,2x,f6.3,2x,3(1pd14.7,2x),i4,i4,2x,6(1pd14.7,2x))

      return
      end
C***********************************************************************
 


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


C***********************************************************************
C           FUNCTION RAND
C***********************************************************************
        FUNCTION RAN(IDUMMY)
C
C     THIS IS AN ADAPTED VERSION OF SUBROUTINE RANECU WRITTEN BY
C     F. JAMES (COMPUT. PHYS. COMMUN. 60 (1990) 329-344, WHICH HAS
C     BEEN MODIFIED TO GIVE A SINGLE RANDOM NUMBER AT EACH CALL.
C     THE 'SEEDS' ISEED1 AND ISEED2 MUST BE INITIALIZED IN THE 
C     MAIN PROGRAM AND TRANSFERRED THROUGH THE NAMED COMMON BLOCK
C     /RSEED/.
        
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      real ran
      PARAMETER (USCALE=1.0D0/2.0D0**31)
      COMMON /RSEED/ ISEED1,ISEED2

      I1=ISEED1/53668
      ISEED1=40014*(ISEED1-I1*53668)-I1*12211
      
      IF (ISEED1.LT.0) ISEED2=ISEED1+2147483563
      
      I2=ISEED2/52774
      ISEED2=40692*(ISEED2-I2*52774)-I2*3791
      
      IF (ISEED2.LT.0) ISEED2=ISEED2+2147483399
      
      IZ=ISEED1-ISEED2
      
      IF(IZ.LT.1) IZ=IZ+2147483562
    
      RAN=IZ*USCALE

       
      RETURN
      END 
C***********************************************************************
        

            
C***********************************************************************
C     SUBROUTINE ODEINT
C     Subroutine for the numerical integration of a system of first  
C     order differential equations
C     ("Numerical recipes in Fortran", Willian H. Press)
C***********************************************************************      

      SUBROUTINE ODEINT(YSTART,NVAR,X1,X2,EPS,H1,HMIN,NOK,
     &           NBAD,DERIVS,RKQC,yscal,y,dydx)
      implicit double precision (a-h,o-z)
      PARAMETER (MAXSTP=10000,NMAX=10,TWO=2.0,ZERO=0.0,TINY=1.E-30)
      COMMON /PATH/ KMAX,KOUNT,DXSAV,XP(200),YP(10,200)
      DIMENSION YSTART(NVAR),YSCAL(Nvar),Y(Nvar),DYDX(Nvar)
      EXTERNAL DERIVS
      EXTERNAL RKQC
      X=X1
      H=dSIGN(H1,X2-X1)
      NOK=0
      NBAD=0
      KOUNT=0
      DO 11 I=1,NVAR
        Y(I)=YSTART(I)
11    CONTINUE
      IF (KMAX.GT.0) XSAV=X-DXSAV*TWO
        DO 16 NSTP=1,MAXSTP
          CALL DERIVS(X,Y,DYDX,xpla,ypla)
          DO 12 I=1,NVAR
            YSCAL(I)=DABS(Y(I))+DABS(H*DYDX(I))+TINY
12        CONTINUE
          IF(KMAX.GT.0)THEN
            IF(DABS(X-XSAV).GT.DABS(DXSAV)) THEN
              IF(KOUNT.LT.KMAX-1)THEN
                KOUNT=KOUNT+1
                XP(KOUNT)=X
                DO 13 I=1,NVAR
                  YP(I,KOUNT)=Y(I)
13              CONTINUE
                XSAV=X
              ENDIF
            ENDIF
          ENDIF
          IF((X+H-X2)*(X+H-X1).GT.ZERO) H=X2-X
          CALL RKQC(Y,DYDX,NVAR,X,H,EPS,YSCAL,HDID,HNEXT,DERIVS)
          IF(HDID.EQ.H)THEN
            NOK=NOK+1
          ELSE
            NBAD=NBAD+1
          ENDIF
          IF((X-X2)*(X2-X1).GE.ZERO)THEN
            DO 14 I=1,NVAR
              YSTART(I)=Y(I)
14          CONTINUE
            IF(KMAX.NE.0)THEN
              KOUNT=KOUNT+1
              XP(KOUNT)=X
              DO 15 I=1,NVAR
                YP(I,KOUNT)=Y(I)
15            CONTINUE
            ENDIF
            RETURN
          ENDIF
          IF(dABS(HNEXT).LT.HMIN) stop 'Stepsize smaller than minimum.'
            H=HNEXT
16      CONTINUE
        STOP 'Too many steps.'

      END
C***********************************************************************

  
C***********************************************************************
C     SUBROUTINE RKQC
C     Runge-Kutta Quality Control
C     Compares the results between one big step and two small steps and
C     decides the size of the next step according to how different the
C     results are
C     ("Numerical recipes in Fortran", Willian H. Press)
C***********************************************************************
      SUBROUTINE RKQC(Y,DYDX,N,X,HTRY,EPS,YSCAL,HDID,HNEXT,DERIVS)
      implicit double precision(a-h,o-z)
      PARAMETER (NMAX=10,FCOR=.0666666667,ONE=1.,SAFETY=0.9,
     &          ERRCON=6.E-4)
      EXTERNAL DERIVS
      DIMENSION Y(N),DYDX(N),YSCAL(N),YTEMP(NMAX),YSAV(NMAX),DYSAV(NMAX)
      PGROW=-0.20
      PSHRNK=-0.25
      XSAV=X
      DO 11 I=1,N
        YSAV(I)=Y(I)
        DYSAV(I)=DYDX(I)
11    CONTINUE
      H=HTRY
1     HH=0.5*H
      CALL RK4(YSAV,DYSAV,N,XSAV,HH,YTEMP,DERIVS)
      X=XSAV+HH
      CALL DERIVS(X,YTEMP,DYDX)
      CALL RK4(YTEMP,DYDX,N,X,HH,Y,DERIVS)
      X=XSAV+H
      IF(X.EQ.XSAV) stop 'Stepsize not significant in RKQC.'
      CALL RK4(YSAV,DYSAV,N,XSAV,H,YTEMP,DERIVS)
      ERRMAX=0.
      DO 12 I=1,N
        YTEMP(I)=Y(I)-YTEMP(I)
        ERRMAX=dMAX1(ERRMAX,DABS(YTEMP(I)/YSCAL(I)))
12    CONTINUE
      ERRMAX=ERRMAX/EPS
      IF(ERRMAX.GT.ONE) THEN
        H=SAFETY*H*(ERRMAX**PSHRNK)
        GOTO 1
      ELSE
        HDID=H
        IF(ERRMAX.GT.ERRCON)THEN
          HNEXT=SAFETY*H*(ERRMAX**PGROW)
        ELSE
          HNEXT=4.*H
        ENDIF
      ENDIF
      DO 13 I=1,N
        Y(I)=Y(I)+YTEMP(I)*FCOR
13    CONTINUE
      
      RETURN
      END
C***********************************************************************


C***********************************************************************
C     SUBROUTINE RK4
C     Runge-Kutta 4
C     Runge-Kutta step on a set of n differential equations. You input  
C     the values of the independent variables, and you get out new  
C     values which are stepped by a stepsize h (which can be positive or 
C     negative). The routine that has three calls to derivs
C     ("Numerical recipes in Fortran", Willian H. Press)
C***********************************************************************        
      SUBROUTINE RK4(Y,DYDX,N,X,H,YOUT,DERIVS)
      implicit double precision (a-h,o-z)
      
      PARAMETER (NMAX=10)
      DIMENSION Y(N),DYDX(N),YOUT(N),YT(NMAX),DYT(NMAX),DYM(NMAX)
      EXTERNAL DERIVS
      HH=H*0.5
      H6=H/6.
      XH=X+HH
      DO 11 I=1,N
        YT(I)=Y(I)+HH*DYDX(I)
11    CONTINUE
      CALL DERIVS(XH,YT,DYT)
      DO 12 I=1,N
        YT(I)=Y(I)+HH*DYT(I)
12    CONTINUE
      CALL DERIVS(XH,YT,DYM)
      DO 13 I=1,N
        YT(I)=Y(I)+H*DYM(I)
        DYM(I)=DYT(I)+DYM(I)
13    CONTINUE
      CALL DERIVS(X+H,YT,DYT)
      DO 14 I=1,N
        YOUT(I)=Y(I)+H6*(DYDX(I)+DYT(I)+2.*DYM(I))
14    CONTINUE
      RETURN
      END
C***********************************************************************        


C***********************************************************************
C     SUBROUTINE derivs
C     Sample derivatives routine for stiff 
C     ("Numerical recipes in Fortran", Willian H. Press)
C***********************************************************************              
      subroutine derivs(x,y,dydx)
C=======================================================================
C     This subroutine provides the values of the derivatives of the 
C     variables y_i
C-------------------------------------------------------------------
      implicit double precision (a-h,o-z)

C     ---   Dimensions   ---      
      dimension y(2),dydx(2)

C     ---   Calculating the derivatives    ---
      dydx(1)=y(2)
      call fuerza(y(1),f)
      dydx(2)=f
      
      return
      end     
C***********************************************************************   

      SUBROUTINE CHSTWO(BINS1,BINS2,NBINS,KNSTRN,DF,CHSQ,PROB)
      implicit double precision(a-h,o-z)
      integer knstrn
      DIMENSION BINS1(NBINS),BINS2(NBINS)
      DF=dfloat(NBINS-1-KNSTRN)
      CHSQ=0.
      DO 11 J=1,NBINS
        IF(BINS1(J).EQ.0..AND.BINS2(J).EQ.0.)THEN
          DF=DF-1.
        ELSE
          CHSQ=CHSQ+(BINS1(J)-BINS2(J))**2/(BINS1(J)+BINS2(J))
        ENDIF
11    CONTINUE
      PROB=GAMMQ(0.5*DF,0.5*CHSQ)
      RETURN
      END

      FUNCTION GAMMQ(A,X)
      implicit double precision (a-h,o-z)
      IF(X.LT.0..OR.A.LE.0.)PAUSE
      IF(X.LT.A+1.)THEN
        CALL GSER(GAMSER,A,X,GLN)
        GAMMQ=1.-GAMSER
      ELSE
        CALL GCF(GAMMQ,A,X,GLN)
      ENDIF
      RETURN
      END
      
      SUBROUTINE GSER(GAMSER,A,X,GLN)
      implicit double precision(a-h,o-z)
      PARAMETER (ITMAX=100,EPS=3.E-7)
      GLN=GAMMLN(A)
      IF(X.LE.0.)THEN
        IF(X.LT.0.)PAUSE
        GAMSER=0.
        RETURN
      ENDIF
      AP=A
      SUM=1./A
      DEL=SUM
      DO 11 N=1,ITMAX
        AP=AP+1.
        DEL=DEL*X/AP
        SUM=SUM+DEL
        IF(dABS(DEL).LT.dABS(SUM)*EPS)GO TO 1
11    CONTINUE
      PAUSE 'A too large, ITMAX too small'
1     GAMSER=SUM*dEXP(-X+A*dLOG(X)-GLN)
      RETURN
      END

      SUBROUTINE GCF(GAMMCF,A,X,GLN)
      implicit double precision (a-h,o-z)
      PARAMETER (ITMAX=100,EPS=3.E-7)
      GLN=GAMMLN(A)
      GOLD=0.
      A0=1.
      A1=X
      B0=0.
      B1=1.
      FAC=1.
      DO 11 N=1,ITMAX
        AN=dFLOAT(N)
        ANA=AN-A
        A0=(A1+A0*ANA)*FAC
        B0=(B1+B0*ANA)*FAC
        ANF=AN*FAC
        A1=X*A0+ANF*A1
        B1=X*B0+ANF*B1
        IF(A1.NE.0.)THEN
          FAC=1./A1
          G=B1*FAC
          IF(dABS((G-GOLD)/G).LT.EPS)GO TO 1
          GOLD=G
        ENDIF
11    CONTINUE
      PAUSE 'A too large, ITMAX too small'
1     GAMMCF=dEXP(-X+A*dLOG(X)-GLN)*G
      RETURN
      END


      FUNCTION GAMMLN(XX)
      implicit double precision (a-h,o-z)
      REAL*8 COF(6),STP,HALF,ONE,FPF,X,TMP,SER
      DATA COF,STP/76.18009173D0,-86.50532033D0,24.01409822D0,
     *    -1.231739516D0,.120858003D-2,-.536382D-5,2.50662827465D0/
      DATA HALF,ONE,FPF/0.5D0,1.0D0,5.5D0/
      X=XX-ONE
      TMP=X+FPF
      TMP=(X+HALF)*dLOG(TMP)-TMP
      SER=ONE
      DO 11 J=1,6
        X=X+ONE
        SER=SER+COF(J)/X
11    CONTINUE
      GAMMLN=TMP+dLOG(STP*SER)
      RETURN
      END