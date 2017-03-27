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

C     NOTE: modules are better than includes

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

      integer i,j,k,ISEED1,ISEED2,iseed,numberOfStarsInSample
      double precision randomNumber,fractionOfDB
      double precision parameterIFMR

      TYPE(FileGroupInfo),DIMENSION(11) :: table

      common /RSEED/ ISEED1,ISEED2
      common /param/ fractionOfDB,galacticDiskAge,parameterIMF,
     &               parameterIFMR,timeOfBurst

C     filling info about groups of files (cooling, color tables)      
      include 'code/tables_linking.f'
   

C     ---  Reading free parameters  ---
C ======================================================================
C      galacticDiskAge (Gyr)
C      parameterIMF (alpha): M^{alpha}
C      Initial-to-Final Mass Relation (IFMR) : 
C         mfinal_new=parameterIFMR*mfinal_old
        
C     Reading parameters line from $temporary_files/grid_set_line.in:
      read (10,*) fractionOfDB,galacticDiskAge,parameterIMF,
     &            parameterIFMR,timeOfBurst

C       Fiducial values (trusted):
C           fractionOfDB=0.20
C           galacticDiskAge=8.9
C           parameterIMF=-2.35
C           parameterIFMR=1.0
C           timeOfBurst=0.6

C     Overwriting parameters (this is not good)
      fractionOfDB=0.20 
      galacticDiskAge=8.9
      parameterIMF=-2.35
      parameterIFMR=1.0
      timeOfBurst=0.6     

      write(6,*) '=========================================='
      write(6,*) ' '
      write(6,*) '            Programa monte.f'
      write(6,*) '          by S.Torres, 14.02.11 '
      write(6,*) ' '
      write(6,*) '            Used parameters:'
      write(6,*) 'numberOfStars=    ',numberOfStars
      write(6,*) 'SFR: parameterOfSFR=',parameterOfSFR,'Gyr'
      write(6,*) 'galacticDiskAge=    ',galacticDiskAge,'Gyr'
      write(6,*) 'minimumSectorRadius=',minimumSectorRadius,'kpc' 
      write(6,*) 'maximumSectorRadius=',maximumSectorRadius,'kpc'
      write(6,*) 'radiusOfSector=     ',radiusOfSector,'kpc'
      write(6,*) ' '
      write(6,*) '=========================================='
      write(6,*) ' '
      write(6,*) '          Start of calculations:'
      write(6,*) ' '
      write(6,*) 'Initializing random number generator and reading the s
     &eeds'
      write(6,*) ' '

      iseed=-9
C     Reading line from $temporary_files/seeds_line.in
      read(72,100) iseed1,iseed2
100   format(I6,2x,I6)  
      write(6,*) 'iseed1=',iseed1
      write(6,*) 'iseed2=',iseed2

C     QUESTION: why do we need this part?      
      do 8123 i=1,10
        randomNumber=ran(iseed)
        write (6,*) i,randomNumber
8123  continue  

      write (157,157) fractionOfDB,galacticDiskAge,parameterIMF,
     &                parameterIFMR,timeOfBurst
 157  format(5(f6.3,2x))



C     ---  Calculating the area of the sector  ---
C ======================================================================
C     This style ensures maximum precision when assigning a value to PI.
      pi=4.0*atan(1.0)
C     QUESTION: what about square function?
      areaOfSector=pi*radiusOfSector**2


C     ---  Program itself  ---
C ======================================================================
      write(6,*) '1. Reading the cooling tables (1/10)'

      write(6,*) '   1.1 Tracks of CO DA WD Z=0.001;0.01;0.03;0.06'
C     Calling the function 'incoolda' for 4 metalicities that we have
      call incoolda(table(1))
      call incoolda(table(2))
      call incoolda(table(3))
      call incoolda(table(4))
      
      write(6,*) '   1.2 Tracks of CO non-DA (DB) WD'
C     Calling the function 'incooldb' for 3 metalicities that we have
      call incooldb(table(5))
      call incooldb(table(6))
      call incooldb(table(7))

      write(6,*) '   1.3 Tracks of ONe DA WD'
      call incoolone
      
      write(6,*) '   1.4 Reading the colors table of Rene(DAs) and Berge
     &ron(DBs)'
      call color(table(11))
      call colordb2(table(10))

      write (6,*) '   1.5 Reading the tables of CO DA with G variable'      
      call incoolea

      write(6,*) '2. Calling the IMF and SFR (2/10)'
      call gen(iseed,parameterOfSFR,areaOfSector,numberOfStarsInSample,
     &     galacticDiskAge,timeOfBurst)
      write(6,*) "numberOfStarsInSample=", numberOfStarsInSample
      
      write(6,*) '3. Calculating luminosities (3/10)'
      call lumx(iseed,numberOfStarsInSample)      

      write(6,*) '4. Calculating polar coordinates (4/10)'
      call polar(iseed,minimumSectorRadius,maximumSectorRadius,
     &     angleCoveringSector,radiusOfSector,
     &     solarGalactocentricDistance,scaleLength)

      write(6,*) '5. Generating heliocentric velocities (5/10)'
      call velh(iseed,numberOfStarsInSample)

C     QUESTION: why are we missing the next step?
      goto 7
C     QUESTION: what does this mean?
C     ---   Calculating the trajectories according to/along z-coordinate ---
      write(6,*) '6. Integrating trajectories (6/10)'
      call traject(galacticDiskAge)

7     write(6,*) '7. Calculating coordinates (7/10)'
      call coor(solarGalactocentricDistance)

      write(6,*) '8. Determinating visual magnitudes (8/10)'
      call magi(fractionOfDB,table) 

C     TODO: give a better description to this step
      write(6,*) '9. Generating the Luminosity Function (9/10)'
      call volum_40pc
      
C     NOTE: here it is useless
      write(6, *) '10. Making vrad to be null (10/10)'
C      call vrado

      write (6,*) 'The end'
 

      stop
      end
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