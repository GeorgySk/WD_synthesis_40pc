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

      integer i,ISEED1,ISEED2,iseed,numberOfStarsInSample
      double precision randomNumber,fractionOfDB

C     QUESTION: what is gamma?    
      double precision parameterIFMR,variationOfGravConst,gamma

C     TODO: make only 11 instances for each type not for each file      
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
      
      common /vargra/ variationOfGravConst,gamma
      common /param/ fractionOfDB,galacticDiskAge,parameterIMF,
     &               parameterIFMR,timeOfBurst
      common /tables/ table
      
      include 'code/tables_linking.f'

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
     &     table(1)%ncol,table(1)%mass,
     &     table(1)%coolingTime,table(1)%prevTime,table(1)%luminosity,
     &     table(1)%effTemp,table(1)%gravAcc)
      call incoolda(table(8)%flag,table(8)%initLink,table(8)%ntrk,
     &     table(8)%ncol,table(8)%mass,table(8)%coolingTime,
     &     table(8)%prevTime,table(8)%luminosity,table(8)%effTemp,
     &     table(8)%gravAcc)
      call incoolda(table(18)%flag,table(18)%initLink,table(18)%ntrk,
     &     table(18)%ncol,table(18)%mass,table(18)%coolingTime,
     &     table(18)%prevTime,table(18)%luminosity,
     &     table(18)%effTemp,
     &     table(18)%gravAcc)
      call incoolda(table(26)%flag,table(26)%initLink,table(26)%ntrk,
     &     table(26)%ncol,table(26)%mass,table(26)%coolingTime,
     &     table(26)%prevTime,table(26)%luminosity,
     &     table(26)%effTemp,
     &     table(26)%gravAcc)
      
      write(6,*) '   1.2 Tracks of CO non-DA (DB) WD'

C     TODO: rename the function 'incooldb'
      call incooldb(table(34)%flag,table(34)%initLink,
     &     table(34)%ncol,table(34)%ntrk,
     &     table(34)%mass,table(34)%coolingTime,
     &     table(34)%prevTime,table(34)%luminosity,
     &     table(34)%effTemp,table(34)%gravAcc)
      call incooldb(table(41)%flag,table(41)%initLink,
     &     table(41)%ncol,table(41)%ntrk,
     &     table(41)%mass,table(41)%coolingTime,
     &     table(41)%prevTime,table(41)%luminosity,
     &     table(41)%effTemp,table(41)%gravAcc)
      call incooldb(table(50)%flag,table(50)%initLink,
     &     table(50)%ncol,table(50)%ntrk,
     &     table(50)%mass,table(50)%coolingTime,
     &     table(50)%prevTime,table(50)%luminosity,
     &     table(50)%effTemp,table(50)%gravAcc)

      write(6,*) '   1.3 Tracks of ONe DA WD'

C     TODO: rename the function 'incoolone'
      call incoolone
      
      write(6,*) '   1.4 Reading the colors table of Rene(DAs) and Berge
     &ron(DBs)'
C     TODO: rename these functions      
      call color(table(77)%ncol,table(77)%ntrk,table(77)%mass,
     &     table(77)%luminosity,
     &     table(77)%color_U,table(77)%color_B,table(77)%color_R,
     &     table(77)%color_V,table(77)%color_I)      
      call colordb(table(70)%ncol,table(70)%ntrk,
     &     table(70)%mass,table(70)%luminosity,table(70)%color_U,
     &     table(70)%color_B,table(70)%color_V,table(70)%color_R,
     &     table(70)%color_I)

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