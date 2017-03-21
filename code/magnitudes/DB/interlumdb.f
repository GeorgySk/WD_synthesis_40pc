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
      use external_types
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

      TYPE(FileInfo),DIMENSION(86) :: table

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

      common /tables/ table
     
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
     &       numberOfSequencesInGroup_1,table(34)%ntrk,
     &       table(34)%coolingTime,table(34)%prevTime,
     &       table(34)%mass,table(34)%luminosity,lum1)
        modlog=1
        call interp(model,modlog,tcool,mass,
     &       numberOfSequencesInGroup_1,table(34)%ntrk,
     &       table(34)%coolingTime,table(34)%prevTime,
     &       table(34)%mass,table(34)%effTemp,teff1)
        modlog=0
        call interp(model,modlog,tcool,mass,
     &       numberOfSequencesInGroup_1,table(34)%ntrk,
     &       table(34)%coolingTime,table(34)%prevTime,
     &       table(34)%mass,table(34)%gravAcc,logg1)

        Z2=zet2
        call interp(model,modlog,tcool,mass,
     &       numberOfSequencesInGroup_2,table(41)%ntrk,
     &       table(41)%coolingTime,table(41)%prevTime,
     &       table(41)%mass,table(41)%luminosity,lum2)
        modlog=1
        call interp(model,modlog,tcool,mass,
     &       numberOfSequencesInGroup_2,table(41)%ntrk,
     &       table(41)%coolingTime,table(41)%prevTime,
     &       table(41)%mass,table(41)%effTemp,teff2)
        modlog=0
        call interp(model,modlog,tcool,mass,
     &       numberOfSequencesInGroup_2,table(41)%ntrk,
     &       table(41)%coolingTime,table(41)%prevTime,
     &       table(41)%mass,table(41)%gravAcc,logg2)
      end if
      
      if(Z.ge.zet2.AND.Z.lt.zet3) then
        Z1=zet2
        call interp(model,modlog,tcool,mass,
     &       numberOfSequencesInGroup_2,table(41)%ntrk,
     &       table(41)%coolingTime,table(41)%prevTime,
     &       table(41)%mass,table(41)%luminosity,lum1)
        modlog=1
        call interp(model,modlog,tcool,mass,
     &       numberOfSequencesInGroup_2,table(41)%ntrk,
     &       table(41)%coolingTime,table(41)%prevTime,
     &       table(41)%mass,table(41)%effTemp,teff1)
        modlog=0
        call interp(model,modlog,tcool,mass,
     &       numberOfSequencesInGroup_2,table(41)%ntrk,
     &       table(41)%coolingTime,table(41)%prevTime,
     &       table(41)%mass,table(41)%gravAcc,logg1)
      Z2=zet3
        call interp(model,modlog,tcool,mass,numberOfSequencesInGroup_3,
     &       table(50)%ntrk,table(50)%coolingTime,
     &       table(50)%prevTime,table(50)%mass,
     &       table(50)%luminosity,lum2)
        modlog=1
        call interp(model,modlog,tcool,mass,numberOfSequencesInGroup_3,
     &       table(50)%ntrk,table(50)%coolingTime,
     &       table(50)%prevTime,table(50)%mass,
     &       table(50)%effTemp,teff2)
        modlog=0
        call interp(model,modlog,tcool,mass,numberOfSequencesInGroup_3,
     &       table(50)%ntrk,table(50)%coolingTime,
     &       table(50)%prevTime,table(50)%mass,table(50)%gravAcc,
     &       logg2)
      end if
      
C     interpolation 3D in function of Z
      lum=lum1+(lum2-lum1)*(Z-Z1)/(Z2-Z1)
      teff=teff1+(teff2-teff1)*(Z-Z1)/(Z2-Z1)
      logg=logg1+(logg2-logg1)*(Z-Z1)/(Z2-Z1)
      
C     colors interpolation, here there is no metallicity
      call intermag(mass,lum,numberOfSequences,
     &     table(70)%ntrk,table(70)%mass,table(70)%luminosity,
     &     table(70)%color_U,table(70)%color_B,table(70)%color_V,
     &     table(70)%color_R,table(70)%color_I,c1,c2,c3,c4,c5)

      return
      end
C***********************************************************************