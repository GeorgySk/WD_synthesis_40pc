#!/bin/ksh
#This is a reformatted version of disk40 by Georgy Skorobogatov


totalNumberOfExecutions=1
#next variable can be used to continue previous series of simulations
#all new output files will have indexes beginning with this: _0001.out, _0002.out etc.
initialIndexForOutputs=1

#making links to folders
DA_Tracks=./seclea
DB_Tracks_z0001=./DB_z=0.001
DB_Tracks_z001=./DB_z=0.01
DB_Tracks_z006=./DB_z=0.06
ONe_Tracks=./data
DB_Colors_byBergeron=./bergeron
Colors_byRenedo=./renedo
Results=./result
BASTI=./basti
G_Variable=./tracks

#TODO what does 'rm' do?
'rm' chidos.out
touch chidos.out
mkdir -p $Results/fl-althaus/
mkdir -p ./bootstrap/R-althaus/



currentNumberOfExecution=1
currentIndexForOutput=$initialIndexForOutputs

while test $currentNumberOfExecution -le $totalNumberOfExecutions
do
   #TODO what does next 4 lines mean?
   'rm' seeds.in
   'rm' parameter.in
   tail -$currentNumberOfExecution seeds.dat  > seeds.in
   tail -$currentNumberOfExecution parameter.dat > parameter.in

   echo "Current execution: " $currentNumberOfExecution "from total of" $totalNumberOfExecutions
   echo "Simulation Nº" $currentIndexForOutput

   #TODO Why is the next line is here? Where does the program runs?
   gfortran disk40.f -o disk40.e

   
   #Parameters of the simulations 
   #TODO I need to place all these lines somewhere. + where to put $currentIndexForOutput? + I should relocate folders

   #Various input parameters
   ln -s ./parameter.in fort.10

   #TRACKS WD CO DA,  Z=0.001
   ln -s $DA_Tracks/wdtracks_z0001/wd0505_z0001.trk fort.11
   ln -s $DA_Tracks/wdtracks_z0001/wd0553_z0001.trk fort.12
   ln -s $DA_Tracks/wdtracks_z0001/wd0593_z0001.trk fort.13
   ln -s $DA_Tracks/wdtracks_z0001/wd0627_z0001.trk fort.14
   ln -s $DA_Tracks/wdtracks_z0001/wd0660_z0001.trk fort.15
   ln -s $DA_Tracks/wdtracks_z0001/wd0692_z0001.trk fort.16
   ln -s $DA_Tracks/wdtracks_z0001/wd0863_z0001.trk fort.17
   #TRACKS WD CO DA,  Z=0.01
   ln -s $DA_Tracks/wdtracks_z001/wd0524_z001.trk fort.21
   ln -s $DA_Tracks/wdtracks_z001/wd0570_z001.trk fort.22
   ln -s $DA_Tracks/wdtracks_z001/wd0593_z001.trk fort.23
   ln -s $DA_Tracks/wdtracks_z001/wd0609_z001.trk fort.24
   ln -s $DA_Tracks/wdtracks_z001/wd0632_z001.trk fort.25
   ln -s $DA_Tracks/wdtracks_z001/wd0659_z001.trk fort.26
   ln -s $DA_Tracks/wdtracks_z001/wd0705_z001.trk fort.27
   ln -s $DA_Tracks/wdtracks_z001/wd0767_z001.trk fort.28
   ln -s $DA_Tracks/wdtracks_z001/wd0837_z001.trk fort.29
   ln -s $DA_Tracks/wdtracks_z001/wd0877_z001.trk fort.30
   #TRACKS WD CO DA,  Z=0.0382
   ln -s $DA_Tracks/diff_003/0524_003_sflhdiff.trk fort.31
   ln -s $DA_Tracks/diff_003/0570_003_sflhdiff.trk fort.32
   ln -s $DA_Tracks/diff_003/0593_003_sflhdiff.trk fort.33
   ln -s $DA_Tracks/diff_003/0610_003_sflhdiff.trk fort.34
   ln -s $DA_Tracks/diff_003/0632_003_sflhdiff.trk fort.35
   ln -s $DA_Tracks/diff_003/0659_003_sflhdiff.trk fort.36
   ln -s $DA_Tracks/diff_003/0705_003_sflhdiff.trk fort.37
   ln -s $DA_Tracks/diff_003/1000_003_sflhdiff.trk fort.38
   #TRACKS WD CO DA,  Z=0.06
   ln -s $DA_Tracks/diff_006/0524_006_sflhdiff.trk fort.41
   ln -s $DA_Tracks/diff_006/0570_006_sflhdiff.trk fort.42
   ln -s $DA_Tracks/diff_006/0593_006_sflhdiff.trk fort.43
   ln -s $DA_Tracks/diff_006/0610_006_sflhdiff.trk fort.44
   ln -s $DA_Tracks/diff_006/0632_006_sflhdiff.trk fort.45
   ln -s $DA_Tracks/diff_006/0659_006_sflhdiff.trk fort.46
   ln -s $DA_Tracks/diff_006/0705_006_sflhdiff.trk fort.47
   ln -s $DA_Tracks/diff_006/1000_006_sflhdiff.trk fort.48

   #TRACKS WD CO DB,  Z=0.001
   ln -s $DB_Tracks_z0001/05047_db_Z=0.001.trk fort.91
   ln -s $DB_Tracks_z0001/05527_db_Z=0.001.trk fort.92
   ln -s $DB_Tracks_z0001/059328_db_Z=0.001.trk fort.93
   ln -s $DB_Tracks_z0001/062738_db_Z=0.001.trk fort.94
   ln -s $DB_Tracks_z0001/06602_db_Z=0.001.trk fort.95
   ln -s $DB_Tracks_z0001/069289_db_Z=0.001.trk fort.96
   ln -s $DB_Tracks_z0001/08637_db_Z=0.001.trk fort.97
   #TRACKS WD CO DB,  Z=0.01
   ln -s $DB_Tracks_z001/db_cool_0514.seq fort.101
   ln -s $DB_Tracks_z001/db_cool_0530.seq fort.102
   ln -s $DB_Tracks_z001/db_cool_0542.seq fort.103
   ln -s $DB_Tracks_z001/db_cool_0565.seq fort.104
   ln -s $DB_Tracks_z001/db_cool_0584.seq fort.105
   ln -s $DB_Tracks_z001/db_cool_0609.seq fort.106
   ln -s $DB_Tracks_z001/db_cool_0664.seq fort.107
   ln -s $DB_Tracks_z001/db_cool_0741.seq fort.108
   ln -s $DB_Tracks_z001/db_cool_0869.seq fort.109
   #TRACKS WD CO DB,  Z=0.06
   ln -s $DB_Tracks_z006/0524db_006_sflhdiff.trk fort.111
   ln -s $DB_Tracks_z006/0570db_006_sflhdiff.trk fort.112
   ln -s $DB_Tracks_z006/0593db_006_sflhdiff.trk fort.113
   ln -s $DB_Tracks_z006/061db_006_sflhdiff.trk fort.114
   ln -s $DB_Tracks_z006/0632db_006_sflhdiff.trk fort.115
   ln -s $DB_Tracks_z006/0659db_006_sflhdiff.trk fort.116
   ln -s $DB_Tracks_z006/070db_006_sflhdiff.trk fort.117
   ln -s $DB_Tracks_z006/076db_006_sflhdiff.trk fort.118
   ln -s $DB_Tracks_z006/087db_006_sflhdiff.trk fort.119

   #TRACKS WD ONe
   ln -s $ONe_Tracks/color_106.out fort.121
   ln -s $ONe_Tracks/color_110.out fort.122
   ln -s $ONe_Tracks/color_116.out fort.123
   ln -s $ONe_Tracks/color_120.out fort.124
   ln -s $ONe_Tracks/color_124.out fort.125
   ln -s $ONe_Tracks/color_128.out fort.126
   ln -s $ONe_Tracks/t106_he.trk fort.127
   ln -s $ONe_Tracks/t110_he.trk fort.128
   ln -s $ONe_Tracks/t116_he.trk fort.129
   ln -s $ONe_Tracks/t120_he.trk fort.130
   ln -s $ONe_Tracks/t128_he.trk fort.131

   #TODO read more about it
   #COLORS WD DB (Bergeron)
   ln -s $DB_Colors_byBergeron/Table_Mass_0.5_DB_sort   fort.132
   ln -s $DB_Colors_byBergeron/Table_Mass_0.6_DB_sort   fort.133
   ln -s $DB_Colors_byBergeron/Table_Mass_0.7_DB_sort   fort.134
   ln -s $DB_Colors_byBergeron/Table_Mass_0.8_DB_sort   fort.135
   ln -s $DB_Colors_byBergeron/Table_Mass_0.9_DB_sort   fort.136
   ln -s $DB_Colors_byBergeron/Table_Mass_1.0_DB_sort   fort.137
   ln -s $DB_Colors_byBergeron/Table_Mass_1.2_DB_sort   fort.138

   #TODO read more about it
   #COLORS (Renedo)
   ln -s $Colors_byRenedo/cox_0524.dat   fort.61
   ln -s $Colors_byRenedo/cox_0570.dat   fort.62
   ln -s $Colors_byRenedo/cox_0593.dat   fort.63
   ln -s $Colors_byRenedo/cox_0609.dat   fort.64
   ln -s $Colors_byRenedo/cox_0632.dat   fort.65
   ln -s $Colors_byRenedo/cox_0659.dat   fort.66
   ln -s $Colors_byRenedo/cox_0705.dat   fort.67
   ln -s $Colors_byRenedo/cox_0767.dat   fort.68
   ln -s $Colors_byRenedo/cox_0837.dat   fort.69
   ln -s $Colors_byRenedo/cox_0877.dat   fort.70

   #TODO why in output?
   #Observational Luminosity Function
   ln -s $Results/fl-althaus/fl_obs_40pc.out fort.71
   
   #SEEDS
   ln -s ./seeds.in fort.72

   #Outputs
   ln -s $Results/fl-althaus/coloresfull.out fort.73
   ln -s $Results/fl-althaus/tburstdata.out  fort.81
   ln -s ./extrapolcolores.out fort.1007
   ln -s ./errors.out fort.1008
   ln -s $Results/fl_prueba.out fort.155
   ln -s ./bootstrap/R-althaus/boot_rowell_thin_$currentNumberOfExecution.out fort.156
   ln -s /home/georgy/Documents/simulation_results/brt/new_file.out fort.1156
   ln -s ./parametersMC.out  fort.157
   #Velocities
   ln -s $Results/fl_velocidades.out  fort.158
   ln -s $Results/fl_velocidades_u.out  fort.400
   ln -s $Results/fl_velocidades_v.out  fort.401
   ln -s $Results/fl_velocidades_w.out  fort.402
   ln -s $Results/fl_velocidadesf.out  fort.159
   ln -s $Results/fl_velocidades_uf.out  fort.403
   ln -s $Results/fl_velocidades_vf.out  fort.404
   ln -s $Results/fl_velocidades_wf.out  fort.405

   ln -s $Results/fl-althaus/mboldis_complete.out fort.160
   ln -s $Results/flmass$currentNumberOfExecution.out fort.161
   ln -s $Results/mass_dis.out fort.162
   ln -s $Results/fl-althaus/mass_tms.out fort.666
   ln -s $Results/mimf.out fort.667
   ln -s $Results/ttg.out fort.181

   #TODO read more about it
   #BASTI files
   ln -s $BASTI/timeline.trk fort.200
   ln -s $BASTI/z0001.basti fort.201
   ln -s $BASTI/z0003.basti fort.202
   ln -s $BASTI/z0006.basti fort.203
   ln -s $BASTI/z0010.basti fort.204
   ln -s $BASTI/z0020.basti fort.205
   ln -s $BASTI/z0040.basti fort.206
   ln -s $BASTI/z0080.basti fort.207
   ln -s $BASTI/z0100.basti fort.208
   ln -s $BASTI/z0198.basti fort.209
   ln -s $BASTI/z0300.basti fort.210
   ln -s $BASTI/z0400.basti fort.211
   ln -s $BASTI/metal.trk fort.212

   #TODO read more about it
   #Variable G files
   ln -s $G_Variable/052_1d11_g1020.trk          fort.300
   ln -s $G_Variable/060_1d11_g1020.trk          fort.301
   ln -s $G_Variable/100_1d11_g1020.trk          fort.302
   ln -s $G_Variable/052_1d11_g1050.trk          fort.303
   ln -s $G_Variable/060_1d11_g1050.trk          fort.304
   ln -s $G_Variable/100_1d11_g1050.trk          fort.305
   ln -s $G_Variable/052_1d11_g1100.trk          fort.306
   ln -s $G_Variable/060_1d11_g1100.trk          fort.307
   ln -s $G_Variable/100_1d11_g1100.trk          fort.308


   #TODO What do the next 2 lines mean?
   time ./disk40.e 
   rm fort.*

   currentNumberOfExecution=`expr $currentNumberOfExecution + 1`
   currentIndexForOutput=`expr $currentIndexForOutput + 1`
done

rm ./disk40.e
