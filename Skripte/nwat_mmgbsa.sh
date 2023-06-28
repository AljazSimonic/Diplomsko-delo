#!/bin/bash

module load NAMD/2.14-foss-2022a-mpi

#NAMD_COMMAND="/d/hpc/home/asimonic/NAMD_2.14_Linux-x86_64-multicore/namd2 +p1 +isomalloc_sync"
NAMD_COMMAND="singularity exec --bind=/d/hpc/projects/FKKT/EpCAM_MD /d/hpc/singularity/namd-2.13.sif namd2 +p1 +isomalloc_sync"
FORCEFIELD_PATH="../00_Supporting_files/toppar_c36_jul21/"
EXTRACTED_PATH=extractnumwat_4th5th
NAME=${1}

FIRST_FRAME=0
STRIDE=100
LAST_FRAME=100000

PREFIX=4th5thns/

cd ${NAME}/EpEX_${NAME}_sim${2}

mkdir 4th5thns

echo "" > nwat_mmgbsa_A.log
echo "" > nwat_mmgbsa_B.log
echo "" > nwat_mmgbsa_AB.log

FRAME=$FIRST_FRAME
mkdir temp_nwat_gbsa
export PATH=$(pwd):$PATH
while [ $FRAME -lt $LAST_FRAME ]
    do
        echo $FRAME
        #Opozorilo - uporablja polje sil CHARMM36. Za CHARMM36m je treba dodati še	parameters ${FORCEFIELD_PATH}/par_all36m_prot.prm
        CONF_FILE="structure ${EXTRACTED_PATH}/${NAME}_frame_%s${FRAME}.psf\ncoordinates ${EXTRACTED_PATH}/${NAME}_frame_%s${FRAME}.pdb\noutputName temp_nwat_gbsa/tempframe_%s${FRAME}\nfirsttimestep ${FRAME}\nparaTypeCharmm on\nparameters ${FORCEFIELD_PATH}/par_all36_prot.prm\nparameters ${FORCEFIELD_PATH}/param19.inp\nparameters ${FORCEFIELD_PATH}/par_all36_carb.prm\nparameters ${FORCEFIELD_PATH}/par_all36_lipid.prm\nparameters ${FORCEFIELD_PATH}/par_all36_na.prm\nparameters ${FORCEFIELD_PATH}/par_all36_cgenff.prm\nparameters ${FORCEFIELD_PATH}/toppar_water_ions_namd.str\ntemperature 310\nGBIS on\nsolventDielectric 78.5\nionConcentration 0.15\nalphaCutoff 15\nsasa on\nsurfaceTension 0.00542\nexclude scaled1-4\n1-4scaling 1.0\ncutoff 18\nswitching on\nswitchdist 16\npairlistdist 20\ntimestep 2\nrigidBonds all\nnonbondedFreq 2\nfullElectFrequency 4\nstepspercycle 20\nrun 0\nexit\n"

        printf "$CONF_FILE" "A" "A" "A" > ${PREFIX}temp_nwatmmgbsa.conf
        $NAMD_COMMAND ${PREFIX}temp_nwatmmgbsa.conf >> ${PREFIX}nwat_mmgbsa_A.log
        printf "$CONF_FILE" "B" "B" "B" > ${PREFIX}temp_nwatmmgbsa.conf
        $NAMD_COMMAND ${PREFIX}temp_nwatmmgbsa.conf >> ${PREFIX}nwat_mmgbsa_B.log
        printf "$CONF_FILE" "AB" "AB" "AB" > ${PREFIX}temp_nwatmmgbsa.conf
        $NAMD_COMMAND ${PREFIX}temp_nwatmmgbsa.conf >> ${PREFIX}nwat_mmgbsa_AB.log
        FRAME=`expr $FRAME + $STRIDE`
        echo "Frame $FRAME finished"
    done

cat ${PREFIX}nwat_mmgbsa_A.log | grep "ENERGY:\ \ " > ${PREFIX}nwat_mmgbsa_A.tsv
cat ${PREFIX}nwat_mmgbsa_B.log | grep "ENERGY:\ \ " > ${PREFIX}nwat_mmgbsa_B.tsv
cat ${PREFIX}nwat_mmgbsa_AB.log | grep "ENERGY:\ \ " > ${PREFIX}nwat_mmgbsa_AB.tsv
export LC_NUMERIC="C"
paste ${PREFIX}nwat_mmgbsa_AB.tsv ${PREFIX}nwat_mmgbsa_A.tsv ${PREFIX}nwat_mmgbsa_B.tsv | awk '{print $14-$30-$46}' > ${PREFIX}${NAME}_nwat_mmgbsa.tsv
rm temp_nwatmmgbsa.conf
rm -r temp_nwat_gbsa
exit
