#!/bin/bash

# 400 ns

#export seed=`date +%s`


charmm=/home/xrliu/program/charmm/c35b2m/exec/gnu/charmm
$charmm run=NVE pdbid=aaqaa -i hyres.inp |tee aaqaa.0.out
#$charmm run=NVE pdbid=aaqaa -i hyres+nb.inp |tee aaqaa.nb.out
#$charmm run=NVE pdbid=aaqaa -i hyres_full.inp |tee aaqaa.all.out

#grep 'DYNA> *' aaqaa.0.out | awk '{print $4}' > charmm_cpu.0.dat
#grep 'DYNA> *' aaqaa.nb.out | awk '{print $4}' > charmm_cpu.nb.dat
#grep 'AVER>  *' aaqaa.all.out | awk '{print $4}' > charmm_cpu.all.dat
#grep 'DYNA EXTERN> *' aaqaa.all.out | awk '{print $5}' > charmm_cpu.hb.dat

#$charmm run=NVE pdbid=aaqaa -i hyres.inp |tee aaqaa.hyres.out
#grep 'DYNA> *' aaqaa.nvenonb.out | awk '{print $2, $4}' > charmm_cpu.hb.totE.dat
