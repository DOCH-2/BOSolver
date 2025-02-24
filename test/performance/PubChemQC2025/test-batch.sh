#!/bin/bash

echo "~= $HOME"

shopt -s expand_aliases
source ~/.bashrc_noninteractive
source ~/BOSolver/test/performance/${TESTSET}/setup.sh

if [ ! -z $SLURM_SUMBIT_DIR ]; then
  cd $SLURM_SUMBIT_DIR
fi

#NPROCS=$(wc -l <$PBS_NODEFILE)
NPROCS=$SLURM_CPUS_ON_NODE
echo "NPROCS: $NPROCS"

SCR=/scratch/leejinwontmp/$SLURM_JOB_ID
TDR=$SCR/tmp
LOG=$SCR/LOG.txt
ERR=$SCR/ERR.txt

DNAME="exp${DPOSTFIX}"

if [ $((XYZ2MOL)) -eq 1 ]; then
  DNAME="${DNAME}/xyz2mol"
  modOpts="--xyz2mol"
elif [ $((OBABEL)) -eq 1 ]; then
  DNAME="${DNAME}/obabel"
  modOpts="--obabel"
elif [ $((INDIGOX)) -eq 1 ]; then
  DNAME="${DNAME}/indigox"
  modOpts="--indigo"
else
  DNAME="${DNAME}/bosolver"
  modOpts="--bosolver"
fi

CHGSTR="neutral"
CHGINT=0
SNAME="${CHGSTR}_sel"

SEL="/home/leejinwon/BOSolver/test/performance/${TESTSET}/test-list/neutral_sel"
EXP="/home/leejinwon/BOSolver/test/performance/${TESTSET}/${DNAME}"

mkdir -p $EXP $TDR

STARTTIME=$(date)

parallel \
  --jobs $NPROCS --tmpdir $TDR \
  -a $SEL --timeout 60 \
  python -u ../test.py {} $CHGINT "$modOpts" \
  1>$LOG 2>$ERR

ENDTIME=$(date)

cp -f $LOG $EXP/LOG.txt
cp -f $ERR $EXP/ERR.txt
rm -rf $SCR

echo "Start time: ${STARTTIME}"
echo "End time: ${ENDTIME}"
