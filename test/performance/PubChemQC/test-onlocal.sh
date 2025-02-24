#!/bin/bash
shopt -s expand_aliases

TESTSET="PubChemQC"
source $PWD/${TESTSET}/setup.sh

NPROCS=8 # depends on CPU of your local machine
echo "NPROCS: $NPROCS"

DNAME="whatever_you_want" # modify this

#XYZ2MOL=1 # uncomment this line if you want to test xyz2mol
#INDIGOX=1 # uncomment this line if you want to test IndigoX
#OBABEL=1 # uncomment this line if you want to test Openbabel

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
SNAME="${CHGSTR}_sel"

SEL="${PWD}/${TESTSET}/test-list/${SNAME}"
EXP="${PWD}/${TESTSET}/${DNAME}"

mkdir -p $EXP

STARTTIME=$(date)

parallel --jobs $NPROCS -a $SEL --timeout 60 \
  python -u ../test.py {} 0 "$modOpts" \
  1>$LOG 2>$ERR

ENDTIME=$(date)

echo "Start time: ${STARTTIME}"
echo "End time: ${ENDTIME}"
