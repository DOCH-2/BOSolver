JNAMEBASE="PubChemQC"
DPOSTFIX="250224"
TESTSET="PubChemQC2025"
#CHG="p"
#CHGODR="1"
#XYZ2MOL=1 # uncomment this line if you want to test xyz2mol
#INDIGOX=1 # uncomment this line if you want to test indigoX
#OBABEL=1 # uncomment this line if you want to test OBABEL
PARTITION="20core"
#PARTITION="16core"
#NODE="cpu1"

if [ $PARTITION == "20core" ]; then
  NCPUS=20
elif [ $PARTITION == "16core" ]; then
  NCPUS=16
fi

if [ ! -v DPOSTFIX ]; then
  JNAMEPOSTFIX=""
  echo "no postfix"
else
  JNAMEPOSTFIX="${DPOSTFIX}"
fi

if [ -v XYZ2MOL ]; then
  JNAMEPOSTFIX="${JNAMEPOSTFIX}_xm"
  echo "xyz2mol"
fi

if [ -v INDIGOX ]; then
  JNAMEPOSTFIX="${JNAMEPOSTFIX}_ix"
  echo "indigox"
fi

if [ -v OBABEL ]; then
  JNAMEPOSTFIX="${JNAMEPOSTFIX}_ob"
  echo "obabel"
fi

echo "${JNAMEBASE}${JNAMEPOSTFIX}"

sbatch --job-name="${JNAMEBASE}${JNAMEPOSTFIX}" \
  -p ${PARTITION} \
  -n 1 -c ${NCPUS} \
  -o "slurm/exp${JNAMEPOSTFIX}.test.out" \
  -e "slurm/exp${JNAMEPOSTFIX}.test.err" \
  --export TESTSET="${TESTSET}",DPOSTFIX="${DPOSTFIX}",XYZ2MOL="${XYZ2MOL}",INDIGOX="${INDIGOX}",OBABEL="${OBABEL}" \
  test-batch.sh #-d afterany:3030 \

unset XYZ2MOL
unset INDIGOX
unset OBABEL
unset JNAMEPOSTFIX
unset CHG
unset CHGODR
unset PARTITION
unset NCPUS
