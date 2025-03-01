echo "BOS"
bash test-onlocal.sh
echo ""

echo "X2M"
export XYZ2MOL=1
bash test-onlocal.sh
unset XYZ2MOL
echo ""

echo "IGX"
export INDIGOX=1
bash test-onlocal.sh
unset INDIGOX
echo ""

echo "OBB"
export OBABEL=1
bash test-onlocal.sh
unset OBABEL
echo ""
