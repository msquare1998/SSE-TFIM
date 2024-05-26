# ====================
#   Params
# ====================
L=8
beta=5
J=1
h=1
nThm=10000
nStat=50000
nBins=10

# ====================
#   Compile & run
# ====================
rm -rf build && mkdir build
rm -rf data && mkdir data
cd build && cmake ..
make
./tfim $L $beta $J $h $nThm $nStat $nBins