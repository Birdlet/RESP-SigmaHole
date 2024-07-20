#!/bin/bash

base=`pwd`
_pwd=`pwd`/example


for TEST in BBr BCl BI C4 CF2 CF3 CH4 F3Cl SC4 SCCO ThiaS ThioS
do
    
    NC=0
    
    echo "Run PSI4 for ESP Calculation"
    cd $_pwd
    mkdir ${TEST}
    cd $_pwd/${TEST}
    mkdir PSI4; cd PSI4
    # conda activate psi4
    obabel -imol2 ../../${TEST}.mol2 -oxyz -O ${TEST}.xyz
    python $base/psi4_esp.py ${TEST}.xyz $NC
    
    echo "Generate Original ITP"
    cd $_pwd
    cd $_pwd/${TEST}
    # conda activate psi4
    acpype -i PSI4/opt.mol2 -n $NC -a gaff2 -k "maxcyc=0"
    
    
    echo "Prepare for RESP"
    cd $_pwd/${TEST}
    # conda activate psi4
    mkdir ATC; cd ATC
    antechamber -i ../PSI4/gas.gesp -fi gesp -o gas.mol2 -fo mol2 -at sybyl -c resp -nc $NC
    ## antechamber -i ../PSI4/sol.gesp -fi gesp -o sol.mol2 -fo mol2 -at sybyl -c resp -nc $NC
    
    
    echo "Add Dummy ATOMS Coordinates"
    cd $_pwd/${TEST}
    # conda activate rdkit
    echo "  Load opt.mol2"
    python $base/rdkit2gmx.py -l PSI4/opt.mol2 -t opt.acpype/opt_GMX.itp -o ${TEST}
    echo "  Load ${TEST}.mol2"
    python $base/rdkit2gmx.py -l ../${TEST}.mol2 -t opt.acpype/opt_GMX.itp -o ${TEST} --mol2_only
    python $base/rdkit2gmx.py -l ATC/gas.mol2 -t opt.acpype/opt_GMX.itp -o gas --mol2_only
    
    
    echo "Run 2 Steps RESP"
    cd $_pwd/${TEST}
    # conda activate psi4
    mkdir RESP; cd RESP
    rm -rf MOD1 MOD2
    mkdir MOD1; mkdir MOD2
    cp ../ATC/ANTECHAMBER_RESP1.IN .
    cp ../ATC/ANTECHAMBER_RESP2.IN .
    cp ../ATC/ANTECHAMBER.ESP .
    python $base/run_resp.py -l ../${TEST}_LP.gro -t ../${TEST}_LP.itp -m ../${TEST}_LP.mol2 -dp ../${TEST}_LP.pairs -n $NC -o ../${TEST}
    # #  cd MOD1; resp -i ANTECHAMBER_RESP1_MOD.IN -o output -p punch -q qin -t qout -e ANTECHAMBER_MOD.ESP -s esout
    obabel -imol2 ../${TEST}_LP_RESP.mol2 -ogro -O ../${TEST}_LP.gro
    
    echo "Generate ESP test files"
    cd $_pwd/${TEST}
    echo "  QM ESP"
    python $base/gesp2mol2.py RESP/ANTECHAMBER.ESP  ${TEST}_esp.mol2
    # conda activate rdkit
    echo "  GAFF ESP"
    python $base/rebuild_esp.py -l ATC/gas.mol2 -g PSI4/grid.dat -o ffESP.mol2
    echo "  GAFF-SH ESP"
    python $base/rebuild_esp.py -l gas_LP.mol2 -g PSI4/grid.dat -o ffESPLP.mol2
    
    cd $base
    
done
