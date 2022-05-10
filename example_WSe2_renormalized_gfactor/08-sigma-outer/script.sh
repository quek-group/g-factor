ln -snf WFN_after_CBM.complex WFN_outer
$MPIRUN $SIGMA  &> ./sigma.out
cp sigma_hp.log CBM_sigma_hp.log
ln -snf WFN_after_VBM.complex WFN_outer
rm vxc.dat x.dat
$MPIRUN $SIGMA  &> ./sigma.out
cp sigma_hp.log VBM_sigma_hp.log
