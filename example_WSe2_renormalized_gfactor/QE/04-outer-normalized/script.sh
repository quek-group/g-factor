./writedudk.py WFN_kkxky.h5 WFN_before.h5 WFN_tempout.h5 WFN_after.h5 > CBM_write.out
rm WFN_tempout.h5
mv WFN_after.h5 WFN_after_CBM.h5
hdf2wfn.x BIN WFN_after_CBM.h5 WFN_after_CBM.complex
./writedudk_V.py WFN_kkxky.h5 WFN_before.h5 WFN_tempout.h5 WFN_after.h5 > VBM_write.out
rm WFN_tempout.h5
mv WFN_after.h5 WFN_after_VBM.h5
hdf2wfn.x BIN WFN_after_VBM.h5 WFN_after_VBM.complex
