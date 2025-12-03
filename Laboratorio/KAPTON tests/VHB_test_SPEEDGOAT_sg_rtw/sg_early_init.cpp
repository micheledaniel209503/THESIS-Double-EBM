#include "VHB_test_SPEEDGOAT.h"
#include "VHB_test_SPEEDGOAT_private.h"
#include "sg_printf.h"
#include "sg_early_init.h"
#include "simstruc.h" // This is required if there are no S-function blocks in the model

static RTWSfcnInfo sgEarlyInitSfcnInfo;
struct _ssBlkInfo2 sgEarlyInitBlkInfo2 = {.rtwSfcnInfo = &sgEarlyInitSfcnInfo};

void sg_init_sfcns(void)
{
    sg_printf(sg_debug, "Initializing [EARLY INIT] S-functions\n");
    
    rtssSetErrorStatusPtr(&sgEarlyInitSfcnInfo, (&rtmGetErrorStatus(VHB_test_SPEEDGOAT_M)));
    sg_early_init_set_blkInfo2((void*)&sgEarlyInitBlkInfo2);
    
    // Level2 S-Function Block: VHB_test_SPEEDGOAT/<Root>/Setup  (sg_IO130_131_setup_s)
    {
        SimStruct *rts = &VHB_test_SPEEDGOAT_M->NonInlinedSFcns.childSFunctions[0];
        rts->blkInfo.blkInfo2 = (struct _ssBlkInfo2*)sg_early_init_get_blkInfo2();
        mxArray **sfcnParams = (mxArray **)
        &VHB_test_SPEEDGOAT_M->NonInlinedSFcns.Sfcn0.params;
        ssSetSFcnParamsCount(rts, 6);
        ssSetSFcnParamsPtr(rts, &sfcnParams[0]);
        ssSetSFcnParam(rts, 0, (mxArray*)VHB_test_SPEEDGOAT_cal->Setup_P1_Size);
        ssSetSFcnParam(rts, 1, (mxArray*)VHB_test_SPEEDGOAT_cal->Setup_P2_Size);
        ssSetSFcnParam(rts, 2, (mxArray*)VHB_test_SPEEDGOAT_cal->Setup_P3_Size);
        ssSetSFcnParam(rts, 3, (mxArray*)VHB_test_SPEEDGOAT_cal->Setup_P4_Size);
        ssSetSFcnParam(rts, 4, (mxArray*)VHB_test_SPEEDGOAT_cal->Setup_P5_Size);
        ssSetSFcnParam(rts, 5, (mxArray*)VHB_test_SPEEDGOAT_cal->Setup_P6_Size);
        sg_IO130_131_setup_s(rts);
        sfcnStart(rts);
        if (ssGetErrorStatus(rts) != (NULL))
            return;
    }
    
    // Level2 S-Function Block: VHB_test_SPEEDGOAT/<Root>/Analog input  (sg_IO130_131_ad_s)
    {
        SimStruct *rts = &VHB_test_SPEEDGOAT_M->NonInlinedSFcns.childSFunctions[1];
        rts->blkInfo.blkInfo2 = (struct _ssBlkInfo2*)sg_early_init_get_blkInfo2();
        mxArray **sfcnParams = (mxArray **)
        &VHB_test_SPEEDGOAT_M->NonInlinedSFcns.Sfcn1.params;
        ssSetSFcnParamsCount(rts, 12);
        ssSetSFcnParamsPtr(rts, &sfcnParams[0]);
        ssSetSFcnParam(rts, 0, (mxArray*)
        VHB_test_SPEEDGOAT_cal->Analoginput_P1_Size);
        ssSetSFcnParam(rts, 1, (mxArray*)
        VHB_test_SPEEDGOAT_cal->Analoginput_P2_Size);
        ssSetSFcnParam(rts, 2, (mxArray*)
        VHB_test_SPEEDGOAT_cal->Analoginput_P3_Size);
        ssSetSFcnParam(rts, 3, (mxArray*)
        VHB_test_SPEEDGOAT_cal->Analoginput_P4_Size);
        ssSetSFcnParam(rts, 4, (mxArray*)
        VHB_test_SPEEDGOAT_cal->Analoginput_P5_Size);
        ssSetSFcnParam(rts, 5, (mxArray*)
        VHB_test_SPEEDGOAT_cal->Analoginput_P6_Size);
        ssSetSFcnParam(rts, 6, (mxArray*)
        VHB_test_SPEEDGOAT_cal->Analoginput_P7_Size);
        ssSetSFcnParam(rts, 7, (mxArray*)
        VHB_test_SPEEDGOAT_cal->Analoginput_P8_Size);
        ssSetSFcnParam(rts, 8, (mxArray*)
        VHB_test_SPEEDGOAT_cal->Analoginput_P9_Size);
        ssSetSFcnParam(rts, 9, (mxArray*)
        VHB_test_SPEEDGOAT_cal->Analoginput_P10_Size);
        ssSetSFcnParam(rts, 10, (mxArray*)
        VHB_test_SPEEDGOAT_cal->Analoginput_P11_Size);
        ssSetSFcnParam(rts, 11, (mxArray*)
        VHB_test_SPEEDGOAT_cal->Analoginput_P12_Size);
        sg_IO130_131_ad_s(rts);
        sfcnStart(rts);
        if (ssGetErrorStatus(rts) != (NULL))
            return;
    }
    
    // Level2 S-Function Block: VHB_test_SPEEDGOAT/<Root>/Analog output  (sg_IO130_131_da_s)
    {
        SimStruct *rts = &VHB_test_SPEEDGOAT_M->NonInlinedSFcns.childSFunctions[2];
        rts->blkInfo.blkInfo2 = (struct _ssBlkInfo2*)sg_early_init_get_blkInfo2();
        mxArray **sfcnParams = (mxArray **)
        &VHB_test_SPEEDGOAT_M->NonInlinedSFcns.Sfcn2.params;
        ssSetSFcnParamsCount(rts, 15);
        ssSetSFcnParamsPtr(rts, &sfcnParams[0]);
        ssSetSFcnParam(rts, 0, (mxArray*)
        VHB_test_SPEEDGOAT_cal->Analogoutput_P1_Size);
        ssSetSFcnParam(rts, 1, (mxArray*)
        VHB_test_SPEEDGOAT_cal->Analogoutput_P2_Size);
        ssSetSFcnParam(rts, 2, (mxArray*)
        VHB_test_SPEEDGOAT_cal->Analogoutput_P3_Size);
        ssSetSFcnParam(rts, 3, (mxArray*)
        VHB_test_SPEEDGOAT_cal->Analogoutput_P4_Size);
        ssSetSFcnParam(rts, 4, (mxArray*)
        VHB_test_SPEEDGOAT_cal->Analogoutput_P5_Size);
        ssSetSFcnParam(rts, 5, (mxArray*)
        VHB_test_SPEEDGOAT_cal->Analogoutput_P6_Size);
        ssSetSFcnParam(rts, 6, (mxArray*)
        VHB_test_SPEEDGOAT_cal->Analogoutput_P7_Size);
        ssSetSFcnParam(rts, 7, (mxArray*)
        VHB_test_SPEEDGOAT_cal->Analogoutput_P8_Size);
        ssSetSFcnParam(rts, 8, (mxArray*)
        VHB_test_SPEEDGOAT_cal->Analogoutput_P9_Size);
        ssSetSFcnParam(rts, 9, (mxArray*)
        VHB_test_SPEEDGOAT_cal->Analogoutput_P10_Size);
        ssSetSFcnParam(rts, 10, (mxArray*)
        VHB_test_SPEEDGOAT_cal->Analogoutput_P11_Size);
        ssSetSFcnParam(rts, 11, (mxArray*)
        VHB_test_SPEEDGOAT_cal->Analogoutput_P12_Size);
        ssSetSFcnParam(rts, 12, (mxArray*)
        VHB_test_SPEEDGOAT_cal->Analogoutput_P13_Size);
        ssSetSFcnParam(rts, 13, (mxArray*)
        VHB_test_SPEEDGOAT_cal->Analogoutput_P14_Size);
        ssSetSFcnParam(rts, 14, (mxArray*)
        VHB_test_SPEEDGOAT_cal->Analogoutput_P15_Size);
        sg_IO130_131_da_s(rts);
        sfcnStart(rts);
        if (ssGetErrorStatus(rts) != (NULL))
            return;
    }
}

__attribute__((constructor)) void early_init_setup(void)
{
    sg_register_early_init_function(sg_init_sfcns);
}
