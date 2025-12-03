/*
 * VHB_test_SPEEDGOAT_private.h
 *
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * Code generation for model "VHB_test_SPEEDGOAT".
 *
 * Model version              : 2.3
 * Simulink Coder version : 25.2 (R2025b) 28-Jul-2025
 * C++ source code generated on : Thu Nov 27 15:57:42 2025
 *
 * Target selection: speedgoat.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: Intel->x86-64 (Linux 64)
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#ifndef VHB_test_SPEEDGOAT_private_h_
#define VHB_test_SPEEDGOAT_private_h_
#include "rtwtypes.h"
#include "multiword_types.h"
#include "VHB_test_SPEEDGOAT_types.h"
#include "VHB_test_SPEEDGOAT.h"

/* Private macros used by the generated code to access rtModel */
#ifndef rtmIsMajorTimeStep
#define rtmIsMajorTimeStep(rtm)        (((rtm)->Timing.simTimeStep) == MAJOR_TIME_STEP)
#endif

#ifndef rtmIsMinorTimeStep
#define rtmIsMinorTimeStep(rtm)        (((rtm)->Timing.simTimeStep) == MINOR_TIME_STEP)
#endif

#ifndef rtmSetTFinal
#define rtmSetTFinal(rtm, val)         ((rtm)->Timing.tFinal = (val))
#endif

#ifndef rtmSetTPtr
#define rtmSetTPtr(rtm, val)           ((rtm)->Timing.t = (val))
#endif

extern void* slrtRegisterSignalToLoggingService(uintptr_t sigAddr);
extern "C" void sg_IO130_131_setup_s(SimStruct *rts);
extern "C" void sg_IO130_131_ad_s(SimStruct *rts);
extern "C" void sg_IO130_131_da_s(SimStruct *rts);

/* private model entry point functions */
extern void VHB_test_SPEEDGOAT_derivatives(void);

#endif                                 /* VHB_test_SPEEDGOAT_private_h_ */
