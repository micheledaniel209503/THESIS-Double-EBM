/*
 * Continues_motor_control_v9_tensiontest_private.h
 *
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * Code generation for model "Continues_motor_control_v9_tensiontest".
 *
 * Model version              : 2.0
 * Simulink Coder version : 25.2 (R2025b) 28-Jul-2025
 * C++ source code generated on : Tue Nov 25 18:53:53 2025
 *
 * Target selection: speedgoat.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: Intel->x86-64 (Linux 64)
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#ifndef Continues_motor_control_v9_tensiontest_private_h_
#define Continues_motor_control_v9_tensiontest_private_h_
#include "rtwtypes.h"
#include "multiword_types.h"
#include "Continues_motor_control_v9_tensiontest_types.h"
#include "Continues_motor_control_v9_tensiontest.h"

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
extern void Continues_motor_control_v9_tensiontest_derivatives(void);

#endif                   /* Continues_motor_control_v9_tensiontest_private_h_ */
