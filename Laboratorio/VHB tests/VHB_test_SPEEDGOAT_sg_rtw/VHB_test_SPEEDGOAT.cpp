/*
 * VHB_test_SPEEDGOAT.cpp
 *
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * Code generation for model "VHB_test_SPEEDGOAT".
 *
 * Model version              : 2.2
 * Simulink Coder version : 25.2 (R2025b) 28-Jul-2025
 * C++ source code generated on : Wed Nov 26 11:30:29 2025
 *
 * Target selection: speedgoat.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: Intel->x86-64 (Linux 64)
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#include "VHB_test_SPEEDGOAT.h"
#include "VHB_test_SPEEDGOAT_cal.h"
#include "rtwtypes.h"
#include "rte_VHB_test_SPEEDGOAT_parameters.h"
#include "VHB_test_SPEEDGOAT_private.h"
#include <cstring>

extern "C"
{

#include "rt_nonfinite.h"

}

const real_T VHB_test_SPEEDGOAT_RGND = 0.0;/* real_T ground */

/* Block signals (default storage) */
B_VHB_test_SPEEDGOAT_T VHB_test_SPEEDGOAT_B;

/* Continuous states */
X_VHB_test_SPEEDGOAT_T VHB_test_SPEEDGOAT_X;

/* Disabled State Vector */
XDis_VHB_test_SPEEDGOAT_T VHB_test_SPEEDGOAT_XDis;

/* Block states (default storage) */
DW_VHB_test_SPEEDGOAT_T VHB_test_SPEEDGOAT_DW;

/* Real-time model */
RT_MODEL_VHB_test_SPEEDGOAT_T VHB_test_SPEEDGOAT_M_ =
  RT_MODEL_VHB_test_SPEEDGOAT_T();
RT_MODEL_VHB_test_SPEEDGOAT_T *const VHB_test_SPEEDGOAT_M =
  &VHB_test_SPEEDGOAT_M_;

/*
 * This function updates continuous states using the ODE3 fixed-step
 * solver algorithm
 */
static void rt_ertODEUpdateContinuousStates(RTWSolverInfo *si )
{
  /* Solver Matrices */
  static const real_T rt_ODE3_A[3] = {
    1.0/2.0, 3.0/4.0, 1.0
  };

  static const real_T rt_ODE3_B[3][3] = {
    { 1.0/2.0, 0.0, 0.0 },

    { 0.0, 3.0/4.0, 0.0 },

    { 2.0/9.0, 1.0/3.0, 4.0/9.0 }
  };

  time_T t = rtsiGetT(si);
  time_T tnew = rtsiGetSolverStopTime(si);
  time_T h = rtsiGetStepSize(si);
  real_T *x = rtsiGetContStates(si);
  ODE3_IntgData *id = static_cast<ODE3_IntgData *>(rtsiGetSolverData(si));
  real_T *y = id->y;
  real_T *f0 = id->f[0];
  real_T *f1 = id->f[1];
  real_T *f2 = id->f[2];
  real_T hB[3];
  int_T i;
  int_T nXc = 1;
  rtsiSetSimTimeStep(si,MINOR_TIME_STEP);

  /* Save the state values at time t in y, we'll use x as ynew. */
  (void) std::memcpy(y, x,
                     static_cast<uint_T>(nXc)*sizeof(real_T));

  /* Assumes that rtsiSetT and ModelOutputs are up-to-date */
  /* f0 = f(t,y) */
  rtsiSetdX(si, f0);
  VHB_test_SPEEDGOAT_derivatives();

  /* f(:,2) = feval(odefile, t + hA(1), y + f*hB(:,1), args(:)(*)); */
  hB[0] = h * rt_ODE3_B[0][0];
  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (f0[i]*hB[0]);
  }

  rtsiSetT(si, t + h*rt_ODE3_A[0]);
  rtsiSetdX(si, f1);
  VHB_test_SPEEDGOAT_step();
  VHB_test_SPEEDGOAT_derivatives();

  /* f(:,3) = feval(odefile, t + hA(2), y + f*hB(:,2), args(:)(*)); */
  for (i = 0; i <= 1; i++) {
    hB[i] = h * rt_ODE3_B[1][i];
  }

  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (f0[i]*hB[0] + f1[i]*hB[1]);
  }

  rtsiSetT(si, t + h*rt_ODE3_A[1]);
  rtsiSetdX(si, f2);
  VHB_test_SPEEDGOAT_step();
  VHB_test_SPEEDGOAT_derivatives();

  /* tnew = t + hA(3);
     ynew = y + f*hB(:,3); */
  for (i = 0; i <= 2; i++) {
    hB[i] = h * rt_ODE3_B[2][i];
  }

  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (f0[i]*hB[0] + f1[i]*hB[1] + f2[i]*hB[2]);
  }

  rtsiSetT(si, tnew);
  rtsiSetSimTimeStep(si,MAJOR_TIME_STEP);
}

/* Model step function */
void VHB_test_SPEEDGOAT_step(void)
{
  real_T currentTime;
  real_T u1;
  real_T u2;
  boolean_T tmp;
  if (rtmIsMajorTimeStep(VHB_test_SPEEDGOAT_M)) {
    /* set solver stop time */
    if (!(VHB_test_SPEEDGOAT_M->Timing.clockTick0+1)) {
      rtsiSetSolverStopTime(&VHB_test_SPEEDGOAT_M->solverInfo,
                            ((VHB_test_SPEEDGOAT_M->Timing.clockTickH0 + 1) *
        VHB_test_SPEEDGOAT_M->Timing.stepSize0 * 4294967296.0));
    } else {
      rtsiSetSolverStopTime(&VHB_test_SPEEDGOAT_M->solverInfo,
                            ((VHB_test_SPEEDGOAT_M->Timing.clockTick0 + 1) *
        VHB_test_SPEEDGOAT_M->Timing.stepSize0 +
        VHB_test_SPEEDGOAT_M->Timing.clockTickH0 *
        VHB_test_SPEEDGOAT_M->Timing.stepSize0 * 4294967296.0));
    }
  }                                    /* end MajorTimeStep */

  /* Update absolute time of base rate at minor time step */
  if (rtmIsMinorTimeStep(VHB_test_SPEEDGOAT_M)) {
    VHB_test_SPEEDGOAT_M->Timing.t[0] = rtsiGetT
      (&VHB_test_SPEEDGOAT_M->solverInfo);
  }

  /* Step: '<S4>/Step' */
  u2 = *get_T();
  tmp = rtmIsMajorTimeStep(VHB_test_SPEEDGOAT_M);
  if (tmp) {
    /* S-Function (sg_IO130_131_setup_s): '<Root>/Setup ' */

    /* Level2 S-Function Block: '<Root>/Setup ' (sg_IO130_131_setup_s) */
    {
      SimStruct *rts = VHB_test_SPEEDGOAT_M->childSfunctions[0];
      sfcnOutputs(rts,0);
    }
  }

  /* Step: '<S4>/Step' */
  currentTime = VHB_test_SPEEDGOAT_M->Timing.t[0];
  u1 = 3.0 * u2;
  if (currentTime < u1) {
    /* Step: '<S4>/Step' */
    VHB_test_SPEEDGOAT_B.Step = VHB_test_SPEEDGOAT_cal->Step_Y0;
  } else {
    /* Step: '<S4>/Step' */
    VHB_test_SPEEDGOAT_B.Step = *get_slope();
  }

  /* Clock: '<S4>/Clock' */
  VHB_test_SPEEDGOAT_B.Clock = VHB_test_SPEEDGOAT_M->Timing.t[0];

  /* Sum: '<S4>/Sum' */
  VHB_test_SPEEDGOAT_B.Sum = VHB_test_SPEEDGOAT_B.Clock - u1;

  /* Product: '<S4>/Product' */
  VHB_test_SPEEDGOAT_B.Product = VHB_test_SPEEDGOAT_B.Step *
    VHB_test_SPEEDGOAT_B.Sum;

  /* Sum: '<S4>/Output' incorporates:
   *  Constant: '<S4>/Constant1'
   */
  VHB_test_SPEEDGOAT_B.Output = VHB_test_SPEEDGOAT_B.Product +
    VHB_test_SPEEDGOAT_cal->Ramp_InitialOutput;
  if (tmp) {
    /* Step: '<Root>/Step6' */
    currentTime = VHB_test_SPEEDGOAT_M->Timing.t[1];
    u1 = *get_N_finish() * u2;
    if (currentTime < u1) {
      /* Step: '<Root>/Step6' */
      VHB_test_SPEEDGOAT_B.Step6 = VHB_test_SPEEDGOAT_cal->Step6_Y0;
    } else {
      /* Step: '<Root>/Step6' */
      VHB_test_SPEEDGOAT_B.Step6 = VHB_test_SPEEDGOAT_cal->Step6_YFinal;
    }

    /* End of Step: '<Root>/Step6' */
  }

  /* Product: '<Root>/Product' */
  VHB_test_SPEEDGOAT_B.Product_a = VHB_test_SPEEDGOAT_B.Output *
    VHB_test_SPEEDGOAT_B.Step6;

  /* Sum: '<Root>/Subtract' incorporates:
   *  Constant: '<Root>/Constant3'
   */
  VHB_test_SPEEDGOAT_B.reference = *get_bias() + VHB_test_SPEEDGOAT_B.Product_a;
  if (tmp) {
    /* S-Function (sg_IO130_131_ad_s): '<Root>/Analog input ' */

    /* Level2 S-Function Block: '<Root>/Analog input ' (sg_IO130_131_ad_s) */
    {
      SimStruct *rts = VHB_test_SPEEDGOAT_M->childSfunctions[1];
      sfcnOutputs(rts,0);
    }

    /* Gain: '<Root>/potent gain' */
    VHB_test_SPEEDGOAT_B.Pott_mm = *get_Potent_gain() *
      VHB_test_SPEEDGOAT_B.Potent;
  }

  /* Sum: '<Root>/Sum' */
  VHB_test_SPEEDGOAT_B.Sum_p = VHB_test_SPEEDGOAT_B.reference -
    VHB_test_SPEEDGOAT_B.Pott_mm;

  /* Gain: '<S3>/Gain' */
  VHB_test_SPEEDGOAT_B.Gain = *get_Kp() * VHB_test_SPEEDGOAT_B.Sum_p;

  /* Integrator: '<S3>/Integrator' */
  VHB_test_SPEEDGOAT_B.Integrator = VHB_test_SPEEDGOAT_X.Integrator_CSTATE;

  /* Gain: '<S3>/Gain1' */
  VHB_test_SPEEDGOAT_B.Gain1 = *get_Ki() * VHB_test_SPEEDGOAT_B.Integrator;

  /* Sum: '<S3>/Sum' */
  VHB_test_SPEEDGOAT_B.Sum_j = VHB_test_SPEEDGOAT_B.Gain +
    VHB_test_SPEEDGOAT_B.Gain1;

  /* Switch: '<Root>/Switch2' */
  if (VHB_test_SPEEDGOAT_B.Sum_j > VHB_test_SPEEDGOAT_cal->Switch2_Threshold) {
    /* Switch: '<Root>/Switch2' */
    VHB_test_SPEEDGOAT_B.Switch2 = VHB_test_SPEEDGOAT_B.Sum_j;
  } else {
    /* Switch: '<Root>/Switch2' incorporates:
     *  Constant: '<Root>/Constant5'
     */
    VHB_test_SPEEDGOAT_B.Switch2 = VHB_test_SPEEDGOAT_cal->Constant5_Value;
  }

  /* End of Switch: '<Root>/Switch2' */

  /* Saturate: '<Root>/Saturation' */
  currentTime = VHB_test_SPEEDGOAT_B.Switch2;
  u1 = VHB_test_SPEEDGOAT_cal->Saturation_LowerSat;
  u2 = VHB_test_SPEEDGOAT_cal->Saturation_UpperSat;
  if (currentTime > u2) {
    /* Saturate: '<Root>/Saturation' */
    VHB_test_SPEEDGOAT_B.Forward = u2;
  } else if (currentTime < u1) {
    /* Saturate: '<Root>/Saturation' */
    VHB_test_SPEEDGOAT_B.Forward = u1;
  } else {
    /* Saturate: '<Root>/Saturation' */
    VHB_test_SPEEDGOAT_B.Forward = currentTime;
  }

  /* End of Saturate: '<Root>/Saturation' */

  /* Switch: '<Root>/Switch3' */
  if (VHB_test_SPEEDGOAT_B.Sum_j >= VHB_test_SPEEDGOAT_cal->Switch3_Threshold) {
    /* Switch: '<Root>/Switch3' incorporates:
     *  Constant: '<Root>/Constant4'
     */
    VHB_test_SPEEDGOAT_B.Switch3 = VHB_test_SPEEDGOAT_cal->Constant4_Value;
  } else {
    /* Gain: '<Root>/Gain1' */
    VHB_test_SPEEDGOAT_B.Gain1_j = VHB_test_SPEEDGOAT_cal->Gain1_Gain *
      VHB_test_SPEEDGOAT_B.Sum_j;

    /* Switch: '<Root>/Switch3' */
    VHB_test_SPEEDGOAT_B.Switch3 = VHB_test_SPEEDGOAT_B.Gain1_j;
  }

  /* End of Switch: '<Root>/Switch3' */

  /* Saturate: '<Root>/Saturation3' */
  currentTime = VHB_test_SPEEDGOAT_B.Switch3;
  u1 = VHB_test_SPEEDGOAT_cal->Saturation3_LowerSat;
  u2 = VHB_test_SPEEDGOAT_cal->Saturation3_UpperSat;
  if (currentTime > u2) {
    /* Saturate: '<Root>/Saturation3' */
    VHB_test_SPEEDGOAT_B.Backward = u2;
  } else if (currentTime < u1) {
    /* Saturate: '<Root>/Saturation3' */
    VHB_test_SPEEDGOAT_B.Backward = u1;
  } else {
    /* Saturate: '<Root>/Saturation3' */
    VHB_test_SPEEDGOAT_B.Backward = currentTime;
  }

  /* End of Saturate: '<Root>/Saturation3' */
  if (tmp) {
    /* S-Function (sg_IO130_131_da_s): '<Root>/Analog output ' */

    /* Level2 S-Function Block: '<Root>/Analog output ' (sg_IO130_131_da_s) */
    {
      SimStruct *rts = VHB_test_SPEEDGOAT_M->childSfunctions[2];
      sfcnOutputs(rts,0);
    }
  }

  /* TransferFcn: '<Root>/Transfer Fcn2' */
  VHB_test_SPEEDGOAT_B.Loadcell = 0.0;
  VHB_test_SPEEDGOAT_B.Loadcell += VHB_test_SPEEDGOAT_cal->TransferFcn2_D *
    VHB_test_SPEEDGOAT_B.loadcell;
  if (tmp) {
    /* ToAsyncQueueBlock generated from: '<S1>/Transfer Fcn2' */
    slrtLogSignal
      (VHB_test_SPEEDGOAT_DW.TAQSigLogging_InsertedFor_Trans.SLRTSigHandles,
       VHB_test_SPEEDGOAT_M->Timing.t[1]);

    /* ToAsyncQueueBlock generated from: '<S2>/potent gain' */
    slrtLogSignal
      (VHB_test_SPEEDGOAT_DW.TAQSigLogging_InsertedFor_poten.SLRTSigHandles,
       VHB_test_SPEEDGOAT_M->Timing.t[1]);

    /* ToAsyncQueueBlock generated from: '<Root>/Analog input ' */
    slrtLogSignal
      (VHB_test_SPEEDGOAT_DW.TAQSigLogging_InsertedFor_Ana_p.SLRTSigHandles,
       VHB_test_SPEEDGOAT_M->Timing.t[1]);
  }

  if (rtmIsMajorTimeStep(VHB_test_SPEEDGOAT_M)) {
    rt_ertODEUpdateContinuousStates(&VHB_test_SPEEDGOAT_M->solverInfo);

    /* Update absolute time for base rate */
    /* The "clockTick0" counts the number of times the code of this task has
     * been executed. The absolute time is the multiplication of "clockTick0"
     * and "Timing.stepSize0". Size of "clockTick0" ensures timer will not
     * overflow during the application lifespan selected.
     * Timer of this task consists of two 32 bit unsigned integers.
     * The two integers represent the low bits Timing.clockTick0 and the high bits
     * Timing.clockTickH0. When the low bit overflows to 0, the high bits increment.
     */
    if (!(++VHB_test_SPEEDGOAT_M->Timing.clockTick0)) {
      ++VHB_test_SPEEDGOAT_M->Timing.clockTickH0;
    }

    VHB_test_SPEEDGOAT_M->Timing.t[0] = rtsiGetSolverStopTime
      (&VHB_test_SPEEDGOAT_M->solverInfo);

    {
      /* Update absolute timer for sample time: [0.0001s, 0.0s] */
      /* The "clockTick1" counts the number of times the code of this task has
       * been executed. The absolute time is the multiplication of "clockTick1"
       * and "Timing.stepSize1". Size of "clockTick1" ensures timer will not
       * overflow during the application lifespan selected.
       * Timer of this task consists of two 32 bit unsigned integers.
       * The two integers represent the low bits Timing.clockTick1 and the high bits
       * Timing.clockTickH1. When the low bit overflows to 0, the high bits increment.
       */
      if (!(++VHB_test_SPEEDGOAT_M->Timing.clockTick1)) {
        ++VHB_test_SPEEDGOAT_M->Timing.clockTickH1;
      }

      VHB_test_SPEEDGOAT_M->Timing.t[1] =
        VHB_test_SPEEDGOAT_M->Timing.clockTick1 *
        VHB_test_SPEEDGOAT_M->Timing.stepSize1 +
        VHB_test_SPEEDGOAT_M->Timing.clockTickH1 *
        VHB_test_SPEEDGOAT_M->Timing.stepSize1 * 4294967296.0;
    }
  }                                    /* end MajorTimeStep */
}

/* Derivatives for root system: '<Root>' */
void VHB_test_SPEEDGOAT_derivatives(void)
{
  XDot_VHB_test_SPEEDGOAT_T *_rtXdot;
  _rtXdot = ((XDot_VHB_test_SPEEDGOAT_T *) VHB_test_SPEEDGOAT_M->derivs);

  /* Derivatives for Integrator: '<S3>/Integrator' */
  _rtXdot->Integrator_CSTATE = VHB_test_SPEEDGOAT_B.Sum_p;
}

/* Model initialize function */
void VHB_test_SPEEDGOAT_initialize(void)
{
  /* Registration code */

  /* initialize non-finites */
  rt_InitInfAndNaN(sizeof(real_T));

  {
    /* Setup solver object */
    rtsiSetSimTimeStepPtr(&VHB_test_SPEEDGOAT_M->solverInfo,
                          &VHB_test_SPEEDGOAT_M->Timing.simTimeStep);
    rtsiSetTPtr(&VHB_test_SPEEDGOAT_M->solverInfo, &rtmGetTPtr
                (VHB_test_SPEEDGOAT_M));
    rtsiSetStepSizePtr(&VHB_test_SPEEDGOAT_M->solverInfo,
                       &VHB_test_SPEEDGOAT_M->Timing.stepSize0);
    rtsiSetdXPtr(&VHB_test_SPEEDGOAT_M->solverInfo,
                 &VHB_test_SPEEDGOAT_M->derivs);
    rtsiSetContStatesPtr(&VHB_test_SPEEDGOAT_M->solverInfo, (real_T **)
                         &VHB_test_SPEEDGOAT_M->contStates);
    rtsiSetNumContStatesPtr(&VHB_test_SPEEDGOAT_M->solverInfo,
      &VHB_test_SPEEDGOAT_M->Sizes.numContStates);
    rtsiSetNumPeriodicContStatesPtr(&VHB_test_SPEEDGOAT_M->solverInfo,
      &VHB_test_SPEEDGOAT_M->Sizes.numPeriodicContStates);
    rtsiSetPeriodicContStateIndicesPtr(&VHB_test_SPEEDGOAT_M->solverInfo,
      &VHB_test_SPEEDGOAT_M->periodicContStateIndices);
    rtsiSetPeriodicContStateRangesPtr(&VHB_test_SPEEDGOAT_M->solverInfo,
      &VHB_test_SPEEDGOAT_M->periodicContStateRanges);
    rtsiSetContStateDisabledPtr(&VHB_test_SPEEDGOAT_M->solverInfo, (boolean_T**)
      &VHB_test_SPEEDGOAT_M->contStateDisabled);
    rtsiSetErrorStatusPtr(&VHB_test_SPEEDGOAT_M->solverInfo, (&rtmGetErrorStatus
      (VHB_test_SPEEDGOAT_M)));
    rtsiSetRTModelPtr(&VHB_test_SPEEDGOAT_M->solverInfo, VHB_test_SPEEDGOAT_M);
  }

  rtsiSetSimTimeStep(&VHB_test_SPEEDGOAT_M->solverInfo, MAJOR_TIME_STEP);
  rtsiSetIsMinorTimeStepWithModeChange(&VHB_test_SPEEDGOAT_M->solverInfo, false);
  rtsiSetIsContModeFrozen(&VHB_test_SPEEDGOAT_M->solverInfo, false);
  VHB_test_SPEEDGOAT_M->intgData.y = VHB_test_SPEEDGOAT_M->odeY;
  VHB_test_SPEEDGOAT_M->intgData.f[0] = VHB_test_SPEEDGOAT_M->odeF[0];
  VHB_test_SPEEDGOAT_M->intgData.f[1] = VHB_test_SPEEDGOAT_M->odeF[1];
  VHB_test_SPEEDGOAT_M->intgData.f[2] = VHB_test_SPEEDGOAT_M->odeF[2];
  VHB_test_SPEEDGOAT_M->contStates = ((X_VHB_test_SPEEDGOAT_T *)
    &VHB_test_SPEEDGOAT_X);
  VHB_test_SPEEDGOAT_M->contStateDisabled = ((XDis_VHB_test_SPEEDGOAT_T *)
    &VHB_test_SPEEDGOAT_XDis);
  VHB_test_SPEEDGOAT_M->Timing.tStart = (0.0);
  rtsiSetSolverData(&VHB_test_SPEEDGOAT_M->solverInfo, static_cast<void *>
                    (&VHB_test_SPEEDGOAT_M->intgData));
  rtsiSetSolverName(&VHB_test_SPEEDGOAT_M->solverInfo,"ode3");
  VHB_test_SPEEDGOAT_M->solverInfoPtr = (&VHB_test_SPEEDGOAT_M->solverInfo);

  /* Initialize timing info */
  {
    int_T *mdlTsMap = VHB_test_SPEEDGOAT_M->Timing.sampleTimeTaskIDArray;
    mdlTsMap[0] = 0;
    mdlTsMap[1] = 1;
    VHB_test_SPEEDGOAT_M->Timing.sampleTimeTaskIDPtr = (&mdlTsMap[0]);
    VHB_test_SPEEDGOAT_M->Timing.sampleTimes =
      (&VHB_test_SPEEDGOAT_M->Timing.sampleTimesArray[0]);
    VHB_test_SPEEDGOAT_M->Timing.offsetTimes =
      (&VHB_test_SPEEDGOAT_M->Timing.offsetTimesArray[0]);

    /* task periods */
    VHB_test_SPEEDGOAT_M->Timing.sampleTimes[0] = (0.0);
    VHB_test_SPEEDGOAT_M->Timing.sampleTimes[1] = (0.0001);

    /* task offsets */
    VHB_test_SPEEDGOAT_M->Timing.offsetTimes[0] = (0.0);
    VHB_test_SPEEDGOAT_M->Timing.offsetTimes[1] = (0.0);
  }

  rtmSetTPtr(VHB_test_SPEEDGOAT_M, &VHB_test_SPEEDGOAT_M->Timing.tArray[0]);

  {
    int_T *mdlSampleHits = VHB_test_SPEEDGOAT_M->Timing.sampleHitArray;
    mdlSampleHits[0] = 1;
    mdlSampleHits[1] = 1;
    VHB_test_SPEEDGOAT_M->Timing.sampleHits = (&mdlSampleHits[0]);
  }

  rtmSetTFinal(VHB_test_SPEEDGOAT_M, -1);
  VHB_test_SPEEDGOAT_M->Timing.stepSize0 = 0.0001;
  VHB_test_SPEEDGOAT_M->Timing.stepSize1 = 0.0001;
  VHB_test_SPEEDGOAT_M->solverInfoPtr = (&VHB_test_SPEEDGOAT_M->solverInfo);
  VHB_test_SPEEDGOAT_M->Timing.stepSize = (0.0001);
  rtsiSetFixedStepSize(&VHB_test_SPEEDGOAT_M->solverInfo, 0.0001);
  rtsiSetSolverMode(&VHB_test_SPEEDGOAT_M->solverInfo, SOLVER_MODE_SINGLETASKING);

  /* block I/O */
  {
    VHB_test_SPEEDGOAT_B.Step = 0.0;
    VHB_test_SPEEDGOAT_B.Clock = 0.0;
    VHB_test_SPEEDGOAT_B.Sum = 0.0;
    VHB_test_SPEEDGOAT_B.Product = 0.0;
    VHB_test_SPEEDGOAT_B.Output = 0.0;
    VHB_test_SPEEDGOAT_B.Step6 = 0.0;
    VHB_test_SPEEDGOAT_B.Product_a = 0.0;
    VHB_test_SPEEDGOAT_B.reference = 0.0;
    VHB_test_SPEEDGOAT_B.Potent = 0.0;
    VHB_test_SPEEDGOAT_B.Laser = 0.0;
    VHB_test_SPEEDGOAT_B.Voltage = 0.0;
    VHB_test_SPEEDGOAT_B.loadcell = 0.0;
    VHB_test_SPEEDGOAT_B.Pott_mm = 0.0;
    VHB_test_SPEEDGOAT_B.Sum_p = 0.0;
    VHB_test_SPEEDGOAT_B.Gain = 0.0;
    VHB_test_SPEEDGOAT_B.Integrator = 0.0;
    VHB_test_SPEEDGOAT_B.Gain1 = 0.0;
    VHB_test_SPEEDGOAT_B.Sum_j = 0.0;
    VHB_test_SPEEDGOAT_B.Switch2 = 0.0;
    VHB_test_SPEEDGOAT_B.Forward = 0.0;
    VHB_test_SPEEDGOAT_B.Switch3 = 0.0;
    VHB_test_SPEEDGOAT_B.Backward = 0.0;
    VHB_test_SPEEDGOAT_B.Loadcell = 0.0;
    VHB_test_SPEEDGOAT_B.Gain1_j = 0.0;
  }

  /* states (continuous) */
  {
    (void) std::memset(static_cast<void *>(&VHB_test_SPEEDGOAT_X), 0,
                       sizeof(X_VHB_test_SPEEDGOAT_T));
  }

  /* disabled states */
  {
    (void) std::memset(static_cast<void *>(&VHB_test_SPEEDGOAT_XDis), 0,
                       sizeof(XDis_VHB_test_SPEEDGOAT_T));
  }

  /* states (dwork) */
  (void) std::memset(static_cast<void *>(&VHB_test_SPEEDGOAT_DW), 0,
                     sizeof(DW_VHB_test_SPEEDGOAT_T));
  VHB_test_SPEEDGOAT_DW.Analoginput_RWORK = 0.0;

  /* child S-Function registration */
  {
    RTWSfcnInfo *sfcnInfo = &VHB_test_SPEEDGOAT_M->NonInlinedSFcns.sfcnInfo;
    VHB_test_SPEEDGOAT_M->sfcnInfo = (sfcnInfo);
    rtssSetErrorStatusPtr(sfcnInfo, (&rtmGetErrorStatus(VHB_test_SPEEDGOAT_M)));
    VHB_test_SPEEDGOAT_M->Sizes.numSampTimes = (2);
    rtssSetNumRootSampTimesPtr(sfcnInfo,
      &VHB_test_SPEEDGOAT_M->Sizes.numSampTimes);
    VHB_test_SPEEDGOAT_M->NonInlinedSFcns.taskTimePtrs[0] = (&rtmGetTPtr
      (VHB_test_SPEEDGOAT_M)[0]);
    VHB_test_SPEEDGOAT_M->NonInlinedSFcns.taskTimePtrs[1] = (&rtmGetTPtr
      (VHB_test_SPEEDGOAT_M)[1]);
    rtssSetTPtrPtr(sfcnInfo,VHB_test_SPEEDGOAT_M->NonInlinedSFcns.taskTimePtrs);
    rtssSetTStartPtr(sfcnInfo, &rtmGetTStart(VHB_test_SPEEDGOAT_M));
    rtssSetTFinalPtr(sfcnInfo, &rtmGetTFinal(VHB_test_SPEEDGOAT_M));
    rtssSetTimeOfLastOutputPtr(sfcnInfo, &rtmGetTimeOfLastOutput
      (VHB_test_SPEEDGOAT_M));
    rtssSetStepSizePtr(sfcnInfo, &VHB_test_SPEEDGOAT_M->Timing.stepSize);
    rtssSetStopRequestedPtr(sfcnInfo, &rtmGetStopRequested(VHB_test_SPEEDGOAT_M));
    rtssSetDerivCacheNeedsResetPtr(sfcnInfo,
      &VHB_test_SPEEDGOAT_M->derivCacheNeedsReset);
    rtssSetZCCacheNeedsResetPtr(sfcnInfo,
      &VHB_test_SPEEDGOAT_M->zCCacheNeedsReset);
    rtssSetContTimeOutputInconsistentWithStateAtMajorStepPtr(sfcnInfo,
      &VHB_test_SPEEDGOAT_M->CTOutputIncnstWithState);
    rtssSetSampleHitsPtr(sfcnInfo, &VHB_test_SPEEDGOAT_M->Timing.sampleHits);
    rtssSetPerTaskSampleHitsPtr(sfcnInfo,
      &VHB_test_SPEEDGOAT_M->Timing.perTaskSampleHits);
    rtssSetSimModePtr(sfcnInfo, &VHB_test_SPEEDGOAT_M->simMode);
    rtssSetSolverInfoPtr(sfcnInfo, &VHB_test_SPEEDGOAT_M->solverInfoPtr);
  }

  VHB_test_SPEEDGOAT_M->Sizes.numSFcns = (3);

  /* register each child */
  {
    (void) std::memset(static_cast<void *>
                       (&VHB_test_SPEEDGOAT_M->NonInlinedSFcns.childSFunctions[0]),
                       0,
                       3*sizeof(SimStruct));
    VHB_test_SPEEDGOAT_M->childSfunctions =
      (&VHB_test_SPEEDGOAT_M->NonInlinedSFcns.childSFunctionPtrs[0]);
    VHB_test_SPEEDGOAT_M->childSfunctions[0] =
      (&VHB_test_SPEEDGOAT_M->NonInlinedSFcns.childSFunctions[0]);
    VHB_test_SPEEDGOAT_M->childSfunctions[1] =
      (&VHB_test_SPEEDGOAT_M->NonInlinedSFcns.childSFunctions[1]);
    VHB_test_SPEEDGOAT_M->childSfunctions[2] =
      (&VHB_test_SPEEDGOAT_M->NonInlinedSFcns.childSFunctions[2]);

    /* Level2 S-Function Block: VHB_test_SPEEDGOAT/<Root>/Setup  (sg_IO130_131_setup_s) */
    {
      SimStruct *rts = VHB_test_SPEEDGOAT_M->childSfunctions[0];

      /* timing info */
      time_T *sfcnPeriod =
        VHB_test_SPEEDGOAT_M->NonInlinedSFcns.Sfcn0.sfcnPeriod;
      time_T *sfcnOffset =
        VHB_test_SPEEDGOAT_M->NonInlinedSFcns.Sfcn0.sfcnOffset;
      int_T *sfcnTsMap = VHB_test_SPEEDGOAT_M->NonInlinedSFcns.Sfcn0.sfcnTsMap;
      (void) std::memset(static_cast<void*>(sfcnPeriod), 0,
                         sizeof(time_T)*1);
      (void) std::memset(static_cast<void*>(sfcnOffset), 0,
                         sizeof(time_T)*1);
      ssSetSampleTimePtr(rts, &sfcnPeriod[0]);
      ssSetOffsetTimePtr(rts, &sfcnOffset[0]);
      ssSetSampleTimeTaskIDPtr(rts, sfcnTsMap);

      {
        ssSetBlkInfo2Ptr(rts, &VHB_test_SPEEDGOAT_M->NonInlinedSFcns.blkInfo2[0]);
      }

      _ssSetBlkInfo2PortInfo2Ptr(rts,
        &VHB_test_SPEEDGOAT_M->NonInlinedSFcns.inputOutputPortInfo2[0]);

      /* Set up the mdlInfo pointer */
      ssSetRTWSfcnInfo(rts, VHB_test_SPEEDGOAT_M->sfcnInfo);

      /* Allocate memory of model methods 2 */
      {
        ssSetModelMethods2(rts, &VHB_test_SPEEDGOAT_M->NonInlinedSFcns.methods2
                           [0]);
      }

      /* Allocate memory of model methods 3 */
      {
        ssSetModelMethods3(rts, &VHB_test_SPEEDGOAT_M->NonInlinedSFcns.methods3
                           [0]);
      }

      /* Allocate memory of model methods 4 */
      {
        ssSetModelMethods4(rts, &VHB_test_SPEEDGOAT_M->NonInlinedSFcns.methods4
                           [0]);
      }

      /* Allocate memory for states auxilliary information */
      {
        ssSetStatesInfo2(rts, &VHB_test_SPEEDGOAT_M->
                         NonInlinedSFcns.statesInfo2[0]);
        ssSetPeriodicStatesInfo(rts,
          &VHB_test_SPEEDGOAT_M->NonInlinedSFcns.periodicStatesInfo[0]);
      }

      /* path info */
      ssSetModelName(rts, "Setup ");
      ssSetPath(rts, "VHB_test_SPEEDGOAT/Setup ");
      ssSetRTModel(rts,VHB_test_SPEEDGOAT_M);
      ssSetParentSS(rts, (NULL));
      ssSetRootSS(rts, rts);
      ssSetVersion(rts, SIMSTRUCT_VERSION_LEVEL2);

      /* parameters */
      {
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
      }

      /* work vectors */
      ssSetPWork(rts, (void **) &VHB_test_SPEEDGOAT_DW.Setup_PWORK);

      {
        struct _ssDWorkRecord *dWorkRecord = (struct _ssDWorkRecord *)
          &VHB_test_SPEEDGOAT_M->NonInlinedSFcns.Sfcn0.dWork;
        struct _ssDWorkAuxRecord *dWorkAuxRecord = (struct _ssDWorkAuxRecord *)
          &VHB_test_SPEEDGOAT_M->NonInlinedSFcns.Sfcn0.dWorkAux;
        ssSetSFcnDWork(rts, dWorkRecord);
        ssSetSFcnDWorkAux(rts, dWorkAuxRecord);
        ssSetNumDWorkAsInt(rts, 1);

        /* PWORK */
        ssSetDWorkWidthAsInt(rts, 0, 1);
        ssSetDWorkDataType(rts, 0,SS_POINTER);
        ssSetDWorkComplexSignal(rts, 0, 0);
        ssSetDWork(rts, 0, &VHB_test_SPEEDGOAT_DW.Setup_PWORK);
      }

      /* registration */
      sg_IO130_131_setup_s(rts);
      sfcnInitializeSizes(rts);
      sfcnInitializeSampleTimes(rts);

      /* adjust sample time */
      ssSetSampleTime(rts, 0, 0.0001);
      ssSetOffsetTime(rts, 0, 0.0);
      sfcnTsMap[0] = 1;

      /* set compiled values of dynamic vector attributes */
      ssSetNumNonsampledZCsAsInt(rts, 0);

      /* Update connectivity flags for each port */
      /* Update the BufferDstPort flags for each input port */
    }

    /* Level2 S-Function Block: VHB_test_SPEEDGOAT/<Root>/Analog input  (sg_IO130_131_ad_s) */
    {
      SimStruct *rts = VHB_test_SPEEDGOAT_M->childSfunctions[1];

      /* timing info */
      time_T *sfcnPeriod =
        VHB_test_SPEEDGOAT_M->NonInlinedSFcns.Sfcn1.sfcnPeriod;
      time_T *sfcnOffset =
        VHB_test_SPEEDGOAT_M->NonInlinedSFcns.Sfcn1.sfcnOffset;
      int_T *sfcnTsMap = VHB_test_SPEEDGOAT_M->NonInlinedSFcns.Sfcn1.sfcnTsMap;
      (void) std::memset(static_cast<void*>(sfcnPeriod), 0,
                         sizeof(time_T)*1);
      (void) std::memset(static_cast<void*>(sfcnOffset), 0,
                         sizeof(time_T)*1);
      ssSetSampleTimePtr(rts, &sfcnPeriod[0]);
      ssSetOffsetTimePtr(rts, &sfcnOffset[0]);
      ssSetSampleTimeTaskIDPtr(rts, sfcnTsMap);

      {
        ssSetBlkInfo2Ptr(rts, &VHB_test_SPEEDGOAT_M->NonInlinedSFcns.blkInfo2[1]);
      }

      _ssSetBlkInfo2PortInfo2Ptr(rts,
        &VHB_test_SPEEDGOAT_M->NonInlinedSFcns.inputOutputPortInfo2[1]);

      /* Set up the mdlInfo pointer */
      ssSetRTWSfcnInfo(rts, VHB_test_SPEEDGOAT_M->sfcnInfo);

      /* Allocate memory of model methods 2 */
      {
        ssSetModelMethods2(rts, &VHB_test_SPEEDGOAT_M->NonInlinedSFcns.methods2
                           [1]);
      }

      /* Allocate memory of model methods 3 */
      {
        ssSetModelMethods3(rts, &VHB_test_SPEEDGOAT_M->NonInlinedSFcns.methods3
                           [1]);
      }

      /* Allocate memory of model methods 4 */
      {
        ssSetModelMethods4(rts, &VHB_test_SPEEDGOAT_M->NonInlinedSFcns.methods4
                           [1]);
      }

      /* Allocate memory for states auxilliary information */
      {
        ssSetStatesInfo2(rts, &VHB_test_SPEEDGOAT_M->
                         NonInlinedSFcns.statesInfo2[1]);
        ssSetPeriodicStatesInfo(rts,
          &VHB_test_SPEEDGOAT_M->NonInlinedSFcns.periodicStatesInfo[1]);
      }

      /* outputs */
      {
        ssSetPortInfoForOutputs(rts,
          &VHB_test_SPEEDGOAT_M->NonInlinedSFcns.Sfcn1.outputPortInfo[0]);
        ssSetPortInfoForOutputs(rts,
          &VHB_test_SPEEDGOAT_M->NonInlinedSFcns.Sfcn1.outputPortInfo[0]);
        _ssSetNumOutputPorts(rts, 4);
        _ssSetPortInfo2ForOutputUnits(rts,
          &VHB_test_SPEEDGOAT_M->NonInlinedSFcns.Sfcn1.outputPortUnits[0]);
        ssSetOutputPortUnit(rts, 0, 0);
        ssSetOutputPortUnit(rts, 1, 0);
        ssSetOutputPortUnit(rts, 2, 0);
        ssSetOutputPortUnit(rts, 3, 0);
        _ssSetPortInfo2ForOutputCoSimAttribute(rts,
          &VHB_test_SPEEDGOAT_M->NonInlinedSFcns.Sfcn1.outputPortCoSimAttribute
          [0]);
        ssSetOutputPortIsContinuousQuantity(rts, 0, 0);
        ssSetOutputPortIsContinuousQuantity(rts, 1, 0);
        ssSetOutputPortIsContinuousQuantity(rts, 2, 0);
        ssSetOutputPortIsContinuousQuantity(rts, 3, 0);

        /* port 0 */
        {
          _ssSetOutputPortNumDimensions(rts, 0, 1);
          ssSetOutputPortWidthAsInt(rts, 0, 1);
          ssSetOutputPortSignal(rts, 0, ((real_T *) &VHB_test_SPEEDGOAT_B.Potent));
        }

        /* port 1 */
        {
          _ssSetOutputPortNumDimensions(rts, 1, 1);
          ssSetOutputPortWidthAsInt(rts, 1, 1);
          ssSetOutputPortSignal(rts, 1, ((real_T *) &VHB_test_SPEEDGOAT_B.Laser));
        }

        /* port 2 */
        {
          _ssSetOutputPortNumDimensions(rts, 2, 1);
          ssSetOutputPortWidthAsInt(rts, 2, 1);
          ssSetOutputPortSignal(rts, 2, ((real_T *)
            &VHB_test_SPEEDGOAT_B.Voltage));
        }

        /* port 3 */
        {
          _ssSetOutputPortNumDimensions(rts, 3, 1);
          ssSetOutputPortWidthAsInt(rts, 3, 1);
          ssSetOutputPortSignal(rts, 3, ((real_T *)
            &VHB_test_SPEEDGOAT_B.loadcell));
        }
      }

      /* path info */
      ssSetModelName(rts, "Analog input ");
      ssSetPath(rts, "VHB_test_SPEEDGOAT/Analog input ");
      ssSetRTModel(rts,VHB_test_SPEEDGOAT_M);
      ssSetParentSS(rts, (NULL));
      ssSetRootSS(rts, rts);
      ssSetVersion(rts, SIMSTRUCT_VERSION_LEVEL2);

      /* parameters */
      {
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
      }

      /* work vectors */
      ssSetRWork(rts, (real_T *) &VHB_test_SPEEDGOAT_DW.Analoginput_RWORK);
      ssSetIWork(rts, (int_T *) &VHB_test_SPEEDGOAT_DW.Analoginput_IWORK);
      ssSetPWork(rts, (void **) &VHB_test_SPEEDGOAT_DW.Analoginput_PWORK[0]);

      {
        struct _ssDWorkRecord *dWorkRecord = (struct _ssDWorkRecord *)
          &VHB_test_SPEEDGOAT_M->NonInlinedSFcns.Sfcn1.dWork;
        struct _ssDWorkAuxRecord *dWorkAuxRecord = (struct _ssDWorkAuxRecord *)
          &VHB_test_SPEEDGOAT_M->NonInlinedSFcns.Sfcn1.dWorkAux;
        ssSetSFcnDWork(rts, dWorkRecord);
        ssSetSFcnDWorkAux(rts, dWorkAuxRecord);
        ssSetNumDWorkAsInt(rts, 3);

        /* RWORK */
        ssSetDWorkWidthAsInt(rts, 0, 1);
        ssSetDWorkDataType(rts, 0,SS_DOUBLE);
        ssSetDWorkComplexSignal(rts, 0, 0);
        ssSetDWork(rts, 0, &VHB_test_SPEEDGOAT_DW.Analoginput_RWORK);

        /* IWORK */
        ssSetDWorkWidthAsInt(rts, 1, 1);
        ssSetDWorkDataType(rts, 1,SS_INTEGER);
        ssSetDWorkComplexSignal(rts, 1, 0);
        ssSetDWork(rts, 1, &VHB_test_SPEEDGOAT_DW.Analoginput_IWORK);

        /* PWORK */
        ssSetDWorkWidthAsInt(rts, 2, 3);
        ssSetDWorkDataType(rts, 2,SS_POINTER);
        ssSetDWorkComplexSignal(rts, 2, 0);
        ssSetDWork(rts, 2, &VHB_test_SPEEDGOAT_DW.Analoginput_PWORK[0]);
      }

      /* registration */
      sg_IO130_131_ad_s(rts);
      sfcnInitializeSizes(rts);
      sfcnInitializeSampleTimes(rts);

      /* adjust sample time */
      ssSetSampleTime(rts, 0, 0.0001);
      ssSetOffsetTime(rts, 0, 0.0);
      sfcnTsMap[0] = 1;

      /* set compiled values of dynamic vector attributes */
      ssSetNumNonsampledZCsAsInt(rts, 0);

      /* Update connectivity flags for each port */
      _ssSetOutputPortConnected(rts, 0, 1);
      _ssSetOutputPortConnected(rts, 1, 0);
      _ssSetOutputPortConnected(rts, 2, 0);
      _ssSetOutputPortConnected(rts, 3, 1);
      _ssSetOutputPortBeingMerged(rts, 0, 0);
      _ssSetOutputPortBeingMerged(rts, 1, 0);
      _ssSetOutputPortBeingMerged(rts, 2, 0);
      _ssSetOutputPortBeingMerged(rts, 3, 0);

      /* Update the BufferDstPort flags for each input port */
    }

    /* Level2 S-Function Block: VHB_test_SPEEDGOAT/<Root>/Analog output  (sg_IO130_131_da_s) */
    {
      SimStruct *rts = VHB_test_SPEEDGOAT_M->childSfunctions[2];

      /* timing info */
      time_T *sfcnPeriod =
        VHB_test_SPEEDGOAT_M->NonInlinedSFcns.Sfcn2.sfcnPeriod;
      time_T *sfcnOffset =
        VHB_test_SPEEDGOAT_M->NonInlinedSFcns.Sfcn2.sfcnOffset;
      int_T *sfcnTsMap = VHB_test_SPEEDGOAT_M->NonInlinedSFcns.Sfcn2.sfcnTsMap;
      (void) std::memset(static_cast<void*>(sfcnPeriod), 0,
                         sizeof(time_T)*1);
      (void) std::memset(static_cast<void*>(sfcnOffset), 0,
                         sizeof(time_T)*1);
      ssSetSampleTimePtr(rts, &sfcnPeriod[0]);
      ssSetOffsetTimePtr(rts, &sfcnOffset[0]);
      ssSetSampleTimeTaskIDPtr(rts, sfcnTsMap);

      {
        ssSetBlkInfo2Ptr(rts, &VHB_test_SPEEDGOAT_M->NonInlinedSFcns.blkInfo2[2]);
      }

      _ssSetBlkInfo2PortInfo2Ptr(rts,
        &VHB_test_SPEEDGOAT_M->NonInlinedSFcns.inputOutputPortInfo2[2]);

      /* Set up the mdlInfo pointer */
      ssSetRTWSfcnInfo(rts, VHB_test_SPEEDGOAT_M->sfcnInfo);

      /* Allocate memory of model methods 2 */
      {
        ssSetModelMethods2(rts, &VHB_test_SPEEDGOAT_M->NonInlinedSFcns.methods2
                           [2]);
      }

      /* Allocate memory of model methods 3 */
      {
        ssSetModelMethods3(rts, &VHB_test_SPEEDGOAT_M->NonInlinedSFcns.methods3
                           [2]);
      }

      /* Allocate memory of model methods 4 */
      {
        ssSetModelMethods4(rts, &VHB_test_SPEEDGOAT_M->NonInlinedSFcns.methods4
                           [2]);
      }

      /* Allocate memory for states auxilliary information */
      {
        ssSetStatesInfo2(rts, &VHB_test_SPEEDGOAT_M->
                         NonInlinedSFcns.statesInfo2[2]);
        ssSetPeriodicStatesInfo(rts,
          &VHB_test_SPEEDGOAT_M->NonInlinedSFcns.periodicStatesInfo[2]);
      }

      /* inputs */
      {
        _ssSetNumInputPorts(rts, 3);
        ssSetPortInfoForInputs(rts,
          &VHB_test_SPEEDGOAT_M->NonInlinedSFcns.Sfcn2.inputPortInfo[0]);
        ssSetPortInfoForInputs(rts,
          &VHB_test_SPEEDGOAT_M->NonInlinedSFcns.Sfcn2.inputPortInfo[0]);
        _ssSetPortInfo2ForInputUnits(rts,
          &VHB_test_SPEEDGOAT_M->NonInlinedSFcns.Sfcn2.inputPortUnits[0]);
        ssSetInputPortUnit(rts, 0, 0);
        ssSetInputPortUnit(rts, 1, 0);
        ssSetInputPortUnit(rts, 2, 0);
        _ssSetPortInfo2ForInputCoSimAttribute(rts,
          &VHB_test_SPEEDGOAT_M->NonInlinedSFcns.Sfcn2.inputPortCoSimAttribute[0]);
        ssSetInputPortIsContinuousQuantity(rts, 0, 0);
        ssSetInputPortIsContinuousQuantity(rts, 1, 0);
        ssSetInputPortIsContinuousQuantity(rts, 2, 0);

        /* port 0 */
        {
          ssSetInputPortRequiredContiguous(rts, 0, 1);
          ssSetInputPortSignal(rts, 0, &VHB_test_SPEEDGOAT_B.Forward);
          _ssSetInputPortNumDimensions(rts, 0, 1);
          ssSetInputPortWidthAsInt(rts, 0, 1);
        }

        /* port 1 */
        {
          ssSetInputPortRequiredContiguous(rts, 1, 1);
          ssSetInputPortSignal(rts, 1, &VHB_test_SPEEDGOAT_B.Backward);
          _ssSetInputPortNumDimensions(rts, 1, 1);
          ssSetInputPortWidthAsInt(rts, 1, 1);
        }

        /* port 2 */
        {
          ssSetInputPortRequiredContiguous(rts, 2, 1);
          ssSetInputPortSignal(rts, 2, (const_cast<real_T*>
            (&VHB_test_SPEEDGOAT_RGND)));
          _ssSetInputPortNumDimensions(rts, 2, 1);
          ssSetInputPortWidthAsInt(rts, 2, 1);
        }
      }

      /* path info */
      ssSetModelName(rts, "Analog output ");
      ssSetPath(rts, "VHB_test_SPEEDGOAT/Analog output ");
      ssSetRTModel(rts,VHB_test_SPEEDGOAT_M);
      ssSetParentSS(rts, (NULL));
      ssSetRootSS(rts, rts);
      ssSetVersion(rts, SIMSTRUCT_VERSION_LEVEL2);

      /* parameters */
      {
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
      }

      /* work vectors */
      ssSetIWork(rts, (int_T *) &VHB_test_SPEEDGOAT_DW.Analogoutput_IWORK[0]);
      ssSetPWork(rts, (void **) &VHB_test_SPEEDGOAT_DW.Analogoutput_PWORK[0]);

      {
        struct _ssDWorkRecord *dWorkRecord = (struct _ssDWorkRecord *)
          &VHB_test_SPEEDGOAT_M->NonInlinedSFcns.Sfcn2.dWork;
        struct _ssDWorkAuxRecord *dWorkAuxRecord = (struct _ssDWorkAuxRecord *)
          &VHB_test_SPEEDGOAT_M->NonInlinedSFcns.Sfcn2.dWorkAux;
        ssSetSFcnDWork(rts, dWorkRecord);
        ssSetSFcnDWorkAux(rts, dWorkAuxRecord);
        ssSetNumDWorkAsInt(rts, 2);

        /* IWORK */
        ssSetDWorkWidthAsInt(rts, 0, 2);
        ssSetDWorkDataType(rts, 0,SS_INTEGER);
        ssSetDWorkComplexSignal(rts, 0, 0);
        ssSetDWork(rts, 0, &VHB_test_SPEEDGOAT_DW.Analogoutput_IWORK[0]);

        /* PWORK */
        ssSetDWorkWidthAsInt(rts, 1, 3);
        ssSetDWorkDataType(rts, 1,SS_POINTER);
        ssSetDWorkComplexSignal(rts, 1, 0);
        ssSetDWork(rts, 1, &VHB_test_SPEEDGOAT_DW.Analogoutput_PWORK[0]);
      }

      /* registration */
      sg_IO130_131_da_s(rts);
      sfcnInitializeSizes(rts);
      sfcnInitializeSampleTimes(rts);

      /* adjust sample time */
      ssSetSampleTime(rts, 0, 0.0001);
      ssSetOffsetTime(rts, 0, 0.0);
      sfcnTsMap[0] = 1;

      /* set compiled values of dynamic vector attributes */
      ssSetNumNonsampledZCsAsInt(rts, 0);

      /* Update connectivity flags for each port */
      _ssSetInputPortConnected(rts, 0, 1);
      _ssSetInputPortConnected(rts, 1, 1);
      _ssSetInputPortConnected(rts, 2, 0);

      /* Update the BufferDstPort flags for each input port */
      ssSetInputPortBufferDstPort(rts, 0, -1);
      ssSetInputPortBufferDstPort(rts, 1, -1);
      ssSetInputPortBufferDstPort(rts, 2, -1);
    }
  }

  /* Start for S-Function (sg_IO130_131_setup_s): '<Root>/Setup ' */
  /* Level2 S-Function Block: '<Root>/Setup ' (sg_IO130_131_setup_s) */
  {
    SimStruct *rts = VHB_test_SPEEDGOAT_M->childSfunctions[0];
    sfcnStart(rts);
    if (ssGetErrorStatus(rts) != (NULL))
      return;
  }

  /* Start for S-Function (sg_IO130_131_ad_s): '<Root>/Analog input ' */
  /* Level2 S-Function Block: '<Root>/Analog input ' (sg_IO130_131_ad_s) */
  {
    SimStruct *rts = VHB_test_SPEEDGOAT_M->childSfunctions[1];
    sfcnStart(rts);
    if (ssGetErrorStatus(rts) != (NULL))
      return;
  }

  /* Start for S-Function (sg_IO130_131_da_s): '<Root>/Analog output ' */
  /* Level2 S-Function Block: '<Root>/Analog output ' (sg_IO130_131_da_s) */
  {
    SimStruct *rts = VHB_test_SPEEDGOAT_M->childSfunctions[2];
    sfcnStart(rts);
    if (ssGetErrorStatus(rts) != (NULL))
      return;
  }

  /* Start for ToAsyncQueueBlock generated from: '<S1>/Transfer Fcn2' */
  VHB_test_SPEEDGOAT_DW.TAQSigLogging_InsertedFor_Trans.SLRTSigHandles =
    slrtRegisterSignalToLoggingService(reinterpret_cast<uintptr_t>
    (&VHB_test_SPEEDGOAT_B.Loadcell));

  /* Start for ToAsyncQueueBlock generated from: '<S2>/potent gain' */
  VHB_test_SPEEDGOAT_DW.TAQSigLogging_InsertedFor_poten.SLRTSigHandles =
    slrtRegisterSignalToLoggingService(reinterpret_cast<uintptr_t>
    (&VHB_test_SPEEDGOAT_B.Pott_mm));

  /* Start for ToAsyncQueueBlock generated from: '<Root>/Analog input ' */
  VHB_test_SPEEDGOAT_DW.TAQSigLogging_InsertedFor_Ana_p.SLRTSigHandles =
    slrtRegisterSignalToLoggingService(reinterpret_cast<uintptr_t>
    (&VHB_test_SPEEDGOAT_B.Potent));

  /* InitializeConditions for Integrator: '<S3>/Integrator' */
  VHB_test_SPEEDGOAT_X.Integrator_CSTATE = VHB_test_SPEEDGOAT_cal->Integrator_IC;
}

/* Model terminate function */
void VHB_test_SPEEDGOAT_terminate(void)
{
  /* Terminate for S-Function (sg_IO130_131_setup_s): '<Root>/Setup ' */
  /* Level2 S-Function Block: '<Root>/Setup ' (sg_IO130_131_setup_s) */
  {
    SimStruct *rts = VHB_test_SPEEDGOAT_M->childSfunctions[0];
    sfcnTerminate(rts);
  }

  /* Terminate for S-Function (sg_IO130_131_ad_s): '<Root>/Analog input ' */
  /* Level2 S-Function Block: '<Root>/Analog input ' (sg_IO130_131_ad_s) */
  {
    SimStruct *rts = VHB_test_SPEEDGOAT_M->childSfunctions[1];
    sfcnTerminate(rts);
  }

  /* Terminate for S-Function (sg_IO130_131_da_s): '<Root>/Analog output ' */
  /* Level2 S-Function Block: '<Root>/Analog output ' (sg_IO130_131_da_s) */
  {
    SimStruct *rts = VHB_test_SPEEDGOAT_M->childSfunctions[2];
    sfcnTerminate(rts);
  }
}
