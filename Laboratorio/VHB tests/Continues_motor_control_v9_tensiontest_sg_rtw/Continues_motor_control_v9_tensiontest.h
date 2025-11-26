/*
 * Continues_motor_control_v9_tensiontest.h
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

#ifndef Continues_motor_control_v9_tensiontest_h_
#define Continues_motor_control_v9_tensiontest_h_
#include <logsrv.h>
#include "rtwtypes.h"
#include "simstruc.h"
#include "fixedpoint.h"
#include "rtw_extmode.h"
#include "sysran_types.h"
#include "rtw_continuous.h"
#include "rtw_solver.h"
#include "Continues_motor_control_v9_tensiontest_types.h"
#include <stddef.h>
#include <cstring>
#include "Continues_motor_control_v9_tensiontest_cal.h"

extern "C"
{

#include "rt_nonfinite.h"

}

/* Macros for accessing real-time model data structure */
#ifndef rtmGetContStateDisabled
#define rtmGetContStateDisabled(rtm)   ((rtm)->contStateDisabled)
#endif

#ifndef rtmSetContStateDisabled
#define rtmSetContStateDisabled(rtm, val) ((rtm)->contStateDisabled = (val))
#endif

#ifndef rtmGetContStates
#define rtmGetContStates(rtm)          ((rtm)->contStates)
#endif

#ifndef rtmSetContStates
#define rtmSetContStates(rtm, val)     ((rtm)->contStates = (val))
#endif

#ifndef rtmGetContTimeOutputInconsistentWithStateAtMajorStepFlag
#define rtmGetContTimeOutputInconsistentWithStateAtMajorStepFlag(rtm) ((rtm)->CTOutputIncnstWithState)
#endif

#ifndef rtmSetContTimeOutputInconsistentWithStateAtMajorStepFlag
#define rtmSetContTimeOutputInconsistentWithStateAtMajorStepFlag(rtm, val) ((rtm)->CTOutputIncnstWithState = (val))
#endif

#ifndef rtmGetDerivCacheNeedsReset
#define rtmGetDerivCacheNeedsReset(rtm) ((rtm)->derivCacheNeedsReset)
#endif

#ifndef rtmSetDerivCacheNeedsReset
#define rtmSetDerivCacheNeedsReset(rtm, val) ((rtm)->derivCacheNeedsReset = (val))
#endif

#ifndef rtmGetFinalTime
#define rtmGetFinalTime(rtm)           ((rtm)->Timing.tFinal)
#endif

#ifndef rtmGetIntgData
#define rtmGetIntgData(rtm)            ((rtm)->intgData)
#endif

#ifndef rtmSetIntgData
#define rtmSetIntgData(rtm, val)       ((rtm)->intgData = (val))
#endif

#ifndef rtmGetOdeF
#define rtmGetOdeF(rtm)                ((rtm)->odeF)
#endif

#ifndef rtmSetOdeF
#define rtmSetOdeF(rtm, val)           ((rtm)->odeF = (val))
#endif

#ifndef rtmGetOdeY
#define rtmGetOdeY(rtm)                ((rtm)->odeY)
#endif

#ifndef rtmSetOdeY
#define rtmSetOdeY(rtm, val)           ((rtm)->odeY = (val))
#endif

#ifndef rtmGetPeriodicContStateIndices
#define rtmGetPeriodicContStateIndices(rtm) ((rtm)->periodicContStateIndices)
#endif

#ifndef rtmSetPeriodicContStateIndices
#define rtmSetPeriodicContStateIndices(rtm, val) ((rtm)->periodicContStateIndices = (val))
#endif

#ifndef rtmGetPeriodicContStateRanges
#define rtmGetPeriodicContStateRanges(rtm) ((rtm)->periodicContStateRanges)
#endif

#ifndef rtmSetPeriodicContStateRanges
#define rtmSetPeriodicContStateRanges(rtm, val) ((rtm)->periodicContStateRanges = (val))
#endif

#ifndef rtmGetSampleHitArray
#define rtmGetSampleHitArray(rtm)      ((rtm)->Timing.sampleHitArray)
#endif

#ifndef rtmGetStepSize
#define rtmGetStepSize(rtm)            ((rtm)->Timing.stepSize)
#endif

#ifndef rtmGetZCCacheNeedsReset
#define rtmGetZCCacheNeedsReset(rtm)   ((rtm)->zCCacheNeedsReset)
#endif

#ifndef rtmSetZCCacheNeedsReset
#define rtmSetZCCacheNeedsReset(rtm, val) ((rtm)->zCCacheNeedsReset = (val))
#endif

#ifndef rtmGet_TimeOfLastOutput
#define rtmGet_TimeOfLastOutput(rtm)   ((rtm)->Timing.timeOfLastOutput)
#endif

#ifndef rtmGetdX
#define rtmGetdX(rtm)                  ((rtm)->derivs)
#endif

#ifndef rtmSetdX
#define rtmSetdX(rtm, val)             ((rtm)->derivs = (val))
#endif

#ifndef rtmGetErrorStatus
#define rtmGetErrorStatus(rtm)         ((rtm)->errorStatus)
#endif

#ifndef rtmSetErrorStatus
#define rtmSetErrorStatus(rtm, val)    ((rtm)->errorStatus = (val))
#endif

#ifndef rtmGetStopRequested
#define rtmGetStopRequested(rtm)       ((rtm)->Timing.stopRequestedFlag)
#endif

#ifndef rtmSetStopRequested
#define rtmSetStopRequested(rtm, val)  ((rtm)->Timing.stopRequestedFlag = (val))
#endif

#ifndef rtmGetStopRequestedPtr
#define rtmGetStopRequestedPtr(rtm)    (&((rtm)->Timing.stopRequestedFlag))
#endif

#ifndef rtmGetT
#define rtmGetT(rtm)                   (rtmGetTPtr((rtm))[0])
#endif

#ifndef rtmGetTFinal
#define rtmGetTFinal(rtm)              ((rtm)->Timing.tFinal)
#endif

#ifndef rtmGetTPtr
#define rtmGetTPtr(rtm)                ((rtm)->Timing.t)
#endif

#ifndef rtmGetTStart
#define rtmGetTStart(rtm)              ((rtm)->Timing.tStart)
#endif

#ifndef rtmGetTimeOfLastOutput
#define rtmGetTimeOfLastOutput(rtm)    ((rtm)->Timing.timeOfLastOutput)
#endif

/* Block signals (default storage) */
struct B_Continues_motor_control_v9_tensiontest_T {
  real_T Step;                         /* '<S5>/Step' */
  real_T Clock;                        /* '<S5>/Clock' */
  real_T Sum;                          /* '<S5>/Sum' */
  real_T Product;                      /* '<S5>/Product' */
  real_T Output;                       /* '<S5>/Output' */
  real_T Step6;                        /* '<Root>/Step6' */
  real_T Product_a;                    /* '<Root>/Product' */
  real_T Subtract;                     /* '<Root>/Subtract' */
  real_T Potent;                       /* '<Root>/Analog input ' */
  real_T Laser;                        /* '<Root>/Analog input ' */
  real_T Voltage;                      /* '<Root>/Analog input ' */
  real_T Analoginput_o4;               /* '<Root>/Analog input ' */
  real_T Pott_mm;                      /* '<Root>/Gain3' */
  real_T Sum_p;                        /* '<Root>/Sum' */
  real_T Gain;                         /* '<S6>/Gain' */
  real_T Integrator;                   /* '<S6>/Integrator' */
  real_T Gain1;                        /* '<S6>/Gain1' */
  real_T Sum_j;                        /* '<S6>/Sum' */
  real_T Switch2;                      /* '<Root>/Switch2' */
  real_T Forward;                      /* '<Root>/Saturation' */
  real_T Switch3;                      /* '<Root>/Switch3' */
  real_T Backward;                     /* '<Root>/Saturation3' */
  real_T Loadcell;                     /* '<Root>/Transfer Fcn2' */
  real_T Gain1_j;                      /* '<Root>/Gain1' */
};

/* Block states (default storage) for system '<Root>' */
struct DW_Continues_motor_control_v9_tensiontest_T {
  real_T Analoginput_RWORK;            /* '<Root>/Analog input ' */
  void *Setup_PWORK;                   /* '<Root>/Setup ' */
  void *Analoginput_PWORK[3];          /* '<Root>/Analog input ' */
  void *Analogoutput_PWORK[3];         /* '<Root>/Analog output ' */
  struct {
    void *AQHandles;
    void *SLRTSigHandles;
  } TAQSigLogging_InsertedFor_Trans;   /* synthesized block */

  struct {
    void *AQHandles;
    void *SLRTSigHandles;
  } TAQSigLogging_InsertedFor_Analo;   /* synthesized block */

  struct {
    void *AQHandles;
    void *SLRTSigHandles;
  } TAQSigLogging_InsertedFor_Ana_a;   /* synthesized block */

  struct {
    void *AQHandles;
    void *SLRTSigHandles;
  } TAQSigLogging_InsertedFor_Gain3;   /* synthesized block */

  struct {
    void *LoggedData;
  } Scope_PWORK;                       /* '<Root>/Scope' */

  struct {
    void *LoggedData;
  } Scope1_PWORK;                      /* '<Root>/Scope1' */

  struct {
    void *LoggedData;
  } Scope2_PWORK;                      /* '<Root>/Scope2' */

  struct {
    void *LoggedData[3];
  } Scope3_PWORK;                      /* '<Root>/Scope3' */

  struct {
    void *LoggedData;
  } Scope4_PWORK;                      /* '<Root>/Scope4' */

  struct {
    void *LoggedData[3];
  } Scope5_PWORK;                      /* '<Root>/Scope5' */

  struct {
    void *LoggedData[2];
  } Scope6_PWORK;                      /* '<Root>/Scope6' */

  struct {
    void *LoggedData;
  } Scope8_PWORK;                      /* '<Root>/Scope8' */

  struct {
    void *LoggedData[2];
  } Scope9_PWORK;                      /* '<Root>/Scope9' */

  struct {
    void *AQHandles;
  } TAQSigLogging_InsertedFor_Ana_l;   /* synthesized block */

  int_T Analoginput_IWORK;             /* '<Root>/Analog input ' */
  int_T Analogoutput_IWORK[2];         /* '<Root>/Analog output ' */
};

/* Continuous states (default storage) */
struct X_Continues_motor_control_v9_tensiontest_T {
  real_T Integrator_CSTATE;            /* '<S6>/Integrator' */
};

/* State derivatives (default storage) */
struct XDot_Continues_motor_control_v9_tensiontest_T {
  real_T Integrator_CSTATE;            /* '<S6>/Integrator' */
};

/* State disabled  */
struct XDis_Continues_motor_control_v9_tensiontest_T {
  boolean_T Integrator_CSTATE;         /* '<S6>/Integrator' */
};

#ifndef ODE3_INTG
#define ODE3_INTG

/* ODE3 Integration Data */
struct ODE3_IntgData {
  real_T *y;                           /* output */
  real_T *f[3];                        /* derivatives */
};

#endif

/* Real-time Model Data Structure */
struct tag_RTM_Continues_motor_control_v9_tensiontest_T {
  struct SimStruct_tag * *childSfunctions;
  const char_T *errorStatus;
  SS_SimMode simMode;
  RTWSolverInfo solverInfo;
  RTWSolverInfo *solverInfoPtr;
  void *sfcnInfo;

  /*
   * NonInlinedSFcns:
   * The following substructure contains information regarding
   * non-inlined s-functions used in the model.
   */
  struct {
    RTWSfcnInfo sfcnInfo;
    time_T *taskTimePtrs[2];
    SimStruct childSFunctions[3];
    SimStruct *childSFunctionPtrs[3];
    struct _ssBlkInfo2 blkInfo2[3];
    struct _ssSFcnModelMethods2 methods2[3];
    struct _ssSFcnModelMethods3 methods3[3];
    struct _ssSFcnModelMethods4 methods4[3];
    struct _ssStatesInfo2 statesInfo2[3];
    ssPeriodicStatesInfo periodicStatesInfo[3];
    struct _ssPortInfo2 inputOutputPortInfo2[3];
    struct {
      time_T sfcnPeriod[1];
      time_T sfcnOffset[1];
      int_T sfcnTsMap[1];
      uint_T attribs[6];
      mxArray *params[6];
      struct _ssDWorkRecord dWork[1];
      struct _ssDWorkAuxRecord dWorkAux[1];
    } Sfcn0;

    struct {
      time_T sfcnPeriod[1];
      time_T sfcnOffset[1];
      int_T sfcnTsMap[1];
      struct _ssPortOutputs outputPortInfo[4];
      struct _ssOutPortUnit outputPortUnits[4];
      struct _ssOutPortCoSimAttribute outputPortCoSimAttribute[4];
      uint_T attribs[12];
      mxArray *params[12];
      struct _ssDWorkRecord dWork[3];
      struct _ssDWorkAuxRecord dWorkAux[3];
    } Sfcn1;

    struct {
      time_T sfcnPeriod[1];
      time_T sfcnOffset[1];
      int_T sfcnTsMap[1];
      struct _ssPortInputs inputPortInfo[3];
      struct _ssInPortUnit inputPortUnits[3];
      struct _ssInPortCoSimAttribute inputPortCoSimAttribute[3];
      uint_T attribs[15];
      mxArray *params[15];
      struct _ssDWorkRecord dWork[2];
      struct _ssDWorkAuxRecord dWorkAux[2];
    } Sfcn2;
  } NonInlinedSFcns;

  X_Continues_motor_control_v9_tensiontest_T *contStates;
  int_T *periodicContStateIndices;
  real_T *periodicContStateRanges;
  real_T *derivs;
  XDis_Continues_motor_control_v9_tensiontest_T *contStateDisabled;
  boolean_T zCCacheNeedsReset;
  boolean_T derivCacheNeedsReset;
  boolean_T CTOutputIncnstWithState;
  real_T odeY[1];
  real_T odeF[3][1];
  ODE3_IntgData intgData;

  /*
   * Sizes:
   * The following substructure contains sizes information
   * for many of the model attributes such as inputs, outputs,
   * dwork, sample times, etc.
   */
  struct {
    uint32_T options;
    int_T numContStates;
    int_T numPeriodicContStates;
    int_T numU;
    int_T numY;
    int_T numSampTimes;
    int_T numBlocks;
    int_T numBlockIO;
    int_T numBlockPrms;
    int_T numDwork;
    int_T numSFcnPrms;
    int_T numSFcns;
    int_T numIports;
    int_T numOports;
    int_T numNonSampZCs;
    int_T sysDirFeedThru;
    int_T rtwGenSfcn;
  } Sizes;

  /*
   * Timing:
   * The following substructure contains information regarding
   * the timing information for the model.
   */
  struct {
    time_T stepSize;
    uint32_T clockTick0;
    uint32_T clockTickH0;
    time_T stepSize0;
    uint32_T clockTick1;
    uint32_T clockTickH1;
    time_T stepSize1;
    time_T tStart;
    time_T tFinal;
    time_T timeOfLastOutput;
    SimTimeStep simTimeStep;
    boolean_T stopRequestedFlag;
    time_T *sampleTimes;
    time_T *offsetTimes;
    int_T *sampleTimeTaskIDPtr;
    int_T *sampleHits;
    int_T *perTaskSampleHits;
    time_T *t;
    time_T sampleTimesArray[2];
    time_T offsetTimesArray[2];
    int_T sampleTimeTaskIDArray[2];
    int_T sampleHitArray[2];
    int_T perTaskSampleHitsArray[4];
    time_T tArray[2];
  } Timing;
};

/* Block signals (default storage) */
#ifdef __cplusplus

extern "C"
{

#endif

  extern struct B_Continues_motor_control_v9_tensiontest_T
    Continues_motor_control_v9_tensiontest_B;

#ifdef __cplusplus

}

#endif

/* Continuous states (default storage) */
extern X_Continues_motor_control_v9_tensiontest_T
  Continues_motor_control_v9_tensiontest_X;

/* Disabled states (default storage) */
extern XDis_Continues_motor_control_v9_tensiontest_T
  Continues_motor_control_v9_tensiontest_XDis;

/* Block states (default storage) */
extern struct DW_Continues_motor_control_v9_tensiontest_T
  Continues_motor_control_v9_tensiontest_DW;

/* External data declarations for dependent source files */
extern const real_T Continues_motor_control_v9_tensiontest_RGND;/* real_T ground */

#ifdef __cplusplus

extern "C"
{

#endif

  /* Model entry point functions */
  extern void Continues_motor_control_v9_tensiontest_initialize(void);
  extern void Continues_motor_control_v9_tensiontest_step(void);
  extern void Continues_motor_control_v9_tensiontest_terminate(void);

#ifdef __cplusplus

}

#endif

/* Real-time Model object */
#ifdef __cplusplus

extern "C"
{

#endif

  extern RT_MODEL_Continues_motor_control_v9_tensiontest_T *const
    Continues_motor_control_v9_tensiontest_M;

#ifdef __cplusplus

}

#endif

/*-
 * The generated code includes comments that allow you to trace directly
 * back to the appropriate location in the model.  The basic format
 * is <system>/block_name, where system is the system number (uniquely
 * assigned by Simulink) and block_name is the name of the block.
 *
 * Use the MATLAB hilite_system command to trace the generated code back
 * to the model.  For example,
 *
 * hilite_system('<S3>')    - opens system 3
 * hilite_system('<S3>/Kp') - opens and selects block Kp which resides in S3
 *
 * Here is the system hierarchy for this model
 *
 * '<Root>' : 'Continues_motor_control_v9_tensiontest'
 * '<S1>'   : 'Continues_motor_control_v9_tensiontest/File Log1'
 * '<S2>'   : 'Continues_motor_control_v9_tensiontest/File Log2'
 * '<S3>'   : 'Continues_motor_control_v9_tensiontest/File Log3'
 * '<S4>'   : 'Continues_motor_control_v9_tensiontest/File Log4'
 * '<S5>'   : 'Continues_motor_control_v9_tensiontest/Ramp'
 * '<S6>'   : 'Continues_motor_control_v9_tensiontest/Subsystem'
 */
#endif                           /* Continues_motor_control_v9_tensiontest_h_ */
