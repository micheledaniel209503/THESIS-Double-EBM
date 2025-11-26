#ifndef Continues_motor_control_v9_tensiontest_cal_h_
#define Continues_motor_control_v9_tensiontest_cal_h_
#include "rtwtypes.h"

/* Storage class 'PageSwitching', for system '<Root>' */
struct Continues_motor_contro_cal_type {
  real_T Ramp_InitialOutput;           /* Mask Parameter: Ramp_InitialOutput
                                        * Referenced by: '<S5>/Constant1'
                                        */
  real_T Ramp_slope;                   /* Mask Parameter: Ramp_slope
                                        * Referenced by: '<S5>/Step'
                                        */
  real_T Gain1_Gain;                   /* Expression: -1
                                        * Referenced by: '<Root>/Gain1'
                                        */
  real_T Setup_P1_Size[2];             /* Computed Parameter: Setup_P1_Size
                                        * Referenced by: '<Root>/Setup '
                                        */
  real_T Setup_P1;                     /* Expression: parPciSlot
                                        * Referenced by: '<Root>/Setup '
                                        */
  real_T Setup_P2_Size[2];             /* Computed Parameter: Setup_P2_Size
                                        * Referenced by: '<Root>/Setup '
                                        */
  real_T Setup_P2;                     /* Expression: parModuleId
                                        * Referenced by: '<Root>/Setup '
                                        */
  real_T Setup_P3_Size[2];             /* Computed Parameter: Setup_P3_Size
                                        * Referenced by: '<Root>/Setup '
                                        */
  real_T Setup_P3;                     /* Expression: parBoardType
                                        * Referenced by: '<Root>/Setup '
                                        */
  real_T Setup_P4_Size[2];             /* Computed Parameter: Setup_P4_Size
                                        * Referenced by: '<Root>/Setup '
                                        */
  real_T Setup_P4;                     /* Expression: parTriggerInitiator
                                        * Referenced by: '<Root>/Setup '
                                        */
  real_T Setup_P5_Size[2];             /* Computed Parameter: Setup_P5_Size
                                        * Referenced by: '<Root>/Setup '
                                        */
  real_T Setup_P5;                     /* Expression: parAdDmaEn
                                        * Referenced by: '<Root>/Setup '
                                        */
  real_T Setup_P6_Size[2];             /* Computed Parameter: Setup_P6_Size
                                        * Referenced by: '<Root>/Setup '
                                        */
  real_T Setup_P6;                     /* Expression: parDaDmaEn
                                        * Referenced by: '<Root>/Setup '
                                        */
  real_T Step_Y0;                      /* Expression: 0
                                        * Referenced by: '<S5>/Step'
                                        */
  real_T Step6_Y0;                     /* Expression: 1
                                        * Referenced by: '<Root>/Step6'
                                        */
  real_T Step6_YFinal;                 /* Expression: 0
                                        * Referenced by: '<Root>/Step6'
                                        */
  real_T Analoginput_P1_Size[2];      /* Computed Parameter: Analoginput_P1_Size
                                       * Referenced by: '<Root>/Analog input '
                                       */
  real_T Analoginput_P1;               /* Expression: parSampleTime
                                        * Referenced by: '<Root>/Analog input '
                                        */
  real_T Analoginput_P2_Size[2];      /* Computed Parameter: Analoginput_P2_Size
                                       * Referenced by: '<Root>/Analog input '
                                       */
  real_T Analoginput_P2;               /* Expression: parPciSlot
                                        * Referenced by: '<Root>/Analog input '
                                        */
  real_T Analoginput_P3_Size[2];      /* Computed Parameter: Analoginput_P3_Size
                                       * Referenced by: '<Root>/Analog input '
                                       */
  real_T Analoginput_P3;               /* Expression: parModuleId
                                        * Referenced by: '<Root>/Analog input '
                                        */
  real_T Analoginput_P4_Size[2];      /* Computed Parameter: Analoginput_P4_Size
                                       * Referenced by: '<Root>/Analog input '
                                       */
  real_T Analoginput_P4;               /* Expression: parBoardType
                                        * Referenced by: '<Root>/Analog input '
                                        */
  real_T Analoginput_P5_Size[2];      /* Computed Parameter: Analoginput_P5_Size
                                       * Referenced by: '<Root>/Analog input '
                                       */
  real_T Analoginput_P5;               /* Expression: parRange
                                        * Referenced by: '<Root>/Analog input '
                                        */
  real_T Analoginput_P6_Size[2];      /* Computed Parameter: Analoginput_P6_Size
                                       * Referenced by: '<Root>/Analog input '
                                       */
  real_T Analoginput_P6;               /* Expression: parOversampling
                                        * Referenced by: '<Root>/Analog input '
                                        */
  real_T Analoginput_P7_Size[2];      /* Computed Parameter: Analoginput_P7_Size
                                       * Referenced by: '<Root>/Analog input '
                                       */
  real_T Analoginput_P7[4];            /* Expression: parChannels
                                        * Referenced by: '<Root>/Analog input '
                                        */
  real_T Analoginput_P8_Size[2];      /* Computed Parameter: Analoginput_P8_Size
                                       * Referenced by: '<Root>/Analog input '
                                       */
  real_T Analoginput_P8;               /* Expression: parContinuousSampling
                                        * Referenced by: '<Root>/Analog input '
                                        */
  real_T Analoginput_P9_Size[2];      /* Computed Parameter: Analoginput_P9_Size
                                       * Referenced by: '<Root>/Analog input '
                                       */
  real_T Analoginput_P9;               /* Expression: parDmaEn
                                        * Referenced by: '<Root>/Analog input '
                                        */
  real_T Analoginput_P10_Size[2];    /* Computed Parameter: Analoginput_P10_Size
                                      * Referenced by: '<Root>/Analog input '
                                      */
  real_T Analoginput_P10;              /* Expression: parDmaChannels
                                        * Referenced by: '<Root>/Analog input '
                                        */
  real_T Analoginput_P11_Size[2];    /* Computed Parameter: Analoginput_P11_Size
                                      * Referenced by: '<Root>/Analog input '
                                      */
  real_T Analoginput_P11;              /* Expression: parDmaConversionTime
                                        * Referenced by: '<Root>/Analog input '
                                        */
  real_T Analoginput_P12_Size[2];    /* Computed Parameter: Analoginput_P12_Size
                                      * Referenced by: '<Root>/Analog input '
                                      */
  real_T Analoginput_P12;              /* Expression: parDmaFrameSize
                                        * Referenced by: '<Root>/Analog input '
                                        */
  real_T Integrator_IC;                /* Expression: 0
                                        * Referenced by: '<S6>/Integrator'
                                        */
  real_T Constant5_Value;              /* Expression: 0
                                        * Referenced by: '<Root>/Constant5'
                                        */
  real_T Switch2_Threshold;            /* Expression: 0
                                        * Referenced by: '<Root>/Switch2'
                                        */
  real_T Saturation_UpperSat;          /* Expression: 5
                                        * Referenced by: '<Root>/Saturation'
                                        */
  real_T Saturation_LowerSat;          /* Expression: 0
                                        * Referenced by: '<Root>/Saturation'
                                        */
  real_T Constant4_Value;              /* Expression: 0
                                        * Referenced by: '<Root>/Constant4'
                                        */
  real_T Switch3_Threshold;            /* Expression: 0
                                        * Referenced by: '<Root>/Switch3'
                                        */
  real_T Saturation3_UpperSat;         /* Expression: 5
                                        * Referenced by: '<Root>/Saturation3'
                                        */
  real_T Saturation3_LowerSat;         /* Expression: 0
                                        * Referenced by: '<Root>/Saturation3'
                                        */
  real_T Analogoutput_P1_Size[2];    /* Computed Parameter: Analogoutput_P1_Size
                                      * Referenced by: '<Root>/Analog output '
                                      */
  real_T Analogoutput_P1;              /* Expression: parSampleTime
                                        * Referenced by: '<Root>/Analog output '
                                        */
  real_T Analogoutput_P2_Size[2];    /* Computed Parameter: Analogoutput_P2_Size
                                      * Referenced by: '<Root>/Analog output '
                                      */
  real_T Analogoutput_P2;              /* Expression: parPciSlot
                                        * Referenced by: '<Root>/Analog output '
                                        */
  real_T Analogoutput_P3_Size[2];    /* Computed Parameter: Analogoutput_P3_Size
                                      * Referenced by: '<Root>/Analog output '
                                      */
  real_T Analogoutput_P3;              /* Expression: parModuleId
                                        * Referenced by: '<Root>/Analog output '
                                        */
  real_T Analogoutput_P4_Size[2];    /* Computed Parameter: Analogoutput_P4_Size
                                      * Referenced by: '<Root>/Analog output '
                                      */
  real_T Analogoutput_P4;              /* Expression: parBoardType
                                        * Referenced by: '<Root>/Analog output '
                                        */
  real_T Analogoutput_P5_Size[2];    /* Computed Parameter: Analogoutput_P5_Size
                                      * Referenced by: '<Root>/Analog output '
                                      */
  real_T Analogoutput_P5;              /* Expression: parRange
                                        * Referenced by: '<Root>/Analog output '
                                        */
  real_T Analogoutput_P6_Size[2];    /* Computed Parameter: Analogoutput_P6_Size
                                      * Referenced by: '<Root>/Analog output '
                                      */
  real_T Analogoutput_P6[3];           /* Expression: parChannels
                                        * Referenced by: '<Root>/Analog output '
                                        */
  real_T Analogoutput_P7_Size[2];    /* Computed Parameter: Analogoutput_P7_Size
                                      * Referenced by: '<Root>/Analog output '
                                      */
  real_T Analogoutput_P7[3];           /* Expression: parInitValues
                                        * Referenced by: '<Root>/Analog output '
                                        */
  real_T Analogoutput_P8_Size[2];    /* Computed Parameter: Analogoutput_P8_Size
                                      * Referenced by: '<Root>/Analog output '
                                      */
  real_T Analogoutput_P8[3];           /* Expression: parResetValues
                                        * Referenced by: '<Root>/Analog output '
                                        */
  real_T Analogoutput_P9_Size[2];    /* Computed Parameter: Analogoutput_P9_Size
                                      * Referenced by: '<Root>/Analog output '
                                      */
  real_T Analogoutput_P9;              /* Expression: parDmaChannels
                                        * Referenced by: '<Root>/Analog output '
                                        */
  real_T Analogoutput_P10_Size[2];  /* Computed Parameter: Analogoutput_P10_Size
                                     * Referenced by: '<Root>/Analog output '
                                     */
  real_T Analogoutput_P10;             /* Expression: parDmaEn
                                        * Referenced by: '<Root>/Analog output '
                                        */
  real_T Analogoutput_P11_Size[2];  /* Computed Parameter: Analogoutput_P11_Size
                                     * Referenced by: '<Root>/Analog output '
                                     */
  real_T Analogoutput_P11;             /* Expression: parDmaClkSrc
                                        * Referenced by: '<Root>/Analog output '
                                        */
  real_T Analogoutput_P12_Size[2];  /* Computed Parameter: Analogoutput_P12_Size
                                     * Referenced by: '<Root>/Analog output '
                                     */
  real_T Analogoutput_P12;             /* Expression: parDmaConversionTime
                                        * Referenced by: '<Root>/Analog output '
                                        */
  real_T Analogoutput_P13_Size[2];  /* Computed Parameter: Analogoutput_P13_Size
                                     * Referenced by: '<Root>/Analog output '
                                     */
  real_T Analogoutput_P13;             /* Expression: parDmaFrameSize
                                        * Referenced by: '<Root>/Analog output '
                                        */
  real_T Analogoutput_P14_Size[2];  /* Computed Parameter: Analogoutput_P14_Size
                                     * Referenced by: '<Root>/Analog output '
                                     */
  real_T Analogoutput_P14;             /* Expression: parDoInterrupt
                                        * Referenced by: '<Root>/Analog output '
                                        */
  real_T Analogoutput_P15_Size[2];  /* Computed Parameter: Analogoutput_P15_Size
                                     * Referenced by: '<Root>/Analog output '
                                     */
  real_T Analogoutput_P15;             /* Expression: parLatency
                                        * Referenced by: '<Root>/Analog output '
                                        */
  real_T TransferFcn2_D;               /* Computed Parameter: TransferFcn2_D
                                        * Referenced by: '<Root>/Transfer Fcn2'
                                        */
};

/* Storage class 'PageSwitching' */
extern Continues_motor_contro_cal_type Continues_motor_contro_cal_impl;
extern Continues_motor_contro_cal_type *Continues_motor_control_v9__cal;

#endif                       /* Continues_motor_control_v9_tensiontest_cal_h_ */
