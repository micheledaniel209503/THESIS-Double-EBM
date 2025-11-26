#ifndef _RTE_CONTINUES_MOTOR_CONTROL_V__TENSIONTEST_PARAMETERS_H
#define _RTE_CONTINUES_MOTOR_CONTROL_V__TENSIONTEST_PARAMETERS_H
#include "rtwtypes.h"
#include "SegmentInfo.hpp"
#include "multiword_types.h"
#include "Continues_motor_control_v9_tensiontest_types.h"

struct RTE_Param_Service_T {
  real_T Ki;
  real_T Kp;
  real_T N_finish;
  real_T Potent_gain;
  real_T T;
  real_T bias;
};

extern RTE_Param_Service_T RTE_Param_Service;
extern RTE_Param_Service_T *RTE_Param_Service_ptr;
real_T* get_Ki(void);
real_T* get_Kp(void);
real_T* get_N_finish(void);
real_T* get_Potent_gain(void);
real_T* get_T(void);
real_T* get_bias(void);
namespace slrealtime
{
  SegmentVector &getSegmentVector(void);
}                                      // slrealtime

#endif
