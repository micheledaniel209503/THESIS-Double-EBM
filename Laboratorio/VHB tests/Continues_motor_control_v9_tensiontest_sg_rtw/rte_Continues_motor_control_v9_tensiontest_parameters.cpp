#include "rte_Continues_motor_control_v9_tensiontest_parameters.h"
#include "Continues_motor_control_v9_tensiontest.h"
#include "Continues_motor_control_v9_tensiontest_cal.h"

RTE_Param_Service_T RTE_Param_Service = {
  1.0,
  15.0,
  11.0,
  27.425,
  4.0,
  20.0
};

RTE_Param_Service_T *RTE_Param_Service_ptr = &RTE_Param_Service;
real_T* get_Ki(void)
{
  return &RTE_Param_Service_ptr->Ki;
}

real_T* get_Kp(void)
{
  return &RTE_Param_Service_ptr->Kp;
}

real_T* get_N_finish(void)
{
  return &RTE_Param_Service_ptr->N_finish;
}

real_T* get_Potent_gain(void)
{
  return &RTE_Param_Service_ptr->Potent_gain;
}

real_T* get_T(void)
{
  return &RTE_Param_Service_ptr->T;
}

real_T* get_bias(void)
{
  return &RTE_Param_Service_ptr->bias;
}

extern Continues_motor_contro_cal_type Continues_motor_contro_cal_impl;
extern RTE_Param_Service_T RTE_Param_Service;
namespace slrealtime
{
  /* Description of SEGMENTS */
  SegmentVector segmentInfo {
    { (void*)&RTE_Param_Service, (void**)&RTE_Param_Service_ptr, sizeof
      (RTE_Param_Service_T), 2 },

    { (void*)&Continues_motor_contro_cal_impl, (void**)
      &Continues_motor_control_v9__cal, sizeof(Continues_motor_contro_cal_type),
      2 }
  };

  SegmentVector &getSegmentVector(void)
  {
    return segmentInfo;
  }
}                                      // slrealtime
