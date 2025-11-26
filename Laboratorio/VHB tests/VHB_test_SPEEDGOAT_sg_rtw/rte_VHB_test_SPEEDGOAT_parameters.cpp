#include "rte_VHB_test_SPEEDGOAT_parameters.h"
#include "VHB_test_SPEEDGOAT.h"
#include "VHB_test_SPEEDGOAT_cal.h"

RTE_Param_Service_T RTE_Param_Service = {
  1.0,
  15.0,
  17.0,
  27.425,
  4.0,
  80.0,
  0.5
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

real_T* get_slope(void)
{
  return &RTE_Param_Service_ptr->slope;
}

extern VHB_test_SPEEDGOAT_cal_type VHB_test_SPEEDGOAT_cal_impl;
extern RTE_Param_Service_T RTE_Param_Service;
namespace slrealtime
{
  /* Description of SEGMENTS */
  SegmentVector segmentInfo {
    { (void*)&RTE_Param_Service, (void**)&RTE_Param_Service_ptr, sizeof
      (RTE_Param_Service_T), 2 },

    { (void*)&VHB_test_SPEEDGOAT_cal_impl, (void**)&VHB_test_SPEEDGOAT_cal,
      sizeof(VHB_test_SPEEDGOAT_cal_type), 2 }
  };

  SegmentVector &getSegmentVector(void)
  {
    return segmentInfo;
  }
}                                      // slrealtime
