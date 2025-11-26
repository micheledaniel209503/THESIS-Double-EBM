/* Main generated for Simulink Real-Time model Continues_motor_control_v9_tensiontest */
#include <ModelInfo.hpp>
#include <utilities.hpp>
#include "Continues_motor_control_v9_tensiontest.h"
#include "rte_Continues_motor_control_v9_tensiontest_parameters.h"

/* Task wrapper function definitions */
void Continues_motor_control_v9_tensiontest_Task1(void)
{ 
    Continues_motor_control_v9_tensiontest_step();
} 
/* Task descriptors */
slrealtime::TaskInfo task_1( 0u, std::bind(Continues_motor_control_v9_tensiontest_Task1), slrealtime::TaskInfo::PERIODIC, 0.0001, 0, 40);

/* Executable base address for XCP */
#ifdef __linux__
extern char __executable_start;
static uintptr_t const base_address = reinterpret_cast<uintptr_t>(&__executable_start);
#else
/* Set 0 as placeholder, to be parsed later from /proc filesystem */
static uintptr_t const base_address = 0;
#endif

/* Model descriptor */
slrealtime::ModelInfo Continues_motor_control_v9_tensiontest_Info =
{
    "Continues_motor_control_v9_tensiontest",
    Continues_motor_control_v9_tensiontest_initialize,
    Continues_motor_control_v9_tensiontest_terminate,
    []()->char const*& { return Continues_motor_control_v9_tensiontest_M->errorStatus; },
    []()->unsigned char& { return Continues_motor_control_v9_tensiontest_M->Timing.stopRequestedFlag; },
    { task_1 },
    slrealtime::getSegmentVector()
};

int main(int argc, char *argv[]) {
    slrealtime::BaseAddress::set(base_address);
    return slrealtime::runModel(argc, argv, Continues_motor_control_v9_tensiontest_Info);
}
