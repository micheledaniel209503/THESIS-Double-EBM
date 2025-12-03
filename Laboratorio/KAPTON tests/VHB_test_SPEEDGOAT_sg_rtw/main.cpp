/* Main generated for Simulink Real-Time model VHB_test_SPEEDGOAT */
#include <ModelInfo.hpp>
#include <utilities.hpp>
#include "VHB_test_SPEEDGOAT.h"
#include "rte_VHB_test_SPEEDGOAT_parameters.h"

/* Task wrapper function definitions */
void VHB_test_SPEEDGOAT_Task1(void)
{ 
    VHB_test_SPEEDGOAT_step();
} 
/* Task descriptors */
slrealtime::TaskInfo task_1( 0u, std::bind(VHB_test_SPEEDGOAT_Task1), slrealtime::TaskInfo::PERIODIC, 0.0001, 0, 40);

/* Executable base address for XCP */
#ifdef __linux__
extern char __executable_start;
static uintptr_t const base_address = reinterpret_cast<uintptr_t>(&__executable_start);
#else
/* Set 0 as placeholder, to be parsed later from /proc filesystem */
static uintptr_t const base_address = 0;
#endif

/* Model descriptor */
slrealtime::ModelInfo VHB_test_SPEEDGOAT_Info =
{
    "VHB_test_SPEEDGOAT",
    VHB_test_SPEEDGOAT_initialize,
    VHB_test_SPEEDGOAT_terminate,
    []()->char const*& { return VHB_test_SPEEDGOAT_M->errorStatus; },
    []()->unsigned char& { return VHB_test_SPEEDGOAT_M->Timing.stopRequestedFlag; },
    { task_1 },
    slrealtime::getSegmentVector()
};

int main(int argc, char *argv[]) {
    slrealtime::BaseAddress::set(base_address);
    return slrealtime::runModel(argc, argv, VHB_test_SPEEDGOAT_Info);
}
