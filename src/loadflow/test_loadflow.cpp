#include <iostream>
#include <iterator>
#include <algorithm>
#include "loadflow_baseline.h"

extern "C" {
 int loadflow(const double *bus, const double *line, const unsigned int noofgens,
              const unsigned int SB,
              const unsigned int noofslack, double *bus_sol);
}

#define PASS "PASS \n"
#define FAIL "\n             *** FAIL ***\n\n"
#define TEST(a) (a ? std::cout << PASS : std::cout << FAIL)

#define TOLERANCE 0.00000001
#define TEST_EQ(A, B) {\
                   TEST(std::equal(std::begin(A), std::end(A), std::begin(B),  \
                   [](double a, double b) -> bool { \
                        return ((a-b) < TOLERANCE) && ((a-b) > -TOLERANCE); \
                   }));}

int main () {
    double bus_sol[390];
    //unsigned int out = loadflow(bus, line, noofgens, SB, noofslack, bus_sol);
    loadflow(bus, line_full, 7, 1, 1, bus_sol);

    // Test the baseline; uncomment the following
//    bus_sol[389] = 4;

    TEST_EQ(bus_sol, bus_sol2)
}
