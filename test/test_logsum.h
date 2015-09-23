#ifndef TEST_LOGSUM
#define TEST_LOGSUM
#include "../logsum.h"
#include <UnitTest++/UnitTest++.h>

TEST(test_logsum)
{
    for (double log_p=-10; log_p<=10; log_p+=0.1)
        for (double log_q=-10; log_q<=10; log_q+=0.1)
            CHECK_CLOSE(log_sum(log_p, log_q), log(exp(log_p) + exp(log_q)), 1e-12);
}


#endif // TEST_LOGSUM

