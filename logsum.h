#ifndef LOGSUM
#define LOGSUM
#include <math.h>
// summation of probability in log space (ineffient implementation, should be replace by hash based fast algorithm)
inline double log_sum(double log_p, double log_q){
    if (log_p >= log_q)
        return log_p + log(1+exp(log_q - log_p));
    else
        return log_q + log(1+exp(log_p - log_q));
}
#endif // LOGSUM

