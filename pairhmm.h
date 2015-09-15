#ifndef PAIRHMM
#define PAIRHMM

#include "stl.h"
struct scoreCell
{
    // current score of match/mismatch, gap in X, gap in Y
    double Vm, Vx, Vy;
    char base;

};

class pairHMM
{
    public:

    protected:
        string seqX;
        string seqY;
        vector < pair<char, int> > cigar;

};



#endif // PAIRHMM

