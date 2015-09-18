#include "pairhmm.h"

void PairHMM::viterbi()
{
    if (!is_setSeq)
        throw runtime_error("sequences have not be set.");
    if (!is_setPar)
        throw runtime_error("parameters have not be set.");
    if (!is_set_Px)
        throw runtime_error("Px have not be set.");
    if (!is_set_Py)
        throw runtime_error("Py have not be set.");
    if (!is_set_Pxy)
        throw runtime_error("Pxy have not be set.");

    // start to fill out scoreMat
    for (int i=0; i<scoreMat.nrow; i++){
        for (int j=0; j<scoreMat.ncol; j++){
            if (i==0 && j==0){
                // init
                scoreMat[i][j].log_Vm = 0;
                scoreMat[i][j].log_Vx = log(TOL);
                scoreMat[i][j].log_Vy = log(TOL);
            }else{
                // init first row
                if (i==0){
                    scoreMat[i][j].log_Vm = log(TOL);
                    scoreMat[i][j].log_Vx = log(TOL);
                    if (j==1)
                        scoreMat[i][j].log_Vy = Py(seqY[j],true) + log_dlt_y + scoreMat[i][j-1].log_Vm;
                    else
                        scoreMat[i][j].log_Vy = Py(seqY[j],true) + log_eps_y + scoreMat[i][j-1].log_Vy;
                    scoreMat[i][j].path = "D";
                    continue;
                }
                // init first col
                if (j==0){
                    scoreMat[i][j].log_Vm = log(TOL);
                    scoreMat[i][j].log_Vy = log(TOL);
                    if (i==1)
                        scoreMat[i][j].log_Vx = Px(seqX[i],true) + log_dlt_x + scoreMat[i-1][j].log_Vm;
                    else
                        scoreMat[i][j].log_Vx = Px(seqX[i],true) + log_eps_x + scoreMat[i-1][j].log_Vx;
                    scoreMat[i][j].path = "I";
                    continue;
                }
                // fill normal cells
                scoreMat[i][j].log_Vm = log(TOL);
                scoreMat[i][j].log_Vx = log(TOL);
                scoreMat[i][j].log_Vy = log(TOL);


            }

        }
    }
}
