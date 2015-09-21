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

    // initialize
    scoreMat[0][0].log_Vm = 0; scoreMat[0][0].log_Vx = log(TOL); scoreMat[0][0].log_Vy = log(TOL);
    for (int i=1; i<scoreMat.nrow; i++){
        if (i==1)
            scoreMat[i][0].log_Vx = Px(seqX[i-1],true) + log_dlt_x + scoreMat[i-1][0].log_Vm;
        else
            scoreMat[i][0].log_Vx = Px(seqX[i-1],true) + log_eps_x + scoreMat[i-1][0].log_Vx;
        scoreMat[i][0].log_Vm = log(TOL);
        scoreMat[i][0].log_Vy = log(TOL);
        scoreMat[i][0].path_M = 'I';
        scoreMat[i][0].path_X = 'I';
        scoreMat[i][0].path_Y = 'I';
    }

    for (int j=1; j<scoreMat.ncol; j++){
        if (j==1)
            scoreMat[0][j].log_Vy = Py(seqY[j-1],true) + log_dlt_x + scoreMat[0][j-1].log_Vm;
        else
            scoreMat[0][j].log_Vy = Py(seqY[j-1],true) + log_eps_x + scoreMat[0][j-1].log_Vy;
        scoreMat[0][j].log_Vm = log(TOL);
        scoreMat[0][j].log_Vx = log(TOL);
        scoreMat[0][j].path_M = 'D';
        scoreMat[0][j].path_X = 'D';
        scoreMat[0][j].path_Y = 'D';
    }
    // start to fill out scoreMat
    for (int i=1; i<scoreMat.nrow; i++){
        for (int j=1; j<scoreMat.ncol; j++){
            // match
            double log_Vm_mm = log_a_mm + scoreMat[i-1][j-1].log_Vm;
            double log_Vm_xm = log_a_im_x + scoreMat[i-1][j-1].log_Vx;
            double log_Vm_ym = log_a_im_y + scoreMat[i-1][j-1].log_Vy;

            scoreMat[i][j].log_Vm = log_Vm_mm;
            scoreMat[i][j].path_M = "M";
            if (log_Vm_xm > scoreMat[i][j].log_Vm){
                scoreMat[i][j].log_Vm = log_Vm_xm;
                scoreMat[i][j].path_M = "I";
            }
            if (log_Vm_ym > scoreMat[i][j].log_Vm){
                scoreMat[i][j].log_Vm = log_Vm_ym;
                scoreMat[i][j].path_M = "D";
            }

            scoreMat[i][j].log_Vm += Pxy(seqX[i-1],seqY[j-1],true);

            // insertion in X
            double log_Vx_mx = log_dlt_x + scoreMat[i-1][j].log_Vm;
            double log_Vx_xx = log_eps_x + scoreMat[i-1][j].log_Vx;

            scoreMat[i][j].log_Vx = log_Vx_mx;
            scoreMat[i][j].path_X = "M";
            if (log_Vx_xx > scoreMat[i][j].log_Vx){
                scoreMat[i][j].log_Vx = log_Vx_xx;
                scoreMat[i][j].path_X = "I";
            }

            scoreMat[i][j].log_Vx += Px(seqX[i-1],true);

            // insertion in Y
            double log_Vy_my = log_dlt_y + scoreMat[i][j-1].log_Vm;
            double log_Vy_yy = log_eps_y + scoreMat[i][j-1].log_Vy;

            scoreMat[i][j].log_Vy = log_Vy_my;
            scoreMat[i][j].path_Y = "M";
            if (log_Vy_yy > scoreMat[i][j].log_Vy){
                scoreMat[i][j].log_Vy = log_Vy_yy;
                scoreMat[i][j].path_Y = "D";
            }
        }
    }
    // trace back
    int i = scoreMat.nrow - 1;
    int j = scoreMat.ncol - 1;
    double cur_max_log_V = scoreMat[i][j].log_Vm;
    char cur_path = 'M';
    if (scoreMat[i][j].log_Vx > cur_max_log_V){
        cur_max_log_V = scoreMat[i][j].log_Vx;
        cur_path = 'I';
    }
    if (scoreMat[i][j].log_Vy > cur_max_log_V){
        cur_max_log_V = scoreMat[i][j].log_Vy;
        cur_path = 'D';
    }
    cigar.push_back(pair<char,int>(cur_path,1));
    if (cur_path=='M'){
        i--;
        j--;
    }
    if (cur_path=='I')
        i--;
    if (cur_path=='D')
        j--;
    char prev_path = cur_path;
    while (i>0 || j>0){
        if (prev_path=='M'){
            cigar.push_back(pair<char,int>(scoreMat[i][j].path_M[0],1));
            prev_path = scoreMat[i][j].path_M[0];
            i--;
            j--;
        }
        if (prev_path=='I'){
            cigar.push_back(pair<char,int>(scoreMat[i][j].path_X[0],1));
            prev_path = scoreMat[i][j].path_X[0];
            i--;
        }
        if (prev_path=='D'){
            cigar.push_back(pair<char,int>(scoreMat[i][j].path_Y[0],1));
            prev_path = scoreMat[i][j].path_Y[0];
            j--;
        }
    }

}
