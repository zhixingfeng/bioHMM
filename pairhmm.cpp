#include "pairhmm.h"
double PairHMM::cal_likelihood_from_cigar(bool is_log)
{
    if (cigar.size()==0)
        throw runtime_error("trying to calculate likelihood from cigar buy cigar is empty");
    // init state is forced to be M, i.e. Pr(S_0) = Pr(S_0|M)
    double log_L = 0;
    int i = -1;
    int j = -1;
    for (int k=0; k < (int) cigar.size(); k++){
        double log_prob_jump;
        if (k == 0)
            log_prob_jump = transProb('M', cigar[k].first, true);
        else
            log_prob_jump = transProb(cigar[k-1].first, cigar[k].first, true);
        if (cigar[k].first == 'M'){
            ++i; ++j;
            log_prob_jump += Pxy(seqX[i], seqY[j], true);
        }
        if (cigar[k].first == 'I'){
            ++i;
            log_prob_jump += Px(seqX[i], true);
        }
        if (cigar[k].first == 'D'){
            ++j;
            log_prob_jump += Py(seqY[j], true);
        }
        log_L += log_prob_jump;


        for (int n=1; n < cigar[k].second; n++){
            double log_prob_run = transProb(cigar[k].first, cigar[k].first, true);
            if (cigar[k].first == 'M'){
                ++i; ++j;
                log_prob_run += Pxy(seqX[i], seqY[j], true);
            }
            if (cigar[k].first == 'I'){
                ++i;
                log_prob_run += Px(seqX[i], true);
            }
            if (cigar[k].first == 'D'){
                ++j;
                log_prob_run += Py(seqY[j], true);
            }
            log_L += log_prob_run;
        }
    }
    if (is_log)
        return log_L;
    else
        return exp(log_L);
}

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
    scoreMat[0][0].log_Vm = 0; scoreMat[0][0].log_Vx = LOG_MIN; scoreMat[0][0].log_Vy = LOG_MIN;
    for (int i=1; i<scoreMat.nrow; i++){
        if (i==1)
            scoreMat[i][0].log_Vx = Px(seqX[i-1],true) + log_dlt_x + scoreMat[i-1][0].log_Vm;
        else
            scoreMat[i][0].log_Vx = Px(seqX[i-1],true) + log_eps_x + scoreMat[i-1][0].log_Vx;
        scoreMat[i][0].log_Vm = LOG_MIN;
        scoreMat[i][0].log_Vy = LOG_MIN;
        scoreMat[i][0].path_M = 'I';
        scoreMat[i][0].path_X = 'I';
        scoreMat[i][0].path_Y = 'I';
    }

    for (int j=1; j<scoreMat.ncol; j++){
        if (j==1)
            scoreMat[0][j].log_Vy = Py(seqY[j-1],true) + log_dlt_x + scoreMat[0][j-1].log_Vm;
        else
            scoreMat[0][j].log_Vy = Py(seqY[j-1],true) + log_eps_x + scoreMat[0][j-1].log_Vy;
        scoreMat[0][j].log_Vm = LOG_MIN;
        scoreMat[0][j].log_Vx = LOG_MIN;
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

    char prev_path;
    if (cur_path=='M'){
        prev_path = scoreMat[i][j].path_M[0];
        i--;
        j--;
    }
    if (cur_path=='I'){
        prev_path = scoreMat[i][j].path_X[0];
        i--;
    }
    if (cur_path=='D'){
        prev_path = scoreMat[i][j].path_Y[0];
        j--;
    }

    while (i>0 || j>0){
        if (cigar.back().first == prev_path)
            cigar.back().second ++;
        else
            cigar.push_back(pair<char,int>(prev_path,1));
        if (prev_path=='M'){
            prev_path = scoreMat[i][j].path_M[0];
            i--;
            j--;
            continue;
        }
        if (prev_path=='I'){
            prev_path = scoreMat[i][j].path_X[0];
            i--;
            continue;
        }
        if (prev_path=='D'){
            prev_path = scoreMat[i][j].path_Y[0];
            j--;
            continue;
        }
    }
    // reverse cigar
    reverse(cigar.begin(), cigar.end());
}

double PairHMM::forward()
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
    fwdMat[0][0].log_Vm = 0; fwdMat[0][0].log_Vx = LOG_MIN; fwdMat[0][0].log_Vy = LOG_MIN;
    for (int i=1; i<fwdMat.nrow; i++){
        if (i==1)
            fwdMat[i][0].log_Vx = Px(seqX[i-1],true) + log_dlt_x + fwdMat[i-1][0].log_Vm;
        else
            fwdMat[i][0].log_Vx = Px(seqX[i-1],true) + log_eps_x + fwdMat[i-1][0].log_Vx;
        fwdMat[i][0].log_Vm = LOG_MIN;
        fwdMat[i][0].log_Vy = LOG_MIN;
    }

    for (int j=1; j<fwdMat.ncol; j++){
        if (j==1)
            fwdMat[0][j].log_Vy = Py(seqY[j-1],true) + log_dlt_x + fwdMat[0][j-1].log_Vm;
        else
            fwdMat[0][j].log_Vy = Py(seqY[j-1],true) + log_eps_x + fwdMat[0][j-1].log_Vy;
        fwdMat[0][j].log_Vm = LOG_MIN;
        fwdMat[0][j].log_Vx = LOG_MIN;
    }

    // fill in fwdMatfwdMat
    for (int i=1; i<fwdMat.nrow; i++){
        for (int j=1; j<fwdMat.ncol; j++){
            // match
            fwdMat[i][j].log_Vm = log_sum(log_a_mm + fwdMat[i-1][j-1].log_Vm, log_a_im_x + fwdMat[i-1][j-1].log_Vx);
            fwdMat[i][j].log_Vm = log_sum(fwdMat[i][j].log_Vm, log_a_im_y + fwdMat[i-1][j-1].log_Vy);
            fwdMat[i][j].log_Vm += Pxy(seqX[i-1],seqY[j-1],true);

            // insertion on X
            fwdMat[i][j].log_Vx = log_sum(log_dlt_x + fwdMat[i-1][j].log_Vm, log_eps_x + fwdMat[i-1][j].log_Vx);
            fwdMat[i][j].log_Vx += Px(seqX[i-1],true);

            // insertion on Y
            fwdMat[i][j].log_Vy = log_sum(log_dlt_y + fwdMat[i][j-1].log_Vm, log_eps_y + fwdMat[i][j-1].log_Vy);
            fwdMat[i][j].log_Vy += Py(seqY[j-1],true);
        }
    }
    double rl = log_sum(fwdMat[fwdMat.nrow-1][fwdMat.ncol-1].log_Vm, fwdMat[fwdMat.nrow-1][fwdMat.ncol-1].log_Vx);
    rl = log_sum(rl, fwdMat[fwdMat.nrow-1][fwdMat.ncol-1].log_Vy);
    return rl;
}


