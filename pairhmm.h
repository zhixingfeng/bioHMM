#ifndef PAIRHMM_H
#define PAIRHMM_H

#include "stl.h"
#include "matrix.h"
#include "math_utils.h"
#include "logsum.h"
#define TOL 1e-16
#define LOG_MIN -1e16

struct QualityScore
{
    QualityScore()
    {
        M = -1;
        I = -1;
        D = -1;
        G = -1;
    }
    QualityScore(int _M, int _I, int _D, int _G)
    {
        M = _M;
        I = _I;
        D = _D;
        G = _G;
    }
    QualityScore(double _M, double _I, double _D, double _G)
    {
        M = _M;
        I = _I;
        D = _D;
        G = _G;
    }
    double M;
    double I;
    double D;
    double G;
};

struct ScoreCell
{
    // init
    ScoreCell()
    {
        bases_X = 'N';
        bases_Y = 'N';
    }

    // current score of match/mismatch, gap in X, gap in Y
    double Vm, Vx, Vy, log_Vm, log_Vx, log_Vy;

    // matched bases in the current cell, for example <A,A> or <-,G>
    char bases_X;
    char bases_Y;

    // quality score
    QualityScore QV;

    // path
    string path_M;
    string path_X;
    string path_Y;
};

typedef Matrix<ScoreCell> ScoreMatrix;




class PairHMM
{
    public:
        PairHMM()
        {
            is_setSeq = false;
            is_setPar = false;
            is_set_Px = false;
            is_set_Py = false;
            is_set_Pxy = false;
        }
        // --------------- parameters --------------//
        // set sequences
        inline void setSeq(string _seqX, string _seqY, bool is_val = true)
        {
            seqX = _seqX;
            seqY = _seqY;
            if (is_val){
                validateSeq(seqX);
                validateSeq(seqY);
            }
            scoreMat.setDim(seqX.size()+1, seqY.size()+1);
            fwdMat.setDim(seqX.size()+1, seqY.size()+1);
            bwdMat.setDim(seqX.size()+1, seqY.size()+1);
            is_setSeq = true;
        }

        // get sequence
        inline string getSeqX(){return seqX;}
        inline string getSeqY(){return seqY;}

        // set parameters eps and dlt
        inline void setPar(double _eps_x, double _dlt_x, double _eps_y, double _dlt_y)
        {
            if (_eps_x <= 0 || _eps_x >= 1 || _dlt_x <= 0 || _dlt_x >= 1 || _eps_y <= 0 || _eps_y >= 1 || _dlt_y <= 0 || _dlt_y >= 1){
                throw runtime_error("in setPar, eps and dlt should be in (0,1)");
            }
            if (_dlt_x + _dlt_y >= 1)
                throw runtime_error("in setPar, dlt_x + dlt_y should < 1");

            // set parameters
            eps_x = _eps_x; eps_y = _eps_y;
            dlt_x = _dlt_x; dlt_y = _dlt_y;
            log_eps_x = log(eps_x); log_eps_y = log(eps_y);
            log_dlt_x = log(dlt_x); log_dlt_y = log(dlt_y);

            a_mm = 1 - dlt_x - dlt_y;
            a_im_x = 1 - eps_x;
            a_im_y = 1 - eps_y;
            log_a_mm = log(a_mm);
            log_a_im_x = log(a_im_x);
            log_a_im_y = log(a_im_y);

            // set transition probability, M = 0, I/X = 1, D/Y = 2
            transProb_value.setDim(3,3);
            transProb_value[0][0] = 1 - dlt_x - dlt_y; transProb_value[0][1] = dlt_x; transProb_value[0][2] = dlt_y;
            transProb_value[1][0] = 1 - eps_x; transProb_value[1][1] = eps_x; transProb_value[1][2] = 0;
            transProb_value[2][0] = 1 - eps_y; transProb_value[2][1] = 0; transProb_value[2][2] = eps_y;
            log_transProb_value.setDim(transProb_value.nrow, transProb_value.ncol);
            for (int i=0; i < transProb_value.nrow; i++)
                for (int j=0; j < transProb_value.ncol; j++)
                    log_transProb_value[i][j] = transProb_value[i][j] > TOL ? log(transProb_value[i][j]) : LOG_MIN;

            is_setPar = true;
        }

        // get parameters
        inline double getPar_eps_x(){return eps_x;}
        inline double getPar_dlt_x(){return dlt_x;}
        inline double getPar_log_eps_x(){return log_eps_x;}
        inline double getPar_log_dlt_x(){return log_dlt_x;}

        inline double getPar_eps_y(){return eps_y;}
        inline double getPar_dlt_y(){return dlt_y;}
        inline double getPar_log_eps_y(){return log_eps_y;}
        inline double getPar_log_dlt_y(){return log_dlt_y;}

        inline double getPar_a_mm(){return a_mm;}
        inline double getPar_a_im_x(){return a_im_x;}
        inline double getPar_a_im_y(){return a_im_y;}
        inline double getPar_log_a_mm(){return log_a_mm;}
        inline double getPar_log_a_im_x(){return log_a_im_x;}
        inline double getPar_log_a_im_y(){return log_a_im_y;}

        // set proablity Px, Py, and Pxy
        inline void set_Px(vector <double> & _Px_value)
        {
            if (_Px_value.size() != 4 && _Px_value.size() != 5)
                throw runtime_error("in PairHMM::set_Px, size of Px_value should be 4 or 5.");

            for (int i=0; i<(int)_Px_value.size(); i++){
                if (_Px_value[i] <= 0 || _Px_value[i] >= 1)
                    throw runtime_error("in PairHMM::set_Px, Px_value should strictly in (0,1).");
            }
            if (fabs(sum(_Px_value) - 1) > TOL)
                throw runtime_error("in PairHMM::set_Px, Px_value is not a probability.");

            Px_value = _Px_value;
            log_Px_value.clear();
            for (int i=0; i<(int)Px_value.size(); i++)
                log_Px_value.push_back(log(Px_value[i]));

            is_set_Px = true;
        }
        inline void set_Py(vector <double> & _Py_value)
        {
            if (_Py_value.size() != 4 && _Py_value.size() != 5)
                throw runtime_error("in PairHMM::set_Py, size of Py_value should be 4 or 5.");

            for (int i=0; i<(int)_Py_value.size(); i++){
                if (_Py_value[i] <= 0 || _Py_value[i] >= 1)
                    throw runtime_error("in PairHMM::set_Py, Py_value should strictly in (0,1).");
            }
            if (fabs(sum(_Py_value) - 1) > TOL )
                throw runtime_error("in PairHMM::set_Py, Py_value is not a probability.");

            Py_value = _Py_value;
            log_Py_value.clear();
            for (int i=0; i<(int)Py_value.size(); i++)
                log_Py_value.push_back(log(Py_value[i]));

            is_set_Py = true;
        }
        inline void set_Pxy(Matrix <double> & _Pxy_value)
        {
            if ((_Pxy_value.nrow != 4 && _Pxy_value.nrow != 5) || (_Pxy_value.ncol != 4 && _Pxy_value.ncol != 5))
                throw runtime_error("in PairHMM::set_Pxy, size of Pxy_value should be 4x4 or 5x5.");

            for (int i=0; i < _Pxy_value.nrow; i++){
                for (int j=0; j < _Pxy_value.ncol; j++){
                    if (_Pxy_value[i][j] <= 0 || _Pxy_value[i][j] >= 1)
                        throw runtime_error("in PairHMM::set_Pxy, Pxy_value should be strictly in (0,1).");
                }
            }
            if (fabs(sum(_Pxy_value) - 1) > TOL)
                throw runtime_error("in PairHMM::set_Pxy, Pxy_value is not a probability.");

            Pxy_value = _Pxy_value;
            log_Pxy_value.setDim(Pxy_value.nrow, Pxy_value.ncol);
            for (int i=0; i<Pxy_value.nrow; i++)
                for (int j=0; j<Pxy_value.ncol; j++)
                    log_Pxy_value[i][j] = log(Pxy_value[i][j]);

            is_set_Pxy = true;
        }

        // get probablity Px, Py, Pxy, index is A=0, C=1, G=2, T=3, N=4
        inline double Px(const char &x, int is_log=false){
            if (Px_value.size() != 4 && Px_value.size() != 5)
                throw runtime_error("in PairHMM::Px, size of Px_value should be 4 or 5.");
            int i;
            switch(x){
                case 'A':
                    i = 0;
                    break;
                case 'C':
                    i = 1;
                    break;
                case 'G':
                    i = 2;
                    break;
                case 'T':
                    i = 3;
                    break;
                case 'N':
                    i = 4;
                    break;
                default:
                    throw runtime_error("PairHMM::Px, input should be A, C, G, T, or N");
            }
            if (is_log)
                return log_Px_value[i];
            else
                return Px_value[i];
        }
        inline double Py(const char &y, int is_log=false){
            if (Py_value.size() != 4 && Py_value.size() != 5)
                throw runtime_error("in PairHMM::Py, size of Py_value should be 4 or 5.");
            int j;
            switch(y){
                case 'A':
                    j = 0;
                    break;
                case 'C':
                    j = 1;
                    break;
                case 'G':
                    j = 2;
                    break;
                case 'T':
                    j = 3;
                    break;
                case 'N':
                    j = 4;
                    break;
                default:
                    throw runtime_error("PairHMM::Py, input should be A, C, G, T, or N");
            }
            if (is_log)
                return log_Py_value[j];
            else
                return Py_value[j];
        }
        inline double Pxy(const char &x, const char &y, int is_log=false){
            if ((Pxy_value.nrow != 4 && Pxy_value.nrow != 5) || (Pxy_value.ncol != 4 && Pxy_value.ncol != 5))
                throw runtime_error("in PairHMM::Pxy, size of Pxy_value should be 4x4 or 5x5.");
            int i, j;
            switch(x){
                case 'A':
                    i = 0;
                    break;
                case 'C':
                    i = 1;
                    break;
                case 'G':
                    i = 2;
                    break;
                case 'T':
                    i = 3;
                    break;
                case 'N':
                    i = 4;
                    break;
                default:
                    throw runtime_error("PairHMM::Pxy, input should be A, C, G, T, or N");
            }
            switch(y){
                case 'A':
                    j = 0;
                    break;
                case 'C':
                    j = 1;
                    break;
                case 'G':
                    j = 2;
                    break;
                case 'T':
                    j = 3;
                    break;
                case 'N':
                    j = 4;
                    break;
                default:
                    throw runtime_error("PairHMM::Pxy, input should be A, C, G, T, or N");
            }
            if (is_log)
                return log_Pxy_value[i][j];
            else
                return Pxy_value[i][j];
        }

        inline double transProb(char S_prev, char S, bool is_log = true)
        {
            if (transProb_value.nrow != 3 || transProb_value.ncol != 3)
                throw runtime_error("transProb_value should be 3x3.");
            int i,j;
            switch(S_prev){
                case 'M':
                    i = 0;
                    break;
                case 'I':
                    i = 1;
                    break;
                case 'D':
                    i = 2;
                    break;
                default:
                    throw runtime_error("state should be M, I or D");
            }
            switch(S){
                case 'M':
                    j = 0;
                    break;
                case 'I':
                    j = 1;
                    break;
                case 'D':
                    j = 2;
                    break;
                default:
                    throw runtime_error("state should be M, I or D");
            }
            if (is_log)
                return log_transProb_value[i][j];
            else
                return transProb_value[i][j];
        }

        // get cigar
        inline vector < pair<char, int> > & get_cigar(){return cigar;}

        // print cigar
        inline void print_cigar(ostream & out = cout)
        {
            for (int i=0; i<(int)cigar.size(); i++)
                out << cigar[i].second << cigar[i].first;
        }

        // get scoreMat, fwdMat and bwdMat
        inline ScoreMatrix & get_scoreMat(){return scoreMat;}
        inline ScoreMatrix & get_fwdMat(){return fwdMat;}
        inline ScoreMatrix & get_bwdMat(){return bwdMat;}

        inline void print_Mat(ScoreMatrix & Mat, ostream & out = cout, bool is_log_V=true, bool is_V=false, bool is_path=false){
            if (is_log_V){
                for (int i=0; i<Mat.nrow; i++){
                    for (int j=0; j<Mat.ncol; j++){
                        out << '(' << Mat[i][j].log_Vm << ',' << Mat[i][j].log_Vx << ',' << Mat[i][j].log_Vy << ")\t";
                    }
                    out << endl;
                }
            }
        }

        // -------------- algorithms ----------- //
        double cal_likelihood_from_cigar(bool is_log = true);
        void viterbi();
        double forward();
        double backward();

    protected:
        // input sequences
        string seqX;
        string seqY;

        // parameters: eps is "prob from indel to indel", dlt is "prob from match to indel"
        double eps_x, dlt_x, log_eps_x, log_dlt_x, eps_y, dlt_y, log_eps_y, log_dlt_y;

        // derived parameters: a_mm = 1 - dlt_x - dlt_y, a_mi_x = 1 - eps_x, a_mi_y = 1 - eps_y
        double a_mm, a_im_x, a_im_y, log_a_mm, log_a_im_x, log_a_im_y;

        //  marginal probability of X
        vector<double> Px_value;
        vector<double> log_Px_value;

        //  marginal probability of Y
        vector<double> Py_value;
        vector<double> log_Py_value;

        //  joint probability of XY
        Matrix<double> Pxy_value;
        Matrix<double> log_Pxy_value;

        // transition probability
        Matrix<double> transProb_value;
        Matrix<double> log_transProb_value;

        // alignment cigar
        vector < pair<char, int> > cigar;

        // score matrix for viterbi, work on log space
        ScoreMatrix scoreMat;

        // score matrix for forward algorithm, work on log space
        ScoreMatrix fwdMat;

        // score matrix for backward algorithm, work on log space
        ScoreMatrix bwdMat;

        // indicator, have to set seqX, seqY, and parameters before running algorithms
        bool is_setSeq;
        bool is_setPar;
        bool is_set_Px;
        bool is_set_Py;
        bool is_set_Pxy;

    protected:
        // validate sequence
        void validateSeq(string & seq)
        {
            for (int i=0; i<(int)seq.size(); i++){
                if (seq[i] != 'A' && seq[i] != 'C' && seq[i] != 'G' && seq[i] != 'T' && seq[i] != 'N'){
                    throw runtime_error("in PairHMM::validateSeq, sequence contains letter other than A, C, G, T or N.");
                }
            }
            if (seq.size() < 1)
                throw runtime_error("empty sequences is not allowed.");
        }

};



#endif // PAIRHMM

