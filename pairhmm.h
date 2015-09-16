#ifndef PAIRHMM_H
#define PAIRHMM_H

#include "stl.h"
#include "matrix.h"
#include "math_utils.h"
#define TOL 0.000001

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
        Vm = -1;
        Vx = -1;
        Vy = -1;
        bases_X = 'N';
        bases_Y = 'N';
    }

    // current score of match/mismatch, gap in X, gap in Y
    double Vm, Vx, Vy;

    // matched bases in the current cell, for example <A,A> or <-,G>
    char bases_X;
    char bases_Y;

    // quality score
    QualityScore QV;
};

typedef Matrix<ScoreCell> ScoreMatrix;




class PairHMM
{
    public:
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
        }

        // get sequence
        inline string getSeqX(){return seqX;}
        inline string getSeqY(){return seqY;}

        // set parameters eps and dlt
        inline void setPar(double _eps, double _dlt)
        {
            if (_eps <= 0 || _eps >= 1 || _dlt <= 0 || _dlt >= 1){
                throw runtime_error("in setPar, eps and dlt should be in (0,1)");
            }
            eps = _eps;
            dlt = _dlt;
            log_eps = log(eps);
            log_dlt = log(dlt);
        }

        // get parameters
        inline double getPar_eps(){return eps;}
        inline double getPar_dlt(){return dlt;}
        inline double getPar_log_eps(){return log_eps;}
        inline double getPar_log_dlt(){return log_dlt;}

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

        // get cigar
        inline vector < pair<char, int> > & get_cigar(){return cigar;}
        inline ScoreMatrix & get_scoreMat(){return scoreMat;}
        // -------------- algorithms ----------- //

        void viterbi();
    protected:
        // input sequences
        string seqX;
        string seqY;

        // parameters: eps is "prob from indel to indel", dlt is "prob from match to indel"
        double eps, dlt, log_eps, log_dlt;

        //  marginal probability of X
        vector<double> Px_value;
        vector<double> log_Px_value;

        //  marginal probability of Y
        vector<double> Py_value;
        vector<double> log_Py_value;

        //  joint probability of XY
        Matrix<double> Pxy_value;
        Matrix<double> log_Pxy_value;


        // alignment cigar
        vector < pair<char, int> > cigar;

        // score matrix, work on log space
        ScoreMatrix scoreMat;

    protected:
        // validate sequence
        void validateSeq(string & seq)
        {
            for (int i=0; i<(int)seq.size(); i++){
                if (seq[i] != 'A' && seq[i] != 'C' && seq[i] != 'G' && seq[i] != 'T' && seq[i] != 'N'){
                    throw runtime_error("in PairHMM::validateSeq, sequence contains letter other than A, C, G, T or N");
                }
            }

        }


};



#endif // PAIRHMM

