#ifndef PAIRHMM_H
#define PAIRHMM_H

#include "stl.h"
#include "matrix.h"
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
        // set sequences
        inline void setSeq(string & _seqX, string & _seqY, bool is_val = true)
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
            eps = _eps;
            dlt = _dlt;
            if (eps < 0 || eps > 1 || dlt < 0 || dlt > 1){
                cerr << "Error: eps and dlt should be in [0,1]" << endl;
                exit(1);
            }
        }

        // get parameters
        double getParEps(){return eps;}
        double getParDlt(){return dlt;}

    protected:
        // input sequences
        string seqX;
        string seqY;

        // parameters: eps is "prob from indel to indel", dlt is "prob from match to indel"
        double eps, dlt;

        // alignment cigar
        vector < pair<char, int> > cigar;

        // score/probability matrix
        ScoreMatrix scoreMat;

    protected:
        // validate sequence
        void validateSeq(string & seq)
        {
            for (int i=0; i<(int)seq.size(); i++){
                if (seq[i] != 'A' && seq[i] != 'C' && seq[i] != 'G' && seq[i] != 'T' && seq[i] != 'N'){
                    cerr << "Error: sequence contains letter other than A, C, G, T or N" << endl;
                    exit(1);
                }
            }

        }


};



#endif // PAIRHMM

