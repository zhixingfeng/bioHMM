#ifndef TEST_PAIRHMM_H
#define TEST_PAIRHMM_H

#include "../pairhmm.h"
#include <UnitTest++/UnitTest++.h>
SUITE(TestPairHMM)
{
    class TestPairHMM
    {
        public:
            PairHMM hmm;
            ScoreMatrix scoremat;
        public:
            void set_std_par(){
                hmm.setPar(0.02, 0.01, 0.02, 0.01);
                vector<double> Px(4, 0.25);
                hmm.set_Px(Px);
                vector<double> Py(4, 0.25);
                hmm.set_Py(Py);
                Matrix <double> Pxy(4,4);
                vector<double> cur_row(4,0.01);
                for (int i=0; i<Pxy.nrow; i++) Pxy[i] = cur_row;
                for (int i=0; i<Pxy.nrow; i++) Pxy[i][i] = 0.22;
                hmm.set_Pxy(Pxy);
            }
    };

    TEST_FIXTURE(TestPairHMM, PairHMM_setSeq)
    {
        // set normal seq
        hmm.setSeq("CGATGC","CGTATTTTTTNNT");
        CHECK_EQUAL(hmm.getSeqX(), "CGATGC");
        CHECK_EQUAL(hmm.getSeqY(), "CGTATTTTTTNNT");

        // set abnormal seq
        CHECK_THROW(hmm.setSeq("CGATGCXX","CGTATTTTTTNNT"), exception);
        CHECK_THROW(hmm.setSeq("CGATGCN","CGTATTTTTTNNTXXX"), exception);
        CHECK_THROW(hmm.setSeq("CGATGC12342  ","CGTATTTTTTNNTXXX"), exception);
        CHECK_THROW(hmm.setSeq("","CGTATTTTTTNNTXXX"), exception);
        CHECK_THROW(hmm.setSeq("CGATGC",""), exception);
    }

    TEST_FIXTURE(TestPairHMM, PairHMM_setPar){
        // set normal par
        hmm.setPar(0.1,0.8, 0.4,0.1);
        CHECK_EQUAL(hmm.getPar_eps_x(), 0.1);
        CHECK_EQUAL(hmm.getPar_dlt_x(), 0.8);
        CHECK_EQUAL(hmm.getPar_log_eps_x(), log(0.1));
        CHECK_EQUAL(hmm.getPar_log_dlt_x(), log(0.8));

        CHECK_EQUAL(hmm.getPar_eps_y(), 0.4);
        CHECK_EQUAL(hmm.getPar_dlt_y(), 0.1);
        CHECK_EQUAL(hmm.getPar_log_eps_y(), log(0.4));
        CHECK_EQUAL(hmm.getPar_log_dlt_y(), log(0.1));

        CHECK_EQUAL(hmm.getPar_a_mm(), 1 - 0.8 - 0.1);
        CHECK_EQUAL(hmm.getPar_a_im_x(), 1 - 0.1);
        CHECK_EQUAL(hmm.getPar_a_im_y(), 1 - 0.4);
        CHECK_EQUAL(hmm.getPar_log_a_mm(), log(1 - 0.8 - 0.1));
        CHECK_EQUAL(hmm.getPar_log_a_im_x(), log(1 - 0.1));
        CHECK_EQUAL(hmm.getPar_log_a_im_y(), log(1 - 0.4));


        // set abnormal par
        CHECK_THROW(hmm.setPar(0, 0.8, 0.4, 0.1), exception);
        CHECK_THROW(hmm.setPar(0.1, 1, 0.4, 0), exception);
        CHECK_THROW(hmm.setPar(0, 1, 0.4, 0), exception);
        CHECK_THROW(hmm.setPar(-10390, 1879371937, 0.4, 0.6), exception);
        CHECK_THROW(hmm.setPar(0.1, 0.8, 0.4, 0.21), exception);
    }

    TEST_FIXTURE(TestPairHMM, PairHMM_set_Px)
    {
        // set_Px, normal, A,C,G,T,N
        vector<double> P_1(5,-1);
        P_1[0]=0.1, P_1[1]=0.2, P_1[2]=0.4, P_1[3]=0.25, P_1[4]=0.05;
        hmm.set_Px(P_1);
        CHECK_EQUAL(hmm.Px('A'), 0.1);CHECK_EQUAL(hmm.Px('C'), 0.2);CHECK_EQUAL(hmm.Px('G'), 0.4);CHECK_EQUAL(hmm.Px('T'), 0.25); CHECK_EQUAL(hmm.Px('N'), 0.05);
        CHECK_CLOSE(exp(hmm.Px('A',true)), 0.1,10e-6);CHECK_CLOSE(exp(hmm.Px('C',true)), 0.2,10e-6);CHECK_CLOSE(exp(hmm.Px('G',true)), 0.4,10e-6);CHECK_CLOSE(exp(hmm.Px('T',true)), 0.25,10e-6); CHECK_CLOSE(exp(hmm.Px('N',true)), 0.05,10e-6);

        // set_Px, normal, A,C,G,T
        vector<double> P_2(4,-1);
        P_2[0]=0.1, P_2[1]=0.2, P_2[2]=0.4, P_2[3]=0.3; hmm.set_Px(P_2);
        CHECK_EQUAL(hmm.Px('A'), 0.1); CHECK_EQUAL(hmm.Px('C'), 0.2); CHECK_EQUAL(hmm.Px('G'), 0.4); CHECK_EQUAL(hmm.Px('T'), 0.3);
        CHECK_CLOSE(exp(hmm.Px('A',true)), 0.1,10e-6);CHECK_CLOSE(exp(hmm.Px('C',true)), 0.2,10e-6);CHECK_CLOSE(exp(hmm.Px('G',true)), 0.4,10e-6);CHECK_CLOSE(exp(hmm.Px('T',true)), 0.3,10e-6);

        // set_Px, abnormal, incorrect size
        vector<double> P_3(2,-1);
        P_3[0]=0.5, P_3[1]=0.5 ;
        CHECK_THROW(hmm.set_Px(P_3), exception);

        // set_Px, abnormal, not a prob
        vector<double> P_4(4,-1);
        P_4[0]=0.1, P_4[1]=0.1; P_4[2]=0.1; P_4[3]=0.1;
        CHECK_THROW(hmm.set_Px(P_4), exception);
    }

    TEST_FIXTURE(TestPairHMM, PairHMM_set_Py)
    {
        // set_Py, normal, A,C,G,T,N
        vector<double> P_1(5,-1);
        P_1[0]=0.1, P_1[1]=0.2, P_1[2]=0.4, P_1[3]=0.25, P_1[4]=0.05;
        hmm.set_Py(P_1);
        CHECK_EQUAL(hmm.Py('A'), 0.1);CHECK_EQUAL(hmm.Py('C'), 0.2);CHECK_EQUAL(hmm.Py('G'), 0.4);CHECK_EQUAL(hmm.Py('T'), 0.25); CHECK_EQUAL(hmm.Py('N'), 0.05);
        CHECK_CLOSE(exp(hmm.Py('A',true)), 0.1,10e-6);CHECK_CLOSE(exp(hmm.Py('C',true)), 0.2,10e-6);CHECK_CLOSE(exp(hmm.Py('G',true)), 0.4,10e-6);CHECK_CLOSE(exp(hmm.Py('T',true)), 0.25,10e-6); CHECK_CLOSE(exp(hmm.Py('N',true)), 0.05,10e-6);
        // set_Py, normal, A,C,G,T
        vector<double> P_2(4,-1);
        P_2[0]=0.1, P_2[1]=0.2, P_2[2]=0.4, P_2[3]=0.3; hmm.set_Py(P_2);
        CHECK_EQUAL(hmm.Py('A'), 0.1); CHECK_EQUAL(hmm.Py('C'), 0.2); CHECK_EQUAL(hmm.Py('G'), 0.4); CHECK_EQUAL(hmm.Py('T'), 0.3);
        CHECK_CLOSE(exp(hmm.Py('A',true)), 0.1,10e-6);CHECK_CLOSE(exp(hmm.Py('C',true)), 0.2,10e-6);CHECK_CLOSE(exp(hmm.Py('G',true)), 0.4,10e-6);CHECK_CLOSE(exp(hmm.Py('T',true)), 0.3,10e-6);
        // set_Py, abnormal, incorrect size
        vector<double> P_3(2,-1);
        P_3[0]=0.5, P_3[1]=0.5 ;
        CHECK_THROW(hmm.set_Py(P_3), exception);

        // set_Py, abnormal, not a prob
        vector<double> P_4(4,-1);
        P_4[0]=0.1, P_4[1]=0.1; P_4[2]=0.1; P_4[3]=0.1;
        CHECK_THROW(hmm.set_Py(P_4), exception);
    }

    TEST_FIXTURE(TestPairHMM, PairHMM_set_Pxy_normal_5_letters)
    {
        // set_Pxy, normal, A,C,G,T,N
        Matrix<double> P_1(5,5);
        for (int i=0; i<P_1.nrow; i++)
            for (int j=0; j<P_1.ncol; j++){
                P_1[i][j] = (i*P_1.ncol + j + 1)/325.0;
            }
        hmm.set_Pxy(P_1);


        Matrix<double> _P_1(5,5);
        _P_1[0][0]=hmm.Pxy('A','A'); _P_1[0][1]=hmm.Pxy('A','C'); _P_1[0][2]=hmm.Pxy('A','G'); _P_1[0][3]=hmm.Pxy('A','T'); _P_1[0][4]=hmm.Pxy('A','N');
        _P_1[1][0]=hmm.Pxy('C','A'); _P_1[1][1]=hmm.Pxy('C','C'); _P_1[1][2]=hmm.Pxy('C','G'); _P_1[1][3]=hmm.Pxy('C','T'); _P_1[1][4]=hmm.Pxy('C','N');
        _P_1[2][0]=hmm.Pxy('G','A'); _P_1[2][1]=hmm.Pxy('G','C'); _P_1[2][2]=hmm.Pxy('G','G'); _P_1[2][3]=hmm.Pxy('G','T'); _P_1[2][4]=hmm.Pxy('G','N');
        _P_1[3][0]=hmm.Pxy('T','A'); _P_1[3][1]=hmm.Pxy('T','C'); _P_1[3][2]=hmm.Pxy('T','G'); _P_1[3][3]=hmm.Pxy('T','T'); _P_1[3][4]=hmm.Pxy('T','N');
        _P_1[4][0]=hmm.Pxy('N','A'); _P_1[4][1]=hmm.Pxy('N','C'); _P_1[4][2]=hmm.Pxy('N','G'); _P_1[4][3]=hmm.Pxy('N','T'); _P_1[4][4]=hmm.Pxy('N','N');
        CHECK_ARRAY_EQUAL(P_1[0], _P_1[0], 5);
        CHECK_ARRAY_EQUAL(P_1[1], _P_1[1], 5);
        CHECK_ARRAY_EQUAL(P_1[2], _P_1[2], 5);
        CHECK_ARRAY_EQUAL(P_1[3], _P_1[3], 5);
        CHECK_ARRAY_EQUAL(P_1[4], _P_1[4], 5);

        Matrix<double> _P_1_log(5,5);
        _P_1_log[0][0]=hmm.Pxy('A','A',true); _P_1_log[0][1]=hmm.Pxy('A','C',true); _P_1_log[0][2]=hmm.Pxy('A','G',true); _P_1_log[0][3]=hmm.Pxy('A','T',true); _P_1_log[0][4]=hmm.Pxy('A','N',true);
        _P_1_log[1][0]=hmm.Pxy('C','A',true); _P_1_log[1][1]=hmm.Pxy('C','C',true); _P_1_log[1][2]=hmm.Pxy('C','G',true); _P_1_log[1][3]=hmm.Pxy('C','T',true); _P_1_log[1][4]=hmm.Pxy('C','N',true);
        _P_1_log[2][0]=hmm.Pxy('G','A',true); _P_1_log[2][1]=hmm.Pxy('G','C',true); _P_1_log[2][2]=hmm.Pxy('G','G',true); _P_1_log[2][3]=hmm.Pxy('G','T',true); _P_1_log[2][4]=hmm.Pxy('G','N',true);
        _P_1_log[3][0]=hmm.Pxy('T','A',true); _P_1_log[3][1]=hmm.Pxy('T','C',true); _P_1_log[3][2]=hmm.Pxy('T','G',true); _P_1_log[3][3]=hmm.Pxy('T','T',true); _P_1_log[3][4]=hmm.Pxy('T','N',true);
        _P_1_log[4][0]=hmm.Pxy('N','A',true); _P_1_log[4][1]=hmm.Pxy('N','C',true); _P_1_log[4][2]=hmm.Pxy('N','G',true); _P_1_log[4][3]=hmm.Pxy('N','T',true); _P_1_log[4][4]=hmm.Pxy('N','N',true);
        for(int i=0; i<5; i++)
            for (int j=0; j<5; j++)
                CHECK_CLOSE(P_1[i][j],exp(_P_1_log[i][j]),10e-6);
    }

    TEST_FIXTURE(TestPairHMM, PairHMM_set_Pxy_normal_4_letters)
    {
        // set_Pxy, normal, A,C,G,T,N
        Matrix<double> P_1(4,4);
        for (int i=0; i<P_1.nrow; i++)
            for (int j=0; j<P_1.ncol; j++){
                P_1[i][j] = (i*P_1.ncol + j + 1)/136.0;
            }
        hmm.set_Pxy(P_1);

        Matrix<double> _P_1(4,4);
        _P_1[0][0]=hmm.Pxy('A','A'); _P_1[0][1]=hmm.Pxy('A','C'); _P_1[0][2]=hmm.Pxy('A','G'); _P_1[0][3]=hmm.Pxy('A','T');
        _P_1[1][0]=hmm.Pxy('C','A'); _P_1[1][1]=hmm.Pxy('C','C'); _P_1[1][2]=hmm.Pxy('C','G'); _P_1[1][3]=hmm.Pxy('C','T');
        _P_1[2][0]=hmm.Pxy('G','A'); _P_1[2][1]=hmm.Pxy('G','C'); _P_1[2][2]=hmm.Pxy('G','G'); _P_1[2][3]=hmm.Pxy('G','T');
        _P_1[3][0]=hmm.Pxy('T','A'); _P_1[3][1]=hmm.Pxy('T','C'); _P_1[3][2]=hmm.Pxy('T','G'); _P_1[3][3]=hmm.Pxy('T','T');
        CHECK_ARRAY_EQUAL(P_1[0], _P_1[0], 4);
        CHECK_ARRAY_EQUAL(P_1[1], _P_1[1], 4);
        CHECK_ARRAY_EQUAL(P_1[2], _P_1[2], 4);
        CHECK_ARRAY_EQUAL(P_1[3], _P_1[3], 4);
    }

    TEST_FIXTURE(TestPairHMM, PairHMM_set_Pxy_not_prob_4_letters)
    {
        // set_Pxy, normal, A,C,G,T,N
        Matrix<double> P_1(4,4);
        for (int i=0; i<P_1.nrow; i++)
            for (int j=0; j<P_1.ncol; j++){
                P_1[i][j] = (i*P_1.ncol + j + 1)/130.0;
            }
        CHECK_THROW(hmm.set_Pxy(P_1), exception);

    }

    TEST_FIXTURE(TestPairHMM, PairHMM_set_Pxy_not_prob_5_letters)
    {
        // set_Pxy, normal, A,C,G,T,N
        Matrix<double> P_1(5,5);
        for (int i=0; i<P_1.nrow; i++)
            for (int j=0; j<P_1.ncol; j++){
                P_1[i][j] = (i*P_1.ncol + j + 1)/230.0;
            }
        CHECK_THROW(hmm.set_Pxy(P_1), exception);
    }

    TEST_FIXTURE(TestPairHMM, PairHMM_set_Pxy_incorrect_size)
    {
        // set_Pxy, normal, A,C,G,T,N
        Matrix<double> P_1(10,10);
        for (int i=0; i<P_1.nrow; i++)
            for (int j=0; j<P_1.ncol; j++){
                P_1[i][j] = (i*P_1.ncol + j + 1)/5050.0;
            }
        CHECK_THROW(hmm.set_Pxy(P_1), exception);
    }

    TEST_FIXTURE(TestPairHMM, PairHMM_viterbi_basic)
    {
        CHECK_THROW(hmm.viterbi(), exception);
        // set seqX and seqY
        hmm.setSeq("AACCTGAGAG", "ACGTGGAG");
        CHECK_THROW(hmm.viterbi(), exception);
        // eps_x = 0.02, dlt_x = 0.01, eps_y=0.02, eps_y=0.01
        hmm.setPar(0.02, 0.01, 0.02, 0.01);
        CHECK_THROW(hmm.viterbi(), exception);
        // Px = [0.25, 0.25, 0.25, 0.25]
        vector<double> Px(4, 0.25);
        hmm.set_Px(Px);
        CHECK_THROW(hmm.viterbi(), exception);
        // Py = [0.25, 0.25, 0.25, 0.25]
        vector<double> Py(4, 0.25);
        hmm.set_Py(Py);
        CHECK_THROW(hmm.viterbi(), exception);

        // Pxy = [0.22, 0.01, 0.01, 0.01
        //        0.01, 0.22, 0.01, 0.01
        //        0.01, 0.01, 0.22, 0.01
        //        0.01, 0.01, 0.01, 0.22]
        Matrix <double> Pxy(4,4);
        vector<double> cur_row(4,0.01);
        for (int i=0; i<Pxy.nrow; i++) Pxy[i] = cur_row;
        for (int i=0; i<Pxy.nrow; i++) Pxy[i][i] = 0.22;
        hmm.set_Pxy(Pxy);

        // run viterbi
        try{
            hmm.viterbi();
        }
        catch (exception &e){
            cout << e.what() << endl;
        }

        // get scoreMat
        Matrix<ScoreCell> outMat = hmm.get_scoreMat();
        double log_max_L = outMat[outMat.nrow - 1][outMat.ncol - 1].log_Vm;
        if (log_max_L < outMat[outMat.nrow - 1][outMat.ncol - 1].log_Vx)
            log_max_L = outMat[outMat.nrow - 1][outMat.ncol - 1].log_Vx;
        if (log_max_L < outMat[outMat.nrow - 1][outMat.ncol - 1].log_Vy)
            log_max_L = outMat[outMat.nrow - 1][outMat.ncol - 1].log_Vy;
        CHECK_CLOSE(hmm.cal_likelihood_from_cigar(), log_max_L, TOL);
        CHECK_CLOSE(hmm.cal_likelihood_from_cigar(false), exp(log_max_L), TOL);

        //hmm.print_Mat(hmm.get_scoreMat());
        /*
        for (int i=0; i<outMat.nrow; i++){
            for (int j=0; j<outMat.ncol; j++){
                cout << '(' << outMat[i][j].log_Vm << ',' << outMat[i][j].log_Vx << ',' << outMat[i][j].log_Vy << ")\t";
            }
            cout << endl;
        }

        for (int i=0; i<outMat.nrow; i++){
            for (int j=0; j<outMat.ncol; j++){
                cout << '(' << outMat[i][j].path_M << ',' << outMat[i][j].path_X << ',' << outMat[i][j].path_Y << ")\t";
            }
            cout << endl;
        }
        // get cigar
        cout << "cigar" << endl;
        vector<pair<char, int> > cigar = hmm.get_cigar();
        for (int i=0; i<(int)cigar.size(); i++)
            cout << '(' << cigar[i].first << ',' << cigar[i].second << "),";
        cout << endl;
        hmm.print_cigar();
        cout << endl;

        // get max log-likelihood
        cout << "max log-likelihood is " << hmm.cal_likelihood_from_cigar() << endl;
        cout << "max likelihood is " << hmm.cal_likelihood_from_cigar(false) << endl;
        */
    }

    TEST_FIXTURE(TestPairHMM, PairHMM_viterbi_random_input)
    {
        set_std_par();
        hmm.setSeq("AACCTGAGAGGTGGGGCGATGCGATCGATC", "ACGTGGAGGTGGGGTCGATGCTAGCGTAGCGACGCTACT");

        // run viterbi
        try{
            hmm.viterbi();
        }
        catch (exception &e){
            cout << e.what() << endl;
        }


    }
    TEST_FIXTURE(TestPairHMM, PairHMM_transProb)
    {
        // eps_x = 0.02, dlt_x = 0.01, eps_y=0.02, eps_y=0.01
        hmm.setPar(0.02, 0.01, 0.02, 0.01);
        CHECK_EQUAL(hmm.transProb('M','M',false), 1 - 0.01 - 0.01);
        CHECK_EQUAL(hmm.transProb('M','I',false), 0.01);
        CHECK_EQUAL(hmm.transProb('M','D',false), 0.01);
        CHECK_EQUAL(hmm.transProb('I','M',false), 1 - 0.02);
        CHECK_EQUAL(hmm.transProb('I','I',false), 0.02);
        CHECK_EQUAL(hmm.transProb('I','D',false), 0);
        CHECK_EQUAL(hmm.transProb('D','M',false), 1 - 0.02);
        CHECK_EQUAL(hmm.transProb('D','D',false), 0.02);
        CHECK_EQUAL(hmm.transProb('D','I',false), 0);

        CHECK_EQUAL(hmm.transProb('M','M',true), log(1 - 0.01 - 0.01));
        CHECK_EQUAL(hmm.transProb('M','I',true), log(0.01));
        CHECK_EQUAL(hmm.transProb('M','D',true), log(0.01));
        CHECK_EQUAL(hmm.transProb('I','M',true), log(1 - 0.02));
        CHECK_EQUAL(hmm.transProb('I','I',true), log(0.02));
        CHECK_EQUAL(hmm.transProb('I','D',true), LOG_MIN);
        CHECK_EQUAL(hmm.transProb('D','M',true), log(1 - 0.02));
        CHECK_EQUAL(hmm.transProb('D','D',true), log(0.02));
        CHECK_EQUAL(hmm.transProb('D','I',true), LOG_MIN);
    }

    TEST_FIXTURE(TestPairHMM, PairHMM_forward)
    {
        set_std_par();
        //hmm.setSeq("AACCTGAGAGGTGGGGCGATGCGATCGATC", "ACGTGGAGGTGGGGTCGATGCTAGCGTAGCGACGCTACT");
        hmm.setSeq("ACCTGAGAG", "ACGTGGAG");
        // run forward algorithm
        //hmm.forward(hmm.getSeqX().size()+1, hmm.getSeqY().size()+1);
        double log_L_fwd = hmm.forward();
        //hmm.print_Mat(hmm.get_fwdMat());
        //cout << log_L_fwd << endl << endl;

        /*hmm.setSeq("ACATGAGAG", "ACGTGGAG");
        cout << hmm.forward() << endl;
        hmm.setSeq("ACCTGAGAG", "ACGTGGAG");
        cout << hmm.forward() << endl;
        hmm.setSeq("ACGTGAGAG", "ACGTGGAG");
        cout << hmm.forward() << endl;
        hmm.setSeq("ACTTGAGAG", "ACGTGGAG");
        cout << hmm.forward() << endl;*/

        hmm.setSeq("GATC", "GATTC");
        double log_A = hmm.forward();
        cout << log_A << endl;
        hmm.print_Mat(hmm.get_fwdMat());
        hmm.setSeq("GCTC", "GATTC");
        double log_C = hmm.forward();
        cout << log_C << endl;
        hmm.setSeq("GGTC", "GATTC");
        double log_G = hmm.forward();
        cout << log_G << endl;
        hmm.setSeq("GTTC", "GATTC");
        double log_T = hmm.forward();
        cout << log_T << endl;
        hmm.print_Mat(hmm.get_fwdMat());
        double norm = exp(log_A) + exp(log_C) + exp(log_G) + exp(log_T);
        cout << "prob:" << endl;
        cout << exp(log_A)/norm << ',' << exp(log_C)/norm << ',' << exp(log_G)/norm << ',' << exp(log_T)/norm << endl;

    }
}



#endif // TEST_PAIRHMM

