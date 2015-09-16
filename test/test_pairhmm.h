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
    }

    TEST_FIXTURE(TestPairHMM, PairHMM_setPar){
        // set normal par
        hmm.setPar(0.1,0.8);
        CHECK_EQUAL(hmm.getPar_eps(), 0.1);
        CHECK_EQUAL(hmm.getPar_dlt(), 0.8);
        CHECK_EQUAL(hmm.getPar_log_eps(), log(0.1));
        CHECK_EQUAL(hmm.getPar_log_dlt(), log(0.8));

        // set abnormal par
        CHECK_THROW(hmm.setPar(0,0.8), exception);
        CHECK_THROW(hmm.setPar(0.1,1), exception);
        CHECK_THROW(hmm.setPar(0,1), exception);
        CHECK_THROW(hmm.setPar(-10390,1879371937), exception);
    }

    TEST_FIXTURE(TestPairHMM, PairHMM_set_Px)
    {
        // set_Px, normal, A,C,G,T,N
        vector<double> P_1(5,-1);
        P_1[0]=0.1, P_1[1]=0.2, P_1[2]=0.4, P_1[3]=0.25, P_1[4]=0.05;
        hmm.set_Px(P_1);
        CHECK_EQUAL(hmm.Px('A'), 0.1);CHECK_EQUAL(hmm.Px('C'), 0.2);CHECK_EQUAL(hmm.Px('G'), 0.4);CHECK_EQUAL(hmm.Px('T'), 0.25); CHECK_EQUAL(hmm.Px('N'), 0.05);

        // set_Px, normal, A,C,G,T
        vector<double> P_2(4,-1);
        P_2[0]=0.1, P_2[1]=0.2, P_2[2]=0.4, P_2[3]=0.3; hmm.set_Px(P_2);
        CHECK_EQUAL(hmm.Px('A'), 0.1); CHECK_EQUAL(hmm.Px('C'), 0.2); CHECK_EQUAL(hmm.Px('G'), 0.4); CHECK_EQUAL(hmm.Px('T'), 0.3);

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

        // set_Py, normal, A,C,G,T
        vector<double> P_2(4,-1);
        P_2[0]=0.1, P_2[1]=0.2, P_2[2]=0.4, P_2[3]=0.3; hmm.set_Py(P_2);
        CHECK_EQUAL(hmm.Py('A'), 0.1); CHECK_EQUAL(hmm.Py('C'), 0.2); CHECK_EQUAL(hmm.Py('G'), 0.4); CHECK_EQUAL(hmm.Py('T'), 0.3);

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
        try{
        hmm.set_Pxy(P_1);
        }
        catch(exception &e){
            cout << e.what() << endl;
        }

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
    }

    TEST_FIXTURE(TestPairHMM, PairHMM_set_Pxy_normal_4_letters)
    {
        // set_Pxy, normal, A,C,G,T,N
        Matrix<double> P_1(4,4);
        for (int i=0; i<P_1.nrow; i++)
            for (int j=0; j<P_1.ncol; j++){
                P_1[i][j] = (i*P_1.ncol + j + 1)/136.0;
            }
        try{
            hmm.set_Pxy(P_1);
        }
        catch(exception &e){
            cout << e.what() << endl;
        }

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


}



#endif // TEST_PAIRHMM

