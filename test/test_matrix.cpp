#include "../matrix.h"
#include <UnitTest++/UnitTest++.h>
SUITE(TestMatrix)
{
    class TestMatrix
    {
        public:
            Matrix <double> mat;
            void setValue(vector<double> &x, int nrow, int ncol)
            {
                if (nrow * ncol != (int)x.size()){
                    throw std::runtime_error("nrow * ncol > x.size()");
                }
                mat.setDim(nrow, ncol);
                for (int i=0; i<nrow; i++){
                    for (int j=0; j<ncol; j++){
                        mat.value[i][j] = x[i*ncol + j];
                    }
                }
            }

    };

    TEST_FIXTURE(TestMatrix, Matrix_setValue)
    {
        vector<double> x;
        for (int i=1; i<=12; i++) x.push_back(i);
        setValue(x, 3, 4);
        vector<double> y(12,-1);
        for (int i=0; i<3; i++){
            for (int j=0; j<4; j++){
                y[i*4 + j] = mat[i][j];
            }
        }
        CHECK_ARRAY_EQUAL(x,y,12);
    }

}
