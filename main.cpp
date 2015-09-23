#include "stl.h"
#include <tclap/CmdLine.h>
#include "pairhmm.h"
#include "./test/test_matrix.h"
#include "./test/test_pairhmm.h"
#include "./test/test_logsum.h"

int main(int argc, char *argv[])
{
    /*try{
    PairHMM a;
    //a.Px('G');
    a.log_Py('G');
    }
    catch(exception& e){
        cerr << "Expection: " << e.what() << endl;
    }*/
    return UnitTest::RunAllTests();
}

