#include "stl.h"
#include <tclap/CmdLine.h>
#include "pairhmm.h"
#include "./test/test_matrix.cpp"
int main(int argc, char *argv[])
{
    Matrix <int> mat(2, 3);
    return UnitTest::RunAllTests();
}

