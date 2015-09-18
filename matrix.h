#ifndef MATRIX_H
#define MATRIX_H

#include "stl.h"

template<typename T> struct Matrix
{
    // init
    Matrix(int _nrow, int _ncol){ setDim(_nrow, _ncol); }
    Matrix(){}

    // set dimension
    void setDim(int _nrow, int _ncol){
        for (int i=0; i < (int)value.size(); i++)
            value[i].clear();
        value.clear();
        nrow = _nrow;
        ncol = _ncol;
        for (int i=0; i<nrow; i++){
            vector<T> cur_value;
            for (int j=0; j<ncol; j++){
                T cur_cell;
                cur_value.push_back(cur_cell);
            }
            value.push_back(cur_value);
        }
    }

    // define operator
    vector<T> & operator[](const int i)
    {
        return value[i];
    }



    // cell value
    vector<vector<T> > value;

    // dimension
    int nrow;
    int ncol;
};

#endif // MATRIX

