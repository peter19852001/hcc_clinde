#ifndef SIMPLE_TAB_SEPARATED_VALUES_READER
#define SIMPLE_TAB_SEPARATED_VALUES_READER

#include <stdio.h>

double* read_TSV_file(FILE* fp, int* out_rows, int* out_cols);
double* read_TSV(char* filename, int* out_rows, int* out_cols);

#endif
