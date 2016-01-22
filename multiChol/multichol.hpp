#pragma once

void multichol_interleave(int num_problem, int m, int n, double* a, int lda, int* info);
void multichol_standard(int num_problem, int m, int n, double* a, int lda, int* info);
