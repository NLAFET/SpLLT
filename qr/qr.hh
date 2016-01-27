#ifndef QR_H
#define QR_H

#include <cmath>

//void qr(int m, int n, double const *a, int lda);

void qr_mgs(int m, int n, double *a, int lda, double *r, int ldr);

#endif // QR_H
