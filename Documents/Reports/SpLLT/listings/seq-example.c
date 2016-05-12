for (i = 1; i < N; i++) {
  x[i] = f(x[i]);
  y[i] = g(x[i], y[i-1]);
}
