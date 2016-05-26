for (i = 1; i < N; i++) {
   submit(f, x[i]:RW);
   submit(g, x[i]:R, y[i-1]:R, y[i]:W);
}
