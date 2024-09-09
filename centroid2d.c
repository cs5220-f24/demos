#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define N 1000000

void sum1(int n, double* restrict xy, double* restrict result)
{
    double x = 0;
    double y = 0;
    for (int i = 0; i < n; ++i) {
        x += xy[2*i+0];
        y += xy[2*i+1];
    }
    result[0] = x/n;
    result[1] = y/n;
}

void sum2(int n, double* restrict xy, double* restrict result)
{
    double x = 0;
    double y = 0;
    for (int i = 0; i < n; ++i)
        x += xy[2*i+0];
    for (int i = 0; i < n; ++i)
        y += xy[2*i+1];
    result[0] = x/n;
    result[1] = y/n;
}

void sum3(int n, double* restrict xy, double* restrict result)
{
    double x = 0;
    double y = 0;
    for (int i = 0; i < n; ++i)
        x += xy[i];
    for (int i = 0; i < n; ++i)
        y += xy[n+i];
    
    result[0] = x/n;
    result[1] = y/n;
}

void init_rand(int n, double* x)
{
    for (int i = 0; i < N; ++i)
        x[i] = drand48();
}

void timeit(const char* s, int ntrials,
            void (*f)(int, double*, double*), int n, double* xy)
{
    struct timespec ts1, ts2;
    double result[2] = {0, 0};
    timespec_get(&ts1, TIME_UTC);
    for (int trial = 0; trial < ntrials; ++trial)
        f(n, xy, result);
    timespec_get(&ts2, TIME_UTC);
    printf("%s: %g\n", s,
           ((ts2.tv_sec-ts1.tv_sec) + (ts2.tv_nsec-ts1.tv_nsec)*1e-9)/ntrials);
}

int main()
{
    double* xy = (double*) malloc(2*N * sizeof(double));
    init_rand(2*N, xy);
    timeit("V1", 100, sum1, N, xy);
    timeit("V2", 100, sum2, N, xy);
    timeit("V3", 100, sum3, N, xy);
    free(xy);
}
