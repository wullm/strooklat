/*******************************************************************************
 * Copyright (c) 2021 Willem Elbers (whe@willemelbers.com)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

/* Standard headers */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <sys/time.h>

/* GSL interpolation library */
#include <gsl/gsl_spline.h>

/* Strooklat interpolation library */
#include "strooklat_inline.h"

int main() {
    /* Allocate memory for some data points */
    int N = 1000;
    double *x1 = malloc(N * sizeof(double));
    double *x2 = malloc(N * sizeof(double));
    double *y1 = malloc(N * sizeof(double));
    double *y2 = malloc(N * sizeof(double));

    /* Generate the data */
    for (int i = 0; i < N; i++) {
        x1[i] = i * 0.1;
        x2[i] = 100 - i * 0.1;
        y1[i] = sin(sqrt(x1[i] + 1)) * sqrt(x1[i]);
        y2[i] = sin(sqrt(x2[i] + 1)) * sqrt(x2[i]);
    }

    /* Initialize the splines */
    struct strooklat spline1 = {x1, N};
    struct strooklat spline2 = {x2, N};
    init_strooklat_spline(&spline1, 100);
    init_strooklat_spline(&spline2, 100);

    /* Initialize GSL splines to compare against */
    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    const gsl_interp_type *t = gsl_interp_linear;
    gsl_spline *gsl_spline = gsl_spline_alloc(t, N);

    gsl_spline_init(gsl_spline, x1, y1, N);

    /* Test the splines */
    for (int i = 1; i < N - 1; i++) {
        double X1 = i * 0.1 + 0.03;
        double X2 = 100 - i * 0.1 - 0.03;

        double Y1 = strooklat_interp(&spline1, y1, X1);
        double Y2 = strooklat_interp(&spline2, y2, X2);

        double gslY1 = gsl_spline_eval (gsl_spline, X1, acc);
        double gslY2 = gsl_spline_eval (gsl_spline, X2, acc);

        double err1 = (Y1 - gslY1)/(Y1 + gslY1);
        double err2 = (Y2 - gslY2)/(Y2 + gslY2);

        assert(fabs(err1) < 1e-6);
        assert(fabs(err2) < 1e-6);

        // printf("%f %f %e %e\n", X1, X2, err1, err2);
    }

    /* Start the timer */
    struct timeval time_stop, time_start;
    gettimeofday(&time_start, NULL);

    int M = 1e8;
    double sum = 0;
    for (int i=0; i<M; i++) {
        double x = (double) 99 * rand()/RAND_MAX;
        double Y = strooklat_interp(&spline1, y1, x);
        sum += Y;
    }

    printf("Strooklat mean: %e\n", sum/M);

    /* End the timer */
    gettimeofday(&time_stop, NULL);
    long unsigned microsec = (time_stop.tv_sec - time_start.tv_sec) * 1000000
                           + time_stop.tv_usec - time_start.tv_usec;
    printf("Strooklat time: %.5f s\n", microsec/1e6);

    /* Reset the timer */
    gettimeofday(&time_start, NULL);

    sum = 0;
    for (int i=0; i<M; i++) {
        double x = (double) 99 * rand()/RAND_MAX;
        double Y = gsl_spline_eval (gsl_spline, x, acc);
        sum += Y;
    }

    printf("GSL mean: %e\n", sum/M);

    /* End the timer */
    gettimeofday(&time_stop, NULL);
    microsec = (time_stop.tv_sec - time_start.tv_sec) * 1000000
                           + time_stop.tv_usec - time_start.tv_usec;
    printf("GSL time: %.5f s\n", microsec/1e6);


    /* Free the splines */
    free_strooklat_spline(&spline1);
    free_strooklat_spline(&spline2);

    /* Free the GSL splines */
    gsl_spline_free(gsl_spline);
    gsl_interp_accel_free(acc);

    /* Free the data */
    free(x1);
    free(x2);
    free(y1);
    free(y2);

    return 0;
}
