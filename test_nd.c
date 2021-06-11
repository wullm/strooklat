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
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

/* GSL interpolation library */
#include <gsl/gsl_spline.h>
#include <gsl/gsl_spline2d.h>

/* Strooklat interpolation library */
#include "strooklat_inline.h"

int main() {
    /* Allocate memory for some data points */
    int N1 = 7000;
    int N2 = 1000;
    double *x1 = malloc(N1 * sizeof(double));
    double *x2 = malloc(N2 * sizeof(double));
    double *y1 = malloc(N1 * sizeof(double));
    double *y2 = malloc(N2 * sizeof(double));
    double *z = malloc(N1 * N2 * sizeof(double));

    /* Generate the first data set */
    for (int i = 0; i < N1; i++) {
        x1[i] = i * 100.0 / N1;
        y1[i] = sin(sqrt(x1[i] + 1)) * sqrt(x1[i]);
    }

    /* Generate the second data set */
    for (int i = 0; i < N2; i++) {
        x2[i] = 100 + i * 100.0 / N2;
        y2[i] = sin(sqrt(x2[i] + 1)) * sqrt(x2[i]);
    }

    /* Compute the exterior product (row major format) */
    for (int i = 0; i < N1; i++) {
        for (int j = 0; j < N2; j++) {
            z[j + N2 * i] = y1[i] * y2[j];
        }
    }

    /* Initialize the splines */
    struct strooklat spline1 = {x1, N1};
    struct strooklat spline2 = {x2, N2};
    init_strooklat_spline(&spline1, 100);
    init_strooklat_spline(&spline2, 100);

    /* Prepare the n-dimensional interpolation */
    struct strooklat *splines[2] = {&spline1, &spline2};

    /* Initialize GSL splines to compare against */
    gsl_interp_accel *acc1 = gsl_interp_accel_alloc();
    gsl_interp_accel *acc2 = gsl_interp_accel_alloc();
    const gsl_interp_type *t = gsl_interp_linear;
    gsl_spline *gsl_spline1 = gsl_spline_alloc(t, N1);
    gsl_spline *gsl_spline2 = gsl_spline_alloc(t, N2);

    gsl_spline_init(gsl_spline1, x1, y1, N1);
    gsl_spline_init(gsl_spline2, x2, y2, N2);

    /* Initialise 2D GSL splines to compare against */
    gsl_spline2d *gsl_spline_2d =
        gsl_spline2d_alloc(gsl_interp2d_bilinear, N2, N1);
    gsl_interp_accel *x1acc = gsl_interp_accel_alloc();
    gsl_interp_accel *x2acc = gsl_interp_accel_alloc();

    /* Flip the axes, because GSL assumes column major format */
    gsl_spline2d_init(gsl_spline_2d, x2, x1, z, N2, N1);

    /* Test the splines */
    int N = 1000;
    for (int i = 1; i < N - 1; i++) {
        double X1 = i * 0.1 + 0.07;
        double X2 = 100 + i * 0.1 - 0.07;
        double X[2] = {X1, X2};

        /* Compute the exact values */
        double Y1_exact = sin(sqrt(X1 + 1)) * sqrt(X1);
        double Y2_exact = sin(sqrt(X2 + 1)) * sqrt(X2);
        double Z_exact = Y1_exact * Y2_exact;

        /* Compute the strooklat interpolated values */
        double Y1 = strooklat_interp(&spline1, y1, X1);
        double Y2 = strooklat_interp(&spline2, y2, X2);
        double Z = strooklat_interp_nd(2, splines, z, X);
        double Z_alt = strooklat_interp_2d(&spline1, &spline2, z, X1, X2);

        /* Compute the GSL interpolated values (flip the axes for 2D) */
        double gslY1 = gsl_spline_eval(gsl_spline1, X1, acc1);
        double gslY2 = gsl_spline_eval(gsl_spline2, X2, acc2);
        double gslZ = gsl_spline2d_eval(gsl_spline_2d, X2, X1, x1acc, x2acc);

        /* Compute the relative errors between strooklat and GSL */
        double err1 = (Y1 - gslY1) / (Y1 + gslY1);
        double err2 = (Y2 - gslY2) / (Y2 + gslY2);
        double err3 = (Z - gslZ) / (Z + gslZ);
        double err4 = (Z_alt - gslZ) / (Z_alt + gslZ);

        /* Compute the relative errors with the exact values */
        double err_exact1 = (Y1 - Y1_exact) / (Y1 + Y1_exact);
        double err_exact2 = (Y2 - Y2_exact) / (Y2 + Y2_exact);
        double err_exact3 = (Z - Z_exact) / (Z + Z_exact);
        double err_exact4 = (Z_alt - Z_exact) / (Z_alt + Z_exact);

        // printf("%f %f %e %e %e %e %e %e\n", X1, X2, err1, err2, err3,
        // err_exact1, err_exact2, err_exact3);

        assert(fabs(err1) < 1e-12);
        assert(fabs(err2) < 1e-12);
        assert(fabs(err3) < 1e-12);
        assert(fabs(err4) < 1e-12);
        assert(fabs(err_exact1) < 1e-3);
        assert(fabs(err_exact2) < 1e-3);
        assert(fabs(err_exact3) < 1e-3);
        assert(fabs(err_exact4) < 1e-3);
    }

    /* Start the timer */
    struct timeval time_stop, time_start;
    gettimeofday(&time_start, NULL);

    int M = 1e8;
    double sum = 0;
    for (int i = 0; i < M; i++) {
        double X1 = (double)99 * rand() / RAND_MAX;
        double X2 = 100.0 + (double)99 * rand() / RAND_MAX;
        // double X[2] = {X1, X2};
        // double Z = strooklat_interp_nd(2, splines, z, X);
        double Z = strooklat_interp_2d(&spline1, &spline2, z, X1, X2);
        sum += Z;
    }

    printf("Strooklat mean: %e\n", sum / M);

    /* End the timer */
    gettimeofday(&time_stop, NULL);
    long unsigned microsec = (time_stop.tv_sec - time_start.tv_sec) * 1000000 +
                             time_stop.tv_usec - time_start.tv_usec;
    printf("Strooklat time: %.5f s\n", microsec / 1e6);

    /* Reset the timer */
    gettimeofday(&time_start, NULL);

    sum = 0;
    for (int i = 0; i < M; i++) {
        double X1 = (double)99 * rand() / RAND_MAX;
        double X2 = 100.0 + (double)99 * rand() / RAND_MAX;
        double Z = gsl_spline2d_eval(gsl_spline_2d, X2, X1, x1acc, x2acc);
        sum += Z;
    }

    printf("GSL mean: %e\n", sum / M);

    /* End the timer */
    gettimeofday(&time_stop, NULL);
    microsec = (time_stop.tv_sec - time_start.tv_sec) * 1000000 +
               time_stop.tv_usec - time_start.tv_usec;
    printf("GSL time: %.5f s\n", microsec / 1e6);

    /* Free the splines */
    free_strooklat_spline(&spline1);
    free_strooklat_spline(&spline2);

    /* Free the GSL splines */
    gsl_spline_free(gsl_spline1);
    gsl_spline_free(gsl_spline2);
    gsl_interp_accel_free(acc1);
    gsl_interp_accel_free(acc2);

    /* Free the 2D GSL splines */
    gsl_spline2d_free(gsl_spline_2d);
    gsl_interp_accel_free(x1acc);
    gsl_interp_accel_free(x2acc);

    /* Free the data */
    free(x1);
    free(x2);
    free(y1);
    free(y2);
    free(z);

    return 0;
}
