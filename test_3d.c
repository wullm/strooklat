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
#include "strooklat.h"

int main() {
    /* Allocate memory for some data points */
    int N1 = 700;
    int N2 = 200;
    int N3 = 400;
    double *x1 = malloc(N1 * sizeof(double));
    double *x2 = malloc(N2 * sizeof(double));
    double *x3 = malloc(N3 * sizeof(double));
    double *y1 = malloc(N1 * sizeof(double));
    double *y2 = malloc(N2 * sizeof(double));
    double *y3 = malloc(N3 * sizeof(double));
    double *z = malloc(N1 * N2 * N3 * sizeof(double));

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

    /* Generate the second data set */
    for (int i = 0; i < N3; i++) {
        x3[i] = 100 + i * 100.0 / N3;
        y3[i] = sin(sqrt(x3[i] + 1)) * sqrt(x3[i]);
    }

    /* Compute the exterior product (row major format) */
    for (int i = 0; i < N1; i++) {
        for (int j = 0; j < N2; j++) {
            for (int k = 0; k < N3; k++) {
                z[k + j * N3 + i * N2 * N3] = y1[i] * y2[j] * y3[k];
            }
        }
    }

    /* Initialize the splines */
    struct strooklat spline1 = {x1, N1};
    struct strooklat spline2 = {x2, N2};
    struct strooklat spline3 = {x3, N3};
    init_strooklat_spline(&spline1, 100);
    init_strooklat_spline(&spline2, 100);
    init_strooklat_spline(&spline3, 100);

    /* Prepare the n-dimensional interpolation */
    struct strooklat const *const splines[3] = {&spline1, &spline2, &spline3};

    /* Test the splines */
    int N = 1000;
    for (int i = 1; i < N - 1; i++) {
        double X1 = i * 0.1 + 0.07;
        double X2 = 100 + i * 0.1 - 0.07;
        double X3 = 100 + i * 0.1 - 0.07;
        double X[3] = {X1, X2, X3};

        /* Compute the strooklat interpolated values */
        double Z = strooklat_interp_nd(3, splines, z, X);
        double Z_alt = strooklat_interp_3d(&spline1, &spline2, &spline3, z, X1, X2, X3);

        /* Compute the relative errors with the exact values */
        double err = (Z_alt - Z) / (Z_alt + Z);

        assert(fabs(err) < 1e-12);
    }

    printf("Check complete.\n");

    /* Start the timer */
    struct timeval time_stop, time_start;
    gettimeofday(&time_start, NULL);

    int M = 1e8;
    double sum = 0;
    for (int i = 0; i < M; i++) {
        double X1 = (double)99 * rand() / RAND_MAX;
        double X2 = 100.0 + (double)99 * rand() / RAND_MAX;
        double X3 = 100.0 + (double)99 * rand() / RAND_MAX;
        double Z = strooklat_interp_3d(&spline1, &spline2, &spline3, z, X1, X2, X3);
        sum += Z;
    }

    printf("Strooklat 3D mean: %e\n", sum / M);

    /* End the timer */
    gettimeofday(&time_stop, NULL);
    long unsigned microsec = (time_stop.tv_sec - time_start.tv_sec) * 1000000 +
                             time_stop.tv_usec - time_start.tv_usec;
    printf("Strooklat 3D time: %.5f s\n", microsec / 1e6);

    /* Reset the timer */
    gettimeofday(&time_start, NULL);

    sum = 0;
    for (int i = 0; i < M; i++) {
        double X1 = (double)99 * rand() / RAND_MAX;
        double X2 = 100.0 + (double)99 * rand() / RAND_MAX;
        double X3 = 100.0 + (double)99 * rand() / RAND_MAX;
        double X[3] = {X1, X2, X3};
        double Z = strooklat_interp_nd(3, splines, z, X);
        sum += Z;
    }

    printf("Strooklat ND mean: %e\n", sum / M);

    /* End the timer */
    gettimeofday(&time_stop, NULL);
    microsec = (time_stop.tv_sec - time_start.tv_sec) * 1000000 +
               time_stop.tv_usec - time_start.tv_usec;
    printf("Strooklat ND time: %.5f s\n", microsec / 1e6);

    /* Free the splines */
    free_strooklat_spline(&spline1);
    free_strooklat_spline(&spline2);
    free_strooklat_spline(&spline3);

    /* Free the data */
    free(x1);
    free(x2);
    free(x3);
    free(y1);
    free(y2);
    free(y3);
    free(z);

    return 0;
}
