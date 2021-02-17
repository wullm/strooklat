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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

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
        x1[i] = 100 - i * 0.1;
        x2[i] = i * 0.1;
        y1[i] = sin(sqrt(x1[i] + 1)) * sqrt(x1[i]);
        y2[i] = sin(sqrt(x2[i] + 1)) * sqrt(x2[i]);
    }

    // /* Initialize the splines */
    struct strooklat spline1 = {x1, N};
    struct strooklat spline2 = {x2, N};
    init_strooklat_spline(&spline1, 100);
    init_strooklat_spline(&spline2, 100);

    /* Test the splines */
    for (int i = 0; i < N * 1; i++) {
        double X1 = 100 - i * 0.1 + 0.03;
        double X2 = i * 0.1 + 0.03;

        double Y1 = strooklat_interp(&spline1, y1, X1);
        double Y2 = strooklat_interp(&spline2, y2, X2);

        printf("%f %f\n", X1, Y1);
    }

    /* Free the splines */
    free_strooklat_spline(&spline1);
    free_strooklat_spline(&spline2);

    /* Free the data */
    free(x1);
    free(x2);
    free(y1);
    free(y2);

    return 0;
}
