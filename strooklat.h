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

#ifndef STROOKLAT_SPLINE_H
#define STROOKLAT_SPLINE_H

/* Spline struct with lookup table */
struct strooklat {
    /* The x values of the data points to be interpolated */
    const double *x;
    /* The number of data points */
    const int size;
    /* Whether the x-array is ascending */
    int ascend;
    /* The last index returned */
    int last_index;

    struct lookup {
        /* A lookup table for faster searches */
        double *lookup_table;
        /* Size of the lookup table */
        int lookup_table_size;
    } lookup;
};

int init_strooklat_spline(struct strooklat *spline, int index_size);
int free_strooklat_spline(struct strooklat *spline);
int strooklat_find_x(struct strooklat *spline, double x, int *ind, double *u);
double strooklat_interp_index(struct strooklat *spline, double *y, int ind,
                              double u);
double strooklat_interp(struct strooklat *spline, double *y, double x);

#endif
