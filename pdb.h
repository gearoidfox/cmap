/* 
 *  Copyright 2018 Gearoid Fox
 * 
 *  This file is part of cmap.
 *
 *  cmap is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  cmap is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with cmap.  If not, see <http://www.gnu.org/licenses/>. 
 *
 */

#ifndef CMAP_PDB_H_
#define CMAP_PDB_H_

struct coords{
        double **coords;
        char *source_filename;
        int nres;
        char source_chain;
};

struct distmat{
        double **mat;
        char *source_filename;
        int nres;
        char source_chain;
};

struct distmat * calculate_distmat(struct coords cs);
double euclid3d(double x1, double y1, double z1, double x2, double y2, double z2);
void freecoords(struct coords *cs);
void freedm(struct distmat *dm);
struct coords * getcoords(char* filename, char chain);
double getdist(struct distmat dm, int i, int j);

#endif // CMAP_PDB_H_
