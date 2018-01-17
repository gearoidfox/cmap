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

#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#include"pdb.h"

/**
 * Find euclidean distance between two atoms located at coordinates (x1, y1, z1)
 * and (x2, y2, z2), respectively.
 */
double
euclid3d(double x1, double y1, double z1, double x2, double y2, double z2)
{
        return sqrt(pow(x1 - x2, 2) + pow(y1 - y2, 2) + pow(z1 - z2, 2));
}


/**
 * Create a distance matrix from atomic coordinates coordinates
 *
 * @coords: a structure containing 3D coordinates of a group of atoms
 *
 * returns a pointer to distmat structure with the pairwise euclidean distances
 * of the atoms. 
 *
 * Allocates memory for this structure
 *
 */
struct distmat *
calculate_distmat(struct coords cs)
{
        double **dist = NULL;
        int i, j;
        struct distmat * dm = NULL;

        if(cs.coords == NULL) return NULL;
        if(cs.nres == 0) return NULL;

        dm = (struct distmat *) malloc(sizeof(struct distmat));
        if(dm == NULL) return NULL;
        dm->mat = NULL;
        dm->nres = cs.nres;

        dist = (double**) malloc((cs.nres - 1) * sizeof(double*));
        if(dist == NULL) goto cdm_error_cleanup;

        for(i = 0; i < cs.nres - 1; i++){
                dist[i] = NULL;
                dist[i] = (double*) malloc((cs.nres - 1 - i) * sizeof(double));
                if(dist[i] == NULL) goto cdm_error_cleanup;
        }
        for(i = 0; i < cs.nres - 1; i++){
                for(j = i + 1; j < cs.nres; j++){
                        if (cs.coords[i] == NULL || cs.coords[j] == NULL){
                            dist[i][j - i - 1] = 999;
                            continue;
                        }
                        dist[i][j - i - 1] = euclid3d(cs.coords[i][0],
                            cs.coords[i][1], cs.coords[i][2], cs.coords[j][0],
                            cs.coords[j][1], cs.coords[j][2]);
                }
        }
        dm->mat = dist;
        return dm;
        /* 
         * Avoid leaking memory by freeing all allocated memory before
         * return NULL. Needed if any call to malloc after the first fails.
         */
        cdm_error_cleanup:
        if(dm != NULL){
                if(dm->mat != NULL) free(dm->mat);
                free(dm);
                dm = NULL;
        }
        if(dist != NULL){
                for(i = 0; i < cs.nres; i++){
                        if(dist[i] != NULL){
                                free(dist[i]);
                                dist[i] = NULL;
                        }
                }
                free(dist);
                dist = NULL;
        }
        return NULL;
}

/**
 *
 */
double
getdist(struct distmat dm, int i, int j)
{
        if(i==j)
                return 0;
        if(i < j)
                return fabs(dm.mat[i][j-i-1]);
        else
                return fabs(dm.mat[j][i-j-1]);
}



/* 
 * freedm: free a distance matrix allocated by calculate_distmat().
 */
void
freedm(struct distmat *dm)
{
        int i;
        for(i = 0; i < dm->nres - 1 ; i++){
                if(dm->mat[i] != NULL)
                        free(dm->mat[i]) ;
                dm->mat[i] = NULL;
        }
        free(dm->mat);
        free(dm);
}

/**
 * getcoords: read alpha carbon co-ordinates from a PDB file
 *
 * @filename: string containing path to Protein Data Bank file
 * @chain: chain identifier
 */
struct coords *
getcoords(char* filename, char target_chain){
        FILE *fp;
        char buffer[1028];
        char recname[7];
        char name[5];
        char x[9], y[9], z[9];
        char resseq[5];
        char chain;
        int i;
        int n;
        int nres = 0;
        double **coords=NULL;
        struct coords *cs;
        
        fp = fopen(filename, "r");
        if (fp == NULL) return NULL;

        cs = malloc(sizeof (struct coords));
        if (cs == NULL) return NULL;

        /* Initial pass through file to count residues in chain of interest */
        while(fgets(buffer, 1028, fp)!= NULL){
                strncpy(recname, buffer, 6); 
                recname[6] = '\0';
                if(strcmp("ATOM  ", recname) == 0) {
                        chain = buffer[21];
                        if(chain == target_chain){
                                strncpy(resseq, buffer+22, 4);
                                resseq[4] =  '\0';
                                n = atoi(resseq);
                                if(n > nres)
                                        nres = n;
                        }
                }
        }
        if(nres == 0) goto gc_error_cleanup;
        
        coords = (double**) malloc(nres * sizeof(double*));
        if(coords == NULL) goto gc_error_cleanup;
        for(i = 0; i < nres; i++ ){
                coords[i] = NULL;
        }

        /* Second pass to record co-ordinates into matrix coords */ 
        rewind(fp);
        while(fgets(buffer, 1028, fp)!= NULL){
                strncpy(recname, buffer, 6); 
                recname[6] = '\0';
                if(strcmp("ATOM  ", recname) == 0) {
                        strncpy(name, buffer+12, 4);
                        chain = buffer[21];
                        name[4] = '\0';
                        if(strcmp(" CA ", name) == 0 && chain == target_chain){
                                strncpy(resseq, buffer+22, 4);
                                resseq[4] =  '\0';
                                strncpy(x, buffer+30, 8);
                                strncpy(y, buffer+38, 8);
                                strncpy(z, buffer+46, 8);
                                x[8] = '\0';
                                y[8] = '\0';
                                z[8] = '\0';

                                n = atoi(resseq);
                                if (n <= 0 || n > nres){ // Bad input data; would segfault
                                        goto gc_error_cleanup;
                                }
                                if(coords[n-1] != NULL){
                                        /* 
                                         * We already stored coordinates for
                                         * this atom. (i.e. the file contains
                                         * an alternate position for this 
                                         * residue. Not handled right now, and
                                         * we simply update the coordinates to 
                                         * the ones we just found.
                                         */
                                }
                                else coords[n-1] = (double*) malloc(3 * sizeof (double));
                                if(coords[n-1] == NULL) goto gc_error_cleanup;
                                coords[n-1][0] = atof(x);
                                coords[n-1][1] = atof(y);
                                coords[n-1][2] = atof(z);
                        }
                }
        }
        cs->nres = nres;
        cs->coords = coords;
        fclose(fp);
        return cs;

        /* 
         * Avoid leaking memory by freeing all allocated memory before
         * returning NULL. Needed if any call to malloc after the first fails.
         */
        gc_error_cleanup:
        if(cs) free(cs);
        if(coords){
                for(i=0; i<nres; i++ ){
                        if(coords[i]){
                                free(coords[i]);
                                coords[i] = NULL;
                        }
                }
                free(coords);
        }
        return NULL;
}

/**
 * freecoords: free a struct coords allocated by etcoords()
 */
void freecoords(struct coords *cs){
        int i;
        if(cs == NULL) return;
        for(i = 0; i < cs->nres; i++){
                if(cs->coords[i] != NULL){
                        free(cs->coords[i]);
                        cs->coords[i] = NULL;
                }
        }
        free(cs->coords);
        cs->coords = NULL;
        free(cs);
}

