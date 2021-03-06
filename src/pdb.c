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
 * euclid3d: Find euclidean distance between two atoms located at coordinates
 * (x1, y1, z1) and (x2, y2, z2), respectively.
 */
double
euclid3d(double x1, double y1, double z1, double x2, double y2, double z2)
{
        return sqrt(pow(x1 - x2, 2) + pow(y1 - y2, 2) + pow(z1 - z2, 2));
}


/**
 * calculate_distmat: Create a distance matrix from atomic coordinates
 *
 * @coords: a structure containing 3D coordinates of a group of atoms
 *
 * Returns a pointer to distmat structure storing the pairwise euclidean
 * distances of the atoms. 
 *
 * Allocates memory for this structure, which should be freed with function
 * freedm()
 *
 * Returns NULL if a distance matrix cannot be calculated/allocated. 
 *
 * Use getdist() on the distance matrix to query distances between residues
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

        dm = malloc(sizeof(*dm));
        if(dm == NULL) return NULL;
        dm->mat = NULL;
        dm->nres = cs.nres;

        dm->source_chain = cs.source_chain;

        dm->source_filename = NULL;
        if(cs.source_filename != NULL){
                i = strlen(cs.source_filename) + 1;
                dm->source_filename = malloc(i * sizeof(*(dm->source_filename)));
                if(dm->source_filename == NULL) goto cdm_error_cleanup;
                strncpy(dm->source_filename, cs.source_filename, i);
        }
        
        dm->sequence = NULL;
        if(cs.sequence != NULL){
                i = strlen(cs.sequence) + 1;
                dm->sequence = malloc(i * sizeof(*(dm->sequence)));
                if(dm->sequence == NULL) goto cdm_error_cleanup;
                strncpy(dm->sequence, cs.sequence, i);
        }

        dist = malloc((cs.nres - 1) * sizeof(*dist));
        if(dist == NULL) goto cdm_error_cleanup;

        /* Allocate "triangular" matrix 
         * Distance matrix has diagonal symmetry, so only store one half
         */
        for(i = 0; i < cs.nres - 1; i++){
                dist[i] = NULL;
                dist[i] = malloc((cs.nres - 1 - i) * sizeof(*dist[i]));
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
         * In the event of an error occurring inside this function,
         * avoid leaking memory by freeing all allocated memory before
         * return NULL. Needed if any call to malloc after the first fails.
         */
        cdm_error_cleanup:
        if(dm != NULL){
                if(dm->mat != NULL) free(dm->mat);
                if(dm->source_filename != NULL) free(dm->source_filename);
                if(dm->sequence != NULL) free(dm->sequence);
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
 * getdist: extract distances from a struct distmat object
 * 
 * @dm: distance matrix
 * @i:  index of a residue in protein chain
 * @j:  index of a second residue
 *
 * Returns the distance between residues @i and @j as stored in @dm
 *
 * dm.mat is stored in a triangular form, which is why we can't just read
 * off dm.mat[i][j] directly.
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
        if(dm->source_filename != NULL) free(dm->source_filename);
        dm->source_filename = NULL;
        if(dm->sequence != NULL) free(dm->sequence);
        dm->sequence = NULL;
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
        /* common record elements */
        char recname[7];
        char chain;
        /* in seqres records */
        char tlc[4];
        char numres[5];
        /* in atom records */
        char name[5];
        char resname[4];
        char x[9], y[9], z[9];
        char resseq[5];

        int i;
        int n;
        int nres = 0;
        double **coords=NULL;
        struct coords *cs;
        char *seq_ptr;
        
        fp = fopen(filename, "r");
        if (fp == NULL) return NULL;

        cs = malloc(sizeof (*cs));
        if (cs == NULL) return NULL;

        cs->source_chain = target_chain;

        /* Allocate memory for/store input filename */
        n = strlen(filename) + 1;
        cs->source_filename = NULL;
        cs->source_filename = malloc(n * sizeof(*(cs->source_filename)));
        if(cs->source_filename == NULL) goto gc_error_cleanup;
        strncpy(cs->source_filename, filename, n);

        cs->sequence = NULL;
        n = 0;
        /* Parse SEQRES records from PDB header to find primary seq info */
        while(fgets(buffer, 1028, fp)!= NULL){
                strncpy(recname, buffer, 6); 
                recname[6] = '\0';
                if(strncmp("SEQRES", recname, 6) == 0) {
                        chain = buffer[11];
                        /* At first SEQRES for our chain -- read # of residues*/
                        if(chain == target_chain && cs->sequence == NULL){
                                strncpy(numres, buffer+13, 4);
                                numres[4] =  '\0';
                                nres = atoi(numres);
                                cs->sequence = malloc((nres + 1) * sizeof(*(cs->sequence)));
                                if(cs->sequence == NULL) goto gc_error_cleanup;
                                memset(cs->sequence, 0, nres + 1);
                                seq_ptr = cs->sequence;
                        }
                        /* At all SEQRES for our chain -- read primary sequence*/
                        if(chain == target_chain){
                            n += read_seqres_line(cs->sequence + n, buffer, nres - n);
                            if(n == nres) break;
                        }
                }
        }
        if(cs->sequence != NULL)
                cs->sequence[nres] = '\0';
        if(nres == 0) goto gc_error_cleanup;
        
        coords = malloc(nres * sizeof(*coords));
        if(coords == NULL) goto gc_error_cleanup;
        for(i = 0; i < nres; i++ ){
                coords[i] = NULL;
        }

        /* Second pass to record co-ordinates into matrix coords */ 
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

                                /* found an atom at a beyond the terminus of
                                 * the chain recorded in the header.
                                 * No memory allocated to store this info
                                 */
                                if (n <= 0 || n > nres){ 
                                        fprintf(stderr, "WARNING: unexpected ATOM records found in chain %c [length %d].\n%s", chain, nres, buffer);
                                        continue;
                                        /* goto gc_error_cleanup; */
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
                                else coords[n-1] = malloc(3 * sizeof(*coords[n-1]));
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

        /* In the event of an error in this function: 
         * Avoid leaking memory by freeing all allocated memory before
         * returning NULL. Needed if any call to malloc after the first fails.
         */
        gc_error_cleanup:
        if(cs){
                if(cs->source_filename != NULL) free(cs->source_filename);
                if(cs->sequence != NULL) free(cs->sequence);
                free(cs);
        }
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
 * freecoords: free a struct coords allocated by getcoords()
 * @cs: pointer to struct coords which is to be freed
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
        if(cs->source_filename != NULL) free(cs->source_filename);
        cs->source_filename = NULL;
        if(cs->sequence != NULL) free(cs->sequence);
        cs->sequence = NULL;
        free(cs);
}


/**
 * one_letter_code
 * convert three letter amino acid codes to one letter codes
 *
 * @three_letter_code: pointer to string containing a three-letter amino acid
 *                     code in UPPERCASE
 *
 * Will also process one and two letter RNA/DNA codes valid in PDB SEQRES 
 * records.
 *
 */
char
one_letter_code(char *three_letter_code)
{
        if(three_letter_code == NULL) return 'X';
        /* Amino acids */
        if(0 == strncmp(three_letter_code, "ALA", 3)) return 'A';
        if(0 == strncmp(three_letter_code, "ASX", 3)) return 'B';
        if(0 == strncmp(three_letter_code, "CYS", 3)) return 'C';
        if(0 == strncmp(three_letter_code, "ASP", 3)) return 'D';
        if(0 == strncmp(three_letter_code, "GLU", 3)) return 'E';
        if(0 == strncmp(three_letter_code, "PHE", 3)) return 'F';
        if(0 == strncmp(three_letter_code, "GLY", 3)) return 'G';
        if(0 == strncmp(three_letter_code, "HIS", 3)) return 'H';
        if(0 == strncmp(three_letter_code, "ILE", 3)) return 'I';
        if(0 == strncmp(three_letter_code, "LYS", 3)) return 'K';
        if(0 == strncmp(three_letter_code, "LEU", 3)) return 'L';
        if(0 == strncmp(three_letter_code, "MET", 3)) return 'M';
        if(0 == strncmp(three_letter_code, "ASN", 3)) return 'N';
        if(0 == strncmp(three_letter_code, "PRO", 3)) return 'P';
        if(0 == strncmp(three_letter_code, "GLN", 3)) return 'Q';
        if(0 == strncmp(three_letter_code, "ARG", 3)) return 'R';
        if(0 == strncmp(three_letter_code, "SER", 3)) return 'S';
        if(0 == strncmp(three_letter_code, "THR", 3)) return 'T';
        if(0 == strncmp(three_letter_code, "VAL", 3)) return 'V';
        if(0 == strncmp(three_letter_code, "TRP", 3)) return 'W';
        if(0 == strncmp(three_letter_code, "XAA", 3)) return 'X';
        if(0 == strncmp(three_letter_code, "TYR", 3)) return 'Y';
        if(0 == strncmp(three_letter_code, "GLX", 3)) return 'Z';
        /* Ribonucleotides */ 
        if(0 == strncmp(three_letter_code, "  A", 3)) return 'A';
        if(0 == strncmp(three_letter_code, "  C", 3)) return 'C';
        if(0 == strncmp(three_letter_code, "  G", 3)) return 'G';
        if(0 == strncmp(three_letter_code, "  U", 3)) return 'U';
        /* Deoxyribonucleotides */ 
        if(0 == strncmp(three_letter_code, " DA", 3)) return 'A';
        if(0 == strncmp(three_letter_code, " DC", 3)) return 'C';
        if(0 == strncmp(three_letter_code, " DG", 3)) return 'G';
        if(0 == strncmp(three_letter_code, " DT", 3)) return 'T';
        /* Default */
        return 'X';
}

/**
 * read_seqres_line: parse a seqres record into a primary sequence string
 *                   fragment
 *
 * @out_buffer: pointer to allocated character array of length >= n
 * @line:       pointer to string containing a PDB seqres record
 * @n:          maximum number of characters to write to out_buffer
 *
 * Returns the number of characters written to @out_buffer.
 *
 * Each seqres record stores part of the primary sequence of a chain in 
 * the PDB file. Read the three-letter code sequence in the record pointed to
 * by @line, and convert to a one-letter-code string which is then stored at
 * @out_buffer.
 *
 * This function does not write a null terminator to out_buffer.
 *
 * EXAMPLE:
 * If @n is 13, and @line points to the string: 
 * "SEQRES   1 A   21  GLY ILE VAL GLU GLN CYS CYS THR SER ILE CYS SER LEU"
 * then this function will write 
 * "GIVEQCCTSICSL"
 * to the location pointed to by @out_buffer, and return the value 13.
 *
 */
int
read_seqres_line(char* out_buffer, char *line, int n)
{
        char recname[7];
        char sernum[4];
        char resname[4];
        char numres[5];

        int i;
        int total_res;
        int serial;
        int recorded_res;
        int offset;

        if(line == NULL) return 0;
        if(out_buffer == NULL) return 0;
        if(n <= 0) return 0;
        
        strncpy(recname, line, 6);
        recname[6] = '\0';
        if(0 != strncmp("SEQRES", recname, 6)) return 0;

        strncpy(sernum, line + 7, 3);
        sernum[3] = '\0';
        serial = atoi(sernum);
       
        strncpy(numres, line+13, 4);
        numres[4] = '\0';
        total_res = atoi(numres);

        offset =  (serial - 1) * 13;        
        recorded_res = total_res - offset > 13 ? 13 : total_res - offset;

        i=0;
        while(i < recorded_res && i < n){
                strncpy(resname, line + 19 + (4 * i), 3);
                resname[3] = '\0';
                out_buffer[i] = one_letter_code(resname);
                i++;
        }
        return i;
}
