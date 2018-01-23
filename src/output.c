#include<pdb.h>
#include<stdio.h>
#include "pdb.h"

/**
 * write_contacts: write a list of contacts to a text file
 *
 * @fp: file pointer open for writing
 * @dm: distance matrix
 * @threshold: distance threshold used to calculate contacts
 */
void
write_contacts(FILE *fp, struct distmat dm, double threshold)
{
        int i, j;

        if(fp == NULL){
                return;
        }
        fprintf(fp, "# cmap v0.1\n");
        if(dm.source_filename != NULL)
                fprintf(fp, "# source file: %s\n", dm.source_filename);
        if(dm.source_chain != '\0')
                fprintf(fp, "# source chain: %c\n", dm.source_chain);
        if(dm.sequence != NULL)
                fprintf(fp, "# sequence: %s\n", dm.sequence);
        fprintf(fp, "# threshold: %f\n", threshold);
        for(i = 0; i < dm.nres - 1; i++){
                for(j = i + 1; j < dm.nres; j++){
                        if(getdist(dm, i, j) < threshold)
                                fprintf(fp, "%d\t%d\n", i+1, j+1);
                }
        }
        return;
}

/**
 * write_eps: write encapsulated postscript diagram of a contact map
 *
 * @fp: file pointer open for writing
 * @dm: distance matrix
 * @threshold: distance threshold used to calculate contacts
 */
void
write_eps(FILE *fp, struct distmat dm, double threshold)
{
        int i, j;
        int xmax;
        int ymax;

        if(fp == NULL){
                return;
        }

        xmax = dm.nres * 5 + 20;
        ymax = dm.nres * 5 + 20;

        fprintf(fp, "%%!PS-Adobe EPSF-3.0\n");
        fprintf(fp, "%%%%BoundingBox: 0 0 %d %d\n", xmax, ymax);
        fprintf(fp, "1 1 scale\n");

        fprintf(fp, "newpath\n"
                    "10 10 moveto\n"
                    "0 %d rlineto\n"
                    "%d 0 rlineto\n"
                    "0 -%d rlineto\n"
                    "-%d 0 rlineto\n"
                    "closepath\n"
                    "1 setlinewidth\n"
                    "stroke\n",
                    dm.nres * 5 + 1, dm.nres * 5 + 1, dm.nres * 5 + 1, 
                    dm.nres * 5 + 1);

        fprintf(fp, "/square {\n"
                    "/y exch def\n"
                    "/x exch def\n"
                    "gsave\n"
                    "newpath\n"
                    "y x moveto\n"
                    "0 5 rlineto\n"
                    "5 0 rlineto\n"
                    "0 -5 rlineto\n"
                    "-5 0 rlineto\n"
                    "closepath\n"
                    "fill\n"
                    "}\n"
                    "def\n");

        for(i = 0; i < dm.nres; i++){
                for(j = 0; j < dm.nres; j++){
                       if(getdist(dm, i, j) < threshold) 
                               fprintf(fp, "%d %d square\n", 
                                               (dm.nres * 5 + 6) - i*5,
                                               11 + j * 5);
                }
        }
        return;
}
