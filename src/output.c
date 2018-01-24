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
        fprintf(fp, "# cmap v%s\n", PACKAGE_VERSION);
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
        float xmax;
        float ymax;
        float oma = 0;
        float boxw;
        float framewidth = 1;
        if(fp == NULL){
                return;
        }
        /* Draw about six sq inches, irrespective of chain length
         * Might be bad for long chains.
         */
        boxw = (72. * 6. / dm.nres ) ;

        xmax = 2 * oma + 2 * framewidth + dm.nres * boxw;
        ymax = xmax + 40;

        fprintf(fp, "%%!PS-Adobe EPSF-3.0\n");
        fprintf(fp, "%%%%BoundingBox: %f %f %f %f\n", oma, oma, xmax, ymax);
        fprintf(fp, "/Courier\n12 selectfont\n");
        fprintf(fp, "5 %f moveto\n", ymax - 15);
        fprintf(fp, "(File: %s) show\n", dm.source_filename);
        fprintf(fp, "5 %f moveto\n", ymax - 30);
        fprintf(fp, "(Threshold: %.2f) show\n", threshold);

        fprintf(fp, "gsave\n.75 .75 .75 setrgbcolor\n");
        fprintf(fp, "newpath\n"
                    "%f %f moveto\n"
                    "0 %f rlineto\n"
                    "%f 0 rlineto\n"
                    "0 -%f rlineto\n"
                    "-%f 0 rlineto\n"
                    "closepath\n"
                    "%f setlinewidth\n"
                    "stroke\n",
                    oma + framewidth/2., oma + framewidth/2.,
                    dm.nres * boxw + 1, dm.nres * boxw + 1, dm.nres * boxw + 1,
                    dm.nres * boxw + 1, framewidth);


        fprintf(fp, "grestore\n");
        fprintf(fp, "/square {\n"
                    "/x exch def\n"
                    "/y exch def\n"
                    "gsave\n"
                    "newpath\n"
                    "x y moveto\n"
                    "0 %f rlineto\n"
                    "%f 0 rlineto\n"
                    "0 -%f rlineto\n"
                    "-%f 0 rlineto\n"
                    "closepath\n"
                    "fill\n"
                    "}\n"
                    "def\n", boxw, boxw, boxw, boxw);

        for(i = 0; i < dm.nres; i++){
                for(j = 0; j < dm.nres; j++){
                       if(getdist(dm, i, j) < threshold)
                               fprintf(fp, "%f %f square\n",
                                               xmax - oma - framewidth -  (i+1)*boxw,
                                               oma + framewidth + j * boxw);
                }
        }
        return;
}
