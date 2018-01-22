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
                fprintf(fp, "# source chain: %c\n", dm.source_chain);
        for(i = 0; i < dm.nres - 1; i++){
                for(j = i + 1; j < dm.nres; j++){
                        if(getdist(dm, i, j) < threshold)
                                fprintf(fp, "%d\t%d\n", i+1, j+1);
                }
        }
        return;
}
