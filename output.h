#ifndef CMAP_OUTPUT_H_
#define CMAP_OUTPUT_H_

#include<stdio.h>
#include "pdb.h"

void write_contacts(FILE *fp, struct distmat dm, double threshold);

#endif // CMAP_OUTPUT_H_
