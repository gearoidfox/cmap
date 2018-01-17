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

/**
 * CHECKB_LIGHT(y, x):
 * Given screen coordinates (y,x), this macro returns TRUE if we should colour
 * the background at that point with a lighter colour--or FALSE if it should 
 * be coloured dark--to achieve a checkerboard effect.
 * The checkboard is 10 columns * 5 rows in screen characters, which
 * translates to 20 * 20 residues in the contact map.
 */

#define CHECKB_LIGHT(y, x) (((  ((x) % 20  < 10) \
                              &&((y) % 10  <  5) \
                             )                   \
                           ||(  ((x) % 20 >= 10) \
                              &&((y) % 10 >= 5)  \
                             )) ? TRUE : FALSE) 
struct distmat;

void def_colours(void);
void draw_bg(WINDOW * scr, unsigned int cols, unsigned int rows);
WINDOW * draw_contacts_pad(struct distmat dist, double threshold);
WINDOW * draw_status_pad(char *filename, char chain, int nres, double threshold);
void init_curses(void);
int main(int argc, char **argv);
char * make_hpos_str(unsigned int x_draw_limit);
char * make_vpos_str(unsigned int y_draw_limit);
