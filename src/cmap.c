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

#define _XOPEN_SOURCE_EXTENDED
#include<getopt.h>
#include<locale.h>
#include<math.h>
#include<ncurses.h>
#include<stdbool.h>
#include<string.h>
#include<stdio.h>
#include<stdlib.h>
#include<wchar.h>

#include"cmap.h"
#include"pdb.h"
#include"output.h"

bool g_has_colours = false;

/**
 * make_hpos_str: create a string marking every 20th position in protein chain
 * @x_draw_limit: maximum number of columns to be printed in unicode
 * representation of contact map.
 *
 * allocates and returns a string of length @x_draw_limit+1 which should be freed
 * manually if necessary.
 */
char *
make_hpos_str(unsigned int x_draw_limit)
{
        char *s = NULL;
        char buf[16];
        unsigned int residue_counter;
        unsigned int screen_counter;
        size_t num_len;
        unsigned int i;

        s = malloc((x_draw_limit + 1) * sizeof(*s));
        if(s == NULL) return NULL;
        for (i = 0; i < x_draw_limit + 1; i++){
                s[i] = ' ';
        }
        residue_counter = 1;
        num_len = 1;
        screen_counter = 0;
        while ( screen_counter < x_draw_limit - (unsigned int)num_len){
                snprintf(buf, 16, "%u", residue_counter);
                strncpy(s + screen_counter, buf, num_len); 
                residue_counter += 20;
                screen_counter += 10;
                num_len = (size_t) floor(log10((double) residue_counter)) + 1; 
        }
        s[x_draw_limit] = '\0';
        return s;
}

/**
 * make_vpos_str: create a string marking every 20th position in protein chain
 * @y_draw_limit: maximum number of rows to be printed in the unicode
 * representation of the contact map.
 *
 * allocates and returns a string of length @y_draw_limit+1 which should be
 * freed manually if necessary.
 */
char *
make_vpos_str(unsigned int y_draw_limit)
{
        char *s = NULL;
        char buf[16];
        unsigned int residue_counter;
        unsigned int screen_counter;
        size_t num_len;
        unsigned int i;
        s = malloc((y_draw_limit + 1) * sizeof(*s));
        if(s == NULL) return NULL;
        for (i = 0; i < y_draw_limit + 1; i++){
                s[i] = ' ';
        }
        residue_counter = 1;
        num_len = 1;
        screen_counter = 0;
        while ( screen_counter < y_draw_limit - (unsigned int)num_len){
                snprintf(buf, 16, "%u", residue_counter);
                strncpy(s + screen_counter, buf, num_len); 
                residue_counter += 20;
                screen_counter += 5;
                num_len = (size_t) floor(log10( (double) residue_counter)) + 1; 
        }
        s[y_draw_limit] = '\0';
        return s;
}

/**
 * init_curses: initialise curses environment
 */
void
init_curses(void)
{
        initscr();			
        cbreak();
        if(has_colors() == true){
                g_has_colours = true;
                start_color();
                /* Pair 1: background
                 * Pair 2: status bar
                 * Pairs 3,4: contact map
                 * Pair 5: position counter bars
                 */
                init_pair(1, COLOR_WHITE, COLOR_BLACK);
                init_pair(2, COLOR_YELLOW, COLOR_BLACK);
                init_pair(3, COLOR_WHITE, COLOR_BLUE);
                init_pair(4, COLOR_WHITE, COLOR_BLACK);
                init_pair(5, COLOR_WHITE, COLOR_BLACK);
        }
        if(can_change_color() == true){
                def_colours();
                init_pair(1, COLOR_WHITE, 13);
                init_pair(2, 10, 13);
                init_pair(3, COLOR_WHITE, 11);
                init_pair(4, COLOR_WHITE, 12);
                init_pair(5, COLOR_WHITE, 13);
        }
        keypad(stdscr, true);
        curs_set(0);
        noecho();
        refresh();
}

/**
 * def_colours: define custom colours for use with curses
 */
void
def_colours(void)
{
        init_color(COLOR_BLACK, 0, 0, 0);
        init_color(COLOR_WHITE, 1000, 1000, 1000);
        init_color(13, 250, 250, 250); // a medium grey
        /* From flatuicolors.com: */
        init_color(10, 945, 769, 59);  // "sun flower"
        init_color(11, 204, 286, 369); // "wet asphalt"
        init_color(12, 173, 243, 314); // "midnight blue"
}

/**
 * draw_bg: fill a window with background colour
 * @scr: pointer to WINDOW to draw on
 * @rows: rows to draw
 * @rows: columns to draw
 *
 * Draws a box in WINDOW, with one corner at (0,0) and the other at
 * (rows-1, cols-1), with the background colour of COLOR_PAIR(1), or the
 * default background color if global variable g_has_colours is false.
 */
void
draw_bg(WINDOW * scr, unsigned int cols, unsigned int rows){
        unsigned int i, j;
        if(g_has_colours) attron(COLOR_PAIR(1));
        for(j = 0; j < cols; j++){
                for(i = 0; i < rows; i++){
                        mvwaddch(scr, j, i, ' '); 
                }
        }
}

/**
 * draw_contacts_pad
 *
 * @dist
 * @threshold
 */
WINDOW *
draw_contacts_pad(struct distmat dist, double threshold)
{

        int i, j;
        int prot_y, prot_x;
        wchar_t ch = ' ';
        wchar_t s[2] = L" ";
        WINDOW *contacts = NULL;
        int nres;

        nres = dist.nres;
        int x_draw_limit = nres % 2 == 0 ? nres/2 : nres/2 + 1;
        int y_draw_limit = nres % 4 == 0 ? nres/4 : nres/4 + 1;
        contacts = newpad(y_draw_limit, x_draw_limit);
        if(contacts == NULL){
                return NULL;
        }

        for(j = 0, prot_y = 0; j < nres/4; j++, prot_y += 4){
                for (i=0, prot_x=0; i < nres/2; i++, prot_x+=2)
                {
                        if(g_has_colours){
                                if(CHECKB_LIGHT(j,i))
                                        wattron(contacts, COLOR_PAIR(3));
                                else wattron(contacts, COLOR_PAIR(4));
                        }
                        ch = 0x2800;
                        if(getdist(dist, prot_y, prot_x) <= threshold)
                                ch+=0x01;
                        if(getdist(dist, prot_y + 1, prot_x) <= threshold)
                                ch+=0x02;
                        if(getdist(dist, prot_y + 2, prot_x) <= threshold)
                                ch+=0x04;
                        if(getdist(dist, prot_y + 3, prot_x) <= threshold)
                                ch+=0x40;
                        if(getdist(dist, prot_y, prot_x + 1) <= threshold)
                                ch+=0x08;
                        if(getdist(dist, prot_y + 1, prot_x + 1) <= threshold)
                                ch+=0x10;
                        if(getdist(dist, prot_y + 2, prot_x + 1) <= threshold)
                                ch+=0x20;
                        if(getdist(dist, prot_y + 3, prot_x + 1) <= threshold)
                                ch+=0x80;
                        s[0] = ch; 
                        waddwstr(contacts, s);
                }
                /* Special case for printing final column */
                if(i < x_draw_limit){
                        if(g_has_colours){
                                if(CHECKB_LIGHT(j,i))
                                        wattron(contacts, COLOR_PAIR(3));
                                else wattron(contacts, COLOR_PAIR(4));
                        }

                        ch = 0x2800;
                        if( getdist(dist, prot_y, prot_x) <= threshold) 
                                ch+=0x01;
                        if( getdist(dist, prot_y + 1, prot_x) <= threshold)
                                ch+=0x02;
                        if( getdist(dist, prot_y + 2, prot_x) <= threshold)
                                ch+=0x04;
                        if( getdist(dist, prot_y + 3, prot_x) <= threshold)
                                ch+=0x40;
                        s[0] = ch; 
                        waddwstr(contacts, s);
                }
                mvwaddnstr(contacts, j + 1, 0, "", 0);
        }
        /* Special case for last row */
        if(j < y_draw_limit){
                for (i = 0, prot_x = 0; i < nres / 2; i++, prot_x += 2)
                {
                        if(g_has_colours){
                                if(CHECKB_LIGHT(j, i))
                                        wattron(contacts, COLOR_PAIR(3));
                                else wattron(contacts, COLOR_PAIR(4));
                        }

                        ch = 0x2800;
                        /* falling through on purpose here */
                        switch(nres - (4 * j)){
                        case 3:
                                if( getdist(dist, prot_y + 2, prot_x) <= threshold)
                                        ch+=0x04;
                                if( getdist(dist, prot_y + 2, prot_x + 1) <= threshold)
                                        ch+=0x20;
                        case 2: 
                                if( getdist(dist, prot_y + 1, prot_x) <= threshold)
                                        ch+=0x02;
                                if(getdist(dist, prot_y + 1, prot_x + 1) <= threshold)
                                        ch+=0x10;
                        default:
                        case 1:
                                if( getdist(dist, prot_y, prot_x)
                                                <= threshold)
                                        ch+=0x01;
                                if(getdist(dist, prot_y, prot_x + 1) <= threshold)
                                        ch+=0x08;
                                break;
                        }
                        s[0] = ch; 
                        waddwstr(contacts, s);
                }
                /* Special special case for printing final LR corner */
                if(i < x_draw_limit){
                        if(g_has_colours){
                                if(CHECKB_LIGHT(j, i))
                                        wattron(contacts, COLOR_PAIR(3));
                                else wattron(contacts, COLOR_PAIR(4));
                        }

                        ch = 0x2800;
                        switch(nres - (4 * j)){
                        case 3:
                                if( getdist(dist, prot_y + 2, prot_x)
                                                <= threshold)
                                        ch+=0x04;
                        case 2: 
                                if( getdist(dist, prot_y + 1, prot_x)
                                                <= threshold)
                                        ch+=0x02;
                        case 1:
                                if( getdist(dist, prot_y, prot_x)
                                                <= threshold)
                                        ch+=0x01;
                        }
                        s[0] = ch; 
                        waddwstr(contacts, s);
                }
        }
        return contacts;
}

WINDOW *
draw_status_pad(char *filename, char chain, int nres, double threshold)
{
        int i;
        WINDOW *status = NULL;
        status = newpad(1, 1025);
        if(status == NULL){
                return NULL;
        }
        if(g_has_colours) wattron(status, COLOR_PAIR(2));
        for(i =0; i<1024; i++)
                wprintw(status, "\u2501");
        mvwaddstr( status, 0, 0, "\u2501");
        wattron(status, A_REVERSE);
        wprintw(status, " %s ", filename);
        wattroff(status, A_REVERSE);
        wprintw(status, "\u2501");
        wattron(status, A_REVERSE);
        wprintw(status, " Chain: %c ", chain);
        wattroff(status, A_REVERSE);
        wprintw(status, "\u2501");
        wattron(status, A_REVERSE);
        wprintw(status, " Residues: %d ", nres);
        wattroff(status, A_REVERSE);
        wprintw(status, "\u2501");
        wattron(status, A_REVERSE);
        wprintw(status, " Threshold: %2.2f ", threshold);
        wattroff(status, A_REVERSE);
        mvwaddstr(status, 0, 1024,  "");
        return status;
}

int
main(int argc, char **argv)
{	
        setlocale(LC_CTYPE, "");
        int j;
        char *filename = NULL;
        char *ofname = NULL;
        char *epsname = NULL;
        FILE *ofp;
        struct distmat * dist = NULL; 
        double threshold = 8;
        char chain = 'A';
        int nres = 0; 
        struct coords *cs = NULL;
        int nrow, ncol;
        char usage_str[1024];

        snprintf(usage_str, 1024, "cmap version 0.1\n\n"
                        "Usage:\n"
                        "  cmap [-c chain] [-e eps_file] [-o output_file] [-t threshold] <pdb_file>\n"
                        "\nOptions:\n"
                        "  -c, --chain      chain to read from input file\n"
                        "  -e, --eps        save EPS image of contact map (experimental)\n"
                        "  -h, --help       show this message\n"
                        "  -o, --output     save contacts to file\n"
                        "  -t, --threshold  distance threshold for contact (Angstroms)\n"
                        "\n");

        /*
         * Parse command line arguments
         */
        static struct option long_options[] =
        {
                {"chain", required_argument, 0, 'c'},
                {"eps", required_argument, 0, 'e'},
                {"help", no_argument, 0, 'h'},
                {"output", required_argument, 0, 'o'},
                {"threshold", required_argument, 0, 't'},
                {0, 0, 0, 0}
        };

        int option_index = 0;
        int opt;
        while(1){
                opt = getopt_long(argc, argv, "c:e:ho:t:", long_options, &option_index);
                if(opt == -1)
                        break;
                if (opt == 'c'){
                        chain = optarg[0];
                }
                if (opt == 't'){
                        threshold = atof(optarg);
                }
                if(opt == 'h'){
                        printf("%s", usage_str);
                        return 0;
                }
                if(opt == 'o'){
                        ofname = optarg;
                }
                if(opt == 'e'){
                        epsname = optarg;
                }
        }
        if( argc < 2){
                fprintf(stderr, "%s", usage_str);
                fprintf(stderr, "FATAL: Must specify a filename.\n");
                return 1;
        }
        else {
                filename = argv[optind];
                if(filename == NULL){
                        fprintf(stderr, "%s", usage_str);
                        fprintf(stderr, "FATAL: Must specify a filename.\n");
                        return 1;
                }
        }


        /*
         * Read PDB coords and calculate distances
         */
        cs = getcoords( filename, chain);
        if (cs == NULL){
                fprintf(stderr, "FATAL: couldn't read coordinates from file [%s].\nTried to read chain [%c].\n", filename, chain);
                return 1;
        }
        nres = cs->nres;
        dist = calculate_distmat(*cs);
        if (dist == NULL){
                fprintf(stderr, "FATAL: couldn't allocate memory.\n");
                return 1;
        }
        freecoords(cs);
        cs = NULL;

        /*
         * Write contacts to file (optional)
         */

        if(ofname != NULL){
                ofp = fopen(ofname, "w");
                if(ofp == NULL){
                        fprintf(stderr, "%s", usage_str);
                        fprintf(stderr, "FATAL: couldn't open output file [%s]\n", ofname);
                        return 1;
                }
                write_contacts(ofp, *dist, threshold);
                printf("Wrote contacts to file [%s].\n", ofname);
                fclose(ofp);
                ofp = NULL;
        }

        /*
         * Write EPS file (optional)
         */

        if(epsname != NULL){
                ofp = fopen(epsname, "w");
                if(ofp == NULL){
                        fprintf(stderr, "%s", usage_str);
                        fprintf(stderr, "FATAL: couldn't open output file [%s]\n", epsname);
                        return 1;
                }
                write_eps(ofp, *dist, threshold);
                printf("Wrote postscript to file [%s].\n", epsname);
                fclose(ofp);
                ofp = NULL;
        }

        /* 
         * Set up curses display
         */

        init_curses();
        getmaxyx(stdscr, nrow, ncol);

        draw_bg(stdscr, nrow, ncol);
        wnoutrefresh(stdscr);

        /*
         * Store how many rows(y) and columns(x) of screen space are needed to
         * display the contact map. One screen character displays 4 rows and 
         * two columns of the contact map.
         */
        int x_draw_limit = nres % 2 == 0 ? nres/2 : nres/2 + 1;
        int y_draw_limit = nres % 4 == 0 ? nres/4 : nres/4 + 1;

        /* 
         * Draw status bar
         */
        WINDOW *status = NULL;
        status = draw_status_pad(filename, chain, nres, threshold);
        if(status == NULL){
                endwin();
                fprintf(stderr, "FATAL: error drawing curses display.");
                return 1;
        }
        pnoutrefresh(status, 0, 0, 0, 0, 1, ncol-1);


        /*
         * Draw horizontal and vertical position bars
         */

        WINDOW *hpos = NULL;
        hpos = newpad(1, x_draw_limit);
        if(hpos == NULL){
                fprintf(stderr, "FATAL: error drawing curses display.\n");
                return 1;
        } 
        char *s1 = NULL;
        s1 = make_hpos_str(x_draw_limit);
        if (s1 == NULL){
                fprintf(stderr, "FATAL: couldn't allocate memory.\n");
                return 1;
        }
        if(g_has_colours) wattron(hpos, COLOR_PAIR(5));
        waddnstr(hpos, s1, x_draw_limit);
        free(s1);
        s1 = NULL;
        
        WINDOW *vpos = NULL;
        vpos = newpad(y_draw_limit, 1);
        if(vpos == NULL){
                fprintf(stderr, "FATAL: error drawing curses display.\n");
                return 1;
        } 
        s1 = make_vpos_str(y_draw_limit);
        if (s1 == NULL){
                fprintf(stderr, "FATAL: couldn't allocate memory.\n");
                return 1;
        }
        if(g_has_colours) wattron(vpos, COLOR_PAIR(5));
        for(j = 0; j < y_draw_limit; j++)
            waddch(vpos, s1[j]);
        free(s1);


        WINDOW *contacts = NULL;
        contacts = draw_contacts_pad(*dist, threshold);        
        if(contacts == NULL){
                fprintf(stderr, "FATAL: error drawing curses display.\n");
                return 1;
        } 
        pnoutrefresh(hpos, 0, 0, 1, 1, 1, ncol - 1);
        pnoutrefresh(vpos, 0, 0, 2, 0, nrow - 1, 1);
        pnoutrefresh(contacts, 0, 0, 2, 1, nrow - 1, ncol - 1);
        doupdate();

        /*
         * Loop for keyboard input
         */
        int x_offset = 0;
        int y_offset = 0;
        int c;
        bool pressed_g = false; /* */
        while(1){
                c =  getch();		
                /* q or Q to quit */
                if (c == 'q' || c == 'Q')
                        break;
                if (c != 'g') 
                        pressed_g = false;
                switch(c){
                        /* Move view with WASD, hjkl, arrow keys*/ 
                        case KEY_LEFT:
                        case 'A':
                        case 'a':
                        case 'h':
                                x_offset-=10;
                                if (x_offset < 0)
                                        x_offset =0;
                                break;
                        case KEY_RIGHT:
                        case 'D':
                        case 'd':
                        case 'l':
                                x_offset+=10;
                                if (x_offset + ncol > x_draw_limit)
                                        x_offset = x_draw_limit - ncol + 1;
                                if (x_draw_limit < ncol)
                                        x_offset = 0;
                                break;
                        case KEY_UP:
                        case 'W':
                        case 'w':
                        case 'k':
                                y_offset -= 5;
                                if (y_offset < 0)
                                        y_offset = 0;
                                break;
                        case KEY_DOWN:
                        case 'S':
                        case 's':
                        case 'j':
                                y_offset += 5;
                                if (y_offset + nrow > y_draw_limit)
                                        y_offset = y_draw_limit - nrow + 2;
                                if (y_draw_limit < nrow)
                                        y_offset = 0;
                                break;
                                /*
                                 * Vim style movement shortcuts with 0,^,$,G,gg
                                 * $: move to end/right
                                 * ^,0: move to beginning/left
                                 * G: move to bottom
                                 * gg: move to top
                                 */
                        case '$': 
                                x_offset = x_draw_limit - ncol + 1;
                                if (x_draw_limit < ncol)
                                        x_offset = 0;
                                break;
                        case '^':
                        case '0':
                                x_offset = 0;
                                break;
                        case 'G':
                                y_offset = y_draw_limit - nrow + 2;
                                break;
                        case 'g':
                                if (pressed_g){
                                        y_offset = 0;
                                        pressed_g = false;
                                } else {
                                        pressed_g = true;
                                }
                                break;
                        /* Change distance threshold with +/- */
                        case '+':
                                threshold += 0.5;
                                delwin(contacts);
                                contacts = NULL;
                                contacts = draw_contacts_pad(*dist, threshold);        
                                if(contacts == NULL){
                                        fprintf(stderr, "FATAL: error drawing curses display.\n");
                                        return 1;
                                } 
                                status = NULL;
                                status = draw_status_pad(filename, chain, nres, threshold);
                                if(status == NULL){
                                        endwin();
                                        fprintf(stderr, "FATAL: error drawing curses display.");
                                        return 1;
                                }
                                pnoutrefresh(status, 0, 0, 0, 0, 1, ncol-1);

                                break;
                        case '-':
                                threshold -= 0.5;
                                if (threshold < 0) threshold = 0;
                                delwin(contacts);
                                contacts = NULL;
                                contacts = draw_contacts_pad(*dist, threshold);        
                                if(contacts == NULL){
                                        fprintf(stderr, "FATAL: error drawing curses display.\n");
                                        return 1;
                                } 
                                status = NULL;
                                status = draw_status_pad(filename, chain, nres, threshold);
                                if(status == NULL){
                                        endwin();
                                        fprintf(stderr, "FATAL: error drawing curses display.");
                                        return 1;
                                }
                                pnoutrefresh(status, 0, 0, 0, 0, 1, ncol-1);

                                break;
                                /* Handle terminal resizing */
                        case KEY_RESIZE:
                                getmaxyx(stdscr, nrow, ncol);
                                draw_bg(stdscr, nrow, ncol);
                                break;
                }
                wnoutrefresh(stdscr);
                pnoutrefresh(hpos, 0, 0 + x_offset, 1, 1, 1, ncol - 1);
                pnoutrefresh(vpos, 0 + y_offset, 0, 2, 0, nrow - 1, 1);
                pnoutrefresh(contacts, 0 + y_offset, 0 + x_offset,
                                2, 1, nrow - 1, ncol - 1);
                pnoutrefresh(status, 0, 0, 0, 0, 1, ncol - 1);
                doupdate();
        }

        /*
         * Clean up and exit
         */
        freedm(dist);
        dist = NULL;

        delwin(status);
        delwin(contacts);
        delwin(hpos);
        delwin(vpos);
        endwin();	

        return 0;
}

