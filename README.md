# cmap - quickly preview protein contact maps in a terminal

## What is it?
cmap is a command line program which reads [Protein Data Bank](http://pdb101.rcsb.org/learn/guide-to-understanding-pdb-data/introduction) coordinate files and uses ncurses to display a [protein contact map](https://en.wikipedia.org/wiki/Protein_contact_map) in a terminal.

![Screenshot](screenshots/screenshot1.png?raw=true)

## Installation

    git clone https://github.com/gearoidfox/cmap.git
    cd cmap
    ./configure
    make
    make install

### Dependencies

cmap requires an ncurses library with wide character support.

## Compatibility

cmap uses unicode characters to draw contact maps. If the output doesn't display correctly, ensure your terminal emulator is using a compatible encoding like UTF-8. 

### Colours

cmap is designed for 256 colour terminals that support redefining colours. If possible, set the TERM environment variable to xterm-256color:

    export TERM=xterm-256color

cmap will also run in terminals with more limited colour support--e.g. xterm or vt220 modes--but it won't look as nice. 

### Fonts

Because it uses unicode characters for drawing, cmap is sensitive to the choice of font used by the terminal emulator. [Ubuntu Mono Regular](https://design.ubuntu.com/font/) is the best looking font I tested.

## Keyboard shortcuts

- q to quit.
- [w, a, s, d], [h, j, k, l], and the arrow keys move around the contact map view.
- vim-style shortcuts gg, G, ^, and $ move to the top, bottom, left edge, and right edge, respectively.
- [+,-] can be used to increase and decrease the distance threshold.

### Screenshots

![Screenshot](screenshots/screenshot2.png?raw=true)
