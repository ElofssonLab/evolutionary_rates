This folder contains:
png figs
converted png figs (.png.tif) with 300 dpi using Imagick: for FNAME in *.png; do convert -units PixelsPerInch $FNAME -density 300 $FNAME.tif; done
tif figs (that are larger than they should be?)
