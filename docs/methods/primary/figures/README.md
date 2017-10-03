Methods Figures
===============

This directory is used for static copies of all figures for the methods writeup.

Where possible, the source code used to regenerate the figure is stored as well.
Figure regeneration does not follow the portability constraints described in the
above section for rendering the methods pdf.

Additional installed software requirements to generate figures are:

- indelErrorPGMFigure
  - PGF v2+ : Generally any latex distribtion from 2008 on will have a sufficiently up to date PGF version. Notably, the Centos6 latex installation does not.
  - pdfcrop : Available in package `texlive-extra-utils` on Ubuntu.

