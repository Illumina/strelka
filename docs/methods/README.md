Strelka Methods
==============

This directory contains documents describing Strelka's primary methods and any major experimental approaches. Methods documentation is maintained in latex.

Each directory should contain a `makepdf.bash` script, which produces a `methods.pdf` file as output (so constrained to facilite automated testing of the latex source build). The build for each document should only require the contents of the standard `texlive` package (as well as `texlive-latex` on RHEL/centos). Any additional latex requirements should be included in-source.

