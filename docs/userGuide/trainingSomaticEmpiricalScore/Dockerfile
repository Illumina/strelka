#
# minimum environment required to run Strelka somatic EVS training procedures
#

FROM ubuntu:14.04

#
# major required packages:
#
# - python
# - numpy 
# - pandas
# - pip
# - scikit-learn
# - bx-python 
#
RUN apt-get update -qq && \
    apt-get install -qq python python-numpy python-pandas && \
    apt-get install -qq python-pip python-dev build-essential && \ 
    pip install -U scikit-learn && \
    apt-get install -qq zlib1g-dev && \ 
    pip install bx-python

#
# Optional packages for visualization discussed
# in user guide but not required for training:
#
# - R
# - R-ggplot2
#
RUN apt-get install -qq r-base r-cran-ggplot2

