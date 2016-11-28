# DROMPA3 README (version 3.2.3)

#1. Overview
DROMPA (DRaw and Observe Multiple enrichment Profiles and Annotation) is a program for user-friendly and flexible ChIP-seq pipelining. DROMPA can be used for quality check, PCRbias filtering, normalization, peak calling, visualization and other multiple analyses of ChIP-seq data. DROMPA is specially designed so that it is easy to handle, and for users without a strong bioinformatics background.

#2. Install
DROMPA requires the following programs and libraries:
* Cairo libraries (http://www.cairographics.org/)
* GTK library (http://www.gtk.org/)
* GNU Scientific Library (http://www.gnu.org/software/gsl/)
* (optional) SAMtools (http://samtools.sourceforge.net/)
* (optional) R (http://www.r-project.org/)

#### 2.1. Install required libraries
for Ubuntu:

     sudo apt-get install git gcc libgtk2.0-dev libgsl-dev samtools r-base
 
for CentOS:

     sudo yum -y install zlib-devel gsl-devel gtk2-devel

#### 2.2. Install cpdf
 DROMPA uses Coherent PDF (http://community.coherentpdf.com/) for merging pdf files.
    wget http://github.com/coherentgraphics/cpdf-binaries/archive/master.zip
    unzip master.zip
    
#### 2.3. Install DROMPA
    git clone https://github.com/rnakato/DROMPA3.git
    cd DROMPA3
    make

#3. Usage
 See Manual.pdf for detail.

#4. Reference
Nakato, R., Itoh T. and Shirahige K.: DROMPA: easy-to-handle peak calling and visualization software for the computational analysis and validation of ChIP-seq data, Genes to Cells, vol.18, issue 7, 2013.

Please direct bug reports and questions about usage to rnakato@iam.u-tokyo.ac.jp.
