# DROMPA3 README

Please direct bug reports and questions about usage to rnakato@iam.u-tokyo.ac.jp.

# 1. Overview
DROMPA (DRaw and Observe Multiple enrichment Profiles and Annotation) is a program for user-friendly and flexible ChIP-seq pipelining. DROMPA can be used for quality check, PCRbias filtering, normalization, peak calling, visualization and other multiple analyses of ChIP-seq data. DROMPA is specially designed so that it is easy to handle, and for users without a strong bioinformatics background.

# 2. Install
DROMPA is written in C and requires the following programs and libraries:
* [Cairo libraries](http://www.cairographics.org/)
* [GTK library](http://www.gtk.org/)
* [GNU Scientific Library](http://www.gnu.org/software/gsl/)
* [zlib](http://www.zlib.net/)
* [SAMtools](http://samtools.sourceforge.net/) (for BAM formatted input)
* [R](http://www.r-project.org/) (for PROFILE command)

#### 2.1. Install required libraries
for Ubuntu:

     sudo apt-get install git gcc libgtk2.0-dev libgsl-dev samtools r-base
 
for CentOS:

     sudo yum -y install zlib-devel gsl-devel gtk2-devel

#### 2.2. Install cpdf
 DROMPA uses Coherent PDF (http://community.coherentpdf.com/) for merging pdf files.
 
     git clone https://github.com/coherentgraphics/cpdf-binaries.git
    
#### 2.3. Install DROMPA
    git clone https://github.com/rnakato/DROMPA3.git
    cd DROMPA3
    make

If you get an installation error, make sure that all required libraries are installed.

#### 2.4. Add the PATH environment variable
For example, if you downloaded DROMPA and cpdf into the $HOME/my_chipseq_exp directory, type:

    export PATH = $PATH:$HOME/my_chipseq_exp/DROMPA3:$HOME/my_chipseq_exp/cpdf-binaries/Linux-Intel-64bit/


# 3. Usage
 See Manual.pdf for detail.

# 4. Reference
1. Nakato R., Shirahige K. Recent advances in ChIP-seq analysis: from quality management to whole-genome annotation, Briefings in Bioinformatics, 2016.

2. Nakato, R., Itoh T. and Shirahige K.: DROMPA: easy-to-handle peak calling and visualization software for the computational analysis and validation of ChIP-seq data, Genes to Cells, vol.18, issue 7, 2013.
