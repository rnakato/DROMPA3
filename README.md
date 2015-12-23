# DROMPA3 README (version 3.1.0)

#0. Recent news
* 2015-12-24
3.1.0 
improve pdf creating and change the tool for marging pdf files from pdftk to cpdf (Coherent PDF) because pdftk is no longer supported on several OS (e.g., CentOS 7).

* 2015-06-10  
3.0.0 (Beta) first commit.

#1. Overview
DROMPA (DRaw and Observe Multiple enrichment Profiles and Annotation) is a program for user-friendly and flexible ChIP-seq pipelining. DROMPA can be used for quality check, PCRbias filtering, normalization, peak calling, visualization and other multiple analyses of ChIP-seq data. DROMPA is specially designed so that it is easy to handle, and for users without a strong bioinformatics background.

#2. Install and Usage
To install DROMPA, simply type "make".

DROMPA requires the following programs and libraries:
* GCC compiler (http://gcc.gnu.org/)
* Cairo libraries (http://www.cairographics.org/)
* GTK library (http://www.gtk.org/)
* GNU Scientific Library (http://www.gnu.org/software/gsl/)
* Coherent PDF (http://community.coherentpdf.com/)
* (optional) SAMtools (http://samtools.sourceforge.net/)
* (optional) R (http://www.r-project.org/)

 See Manual.pdf for detail.


#3. Reference
Nakato, R., Itoh T. and Shirahige K.: DROMPA: easy-to-handle peak calling and visualization software for the computational analysis and validation of ChIP-seq data, Genes to Cells, vol.18, issue 7, 2013.

Please direct bug reports and questions about usage to rnakato@iam.u-tokyo.ac.jp.
