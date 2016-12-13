Mathematica Notebook for Analytic Sommerfeld Corrections
===========
This folder contains a Mathematica notebook and corresponding files to calculate the Sommerfeld corrections for the annihilation of colored particles for any partial wave as described in [arXiv:1612.02825](https://arxiv.org/abs/1612.02825). Each of the other files are needed by the Mathematica notebook and should be in the same folder.

The repository contains the following files & folders:
* sommerfeld.nb         - Mathematica notebook with all the calculations.
* PlotUtility.m         - Mathematica package used to make the plots.
* math_store/           - Folder which has intermediate Mathematica results stored.
* main_default.c        - Default main file for micrOMEGAs with placeholders for Mathematica.
* main_micromegas.c     - Main file of micrOMEGAs in which the Mathematica notebook inserts the (Sommerfeld-corrected) annihilation cross sections.
* sommerfeld_default.py - Default python file for Sommerfeld calculations with placeholders for Mathematica.
* sommerfeld.py         - Python script that calculates the (Sommerfeld-corrected) annihilation cross sections.

Version: 1.0

License: <a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons License" style="border-width:0"  height=20 src="http://i.creativecommons.org/l/by/4.0/88x31.png"/></a>
