
****** DDE-BIFTOOL v. 3.1.1 ******
    * Installation
    * Reference_and_documentation
    * Contributors
    * Citation
    * Copyright,_License_and_No-warranty_notice
===============================================================================
***** Installation *****
    * Unzipping ddebiftool.zip creates a "dde_biftool" directory (named
      "dde_biftool") containing the subfolders:
          o ddebiftool (basic DDE-BIFTOOL routines),
          o demos (example scripts illustrating the use of DDE-BIFTOOL),
          o ddebiftool_extra_psol (extension for local bifurcations of periodic
            orbits),
          o ddebiftool_extra_nmfm (extension for normal form coefficient
            computations at local bifurcations of equilibria in DDEs with
            constant delay),
          o ddebiftool_utilities (auxiliary functions),
          o ddebiftool_extra_rotsym (extension for systems with rotational
            symmetry).
          o external_tools (support scripts, such as a Mathematica and a Maple
            script for generating derivative functions used in DDE-BIFTOOL).
    * To test the tutorial demo "neuron" (the instructions below assume
      familiarity with Matlab or octave):
          o Start Matlab (version 7.0 or higher) or octave (tested with version
            3.8.1)
          o Inside Matlab or octave change working directory to demos/neuron
            using the "cd" command
          o Execute script "rundemo" to perform all steps of the tutorial demo
          o Compare the outputs on screen and in figure windows with the
            published output in demos/neuron/html/demo1.html.
===============================================================================
***** Reference and documentation *****
  Current download URL on Sourceforge (including access to versions from 3.1
  onward)
      https://sourceforge.net/projects/ddebiftool/
  URL of original DDE-BIFTOOL website (including access to versions up to 3.0)
      http://twr.cs.kuleuven.be/research/software/delay/ddebiftool.shtml
  Contact (bug reports, questions etc)
      https://sourceforge.net/projects/ddebiftool/support
  Manual for version 2.0x
      K. Engelborghs, T. Luzyanina, G. Samaey. DDE-BIFTOOL v. 2.00: a Matlab
      package for bifurcation analysis of delay differential equations.
      Technical Report TW-330
  Manual for current version
      manual.pdf (v. 3.1.1), stored on arxiv: arxiv.org/abs/1406.7144
  Changes for v. 2.03
      Addendum_Manual_DDE-BIFTOOL_2_03.pdf (by K. Verheyden)
  Changes for v. 3.0
      Changes-v3.pdf (by J. Sieber)
  Description of extensions ddebiftool_extra_psol and ddebiftool_extra_rotsym
      Extra_psol_extension.pdf (by J. Sieber)
  Description of the extention ddebiftool_extra_nmfm
      nmfm_extension_desctiption.pdf (by M. Bosschaert, B. Wage,  Yu.A.
      Kuznetsov)
  Overview of documented demos
      demos/index.html
===============================================================================
***** Contributors *****
** Original code and documentation (v. 2.03) **
K. Engelborghs, T. Luzyanina, G. Samaey. D. Roose, K. Verheyden
K.U.Leuven
Department of Computer Science
Celestijnenlaan 200A
B-3001 Leuven
Belgium
** Revision for v. 3.0, 3.1.x **
** Bifurcations of periodic orbits **
J. Sieber
College for Engineering, Mathematics and Physical Sciences, University of
Exeter (UK),
emps.exeter.ac.uk/mathematics/staff/js543
** Normal form coefficients for bifurcations of equilibria **
S. Janssens, B. Wage, M. Bosschaert, Yu.A. Kuznetsov
Utrecht University
Department of Mathematics
Budapestlaan 6
3584 CD Utrecht
The Netherlands
www.staff.science.uu.nl/~kouzn101/_((Y.A._Kuznetsov)
** Automatic generation of right-hand sides and derivatives in Mathematica **
D. Pieroux
Universite Libre de Bruxelles (ULB, Belgium)
** Demo for phase oscillator **
A. Yeldesbay
Potsdam University (Germany)
===============================================================================
***** Citation *****
Scientific publications, for which the package DDE-BIFTOOL has been used, shall
mention usage of the package DDE-BIFTOOL, and shall cite the following
publications to ensure proper attribution and reproducibility:
    * K. Engelborghs, T. Luzyanina, and D. Roose, Numerical bifurcation
      analysis of delay differential equations using DDE-BIFTOOL, ACM Trans.
      Math. Softw. 28 (1), pp. 1-21, 2002.
    * K. Engelborghs, T. Luzyanina, G. Samaey. DDE-BIFTOOL v. 2.00: a Matlab
      package for bifurcation analysis of delay differential equations.
      Technical Report TW-330, Department of Computer Science, K.U.Leuven,
      Leuven, Belgium, 2001.
    * [Manual of current version, permanent link]
      J. Sieber, K. Engelborghs, T. Luzyanina, G. Samaey, D. Roose: DDE-BIFTOOL
      Manual - Bifurcation analysis of delay differential equations. arxiv.org/
      abs/1406.7144.
    * [Theoretical background for computation of normal form coefficients,
      permanent link]
      Sebastiaan Janssens: On a Normalization Technique for Codimension Two
      Bifurcations of Equilibria of Delay Differential Equations. Master
      Thesis, Utrecht University (NL), supervised by Yu.A. Kuznetsov and O.
      Diekmann, dspace.library.uu.nl/handle/1874/312252, 2010.
    * [Normal form implementation for Hopf-related cases, permanent link]
      Bram Wage: Normal form computations for Delay Differential Equations in
      DDE-BIFTOOL. Master Thesis, Utrecht University (NL), supervised by Y.A.
      Kuznetsov, dspace.library.uu.nl/handle/1874/296912, 2014.
      M. M. Bosschaert: Switching from codimension 2 bifurcations of equilibria
      in delay differential equations. Master Thesis, Utrecht University (NL),
      supervised by Y.A. Kuznetsov, dspace.library.uu.nl/handle/1874/334792,
      2016.
===============================================================================
***** Copyright, License and No-warranty Notice *****
BSD 2-Clause license

Copyright (c) 2017, K.U. Leuven, Department of Computer Science, K.
Engelborghs, T. Luzyanina, G. Samaey. D. Roose, K. Verheyden, J. Sieber, B.
Wage, D. Pieroux

All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation and/
or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
