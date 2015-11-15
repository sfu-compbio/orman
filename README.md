# ORMAN
a.k.a. <i><b>O</b>ptimal <b>R</b>esolution of <b>M</b>ultimapping <b>A</b>mbiguity of R<b>N</b>A-Seq Reads</i>
<hr />

# Introduction

### Eh... ORMAN?

ORMAN is a tool for resolving multi-mappings within an RNA-Seq SAM file. ORMAN was publised in [Bioinformatics in October 2013][3].

### How do I get ORMAN?

Just clone our repository and issue `make` command:

	git clone https://github.com/sfu-compbio/orman.git
	cd orman && make -j

This will work only if you have your CPLEX installed in `/opt/ibm/cplex12_5_1`. If not, keep reading, or just try out our x86-64 [binary](https://github.com/sfu-compbio/orman/tarball/binary) (compiled on CentOS 6.x).

### How do I get CPLEX?

You will need to install IBM CPLEX Optimization Studio (preferably versions >= 12.5.1) before compiling ORMAN. You can get the trial version at [IBM's website](http://www-01.ibm.com/software/websphere/products/optimization/cplex-studio-preview-edition/).

If you are student or a faculty member, you can get CPLEX for free (woohoo!). Just sign up for [Academic Initiative Program](http://www-03.ibm.com/ibm/university/academic/pub/page/academic_initiative).

After you are done with installing your CPLEX, please set the `CPLEXDIR` variable to point to the CPLEX installation directory. For example, if you have installed CPLEX in `~/apps/mycplex`, run `make` as follows:

	CPLEXDIR=~/apps/mycplex make

If you are not using x86-64 system, you will also need to set up `CPLEXARCH` variable. Please refer to the CPLEX documentation to discover its proper value (but please note that we cannot guarantee that ORMAN will run properly on such systems).

### How do I run ORMAN?

ORMAN is invoked as following:

	orman -t [thread count] -g [ORMAN gtf file] -s [sam file] [output file]

Example invocation (results will be stored in `orman_awesome.sam` file):

	orman -t 10 -g super.gtf -s awesome.sam orman_awesome.sam

**Please note that input SAM file should be sorted by the read name, and it should have valid headers.**

#### ormanGTF

ORMAN supports only GTF files containing ORMAN-specific additional information (namely `partial_ex` field). If you don't have such file, use the provided `ormanGTF` script to create compatinge GTF file as follows:
	
	ormanGTF < [input gtf file] > [new gtf file]

Example:

	ormanGTF < almost_super.gtf > super.gtf

# Usage

### Parameter explanation

- `--help, -h`
  
  Display not-so-cute help message and exit gracefully.

- `--threads, -t [number]`

  Set up the maximum number of threads ORMAN may use.  ORMAN uses multi-threading to speed up CPLEX smoothing step.

  Beware that CPLEX may require a significant amount of memory to complete, and that large number of threads may result in large memory consumption.

  Default value: **number of available processors**

- `--sam, -s [file]`

  Specify input SAM file here.
  Requirements for SAM file are:
  
  * Header must be present (otherwise ORMAN may behave strange ...)
  * All reads with the same read name have to be clustered together. 

  Or simply: SAM file should be sorted by the read name.
  If not, ORMAN WILL NOT report any error, but the results might (and will) be interesting.

- `--gtf, -g [file]`

  Specify input GTF file here.
  ORMAN supports only GTF files containing ORMAN-specific additional information (namely `partial_ex` field). If you don't have such file, please use the provided ormanGTF script to create one.

- `--single-end, -1`

  Let ORMAN treat input SAM file as a single-end mapping (i.e. assume that no paired-end mapping was performed).


### Some additional information

ORMAN adds `YG` tag to each read. `YG` denotes the gene to which the read was assigned by ORMAN. If ORMAN couldn't handle read (some weird fusion involved, or it is just crappy read), `YG` tag will be assigned to `_`.

Every resolved read will have its `NH` tag set to 1. Also, for every resolved read, `CC`, `CP` and `HI` tags will be removed.

If ORMAN complains about too many crappy reads, please make sure that GTF and SAM chromosome names march. Common cause is the mismatch caused by `chr` prefix (e.g. `1` vs `chr1`).

> **A friendly advice**: ORMAN might use easily >15G memory if the input is quite large SAM file.

# Support

### Contact & Support

Feel free to drop any inquiry to [inumanag at sfu dot oh canada](mailto:inumanag@).

### Authors

ORMAN has been brought to you by:

- [Dr. Phuong Dao](http://www.cs.sfu.ca/~pdao/personal)
- [Ibrahim Numanagić](http://www.sfu.ca/~inumanag)
- [Yen-Yi Lin](http://www.sfu.ca/~yenyil/)
- [Dr. Faraz Hach](http://www.cs.sfu.ca/~fhach/personal/)
- [Dr. Nilgun Donmez](http://www.cs.toronto.edu/~nild)

from the [Lab for Computational Biology](http://compbio.cs.sfu.ca) at [Simon Fraser University](http://www.sfu.ca) and [Eicher Lab](http://eichlerlab.gs.washington.edu/) at [University of Washington](http://www.washington.edu).

### Licence

Copyright (c) 2012--2013, Simon Fraser University. All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

- Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
- Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
- Neither the name of the Simon Fraser University nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT	LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR	A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR	CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,	EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,	PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR	PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF	LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING	NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS	SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

# Release notes

-	*(15-Nov-2015)* **ORMAN version 1.3 release fix**
	-	Added single-end mode (`-1`) for samples which do not contain paired-end reads.

-	*(12-Jun-2014)* **ORMAN version 1.2 release fix**
	-	Additional headers are added so that ORMAN can be compiled on systems without boost libraries. 
	-	ormanGTF script is now working on unix-systems without any modifications.

-	*(06-Mar-2014)* **ORMAN version 1.2 release**
	-	NH tags are properly updated. 
	-	CC, CP and HI tags are removed if the read is resolved.

-	*(14-Jan-2014)* **ORMAN version 1.1.1 bugfix release**
	-	A serious bug was found in the ILP formulation function, which was fixed in this release. Please upgrade to v1.1.1 as soon as possible.

-	*(24-Oct-2013)* **Compile error fixes**
	-	The previously posted source code archive was missing few header files necessary for the compilation process. The missing files are added back now. Please re-download the source code if this bug affected you.

-	*(30-Sep-2013)* **ORMAN version 1.1 release**
	-	Initial public release

[3]: http://bioinformatics.oxfordjournals.org/content/30/5/644
