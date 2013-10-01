===================================================================================================
= ORMAN: Optimal Resolution of Ambiguous RNA-Seq Multi-mappings in the Presence of Novel Isoforms = 
=                              Serious README for ORMAN version 1.1                               =
===================================================================================================

* Description *
	ORMAN is a tool for resolving multi-mappings within an RNA-Seq SAM file.

* Installation *
	Extract archive and run "make".  
	You will also need to install IBM CPLEX software and to set CPLEXDIR variable to point to the 
	CPLEX installation directory before invoking make. For example:
		CPLEXDIR=/opt/mycplex make

* Running *
	ORMAN is invoked as following:
		orman -t [thread count] -g [gtf file] -s [sam file] [output file]

	Example:
		orman -t 10 -g super.gtf -s awesome.sam orman_awesome.sam
		Results will be stored in orman_awesome.sam file.

* Parameter explanation *
	--help, -h
		Display not-so-cute help message and exit gracefully.

	--threads, -t [number]
		Set up the maximum number of threads ORMAN may use.
		ORMAN uses multi-threading to speed up CPLEX smoothing step.
		Beware that CPLEX may require a significant amount of memory to complete, 
		and that large number of threads may result in large memory consumption.

		Default value: number of available processors - 1

	--sam, -s [file]
		Specify input SAM file here.
		Requirements for SAM file:
			1. Header must be present (otherwise ORMAN may behave strange ...)
			2. All reads with the same read name have to be clustered together.
			   Or simply: SAM file should be sorted by the read name.
			   If not, ORMAN WILL NOT report any error, but the results 
			   might (and will) be interesting.

	--gtf, -g [file]
		Specify input GTF file here.
		ORMAN supports only GTF files containing ORMAN-specific additional information
		(namely partial_ex field). If you don't have such file, use the provided ormanGTF 
		script to create one. Details will follow later.

* ormanGTF script *
	Use it as follows:
		ormanGTF < [input gtf file] > [new gtf file]

	Example:
		ormanGTF < almost_super.gtf > super.gtf

* Some additional information *
	ORMAN adds YG tag to each read. YG denotes the gene to which the read was assigned by ORMAN.
	If ORMAN couldn't handle read (some weird fusion involved, or it is just crappy read), 
	YG tag will be assigned to "_".

	A friendly advice: ORMAN might use easily >15G memory if the input is quite large SAM file.

* Contact & support *
	Feel free to drop any inquiry to <inumanag at sfu dot ca>.

* Licence *
	Copyright (c) 2012, 2013, Simon Fraser University
	All rights reserved.
	Redistribution and use in source and binary forms, with or without modification,
	are permitted provided that the following conditions are met:
	Redistributions of source code must retain the above copyright notice, this list
	of conditions and the following disclaimer.
	- Redistributions in binary form must reproduce the above copyright notice, this
	  list of conditions and the following disclaimer in the documentation and/or other
	  materials provided with the distribution.
	- Neither the name of the Simon Fraser University nor the names of its contributors may be
	  used to endorse or promote products derived from this software without specific
	  prior written permission.
	THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
	"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
	LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
	A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
	CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
	EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
	PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
	PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
	LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
	NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
	SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

