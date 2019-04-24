# local_assembly_nonreferenceTE
This script extracts mapped reads including mates of the discordant from a window of 500bp from a given breakpoint and performs a local assembly of the reads. 

Usage:

    perl $scriptname -t <table> -f <file with ltr cordinates> -l <location of bamfiles>[-p <path of the outputdirectory>][-te <te sequences>][-v] [-c] [-h] [s] 
	
	
	
    MANDATORY ARGUMENT:	
    -t,--table (STRING) file contain accession information first column needs to be the IDs, second column BAMIDs
    -f,--file  (STRING) file containing accession information 
    -l,--bamloc(STRING)	location of bam files
    -sq,--seqtk(STRING)  path to seqtk ((https://github.com/lh3/seqtk))
    -pc,--picard(STRING) path to picard tools (https://github.com/broadinstitute/picard)
    -cp,--cap3  (STRING) path to cap3 (http://seq.cs.iastate.edu/cap3.html) 
    -bp,--blast (STRING) path to blast	(ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) 
    
    samtools has to be installed (http://www.htslib.org/) and has to be in the path
    	  
    OPTIONAL ARGUMENTS:
    -p,--path   (STRING) output directory name (path)
                         Default = <current working directory>
    -te,--TEfile(STRING) path to transposable element file (performs blast of assembly sequence against TE sequence provided)
    -c,--chlog  (BOOL)   Print updates
    -v,--v      (BOOL)   Print version if only option
    -s,--verbose(BOOL)   The script will talk to you
    -h,--help>  (BOOL)   Print this usage\n\n";
    
     
