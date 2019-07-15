# local_assembly_nonreferenceTE
This script extracts mapped reads including mates of the discordant from a window of 500bp from a given breakpoint and performs a local assembly of the reads. 

Additional options: If TE sequences and path to blast are provided,the mates of the discordant reads and assembled contigs are blasted against the provided TE sequence.

Usage:

    perl $scriptname -t <table> -f <file with TE breakpoints> -l <location of bamfiles> -sq <path> -pc <path> -cp <cap3> [-p <path of the outputdirectory>][-te <te sequences> -bp <blast>][-v] [-c] [-h] [s] 
	
	
	
    MANDATORY ARGUMENT:	
    -t,--table (STRING) file containing information of the bam files separated by tab (sampleIDs, name of the bamfile e.g (8463 8463.bam) 
    -f,--file  (STRING) file containing the break point information separted by tab (sampleID,chr,genomic location,expected TE(optional))  (e.g.8463    4       140408656 Alu) 
    -l,--bamloc(STRING) path to bam files
    -sq,--seqtk(STRING)  path to seqtk (seqtk can be downloaded from: (https://github.com/lh3/seqtk))
    -pc,--picard(STRING) path to picard tools (picard tools can be downloaded from: https://github.com/broadinstitute/picard))
    -cp,--cap3  (STRING) path to cap3 (cap3 can be downloaded from: (http://seq.cs.iastate.edu/cap3.html)) 

    
    samtools has to be installed (http://www.htslib.org/) and has to be in the path
    	  
    OPTIONAL ARGUMENTS:
    -p,--path   (STRING) output directory name (path)
                         Default = <current working directory>
    -te,--TEfile(STRING) path to transposable element file (performs blast of assembly sequence against TE sequence provided)
    -bp,--blast (STRING) path to blast	(Blast can be downloaded from: (ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)) needed with -te option 
    -c,--chlog  (BOOL)   Print updates
    -v,--v      (BOOL)   Print version if only option
    -s,--verbose(BOOL)   The script will talk to you
    -h,--help>  (BOOL)   Print this usage\n\n";
    
     
