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
 
 
 
Total run time of the above script depends on the number bams, as the number of bams to analyse increases the run also increases. If you are dealing multiple bam files please follow the steps described below to speed up the process.

Step 1:
 The  --file (the file containing coordinates and individual)  is split to multiple files using the script called 'splitfile_for_parallel_individuals.pl'. Please type 'perl scriptname -h' to see the usage of the script.
A folder called 'splitbyindividuals' is created containing the splitfiles based on the number of individuals requested for split.

	perl splitfile_for_parallel_individuals.pl -f <file that needs to split i.e file with TE breakpoints> -s yes -n <number in individuals>

Step 2:
Then do 'cd splitbyindividuals' and then launch multiple jobs using 'parallel' (https://www.gnu.org/software/parallel/parallel_tutorial.html) .(see example commandline below). Max no of jobs submitted for the parellel (-j) depends on the total number of individuals/number of individuals given with -n option in step 1.

	 cat ../listof_files.txt | nohup /usr/bin/parallel -j X --results path_of_directory_forstderr 'perl $scriptname -t <table> -f <file with TE breakpoints> -l <location of bamfiles> -sq <path> -pc <path> -cp <cap3> [-p <path of the outputdirectory>][-te <te sequences> -bp <blast>]' &





 Please contact jainythomas1@gmail.com for questions or support.
