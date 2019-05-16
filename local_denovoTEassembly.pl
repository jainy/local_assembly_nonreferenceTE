#!/usr/bin/perl -w
#######################################################
# Author :  Jainy Thomas
# date   :  Jan 2017
# email  :  jainythomas1@gmail.com
# Pupose :  to assemble a non-reference TE from a bam file
#           
#####################################################
use warnings;
use strict;
use Carp;
use Getopt::Long;
use File::Path qw(make_path remove_tree);
use Cwd;
use Bio::SearchIO; 
use Bio::SeqIO;
use Bio::DB::Fasta;
use Data::Dumper;
use File::Copy;
use List::MoreUtils qw(uniq);


my $version = "1.0";
my $scriptname = "local_denovoTEassembly.pl";
my $changelog = "
#   - v1.0 = does assembly of the all the reads, assembly of mates of discordant reads only 
#			with -te option, blast the contigs and singlets file  with query TE, also with the assembly made with matesonly with the query TE
#		
\n";

my $usage = "\nUsage [$version]: 
    perl $scriptname -t <table> -f <file with cordinates> -l <location of bamfiles>[-p <path of the outputdirectory>][-te <te sequences>][-v] [-c] [-h] [s] 
	
	
	
    MANDATORY ARGUMENT:	
    -t,--table (STRING) file contain accession information first column needs to be the IDs, second column BAMIDs
    -f,--file  (STRING) file containing accession information 
    -l,--bamloc(STRING)	location of bam files
    -sq,--seqtk(STRING)  path to seqtk
    -pc,--picard(STRING) path to picard tools
    -cp,--cap3  (STRING) path to cap3
    -bp,--blast (STRING) path to blast	
    	  
    OPTIONAL ARGUMENTS:
    -p,--path   (STRING) output directory name (path)
                         Default = <current working directory>
    -te,--TEfile(STRING) path to transposable element file (performs blast of assembly sequence against TE sequence provided)
    -c,--chlog  (BOOL)   Print updates
    -v,--v      (BOOL)   Print version if only option
    -s,--verbose(BOOL)   The script will talk to you
    -h,--help>  (BOOL)   Print this usage\n\n";
   


#-----------------------------------------------------------------------------
#------------------------------ LOAD AND CHECK -------------------------------
#-----------------------------------------------------------------------------
my ($file,$path,$table,$teseq,$bamlocation,$seqtkpro,$picardpro,$CAP3pro,$BLASTpro,$verbose,$help,$v,$chlog);
GetOptions ('f=s' => \$file,
            'p=s' => \$path,
            't=s' => \$table,
           'te=s' => \$teseq,
           	'l=s' => \$bamlocation,
           'sq=s' => \$seqtkpro,
           'pc=s' => \$picardpro,
           'cp=s' => \$CAP3pro,
           'bp=s' => \$BLASTpro,	
            'c'   => \$chlog, 
            'h'   => \$help,
            's'   => \$verbose, 
            'v'   => \$v);

#check step to see if mandatory argument is provided + if help/changelog
die "\n Script $scriptname version $version\n\n" if ((! $table) && (! $file) && (! $help) && (! $chlog) && ($v));
die $changelog if ($chlog);
die $usage if ((! $table) || (! $file) ||  ($help));
my $cwd = getcwd();
$path = $cwd if (!$path) ;

#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------


my %bamfile = ();
my $bamid ;
my $extractgenomicseqout;
my %tenames = ();
my $uniqueid;
my $dispath;
my $individual;
my $genomeloc;
my @alldisreads =();
my %hashindividual =();
my %discordantmatelist =();
my $discordmatefile;
my $indiallreads;
my $dismateseqfilepath;
my $dismateIDlistpath;
my $outfile = "$file.withnodisreads.txt";

%bamfile = load_file ($table);
#print "\n",%bamfile, "\n";
#%tenames = &load_tenames ($teseq);

open (my $fh, "<", $file) or confess "\n ERROR (main): could not open to read $file $!\n";
	while(<$fh>) {
		chomp (my $line = $_);
		my @col = split(/\s+/,$line);
		$individual = $col[0];
		my $chr = $col[1] ;
		my $start = $col[2] - 250;
		my $end = $col[2] + 250;
		$genomeloc = $chr.":".$start."-".$end;
		my $tetype = $col[3];
		&find_bamid();
		$uniqueid = $individual.".".$genomeloc;

		
		##Extracting discordant reads
		make_path  ("$path/Discordantreads/$individual");
		system ("samtools view -b -q 20 -F 3854 $bamlocation/$bamid $genomeloc > $path/Discordantreads/$individual/$uniqueid.discordantF3854.outfile.bam") == 0 or die ("unable to extract readsby F 3854 flag $!");
		system ("samtools bam2fq $path/Discordantreads/$individual/$uniqueid.discordantF3854.outfile.bam | $seqtkpro/seqtk seq -A -q20 > $path/Discordantreads/$individual/$uniqueid.dismapped.reads.fasta") == 0 or die ("unable to convert discordant bam file  to fasta $uniqueid \n");
		#identifying mates of discordant reads
		$dispath = "$path/Discordantreads/$individual";
		my $mappedisreads = "$uniqueid.dismapped.reads.fasta";
		%discordantmatelist = &load_readIDs ($dispath,$mappedisreads);
		#print discordant read mates to a file
		make_path  ("$path/Discordantreads/$file/dismateIDLists");	
		$discordmatefile = "$path/Discordantreads/$file/dismateIDLists/$uniqueid.discordantmatesreadIDlist.txt";
		&print_hash(%discordantmatelist);
		#loading the all the discordant reads that needs to be extracted to a hash
		&collectreads_indi();
		
		
	}
#print Dumper %discordantmatelist, "\n";	
#print Dumper %hashindividual, "\n";			
close $fh;


make_path  ("$path/Allreadsforbam/$file");
$indiallreads = "$path/Allreadsforbam/$file"; 
&printreads_indi();
&extractreads_bam();
make_path ("$path/Discordantreads/$file/discordantmatesonly");
$dismateseqfilepath = "$path/Discordantreads/$file/discordantmatesonly"; 
$dismateIDlistpath = "$path/Discordantreads/$file/dismateIDLists";
&extractdiscordmates();


my $tetype;

my $allreadrenamedfile;
my $disrenamedfile;

my $allreadnosoftrenamedfile;
my @disreadszerofiles;

open ($fh, "<", $file) or confess "\n ERROR (main): could not open to read $file $!\n";
	while(<$fh>) {
		chomp (my $line = $_);
		my @col = split(/\s+/,$line);
		my $numcol = @col;
		#print "the number of colums is $numcol\n";
		$individual = $col[0];
		my $chr = $col[1] ;
		my $start = $col[2] - 250;
		my $end = $col[2] + 250;
		$genomeloc = $chr.":".$start."-".$end;
		$tetype = $col[3];
		$tetype = "TE" if (not defined ($col[3]));
		&find_bamid();
		$uniqueid = $individual.".".$genomeloc;
		
		
			
		#Extracting reads for the assembling the reads
		
		make_path  ("$path/ExtractedReads/$tetype") ;
		make_path("$path/Assembly/allreads");		
		#extracting only the mapped reads
		unless (-e "$path/ExtractedReads/$tetype/$uniqueid.onlymappedreadIDs.fasta") {
			system("samtools view -b -q 20 -F 4 $bamlocation/$bamid $genomeloc > $path/ExtractedReads/$tetype/$uniqueid.onlymappedreadIDs.bam ") == 0 or die ("unable to extract bam at $uniqueid \n");
			system("samtools bam2fq $path/ExtractedReads/$tetype/$uniqueid.onlymappedreadIDs.bam | $seqtkpro/seqtk seq -A -q20  > $path/ExtractedReads/$tetype/$uniqueid.onlymappedreadIDs.fasta") == 0 or die ("unable to convert to fasta $uniqueid \n");
		}
	
		if (-e "$path/Discordantreads/$file/discordantmatesonly/$uniqueid.discordantmatesreadIDlist.txt.fasta") {
			#copy the corresponding the discordant mate reads
			copy("$path/Discordantreads/$file/discordantmatesonly/$uniqueid.discordantmatesreadIDlist.txt.fasta", "$path/ExtractedReads/$tetype/$uniqueid.discordantmatesreadIDlist.txt.fasta") or die "Copy failed:$!"; 
			#concatenate mapped and discordant reads
			system ("cat $path/ExtractedReads/$tetype/$uniqueid.onlymappedreadIDs.fasta $path/ExtractedReads/$tetype/$uniqueid.discordantmatesreadIDlist.txt.fasta > $path/Assembly/allreads/$uniqueid.allpairedendread.fasta") == 0 or die ("unable to concatenate allreads for assembly $uniqueid $! \n");
		} else {
			copy("$path/ExtractedReads/$tetype/$uniqueid.onlymappedreadIDs.fasta", "$path/Assembly/allreads/$uniqueid.allpairedendread.fasta") or die "Copy failed:$!"; 
		}
		#Assemble the reads to contigs
		system ("$CAP3pro/cap3 $path/Assembly/allreads/$uniqueid.allpairedendread.fasta > $path/Assembly/allreads/$uniqueid.allpairedendread.cap3.assembld.fasta") == 0 or die ("unable to assemble fasta $uniqueid \n");
		print STDERR " Assembling fasta done \n" if ($verbose);
		#concatenate the contigs and singlets
		system ("cat $path/Assembly/allreads/$uniqueid.allpairedendread.fasta.cap.contigs $path/Assembly/allreads/$uniqueid.allpairedendread.fasta.cap.singlets > $path/Assembly/allreads/$uniqueid.allpairedendread.fasta.cap.contigs_singlets.fasta") == 0 or die ("unable to concatenate contigs and singlets $uniqueid \n");
		print STDERR " concatenation done \n" if ($verbose);
		
		#rename the fasta sequence with its filename
		my $allpath = "$path/Assembly/allreads";
		my $assembledfile = "$uniqueid.allpairedendread.fasta.cap.contigs_singlets.fasta";
		&renameseq_filename ($assembledfile,$allpath); 
		#Assembling the discordant reads
		make_path  ("$path/Discordantreads/$tetype");
		make_path("$path/Assembly/disreadmatesonly");
		#copy the corresponding the discordant mate reads
		if (-e "$path/Discordantreads/$file/discordantmatesonly/$uniqueid.discordantmatesreadIDlist.txt.fasta") {
			copy("$path/Discordantreads/$file/discordantmatesonly/$uniqueid.discordantmatesreadIDlist.txt.fasta", "$path/Assembly/disreadmatesonly/$uniqueid.discordantmatesreadIDlist.txt.fasta") or die "Copy failed:$!"; 
			my $filesize = (-s "$path/Discordantreads/$file/discordantmatesonly/$uniqueid.discordantmatesreadIDlist.txt.fasta");
			if ($filesize > 0) {
				system ("$CAP3pro/cap3 $path/Assembly/disreadmatesonly/$uniqueid.discordantmatesreadIDlist.txt.fasta > $path/Assembly/disreadmatesonly/$uniqueid.discordant.reads.assembld.fasta") == 0 or die ("unable to assemble fasta $uniqueid \n");
				system ("cat $path/Assembly/disreadmatesonly/$uniqueid.discordantmatesreadIDlist.txt.fasta.cap.contigs $path/Assembly/disreadmatesonly/$uniqueid.discordantmatesreadIDlist.txt.fasta.cap.singlets > $path/Assembly/disreadmatesonly/$uniqueid.discordant.reads.cap.contigs_singlets.fasta") == 0 or die ("unable to concatenate contigs and singlets $uniqueid \n");
				#discordantreads-renaming
				my $dismatpath = "$path/Assembly/disreadmatesonly";
				my $disassembdfile = "$uniqueid.discordant.reads.cap.contigs_singlets.fasta";
				&renameseq_filename ($disassembdfile,$dismatpath); 
			} else {
				push (@disreadszerofiles, $uniqueid);
			}
		}	
		
		print  " Please find the assembly of sequences in the $path/Assembly/allreads/.......contigs_singlets.fasta  \n" ;	
		
		#blast the discordant read assembly with TEseq
		
		if (defined ($teseq)) {
			make_path  ("$path/contigTEblast/$tetype");
			$allreadrenamedfile= "$path/Assembly/allreads/Renamedcontigs/$uniqueid.rename.fasta";
			print STDERR "performing BLAST with contigs on TEs......\n";
			my $contigTEblastout = "$path/contigTEblast/$tetype/$uniqueid.rename.te.blast.out";
			unless (-e "$allreadrenamedfile.nhr") {
				system ("$BLASTpro/makeblastdb -in $allreadrenamedfile -dbtype nucl") == 0 or die ("unable to makeblastdb on $allreadrenamedfile \n");
			}
			#system ("$BLASTpro/blastn -db $teseq -query $allreadrenamedfile -evalue 0.0001 -out $contigTEblastout") == 0 or die ("unable to to perform blast $uniqueid \n");
			system ("$BLASTpro/blastn -db $allreadrenamedfile -query $teseq -evalue 0.001 -outfmt 6 -out $contigTEblastout.tabular.out") == 0 or die ("unable to to perform tabular blast $uniqueid \n");#fortableoutput
			$disrenamedfile = "$path/Assembly/disreadmatesonly/Renamedcontigs/$uniqueid.rename.fasta";
			if (-e $disrenamedfile) {
				unless (-z $disrenamedfile) {
					make_path  ("$path/DisreadsTEblast/$tetype");
					unless (-e "$disrenamedfile.nhr") {
						system ("$BLASTpro/makeblastdb -in $disrenamedfile -dbtype nucl") == 0 or die ("unable to makeblastdb on $disrenamedfile \n");
					}
					print STDERR "performing BLAST with disreads on TEs......\n";
					my $disTEblastout = "$path/DisreadsTEblast/$tetype/$uniqueid.disreads.te.blast.out";
					#system ("$BLASTpro/blastn -db $teseq -query $disrenamedfile -evalue 0.0001 -out $disTEblastout") == 0 or die ("unable to to perform blast $uniqueid \n");
					system ("$BLASTpro/blastn -db $teseq -query $disrenamedfile -evalue 0.001 -outfmt 6 -out $disTEblastout.tabular.out") == 0 or die ("unable to to perform tabular blast $uniqueid \n");#fortableoutput
					
				} else {
					print STDERR "No discordant reads identified for $uniqueid\n";
				}
			}
		}
	}
#&email;
close $fh;
print_array(@disreadszerofiles);#if need to print out candidates that does not have any discordant reads identified

exit;
#-----------------------------------------------------------------------------
#----------------------------------- SUB -------------------------------------
#-----------------------------------------------------------------------------
sub load_file {
	my ($file1) = @_;
	my %sbamfile;
	open (my $th, "<", $file1) or confess "\n ERROR (main): could not open to read $file1 $!\n";
		while (my $data = <$th>) { #reading the table
			chomp $data; #storing the table values to the $data and removing the extraline
			my @namebam = split(/\s+/,$data); # splitting the data based on tab and storing into the arrray
			my $name = $namebam[0];
			$sbamfile{$name} = $namebam[1]; #loading to the hash
		}
	return (%sbamfile);
	close $th;	
}
sub find_bamid {
	if (exists ($bamfile { $individual } ))  {
		$bamid = $bamfile{$individual};
		print STDERR " the file now analysing is $bamid \n" if ($verbose);
	}
	else {
		die " Not able to identify the bamid of the individual. Please check the input files $! \n" ;
	}
}
sub load_readIDs {
	my ($dpath,$mappedid) = @_;
	open (my $idh,"<","$dpath/$mappedid") || die ("unable to open $mappedid $! \n ");
	my %idlist;
	while (<$idh>) {
		 my $head = $_;
		 chomp $head;
		if ($head =~ /^\>(.*)\/(\d)/) {
			$head = "$1\/2" if ($2 == 1);
			$head = "$1\/1" if ($2 == 2);
			$idlist{$head} =1;
		} elsif ($head =~ /^\>(.*)/)  {
		 	$head = "$1\/0";
		 	$idlist{$head} =1;
		} else { 
			next;
		}
	}
	return (%idlist);
	close $idh;
}
sub print_hash {
	my %hashtoprint = @_;
	open (my $hp,">","$discordmatefile") || die ("failed to open file to write discordant mate list $!\n");
	foreach my $mate (sort keys %hashtoprint) {
		print $hp "$mate\n";
	}
	close $hp;
}
sub collectreads_indi {#collect all reads for an individual
	@alldisreads =();
	foreach my $read (sort keys %discordantmatelist) {
		$read = substr $read, 0,-2;#remove /1or2 from the file
		push (@{$hashindividual{$individual}},$read);
		(@{$hashindividual{$individual}}) = uniq (@{$hashindividual{$individual}});#for making the array unique
	}
	return (%hashindividual);
}
sub printreads_indi {
	foreach my $indi (sort keys %hashindividual) {
		open (my $ih, ">","$indiallreads/$indi.allreadIDs.txt" ) or die ("cannot write file $indi.allreadIDs.txt $!\n");
		foreach my $reads (@{$hashindividual{$indi}}) {
			print $ih "$reads\n";
		}
		close $ih;
	}
}
sub extractreads_bam {
	my @indifiles = `ls $path/Allreadsforbam/$file`;
	my $bam_id;
	foreach my $indivi (@indifiles) {
		my @allreadsfile_name = split (/\./,$indivi);
		my $allreadindi = $allreadsfile_name[0];
		if (exists ($bamfile{$allreadindi})) {
			$bam_id = $bamfile{$allreadindi};
		} else {
			print STDERR "bam_id cannot be identified for $allreadindi\n";
		}
		#extract reads using picard tools
		unless (-e "$path/Allreadsforbam/$file/$allreadindi.allreadIDs.bam") {
			system ("java -jar $picardpro/picard.jar FilterSamReads INPUT=$bamlocation/$bam_id VALIDATION_STRINGENCY=LENIENT FILTER=includeReadList READ_LIST_FILE=$path/Allreadsforbam/$file/$allreadindi.allreadIDs.txt WRITE_READS_FILES=false OUTPUT=$path/Allreadsforbam/$file/$allreadindi.allreadIDs.bam") == 0 or die ("unable to run picard tools in $file on $allreadindi \n");
			system ("samtools bam2fq $path/Allreadsforbam/$file/$allreadindi.allreadIDs.bam | $seqtkpro/seqtk seq -A -q20 > $path/Allreadsforbam/$file/$allreadindi.allreadIDs.fasta") == 0 or die ("unable to convert to fasta $allreadindi \n");

		}
	}
}
sub extractdiscordmates {
	my @dismates = `ls $dismateIDlistpath`;
	foreach my $matefile (@dismates) {
		chomp $matefile;
		my @matefilename = split (/\./,$matefile);
		my $indivifilename = $matefilename[0];
		my %mates = ();
		open (my $mh, "<", "$dismateIDlistpath/$matefile") or confess "\n ERROR (main): could not open to read $matefile $!\n";
		while (my $dataline = <$mh>) { 
			chomp $dataline; 
			$mates{$dataline} = 1; 
		}
		#print Dumper %mates, "\n";
		
		if (-e "$path/Allreadsforbam/$file/$indivifilename.allreadIDs.fasta") {
			my $readio_obj = Bio::SeqIO->new(-file 	 => "$path/Allreadsforbam/$file/$indivifilename.allreadIDs.fasta", 
											 -format => 'fasta') 
										 or die "\t    ERROR - Failed to create SeqIO FH object from $indivifilename.allreadIDs.fasta $!\n";  
			my $outreadio_obj = Bio::SeqIO->new(-file   => ">$dismateseqfilepath/$matefile.fasta",
												-format => 'fasta') 
											 or die "\t    ERROR - Failed to create SeqIO FH object from $dismateseqfilepath/$matefile.fasta $!\n";  
			while (my $seq = $readio_obj->next_seq() ){
				my $header = $seq->display_id;
				if (exists $mates{$header}) {
					$outreadio_obj->write_seq($seq);
				} else {
					next;
				}
			}
		}
	}
}
sub print_array {
	my (@array)= @_;
	open (my $ah,">","$path/$outfile") || die ("cannot open file $outfile to write $!\n");
		foreach my $dataline (@array) {
			print $ah "$dataline\n";
		}
	close ($ah);
}

sub renameseq_filename {
	my ($contigfile,$fpath) = @_;
	make_path  ("$fpath/Renamedcontigs");
	open (my $bhout, ">","$fpath/Renamedcontigs/$uniqueid.rename.fasta") or die "\n ERROR (main): could not open to read  $!\n";
	open (my $bh, "<", "$fpath/$contigfile") or confess "\n ERROR (main): could not open to read $contigfile $!\n";
		while(my $dataline = <$bh>) {
			chomp($dataline);
			#print STDERR "$line\n";
				if ($dataline =~ m/^\>\w+\d+/) {
				$dataline =~ s/^\>(\w+\d+)/\>$uniqueid\.$1/;
				print $bhout "$dataline\n";
				}
				else {
				print $bhout "$dataline\n";
				}
		}
	close $bh;
	close $bhout;
}

