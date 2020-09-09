$HOSTNAME = ""
params.outdir = 'results'  

//pipeline defaults
params.nucleicAcidType = "rna"


if (!params.mate){params.mate = ""} 
if (!params.reads){params.reads = ""} 

Channel.value(params.mate).into{g_1_mate_g_4;g_1_mate_g_41;g_1_mate_g0_3;g_1_mate_g0_11;g_1_mate_g0_16;g_1_mate_g0_18;g_1_mate_g0_19;g_1_mate_g0_20;g_1_mate_g0_21;g_1_mate_g42_26;g_1_mate_g42_30;g_1_mate_g42_32}
Channel
	.fromFilePairs( params.reads , size: (params.mate != "pair") ? 1 : 2 )
	.ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
	.into{g_2_reads_g0_3;g_2_reads_g0_18}


params.run_Adapter_Removal =   "no"   //* @dropdown @options:"yes","no" @show_settings:"Adapter_Removal"
//* @style @multicolumn:{seed_mismatches, palindrome_clip_threshold, simple_clip_threshold} @condition:{Tool_for_Adapter_Removal="trimmomatic", seed_mismatches, palindrome_clip_threshold, simple_clip_threshold}, {Tool_for_Adapter_Removal="fastx_clipper", discard_non_clipped}

//* autofill
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 1000
    $CPU  = 1
    $MEMORY = 24
    $QUEUE = "long"
}
//* platform
//* autofill
if (!((params.run_Adapter_Removal && (params.run_Adapter_Removal == "yes")) || !params.run_Adapter_Removal)){
g_2_reads_g0_18.into{g0_18_reads_g0_19}
g0_18_log_file_g0_11 = Channel.empty()
} else {


process Adapter_Trimmer_Quality_Module_Adapter_Removal {

input:
 set val(name), file(reads) from g_2_reads_g0_18
 val mate from g_1_mate_g0_18

output:
 set val(name), file("reads/*.fastq")  into g0_18_reads_g0_19
 file "*.{fastx,trimmomatic}.log"  into g0_18_log_file_g0_11

errorStrategy 'retry'

when:
(params.run_Adapter_Removal && (params.run_Adapter_Removal == "yes")) || !params.run_Adapter_Removal

shell:
phred = params.Adapter_Trimmer_Quality_Module_Adapter_Removal.phred
Tool_for_Adapter_Removal = params.Adapter_Trimmer_Quality_Module_Adapter_Removal.Tool_for_Adapter_Removal
Adapter_Sequence = params.Adapter_Trimmer_Quality_Module_Adapter_Removal.Adapter_Sequence
//trimmomatic_inputs
min_length = params.Adapter_Trimmer_Quality_Module_Adapter_Removal.min_length
seed_mismatches = params.Adapter_Trimmer_Quality_Module_Adapter_Removal.seed_mismatches
palindrome_clip_threshold = params.Adapter_Trimmer_Quality_Module_Adapter_Removal.palindrome_clip_threshold
simple_clip_threshold = params.Adapter_Trimmer_Quality_Module_Adapter_Removal.simple_clip_threshold

//fastx_clipper_inputs
discard_non_clipped = params.Adapter_Trimmer_Quality_Module_Adapter_Removal.discard_non_clipped
    
remove_previous_reads = params.Adapter_Trimmer_Quality_Module_Adapter_Removal.remove_previous_reads
discard_non_clipped_text = ""
if (discard_non_clipped == "yes") {discard_non_clipped_text = "-c"}
nameAll = reads.toString()
nameArray = nameAll.split(' ')
file2 = ""
if (nameAll.contains('.gz')) {
    newName =  nameArray[0] - ~/(\.fastq.gz)?(\.fq.gz)?$/
    file1 =  nameArray[0] - '.gz' 
    if (mate == "pair") {file2 =  nameArray[1] - '.gz'}
    runGzip = "ls *.gz | xargs -i echo gzip -df {} | sh"
} else {
    newName =  nameArray[0] - ~/(\.fastq)?(\.fq)?$/
    file1 =  nameArray[0]
    if (mate == "pair") {file2 =  nameArray[1]}
    runGzip = ''
}
'''
#!/usr/bin/env perl
 use List::Util qw[min max];
 use strict;
 use File::Basename;
 use Getopt::Long;
 use Pod::Usage;
 use Cwd qw();
 
runCmd("mkdir reads adapter unpaired");

open(OUT, ">adapter/adapter.fa");
my @adaps=split(/\n/,"!{Adapter_Sequence}");
my $i=1;
foreach my $adap (@adaps)
{
 print OUT ">adapter$i\\n$adap\\n";
 $i++;
}
close(OUT);

system("!{runGzip}");
my $quality="!{phred}";
print "fastq quality: $quality\\n";
print "tool: !{Tool_for_Adapter_Removal}\\n";

if ("!{mate}" eq "pair") {
    if ("!{Tool_for_Adapter_Removal}" eq "trimmomatic") {
        runCmd("trimmomatic PE -threads 1 -phred${quality} !{file1} !{file2} reads/!{name}.1.fastq unpaired/!{name}.1.fastq.unpaired reads/!{name}.2.fastq unpaired/!{name}.2.fastq.unpaired ILLUMINACLIP:adapter/adapter.fa:!{seed_mismatches}:!{palindrome_clip_threshold}:!{simple_clip_threshold} MINLEN:!{min_length} 2> !{name}.trimmomatic.log");
    } elsif ("!{Tool_for_Adapter_Removal}" eq "fastx_clipper") {
        print "Fastx_clipper is not suitable for paired reads.";
    }
} else {
    if ("!{Tool_for_Adapter_Removal}" eq "trimmomatic") {
        runCmd("trimmomatic SE -threads 1  -phred${quality} !{file1} reads/!{name}.fastq ILLUMINACLIP:adapter/adapter.fa:!{seed_mismatches}:!{palindrome_clip_threshold}:!{simple_clip_threshold} MINLEN:!{min_length} 2> !{name}.trimmomatic.log");
    } elsif ("!{Tool_for_Adapter_Removal}" eq "fastx_clipper") {
        runCmd("fastx_clipper  -Q $quality -a !{Adapter_Sequence} -l !{min_length} !{discard_non_clipped_text} -v -i !{file1} -o reads/!{name}.fastq > !{name}.fastx.log");
    }
}
if ("!{remove_previous_reads}" eq "true") {
    my $currpath = Cwd::cwd();
    my @paths = (split '/', $currpath);
    splice(@paths, -2);
    my $workdir= join '/', @paths;
    splice(@paths, -1);
    my $inputsdir = join '/', @paths;
    $inputsdir .= "/inputs";
    print "INFO: inputs reads will be removed if they are located in the workdir inputsdir\\n";
    my @listOfFiles = `readlink -e !{file1} !{file2}`;
    foreach my $targetFile (@listOfFiles){
        if (index($targetFile, $workdir) != -1 || index($targetFile, $inputsdir) != -1) {
            system("rm -f $targetFile");
            print "INFO: $targetFile deleted.\\n";
        }
    }
}


##Subroutines
sub runCmd {
    my ($com) = @_;
    my $error = system($com);
    if   ($error) { die "Command failed: $error $com\\n"; }
    else          { print "Command successful: $com\\n"; }
}
'''

}
}


params.run_Trimmer =   "no"   //* @dropdown @options:"yes","no" @show_settings:"Trimmer"
//* @style @multicolumn:{trim_length_5prime,trim_length_3prime}, {trim_length_5prime_R1,trim_length_3prime_R1}, {trim_length_5prime_R2,trim_length_3prime_R2} @condition:{single_or_paired_end_reads="single", trim_length_5prime,trim_length_3prime}, {single_or_paired_end_reads="pair", trim_length_5prime_R1,trim_length_3prime_R1,trim_length_5prime_R2,trim_length_3prime_R2}

//* autofill
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 500
    $CPU  = 1
    $MEMORY = 8
    $QUEUE = "long"
}
//* platform
//* autofill
if (!((params.run_Trimmer && (params.run_Trimmer == "yes")) || !params.run_Trimmer)){
g0_18_reads_g0_19.into{g0_19_reads_g0_20}
g0_19_log_file_g0_21 = Channel.empty()
} else {


process Adapter_Trimmer_Quality_Module_Trimmer {

input:
 set val(name), file(reads) from g0_18_reads_g0_19
 val mate from g_1_mate_g0_19

output:
 set val(name), file("reads/*q")  into g0_19_reads_g0_20
 file "*.log" optional true  into g0_19_log_file_g0_21

errorStrategy 'retry'

when:
(params.run_Trimmer && (params.run_Trimmer == "yes")) || !params.run_Trimmer

shell:
phred = params.Adapter_Trimmer_Quality_Module_Trimmer.phred
single_or_paired_end_reads = params.Adapter_Trimmer_Quality_Module_Trimmer.single_or_paired_end_reads
trim_length_5prime = params.Adapter_Trimmer_Quality_Module_Trimmer.trim_length_5prime
trim_length_3prime = params.Adapter_Trimmer_Quality_Module_Trimmer.trim_length_3prime
trim_length_5prime_R1 = params.Adapter_Trimmer_Quality_Module_Trimmer.trim_length_5prime_R1
trim_length_3prime_R1 = params.Adapter_Trimmer_Quality_Module_Trimmer.trim_length_3prime_R1
trim_length_5prime_R2 = params.Adapter_Trimmer_Quality_Module_Trimmer.trim_length_5prime_R2
trim_length_3prime_R2 = params.Adapter_Trimmer_Quality_Module_Trimmer.trim_length_3prime_R2
remove_previous_reads = params.Adapter_Trimmer_Quality_Module_Trimmer.remove_previous_reads

nameAll = reads.toString()
nameArray = nameAll.split(' ')
file2 = ""
if (nameAll.contains('.gz')) {
    newName =  nameArray[0] - ~/(\.fastq.gz)?(\.fq.gz)?$/
    file1 =  nameArray[0] - '.gz' 
    if (mate == "pair") {file2 =  nameArray[1] - '.gz'}
    runGzip = "ls *.gz | xargs -i echo gzip -df {} | sh"
} else {
    newName =  nameArray[0] - ~/(\.fastq)?(\.fq)?$/
    file1 =  nameArray[0]
    if (mate == "pair") {file2 =  nameArray[1]}
    runGzip = ''
}
'''
#!/usr/bin/env perl
 use List::Util qw[min max];
 use strict;
 use File::Basename;
 use Getopt::Long;
 use Pod::Usage; 
 use Cwd qw();
 
system("mkdir reads");
system("!{runGzip}");
my $file1 = "";
my $file2 = "";
if ("!{mate}" eq "pair") {
    $file1 = "!{file1}";
    $file2 = "!{file2}";
    my $trim1 = "!{trim_length_5prime_R1}:!{trim_length_3prime_R1}";
    my $trim2 = "!{trim_length_5prime_R2}:!{trim_length_3prime_R2}";
    my $len=getLength($file1);
    print "length of $file1: $len\\n";
    trimFiles($file1, $trim1, $len);
    my $len=getLength($file2);
    print "INFO: length of $file2: $len\\n";
    trimFiles($file2, $trim2, $len);
} else {
    $file1 = "!{file1}";
    my $trim1 = "!{trim_length_5prime}:!{trim_length_3prime}";
    my $len=getLength($file1);
    print "INFO: length of file1: $len\\n";
    trimFiles($file1, $trim1, $len);
}
if ("!{remove_previous_reads}" eq "true") {
    my $currpath = Cwd::cwd();
    my @paths = (split '/', $currpath);
    splice(@paths, -2);
    my $workdir= join '/', @paths;
    splice(@paths, -1);
    my $inputsdir= join '/', @paths;
    $inputsdir .= "/inputs";
    print "INFO: inputs reads will be removed if they are located in the workdir inputsdir\\n";
    my @listOfFiles = `readlink -e !{file1} !{file2}`;
    foreach my $targetFile (@listOfFiles){
        if (index($targetFile, $workdir) != -1 || index($targetFile, $inputsdir) != -1) {
            system("rm -f $targetFile");
            print "INFO: $targetFile deleted.\\n";
        }
    }
}



sub trimFiles
{
  my ($file, $trim, $len)=@_;
    my @nts=split(/[,:\\s\\t]+/,$trim);
    my $inpfile="";
    my $com="";
    my $i=1;
    my $outfile="";
    my $param="";
    my $quality="-Q!{phred}";

    if (scalar(@nts)==2)
    {
      $param = "-f ".($nts[0]+1) if (exists($nts[0]) && $nts[0] >= 0 );
      $param .= " -l ".($len-$nts[1]) if (exists($nts[0]) && $nts[1] > 0 );
      $outfile="reads/$file";  
      $com="fastx_trimmer $quality -v $param -o $outfile -i $file > !{name}.fastx_trimmer.log" if ((exists($nts[0]) && $nts[0] > 0) || (exists($nts[0]) && $nts[1] > 0 ));
      print "INFO: $com\\n";
      if ($com eq ""){
          print "INFO: Trimmer skipped for $file \\n";
          system("mv $file reads/.");
      } else {
          runCmd("$com");
          print "INFO: Trimmer executed for $file \\n";
      }
    }

    
}


sub getLength
{
   my ($filename)=@_;
   open (IN, $filename);
   my $j=1;
   my $len=0;
   while(my $line=<IN>)
   {
     chomp($line);
     if ($j >50) { last;}
     if ($j%4==0)
     {
        $len=length($line);
     }
     $j++;
   }
   close(IN);
   return $len;
}

sub runCmd {
    my ($com) = @_;
    my $error = system($com);
    if   ($error) { die "Command failed: $error $com\\n"; }
    else          { print "Command successful: $com\\n"; }
}

'''

}
}



process Adapter_Trimmer_Quality_Module_Trimmer_Removal_Summary {

input:
 file logfile from g0_19_log_file_g0_21.collect()
 val mate from g_1_mate_g0_21

output:
 file "trimmer_summary.tsv"  into g0_21_outputFileTSV_g_26

shell:
'''
#!/usr/bin/env perl
use List::Util qw[min max];
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;

my @header;
my %all_files;
my %tsv;
my %tsvDetail;
my %headerHash;
my %headerText;
my %headerTextDetail;

my $i = 0;
chomp( my $contents = `ls *.log` );

my @files = split( /[\\n]+/, $contents );
foreach my $file (@files) {
    $i++;
    my $mapOrder = "1";
    if ($file =~ /(.*)\\.fastx_trimmer\\.log/){
        $file =~ /(.*)\\.fastx_trimmer\\.log/;
        my $mapper   = "fastx_trimmer";
        my $name = $1;    ##sample name
        push( @header, $mapper );
        my $in;
        my $out;
        chomp( $in =`cat $file | grep 'Input:' | awk '{sum+=\\$2} END {print sum}'` );
        chomp( $out =`cat $file | grep 'Output:' | awk '{sum+=\\$2} END {print sum}'` );

        $tsv{$name}{$mapper} = [ $in, $out ];
        $headerHash{$mapOrder} = $mapper;
        $headerText{$mapOrder} = [ "Total Reads", "Reads After Trimmer" ];
    }
}

my @mapOrderArray = ( keys %headerHash );
my @sortedOrderArray = sort { $a <=> $b } @mapOrderArray;

my $summary          = "trimmer_summary.tsv";
writeFile( $summary,          \\%headerText,       \\%tsv );

sub writeFile {
    my $summary    = $_[0];
    my %headerText = %{ $_[1] };
    my %tsv        = %{ $_[2] };
    open( OUT, ">$summary" );
    print OUT "Sample\\t";
    my @headArr = ();
    for my $mapOrder (@sortedOrderArray) {
        push( @headArr, @{ $headerText{$mapOrder} } );
    }
    my $headArrAll = join( "\\t", @headArr );
    print OUT "$headArrAll\\n";

    foreach my $name ( keys %tsv ) {
        my @rowArr = ();
        for my $mapOrder (@sortedOrderArray) {
            push( @rowArr, @{ $tsv{$name}{ $headerHash{$mapOrder} } } );
        }
        my $rowArrAll = join( "\\t", @rowArr );
        print OUT "$name\\t$rowArrAll\\n";
    }
    close(OUT);
}

'''
}

params.run_Quality_Filtering =   "no"   //* @dropdown @options:"yes","no" @show_settings:"Quality_Filtering"
//* @style @multicolumn:{window_size,required_quality}, {leading,trailing,minlen}, {minQuality,minPercent} @condition:{tool="trimmomatic", minlen, trailing, leading, required_quality_for_window_trimming, window_size}, {tool="fastx", minQuality, minPercent}

//* autofill
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 240
    $CPU  = 1
    $MEMORY = 8
    $QUEUE = "short"
}
//* platform
//* autofill
if (!((params.run_Quality_Filtering && (params.run_Quality_Filtering == "yes")) || !params.run_Quality_Filtering)){
g0_19_reads_g0_20.into{g0_20_reads_g_41}
g0_20_log_file_g0_16 = Channel.empty()
} else {


process Adapter_Trimmer_Quality_Module_Quality_Filtering {

input:
 set val(name), file(reads) from g0_19_reads_g0_20
 val mate from g_1_mate_g0_20

output:
 set val(name), file("reads/*q")  into g0_20_reads_g_41
 file "*.{fastx,trimmomatic}_quality.log" optional true  into g0_20_log_file_g0_16

errorStrategy 'retry'

when:
(params.run_Quality_Filtering && (params.run_Quality_Filtering == "yes")) || !params.run_Quality_Filtering    

shell:
tool = params.Adapter_Trimmer_Quality_Module_Quality_Filtering.tool
phred = params.Adapter_Trimmer_Quality_Module_Quality_Filtering.phred
window_size = params.Adapter_Trimmer_Quality_Module_Quality_Filtering.window_size
required_quality_for_window_trimming = params.Adapter_Trimmer_Quality_Module_Quality_Filtering.required_quality_for_window_trimming
leading = params.Adapter_Trimmer_Quality_Module_Quality_Filtering.leading
trailing = params.Adapter_Trimmer_Quality_Module_Quality_Filtering.trailing
minlen = params.Adapter_Trimmer_Quality_Module_Quality_Filtering.minlen


// fastx parameters
minQuality = params.Adapter_Trimmer_Quality_Module_Quality_Filtering.minQuality
minPercent = params.Adapter_Trimmer_Quality_Module_Quality_Filtering.minPercent

remove_previous_reads = params.Adapter_Trimmer_Quality_Module_Quality_Filtering.remove_previous_reads

nameAll = reads.toString()
nameArray = nameAll.split(' ')
file2 ="";
if (nameAll.contains('.gz')) {
    newName =  nameArray[0] - ~/(\.fastq.gz)?(\.fq.gz)?$/
    file1 =  nameArray[0] - '.gz' 
    if (mate == "pair") {file2 =  nameArray[1] - '.gz'}
    runGzip = "ls *.gz | xargs -i echo gzip -df {} | sh"
} else {
    newName =  nameArray[0] - ~/(\.fastq)?(\.fq)?$/
    file1 =  nameArray[0]
    if (mate == "pair") {file2 =  nameArray[1]}
    runGzip = ''
}
'''
#!/usr/bin/env perl
 use List::Util qw[min max];
 use strict;
 use File::Basename;
 use Getopt::Long;
 use Pod::Usage; 
 use Cwd qw();
 
system("mkdir reads unpaired");
system("!{runGzip}");
my $param = "SLIDINGWINDOW:"."!{window_size}".":"."!{required_quality_for_window_trimming}";
$param.=" LEADING:"."!{leading}";
$param.=" TRAILING:"."!{trailing}";
$param.=" MINLEN:"."!{minlen}";

my $quality="!{phred}";

print "INFO: fastq quality: $quality\\n";
     
if ("!{tool}" eq "trimmomatic") {
    if ("!{mate}" eq "pair") {
        runCmd("trimmomatic PE -phred${quality} !{file1} !{file2} reads/!{name}.1.fastq unpaired/!{name}.1.fastq.unpaired reads/!{name}.2.fastq unpaired/!{name}.1.fastq.unpaired $param 2> !{name}.trimmomatic_quality.log");
    } else {
        runCmd("trimmomatic SE -phred${quality} !{file1} reads/!{name}.fastq $param 2> !{name}.trimmomatic_quality.log");
    }
} elsif ("!{tool}" eq "fastx") {
    if ("!{mate}" eq "pair") {
        print("WARNING: Fastx option is not suitable for paired reads. This step will be skipped.");
        system("mv !{file1} !{file2} reads/.");
    } else {
        runCmd("fastq_quality_filter  -Q $quality -q !{minQuality} -p !{minPercent} -v -i !{file1} -o reads/!{name}.fastq > !{name}.fastx_quality.log");
    }
}
if ("!{remove_previous_reads}" eq "true") {
    my $currpath = Cwd::cwd();
    my @paths = (split '/', $currpath);
    splice(@paths, -2);
    my $workdir= join '/', @paths;
    splice(@paths, -1);
    my $inputsdir= join '/', @paths;
    $inputsdir .= "/inputs";
    print "INFO: inputs reads will be removed if they are located in the workdir inputsdir\\n";
    my @listOfFiles = `readlink -e !{file1} !{file2}`;
    foreach my $targetFile (@listOfFiles){
        if (index($targetFile, $workdir) != -1 || index($targetFile, $inputsdir) != -1) {
            system("rm -f $targetFile");
            print "INFO: $targetFile deleted.\\n";
        }
    }
}

##Subroutines
sub runCmd {
    my ($com) = @_;
    my $error = system($com);
    if   ($error) { die "Command failed: $error $com\\n"; }
    else          { print "Command successful: $com\\n"; }
}


'''

}
}



process Adapter_Trimmer_Quality_Module_Quality_Filtering_Summary {

input:
 file logfile from g0_20_log_file_g0_16.collect()
 val mate from g_1_mate_g0_16

output:
 file "quality_filter_summary.tsv"  into g0_16_outputFileTSV_g_26

shell:
'''
#!/usr/bin/env perl
use List::Util qw[min max];
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;

my @header;
my %all_files;
my %tsv;
my %headerHash;
my %headerText;

my $i = 0;
chomp( my $contents = `ls *.log` );
my @files = split( /[\\n]+/, $contents );
foreach my $file (@files) {
    $i++;
    my $mapper   = "";
    my $mapOrder = "1";
    if ($file =~ /(.*)\\.fastx_quality\\.log/){
        $mapper   = "fastx";
        $file =~ /(.*)\\.fastx_quality\\.log/;
        my $name = $1;    ##sample name
        push( @header, $mapper );
        my $in;
        my $out;
        chomp( $in =`cat $file | grep 'Input:' | awk '{sum+=\\$2} END {print sum}'` );
        chomp( $out =`cat $file | grep 'Output:' | awk '{sum+=\\$2} END {print sum}'` );
        $tsv{$name}{$mapper} = [ $in, $out ];
        $headerHash{$mapOrder} = $mapper;
        $headerText{$mapOrder} = [ "Total Reads", "Reads After Quality Filtering" ];
    } elsif ($file =~ /(.*)\\.trimmomatic_quality\\.log/){
        $mapper   = "trimmomatic";
        $file =~ /(.*)\\.trimmomatic_quality\\.log/;
        my $name = $1;    ##sample name
        push( @header, $mapper );
        my $in;
        my $out;
        if ( "!{mate}" eq "pair"){
            chomp( $in =`cat $file | grep 'Input Read Pairs:' | awk '{sum+=\\$4} END {print sum}'` );
            chomp( $out =`cat $file | grep 'Input Read Pairs:' | awk '{sum+=\\$7} END {print sum}'` );
        } else {
            chomp( $in =`cat $file | grep 'Input Reads:' | awk '{sum+=\\$3} END {print sum}'` );
            chomp( $out =`cat $file | grep 'Input Reads:' | awk '{sum+=\\$5} END {print sum}'` );
        }
        $tsv{$name}{$mapper} = [ $in, $out ];
        $headerHash{$mapOrder} = $mapper;
        $headerText{$mapOrder} = [ "Total Reads", "Reads After Quality Filtering" ];
    }
    
}

my @mapOrderArray = ( keys %headerHash );
my @sortedOrderArray = sort { $a <=> $b } @mapOrderArray;

my $summary          = "quality_filter_summary.tsv";
writeFile( $summary,          \\%headerText,       \\%tsv );

sub writeFile {
    my $summary    = $_[0];
    my %headerText = %{ $_[1] };
    my %tsv        = %{ $_[2] };
    open( OUT, ">$summary" );
    print OUT "Sample\\t";
    my @headArr = ();
    for my $mapOrder (@sortedOrderArray) {
        push( @headArr, @{ $headerText{$mapOrder} } );
    }
    my $headArrAll = join( "\\t", @headArr );
    print OUT "$headArrAll\\n";

    foreach my $name ( keys %tsv ) {
        my @rowArr = ();
        for my $mapOrder (@sortedOrderArray) {
            push( @rowArr, @{ $tsv{$name}{ $headerHash{$mapOrder} } } );
        }
        my $rowArrAll = join( "\\t", @rowArr );
        print OUT "$name\\t$rowArrAll\\n";
    }
    close(OUT);
}

'''
}


process Adapter_Trimmer_Quality_Module_Adapter_Removal_Summary {

input:
 file logfile from g0_18_log_file_g0_11.collect()
 val mate from g_1_mate_g0_11

output:
 file "adapter_removal_summary.tsv"  into g0_11_outputFileTSV_g_26
 file "adapter_removal_detailed_summary.tsv" optional true  into g0_11_outputFile

shell:
'''
#!/usr/bin/env perl
use List::Util qw[min max];
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;

my @header;
my %all_files;
my %tsv;
my %tsvDetail;
my %headerHash;
my %headerText;
my %headerTextDetail;

my $i = 0;
chomp( my $contents = `ls *.log` );

my @files = split( /[\\n]+/, $contents );
foreach my $file (@files) {
    $i++;
    my $mapOrder = "1";
    if ($file =~ /(.*)\\.fastx\\.log/){
        $file =~ /(.*)\\.fastx\\.log/;
        my $mapper   = "fastx";
        my $name = $1;    ##sample name
        push( @header, $mapper );

        my $in;
        my $out;
        my $tooshort;
        my $adapteronly;
        my $noncliped;
        my $Nreads;

        chomp( $in =`cat $file | grep 'Input:' | awk '{sum+=\\$2} END {print sum}'` );
        chomp( $out =`cat $file | grep 'Output:' | awk '{sum+=\\$2} END {print sum}'` );
        chomp( $tooshort =`cat $file | grep 'too-short reads' | awk '{sum+=\\$2} END {print sum}'`);
        chomp( $adapteronly =`cat $file | grep 'adapter-only reads' | awk '{sum+=\\$2} END {print sum}'`);
        chomp( $noncliped =`cat $file | grep 'non-clipped reads.' | awk '{sum+=\\$2} END {print sum}'`);
        chomp( $Nreads =`cat $file | grep 'N reads.' | awk '{sum+=\\$2} END {print sum}'` );

        $tsv{$name}{$mapper} = [ $in, $out ];
        $headerHash{$mapOrder} = $mapper;
        $headerText{$mapOrder} = [ "Total Reads", "Reads After Adapter Removal" ];
        $tsvDetail{$name}{$mapper} = [ $in, $tooshort, $adapteronly, $noncliped, $Nreads, $out ];
        $headerTextDetail{$mapOrder} = ["Total Reads","Too-short reads","Adapter-only reads","Non-clipped reads","N reads","Reads After Adapter Removal"];
    } elsif ($file =~ /(.*)\\.trimmomatic\\.log/){
        $file =~ /(.*)\\.trimmomatic\\.log/;
        my $mapper   = "trimmomatic";
        my $name = $1;    ##sample name
        push( @header, $mapper );
        
        my $in;
        my $out;

        if ( "!{mate}" eq "pair"){
            chomp( $in =`cat $file | grep 'Input Read Pairs:' | awk '{sum+=\\$4} END {print sum}'` );
            chomp( $out =`cat $file | grep 'Input Read Pairs:' | awk '{sum+=\\$7} END {print sum}'` );
        } else {
            chomp( $in =`cat $file | grep 'Input Reads:' | awk '{sum+=\\$3} END {print sum}'` );
            chomp( $out =`cat $file | grep 'Input Reads:' | awk '{sum+=\\$5} END {print sum}'` );
        }
        


        $tsv{$name}{$mapper} = [ $in, $out ];
        $headerHash{$mapOrder} = $mapper;
        $headerText{$mapOrder} = [ "Total Reads", "Reads After Adapter Removal" ];
        
    }
    
}

my @mapOrderArray = ( keys %headerHash );
my @sortedOrderArray = sort { $a <=> $b } @mapOrderArray;

my $summary          = "adapter_removal_summary.tsv";
my $detailed_summary = "adapter_removal_detailed_summary.tsv";
writeFile( $summary,          \\%headerText,       \\%tsv );
if (%headerTextDetail){
    writeFile( $detailed_summary, \\%headerTextDetail, \\%tsvDetail );  
}

sub writeFile {
    my $summary    = $_[0];
    my %headerText = %{ $_[1] };
    my %tsv        = %{ $_[2] };
    open( OUT, ">$summary" );
    print OUT "Sample\\t";
    my @headArr = ();
    for my $mapOrder (@sortedOrderArray) {
        push( @headArr, @{ $headerText{$mapOrder} } );
    }
    my $headArrAll = join( "\\t", @headArr );
    print OUT "$headArrAll\\n";

    foreach my $name ( keys %tsv ) {
        my @rowArr = ();
        for my $mapOrder (@sortedOrderArray) {
            push( @rowArr, @{ $tsv{$name}{ $headerHash{$mapOrder} } } );
        }
        my $rowArrAll = join( "\\t", @rowArr );
        print OUT "$name\\t$rowArrAll\\n";
    }
    close(OUT);
}

'''
}

params.gtf =  ""  //* @input
params.trna_gtf =  ""  //* @input
params.rrna_fa =  ""  //* @input
params.genome =  ""  //* @input
params.genome2bit =  ""  //* @input
params.genome_origin =  ""  //* @input
params.gtf_origin =  ""  //* @input
params.trna_gtf_origin =  ""  //* @input
params.rrna_fa_origin =  ""  //* @input

def downFile(path){
    if (path.take(1).indexOf("/") == 0){
      target=path
    } else {
      a=file(path)
      fname = a.getName().toString()
      target = "${baseDir}/work/${fname}"
      a.copyTo(target) 
    }
    return target
}

process Ribosome_Profiling_Build_Index_Check_Genome_Gtf {


output:
 val "${params.genome}"  into g40_0_genomePath
 val "${params.gtf}"  into g40_0_gtfPath_g40_1, g40_0_gtfPath_g40_3, g40_0_gtfPath_g40_5
 val "${params.trna_gtf}"  into g40_0_trnaGtfPath_g40_5
 val "${params.rrna_fa}"  into g40_0_rRnaFastaPath_g40_5

when:
params.run_checkAndBuild == "yes"

script:
gtf_dir  = params.gtf.substring(0, params.gtf.lastIndexOf('/')) 
genome_dir  = params.genome.substring(0, params.genome.lastIndexOf('/')) 
basename_and_path  = params.genome.substring(0, params.genome.lastIndexOf('.'))

downGenomePath = ""
downGtfPath = ""
downTrnaGtfPath = ""
downRrnaFastaPath = ""
if ( !file("${params.genome}").exists() ) {
	downGenomePath=downFile(params.genome_origin)
}
if ( !file("${params.gtf}").exists() ) {
	downGtfPath=downFile(params.gtf_origin)
}
if ( !file("${params.trna_gtf}").exists() ) {
	downTrnaGtfPath=downFile(params.trna_gtf_origin)
}
if ( !file("${params.rrna_fa}").exists() ) {
	downRrnaFastaPath=downFile(params.rrna_fa_origin)
}


"""
if [ ! -e "${params.gtf}" ] ; then
    echo "${params.gtf} not found"
    mkdir -p ${gtf_dir}
    cp -n ${downGtfPath} ${params.gtf}
fi
if [ ! -e "${params.trna_gtf}" ] ; then
    echo "${params.trna_gtf} not found"
    mkdir -p ${gtf_dir}
    cp -n $downTrnaGtfPath ${params.trna_gtf}
fi
if [ ! -e "${params.rrna_fa}" ] ; then
    echo "${params.rrna_fa} not found"
    mkdir -p ${gtf_dir}
    cp -n $downRrnaFastaPath ${params.rrna_fa}
fi
if [ ! -e "${params.genome}" ] ; then
    echo "${params.genome} not found"
    mkdir -p ${genome_dir}
    cp -n $downGenomePath ${params.genome}
fi
if [ ! -e "${params.genome}.fai" ] ; then
    echo "${params.genome}.fai not found"
    cd $genome_dir
    samtools faidx ${params.genome}
fi
if [ ! -e "${params.genome2bit}" ] ; then
    echo "${params.genome2bit} not found"
    cd $genome_dir
    faToTwoBit ${params.genome} ${params.genome2bit}
fi

"""




}

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 10
    $MEMORY = 30
}
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 500
    $CPU  = 10
    $MEMORY = 30
    $QUEUE = "long"
} 
//* platform
//* autofill

process Ribosome_Profiling_Build_Index_build_ncRNA_fasta {

input:
 val valid_gtf_path from g40_0_gtfPath_g40_5
 val trna_gtf_path from g40_0_trnaGtfPath_g40_5
 val rrna_fa_path from g40_0_rRnaFastaPath_g40_5

output:
 val ncRNA_output_fasta  into g40_5_ncRNAFastaPath_g40_6

script:
gtfInFilePrefix = valid_gtf_path.substring(0, valid_gtf_path.lastIndexOf('.')) 
gtfDir = valid_gtf_path.substring(0, valid_gtf_path.lastIndexOf('/')) 
ncRNA_output_fasta = gtfDir + "/genes_ncRNA_all.fa"
"""
#!/usr/bin/env python 
import sys, os
import pandas as pd 
from collections import OrderedDict
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import twobitreader
import argparse


### specify annotation files here
gtf_outdir = "${gtfDir}"
gtfInFilePrefix = "${gtfInFilePrefix}"
twobitfile = "${params.genome2bit}"
genome = twobitreader.TwoBitFile(twobitfile) # create handler to open the 2bit file

### ncRNA gtf file
gtf_file_path = "${valid_gtf_path}"
gtf_outfile = "%s/genes_noncodingRNAs.gtf.gz" % (gtf_outdir)
fasta_outfile = "%s/genes_noncodingRNAs.fa" % (gtf_outdir)

### tRNA gtf file
tRNA_gtf_file_path = "${trna_gtf_path}"
tRNA_fasta_outfile = "%s/genes_tRNAs.fa" % (gtf_outdir)

### rRNA fasta file
rRNA_fasta_file = "${rrna_fa_path}"

### ncRNA output fasta
ncRNA_output_fasta = "${ncRNA_output_fasta}"

### functions

def read_in_gtf(infile_path, gtf_rows_to_skip):
	# read in a gtf file
	df = pd.read_csv(infile_path, sep='\t', dtype=str, header=None, skiprows=range(gtf_rows_to_skip))
	cols = ['#chrom', 'source', 'feature', 'chromStart', 'chromEnd', 'score', 'strand', 'frame', 'transcript_id']
	df.columns = cols
	return df


def parse_entry(tr):
	# given a single row of a gtf data.frame, parse column 8 into an ordered dictionary
	# this works for default gencode gtf file formats
	tr = tr.replace('"', '')
	trl = tr.split("; ")
	trdict = OrderedDict()

	for j in trl:
		k = j.split(" ")

		if k[0] in trdict:
			trdict[k[0]].append(k[1])
		else:       
			trdict[k[0]]=[k[1]]
	return trdict

def parse_mod_entry(tr):
	# given a single row of a gtf data.frame, parse column 8 into an ordered dictionary
	#	*** Use this to parse gtf files that have modified column 8
	#	*** this is for files compatible with downstream GFF package
	trl = tr.split(";")
	trdict = OrderedDict()

	for j in trl:
		k = j.split("=")

		if k[0] in trdict:
			trdict[k[0]].append(k[1])
		else:       
			trdict[k[0]]=[k[1]]
	return trdict


def build_gene_indexes(df):
	# Take an input dataframe and output an ordered dict with the index slices for every gene
	geneDict = OrderedDict()

	geneCount = 0
	previousGeneIndex = 0

	current_id=""
	current_gene=""

	for i in range(len(df)):

		if df.loc[i,'feature'] == 'gene':
			trdict = parse_entry(df.loc[i,'transcript_id'])

			curGeneID = trdict['gene_id'][0]
		
			if geneCount != 0:
				newGeneIndex = i
				geneDict[current_id] = [previousGeneIndex,newGeneIndex]
				previousGeneIndex = i
				current_id = trdict['gene_id'][0]
				geneCount += 1

			else:
				newgeneIndex = 0
				geneCount +=1
				current_id = trdict['gene_id'][0]
		if i == (len(df)-1):
			newGeneIndex = i+1
			current_id = trdict['gene_id'][0]
			geneDict[current_id] = [previousGeneIndex,newGeneIndex]
	return geneDict


def find_noncoding_genes(geneDict, df):
	# 3) filter geneDict on noncoding genes

	### tRNA's not part of the annotation, separate GTF file for these... 
	noncoding_types = ['Mt_rRNA', 
					   'Mt_tRNA',
					   'rRNA',
					   'miRNA',
					   'scRNA',
					   'scaRNA',
					   'snoRNA',
					   'snRNA',
					   'sRNA',
					   'vaultRNA'
					  ]

	geneDictNonCoding = OrderedDict() ### store noncoding transcripts only

	for gene in geneDict:
		tr = df.loc[geneDict[gene][0], 'transcript_id']

		trdict = parse_entry(tr)

		if trdict['gene_type'][0] in noncoding_types:
			geneDictNonCoding[gene] = geneDict[gene]

	return geneDictNonCoding


def build_transcript_indexes(geneDictCoding, df):
	# take an input geneDict and find the indexes for all transcripts associated with each gene

	TR_index_dict = OrderedDict()

	for gene in geneDictCoding:

		trDF = df.iloc[geneDictCoding[gene][0]:geneDictCoding[gene][1]]
	
		trPrev = -1
		trNamePrev = ""
		
		### iterate through a slice of the data frame for each gene
		### search for transcripts over that slice
		### find transcript slices
		for i in range(geneDictCoding[gene][0], geneDictCoding[gene][1]):
			if trDF.loc[i,'feature'] == 'transcript':
				trdict = parse_entry(trDF.loc[i,'transcript_id'])
				trCur = i
				trNameCur = trdict['transcript_id'][0]
				
				if trPrev != -1: # do not make an entry for the first transcript
					TR_index_dict[trNamePrev] = [trPrev, trCur]

				trPrev = trCur
				trNamePrev = trNameCur
			
			### for the final transcript
			if i == geneDictCoding[gene][1]-1:
				trdict = parse_entry(trDF.loc[i,'transcript_id'])
				TR_index_dict[trdict['transcript_id'][0]] = [trCur, i+1]
	return TR_index_dict


def convert_trsp_index(geneDictNonCoding, df, TR_index_dict):
	# take input geneDict and output a single transcript for each gene use transcript slices defined in TR_index_dict
	# for noncoding genes, take the longest transcripts
	# There are ~7 snoRNAs that have multiple noncoding transcripts, but the longest of these
	# always encompasses the entire shorter transcripts
	
	geneDictCanon = OrderedDict()
	
	for gene in geneDictNonCoding:
		trDF = df.iloc[geneDictNonCoding[gene][0]:geneDictNonCoding[gene][1]]
		trDFz = trDF.reset_index(drop=True)
		
		trCount = 0
		trDictLoc = OrderedDict()
		
		for i in range(len(trDFz)):
			if trDFz.loc[i, 'feature'] == 'transcript':
				tr = trDFz.loc[i, 'transcript_id']
				trdict = parse_entry(tr)
				trName = trdict['transcript_id'][0]
				trDictLoc[trName] = [trDFz.loc[i, 'chromStart'], trDFz.loc[i, 'chromEnd']]
				trCount += 1
		
		if trCount > 1:
			
			rangeDict = OrderedDict() ## store the ranges, and take the longest
			for key in trDictLoc:
				trRange = len(range(int(trDictLoc[key][0]),int(trDictLoc[key][1])))
				rangeDict[key] = trRange
				
			v=list(rangeDict.values())
			k=list(rangeDict.keys())
			trOut = k[v.index(max(v))]
			geneDictCanon[trOut] = [gene, TR_index_dict[trOut]]

		else: ## for genes with single transcripts
			trOut = trDictLoc.keys()[0]
			geneDictCanon[trOut] = [gene, TR_index_dict[trOut]]
	return geneDictCanon


def build_df_dict(df, geneDictCanon):
	# for each transcript selected in geneDictCanon, take a dataframe slice from df that has all entries for the transcript
	# currently keying the output dictionary on 'geneIDs', could be switched to transcript id's later if desired

	outDict = OrderedDict()

	for tr in geneDictCanon:
		outDict[geneDictCanon[tr][0]] = df.iloc[geneDictCanon[tr][1][0]:geneDictCanon[tr][1][1]]

	return outDict


def edit_col8(dfIn):
	# Edit the format of column 8 to remove spaces and add:
	# ID= ;Name= ; for transcript
	# Parent= ;Name= ; for other entires
	# Sort the entries based on the value of chromStart
	# This will be compaitble with downstream GTF parser from GFF package

	for i in dfIn.index:
		tr = parse_entry(dfIn.loc[i,'transcript_id'])

		### if line is a transcript, add Id= Name=
		if dfIn.loc[i,"feature"] == 'transcript':

			### make entry of new identifies to work with GFF parser 
			newline = "ID=%s;Name=%s" % ("".join(tr['transcript_id']), "".join(tr['gene_name']))
			### remove spaces and add semicolon to items in col 8
			line8 = ";".join("=".join((str(k),str(",".join(v)))) for k,v in tr.items())
			### Merge these into a single string
			outline = newline+";"+line8 
			outline = outline[0:-1] ## strip off last ';'

			### set the value of the cell to the edited line
			dfIn.at[i,"transcript_id"] = outline


		### for all other lines, add Parent= Name= 
		else:
			newline = "Parent=%s;Name=%s" % ("".join(tr['transcript_id']), "".join(tr['gene_name']))
			line8 = ";".join("=".join((str(k),str(",".join(v)))) for k,v in tr.items())
			outline = newline+";"+line8
			outline = outline[0:-1]
			dfIn.at[i,"transcript_id"] = outline

	dfIn['sort_vals'] = dfIn['chromStart'].astype(int)
	dfIn = dfIn.sort_values(by=['sort_vals'], axis=0, ascending=True)
	dfIn = dfIn.drop(['sort_vals'], axis=1)
	dfOut = dfIn.copy()

	return dfOut


def mod_last_column(outDictExclu):
	# Edit the final column of each dataframe
	
	outDictMod = OrderedDict()

	for key in outDictExclu:
		k = key
		dfIn = outDictExclu[key].copy()
		dfMod = edit_col8(dfIn)
		outDictMod[k] = dfMod
	
	return outDictMod


def output_df(outdict, out_file):
	# compile all modified dataframes for each transcript into a master dataframe
	# build the output dataframe from the modified dictionary and write this to a file:

	cols = ['#chrom', 'source', 'feature', 'chromStart', 'chromEnd', 'score', 'strand', 'frame', 'transcript_id']
	colOut = ['#chrom', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'transcript_id']
	gtfDF = pd.DataFrame(columns=cols)

	for trsp in outdict:
		gtfDF = gtfDF.append(outdict[trsp], ignore_index=True)
		
	gtfDF.columns = colOut
	# print gtfDF.head()
	gtfDF.to_csv(out_file, compression='gzip', sep='\t', index=False)


def tr_nc_dict(dfin):
	# similar to build_gene_indexes()
	# use this to define transcript indexes from GTF file written for ncRNA

	tr_nc_index_dict = OrderedDict()
	
	trCount = 0
	previousTrIndex = 0

	current_id=""
	current_tr=""

	for i in range(len(dfin)):

		if dfin.loc[i,'feature'] == 'transcript':
			trdict = parse_mod_entry(dfin.loc[i,'transcript_id'])

			curGeneID = trdict['gene_id'][0]

			if trCount != 0:
				newTrIndex = i
				tr_nc_index_dict[current_id] = [previousTrIndex,newTrIndex]
				previousTrIndex = i
				current_id = trdict['gene_id'][0]
				trCount += 1

			else:
				newTrIndex = 0
				trCount +=1
				current_id = trdict['gene_id'][0]
				
		if i == (len(dfin)-1):
			newTrIndex = i+1
			current_id = trdict['gene_id'][0]
			tr_nc_index_dict[current_id] = [previousTrIndex,newTrIndex]
			
	return tr_nc_index_dict


def check_multi_exon(tr_nc_index_dict, ncdf):
	# simple function that checks for multiple exons for each transcript in tr_nc_index_dict

	for gene in tr_nc_index_dict:
	
		tempdf = ncdf.iloc[tr_nc_index_dict[gene][0]:tr_nc_index_dict[gene][1]]
		exon_count = 0
		
		for i in tempdf.index:
			if tempdf.loc[i,'feature'] == 'exon':
				exon_count += 1
		if exon_count >1 :
			print " more than one exon for %s" % gene
			sys.exit()	# prevent writing fasta if there is multi exon transcript


def make_fasta_dict(ncdf):
	
	fasta_outdict = OrderedDict() 
	for i in ncdf.index:
		if ncdf.loc[i,'feature'] == 'transcript':
			chrom = ncdf.loc[i,'#chrom']
			chrStart = int(ncdf.loc[i,'chromStart'])
			chrEnd = int(ncdf.loc[i,'chromEnd'])
			strand = ncdf.loc[i,'strand']
			
			if strand == "+":
				chrStart = chrStart-1 ## gtf files are 1 based, convert to 0-based for python
				trSeq = SeqIO.Seq(genome[chrom][chrStart:chrEnd])
				trdict = parse_mod_entry(ncdf.loc[i,'transcript_id'])
			
			else: # for neg strand
				chrStart = chrStart-1
				trSeq = SeqIO.Seq(genome[chrom][chrStart:chrEnd])
				trSeq = trSeq.reverse_complement() # negative strand
				trdict = parse_mod_entry(ncdf.loc[i,'transcript_id'])

			### add output annotation line features
			trID = trdict['ID'][0]
			desc = "| "+trdict['gene_type'][0]+" | "+trdict['gene_name'][0]+ " | %s; %s; %s:%s" % (chrom, strand, chrStart, chrEnd)

			trSeqRec = SeqRecord(trSeq, id=trID, name=trdict['gene_name'][0], description=desc)
			fasta_outdict[trID] = trSeqRec
	
	return fasta_outdict


def write_output_fasta_ncRNA(fasta_outdict, fasta_outfile):
	# takes a dictionary of all fasta entries and write an output fasta file

	ncFasta_iterator = (record for record in fasta_outdict.values())
	SeqIO.write(ncFasta_iterator, fasta_outfile, "fasta")


### tRNA database specific functions:

def make_tRNA_fasta_dict(tRNAdf):
	# similar to make_fasta_dict, but for the tRNA database

	tRNA_fasta_outdict = OrderedDict()

	for i in tRNAdf.index:

		if tRNAdf.loc[i,'feature'] == 'tRNA':
			chrom = tRNAdf.loc[i,'#chrom']
			chrStart = int(tRNAdf.loc[i,'chromStart'])
			chrEnd = int(tRNAdf.loc[i,'chromEnd'])
			strand = tRNAdf.loc[i,'strand']
			
			if strand == "+":
				chrStart = chrStart-1 ### gtf files are 1-based, convert to 0-based
				trSeq = SeqIO.Seq(genome[chrom][chrStart:chrEnd])
				trdict = parse_entry(tRNAdf.loc[i,'transcript_id'])
			
			else: # for neg strand
				chrStart = chrStart-1
				trSeq = SeqIO.Seq(genome[chrom][chrStart:chrEnd])
				trSeq = trSeq.reverse_complement()
				trdict = parse_entry(tRNAdf.loc[i,'transcript_id'])

			trID = "tRNA_"+trdict['gene_id'][0]
			desc = "| tRNA | "+trdict['gene_type'][0] + " | %s; %s; %s:%s" % (chrom, strand, chrStart, chrEnd)

			trSeqRec = SeqRecord(trSeq, id=trID, name=trdict['gene_name'][0], description=desc)
			tRNA_fasta_outdict[trID] = trSeqRec
	
	return tRNA_fasta_outdict
		
def write_output_fasta_tRNA(tRNA_fasta_outdict, tRNA_fasta_outfile):

	tRNA_ncFasta_iterator = (record for record in tRNA_fasta_outdict.values())
	SeqIO.write(tRNA_ncFasta_iterator, tRNA_fasta_outfile, "fasta")

def merge_ncRNA_fastas():
	# command to merge all fasta files into a single fasta
	merge_cmnd = "cat %s %s %s > %s" % (rRNA_fasta_file, fasta_outfile, tRNA_fasta_outfile, ncRNA_output_fasta)
	print merge_cmnd
	os.system(merge_cmnd)

	index_cmnd = "samtools faidx %s" % (ncRNA_output_fasta)
	print index_cmnd
	os.system(index_cmnd)

def main():
	# if not os.path.isfile(outfilestring+".csv"):
	### make the genecode ncRNA gtf file
	df = read_in_gtf(gtf_file_path, gtf_rows_to_skip=5)
	geneDict = build_gene_indexes(df)
	geneDictNonCoding = find_noncoding_genes(geneDict, df)
	TR_index_dict = build_transcript_indexes(geneDictNonCoding, df)
	geneDictCanon = convert_trsp_index(geneDictNonCoding, df, TR_index_dict)
	outDict = build_df_dict(df, geneDictCanon)
	outDictMod = mod_last_column(outDictExclu=outDict)
	output_df(outDictMod, gtf_outfile)

	### make fasta file
	ncdf = read_in_gtf(gtf_outfile, gtf_rows_to_skip=1)
	print ncdf.head()
	tr_nc_index_dict = tr_nc_dict(ncdf)
	check_multi_exon(tr_nc_index_dict, ncdf)
	fasta_outdict = make_fasta_dict(ncdf)
	write_output_fasta_ncRNA(fasta_outdict, fasta_outfile)

	### for tRNA fasta file
	tRNAdf = read_in_gtf(tRNA_gtf_file_path, gtf_rows_to_skip=5)
	tRNA_fasta_outdict = make_tRNA_fasta_dict(tRNAdf)
	write_output_fasta_tRNA(tRNA_fasta_outdict, tRNA_fasta_outfile)

	### write the output file:
	merge_ncRNA_fastas()

	

if __name__ == '__main__':
	main()
"""
}

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 1
    $MEMORY = 48
}
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 700
    $CPU  = 1
    $MEMORY = 48
    $QUEUE = "long"
} 
//* platform
//* autofill

process Ribosome_Profiling_Build_Index_build_Star_Index {

input:
 val ncRNAFastaPath from g40_5_ncRNAFastaPath_g40_6

output:
 val ncRNAFastaPath  into g40_6_ncRNAFastaPath_g40_7

script:
gtfInFilePrefix = ncRNAFastaPath.substring(0, ncRNAFastaPath.lastIndexOf('.')) 
indexDir = ncRNAFastaPath.substring(0, ncRNAFastaPath.lastIndexOf('/')) 
"""
#!/usr/bin/env python 
import sys, os

def main():
	
	### ncRNA
	ncSparsity = 1
	ncGenomeDir = "${indexDir}" + "/STARIndex_ncRNA"
	ncGenomeFasta = "${ncRNAFastaPath}"
	ncSAindexNbases = 9
	threadNumb = 4
	
	if not os.path.exists(ncGenomeDir):
		os.makedirs(ncGenomeDir)
		CMMD= 'STAR  --runThreadN %s  --runMode genomeGenerate  --genomeSAsparseD %s  --genomeDir %s  --genomeFastaFiles %s  --genomeSAindexNbases %s  ' % (threadNumb, ncSparsity, ncGenomeDir, ncGenomeFasta, ncSAindexNbases)
		print CMMD
		if os.system(CMMD) != 0:
			raise Exception(CMMD + " failed.")
			
	hgSparsity = 1
	hgGenomeDir = "${indexDir}" + "/STARIndex"
	hgGenomeFasta = "${params.genome}"
	hgSjdbGTF = "${params.gtf}"
	hgSAindexNbases = 14
	
	if not os.path.exists(hgGenomeDir):
		os.makedirs(hgGenomeDir)
		CMMD= 'STAR  --runThreadN %s  --runMode genomeGenerate  --genomeSAsparseD %s  --genomeDir %s  --genomeFastaFiles %s  --genomeSAindexNbases %s --sjdbGTFfile %s' % (threadNumb, hgSparsity, hgGenomeDir, hgGenomeFasta, hgSAindexNbases, hgSjdbGTF)
		print CMMD
		if os.system(CMMD) != 0:
			raise Exception(CMMD + " failed.")
	

if __name__ == '__main__':
	main()
"""
}

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 1
    $MEMORY = 20
}
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 1000
    $CPU  = 1
    $MEMORY = 20
    $QUEUE = "long"
} 
//* platform
//* autofill

process Ribosome_Profiling_Build_Index_gtfToAllTr {

input:
 val gtf from g40_0_gtfPath_g40_3

output:
 val gtf_outfile  into g40_3_valid_gtf_g40_4

script:
gtf_outdir  = params.gtf.substring(0, params.gtf.lastIndexOf('/')) 
gtf_outfile = gtf_outdir+"/genes_all_tr.gtf"
"""
#!/usr/bin/env python 
import sys, os
import pandas as pd 
import time
from collections import OrderedDict
import argparse

gtf_file_path = "${params.gtf}" 
gtf_rows_to_skip = 5

gtf_outdir = "${gtf_outdir}"
gtf_outfile = "${gtf_outfile}"

def read_in_gtf():
	# 1) read in the gtf file here
	df = pd.read_csv(gtf_file_path, sep='\t', dtype=str, header=None, skiprows=range(gtf_rows_to_skip))
	cols = ['#chrom', 'source', 'feature', 'chromStart', 'chromEnd', 'score', 'strand', 'frame', 'transcript_id']
	df.columns = cols
	return df


def parse_entry(tr):
	# given a single row of a gtf data.frame, parse column 8 into an ordered dictionary
	tr = tr.replace('"', '')
	trl = tr.split("; ")
	trdict = OrderedDict()

	for j in trl:
		k = j.split(" ")

		if k[0] in trdict:
			trdict[k[0]].append(k[1])
		else:       
			trdict[k[0]]=[k[1]]
	return trdict


def build_gene_indexes(df):
	# 2) take an input dataframe and output an ordered dict with the index slices for every gene
	geneDict = OrderedDict()

	geneCount = 0
	previousGeneIndex = 0

	current_id=""
	current_gene=""

	for i in range(len(df)):

		if df.loc[i,'feature'] == 'gene':
			trdict = parse_entry(df.loc[i,'transcript_id'])

			curGeneID = trdict['gene_id'][0]
		
			if geneCount != 0:
				newGeneIndex = i
				geneDict[current_id] = [previousGeneIndex,newGeneIndex]
				previousGeneIndex = i
				current_id = trdict['gene_id'][0]
				geneCount += 1

			else:
				newgeneIndex = 0
				geneCount +=1
				current_id = trdict['gene_id'][0]
		if i == (len(df)-1):
			newGeneIndex = i+1
			current_id = trdict['gene_id'][0]
			geneDict[current_id] = [previousGeneIndex,newGeneIndex]
	return geneDict

### -------- ###
			
def find_valid_genes(geneDict, df):
	# 3) filter geneDict on coding genes only

	geneDictOut = OrderedDict() ### store coding transcripts only

	geneTypeList = []

	valid_gene_types = [
						'protein_coding', 
						'unitary_pseudogene',
						'lincRNA',
						'translated_processed_pseudogene',
						'sense_intronic',
						'pseudogene',
						'polymorphic_pseudogene',
						'bidirectional_promoter_lncRNA',
						'sense_overlapping',
						'antisense',
						'unprocessed_pseudogene',
						'transcribed_unprocessed_pseudogene',
						'3prime_overlapping_ncRNA',
						'non_coding',
						'processed_transcript',
						'transcribed_processed_pseudogene',
						'processed_pseudogene',
						'transcribed_unitary_pseudogene'
						]

	for gene in geneDict:
		tr = df.loc[geneDict[gene][0], 'transcript_id']
		
		trdict = parse_entry(tr)

		
		if trdict['gene_type'][0] in valid_gene_types:
			geneDictOut[gene] = geneDict[gene]
			geneTypeList.append(trdict['gene_type'][0])

	print set(geneTypeList)
	return geneDictOut



def build_transcript_indexes(geneDictCoding, df):
	# take an input geneDict and find the indexes for all transcripts associated with each gene

	TR_index_dict = OrderedDict()

	for gene in geneDictCoding:

		trDF = df.iloc[geneDictCoding[gene][0]:geneDictCoding[gene][1]]
	
		trPrev = -1
		trNamePrev = ""
		
		### iterate through a slice of the data frame for each gene
		### search for transcripts ofver that slice
		### find transcript slices
		for i in range(geneDictCoding[gene][0], geneDictCoding[gene][1]):
			if trDF.loc[i,'feature'] == 'transcript':
				trdict = parse_entry(trDF.loc[i,'transcript_id'])
				trCur = i
				trNameCur = trdict['transcript_id'][0]
				
				if trPrev != -1: # do not make an entry for the first transcript
					TR_index_dict[trNamePrev] = [trPrev, trCur]

				trPrev = trCur
				trNamePrev = trNameCur
			
			### for the final transcript
			if i == geneDictCoding[gene][1]-1:
				trdict = parse_entry(trDF.loc[i,'transcript_id'])
				TR_index_dict[trdict['transcript_id'][0]] = [trCur, i+1]
	return TR_index_dict



def find_all_valid_transcripts(geneDictCoding, df, TR_index_dict):
	# take gene and transcript dictionaries, and choose the best transcript accoding these criteria:
	# re-working this to include all protein coding, and all ccds genes

	geneDictCanon = OrderedDict()
	geneDictChrLoc = OrderedDict() # [trxStart, trxEnd, tr_for_start, tr_for_stop]

	ap1 = 0
	ap2 = 0
	ap3 = 0
	ap4 = 0
	ap5 = 0

	ap1_alt = 0
	ap2_alt = 0

	single_iso_pri = 0
	single_iso_alt = 0
	multi_iso = 0
	total_iso = 0
	total_iso_pri = 0
	total_iso_alt = 0

	identical_trsps = 0
	noAppris = 0
	
	noCCDS = 0
	singleCCDS = 0
	multiCCDS = 0

	noProtCode = 0
	singleProtCode = 0
	multiProtCode = 0

	protCodeCount = 0
	cds_NF_count = 0
	noValidUTRsRemaining = 0
	invalidUTRtotal = 0

	### now we have all genes that are coding, with start and stop indexes:
	for gene in geneDictCoding:


		trDF = df.iloc[geneDictCoding[gene][0]:geneDictCoding[gene][1]]
		trDFz = trDF.reset_index(drop=True) # z is for zero based conversion here   

		### Now we are finding the CCDS transcript only
		### and retrieving the 'best' transcript according to appris_principal
		
		protCodeTrDict = OrderedDict()
		ccdsTrDict = OrderedDict() ### store all transcripts that have CCDS annotation
		trCount = 0
		# protCodeCount = 0
		ccdsCount = 0
		
		### find all protein coding transcripts
		for i in range(len(trDFz)):
			if trDFz.loc[i, 'feature'] == 'transcript':
				tr = trDFz.loc[i, 'transcript_id']
				trdict = parse_entry(tr)
				trOut = trdict['transcript_id'][0]
				geneDictCanon[trOut] = [gene, TR_index_dict[trOut]] 
	print "transcripts in geneDictCanon == %s" % (len(geneDictCanon))

	return geneDictCanon

def build_df_dict_byTr(df, geneDictCanon):
	# for each transcript selected in geneDictCanon, take a dataframe slice from df that has all entries for the transcript
	# currently keying the output dictionary on 'geneIDs', could be switched to transcript id's later if desired

	outDict = OrderedDict()

	for tr in geneDictCanon:
		outDict[tr] = df.iloc[geneDictCanon[tr][1][0]:geneDictCanon[tr][1][1]]

	return outDict



def build_df_dict(df, geneDictCanon):
	# for each transcript selected in geneDictCanon, take a dataframe slice from df that has all entries for the transcript
	# currently keying the output dictionary on 'geneIDs', could be switched to transcript id's later if desired

	outDict = OrderedDict()

	for tr in geneDictCanon:
		outDict[geneDictCanon[tr][0]] = df.iloc[geneDictCanon[tr][1][0]:geneDictCanon[tr][1][1]]

	return outDict

def overlap_features(genedf):
	# for a given data.frame of a transcript, extract info to check for overlaps Used in find_transcript_overlaps()
	# output: [chrStart, chrEnd, strand, geneName, chrom]
	chrStart = genedf.loc[0,'chromStart']
	chrEnd = genedf.loc[0, 'chromEnd']
	strand = genedf.loc[0, 'strand']
	trdict = parse_entry(genedf.loc[0, 'transcript_id'])
	geneName = trdict['gene_id'][0]
	trName = trdict['transcript_id'][0]
	chrom = genedf.loc[0, '#chrom']
	outlist = [chrStart, chrEnd, strand, geneName, chrom, trName]
	return outlist

def find_transcript_overlaps(outDict, geneDictChrLoc):
	# For each transcript, look at nearest neighbor on same strand and check if transcripts have overlaps
	# add all transcripts with overlaps to exclusion_overlaps dictionary
	
	exclusion_overlaps = OrderedDict()
	excluded_trsp_count = 0
	gene_key_list = outDict.keys() ### ordered list of gene names in the dictionary

	trsp_entry = -1
	for gene in outDict:

		trsp_entry += 1
		genedf = outDict[gene]
		genedf = genedf.reset_index(drop=True)
		
		# chrStart = genedf.loc[0,'chromStart']
		# chrEnd = genedf.loc[0, 'chromEnd']
		# strand = genedf.loc[0, 'strand']
		# trdict = parse_entry(genedf.loc[0, 'transcript_id'])
		# geneName = trdict['gene_id'][0]
		
		overlap_feats = overlap_features(genedf)
		
		cur_strand = overlap_feats[2]
		cur_chrom = overlap_feats[4]
		
		##### DOWNSTREAM Overlaps 
		next_strand = "0"
		search_index = 1
		break_signal = 0
		
		while next_strand != cur_strand:
		
			next_tr = trsp_entry+search_index ### only for the end of the list
			if next_tr == len(gene_key_list):
				print "reached end of downstream list for %s " % gene
				break_signal = 1
				break
			
			next_gene = gene_key_list[next_tr] ### adding this
			nextdf = outDict[gene_key_list[next_tr]]
			nextdf = nextdf.reset_index(drop=True)
			next_over = overlap_features(nextdf)

			### adding genome position
			
			next_chrom = next_over[4]
			
			### check to make sure that this is still the same chr
			if next_chrom != cur_chrom:
				print "search downstream to next chr for %s" % gene
				break_signal = 1
				break

			next_strand = next_over[2]
			search_index +=1
			
		next_gene_start = geneDictChrLoc[next_gene][0]

		# if int(overlap_feats[1]) > int(next_over[0]) and break_signal !=1:
		if int(overlap_feats[1]) > int(next_gene_start) and break_signal !=1:
			# print "******"
			# print "OVERLAPING DOWNSTREAM TRSP %s, %s" % (gene, overlap_feats[5])
			# print "******"
			# print overlap_feats
			# # print next_over
			# print next_gene

			exclusion_overlaps[gene] = overlap_feats
			excluded_trsp_count += 1
			
			
		##### UPSTREAM Overlaps
		next_strand = "0"
		search_index = 1
		break_signal = 0
		
		while next_strand != cur_strand:
			next_tr = trsp_entry-search_index ### this time subtract search index
			if next_tr == -1:
				print "end of upstream list for %s" % gene
				break_signal = 1
				break
				
			next_gene = gene_key_list[next_tr] ### adding this
			nextdf = outDict[gene_key_list[next_tr]]
			nextdf = nextdf.reset_index(drop=True)
			next_over = overlap_features(nextdf)
			
			next_chrom = next_over[4]

			### check to make sure that this is still the same chr
			if next_chrom != cur_chrom:
				print "search upstream to prev chr for %s" % gene
				break_signal = 1
				break

			next_strand = next_over[2]
			search_index +=1
			
		next_gene_end = geneDictChrLoc[next_gene][1]

		# if int(overlap_feats[0]) < int(next_over[1]) and break_signal !=1:
		if int(overlap_feats[0]) < int(next_gene_end) and break_signal !=1:
			# print "******"
			# print "OVERLAPING UPSTREAM TRSP %s, %s" % (gene, overlap_feats[5])
			# print "******"
			# print overlap_feats
			# # print next_over
			# print next_gene
			
			exclusion_overlaps[gene] = overlap_feats
			excluded_trsp_count += 1

	print "Total Overlapping transcripts revised: %s" % excluded_trsp_count
	return exclusion_overlaps


def remove_overlapping_transcripts(outDict, exclusion_overlaps):
	# simply remove all of the transcripts in exclusion_overlaps from outDict

	outDictExclu = outDict

	for exclu in exclusion_overlaps:
	#     print exclu
		outDictExclu.pop(exclu)
		
	print len(outDictExclu)
	return outDictExclu




def edit_col8(dfIn):
	# Edit the format of column 8 to remove spaces and add:
	# ID= Name= for transcript
	# Parent= Name= for other entires
	# Sort the entries based on the value of chromStart
	# This will be compaitble with downstream GTF parser
	
	for i in dfIn.index:
		tr = parse_entry(dfIn.loc[i,'transcript_id'])

		### if line is a transcript, add Id= Name=
		if dfIn.loc[i,"feature"] == 'transcript':

			### make entry of new identifies to work with GFF parser 
			newline = "ID=%s;Name=%s" % ("".join(tr['transcript_id']), "".join(tr['gene_name']))
			### remove spaces and add semicolon to items in col 8
			line8 = ";".join("=".join((str(k),str(",".join(v)))) for k,v in tr.items())
			### Merge these into a single string
			outline = newline+";"+line8 
			outline = outline[0:-1] ## strip off last ';'

			### set the value of the cell to the edited line
			dfIn.at[i,"transcript_id"] = outline


		### for all other lines, add Parent= Name= 
		else:
			newline = "Parent=%s;Name=%s" % ("".join(tr['transcript_id']), "".join(tr['gene_name']))
			line8 = ";".join("=".join((str(k),str(",".join(v)))) for k,v in tr.items())
			outline = newline+";"+line8
			outline = outline[0:-1]
			dfIn.at[i,"transcript_id"] = outline
			
	dfIn['sort_vals'] = dfIn['chromStart'].astype(int)
	dfIn = dfIn.sort_values(by=['sort_vals'], axis=0, ascending=True)
	dfIn = dfIn.drop(['sort_vals'], axis=1)
	dfOut = dfIn.copy()
	
	return dfOut


def define_pseudogene_positions(pseudoGeneDict, df):

	pseudoChromDict = OrderedDict()

	for gene in pseudoGeneDict:

		# print gene
		tempdf = df.iloc[pseudoGeneDict[gene][0]:pseudoGeneDict[gene][1]]
		tempdf = tempdf.loc[tempdf['feature']=='transcript']
		tempdf = tempdf.reset_index(drop=True)

		# print tempdf
		overlap_feats = overlap_features(tempdf)
		# print overlap_feats

		over_range = overlap_feats[0:2]
		over_range = [int(x) for x in over_range]
		# over_range = map(int, over_range)
		# cur_chrom = overlap_feats[4]
		cur_chrom = overlap_feats[4]+"_"+overlap_feats[2] ## account for strandedness
		# print cur_chrom

		if cur_chrom in pseudoChromDict:
			pseudoChromDict[cur_chrom].append(over_range)
		else:
			pseudoChromDict[cur_chrom] = [over_range]

			# pseudoChromDict[cur_chrom] = pseudoChromDict[cur_chrom].append(over_range)

	# print pseudoChromDict
	return pseudoChromDict

	# print pseudoChromDict


def check_pseudo_overlap(x1, x2, y1, y2):
	# check for overlaps [x1, x2] and [y1, y2]
	# should only be true if these overlap
	return max(x1, y1) <= min(x2, y2) 


def find_pseudo_overlaps(outDictExclu, pseudoChromDict):

	pseudo_exclude = OrderedDict()
	pseudo_exclude_count = 0

	for gene in outDictExclu:
		# print gene
		# print outDictExclu[gene]

		genedf = outDictExclu[gene]
		genedf = genedf.reset_index(drop=True)

		for feat in genedf.index: ### only check UTR entries, not perfect, but can be modified, maybe to exons?
			if genedf.loc[feat, 'feature'] == 'UTR':
				chrStart = genedf.loc[feat,'chromStart']
				chrEnd = genedf.loc[feat, 'chromEnd']
				strand = genedf.loc[feat, 'strand']
				trdict = parse_entry(genedf.loc[feat, 'transcript_id'])
				geneName = trdict['gene_id'][0]
				trName = trdict['transcript_id'][0]
				chrom = genedf.loc[feat, '#chrom']
				outlist = [chrStart, chrEnd, strand, geneName, chrom, trName]

				overlap_feats = outlist
				# overlap_feats = overlap_features(genedf)
				
				gene_over_range = overlap_feats[0:2]
				gene_over_range = [int(x) for x in gene_over_range]
				# gene_over_range = map(int, gene_over_range)
				cur_chrom = overlap_feats[4]+"_"+overlap_feats[2]

				# print gene_over_range
				# print gene_over_range[0]
				# print gene_over_range[1]

				try: # 'chrM_+' is not in pseudoChromDict ... 
					chr_ref_list = pseudoChromDict[cur_chrom]
				except KeyError as e:
					# pseudoChromDict[cur_chrom]
					print(e.message)
					print "no key in pseudoChromDict"
					continue

				for pseudo in chr_ref_list:
					# print pseudo
					overTest = check_pseudo_overlap(gene_over_range[0], gene_over_range[1], pseudo[0], pseudo[1])

					if overTest == True:
						print "overlap found here for %s with pseudo %s" % (gene, pseudo)
						pseudo_exclude[gene] = [cur_chrom, pseudo]
						pseudo_exclude_count += 1

	print pseudo_exclude_count, "excluded pseudogene count"

	return pseudo_exclude


def mod_last_column(outDictExclu):
	# Edit the final column of each dataframe
	
	outDictMod = OrderedDict()

	for key in outDictExclu:
		k = key
		dfIn = outDictExclu[key].copy()
		dfMod = edit_col8(dfIn)
		outDictMod[k] = dfMod
	
	return outDictMod

def output_df(outdict, out_file):
	# compile all modified dataframes for each transcript into a master dataframe
	# build the output dataframe from the modified dictionary and write this to a file:
	
	cols = ['#chrom', 'source', 'feature', 'chromStart', 'chromEnd', 'score', 'strand', 'frame', 'transcript_id']
	colOut = ['#chrom', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'transcript_id']
	gtfDF = pd.DataFrame(columns=cols)

	for trsp in outdict:
		gtfDF = gtfDF.append(outdict[trsp], ignore_index=True)
		
	gtfDF.columns = colOut
	print gtfDF.head(5)
	gtfDF.to_csv(out_file, sep='\t', index=False)


def main():
	if not os.path.isfile(gtf_outfile):
		time_start = time.time()
		df = read_in_gtf()
		print df.head()
		geneDict = build_gene_indexes(df)
		geneDictCoding = find_valid_genes(geneDict, df)
		TR_index_dict = build_transcript_indexes(geneDictCoding, df)
		geneDictCanon= find_all_valid_transcripts(geneDictCoding, df, TR_index_dict)
		outDict = build_df_dict_byTr(df, geneDictCanon)
		outDictMod = mod_last_column(outDict)
		output_df(outDictMod, gtf_outfile)
		time_end = time.time()
		print "total time == ", (time_end - time_start)

if __name__ == '__main__':
	main()
"""
}

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 10
    $MEMORY = 30
}
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 500
    $CPU  = 10
    $MEMORY = 30
    $QUEUE = "long"
} 
//* platform
//* autofill

process Ribosome_Profiling_Build_Index_build_annotation_allTR {

input:
 val valid_gtf_path from g40_3_valid_gtf_g40_4

output:
 val valid_gtf_path  into g40_4_valid_gtf_g40_7

script:
gtfInFilePrefix = valid_gtf_path.substring(0, valid_gtf_path.lastIndexOf('.')) 
gtfDir = valid_gtf_path.substring(0, valid_gtf_path.lastIndexOf('/')) 
"""
#!/usr/bin/env python 
import sys, os
import GFF
import twobitreader
import argparse
from Bio import SeqIO
import csv
from Bio.Seq import Seq
import pandas as pd
from collections import OrderedDict
from pathos.multiprocessing import ProcessingPool as Pool


# Add list of acceptable chromosomes that will be output to the table
validChrs = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 
			'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15',
			'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 
			'chrX', 'chrY', 'chrM', 'chrSinV', 'chrLUC','egfp', 'mcherry']


### specify annotation files here
GenAnPath = "${gtfDir}"
GTFfile = "${valid_gtf_path}"
gtfInFilePrefix = "${gtfInFilePrefix}"
twobitfile = "${params.genome2bit}"
genome = twobitreader.TwoBitFile(twobitfile) # create handler to open the 2bit file


### functions imported from Colin's densbuilder
def makeGFFlist(GTFinput):
	# Create a dictionary with a key for each chromosome in the GFF file
	GTFlist={}
	for chr in GTFinput:
		GTFlist[chr.id]=chr
	return GTFlist

def chrpostomrnapos(chrpos,chrom,featnum,GFFlist):
	# This funciton takes a genomic query position (chrpos) and a chromosome number (ex 'chr6')
	#	along with the feature number (the entry in the gff file for that transcript) defined by build_utr_table()
	#	and the dictionary of transcript from the GFFlist
	# The output is mrnapos which is the transcript relative position (position along the mRNA)
	#	of the original genomic query postion (chrpos)
	
	trsp_id= GFFlist[chrom].features[featnum].id
	trsp_strand= GFFlist[chrom].features[featnum].strand
	trsp_chromstart= int(GFFlist[chrom].features[featnum].location.start.position)  # 0-based
	trsp_chromend= int(GFFlist[chrom].features[featnum].location.end.position)
	sublist=[]

	for subfeature in GFFlist[chrom].features[featnum].sub_features:     # Make list of features
		if subfeature.type== 'exon':
			start= subfeature.location.start.position
			end= subfeature.location.end.position
			sublist.append([start,end])

	if trsp_strand== -1:    
		sublist.reverse()
	assert len(sublist)!= 0, ("transcript %s has a sublist length of zero!" % trsp_id)

	prevexonlen= 0 
	for item in sublist:
		exonstart= item[0]
		exonend= item[1]
		exonlen= exonend- exonstart

		if trsp_strand== 1:
			if chrpos>= exonstart and chrpos< exonend:      
				mrnapos= prevexonlen+ chrpos- exonstart
				return mrnapos
			else:   prevexonlen+= exonlen
		else:
			if chrpos< exonend and chrpos>= exonstart:
				mrnapos= prevexonlen+ (exonend-1)- chrpos       # Need -1 because end is not actual end, it is 1 beyond end.
				return mrnapos
			else:   prevexonlen+= exonlen 

def build_utr_table(GFFlist, include_noncanon_start, include_noncanon_stop):
	# This is a function to get the cds and utr sizes for an mRNA from a GFF file
	# returns a list with: #transcript,chrom,featnum,strand,mrna_len,cds_len,5utr_len,3utr_len,gene_name
	# Includes most of the functions from densebuilder_main but does not return counts
	# GFFlist = GFFinput

	transcriptdict={}
	ucscIDlist = []
	total_transcripts = 0 
	nonvalidchorms = 0
	nonATGstart = 0
	wrongstopcodon = 0
	validchroms = 0
	excluded_chroms = []
	included_chroms = []
	for chrom in GFFlist:
		if not chrom in validChrs:
			excluded_chroms.append(chrom)
			nonvalidchorms += 1
			print chrom
			continue	# check that only valid choromosomes are used
		validchroms+=1
		included_chroms.append(chrom)
		transcriptnum= -1 # set to negative one so first transcript is == to 0
		for transcript in GFFlist[chrom].features:	# this is where the SeqFeatures are actually stored
			tr_attribute_list = []
			transcriptnum+=1
			trsp_id= transcript.id # it is a number 
			trsp_strand= transcript.strand
			### changing this to be compatible with new hg38 annotation
			# print transcript.qualifiers ### these are all of the fields parsed by the GTF parser from column 8, output is a dictionary {'key':['item1', 'item2', 'ect']}
			trsp_genename= transcript.qualifiers['Name'][0]
			trsp_chromstart= int(transcript.location.start.position)  # 0-based
			trsp_chromend= int(transcript.location.end.position) 
			transcriptlist= [0.0 for x in range(abs(trsp_chromend- trsp_chromstart))] # a list for transcript (pre-mRNA), not CDS
			
			exonsplicedseq= SeqIO.Seq('')
			transcriptseq= SeqIO.Seq(genome[chrom][trsp_chromstart: trsp_chromend])
			
			### use lists to handle transcripts with multiple start and stop codons
			startCodonMrnaList = []
			stopCodonMrnaList = []

			for item in GFFlist[chrom].features[transcriptnum].sub_features:
				


				if trsp_strand== 1:

					### dealing with transcripts having multiple start or stop codon entries, if spaning splice junctions


					if item.type== 'exon': # or item.type== 'CDS':	# For yeast, use 'CDS'
						exonstart= int(item.location.start.position)  # 0-based position
						exonend= int(item.location.end.position) # not 0-based
						exonstart_feat= exonstart- trsp_chromstart
						exonend_feat= exonend- trsp_chromstart # Not 0-based, it is fine for length....next line. 
						exonsplicedseq+= transcriptseq[exonstart_feat:exonend_feat] # takes from exonstart to exonend-1
					if item.type== 'start_codon':
						startcodonpos= item.location.start.position # 0-based position
						# startcodonmrnapos=  chrpostomrnapos(startcodonpos,chrom,transcriptnum,GFFlist)	# spliced mRNA position
						startCodonMrnaList.append(chrpostomrnapos(startcodonpos,chrom,transcriptnum,GFFlist))	# spliced mRNA position
						# print startcodonmrnapos
					if item.type== 'stop_codon':
						stopcodonpos= item.location.end.position- 1 # 0-based position
						# stopcodonmrnapos= chrpostomrnapos(stopcodonpos,chrom,transcriptnum,GFFlist)
						stopCodonMrnaList.append(chrpostomrnapos(stopcodonpos,chrom,transcriptnum,GFFlist))
						# print stopcodonmrnapos

				if trsp_strand== -1:
					print 'neg_strand'
					# reverse_complement() # this comes from seqIO
					transcriptseq_rev= transcriptseq.reverse_complement() 

					if item.type== 'exon': # or item.type== 'CDS':	# For yeast, use 'CDS'
						exonstart= int(item.location.start.position)  # 0-based position
						exonend= int(item.location.end.position)	# not 0-based 
						exonstart_feat= (trsp_chromend-1)- (exonend- 1) 		# 0-based
						exonend_feat= (trsp_chromend-1)- exonstart 		# 0-based
						exonseq= transcriptseq_rev[exonstart_feat:exonend_feat+ 1] 
						exonsplicedseq= exonseq+ exonsplicedseq
					if item.type== 'start_codon':
						startcodonpos= item.location.end.position- 1 # Need to -1 to be 0-based.
						print startcodonpos
						# startcodonmrnapos= chrpostomrnapos(startcodonpos,chrom,transcriptnum,GFFlist)
						startCodonMrnaList.append(chrpostomrnapos(startcodonpos,chrom,transcriptnum,GFFlist))
						# print "start codon: ", startcodonmrnapos
					if item.type== 'stop_codon':
						stopcodonpos= item.location.start.position	# start.position is 0-based already. 
						print stopcodonpos
						# stopcodonmrnapos= chrpostomrnapos(stopcodonpos,chrom,transcriptnum,GFFlist)
						stopCodonMrnaList.append(chrpostomrnapos(stopcodonpos,chrom,transcriptnum,GFFlist))
						# print "stop codon: ", stopcodonmrnapos
			
			if len(startCodonMrnaList) > 0:
				print "MORE THAN 1 START", startCodonMrnaList
				startcodonmrnapos = min(startCodonMrnaList)
			else:
				print "!!! no start codon for %s" % (trsp_id)
				startcodonmrnapos = 0 ### adding for transcripts without start codon
			# if len(stopCodonMrnaList)
			if len(stopCodonMrnaList) > 0:
				stopcodonmrnapos = max(stopCodonMrnaList)
			else:
				print "!!! no stop codon for %s" % (trsp_id)
				stopcodonmrnapos = len(exonsplicedseq) - 3 ### leave 3nt's in "3'UTR"



			cdsseq= exonsplicedseq[startcodonmrnapos: stopcodonmrnapos+ 1] # take from startcodonmrnapos to stopcodonmrnapos
			utr5seq = exonsplicedseq[:startcodonmrnapos]
			utr3seq = exonsplicedseq[stopcodonmrnapos+1:]
			
			print trsp_id
			print transcript.qualifiers['transcript_name']
			print trsp_strand
			print utr5seq
			print " - - - "
			print cdsseq
			print " - - - "
			print utr3seq
			print utr5seq+cdsseq+utr3seq
			print ""
			print transcriptseq

			
			if include_noncanon_start == False:
				if str(cdsseq[:3].upper())!= "ATG":	
					nonATGstart += 1
					print "non canon start"
					print trsp_id
					print cdsseq
					print ""
					continue	# ignore non-AUG start codons

			stopcodon= str(cdsseq[-3:].upper())
			if len(utr3seq) > 0:
				stop4nt = stopcodon +str(utr3seq[0].upper())
			elif len(utr3seq) == 0: 
				stop4nt = '0'
			else:
				print "there is a 3'UTR with negative length..."
				sys.exit()
			

			if include_noncanon_stop == False:
				if stopcodon!= "TGA" and stopcodon!= "TAG" and stopcodon!= "TAA":	
					wrongstopcodon += 1
					print "wrong stop!"
					print trsp_id
					print cdsseq
					print ""
					continue	# ignore weird stop codons

			# build itmes in transcript attribute list
			mRNAlen = len(exonsplicedseq)
			cdslen = len(cdsseq)
			utr5len = len(utr5seq)
			utr3len = len(utr3seq)
			assert mRNAlen == utr3len+cdslen+utr5len # check that sum of features equals mRNA length

			trsp_attr_list = [trsp_id, chrom, transcriptnum, trsp_strand, mRNAlen, cdslen, utr5len, utr3len, trsp_genename, stopcodon, stop4nt]
			ucscIDlist.append(trsp_attr_list[0])
			transcriptdict[trsp_id] = trsp_attr_list
			total_transcripts += 1
			#transcript,chrom,featnum,strand,mrna_len,cds_len,5utr_len,3utr_len,gene_name,stopcodon,stop4nt 
	print "total number of transcripts in data table: %s" % total_transcripts
	print "Number of included chromosomes chr: %s" % validchroms
	print "Number of excluded chromosomes chr: %s" % nonvalidchorms
	print "included chroms: ", included_chroms
	print "excluded chroms: ", excluded_chroms
	print "transcripts discarded due to non-AUG start codon %s" % nonATGstart
	print "transcripts discarded due to noncanonical stop codon %s" % wrongstopcodon
	return ucscIDlist, transcriptdict
	


def write_utr_stopcodon_csvfile(ucscIDlist, transcriptdict, outfilestring, headers):
	# For writing dictionary to csv file first headers are written Then one line at a time is added to t, This is a list that will hold all of the trsp_attr_lists 
	# 	Finally these are written to each line of a csv file
	t=[]
	t.append(headers)

	for i in ucscIDlist: # position starts at 0
		newline= transcriptdict[i]
		t.append(newline)

	fa = open(outfilestring+".csv", "w")
	writer = csv.writer(fa)
	writer.writerows(t)
	fa.close()

def main():
	GTFgen = GFF.parse(GTFfile)
	GFFlist = makeGFFlist(GTFgen)
	outfilestring = '%s_UTRs' % (gtfInFilePrefix)
	if not os.path.isfile(outfilestring+".csv"):
		include_noncanon_start = True
		include_noncanon_stop = True
		headers= ['#transcript','chrom','featnum','strand','mrna_len','cds_len','5utr_len','3utr_len','gene_name','stopcodon','stop4nt']
		ucscIDlist, transcriptdict = build_utr_table(GFFlist, include_noncanon_start, include_noncanon_stop)
		write_utr_stopcodon_csvfile(ucscIDlist, transcriptdict, outfilestring, headers)

	

if __name__ == '__main__':
	# execute only if run as a script
	main()
"""
}


process Ribosome_Profiling_Build_Index_gtfToSingleApprisTranscript {

input:
 val gtf from g40_0_gtfPath_g40_1

output:
 val gtf_outfile  into g40_1_valid_gtf_g40_2

script:
gtf_outdir  = params.gtf.substring(0, params.gtf.lastIndexOf('/')) 
gtf_outfile = gtf_outdir+"/genes_protCode_TermStopCodon_validUTRs.gtf"
"""
#!/usr/bin/env python 
import sys, os
import pandas as pd 
import time
import argparse
from collections import OrderedDict

gtf_file_path = "${params.gtf}" 
gtf_rows_to_skip = 5

gtf_outdir = "${gtf_outdir}"
gtf_outfile = "${gtf_outfile}"

def read_in_gtf():
	# 1) read in the gtf file here
	df = pd.read_csv(gtf_file_path, sep='\t', dtype=str, header=None, skiprows=range(gtf_rows_to_skip))
	cols = ['#chrom', 'source', 'feature', 'chromStart', 'chromEnd', 'score', 'strand', 'frame', 'transcript_id']
	df.columns = cols
	return df


def parse_entry(tr):
	# given a single row of a gtf data.frame, parse column 8 into an ordered dictionary
	tr = tr.replace('"', '')
	trl = tr.split("; ")
	trdict = OrderedDict()

	for j in trl:
		k = j.split(" ")
		if k[0] in trdict:
			trdict[k[0]].append(k[1])
		else:       
			trdict[k[0]]=[k[1]]
	return trdict


def build_gene_indexes(df):
	# 2) take an input dataframe and output an ordered dict with the index slices for every gene
	geneDict = OrderedDict()

	geneCount = 0
	previousGeneIndex = 0

	current_id=""
	current_gene=""

	for i in range(len(df)):

		if df.loc[i,'feature'] == 'gene':
			trdict = parse_entry(df.loc[i,'transcript_id'])

			curGeneID = trdict['gene_id'][0]
		
			if geneCount != 0:
				newGeneIndex = i
				geneDict[current_id] = [previousGeneIndex,newGeneIndex]
				previousGeneIndex = i
				current_id = trdict['gene_id'][0]
				geneCount += 1

			else:
				newgeneIndex = 0
				geneCount +=1
				current_id = trdict['gene_id'][0]
		if i == (len(df)-1):
			newGeneIndex = i+1
			current_id = trdict['gene_id'][0]
			geneDict[current_id] = [previousGeneIndex,newGeneIndex]
	return geneDict

### -------- ###
			
def find_coding_genes(geneDict, df):
	## 3) filter geneDict on coding genes only

	geneDictCoding = OrderedDict() ### store coding transcripts only

	for gene in geneDict:
		tr = df.loc[geneDict[gene][0], 'transcript_id']
		
		trdict = parse_entry(tr)
		
		if trdict['gene_type'][0] == 'protein_coding':
			geneDictCoding[gene] = geneDict[gene]

	return geneDictCoding

def find_processed_pseudogenes(geneDict, df):
	# locate all of the processes pseudogenes

	pseudoGeneDict = OrderedDict()

	for gene in geneDict:

		tr = df.loc[geneDict[gene][0], 'transcript_id']
		
		trdict = parse_entry(tr)
		
		if trdict['gene_type'][0] == 'processed_pseudogene':
			pseudoGeneDict[gene] = geneDict[gene]
	return pseudoGeneDict



def find_overlapping_loci(geneDictCoding, df):
	# a huge fraction of genes have 'overlaps' this is not super useful

	overlapLociDict = OrderedDict()

	for gene in geneDictCoding:
		tr = df.loc[geneDictCoding[gene][0], 'transcript_id']

		trdict = parse_entry(tr)

		if 'tag' in trdict:
			if 'overlapping_locus' in trdict['tag']:
				overlapLociDict[gene] = trdict
				# overlapLociDict[trdict['transcript_id'][0]]= trdict['tag']

	print overlapLociDict
	print len(geneDictCoding)
	print len(overlapLociDict)



def build_transcript_indexes(geneDictCoding, df):
	# take an input geneDict and find the indexes for all transcripts associated with each gene

	TR_index_dict = OrderedDict()

	for gene in geneDictCoding:

		trDF = df.iloc[geneDictCoding[gene][0]:geneDictCoding[gene][1]]
	
		trPrev = -1
		trNamePrev = ""
		
		### iterate through a slice of the data frame for each gene
		### search for transcripts ofver that slice
		### find transcript slices
		for i in range(geneDictCoding[gene][0], geneDictCoding[gene][1]):
			if trDF.loc[i,'feature'] == 'transcript':
				trdict = parse_entry(trDF.loc[i,'transcript_id'])
				trCur = i
				trNameCur = trdict['transcript_id'][0]
				
				if trPrev != -1: # do not make an entry for the first transcript
					TR_index_dict[trNamePrev] = [trPrev, trCur]

				trPrev = trCur
				trNamePrev = trNameCur
			
			### for the final transcript
			if i == geneDictCoding[gene][1]-1:
				trdict = parse_entry(trDF.loc[i,'transcript_id'])
				TR_index_dict[trdict['transcript_id'][0]] = [trCur, i+1]
	return TR_index_dict



def choose_appris_canonical(geneDictCoding, df, TR_index_dict):

	#take gene and transcript dictionaries, and choose the best transcript accoding these criteria:
	#re-working this to include all protein coding, and all ccds genes

	geneDictCanon = OrderedDict()
	geneDictChrLoc = OrderedDict() # [trxStart, trxEnd, tr_for_start, tr_for_stop]

	ap1 = 0
	ap2 = 0
	ap3 = 0
	ap4 = 0
	ap5 = 0

	ap1_alt = 0
	ap2_alt = 0

	single_iso_pri = 0
	single_iso_alt = 0
	multi_iso = 0
	total_iso = 0
	total_iso_pri = 0
	total_iso_alt = 0

	identical_trsps = 0
	noAppris = 0
	
	noCCDS = 0
	singleCCDS = 0
	multiCCDS = 0

	noProtCode = 0
	singleProtCode = 0
	multiProtCode = 0

	protCodeCount = 0
	cds_NF_count = 0
	noValidUTRsRemaining = 0
	invalidUTRtotal = 0

	### now we have all genes that are coding, with start and stop indexes:
	for gene in geneDictCoding:


		trDF = df.iloc[geneDictCoding[gene][0]:geneDictCoding[gene][1]]
		trDFz = trDF.reset_index(drop=True) # z is for zero based conversion here   

		### Now we are finding the CCDS transcript only
		### and retrieving the 'best' transcript according to appris_principal
		
		protCodeTrDict = OrderedDict()
		ccdsTrDict = OrderedDict() ### store all transcripts that have CCDS annotation
		trCount = 0
		# protCodeCount = 0
		ccdsCount = 0
		
		### find all protein coding transcripts
		for i in range(len(trDFz)):
			if trDFz.loc[i, 'feature'] == 'transcript':
				tr = trDFz.loc[i, 'transcript_id']
				trdict = parse_entry(tr)
				if trdict['transcript_type'][0] == 'protein_coding':

					### discard transcripts without completed coding regions:
					# if 'cds_start_NF' in trdict['tag'] or 'cds_end_NF' in trdict['tag']:
					if 'cds_start_NF' in trdict['tag'] or 'cds_end_NF' in trdict['tag'] or 'mRNA_start_NF' in trdict['tag'] or 'mRNA_end_NF' in trdict['tag']:
						cds_NF_count +=1
						# print "cds_NF for %s" % (trdict['transcript_id'])
						# print trdict['tag']
					else:
						protCodeTrDict[trdict['transcript_id'][0]] = trdict['tag'] 
						protCodeCount += 1

		# print protCodeCount
		# print protCodeTrDict.keys()


		### find all CCDS transcripts
		for i in range(len(trDFz)):
			if trDFz.loc[i, 'feature'] == 'transcript':
				trCount +=1
				### check for ccds
				tr = trDFz.loc[i, 'transcript_id']
				trdict = parse_entry(tr)
				if 'tag' in trdict:
					if 'CCDS' in trdict['tag']:
						ccdsTrDict[trdict['transcript_id'][0]]= trdict['tag']
						ccdsCount += 1

		# print ccdsCount
		# print ccdsTrDict.keys()
		# print ccdsTrDict

		### find all appris transcripts for the gene
		principal_iso_count = 0 
		alt_iso_count = 0
		
		ap1dict = OrderedDict()
		ap2dict = OrderedDict()
		ap3dict = OrderedDict()
		ap4dict = OrderedDict()
		ap5dict = OrderedDict()

		ap1dict_alt = OrderedDict()
		ap2dict_alt = OrderedDict()
		
		apDictAll_pri = OrderedDict()
		apDictAll_alt = OrderedDict()
		
		for tran in ccdsTrDict:

			if 'appris_principal_1' in ccdsTrDict[tran]:
				ap1dict[tran] = TR_index_dict[tran]
				apDictAll_pri[tran] = 1
				ap1 += 1
				principal_iso_count += 1

			elif 'appris_principal_2' in ccdsTrDict[tran]:
				ap2dict[tran] = TR_index_dict[tran]
				apDictAll_pri[tran] = 2
				ap2 += 1
				principal_iso_count += 1

			elif 'appris_principal_3' in ccdsTrDict[tran]:
				ap3dict[tran] = TR_index_dict[tran]
				apDictAll_pri[tran] = 3
				ap3 += 1
				principal_iso_count += 1

			elif 'appris_principal_4' in ccdsTrDict[tran]:
				ap4dict[tran] = TR_index_dict[tran]
				apDictAll_pri[tran] = 4
				ap4 += 1
				principal_iso_count += 1

			elif 'appris_principal_5' in ccdsTrDict[tran]:
				ap5dict[tran] = TR_index_dict[tran]
				apDictAll_pri[tran] = 5
				ap5 += 1
				principal_iso_count += 1

			elif 'appris_alternative_1' in ccdsTrDict[tran]:
				ap1dict_alt[tran] = TR_index_dict[tran]
				apDictAll_alt[tran] = 6
				ap1_alt += 1
				alt_iso_count += 1

			elif 'appris_alternative_2' in ccdsTrDict[tran]:
				ap2dict_alt[tran] = TR_index_dict[tran]
				apDictAll_alt[tran] = 7
				ap2_alt += 1
				alt_iso_count += 1

		apDictAll = OrderedDict(apDictAll_pri.items()+apDictAll_alt.items())

		# ### add a dictionary to output all appri transcripts with start and stop indexes from df
		# for trnscrpt in apDictAll:
		# 	geneDictAllAppris[trnscrpt] = [gene, TR_index_dict[trnscrpt]] 

		total_iso_pri += principal_iso_count
		total_iso_alt += alt_iso_count ### ALT work in progress
		total_iso += principal_iso_count+alt_iso_count
		

		###############
		## Prot Code ##
		###############


		if len(protCodeTrDict) == 0: 
			# print "no protein coding genes here!"
			noProtCode +=1 ## these should have already been filtered out
			continue



		###
		### Build dfcomp first
		### first build dfcomp for comparison features of each transcript for comparison
		dfcomp = pd.DataFrame(columns=['geneID', 'trsp', 'strand', 'trxStart', 'trxEnd',
									  'scStart', 'scEnd', 'startCodonStart', 'startCodonEnd',
									   'cdsLen', 'utr3Len', 'utr5Len', 'apri_val', 'apri_class', 'CCDS'])
		
		### define whether transcript is a primary or alternative 
		pri_vals = [1,2,3,4,5]
		alt_vals = [6,7]


		for tr in protCodeTrDict.keys():
			# print tr

			### check for appris entry
			if tr in apDictAll.keys():
				# print apDictAll[tr]
				apri_val = apDictAll[tr]

				if apri_val in pri_vals:
					apri_class = 'pri'
				elif apri_val in alt_vals:
					apri_class = 'alt'
				else:
					apri_class = 'unknown'
			else:
				apri_val = 0
				apri_class = "none"

			### check for CCDS entry
			if tr in ccdsTrDict.keys():
				ccds_present = 'yes'
			else:
				ccds_present = 'no'


			### make a df slice of each transcript to retrieve data
			tempdf = df.iloc[TR_index_dict[tr][0]:TR_index_dict[tr][1]]
			tempdf.reset_index(drop=True, inplace=True) ## set to 0-based index
			
			strand = tempdf.loc[0]['strand']
			cdsLen = 3 ## account for inclusion of stop codon, this is not included in gtf file
			utr3Len = 0
			utr5Len = 0
			
			utrVals = []
		
			for row in range(len(tempdf)):
				## define transcript start and end positions
				if tempdf.loc[row]['feature'] == 'transcript':
					trxStart = tempdf.loc[row]['chromStart']
					trxEnd = tempdf.loc[row]['chromEnd']

				## define stop_codon start and end positions
				if tempdf.loc[row]['feature'] == 'stop_codon':
					scStart = tempdf.loc[row]['chromStart']
					scEnd = tempdf.loc[row]['chromEnd']

				## define start_codon start and end positions
				if tempdf.loc[row]['feature'] == 'start_codon':
					startCodonStart = tempdf.loc[row]['chromStart']
					startCodonEnd = tempdf.loc[row]['chromEnd']
			
				## define length of coding sequence
				if tempdf.loc[row]['feature'] == 'CDS':
					cdsLen += int(tempdf.loc[row]['chromEnd']) - int(tempdf.loc[row]['chromStart']) + 1
				
				## build list of all UTR values
				if tempdf.loc[row]['feature'] == 'UTR':
					### 3'UTR
					if strand == '+' and tempdf.loc[row]['chromStart'] >= scStart:
						utr3Len += int(tempdf.loc[row]['chromEnd']) - int(tempdf.loc[row]['chromStart']) + 1
					if strand == '-' and tempdf.loc[row]['chromEnd'] <= scEnd:
						utr3Len += int(tempdf.loc[row]['chromEnd']) - int(tempdf.loc[row]['chromStart']) + 1
					### 5'UTR
					if strand == '+' and tempdf.loc[row]['chromEnd'] < startCodonStart:
						utr5Len += int(tempdf.loc[row]['chromEnd']) - int(tempdf.loc[row]['chromStart']) + 1
					if strand == '-' and tempdf.loc[row]['chromStart'] > startCodonEnd:
						utr5Len += int(tempdf.loc[row]['chromEnd']) - int(tempdf.loc[row]['chromStart']) + 1
			### assemble all features into a list, in order of columns
			df_tr_entry = [[gene, tr, strand, trxStart, trxEnd, scStart, scEnd, 
							startCodonStart, startCodonEnd, cdsLen, utr3Len, utr5Len, 
							apri_val, apri_class, ccds_present]]

			 ## build a temparary df to hold each transcript
			df_to_add = pd.DataFrame(df_tr_entry, columns = dfcomp.columns)
			
			## add transcript df to gene df
			dfcomp = dfcomp.append(df_to_add, ignore_index=True)

			# print dfcomp

		### output the minimum and maximum chromosome start and end positions for all coding transcripts
		# print dfcomp
		trxStartmin = dfcomp['trxStart'].min()
		trxEndmax = dfcomp['trxEnd'].max()
		tr_min = dfcomp.loc[pd.to_numeric(dfcomp['trxStart']).idxmin()]['trsp']
		tr_max = dfcomp.loc[pd.to_numeric(dfcomp['trxEnd']).idxmax()]['trsp']

		geneDictChrLoc[gene] = [trxStartmin, trxEndmax, tr_min, tr_max]

		### Check for valid UTR lengths
			## throw out transcripts with UTR lengths of zero
		pc_tr_count = len(dfcomp)

		dfcomp = dfcomp.loc[dfcomp['utr3Len']>3] # 3 is minmum UTR length here
		dfcomp = dfcomp.loc[dfcomp['utr5Len']>0]

		valUTR_tr_count = len(dfcomp)

		invalid_UTR_count = pc_tr_count - valUTR_tr_count
		invalidUTRtotal += invalid_UTR_count
		# print invalid_UTR_count, "INVAL UTR"

		### discard this gene if no valid transcripts remain
		if len(dfcomp) == 0:
			noValidUTRsRemaining +=1
			continue

		### selecting isoform: 
		if len(dfcomp) == 1:
			singleProtCode +=1
			# tr = protCodeTrDict.keys()[0] ### I think this is wrong here
			dfcomp = dfcomp.reset_index(drop=True)
			tr = dfcomp.loc[0, 'trsp']
			# print tr, "single ccds"

			### !!!! output transcript here
			geneDictCanon[tr] = [gene, TR_index_dict[tr]] 

			### build min and max genome locations for this gene
			tempdf = df.iloc[TR_index_dict[tr][0]:TR_index_dict[tr][1]]
			tempdf.reset_index(drop=True, inplace=True) ## set to 0-based index

			tempdf = tempdf.loc[tempdf['feature'] == 'transcript'] ## should only be one entry

			if len(tempdf) != 1:
				print "not single transcript for %s" % gene
				sys.exit()

			trxStart = tempdf.loc[0]['chromStart']
			trxEnd = tempdf.loc[0]['chromEnd']

			# for row in range(len(tempdf)): ## safest way to make sure that it is the whole transcript entry
			# 		## define transcript start and end positions
			# 		if tempdf.loc[row]['feature'] == 'transcript':
			# 			trxStart = tempdf.loc[row]['chromStart']
			# 			trxEnd = tempdf.loc[row]['chromEnd']

			geneDictChrLoc[gene] = [trxStart, trxEnd, tr, tr]

		elif len(protCodeTrDict) > 1:
			
			multiProtCode += 1

			### QC1 - verify that transcripts map to one strand, these have overlaping transcripts:
			strnds = dfcomp['strand'].unique()
			if len(strnds) > 1:
				print "MULTI STRAND for %s" % gene ## There are none in genecode annotation file
				sys.exit()


			# if gene == 'ENSG00000189195.13':
			# 	print dfcomp

			# ### QC2 - Check for overlaps between transcripts in dfcomp ###
			# overlaps = []
			# for j in range(len(dfcomp)):
			# 	trxspan = set(range(int(dfcomp.loc[j]['trxStart']), int(dfcomp.loc[j]['trxEnd'])))
			# 	overlaps.append(trxspan)

			# interIndex = []
			# maxcount = len(overlaps)
			# ## iteratively check each transcript for overlap with all other transcripts
			# for x in range(len(overlaps)):

			# 	y = x +1
			# 	while y < maxcount:
			# 		intersect = overlaps[x].intersection(overlaps[y])
			# 		if len(intersect) > 0:
			# 			interIndex.append([x,y])
			# 		y = y+ 1
			# if len(interIndex) == 0:
			# 	print "no overlaps for %s" % gene ### all genes in annotation have overlaps, must be a requirement for appris principals...
			# 	### take only transcripts from principal isoform
			# 	dfcomp = dfcomp.loc[dfcomp['apri_class'] == 'pri'] ## chane this to only the principal values and continue with procesing
			# 	print dfcomp
			#### -- end overlaps -- ###



			### QC3 - Take 3' most stop codon
			unique_stops = dfcomp['scEnd'].unique()

			### check for multiple stop codons, there is one case of this for ENSG00000108395.14 (TRIM37), dif exons that give same cds
			if len(unique_stops) > 1:
				# print "MULTIPLE STOP CODONS for %s" % gene 
				# print dfcomp
				# print unique_stops
				# print strnds 

				### take subset of transcripts with terminal stop codon
				if strnds == "+":
					dfcomp = dfcomp.loc[dfcomp['scEnd'] == unique_stops.max()]
				if strnds == "-":
					dfcomp = dfcomp.loc[dfcomp['scEnd'] == unique_stops.min()]
				# print dfcomp

			# if dfcomp.loc[0]['geneID'] == 'ENSG00000108395.14':
			# 	print dfcomp
			
			remaining_classes = dfcomp['apri_class'].unique()
			# print remaining_classes

			### QC4 - filter by appris classification

			if 'pri' in remaining_classes:
				dfcomp = dfcomp.loc[dfcomp['apri_class'] == 'pri']
			elif 'alt' in remaining_classes:
				dfcomp = dfcomp.loc[dfcomp['apri_class'] ==  'alt']
			elif 'none' in remaining_classes:
				dfcomp = dfcomp.loc[dfcomp['apri_class'] == 'none']
			else:
				print 'no valid classes for %s' % gene
				sys.exit()

			### filter by CCDS classification

			ccds_remaining = dfcomp['CCDS'].unique()

			if 'yes' in ccds_remaining:
				dfcomp =dfcomp.loc[dfcomp['CCDS'] == 'yes']
			elif 'no' in ccds_remaining:
				dfcomp = dfcomp.loc[dfcomp['CCDS'] == 'no']
			else:
				print 'ambiguous CCDS classification for %s' % gene
				sys.exit()


			### QC5 - take longest coding sequences

			if len(dfcomp) > 1:
				cdsLenUnique = dfcomp['cdsLen'].unique()
				# print cdsLenUniqueAlt.max()
				dfcomp = dfcomp.loc[dfcomp['cdsLen'] == cdsLenUnique.max()]

			# print dfcomp

			### >>
			###  Select transcripts with MIN 3'UTR length ###
			### >>

			### take minimum utr3 if longer than 3 nt's, unless only 1 transcript and it is that length

			min_utr3_thresh = 3

			utr3lengths = dfcomp['utr3Len'].unique() ## list of 

			if utr3lengths.min() > min_utr3_thresh: # if longer than threshold, simply take the minimum
				utr3min = utr3lengths.min()
			elif utr3lengths.min() <= min_utr3_thresh and len(utr3lengths) == 1: # if only one utr3 value and below the min, just use this
				utr3min = utr3lengths.min()
			else:
				print utr3lengths
				utr3len_thresh = [x for x in utr3lengths if x > min_utr3_thresh]
				print utr3len_thresh
				if len(utr3len_thresh) > 0:
					utr3min = min(utr3len_thresh)
				else:
					utr3min = utr3lengths.min()
				print utr3min

			# utr3min = dfcomp['utr3Len'].min() ### set minimum value for 3'UTR
			# print utr3min
			utr3valCounts = dfcomp['utr3Len'].value_counts() ### get number of occarnaces of each 3'utr length

			utr3minOccurances = utr3valCounts[utr3min]
			
			### single transcript with this 3'UTR
			if utr3minOccurances == 1:
				minInd = pd.to_numeric(dfcomp['utr3Len']).idxmin()
				outTr = dfcomp.loc[minInd,'trsp']
				
				### !!!!! output this transcript here
				geneDictCanon[outTr] = [gene, TR_index_dict[outTr]]
			
			if utr3minOccurances > 1:
				# print "multiple 3'UTRs"


				dfcomp = dfcomp.loc[dfcomp['utr3Len'] == utr3min] ### reset dfcomp to only include min utr3's
				# print dfcomp
				### next take transcript with shortest 5'UTR
				utr5min = dfcomp['utr5Len'].min()
				# if utr5min == 3:
				# 	print "no utr for %s " % gene 
				utr5valCounts = dfcomp['utr5Len'].value_counts()
				utr5minOccurances = utr5valCounts[utr5min]
				
				if utr5minOccurances == 1:
					minUtr5Ind = pd.to_numeric(dfcomp['utr5Len']).idxmin()
					outTr = dfcomp.loc[minUtr5Ind, 'trsp'] ### selected transcript
					
					### !!!!! output a transcript here
					geneDictCanon[outTr] = [gene, TR_index_dict[outTr]]
					
				if utr5minOccurances > 1:
					### some transcripts are completely identical in genecode annotation, must be historical reason
					### only difference is transcript name, and some of the annotations in last column
					# print "MULTIPLE VALID trsps for %s" % gene
					dfcomp = dfcomp.reset_index(drop=True)
					print dfcomp
					identical_trsps += 1
					outTr = dfcomp.loc[0, 'trsp'] # take the first one, since these are identical in terms of genomic features
					### !!!!! output transcript here
					geneDictCanon[outTr] = [gene, TR_index_dict[outTr]] 


	print "" # print out some summaries 
	print "Number of identical transcripts: ", identical_trsps #20 of these in the genome that fall out at the end
	print "Total number of appris isoforms: ", total_iso
	print "Number of genes with multiple appris isoforms: ", multi_iso
	# print "Number of genes with single appris isofrom: ", single_iso
	# print geneDictCanon   
	print "Count for appris principals 1/2/3/4/5/6/7 isoforms: ", ap1, ap2, ap3, ap4, ap5, ap1_alt, ap2_alt
	# print geneDictCanon
	print "CCDS genes without appris isofomrs: ", noAppris
	print "genes included == %s" % (len(geneDictCanon))
	print "genes in protCodeTrDict == %s" % (protCodeCount)
	print "cds_NF count %s" % (cds_NF_count)
	print "protein coding transcripts == %s" % (protCodeCount)
	print "genes discarded due to zero length UTRs == %s" % (noValidUTRsRemaining)
	print "total transcripts with invalid UTRs == %s" % (invalidUTRtotal)
	# print "length of geneDictAllAppris == %s" % (len(geneDictAllAppris))

	return geneDictCanon, geneDictChrLoc


def build_df_dict(df, geneDictCanon):
	## for each transcript selected in geneDictCanon, take a dataframe slice from df that has all entries for the transcript
	## currently keying the output dictionary on 'geneIDs', could be switched to transcript id's later if desired

	outDict = OrderedDict()

	for tr in geneDictCanon:
		outDict[geneDictCanon[tr][0]] = df.iloc[geneDictCanon[tr][1][0]:geneDictCanon[tr][1][1]]

	return outDict


def overlap_features(genedf):
	# for a given data.frame of a transcript, extract info to check for overlaps
	# Used in find_transcript_overlaps()
	# output: [chrStart, chrEnd, strand, geneName, chrom]
	
	chrStart = genedf.loc[0,'chromStart']
	chrEnd = genedf.loc[0, 'chromEnd']
	strand = genedf.loc[0, 'strand']
	trdict = parse_entry(genedf.loc[0, 'transcript_id'])
	geneName = trdict['gene_id'][0]
	trName = trdict['transcript_id'][0]
	chrom = genedf.loc[0, '#chrom']
	outlist = [chrStart, chrEnd, strand, geneName, chrom, trName]
	return outlist

def find_transcript_overlaps(outDict, geneDictChrLoc):
	# For each transcript, look at nearest neighbor on same strand and check if transcripts have overlaps
	# add all transcripts with overlaps to exclusion_overlaps dictionary
	
	exclusion_overlaps = OrderedDict()
	excluded_trsp_count = 0
	gene_key_list = outDict.keys() ### ordered list of gene names in the dictionary

	trsp_entry = -1
	for gene in outDict:

		trsp_entry += 1
		genedf = outDict[gene]
		genedf = genedf.reset_index(drop=True)
		
		# chrStart = genedf.loc[0,'chromStart']
		# chrEnd = genedf.loc[0, 'chromEnd']
		# strand = genedf.loc[0, 'strand']
		# trdict = parse_entry(genedf.loc[0, 'transcript_id'])
		# geneName = trdict['gene_id'][0]
		
		overlap_feats = overlap_features(genedf)
		
		cur_strand = overlap_feats[2]
		cur_chrom = overlap_feats[4]
		
		##### DOWNSTREAM Overlaps 
		next_strand = "0"
		search_index = 1
		break_signal = 0
		
		while next_strand != cur_strand:
		
			next_tr = trsp_entry+search_index ### only for the end of the list
			if next_tr == len(gene_key_list):
				print "reached end of downstream list for %s " % gene
				break_signal = 1
				break
			
			next_gene = gene_key_list[next_tr] ### adding this
			nextdf = outDict[gene_key_list[next_tr]]
			nextdf = nextdf.reset_index(drop=True)
			next_over = overlap_features(nextdf)

			### adding genome position
			
			next_chrom = next_over[4]
			
			### check to make sure that this is still the same chr
			if next_chrom != cur_chrom:
				print "search downstream to next chr for %s" % gene
				break_signal = 1
				break

			next_strand = next_over[2]
			search_index +=1
			
		next_gene_start = geneDictChrLoc[next_gene][0]

		# if int(overlap_feats[1]) > int(next_over[0]) and break_signal !=1:
		if int(overlap_feats[1]) > int(next_gene_start) and break_signal !=1:
			# print "******"
			# print "OVERLAPING DOWNSTREAM TRSP %s, %s" % (gene, overlap_feats[5])
			# print "******"
			# print overlap_feats
			# # print next_over
			# print next_gene

			exclusion_overlaps[gene] = overlap_feats
			excluded_trsp_count += 1
			
			
		##### UPSTREAM Overlaps
		next_strand = "0"
		search_index = 1
		break_signal = 0
		
		while next_strand != cur_strand:
			next_tr = trsp_entry-search_index ### this time subtract search index
			if next_tr == -1:
				print "end of upstream list for %s" % gene
				break_signal = 1
				break
				
			next_gene = gene_key_list[next_tr] ### adding this
			nextdf = outDict[gene_key_list[next_tr]]
			nextdf = nextdf.reset_index(drop=True)
			next_over = overlap_features(nextdf)
			
			next_chrom = next_over[4]

			### check to make sure that this is still the same chr
			if next_chrom != cur_chrom:
				print "search upstream to prev chr for %s" % gene
				break_signal = 1
				break

			next_strand = next_over[2]
			search_index +=1
			
		next_gene_end = geneDictChrLoc[next_gene][1]

		# if int(overlap_feats[0]) < int(next_over[1]) and break_signal !=1:
		if int(overlap_feats[0]) < int(next_gene_end) and break_signal !=1:
			# print "******"
			# print "OVERLAPING UPSTREAM TRSP %s, %s" % (gene, overlap_feats[5])
			# print "******"
			# print overlap_feats
			# # print next_over
			# print next_gene
			
			exclusion_overlaps[gene] = overlap_feats
			excluded_trsp_count += 1

	print "Total Overlapping transcripts revised: %s" % excluded_trsp_count
	return exclusion_overlaps


def remove_overlapping_transcripts(outDict, exclusion_overlaps):
	# simply remove all of the transcripts in exclusion_overlaps from outDict

	outDictExclu = outDict

	for exclu in exclusion_overlaps:
		outDictExclu.pop(exclu)
		
	print len(outDictExclu)
	return outDictExclu


def edit_col8(dfIn):
	# Edit the format of column 8 to remove spaces and add:
	# 	ID= Name= for transcript
	#	Parent= Name= for other entires
	#Sort the entries based on the value of chromStart

	#	This will be compaitble with downstream GTF parser
	
	for i in dfIn.index:
		tr = parse_entry(dfIn.loc[i,'transcript_id'])

		### if line is a transcript, add Id= Name=
		if dfIn.loc[i,"feature"] == 'transcript':

			### make entry of new identifies to work with GFF parser 
			newline = "ID=%s;Name=%s" % ("".join(tr['transcript_id']), "".join(tr['gene_name']))
			### remove spaces and add semicolon to items in col 8
			line8 = ";".join("=".join((str(k),str(",".join(v)))) for k,v in tr.items())
			### Merge these into a single string
			outline = newline+";"+line8 
			outline = outline[0:-1] ## strip off last ';'

			### set the value of the cell to the edited line
			dfIn.at[i,"transcript_id"] = outline


		### for all other lines, add Parent= Name= 
		else:
			newline = "Parent=%s;Name=%s" % ("".join(tr['transcript_id']), "".join(tr['gene_name']))
			line8 = ";".join("=".join((str(k),str(",".join(v)))) for k,v in tr.items())
			outline = newline+";"+line8
			outline = outline[0:-1]
			dfIn.at[i,"transcript_id"] = outline
			
	dfIn['sort_vals'] = dfIn['chromStart'].astype(int)
	dfIn = dfIn.sort_values(by=['sort_vals'], axis=0, ascending=True)
	dfIn = dfIn.drop(['sort_vals'], axis=1)
	dfOut = dfIn.copy()
	
	return dfOut


def define_pseudogene_positions(pseudoGeneDict, df):

	pseudoChromDict = OrderedDict()

	for gene in pseudoGeneDict:

		# print gene
		tempdf = df.iloc[pseudoGeneDict[gene][0]:pseudoGeneDict[gene][1]]
		tempdf = tempdf.loc[tempdf['feature']=='transcript']
		tempdf = tempdf.reset_index(drop=True)

		# print tempdf
		overlap_feats = overlap_features(tempdf)
		# print overlap_feats

		over_range = overlap_feats[0:2]
		over_range = [int(x) for x in over_range]
		# over_range = map(int, over_range)
		# cur_chrom = overlap_feats[4]
		cur_chrom = overlap_feats[4]+"_"+overlap_feats[2] ## account for strandedness
		# print cur_chrom

		if cur_chrom in pseudoChromDict:
			pseudoChromDict[cur_chrom].append(over_range)
		else:
			pseudoChromDict[cur_chrom] = [over_range]

			# pseudoChromDict[cur_chrom] = pseudoChromDict[cur_chrom].append(over_range)

	# print pseudoChromDict
	return pseudoChromDict

	# print pseudoChromDict


def check_pseudo_overlap(x1, x2, y1, y2):
	# check for overlaps [x1, x2] and [y1, y2]
	# should only be true if these overlap
	return max(x1, y1) <= min(x2, y2) 


def find_pseudo_overlaps(outDictExclu, pseudoChromDict):

	pseudo_exclude = OrderedDict()
	pseudo_exclude_count = 0

	for gene in outDictExclu:
		# print gene
		# print outDictExclu[gene]

		genedf = outDictExclu[gene]
		genedf = genedf.reset_index(drop=True)

		for feat in genedf.index: ### only check UTR entries, not perfect, but can be modified, maybe to exons?
			if genedf.loc[feat, 'feature'] == 'UTR':
				chrStart = genedf.loc[feat,'chromStart']
				chrEnd = genedf.loc[feat, 'chromEnd']
				strand = genedf.loc[feat, 'strand']
				trdict = parse_entry(genedf.loc[feat, 'transcript_id'])
				geneName = trdict['gene_id'][0]
				trName = trdict['transcript_id'][0]
				chrom = genedf.loc[feat, '#chrom']
				outlist = [chrStart, chrEnd, strand, geneName, chrom, trName]

				overlap_feats = outlist
				# overlap_feats = overlap_features(genedf)
				
				gene_over_range = overlap_feats[0:2]
				gene_over_range = [int(x) for x in gene_over_range]
				# gene_over_range = map(int, gene_over_range)
				cur_chrom = overlap_feats[4]+"_"+overlap_feats[2]

				# print gene_over_range
				# print gene_over_range[0]
				# print gene_over_range[1]

				try: # 'chrM_+' is not in pseudoChromDict ... 
					chr_ref_list = pseudoChromDict[cur_chrom]
				except KeyError as e:
					# pseudoChromDict[cur_chrom]
					print(e.message)
					print "no key in pseudoChromDict"
					continue

				for pseudo in chr_ref_list:
					# print pseudo
					overTest = check_pseudo_overlap(gene_over_range[0], gene_over_range[1], pseudo[0], pseudo[1])

					if overTest == True:
						print "overlap found here for %s with pseudo %s" % (gene, pseudo)
						pseudo_exclude[gene] = [cur_chrom, pseudo]
						pseudo_exclude_count += 1

	print pseudo_exclude_count, "excluded pseudogene count"

	return pseudo_exclude

def mod_last_column(outDictExclu):
	# Edit the final column of each dataframe
	
	outDictMod = OrderedDict()

	for key in outDictExclu:
		k = key
		dfIn = outDictExclu[key].copy()
		dfMod = edit_col8(dfIn)
		outDictMod[k] = dfMod
	
	return outDictMod

def output_df(outdict, out_file):
	# compile all modified dataframes for each transcript into a master dataframe
	# build the output dataframe from the modified dictionary and write this to a file:
	
	cols = ['#chrom', 'source', 'feature', 'chromStart', 'chromEnd', 'score', 'strand', 'frame', 'transcript_id']
	colOut = ['#chrom', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'transcript_id']
	gtfDF = pd.DataFrame(columns=cols)

	for trsp in outdict:
		gtfDF = gtfDF.append(outdict[trsp], ignore_index=True)
		
	gtfDF.columns = colOut
	print gtfDF.head(5)
	gtfDF.to_csv(out_file, sep='\t', index=False)


def main():
	import os.path
	if not os.path.isfile(gtf_outfile):
		time_start = time.time()
		df = read_in_gtf()
		print df.head()
		geneDict = build_gene_indexes(df)
		geneDictCoding = find_coding_genes(geneDict, df)
		pseudoGeneDict = find_processed_pseudogenes(geneDict, df)
		pseudoChromDict = define_pseudogene_positions(pseudoGeneDict, df)

		TR_index_dict = build_transcript_indexes(geneDictCoding, df)
		geneDictCanon, geneDictChrLoc = choose_appris_canonical(geneDictCoding, df, TR_index_dict)

		outDict = build_df_dict(df, geneDictCanon)
		exclusion_overlaps = find_transcript_overlaps(outDict, geneDictChrLoc)
		outDictExclu = remove_overlapping_transcripts(outDict, exclusion_overlaps)
		pseudo_exclude = find_pseudo_overlaps(outDictExclu, pseudoChromDict)

		outDictExcluPseudo = remove_overlapping_transcripts(outDict=outDictExclu, exclusion_overlaps=pseudo_exclude)

		outDictMod = mod_last_column(outDictExcluPseudo)
		output_df(outDictMod, gtf_outfile)
		time_end = time.time()
		print "total time == ", (time_end - time_start)
	else:
		print "protCode_TermStopCodon_validUTRs.gtf found"

if __name__ == '__main__':
	main()
"""
}

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 10
    $MEMORY = 30
}
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 3000
    $CPU  = 10
    $MEMORY = 25
    $QUEUE = "long"
} 
//* platform
//* autofill

process Ribosome_Profiling_Build_Index_build_annotation_files {

input:
 val valid_gtf_path from g40_1_valid_gtf_g40_2

output:
 val valid_gtf_path  into g40_2_valid_gtf_g40_7

script:
gtfInFilePrefix = valid_gtf_path.substring(0, valid_gtf_path.lastIndexOf('.')) 
gtfDir = valid_gtf_path.substring(0, valid_gtf_path.lastIndexOf('/')) 
"""
#!/usr/bin/env python 
import sys, os
import GFF
import twobitreader
import argparse
from Bio import SeqIO
import csv
from Bio.Seq import Seq
import pandas as pd
from collections import OrderedDict
from pathos.multiprocessing import ProcessingPool as Pool


# Add list of acceptable chromosomes that will be output to the table
validChrs = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 
			'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15',
			'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 
			'chrX', 'chrY', 'chrM', 'chrSinV', 'chrLUC', 'egfp', 'mcherry']

codonList = [
	'AAA', 'AAC', 'AAG', 'AAT',
	'ACA', 'ACC', 'ACG', 'ACT',
	'AGA', 'AGC', 'AGG', 'AGT',
	'ATA', 'ATC', 'ATG', 'ATT',
	'CAA', 'CAC', 'CAG', 'CAT',
	'CCA', 'CCC', 'CCG', 'CCT',
	'CGA', 'CGC', 'CGG', 'CGT',
	'CTA', 'CTC', 'CTG', 'CTT',
	'GAA', 'GAC', 'GAG', 'GAT',
	'GCA', 'GCC', 'GCG', 'GCT',
	'GGA', 'GGC', 'GGG', 'GGT',
	'GTA', 'GTC', 'GTG', 'GTT',
	'TAA', 'TAC', 'TAG', 'TAT',
	'TCA', 'TCC', 'TCG', 'TCT',
	'TGA', 'TGC', 'TGG', 'TGT',
	'TTA', 'TTC', 'TTG', 'TTT'  
]


### specify annotation files here
GenAnPath = "${gtfDir}"
GTFfile = "${valid_gtf_path}"
gtfInFilePrefix = "${gtfInFilePrefix}"
twobitfile = "${params.genome2bit}"
genome = twobitreader.TwoBitFile(twobitfile) # create handler to open the 2bit file
csvOutDir = "%s_codons" % gtfInFilePrefix
genAnName = gtfInFilePrefix.split('/')[-1]


### functions imported from Colin's densbuilder
def makeGFFlist(GTFinput):
	# Create a dictionary with a key for each chromosome in the GFF file
	GTFlist={}
	for chr in GTFinput:
		GTFlist[chr.id]=chr
	return GTFlist

def chrpostomrnapos(chrpos,chrom,featnum,GFFlist):
	# This funciton takes a genomic query position (chrpos) and a chromosome number (ex 'chr6')
	#	along with the feature number (the entry in the gff file for that transcript) defined by build_utr_table()
	#	and the dictionary of transcript from the GFFlist
	# The output is mrnapos which is the transcript relative position (position along the mRNA)
	#	of the original genomic query postion (chrpos)
	
	trsp_id= GFFlist[chrom].features[featnum].id
	trsp_strand= GFFlist[chrom].features[featnum].strand
	trsp_chromstart= int(GFFlist[chrom].features[featnum].location.start.position)  # 0-based
	trsp_chromend= int(GFFlist[chrom].features[featnum].location.end.position)
	sublist=[]

	for subfeature in GFFlist[chrom].features[featnum].sub_features:     # Make list of features
		if subfeature.type== 'exon':
			start= subfeature.location.start.position
			end= subfeature.location.end.position
			sublist.append([start,end])

	if trsp_strand== -1:    
		sublist.reverse()
	assert len(sublist)!= 0, ("transcript %s has a sublist length of zero!" % trsp_id)

	prevexonlen= 0 
	for item in sublist:
		exonstart= item[0]
		exonend= item[1]
		exonlen= exonend- exonstart

		if trsp_strand== 1:
			if chrpos>= exonstart and chrpos< exonend:      
				mrnapos= prevexonlen+ chrpos- exonstart
				return mrnapos
			else:   prevexonlen+= exonlen
		else:
			if chrpos< exonend and chrpos>= exonstart:
				mrnapos= prevexonlen+ (exonend-1)- chrpos       # Need -1 because end is not actual end, it is 1 beyond end.
				return mrnapos
			else:   prevexonlen+= exonlen 

def get_mRNA_sequence(GFFlist):

	transcriptdict={}
	ucscIDlist = []
	total_transcripts = 0 
	nonvalidchorms = 0
	nonATGstart = 0
	wrongstopcodon = 0
	validchroms = 0
	excluded_chroms = []
	included_chroms = []
	for chrom in GFFlist:
		if not chrom in validChrs:
			excluded_chroms.append(chrom)
			nonvalidchorms += 1
			# print chrom
			continue    # check that only valid choromosomes are used
		validchroms+=1
		included_chroms.append(chrom)
		transcriptnum= -1 # set to negative one so first transcript is == to 0
		for transcript in GFFlist[chrom].features:  # this is where the SeqFeatures are actually stored
			tr_attribute_list = []
			transcriptnum+=1
			trsp_id= transcript.id # it is a number 
			trsp_strand= transcript.strand
			trsp_genename= transcript.qualifiers['Name'][0]
			trsp_chromstart= int(transcript.location.start.position)  # 0-based
			trsp_chromend= int(transcript.location.end.position) 
			transcriptlist= [0.0 for x in range(abs(trsp_chromend- trsp_chromstart))] # a list for transcript (pre-mRNA), not CDS
			
			exonsplicedseq= SeqIO.Seq('')
			transcriptseq= SeqIO.Seq(genome[chrom][trsp_chromstart: trsp_chromend])

			### use lists to handle transcripts with multiple start and stop codons
			startCodonMrnaList = []
			stopCodonMrnaList = []
			
			for item in GFFlist[chrom].features[transcriptnum].sub_features:
				if trsp_strand== 1:
					if item.type== 'exon': # or item.type== 'CDS':  # For yeast, use 'CDS'
						exonstart= int(item.location.start.position)  # 0-based position
						exonend= int(item.location.end.position) # not 0-based
						exonstart_feat= exonstart- trsp_chromstart
						exonend_feat= exonend- trsp_chromstart # Not 0-based, it is fine for length....next line. 
						exonsplicedseq+= transcriptseq[exonstart_feat:exonend_feat] # takes from exonstart to exonend-1
					if item.type== 'start_codon':
						startcodonpos= item.location.start.position # 0-based position
						# startcodonmrnapos=  chrpostomrnapos(startcodonpos,chrom,transcriptnum,GFFlist)  # spliced mRNA position
						startCodonMrnaList.append(chrpostomrnapos(startcodonpos,chrom,transcriptnum,GFFlist))	# spliced mRNA position
					if item.type== 'stop_codon':
						stopcodonpos= item.location.end.position- 1 # 0-based position
						# stopcodonmrnapos= chrpostomrnapos(stopcodonpos,chrom,transcriptnum,GFFlist)
						stopCodonMrnaList.append(chrpostomrnapos(stopcodonpos,chrom,transcriptnum,GFFlist))

				if trsp_strand== -1:
					# reverse_complement() # this comes from seqIO
					transcriptseq_rev= transcriptseq.reverse_complement() 
					if item.type== 'exon': # or item.type== 'CDS':  # For yeast, use 'CDS'
						exonstart= int(item.location.start.position)  # 0-based position
						exonend= int(item.location.end.position)    # not 0-based 
						exonstart_feat= (trsp_chromend-1)- (exonend- 1)         # 0-based
						exonend_feat= (trsp_chromend-1)- exonstart      # 0-based
						exonseq= transcriptseq_rev[exonstart_feat:exonend_feat+ 1] 
						exonsplicedseq= exonseq+ exonsplicedseq
					if item.type== 'start_codon':
						startcodonpos= item.location.end.position- 1 # Need to -1 to be 0-based.
						# startcodonmrnapos= chrpostomrnapos(startcodonpos,chrom,transcriptnum,GFFlist)
						startCodonMrnaList.append(chrpostomrnapos(startcodonpos,chrom,transcriptnum,GFFlist))
					if item.type== 'stop_codon':
						stopcodonpos= item.location.start.position  # start.position is 0-based already. 
						# stopcodonmrnapos= chrpostomrnapos(stopcodonpos,chrom,transcriptnum,GFFlist)
						stopCodonMrnaList.append(chrpostomrnapos(stopcodonpos,chrom,transcriptnum,GFFlist))
			
			### choose start and stop codons			
			if len(startCodonMrnaList) > 0:
				# print "MORE THAN 1 START", startCodonMrnaList
				startcodonmrnapos = min(startCodonMrnaList)
			else:
				print "!!! no start codon for %s" % (trsp_id)
				startcodonmrnapos = 0 ### adding for transcripts without start codon
			# if len(stopCodonMrnaList)
			if len(stopCodonMrnaList) > 0:
				stopcodonmrnapos = max(stopCodonMrnaList)
			else:
				print "!!! no stop codon for %s" % (trsp_id)
				stopcodonmrnapos = len(exonsplicedseq) - 3 ### leave 3nt's in "3'UTR"

			mRNAseq = exonsplicedseq
			cdsseq= exonsplicedseq[startcodonmrnapos: stopcodonmrnapos+ 1] # take from startcodonmrnapos to stopcodonmrnapos
			utr5seq = exonsplicedseq[:startcodonmrnapos]
			utr3seq = exonsplicedseq[stopcodonmrnapos+1:]

			outseq = utr5seq.lower()+cdsseq.upper()+utr3seq.lower()

			
			if str(cdsseq[:3].upper())!= "ATG": 
				nonATGstart += 1
				# continue    # ignore non-AUG start codons
			
			### stopcodon is included in cdsseq, represnted by the last 3nt's
			stopcodon= str(cdsseq[-3:].upper())

			if stopcodon!= "TGA" and stopcodon!= "TAG" and stopcodon!= "TAA":   
				wrongstopcodon += 1
				# continue    # ignore weird stop codons

			# build itmes in transcript attribute list
			mRNAlen = len(exonsplicedseq)
			cdslen = len(cdsseq)
			utr5len = len(utr5seq)
			utr3len = len(utr3seq)
			assert mRNAlen == utr3len+cdslen+utr5len # check that sum of features equals mRNA length

			trsp_attr_list = [trsp_id, trsp_genename, outseq]
			ucscIDlist.append(trsp_attr_list[0])
			transcriptdict[trsp_id] = trsp_attr_list
			total_transcripts += 1

	print "total number of transcripts in data table: %s" % total_transcripts
	print "Number of included chromosomes chr: %s" % validchroms
	print "Number of excluded chromosomes chr: %s" % nonvalidchorms
	print "included chroms: ", included_chroms
	print "excluded chroms: ", excluded_chroms
	print "transcripts discarded due to non-AUG start codon %s" % nonATGstart
	print "transcripts discarded due to noncanonical stop codon %s" % wrongstopcodon
	return ucscIDlist, transcriptdict

def get_Prot_sequence(GFFlist):

	transcriptdict={}
	ucscIDlist = []
	total_transcripts = 0 
	nonvalidchorms = 0
	nonATGstart = 0
	wrongstopcodon = 0
	validchroms = 0
	excluded_chroms = []
	included_chroms = []
	for chrom in GFFlist:
		if not chrom in validChrs:
			excluded_chroms.append(chrom)
			nonvalidchorms += 1
			# print chrom
			continue    # check that only valid choromosomes are used
		validchroms+=1
		included_chroms.append(chrom)
		transcriptnum= -1 # set to negative one so first transcript is == to 0
		for transcript in GFFlist[chrom].features:  # this is where the SeqFeatures are actually stored
			tr_attribute_list = []
			transcriptnum+=1
			trsp_id= transcript.id # it is a number 
			trsp_strand= transcript.strand
			trsp_genename= transcript.qualifiers['Name'][0]
			trsp_chromstart= int(transcript.location.start.position)  # 0-based
			trsp_chromend= int(transcript.location.end.position) 
			transcriptlist= [0.0 for x in range(abs(trsp_chromend- trsp_chromstart))] # a list for transcript (pre-mRNA), not CDS
			
			exonsplicedseq= SeqIO.Seq('')
			transcriptseq= SeqIO.Seq(genome[chrom][trsp_chromstart: trsp_chromend])

			### use lists to handle transcripts with multiple start and stop codons
			startCodonMrnaList = []
			stopCodonMrnaList = []
			
			for item in GFFlist[chrom].features[transcriptnum].sub_features:
				if trsp_strand== 1:
					if item.type== 'exon': # or item.type== 'CDS':  # For yeast, use 'CDS'
						exonstart= int(item.location.start.position)  # 0-based position
						exonend= int(item.location.end.position) # not 0-based
						exonstart_feat= exonstart- trsp_chromstart
						exonend_feat= exonend- trsp_chromstart # Not 0-based, it is fine for length....next line. 
						exonsplicedseq+= transcriptseq[exonstart_feat:exonend_feat] # takes from exonstart to exonend-1
					if item.type== 'start_codon':
						startcodonpos= item.location.start.position # 0-based position
						# startcodonmrnapos=  chrpostomrnapos(startcodonpos,chrom,transcriptnum,GFFlist)  # spliced mRNA position
						startCodonMrnaList.append(chrpostomrnapos(startcodonpos,chrom,transcriptnum,GFFlist))	# spliced mRNA position
					if item.type== 'stop_codon':
						stopcodonpos= item.location.end.position- 1 # 0-based position
						# stopcodonmrnapos= chrpostomrnapos(stopcodonpos,chrom,transcriptnum,GFFlist)
						stopCodonMrnaList.append(chrpostomrnapos(stopcodonpos,chrom,transcriptnum,GFFlist))

				if trsp_strand== -1:
					# reverse_complement() # this comes from seqIO
					transcriptseq_rev= transcriptseq.reverse_complement() 
					if item.type== 'exon': # or item.type== 'CDS':  # For yeast, use 'CDS'
						exonstart= int(item.location.start.position)  # 0-based position
						exonend= int(item.location.end.position)    # not 0-based 
						exonstart_feat= (trsp_chromend-1)- (exonend- 1)         # 0-based
						exonend_feat= (trsp_chromend-1)- exonstart      # 0-based
						exonseq= transcriptseq_rev[exonstart_feat:exonend_feat+ 1] 
						exonsplicedseq= exonseq+ exonsplicedseq
					if item.type== 'start_codon':
						startcodonpos= item.location.end.position- 1 # Need to -1 to be 0-based.
						# startcodonmrnapos= chrpostomrnapos(startcodonpos,chrom,transcriptnum,GFFlist)
						startCodonMrnaList.append(chrpostomrnapos(startcodonpos,chrom,transcriptnum,GFFlist))
					if item.type== 'stop_codon':
						stopcodonpos= item.location.start.position  # start.position is 0-based already. 
						# stopcodonmrnapos= chrpostomrnapos(stopcodonpos,chrom,transcriptnum,GFFlist)
						stopCodonMrnaList.append(chrpostomrnapos(stopcodonpos,chrom,transcriptnum,GFFlist))

			### choose start and stop codons			
			if len(startCodonMrnaList) > 0:
				# print "MORE THAN 1 START", startCodonMrnaList
				startcodonmrnapos = min(startCodonMrnaList)
			else:
				print "!!! no start codon for %s" % (trsp_id)
			# if len(stopCodonMrnaList)
			if len(stopCodonMrnaList) > 0:
				stopcodonmrnapos = max(stopCodonMrnaList)
			else:
				print "!!! no stop codon for %s" % (trsp_id)
						
			mRNAseq = exonsplicedseq
			cdsseq= exonsplicedseq[startcodonmrnapos: stopcodonmrnapos+ 1] # take from startcodonmrnapos to stopcodonmrnapos
			utr5seq = exonsplicedseq[:startcodonmrnapos]
			utr3seq = exonsplicedseq[stopcodonmrnapos+1:]

			cdsProt = cdsseq.translate()

			# outseq = utr5seq.lower()+cdsseq.upper()+utr3seq.lower()

			
			if str(cdsseq[:3].upper())!= "ATG": 
				nonATGstart += 1
				continue    # ignore non-AUG start codons
			
			### stopcodon is included in cdsseq, represnted by the last 3nt's
			stopcodon= str(cdsseq[-3:].upper())

			if stopcodon!= "TGA" and stopcodon!= "TAG" and stopcodon!= "TAA":   
				wrongstopcodon += 1
				continue    # ignore weird stop codons

			# build itmes in transcript attribute list
			mRNAlen = len(exonsplicedseq)
			cdslen = len(cdsseq)
			utr5len = len(utr5seq)
			utr3len = len(utr3seq)
			assert mRNAlen == utr3len+cdslen+utr5len # check that sum of features equals mRNA length

			trsp_attr_list = [trsp_id, trsp_genename, cdsProt]
			ucscIDlist.append(trsp_attr_list[0])
			transcriptdict[trsp_id] = trsp_attr_list
			total_transcripts += 1
	print "total number of transcripts in data table: %s" % total_transcripts
	print "Number of included chromosomes chr: %s" % validchroms
	print "Number of excluded chromosomes chr: %s" % nonvalidchorms
	print "included chroms: ", included_chroms
	print "excluded chroms: ", excluded_chroms
	print "transcripts discarded due to non-AUG start codon %s" % nonATGstart
	print "transcripts discarded due to noncanonical stop codon %s" % wrongstopcodon
	return ucscIDlist, transcriptdict

def build_utr_table(GFFlist, include_noncanon_start, include_noncanon_stop):
	# This is a function to get the cds and utr sizes for an mRNA from a GFF file
	# returns a list with: #transcript,chrom,featnum,strand,mrna_len,cds_len,5utr_len,3utr_len,gene_name
	# Includes most of the functions from densebuilder_main but does not return counts
	# GFFlist = GFFinput

	transcriptdict={}
	ucscIDlist = []
	total_transcripts = 0 
	nonvalidchorms = 0
	nonATGstart = 0
	wrongstopcodon = 0
	validchroms = 0
	excluded_chroms = []
	included_chroms = []
	for chrom in GFFlist:
		if not chrom in validChrs:
			excluded_chroms.append(chrom)
			nonvalidchorms += 1
			# print chrom
			continue	# check that only valid choromosomes are used
		validchroms+=1
		included_chroms.append(chrom)
		transcriptnum= -1 # set to negative one so first transcript is == to 0
		for transcript in GFFlist[chrom].features:	# this is where the SeqFeatures are actually stored
			tr_attribute_list = []
			transcriptnum+=1
			trsp_id= transcript.id # it is a number 
			trsp_strand= transcript.strand
			### changing this to be compatible with new hg38 annotation
			# print transcript.qualifiers ### these are all of the fields parsed by the GTF parser from column 8, output is a dictionary {'key':['item1', 'item2', 'ect']}
			trsp_genename= transcript.qualifiers['Name'][0]
			trsp_chromstart= int(transcript.location.start.position)  # 0-based
			trsp_chromend= int(transcript.location.end.position) 
			transcriptlist= [0.0 for x in range(abs(trsp_chromend- trsp_chromstart))] # a list for transcript (pre-mRNA), not CDS
			
			exonsplicedseq= SeqIO.Seq('')
			transcriptseq= SeqIO.Seq(genome[chrom][trsp_chromstart: trsp_chromend])
			
			### use lists to handle transcripts with multiple start and stop codons
			startCodonMrnaList = []
			stopCodonMrnaList = []

			for item in GFFlist[chrom].features[transcriptnum].sub_features:
				


				if trsp_strand== 1:

					### dealing with transcripts having multiple start or stop codon entries, if spaning splice junctions


					if item.type== 'exon': # or item.type== 'CDS':	# For yeast, use 'CDS'
						exonstart= int(item.location.start.position)  # 0-based position
						exonend= int(item.location.end.position) # not 0-based
						exonstart_feat= exonstart- trsp_chromstart
						exonend_feat= exonend- trsp_chromstart # Not 0-based, it is fine for length....next line. 
						exonsplicedseq+= transcriptseq[exonstart_feat:exonend_feat] # takes from exonstart to exonend-1
					if item.type== 'start_codon':
						startcodonpos= item.location.start.position # 0-based position
						# startcodonmrnapos=  chrpostomrnapos(startcodonpos,chrom,transcriptnum,GFFlist)	# spliced mRNA position
						startCodonMrnaList.append(chrpostomrnapos(startcodonpos,chrom,transcriptnum,GFFlist))	# spliced mRNA position
						# print startcodonmrnapos
					if item.type== 'stop_codon':
						stopcodonpos= item.location.end.position- 1 # 0-based position
						# stopcodonmrnapos= chrpostomrnapos(stopcodonpos,chrom,transcriptnum,GFFlist)
						stopCodonMrnaList.append(chrpostomrnapos(stopcodonpos,chrom,transcriptnum,GFFlist))
						# print stopcodonmrnapos

				if trsp_strand== -1:
					# print 'neg_strand'
					# reverse_complement() # this comes from seqIO
					transcriptseq_rev= transcriptseq.reverse_complement() 

					if item.type== 'exon': # or item.type== 'CDS':	# For yeast, use 'CDS'
						exonstart= int(item.location.start.position)  # 0-based position
						exonend= int(item.location.end.position)	# not 0-based 
						exonstart_feat= (trsp_chromend-1)- (exonend- 1) 		# 0-based
						exonend_feat= (trsp_chromend-1)- exonstart 		# 0-based
						exonseq= transcriptseq_rev[exonstart_feat:exonend_feat+ 1] 
						exonsplicedseq= exonseq+ exonsplicedseq
					if item.type== 'start_codon':
						startcodonpos= item.location.end.position- 1 # Need to -1 to be 0-based.
						# print startcodonpos
						# startcodonmrnapos= chrpostomrnapos(startcodonpos,chrom,transcriptnum,GFFlist)
						startCodonMrnaList.append(chrpostomrnapos(startcodonpos,chrom,transcriptnum,GFFlist))
						# print "start codon: ", startcodonmrnapos
					if item.type== 'stop_codon':
						stopcodonpos= item.location.start.position	# start.position is 0-based already. 
						# print stopcodonpos
						# stopcodonmrnapos= chrpostomrnapos(stopcodonpos,chrom,transcriptnum,GFFlist)
						stopCodonMrnaList.append(chrpostomrnapos(stopcodonpos,chrom,transcriptnum,GFFlist))
						# print "stop codon: ", stopcodonmrnapos
			
			if len(startCodonMrnaList) > 0:
				# print "MORE THAN 1 START", startCodonMrnaList
				startcodonmrnapos = min(startCodonMrnaList)
			else:
				print "!!! no start codon for %s" % (trsp_id)
				startcodonmrnapos = 0 ### adding for transcripts without start codon
			# if len(stopCodonMrnaList)
			if len(stopCodonMrnaList) > 0:
				stopcodonmrnapos = max(stopCodonMrnaList)
			else:
				print "!!! no stop codon for %s" % (trsp_id)
				stopcodonmrnapos = len(exonsplicedseq) - 3 ### leave 3nt's in "3'UTR"



			cdsseq= exonsplicedseq[startcodonmrnapos: stopcodonmrnapos+ 1] # take from startcodonmrnapos to stopcodonmrnapos
			utr5seq = exonsplicedseq[:startcodonmrnapos]
			utr3seq = exonsplicedseq[stopcodonmrnapos+1:]
			
			# print trsp_id
			# # print transcript.qualifiers['transcript_name']
			# print trsp_strand
			# print utr5seq
			# print " - - - "
			# print cdsseq
			# print " - - - "
			# print utr3seq
			# # print utr5seq+cdsseq+utr3seq
			# print ""
			# # print transcriptseq

			
			if include_noncanon_start == False:
				if str(cdsseq[:3].upper())!= "ATG":	
					nonATGstart += 1
					print "non canon start"
					print trsp_id
					print cdsseq
					print ""
					continue	# ignore non-AUG start codons

			stopcodon= str(cdsseq[-3:].upper())
			if len(utr3seq) > 0:
				stop4nt = stopcodon +str(utr3seq[0].upper())
			elif len(utr3seq) == 0: 
				stop4nt = '0'
			else:
				print "there is a 3'UTR with negative length..."
				sys.exit()
			

			if include_noncanon_stop == False:
				if stopcodon!= "TGA" and stopcodon!= "TAG" and stopcodon!= "TAA":	
					wrongstopcodon += 1
					print "wrong stop!"
					print trsp_id
					print cdsseq
					print ""
					continue	# ignore weird stop codons

			# build itmes in transcript attribute list
			mRNAlen = len(exonsplicedseq)
			cdslen = len(cdsseq)
			utr5len = len(utr5seq)
			utr3len = len(utr3seq)
			assert mRNAlen == utr3len+cdslen+utr5len # check that sum of features equals mRNA length

			trsp_attr_list = [trsp_id, chrom, transcriptnum, trsp_strand, mRNAlen, cdslen, utr5len, utr3len, trsp_genename, stopcodon, stop4nt]
			ucscIDlist.append(trsp_attr_list[0])
			transcriptdict[trsp_id] = trsp_attr_list
			total_transcripts += 1
			#transcript,chrom,featnum,strand,mrna_len,cds_len,5utr_len,3utr_len,gene_name,stopcodon,stop4nt 
	print "total number of transcripts in data table: %s" % total_transcripts
	print "Number of included chromosomes chr: %s" % validchroms
	print "Number of excluded chromosomes chr: %s" % nonvalidchorms
	print "included chroms: ", included_chroms
	print "excluded chroms: ", excluded_chroms
	print "transcripts discarded due to non-AUG start codon %s" % nonATGstart
	print "transcripts discarded due to noncanonical stop codon %s" % wrongstopcodon
	return ucscIDlist, transcriptdict
	
def build_stopcodon_table(GFFlist, include_noncanon_start, include_noncanon_stop):
	# This is a function to get the cds and utr sizes for an mRNA from a GFF file
	# returns a list with: #transcript,chrom,featnum,strand,mrna_len,cds_len,5utr_len,3utr_len,gene_name
	# Includes most of the functions from densebuilder_main but does not return counts
	# GFFlist = GFFinput

	transcriptdict={}
	ucscIDlist = []
	total_transcripts = 0 
	nonvalidchorms = 0
	nonATGstart = 0
	wrongstopcodon = 0
	validchroms = 0
	excluded_chroms = []
	included_chroms = []
	for chrom in GFFlist:
		if not chrom in validChrs:
			excluded_chroms.append(chrom)
			nonvalidchorms += 1
			# print chrom
			continue    # check that only valid choromosomes are used
		validchroms+=1
		included_chroms.append(chrom)
		transcriptnum= -1 # set to negative one so first transcript is == to 0
		for transcript in GFFlist[chrom].features:  # this is where the SeqFeatures are actually stored
			### testing
			# if transcript == 'uc002zku.3':	print transcript
			###
			tr_attribute_list = []
			transcriptnum+=1
			trsp_id= transcript.id # it is a number 
			trsp_strand= transcript.strand
			trsp_genename= transcript.qualifiers['Name'][0]
			trsp_chromstart= int(transcript.location.start.position)  # 0-based
			trsp_chromend= int(transcript.location.end.position) 
			transcriptlist= [0.0 for x in range(abs(trsp_chromend- trsp_chromstart))] # a list for transcript (pre-mRNA), not CDS
			
			exonsplicedseq= SeqIO.Seq('')
			transcriptseq= SeqIO.Seq(genome[chrom][trsp_chromstart: trsp_chromend])

			### use lists to handle transcripts with multiple start and stop codons
			startCodonMrnaList = []
			stopCodonMrnaList = []
			
			for item in GFFlist[chrom].features[transcriptnum].sub_features:
				if trsp_strand== 1:
					if item.type== 'exon': # or item.type== 'CDS':  # For yeast, use 'CDS'
						exonstart= int(item.location.start.position)  # 0-based position
						exonend= int(item.location.end.position) # not 0-based
						exonstart_feat= exonstart- trsp_chromstart
						exonend_feat= exonend- trsp_chromstart # Not 0-based, it is fine for length....next line. 
						exonsplicedseq+= transcriptseq[exonstart_feat:exonend_feat] # takes from exonstart to exonend-1
					if item.type== 'start_codon':
						startcodonpos= item.location.start.position # 0-based position
						# startcodonmrnapos=  chrpostomrnapos(startcodonpos,chrom,transcriptnum,GFFlist)  # spliced mRNA position
						startCodonMrnaList.append(chrpostomrnapos(startcodonpos,chrom,transcriptnum,GFFlist))	# spliced mRNA position
					if item.type== 'stop_codon':
						stopcodonpos= item.location.end.position- 1 # 0-based position
						# stopcodonmrnapos= chrpostomrnapos(stopcodonpos,chrom,transcriptnum,GFFlist)
						stopCodonMrnaList.append(chrpostomrnapos(stopcodonpos,chrom,transcriptnum,GFFlist))

				if trsp_strand== -1:
					# reverse_complement() # this comes from seqIO
					transcriptseq_rev= transcriptseq.reverse_complement() 
					if item.type== 'exon': # or item.type== 'CDS':  # For yeast, use 'CDS'
						exonstart= int(item.location.start.position)  # 0-based position
						exonend= int(item.location.end.position)    # not 0-based 
						exonstart_feat= (trsp_chromend-1)- (exonend- 1)         # 0-based
						exonend_feat= (trsp_chromend-1)- exonstart      # 0-based
						exonseq= transcriptseq_rev[exonstart_feat:exonend_feat+ 1] 
						exonsplicedseq= exonseq+ exonsplicedseq
					if item.type== 'start_codon':
						startcodonpos= item.location.end.position- 1 # Need to -1 to be 0-based.
						# startcodonmrnapos= chrpostomrnapos(startcodonpos,chrom,transcriptnum,GFFlist)
						startCodonMrnaList.append(chrpostomrnapos(startcodonpos,chrom,transcriptnum,GFFlist))
					if item.type== 'stop_codon':
						stopcodonpos= item.location.start.position  # start.position is 0-based already. 
						# stopcodonmrnapos= chrpostomrnapos(stopcodonpos,chrom,transcriptnum,GFFlist)
						stopCodonMrnaList.append(chrpostomrnapos(stopcodonpos,chrom,transcriptnum,GFFlist))


			if len(startCodonMrnaList) > 0:
				# print "MORE THAN 1 START", startCodonMrnaList
				startcodonmrnapos = min(startCodonMrnaList)
			else:
				print "!!! no start codon for %s" % (trsp_id)
				startcodonmrnapos = 0 ### adding for transcripts without start codon
			# if len(stopCodonMrnaList)
			if len(stopCodonMrnaList) > 0:
				stopcodonmrnapos = max(stopCodonMrnaList)
			else:
				print "!!! no stop codon for %s" % (trsp_id)
				stopcodonmrnapos = len(exonsplicedseq) - 3 ### leave 3nt's in "3'UTR"

			cdsseq= exonsplicedseq[startcodonmrnapos: stopcodonmrnapos+ 1] # take from startcodonmrnapos to stopcodonmrnapos
			utr5seq = exonsplicedseq[:startcodonmrnapos]
			utr3seq = exonsplicedseq[stopcodonmrnapos+1:] ## start utr3 after the stop codon 
			
			if include_noncanon_start == False:
				if str(cdsseq[:3].upper())!= "ATG": 
					nonATGstart += 1
					continue    # ignore non-AUG start codons
			
			stopcodon= str(cdsseq[-3:].upper())
			if len(utr3seq) > 0:
				stop4nt = stopcodon +str(utr3seq[0].upper())
			elif len(utr3seq) == 0: 
				stop4nt = '0'
			else:
				print "there is a 3'UTR with negative length..."
				sys.exit()
			
			if include_noncanon_stop == False:
				if stopcodon!= "TGA" and stopcodon!= "TAG" and stopcodon!= "TAA":   
					wrongstopcodon += 1
					continue    # ignore weird stop codons

			# build itmes in transcript attribute list
			mRNAlen = len(exonsplicedseq)
			cdslen = len(cdsseq)
			utr5len = len(utr5seq)
			utr3len = len(utr3seq)
			assert mRNAlen == utr3len+cdslen+utr5len # check that sum of features equals mRNA length


			###### Finding inframe stop codons ######

			### Frame zero for loop, 
				### count each codon into 3'UTR using 0-based counting
				### With zero-based counting, next stopcodon * 3 == adjusted 3'UTR length
			frameZeroTrans = utr3seq.translate()
			frameZeroPos = -1
			frameZeroStopCounter = 0
			frameZeroUtr3LenAdj = 0
			for codon in frameZeroTrans:
				frameZeroPos +=1 
				if codon == '*':
					frameZeroStopCounter +=1
				if codon == '*' and frameZeroStopCounter == 1:
					frameZeroUtr3LenAdj = frameZeroPos*3
			if frameZeroUtr3LenAdj == 0 and frameZeroStopCounter == 0:
				frameZeroUtr3LenAdj = len(utr3seq)

			### Frame +1 for loop, 
			framePlusOneTrans = utr3seq[1:].translate() # start one nucleotide into 3'UTR for +1 frameshift
			framePlusOnePos = -1
			framePlusOneStopCounter = 0
			framePlusOneUtr3LenAdj = 0
			for codon in framePlusOneTrans:
				framePlusOnePos +=1 
				if codon == '*':
					framePlusOneStopCounter +=1
				if codon == '*' and framePlusOneStopCounter == 1:
					framePlusOneUtr3LenAdj = (framePlusOnePos*3)+1
			if framePlusOneUtr3LenAdj == 0 and framePlusOneStopCounter == 0:
				# framePlusOneUtr3LenAdj = len(utr3seq[1:])
				framePlusOneUtr3LenAdj = len(utr3seq)-1
				if framePlusOneUtr3LenAdj == -1:
					framePlusOneUtr3LenAdj = 0 ## avoid negative length utr's
			### checking for bugs with 3'utr length adjustments in the +1 frame
			if utr3len !=0 and framePlusOneStopCounter == 0 and framePlusOneUtr3LenAdj == 0:
				if utr3len !=1:
					print "oh shit, lenadj is wrong in +1 frame for %s" % trsp_id
					print "utr3len: %s" % utr3len
					print "framePlusOneStopCounter: %s" % framePlusOneStopCounter
					print "framePlusOneUtr3LenAdj: %s" % framePlusOneUtr3LenAdj
					sys.exit()

			### Frame -1 for loop, 
			frameMinusOneTrans = (cdsseq[-1]+utr3seq).translate() # include last nucleotide of cds for -1 frameshift
			frameMinusOnePos = -1
			frameMinusOneStopCounter = 0
			frameMinusOneUtr3LenAdj = 0
			for codon in frameMinusOneTrans:
				frameMinusOnePos +=1 
				if codon == '*':
					frameMinusOneStopCounter +=1
				if codon == '*' and frameMinusOneStopCounter == 1:
					frameMinusOneUtr3LenAdj = (frameMinusOnePos*3)-1
			if frameMinusOneUtr3LenAdj == 0 and frameMinusOneStopCounter == 0:
				# frameMinusOneUtr3LenAdj = len(cdsseq[-1]+utr3seq)
				frameMinusOneUtr3LenAdj = len(utr3seq)+1

			trsp_attr_list = [trsp_id, chrom, transcriptnum, trsp_strand, 
							  mRNAlen, cdslen, utr5len, utr3len, trsp_genename, stopcodon, stop4nt,
							 frameZeroStopCounter, frameZeroUtr3LenAdj,
							 framePlusOneStopCounter, framePlusOneUtr3LenAdj,
							 frameMinusOneStopCounter, frameMinusOneUtr3LenAdj]

			ucscIDlist.append(trsp_attr_list[0])
			transcriptdict[trsp_id] = trsp_attr_list
			total_transcripts += 1
			#transcript,chrom,featnum,strand,mrna_len,cds_len,5utr_len,3utr_len,gene_name,stopcodon,stop4nt 
	print "total number of transcripts in data table: %s" % total_transcripts
	print "Number of included chromosomes chr: %s" % validchroms
	print "Number of excluded chromosomes chr: %s" % nonvalidchorms
	print "included chroms: ", included_chroms
	print "excluded chroms: ", excluded_chroms
	print "transcripts discarded due to non-AUG start codon %s" % nonATGstart
	print "transcripts discarded due to noncanonical stop codon %s" % wrongstopcodon
	return ucscIDlist, transcriptdict

def build_utr3_stop_positions(GFFlist):
	# This is a function to get the cds and utr sizes for an mRNA from a GFF file
	# returns a list with: #transcript,chrom,featnum,strand,mrna_len,cds_len,5utr_len,3utr_len,gene_name
	# Includes most of the functions from densebuilder_main but does not return counts
	# GFFlist = GFFinput

	transcriptdict={}
	ucscIDlist = []
	total_transcripts = 0 
	nonvalidchorms = 0
	nonATGstart = 0
	wrongstopcodon = 0
	shortcontext = 0
	validchroms = 0
	excluded_chroms = []
	included_chroms = []
	for chrom in GFFlist:
		if not chrom in validChrs:
			excluded_chroms.append(chrom)
			nonvalidchorms += 1
			# print chrom
			continue    # check that only valid choromosomes are used
		validchroms+=1
		included_chroms.append(chrom)
		transcriptnum= -1 # set to negative one so first transcript is == to 0
		for transcript in GFFlist[chrom].features:  # this is where the SeqFeatures are actually stored
			tr_attribute_list = []
			transcriptnum+=1
			trsp_id= transcript.id # it is a number 
			trsp_strand= transcript.strand
			trsp_genename= transcript.qualifiers['Name'][0]
			trsp_chromstart= int(transcript.location.start.position)  # 0-based
			trsp_chromend= int(transcript.location.end.position) 
			transcriptlist= [0.0 for x in range(abs(trsp_chromend- trsp_chromstart))] # a list for transcript (pre-mRNA), not CDS
			
			exonsplicedseq= SeqIO.Seq('')
			transcriptseq= SeqIO.Seq(genome[chrom][trsp_chromstart: trsp_chromend])

			startCodonMrnaList = []
			stopCodonMrnaList = []
			
			for item in GFFlist[chrom].features[transcriptnum].sub_features:
				if trsp_strand== 1:
					if item.type== 'exon': # or item.type== 'CDS':  # For yeast, use 'CDS'
						exonstart= int(item.location.start.position)  # 0-based position
						exonend= int(item.location.end.position) # not 0-based
						exonstart_feat= exonstart- trsp_chromstart
						exonend_feat= exonend- trsp_chromstart # Not 0-based, it is fine for length....next line. 
						exonsplicedseq+= transcriptseq[exonstart_feat:exonend_feat] # takes from exonstart to exonend-1
					if item.type== 'start_codon':
						startcodonpos= item.location.start.position # 0-based position
						# startcodonmrnapos=  chrpostomrnapos(startcodonpos,chrom,transcriptnum,GFFlist)  # spliced mRNA position
						startCodonMrnaList.append(chrpostomrnapos(startcodonpos,chrom,transcriptnum,GFFlist))	# spliced mRNA position
					if item.type== 'stop_codon':
						stopcodonpos= item.location.end.position- 1 # 0-based position
						# stopcodonmrnapos= chrpostomrnapos(stopcodonpos,chrom,transcriptnum,GFFlist)
						stopCodonMrnaList.append(chrpostomrnapos(stopcodonpos,chrom,transcriptnum,GFFlist))

				if trsp_strand== -1:
					# reverse_complement() # this comes from seqIO
					transcriptseq_rev= transcriptseq.reverse_complement() 
					if item.type== 'exon': # or item.type== 'CDS':  # For yeast, use 'CDS'
						exonstart= int(item.location.start.position)  # 0-based position
						exonend= int(item.location.end.position)    # not 0-based 
						exonstart_feat= (trsp_chromend-1)- (exonend- 1)         # 0-based
						exonend_feat= (trsp_chromend-1)- exonstart      # 0-based
						exonseq= transcriptseq_rev[exonstart_feat:exonend_feat+ 1] 
						exonsplicedseq= exonseq+ exonsplicedseq
					if item.type== 'start_codon':
						startcodonpos= item.location.end.position- 1 # Need to -1 to be 0-based.
						# startcodonmrnapos= chrpostomrnapos(startcodonpos,chrom,transcriptnum,GFFlist)
						startCodonMrnaList.append(chrpostomrnapos(startcodonpos,chrom,transcriptnum,GFFlist))
					if item.type== 'stop_codon':
						stopcodonpos= item.location.start.position  # start.position is 0-based already. 
						# stopcodonmrnapos= chrpostomrnapos(stopcodonpos,chrom,transcriptnum,GFFlist)
						stopCodonMrnaList.append(chrpostomrnapos(stopcodonpos,chrom,transcriptnum,GFFlist))
			
			### choose start and stop codons			
			if len(startCodonMrnaList) > 0:
				# print "MORE THAN 1 START", startCodonMrnaList
				startcodonmrnapos = min(startCodonMrnaList)
			else:
				print "!!! no start codon for %s" % (trsp_id)
			# if len(stopCodonMrnaList)
			if len(stopCodonMrnaList) > 0:
				stopcodonmrnapos = max(stopCodonMrnaList)
			else:
				print "!!! no stop codon for %s" % (trsp_id)

			mRNAseq = exonsplicedseq
			cdsseq= exonsplicedseq[startcodonmrnapos: stopcodonmrnapos+ 1] # take from startcodonmrnapos to stopcodonmrnapos
			utr5seq = exonsplicedseq[:startcodonmrnapos]
			utr3seq = exonsplicedseq[stopcodonmrnapos+1:]
			
			if str(cdsseq[:3].upper())!= "ATG": 
				nonATGstart += 1
				continue    # ignore non-AUG start codons
			
			### stopcodon is included in cdsseq, represnted by the last 3nt's
			stopcodon= str(cdsseq[-3:].upper())

			if stopcodon!= "TGA" and stopcodon!= "TAG" and stopcodon!= "TAA":   
				wrongstopcodon += 1
				continue    # ignore weird stop codons

			# build itmes in transcript attribute list
			mRNAlen = len(exonsplicedseq)
			cdslen = len(cdsseq)
			utr5len = len(utr5seq)
			utr3len = len(utr3seq)
			assert mRNAlen == utr3len+cdslen+utr5len # check that sum of features equals mRNA length


			###### Finding inframe stop codons ######

			### Frame zero for loop, 
				### count each codon into 3'UTR using 0-based counting
				### With zero-based counting, next stopcodon * 3 == adjusted 3'UTR length
			frameZeroTrans = utr3seq.translate()
			frameZeroStopPositions = []
			frameZeroStopPositionsMRNA = []
			frameZeroPos = -1
			frameZeroStopCounter = 0
			frameZeroUtr3LenAdj = 0
			for codon in frameZeroTrans:
				frameZeroPos +=1 
				if codon == '*':
					frameZeroStopPositions.append(frameZeroPos*3) # get utr3position in nucleotides
					frameZeroStopPositionsMRNA.append((utr5len+cdslen)+(frameZeroPos*3))
					### check mRNA position to make sure stop codons are all valid
					sc = str(mRNAseq[(utr5len+cdslen)+(frameZeroPos*3):(utr5len+cdslen)+(frameZeroPos*3)+3].upper())
					if sc != "TAA" and sc != "TAG" and sc != "TGA":
						print "stop codon in frame 0 for %s is non correct!" % trsp_id
						print "stopcodon is: %s" % sc
						sys.exit()
					frameZeroStopCounter +=1
				if codon == '*' and frameZeroStopCounter == 1:
					frameZeroUtr3LenAdj = frameZeroPos*3
			if frameZeroUtr3LenAdj == 0 and frameZeroStopCounter == 0:
				frameZeroUtr3LenAdj = len(utr3seq)

			### Frame +1 for loop, 
			framePlusOneTrans = utr3seq[1:].translate() # start one nucleotide into 3'UTR for +1 frameshift
			framePlusOneStopPositions = []
			framePlusOneStopPositionsMRNA = []
			framePlusOnePos = -1
			framePlusOneStopCounter = 0
			framePlusOneUtr3LenAdj = 0
			for codon in framePlusOneTrans:
				framePlusOnePos +=1 
				if codon == '*':
					framePlusOneStopPositions.append((framePlusOnePos*3)+1) # get utr3position in nucleotides
					framePlusOneStopPositionsMRNA.append((utr5len+cdslen)+(framePlusOnePos*3)+1)
					### check mRNA position to make sure stop codons are all valid
					sc = str(mRNAseq[((utr5len+cdslen)+(framePlusOnePos*3)+1):((utr5len+cdslen)+(framePlusOnePos*3)+1)+3].upper())
					if sc != "TAA" and sc != "TAG" and sc != "TGA":
						print "stop codon in frame +1 for %s is non correct!" % trsp_id
						print "stopcodon is: %s" % sc
						sys.exit()
					framePlusOneStopCounter +=1
				if codon == '*' and framePlusOneStopCounter == 1:
					framePlusOneUtr3LenAdj = (framePlusOnePos*3)+1
			if framePlusOneUtr3LenAdj == 0 and frameZeroStopCounter == 0:
				framePlusOneUtr3LenAdj = len(utr3seq[1:])

			### Frame -1 for loop, 
			frameMinusOneTrans = (cdsseq[-1]+utr3seq).translate() # include last nucleotide of cds for -1 frameshift
			frameMinusOneStopPositions = []
			frameMinusOneStopPositionsMRNA = []
			frameMinusOnePos = -1
			frameMinusOneStopCounter = 0
			frameMinusOneUtr3LenAdj = 0
			for codon in frameMinusOneTrans:
				frameMinusOnePos +=1 
				if codon == '*':
					frameMinusOneStopPositions.append((frameMinusOnePos*3)-1) # get utr3position in nucleotides
					frameMinusOneStopPositionsMRNA.append((utr5len+cdslen)+(frameMinusOnePos*3)-1)
					### check mRNA position to make sure stop codons are all valid
					sc = str(mRNAseq[((utr5len+cdslen)+(frameMinusOnePos*3)-1):((utr5len+cdslen)+(frameMinusOnePos*3)-1)+3].upper())
					if sc != "TAA" and sc != "TAG" and sc != "TGA":
						print "stop codon in frame -1 for %s is non correct!" % trsp_id
						print "stopcodon is: %s" % sc
						sys.exit()
					frameMinusOneStopCounter +=1
				if codon == '*' and frameMinusOneStopCounter == 1:
					frameMinusOneUtr3LenAdj = (frameMinusOnePos*3)-1
			if frameMinusOneUtr3LenAdj == 0 and frameZeroStopCounter == 0:
				frameMinusOneUtr3LenAdj = len(cdsseq[-1]+utr3seq)




			####

			trsp_attr_list = [trsp_id, trsp_genename, 
							  frameZeroStopPositions, frameZeroStopPositionsMRNA,
							  framePlusOneStopPositions, framePlusOneStopPositionsMRNA,
							  frameMinusOneStopPositions, frameMinusOneStopPositionsMRNA]

			ucscIDlist.append(trsp_attr_list[0])
			transcriptdict[trsp_id] = trsp_attr_list
			total_transcripts += 1
			#transcript,chrom,featnum,strand,mrna_len,cds_len,5utr_len,3utr_len,gene_name,stopcodon,stop4nt 
	print "total number of transcripts in data table: %s" % total_transcripts
	print "Number of included chromosomes chr: %s" % validchroms
	print "Number of excluded chromosomes chr: %s" % nonvalidchorms
	print "included chroms: ", included_chroms
	print "excluded chroms: ", excluded_chroms
	print "transcripts discarded due to non-AUG start codon %s" % nonATGstart
	print "transcripts discarded due to noncanonical stop codon %s" % wrongstopcodon
	return ucscIDlist, transcriptdict

def find_uORFs(GFFlist, uORFtableOutfile,uORFsummaryOutfile):
	# using the same basic structure as denesbuilder_main, this function identifies all uORFs and write csv files

	startCodon = Seq('ATG')

	### build empty data frames, rows will be appended as function iterates over transcripts
	dfCols = ['trxname', 'symbol', 'strand', 'uORFCounter', 
			  'startPosition', 'cdsExtension',
			 'utr5len', 'cdslen', 'utr3len',
			 'uORFlen', 'uORFseq', 'uORFaa']
	uORFdf = pd.DataFrame(columns=dfCols)

	summaryCols = ['trxname', 'symbol', 'chr', 'tr_number', 
				   'strand', 'uORFCounter', 'cdsExtension']
	summarydf = pd.DataFrame(columns=summaryCols)

	####

	total_transcripts = 0 
	nonvalidchorms = 0
	nonATGstart = 0
	wrongstopcodon = 0
	validchroms = 0
	excluded_chroms = []
	included_chroms = []
	for chrom in GFFlist:
		if not chrom in validChrs:
			excluded_chroms.append(chrom)
			nonvalidchorms += 1
	#         print chrom
			continue    # check that only valid choromosomes are used
		validchroms+=1
		included_chroms.append(chrom)
		transcriptnum= -1 # set to negative one so first transcript is == to 0
		# print chrom
		for transcript in GFFlist[chrom].features:  # this is where the SeqFeatures are actually stored
	#         print transcript
			tr_attribute_list = []
			transcriptnum+=1
			trsp_id= transcript.id # it is a number 
			trsp_strand= transcript.strand
			trsp_genename= transcript.qualifiers['Name'][0]
			trsp_chromstart= int(transcript.location.start.position)  # 0-based
			trsp_chromend= int(transcript.location.end.position) 
			transcriptlist= [0.0 for x in range(abs(trsp_chromend- trsp_chromstart))] # a list for transcript (pre-mRNA), not CDS

			exonsplicedseq= SeqIO.Seq('')
			transcriptseq= SeqIO.Seq(genome[chrom][trsp_chromstart: trsp_chromend])

			### handling transcripts with no start or stop codon:
			# startcodonmrnapos = 'absent'
			# stopcodonmrnapos = 'absent'

			startCodonMrnaList = []
			stopCodonMrnaList = []
			
			
			for item in GFFlist[chrom].features[transcriptnum].sub_features:
				if trsp_strand== 1:
					if item.type== 'exon': # or item.type== 'CDS':  # For yeast, use 'CDS'
						exonstart= int(item.location.start.position)  # 0-based position
						exonend= int(item.location.end.position) # not 0-based
						exonstart_feat= exonstart- trsp_chromstart
						exonend_feat= exonend- trsp_chromstart # Not 0-based, it is fine for length....next line. 
						exonsplicedseq+= transcriptseq[exonstart_feat:exonend_feat] # takes from exonstart to exonend-1
					if item.type== 'start_codon':
						startcodonpos= item.location.start.position # 0-based position
						# startcodonmrnapos=  chrpostomrnapos(startcodonpos,chrom,transcriptnum,GFFlist)  # spliced mRNA position
						startCodonMrnaList.append(chrpostomrnapos(startcodonpos,chrom,transcriptnum,GFFlist))
					if item.type== 'stop_codon':
						stopcodonpos= item.location.end.position- 1 # 0-based position
						# stopcodonmrnapos= chrpostomrnapos(stopcodonpos,chrom,transcriptnum,GFFlist)
						stopCodonMrnaList.append(chrpostomrnapos(stopcodonpos,chrom,transcriptnum,GFFlist))

				if trsp_strand== -1:
					# reverse_complement() # this comes from seqIO
					transcriptseq_rev= transcriptseq.reverse_complement() 
					if item.type== 'exon': # or item.type== 'CDS':  # For yeast, use 'CDS'
						exonstart= int(item.location.start.position)  # 0-based position
						exonend= int(item.location.end.position)    # not 0-based 
						exonstart_feat= (trsp_chromend-1)- (exonend- 1)         # 0-based
						exonend_feat= (trsp_chromend-1)- exonstart      # 0-based
						exonseq= transcriptseq_rev[exonstart_feat:exonend_feat+ 1] 
						exonsplicedseq= exonseq+ exonsplicedseq
					if item.type== 'start_codon':
						startcodonpos= item.location.end.position- 1 # Need to -1 to be 0-based.
						# startcodonmrnapos= chrpostomrnapos(startcodonpos,chrom,transcriptnum,GFFlist)
						startCodonMrnaList.append(chrpostomrnapos(startcodonpos,chrom,transcriptnum,GFFlist))
					if item.type== 'stop_codon':
						stopcodonpos= item.location.start.position  # start.position is 0-based already. 
						# stopcodonmrnapos= chrpostomrnapos(stopcodonpos,chrom,transcriptnum,GFFlist)
						stopCodonMrnaList.append(chrpostomrnapos(stopcodonpos,chrom,transcriptnum,GFFlist))
			

			if len(startCodonMrnaList) > 0:
				# print "MORE THAN 1 START", startCodonMrnaList
				startcodonmrnapos = min(startCodonMrnaList)
			else:
				print "!!! no start codon for %s" % (trsp_id)
			# if len(stopCodonMrnaList)
			if len(stopCodonMrnaList) > 0:
				stopcodonmrnapos = max(stopCodonMrnaList)
			else:
				print "!!! no stop codon for %s" % (trsp_id)


	# 		if startcodonmrnapos == 'absent' or stopcodonmrnapos == 'absent':
	# #             print "no start of stop for trsp  %s" % trsp_id
	# 			continue
			
			cdsseq= exonsplicedseq[startcodonmrnapos: stopcodonmrnapos+ 1] # take from startcodonmrnapos to stopcodonmrnapos
			utr5seq = exonsplicedseq[:startcodonmrnapos]
			utr3seq = exonsplicedseq[stopcodonmrnapos+1:]

			if str(cdsseq[:3].upper())!= "ATG": 
				nonATGstart += 1
				continue    # ignore non-AUG start codons
			stopcodon= str(cdsseq[-3:].upper())
			# if len(utr3seq) > 0:
			# 	stop4nt = stopcodon +str(utr3seq[0].upper())
			# elif len(utr3seq) == 0: 
			# 	stop4nt = '0'
			# else:
			# 	print "there is a 3'UTR with negative length..."
			# 	sys.exit()
			if stopcodon!= "TGA" and stopcodon!= "TAG" and stopcodon!= "TAA":   
				wrongstopcodon += 1
				continue    # ignore weird stop codons

			# build itmes in transcript attribute list
			mRNAlen = len(exonsplicedseq)
			cdslen = len(cdsseq)
			utr5len = len(utr5seq)
			utr3len = len(utr3seq)
			assert mRNAlen == utr3len+cdslen+utr5len # check that sum of features equals mRNA length

			#### Counting of uORFs ####
			uORFcounter = 0
			cdsExtension = 0

			for i in range(len(utr5seq)):
				### iterate over every nucleotide in the 5'UTR

				codon = utr5seq[i:i+3] # define the codon at each position
				if str(codon) == str(startCodon): # check if it is a start codon
					uORFcounter += 1

					startPosition = i
					seqIndex = i
					uORFaa = []
					uORFseq = []
					# print "found start codon at pos %s" % startPosition
					aminoAcid = codon.translate()
					uORFseq.append(str(codon))
					uORFaa.append(str(aminoAcid))
					
					while str(aminoAcid) != "*": # continue this loop until a stop codon is encoutered
						seqIndex +=3 # advance by 3 nt's each time (1 codon)
						nextCodon = utr5seq[seqIndex:seqIndex+3]
						aminoAcid = nextCodon.translate()
						if len(nextCodon) == 3: # ensure that a full codon is still present, do not want 1 or 2 nts
							uORFseq.append(str(nextCodon))
							uORFaa.append(str(aminoAcid))

						if seqIndex > len(utr5seq)-2: # if uORF continues into cds, retreive sequences from here
							# -2 is because this will not yeild a full codon (only 2 nt's)
							# print "end of UTR"
							cdsExtension = 1

							utrCdsSeq = utr5seq+cdsseq
							nextCodon = utrCdsSeq[seqIndex:seqIndex+3]
							aminoAcid = nextCodon.translate()
							uORFseq.append(str(nextCodon))
							uORFaa.append(str(aminoAcid))
							# print nextCodon, aminoAcid

							if seqIndex > len(utrCdsSeq): ## if uORF exceeds coding region, stop counting this,
								### could eventually extend to the 3'UTR if any transcript exists here
								print 'end of CDS for trsp %s' % trsp_id
								break

					uORFseqCat = "".join(uORFseq) # remove seperate list entries and concat to a string
					uORFaaCat = "".join(uORFaa)

					### save all uORF features to a list, and build into a dataframe
					uORF_features = [trsp_id, trsp_genename, trsp_strand, uORFcounter, startPosition, 
									 cdsExtension, len(utr5seq), len(cdsseq), len(utr3seq),
									 len(uORFseqCat), uORFseqCat, uORFaaCat]

					dftemp = pd.DataFrame([uORF_features], columns=dfCols) ## 
					# print dftemp

					uORFdf = pd.concat([uORFdf, dftemp], ignore_index=True)

				if i == (len(utr5seq)-1): # at the end of the 5'UTR, do this ...
					# print i
					uORFsummary = [trsp_id, trsp_genename, chrom, transcriptnum, trsp_strand, uORFcounter, cdsExtension]
					# print uORFsummary
					dfSummaryTemp = pd.DataFrame([uORFsummary], columns=summaryCols)
					# print dfSummaryTemp
					summarydf = pd.concat([summarydf, dfSummaryTemp], ignore_index=True)
		
	uORFdf.to_csv(uORFtableOutfile)
	summarydf.to_csv(uORFsummaryOutfile)
	print summarydf.head()

def find_cds_seq(dfin):

	df = dfin.copy()
	
	for tr in df.index:
		mrnaseq = df.loc[tr]['mRNAseqs']

		utr5len = df.loc[tr]['5utr_len']
		cdslen = df.loc[tr]['cds_len']

		cdsStart = utr5len
		cdsEnd = utr5len + cdslen

		cdsseq = mrnaseq[cdsStart:cdsEnd]

		df.loc[tr, 'cdsSeq'] = cdsseq

	return df
	

def find_codon_positions(df, codonList):
	### iterate through codons one at a time

	for cod in codonList:
		
		print cod
		csv_outpath = "%s/%s_%s_1.csv" % (csvOutDir, genAnName, cod)
		print csv_outpath
		dfcolumns = ['headers', 'gene', 'chrom', 'trsp_num', 'cds_pos']
		tempdf = pd.DataFrame(columns = dfcolumns)
		
		for tr in df.index:
			seq = df.loc[tr, 'cdsSeq']
			
			counter = 0
			while counter < len(seq):
				codon = seq[counter:counter+3]

				if codon == cod:
					
					hdr = df.loc[tr, '#transcript']
					gene = df.loc[tr, 'genename']
					chrom = df.loc[tr, 'chrom']
					trspNum = df.loc[tr, 'featnum']
					pos = counter
					
					entryList = [[hdr, gene, chrom, trspNum, pos]]
					entryDf = pd.DataFrame(entryList, columns = dfcolumns)
					
					tempdf = tempdf.append(entryDf, ignore_index=True)


				counter +=3
				
		print tempdf
		tempdf.to_csv(csv_outpath, index=False)
		
def find_codon_positions_multi(cod, df):
	### iterate through codons one at a time

	
	print cod
	csv_outpath = "%s/%s_%s_1.csv" % (csvOutDir, genAnName, cod)
	print csv_outpath
	dfcolumns = ['headers', 'gene', 'chrom', 'trsp_num', 'cds_pos']
	tempdf = pd.DataFrame(columns = dfcolumns)
	
	for tr in df.index:
		seq = df.loc[tr, 'cdsSeq']
		
		counter = 0
		while counter < len(seq):
			codon = seq[counter:counter+3]

			if codon == cod:
				
				hdr = df.loc[tr, '#transcript']
				gene = df.loc[tr, 'genename']
				chrom = df.loc[tr, 'chrom']
				trspNum = df.loc[tr, 'featnum']
				pos = counter
				
				entryList = [[hdr, gene, chrom, trspNum, pos]]
				entryDf = pd.DataFrame(entryList, columns = dfcolumns)
				
				tempdf = tempdf.append(entryDf, ignore_index=True)


			counter +=3
			
	# print tempdf
	tempdf.to_csv(csv_outpath, index=False)

def write_utr_stopcodon_csvfile(ucscIDlist, transcriptdict, outfilestring, headers):
	# For writing dictionary to csv file first headers are written Then one line at a time is added to t, This is a list that will hold all of the trsp_attr_lists 
	# 	Finally these are written to each line of a csv file
	t=[]
	t.append(headers)

	for i in ucscIDlist: # position starts at 0
		newline= transcriptdict[i]
		t.append(newline)

	fa = open(outfilestring+".csv", "w")
	writer = csv.writer(fa)
	writer.writerows(t)
	fa.close()

def main():
	GTFgen = GFF.parse(GTFfile)
	GFFlist = makeGFFlist(GTFgen)
	outfilestring = '%s_mRNAseqs' % (gtfInFilePrefix)
	if not os.path.isfile(outfilestring+".csv"):
		ucscIDlist, transcriptdict = get_mRNA_sequence(GFFlist)
		headers= ['#transcript', 'genename', 'mRNAseqs']
		write_utr_stopcodon_csvfile(ucscIDlist, transcriptdict, outfilestring, headers)
	
	outfilestring = '%s_Protseqs' % (gtfInFilePrefix)
	if not os.path.isfile(outfilestring+".csv"):	
		ucscIDlist, transcriptdict = get_Prot_sequence(GFFlist)
		headers= ['#transcript', 'genename', 'cdsProt']
		write_utr_stopcodon_csvfile(ucscIDlist, transcriptdict, outfilestring, headers)
	
	outfilestring = '%s_UTRs' % (gtfInFilePrefix)
	if not os.path.isfile(outfilestring+".csv"):
		include_noncanon_start = True
		include_noncanon_stop = True
		headers= ['#transcript','chrom','featnum','strand','mrna_len','cds_len','5utr_len','3utr_len','gene_name','stopcodon','stop4nt']
		ucscIDlist, transcriptdict = build_utr_table(GFFlist, include_noncanon_start, include_noncanon_stop)
		write_utr_stopcodon_csvfile(ucscIDlist, transcriptdict, outfilestring, headers)
		
	outfilestring = '%s_stopcodons' % (gtfInFilePrefix)
	if not os.path.isfile(outfilestring+".csv"):
		include_noncanon_start = True
		include_noncanon_stop = True
		headers= ['#transcript','chrom','featnum','strand','mrna_len','cds_len','5utr_len','3utr_len','gene_name','stopcodon','stop4nt','frameZeroStopCount', 'frameZeroUtr3LenAdj','framePlusOneStopCount', 'framePlusOneUtr3LenAdj','frameMinusOneStopCount', 'frameMinusOneUtr3LenAdj']
		ucscIDlist, transcriptdict = build_stopcodon_table(GFFlist, include_noncanon_start, include_noncanon_stop)
		write_utr_stopcodon_csvfile(ucscIDlist, transcriptdict, outfilestring, headers)
	
	outfilestring = '%s_utr3StopPositions' % (gtfInFilePrefix)
	if not os.path.isfile(outfilestring+".csv"):
		headers= ['#transcript', 'genename', 'frameZeroStopPosUtr3', 'frameZeroStopPosMRNA', 'framePlusOneStopPosUtr3', 'framePlusOneStopPosMRNA', 'frameMinusOneStopPosUtr3', 'frameMinusOneStopPosMRNA']
		ucscIDlist, transcriptdict = build_utr3_stop_positions(GFFlist)
		write_utr_stopcodon_csvfile(ucscIDlist, transcriptdict, outfilestring, headers)
		
	uORFtableOutfile = '%s_uORFtable.csv' % (gtfInFilePrefix)
	uORFsummaryOutfile = '%s_uORFsummary.csv' % (gtfInFilePrefix)
	if not os.path.isfile(uORFtableOutfile):
		find_uORFs(GFFlist,uORFtableOutfile,uORFsummaryOutfile)
	
	if not os.path.exists(csvOutDir):
		os.makedirs(csvOutDir)
		UTRfilepath = "%s_UTRs.csv" % gtfInFilePrefix
		mRNAfilepath = "%s_mRNAseqs.csv" % gtfInFilePrefix
		utrdf = pd.read_csv(UTRfilepath)
		mRNAdf = pd.read_csv(mRNAfilepath)
		
		dfref = utrdf.merge(right = mRNAdf, on='#transcript')
		df = find_cds_seq(dfref)
		dfl = [df]*len(codonList) ### replicate the dateframe as a list for the length of the input codon list
		threadNumb = 40
		p = Pool(nodes=int(threadNumb))
		p.map(find_codon_positions_multi, codonList, dfl)

	

if __name__ == '__main__':
	# execute only if run as a script
	main()
"""
}


process Ribosome_Profiling_Build_Index_Validate_Genome_Gtf {

input:
 val "*" from g40_2_valid_gtf_g40_7
 val "*" from g40_4_valid_gtf_g40_7
 val "*" from g40_6_ncRNAFastaPath_g40_7

output:
 val "complete"  into g40_7_complete_g_41

"""
echo "Building Completed."
"""
}

g40_7_complete_g_41= g40_7_complete_g_41.ifEmpty([""]) 

params.ncRNA_star_index =  ""  //* @input
if (!((params.run_ncRNA_Removal && (params.run_ncRNA_Removal == "yes")) || !params.run_ncRNA_Removal)){
g0_20_reads_g_41.into{g_41_reads_g42_32}
g_41_logOut_g_29 = Channel.empty()
} else {


process ncRNA_Removal {

input:
 set val(name), file(reads) from g0_20_reads_g_41
 val mate from g_1_mate_g_41
 val "*" from g40_7_complete_g_41

output:
 set val(name), file("*_no_ncRNA.fastq")  into g_41_reads_g42_32
 set val(name), file("${newName}Log.final.out")  into g_41_logOut_g_29

errorStrategy 'retry'

when:
(params.run_ncRNA_Removal && (params.run_ncRNA_Removal == "yes")) || !params.run_ncRNA_Removal

script:
params_STAR = params.ncRNA_Removal.params_STAR
nameAll = reads.toString()
nameArray = nameAll.split(' ')

if (nameAll.contains('.gz')) {
    newName =  nameArray[0] - ~/(\.fastq.gz)?(\.fq.gz)?$/
    file =  nameAll - '.gz' - '.gz'
    runGzip = "ls *.gz | xargs -i echo gzip -df {} | sh"
} else {
    newName =  nameArray[0] - ~/(\.fastq)?(\.fq)?$/
    file =  nameAll 
    runGzip = ''
}

"""
$runGzip
STAR ${params_STAR}  --genomeDir ${params.ncRNA_star_index} --readFilesIn $file  --outReadsUnmapped Fastx  --outSAMtype BAM SortedByCoordinate     --outFileNamePrefix ${newName}
echo "Alignment completed."
mv ${newName}Unmapped.out.mate1 ${newName}_no_ncRNA.fastq
"""


}
}



process riboseq_ncRNA_Removal_Summary {

input:
 set val(name), file(alignSum) from g_41_logOut_g_29.groupTuple()

output:
 file "*.tsv"  into g_29_outputFile_g_28
 val "ncRNA_removal_sum"  into g_29_name_g_28

shell:
'''
#!/usr/bin/env perl
use List::Util qw[min max];
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage; 
use Data::Dumper;

my %tsv;
my @headers = ();
my $name = "!{name}";

alteredAligned();

my @keys = keys %tsv;
my $summary = "$name"."_ncRNA_removal_sum.tsv";
my $header_string = join("\\t", @headers);
`echo "$header_string" > $summary`;
foreach my $key (@keys){
	my $values = join("\\t", @{ $tsv{$key} });
	`echo "$values" >> $summary`;
}


sub alteredAligned
{
	my @files = qw(!{alignSum});
	my $multimappedSum;
	my $alignedSum;
	my $inputCountSum;
	my $afterFiltering;
	push(@headers, "Sample");
    push(@headers, "Total Reads");
	push(@headers, "Reads After ncRNA Removal");
	foreach my $file (@files){
		my $multimapped;
		my $aligned;
		my $inputCount;
		chomp($inputCount = `cat $file | grep 'Number of input reads' | awk '{sum+=\\$6} END {print sum}'`);
		chomp($aligned = `cat $file | grep 'Uniquely mapped reads number' | awk '{sum+=\\$6} END {print sum}'`);
		chomp($multimapped = `cat $file | grep 'Number of reads mapped to multiple loci' | awk '{sum+=\\$9} END {print sum}'`);
		$multimappedSum += int($multimapped);
        $alignedSum += int($aligned);
        $inputCountSum += int($inputCount);
	}
	$afterFiltering = $inputCountSum - $alignedSum - $multimappedSum;
	$tsv{$name} = [$name, $inputCountSum];
	push(@{$tsv{$name}}, $afterFiltering);
}
'''

}


process Merge_TSV_Files_ncRNA_removal {

input:
 file tsv from g_29_outputFile_g_28.collect()
 val outputFileName from g_29_name_g_28.collect()

output:
 file "${name}.tsv"  into g_28_outputFileTSV_g_26

script:
name = outputFileName[0]
"""    
awk 'FNR==1 && NR!=1 {  getline; } 1 {print} ' *.tsv > ${name}.tsv
"""
}

params.run_Sequential_Mapping =   "yes"   //* @dropdown @options:"yes","no" @show_settings:"Sequential_Mapping" @description:"Filters out or quantify given sequence sets."
params.bowtieInd_rRNA =  ""  //* @input
params.bowtieInd_ercc =  ""  //* @input
params.bowtieInd_miRNA =  ""  //* @input
params.bowtieInd_tRNA =  ""  //* @input
params.bowtieInd_piRNA =  ""  //* @input
params.bowtieInd_snRNA =  ""  //* @input
params.bowtieInd_rmsk =  ""  //* @input
params.bowtie_index =  ""  //* @input
params.bowtie2_index =  ""  //* @input
params.star_index =  ""  //* @input

//both bowtie and bowtie2 indexes located in same path
bowtieIndexes = [rRNA: params.bowtieInd_rRNA, 
                 ercc: params.bowtieInd_ercc,
                 miRNA: params.bowtieInd_miRNA,
                 tRNA: params.bowtieInd_tRNA,
                 piRNA: params.bowtieInd_piRNA,
                 snRNA: params.bowtieInd_snRNA,
                 rmsk: params.bowtieInd_rmsk]
                 
genomeIndexes = [bowtie: params.bowtie_index,
                 bowtie2: params.bowtie2_index,
                 STAR: params.star_index+"/genome"]


//_nucleicAcidType="dna" should be defined in the autofill section of pipeline header in case dna is used.
_select_sequence = params.Sequential_Mapping_Module_Sequential_Mapping._select_sequence
index_directory = params.Sequential_Mapping_Module_Sequential_Mapping.index_directory
name_of_the_index_file = params.Sequential_Mapping_Module_Sequential_Mapping.name_of_the_index_file
_aligner = params.Sequential_Mapping_Module_Sequential_Mapping._aligner
aligner_Parameters = params.Sequential_Mapping_Module_Sequential_Mapping.aligner_Parameters
description = params.Sequential_Mapping_Module_Sequential_Mapping.description
filter_Out = params.Sequential_Mapping_Module_Sequential_Mapping.filter_Out

desc_all=[]
description.eachWithIndex() {param,i -> 
    if (param.isEmpty()){
        desc_all[i] = name_of_the_index_file[i]
    }  else {
        desc_all[i] = param.replaceAll("[ |.|;]", "_")
    }
}
custom_index=[]
index_directory.eachWithIndex() {param,i -> 
    if (_select_sequence[i] == "genome"){
        custom_index[i] = genomeIndexes[_aligner[i]]
    }else if (_select_sequence[i] == "custom"){
        custom_index[i] = param+"/"+name_of_the_index_file[i]
    }else {
        custom_index[i] = bowtieIndexes[_select_sequence[i]]
    }
}

mapList = []
paramList = []
alignerList = []
filterList = []
indexList = []

//concat default mapping and custom mapping
mapList = (desc_all) 
paramList = (aligner_Parameters)
alignerList = (_aligner)
filterList = (filter_Out)
indexList = (custom_index)

mappingList = mapList.join(" ") // convert into space separated format in order to use in bash for loop
paramsList = paramList.join(",") // convert into comma separated format in order to use in as array in bash
alignersList = alignerList.join(",") 
filtersList = filterList.join(",") 
indexesList = indexList.join(",") 
//* @style @condition:{remove_duplicates="yes",remove_duplicates_based_on_UMI_after_mapping},{remove_duplicates="no"},{_select_sequence="custom", index_directory,name_of_the_index_file,description,_aligner,aligner_Parameters,filter_Out},{_select_sequence=("rRNA","ercc","miRNA","tRNA","piRNA","snRNA","rmsk","genome"),_aligner,aligner_Parameters,filter_Out}  @array:{_select_sequence,_select_sequence, index_directory,name_of_the_index_file,_aligner,aligner_Parameters,filter_Out,description} @multicolumn:{_select_sequence,_select_sequence,index_directory,name_of_the_index_file,_aligner,aligner_Parameters,filter_Out, description},{remove_duplicates,remove_duplicates_based_on_UMI_after_mapping}


//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 4
    $MEMORY = 20
}
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 2000
    $CPU  = 4
    $MEMORY = 20
    $QUEUE = "long"
}
//* platform
//* autofill
if (!(params.run_Sequential_Mapping == "yes")){
g_41_reads_g42_32.into{g42_32_reads_g_4}
g42_32_bowfiles_g42_26 = Channel.empty()
g42_32_bam_file_g42_23 = Channel.empty()
g42_32_bam_file_g42_27 = Channel.empty()
g42_32_bam_index_g42_23 = Channel.empty()
g42_32_bam_index_g42_27 = Channel.empty()
g42_32_filter_g42_26 = Channel.empty()
g42_32_log_file_g42_30 = Channel.empty()
} else {


process Sequential_Mapping_Module_Sequential_Mapping {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*\/.*_sorted.bam$/) "sequential_mapping/$filename"
	else if (filename =~ /.*\/.*_sorted.bam.bai$/) "sequential_mapping/$filename"
}

input:
 set val(name), file(reads) from g_41_reads_g42_32
 val mate from g_1_mate_g42_32

output:
 set val(name), file("final_reads/*q")  into g42_32_reads_g_4
 set val(name), file("bowfiles/?*") optional true  into g42_32_bowfiles_g42_26
 file "*/*_sorted.bam" optional true  into g42_32_bam_file_g42_23
 file "*/*_sorted.bam.bai" optional true  into g42_32_bam_index_g42_23
 val filtersList  into g42_32_filter_g42_26
 file "*/*_sorted.dedup.bam" optional true  into g42_32_bam_file_g42_27
 file "*/*_sorted.dedup.bam.bai" optional true  into g42_32_bam_index_g42_27
 file "*/*_duplicates_stats.log" optional true  into g42_32_log_file_g42_30

errorStrategy 'retry'

when:
params.run_Sequential_Mapping == "yes"

script:
nameAll = reads.toString()
nameArray = nameAll.split(' ')
def file2;

if (nameAll.contains('.gz')) {
    newName =  nameArray[0] - ~/(\.fastq.gz)?(\.fq.gz)?$/
    file1 =  nameArray[0] - '.gz' 
    if (mate == "pair") {file2 =  nameArray[1] - '.gz'}
    runGzip = "ls *.gz | xargs -i echo gzip -df {} | sh"
} else {
    newName =  nameArray[0] - ~/(\.fastq)?(\.fq)?$/
    file1 =  nameArray[0]
    if (mate == "pair") {file2 =  nameArray[1]}
    runGzip = ''
}

remove_duplicates = params.Sequential_Mapping_Module_Sequential_Mapping.remove_duplicates
remove_duplicates_based_on_UMI_after_mapping = params.Sequential_Mapping_Module_Sequential_Mapping.remove_duplicates_based_on_UMI_after_mapping
remove_previous_reads = params.Sequential_Mapping_Module_Sequential_Mapping.remove_previous_reads

"""
#!/bin/bash
mkdir reads final_reads bowfiles
workflowWorkDir=\$(cd ../../ && pwd)
if [ -n "${mappingList}" ]; then
    $runGzip
    #rename files to standart format
    if [ "${mate}" == "pair" ]; then
        mv $file1 ${name}.1.fastq 2>/dev/null
        mv $file2 ${name}.2.fastq 2>/dev/null
        mv ${name}.1.fastq ${name}.2.fastq reads/.
    else
        mv $file1 ${name}.fastq 2>/dev/null
        mv ${name}.fastq reads/.
    fi
    #sequential mapping
    k=0
    prev="reads"
    IFS=',' read -r -a paramsListAr <<< "${paramsList}" #create comma separated array 
    IFS=',' read -r -a filtersListAr <<< "${filtersList}"
    IFS=',' read -r -a indexesListAr <<< "${indexesList}"
    IFS=',' read -r -a alignersListAr <<< "${alignersList}"
    wrkDir=\$(pwd)
    for rna_set in ${mappingList}
    do
        ((k++))
        printf -v k2 "%02d" "\$k" #turn into two digit format
        mkdir -p \${rna_set}/unmapped
        cd \$rna_set
        ## create link of the target file to prevent "too many symlinks error"
        for r in \${wrkDir}/\${prev}/*; do
            targetRead=\$(readlink -e \$r)
            rname=\$(basename \$r)
            echo "INFO: ln -s \$targetRead \$rname"
            ln -s \$targetRead \$rname
        done
        genomeDir=`dirname "\${indexesListAr[\$k-1]}"`
        echo "INFO: genomeDir: \$genomeDir"
        if [ -e "\${indexesListAr[\$k-1]}.1.bt2" -o  -e "\${indexesListAr[\$k-1]}.fa"  -o  -e "\${indexesListAr[\$k-1]}.fasta"  -o  -e "\$genomeDir/SAindex" ]; then
            if [ -e "\${indexesListAr[\$k-1]}.fa" ] ; then
                fasta=\${indexesListAr[\$k-1]}.fa
            elif [ -e "\${indexesListAr[\$k-1]}.fasta" ] ; then
                fasta=\${indexesListAr[\$k-1]}.fasta
            fi
            echo "INFO: fasta: \$fasta"
            if [ -e "\${indexesListAr[\$k-1]}.1.bt2" -a "\${alignersListAr[\$k-1]}" == "bowtie2" ] ; then
                echo "INFO: \${indexesListAr[\$k-1]}.1.bt2 Bowtie2 index found."
            elif [ -e "\${indexesListAr[\$k-1]}.1.ebwt" -a "\${alignersListAr[\$k-1]}" == "bowtie" ] ; then
                echo "INFO: \${indexesListAr[\$k-1]}.1.ebwt Bowtie index found."
            elif [ -e "\$genomeDir/SAindex" -a "\${alignersListAr[\$k-1]}" == "STAR" ] ; then
                echo "INFO: \$genomeDir/SAindex STAR index found."
            elif [ -e "\${indexesListAr[\$k-1]}.fa" -o  -e "\${indexesListAr[\$k-1]}.fasta" ] ; then
                if [ "\${alignersListAr[\$k-1]}" == "bowtie2" ]; then
                    bowtie2-build \$fasta \${indexesListAr[\$k-1]}
                elif [ "\${alignersListAr[\$k-1]}" == "STAR" ]; then
                    if [ -e "\${indexesListAr[\$k-1]}.gtf" ]; then
                        STAR --runMode genomeGenerate --genomeDir \$genomeDir --genomeFastaFiles \$fasta --sjdbGTFfile \${indexesListAr[\$k-1]}.gtf --genomeSAindexNbases 5
                    else
                        echo "WARNING: \${indexesListAr[\$k-1]}.gtf not found. STAR index is not generated."
                    fi
                elif [ "\${alignersListAr[\$k-1]}" == "bowtie" ]; then
                    bowtie-build \$fasta \${indexesListAr[\$k-1]}
                fi
            fi
                
            if [ "${mate}" == "pair" ]; then
                if [ "\${alignersListAr[\$k-1]}" == "bowtie2" ]; then
                    bowtie2 \${paramsListAr[\$k-1]} -x \${indexesListAr[\$k-1]} --no-unal --un-conc unmapped/${name}.unmapped.fastq -1 ${name}.1.fastq -2 ${name}.2.fastq --al-conc ${name}.fq.mapped -S \${rna_set}_${name}_alignment.sam 2>&1 | tee \${k2}_${name}.bow_\${rna_set}
                elif [ "\${alignersListAr[\$k-1]}" == "STAR" ]; then
                    STAR \${paramsListAr[\$k-1]}  --genomeDir \$genomeDir --readFilesIn ${name}.1.fastq ${name}.2.fastq --outSAMtype SAM  --outFileNamePrefix ${name}.star --outReadsUnmapped Fastx
                    mv ${name}.starAligned.out.sam \${rna_set}_${name}_alignment.sam
                    mv ${name}.starUnmapped.out.mate1 unmapped/${name}.unmapped.1.fastq
                    mv ${name}.starUnmapped.out.mate2 unmapped/${name}.unmapped.2.fastq
                    mv ${name}.starLog.final.out \${k2}_${name}.star_\${rna_set}
                elif [ "\${alignersListAr[\$k-1]}" == "bowtie" ]; then
                    bowtie \${paramsListAr[\$k-1]}   \${indexesListAr[\$k-1]}  --un  unmapped/${name}.unmapped.fastq -1 ${name}.1.fastq -2 ${name}.2.fastq -S  \${rna_set}_${name}_alignment.sam 2>&1 | tee \${k2}_${name}.bow1_\${rna_set}  
                    mv unmapped/${name}.unmapped_1.fastq unmapped/${name}.unmapped.1.fastq
                    mv unmapped/${name}.unmapped_2.fastq unmapped/${name}.unmapped.2.fastq
                fi
            else
                if [ "\${alignersListAr[\$k-1]}" == "bowtie2" ]; then
                    bowtie2 \${paramsListAr[\$k-1]} -x \${indexesListAr[\$k-1]} --no-unal --un  unmapped/${name}.unmapped.fastq -U ${name}.fastq --al ${name}.fq.mapped -S \${rna_set}_${name}_alignment.sam 2>&1 | tee \${k2}_${name}.bow_\${rna_set}  
                elif [ "\${alignersListAr[\$k-1]}" == "STAR" ]; then
                    STAR \${paramsListAr[\$k-1]}  --genomeDir \$genomeDir --readFilesIn ${name}.fastq --outSAMtype SAM  --outFileNamePrefix ${name}.star --outReadsUnmapped Fastx
                    mv ${name}.starAligned.out.sam \${rna_set}_${name}_alignment.sam
                    mv ${name}.starUnmapped.out.mate1 unmapped/${name}.unmapped.fastq
                    mv ${name}.starLog.final.out \${k2}_${name}.star_\${rna_set}
                elif [ "\${alignersListAr[\$k-1]}" == "bowtie" ]; then
                    bowtie \${paramsListAr[\$k-1]}  \${indexesListAr[\$k-1]}  --un  unmapped/${name}.unmapped.fastq  ${name}.fastq  -S \${rna_set}_${name}_alignment.sam 2>&1 | tee \${k2}_${name}.bow1_\${rna_set}  
                    
                fi
            fi
            echo "INFO: samtools view -bT \${fasta} \${rna_set}_${name}_alignment.sam > \${rna_set}_${name}_alignment.bam"
            samtools view -bT \${fasta} \${rna_set}_${name}_alignment.sam > \${rna_set}_${name}_alignment.bam
            rm -f \${rna_set}_${name}_alignment.sam
            if [ "\${alignersListAr[\$k-1]}" == "bowtie" ]; then
                mv \${rna_set}_${name}_alignment.bam \${rna_set}_${name}_tmp0.bam
                echo "INFO: samtools view -F 0x04 -b \${rna_set}_${name}_tmp0.bam > \${rna_set}_${name}_alignment.bam"
                samtools view -F 0x04 -b \${rna_set}_${name}_tmp0.bam > \${rna_set}_${name}_alignment.bam  # Remove unmapped reads
                if [ "${mate}" == "pair" ]; then
                    echo "# unique mapped reads: \$(samtools view -f 0x40 -F 0x4 -q 255 \${rna_set}_${name}_alignment.bam | cut -f 1 | sort | uniq | wc -l)" >> \${k2}_${name}.bow1_\${rna_set}
                else
                    echo "# unique mapped reads: \$(samtools view -F 0x40 -q 255 \${rna_set}_${name}_alignment.bam | cut -f 1 | sort | uniq | wc -l)" >> \${k2}_${name}.bow1_\${rna_set}
                fi
            fi
            if [ "${mate}" == "pair" ]; then
                mv \${rna_set}_${name}_alignment.bam \${rna_set}_${name}_alignment.tmp1.bam
                echo "INFO: samtools sort -n -o \${rna_set}_${name}_alignment.tmp2 \${rna_set}_${name}_alignment.tmp1.bam"
                samtools sort -n -o \${rna_set}_${name}_alignment.tmp2.bam \${rna_set}_${name}_alignment.tmp1.bam 
                echo "INFO: samtools view -bf 0x02 \${rna_set}_${name}_alignment.tmp2.bam >\${rna_set}_${name}_alignment.bam"
                samtools view -bf 0x02 \${rna_set}_${name}_alignment.tmp2.bam >\${rna_set}_${name}_alignment.bam
                rm \${rna_set}_${name}_alignment.tmp1.bam \${rna_set}_${name}_alignment.tmp2.bam
            fi
            echo "INFO: samtools sort -o \${rna_set}@${name}_sorted.bam \${rna_set}_${name}_alignment.bam"
            samtools sort -o \${rna_set}@${name}_sorted.bam \${rna_set}_${name}_alignment.bam 
            echo "INFO: samtools index \${rna_set}@${name}_sorted.bam"
            samtools index \${rna_set}@${name}_sorted.bam
            
            if [ "${remove_duplicates}" == "yes" ]; then
                ## check read header whether they have UMI tags which are separated with underscore.(eg. NS5HGY:2:11_GTATAACCTT)
                umiCheck=\$(samtools view \${rna_set}@${name}_sorted.bam |head -n 1 | awk 'BEGIN {FS="\\t"}; {print \$1}' | awk 'BEGIN {FS=":"}; \$NF ~ /_/ {print \$NF}')
                
                # based on remove_duplicates_based_on_UMI_after_mapping
                if [ "${remove_duplicates_based_on_UMI_after_mapping}" == "yes" -a ! -z "\$umiCheck" ]; then
                    echo "INFO: umi_mark_duplicates.py will be executed for removing duplicates from bam file"
                    echo "python umi_mark_duplicates.py -f \${rna_set}@${name}_sorted.bam -p 4"
                    python umi_mark_duplicates.py -f \${rna_set}@${name}_sorted.bam -p 4
                else
                    echo "INFO: Picard MarkDuplicates will be executed for removing duplicates from bam file"
                    if [ "${remove_duplicates_based_on_UMI_after_mapping}" == "yes"  ]; then
                        echo "WARNING: Read header have no UMI tags which are separated with underscore. Picard MarkDuplicates will be executed to remove duplicates from alignment file (bam) instead of remove_duplicates_based_on_UMI_after_mapping."
                    fi
                    echo "INFO: picard MarkDuplicates OUTPUT=\${rna_set}@${name}_sorted.deumi.sorted.bam METRICS_FILE=${name}_picard_PCR_duplicates.log  VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=false INPUT=\${rna_set}@${name}_sorted.bam"
                    picard MarkDuplicates OUTPUT=\${rna_set}@${name}_sorted.deumi.sorted.bam METRICS_FILE=${name}_picard_PCR_duplicates.log  VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=false INPUT=\${rna_set}@${name}_sorted.bam 
                fi
                #get duplicates stats (read the sam flags)
                samtools flagstat \${rna_set}@${name}_sorted.deumi.sorted.bam > \${k2}@\${rna_set}@${name}_duplicates_stats.log
                #remove alignments marked as duplicates
                samtools view -b -F 0x400 \${rna_set}@${name}_sorted.deumi.sorted.bam > \${rna_set}@${name}_sorted.deumi.sorted.bam.x_dup
                #sort deduplicated files by chrom pos
                echo "INFO: samtools sort -o \${rna_set}@${name}_sorted.dedup.bam \${rna_set}@${name}_sorted.deumi.sorted.bam.x_dup"
                samtools sort -o \${rna_set}@${name}_sorted.dedup.bam \${rna_set}@${name}_sorted.deumi.sorted.bam.x_dup 
                samtools index \${rna_set}@${name}_sorted.dedup.bam
                #get flagstat after dedup
                echo "##After Deduplication##" >> \${k2}@\${rna_set}@${name}_duplicates_stats.log
                samtools flagstat \${rna_set}@${name}_sorted.dedup.bam >> \${k2}@\${rna_set}@${name}_duplicates_stats.log
            fi
            
        
            for file in unmapped/*; do mv \$file \${file/.unmapped/}; done ##remove .unmapped from filename
            if [ "\${alignersListAr[\$k-1]}" == "bowtie2" ]; then
                grep -v Warning \${k2}_${name}.bow_\${rna_set} > ${name}.tmp
                mv ${name}.tmp \${k2}_${name}.bow_\${rna_set}
                cp \${k2}_${name}.bow_\${rna_set} ./../bowfiles/.
            elif [ "\${alignersListAr[\$k-1]}" == "bowtie" ]; then
                cp \${k2}_${name}.bow1_\${rna_set} ./../bowfiles/.
            elif [ "\${alignersListAr[\$k-1]}" == "STAR" ]; then
                cp \${k2}_${name}.star_\${rna_set} ./../bowfiles/.
            fi
            cd ..
            # if filter is on, remove previously created unmapped fastq. 
            if [ "\${filtersListAr[\$k-1]}" == "Yes" ]; then
                if [ "\${prev}" != "reads" ]; then
                    echo "INFO: remove prev: \${prev}/*"
                    rm -rf \${prev}/*
                elif  [ "${remove_previous_reads}" == "true" ]; then
                    echo "INFO: inputs reads will be removed if they are located in the workdir"
                    for f in \${prev}/*; do
                        targetFile=\$(readlink -e \$f)
                        echo "INFO: targetFile: \$targetFile"
                        if [[ \$targetFile == *"\${workflowWorkDir}"* ]]; then
                            rm -f \$targetFile
                            echo "INFO: \$targetFile located in workdir and deleted."
                        fi
                    done
                fi
            # if filter is off remove current unmapped fastq
            else
                echo "INFO: remove \${rna_set}/unmapped/*"
                rm -rf \${rna_set}/unmapped/*
            fi
        else
            echo "WARNING: \${indexesListAr[\$k-1]} Mapping skipped. File not found."
            cd unmapped 
            ln -s \${wrkDir}/\${rna_set}/*fastq .
            cd ..
            cd ..
        fi
        
        if [ "\${filtersListAr[\$k-1]}" == "Yes" ]; then
            prev=\${rna_set}/unmapped
        fi
    done
    cd final_reads && ln -s \${wrkDir}/\${prev}/* .
else 
    mv ${reads} final_reads/.
fi
"""

}
}


params.star_index =  ""  //* @input

process STAR_align_Single_Best_Multimapper {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /${newName}.(bam|bam.bai)$/) "star_main_alignment/$filename"
	else if (filename =~ /${newName}Log.final.out$/) "star_main_alignment/$filename"
}

input:
 set val(name), file(reads) from g42_32_reads_g_4
 val mate from g_1_mate_g_4

output:
 set val(name), file("${newName}.{bam,bam.bai}")  into g_4_mapped_reads_g_5
 set val(name), file("${newName}_transcriptome.bam")  into g_4_transcriptome_bam
 set val(name), file("${newName}Log.final.out")  into g_4_logOut_g_21

errorStrategy 'retry'

when:
(params.run_STAR && (params.run_STAR == "yes")) || !params.run_STAR

script:
    
params_STAR = params.STAR_align_Single_Best_Multimapper.params_STAR
nameAll = reads.toString()
nameArray = nameAll.split(' ')

if (nameAll.contains('.gz')) {
    newName =  nameArray[0] - ~/(\.fastq.gz)?(\.fq.gz)?$/
    file =  nameAll - '.gz' - '.gz'
    runGzip = "ls *.gz | xargs -i echo gzip -df {} | sh"
} else {
    newName =  nameArray[0] - ~/(\.fastq)?(\.fq)?$/
    file =  nameAll 
    runGzip = ''
}

"""
$runGzip
STAR ${params_STAR}  \
    --genomeDir ${params.star_index} \
    --readFilesIn $file \
    --outSAMtype BAM SortedByCoordinate \
	--outFileNamePrefix ${newName} \
	--quantMode TranscriptomeSAM \
	--outReadsUnmapped Fastx

echo "Alignment completed."

if [ -e "${newName}Aligned.sortedByCoord.out.bam" ] ; then
    mv ${newName}Aligned.sortedByCoord.out.bam ${newName}.bam
fi
samtools index ${newName}.bam

if [ -e "${newName}Aligned.toTranscriptome.out.bam" ] ; then
    mv ${newName}Aligned.toTranscriptome.out.bam ${newName}_transcriptome.bam
fi


"""


}


process STAR_Summary {

input:
 set val(name), file(alignSum) from g_4_logOut_g_21.groupTuple()

output:
 file "*.tsv"  into g_21_outputFile_g_22
 val "star_alignment_sum"  into g_21_name_g_22

shell:
'''
#!/usr/bin/env perl
use List::Util qw[min max];
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage; 
use Data::Dumper;

my %tsv;
my @headers = ();
my $name = "!{name}";

alteredAligned();

my @keys = keys %tsv;
my $summary = "$name"."_star_sum.tsv";
my $header_string = join("\\t", @headers);
`echo "$header_string" > $summary`;
foreach my $key (@keys){
	my $values = join("\\t", @{ $tsv{$key} });
	`echo "$values" >> $summary`;
}


sub alteredAligned
{
	my @files = qw(!{alignSum});
	my $multimappedSum;
	my $alignedSum;
	my $inputCountSum;
	push(@headers, "Sample");
    push(@headers, "Total Reads");
	push(@headers, "Multimapped Reads Aligned (STAR)");
	push(@headers, "Unique Reads Aligned (STAR)");
	foreach my $file (@files){
		my $multimapped;
		my $aligned;
		my $inputCount;
		chomp($inputCount = `cat $file | grep 'Number of input reads' | awk '{sum+=\\$6} END {print sum}'`);
		chomp($aligned = `cat $file | grep 'Uniquely mapped reads number' | awk '{sum+=\\$6} END {print sum}'`);
		chomp($multimapped = `cat $file | grep 'Number of reads mapped to multiple loci' | awk '{sum+=\\$9} END {print sum}'`);
		$multimappedSum += int($multimapped);
        $alignedSum += int($aligned);
        $inputCountSum += int($inputCount);
	}
	$tsv{$name} = [$name, $inputCountSum];
	push(@{$tsv{$name}}, $multimappedSum);
	push(@{$tsv{$name}}, $alignedSum);
}
'''

}


process Merge_TSV_Files_star {

input:
 file tsv from g_21_outputFile_g_22.collect()
 val outputFileName from g_21_name_g_22.collect()

output:
 file "${name}.tsv"  into g_22_outputFileTSV_g_26

script:
name = outputFileName[0]
"""    
awk 'FNR==1 && NR!=1 {  getline; } 1 {print} ' *.tsv > ${name}.tsv
"""
}

params.single_transcript_gtf =  ""  //* @input
params.genome2bit =  ""  //* @input
params.makeavggene_main_path =  ""  //* @input
params.single_transcript_UTR_csv =  ""  //* @input
params.single_transcript_stopcodons_csv =  ""  //* @input
params.densebuilder_main_path =  ""  //* @input
params.riboseq_buildDenseTables_rpkm_path =  ""  //* @input
params.riboseq_buildDenseTables_rpkm_utr3adj_path =  ""  //* @input

thread = "10" 

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 2
    $MEMORY = 50
}
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 1300
    $CPU  = 2
    $MEMORY = 50
    $QUEUE = "long"
} 
//* platform
//* autofill

process riboseq_Densebuilder {

input:
 set val(name), file(bambai) from g_4_mapped_reads_g_5

output:
 set val(name), file("${dirName}")  into g_5_outputDir_g_9

script:
bambai = bambai.toString()
bam = bambai.split(" ")[0]
dirName = name + "Densebuilder"
"""
#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg') ### set backend
import matplotlib.pyplot as plt
plt.rcParams['pdf.fonttype'] = 42 # this keeps most text as actual text in PDFs, not outlines

import sys
import os
import argparse
import subprocess
import collections
from collections import OrderedDict
import pysam
import glob
import pandas as pd
import seaborn as sns
import numpy as np
from pathos.multiprocessing import ProcessingPool as Pool
import multiprocessing
import importlib

ftsize = [15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40]
#ftsize = [35]
### Footprint assignment
norm_type = "rpm" # either 'raw' for no normalization or 'rpm' for normalization to all mapped reads in alignment BAM file
threshold = '0'

flmin = 28
flmax = 35 # set to exact readsize that is desired
eAmin = 21
eAmax = 24
eEmin = 18
eEmax = 19

### Count tables:
cdsDenThresh = True
raw_dense_thresh = 0.01
rpm_dense_thresh = 0.001

# select insets to use to avoid start and stop peaks
inset_choice = 'default' # defaultInsets, zeroInsets, or customInsets
UTRfilestring = "${params.single_transcript_UTR_csv}"

def densebuilder_multi_raw(fiveorthreeprime, ribosome_site="0"):

    assignment = fiveorthreeprime
    totreads = 1E6
    ribosome_shift = "0"
    bampath = "${name}/${name}_star_default/"
    
    def densebuilder_one(readsize):
        readsize = str(readsize)
        riboshiftdict = "{"+readsize+':[0]'+"}"
        densitypath = "${name}/DensityUnnormalized/density"+assignment+"p_"+ribosome_shift+"shift_"+readsize+"/${name}_"+readsize+"f/"
        densityoutput = densitypath + "${name}_"+readsize+"f_"
        bamfileoutput = bampath + "${name}_"+readsize+"_match.bam"
        bamfilesorted = bampath + "${name}_"+readsize+"_match.sorted.bam"
        try:
            if not os.path.exists(densitypath): os.makedirs(densitypath)
        except OSError as err:
            print(err)
        try:
            if not os.path.exists(bampath): os.makedirs(bampath)
        except OSError as err:
            print(err)
        
        cmd = ("python ${params.densebuilder_main_path} --bamfileinput " + "${bam}" + 
            " --GTFfile " + "${params.single_transcript_gtf}" + 
            " --twobitfile " + "${params.genome2bit}" + 
            " --assignment " + assignment + 
            " --riboshiftdict " +  riboshiftdict +
            " --threshold " + str(threshold) + 
            " --totreads " + str(totreads) + 
            " --outputdata " + densityoutput +
            " --bamfileoutput " + bamfileoutput)
        print cmd
        os.system(cmd)
        ### sorting and indexing output bamfile
        sort_bamfile = 'samtools sort -o ' + bamfilesorted + " " + bamfileoutput 
        os.system(sort_bamfile)
        index_bamfile = 'samtools index %s' % (bamfilesorted)
        os.system(index_bamfile)
        clean_bamfile = 'rm %s' % (bamfileoutput)
        os.system(clean_bamfile)
    p = Pool(nodes=int(${thread}))
    p.map(densebuilder_one, ftsize)
    
    # now calculate the total number of reads from the bamfiles in ftsize
    # use pysam to calculate the number of reads present for all read lengths
    readcounter = 0
    for readlen in ftsize:
        bamfilesorted = bampath + "${name}_"+str(readlen)+"_match.sorted.bam"
        bamfile = pysam.AlignmentFile(bamfilesorted, 'rb')
        readnumber = bamfile.count()
        readcounter += readnumber
    
    readcount_outfile = "${name}/${name}_FPassigned_counts.txt"
    readcountout = open(readcount_outfile, "w")
    readcountout.write(str(readcounter))
    readcountout.close()

def densebuilder_multi_norm(self, ribosome_site="0"):
    assignment = "5" # 5' or 3' alignment
    ribosome_shift = "0"
    bampath = "${name}/${name}_star_default/"
    try:
        if not os.path.exists(bampath): os.makedirs(bampath)
    except OSError as err:
        print(err)
    
    totreads_countfile = "${name}/${name}_FPassigned_counts.txt"
    totreadcountf = open(totreads_countfile, "r")
    totreads = int(totreadcountf.read())
    totreadcountf.close()
    
    def densebuilder_one(readsize):
        readsize = str(readsize)
        riboshiftdict = "{"+readsize+':[0]'+"}"
        densitypath = "${name}/Density_rpm/density"+assignment+"p_"+ribosome_shift+"shift_"+readsize+"/${name}_"+readsize+"f/"
        densityoutput = densitypath + "${name}_"+readsize+"f_"
        bamfileoutput = bampath + "${name}_"+readsize+"_match.bam"
        bamfilesorted = bampath + "${name}_"+readsize+"_match.sorted.bam"
        try:
            if not os.path.exists(densitypath): os.makedirs(densitypath)
        except OSError as err:
            print(err)

        
        cmd = ("python ${params.densebuilder_main_path} --bamfileinput " + "${bam}" + 
            " --GTFfile " + "${params.single_transcript_gtf}" + 
            " --twobitfile " + "${params.genome2bit}" + 
            " --assignment " + assignment + 
            " --riboshiftdict " +  riboshiftdict +
            " --threshold " + str(threshold) + 
            " --totreads " + str(totreads) + 
            " --outputdata " + densityoutput +
            " --bamfileoutput " + bamfileoutput)
        print cmd
        os.system(cmd)
        ### sorting and indexing output bamfile
        sort_bamfile = 'samtools sort -o ' + bamfilesorted + " " + bamfileoutput 
        os.system(sort_bamfile)
        index_bamfile = 'samtools index %s' % (bamfilesorted)
        os.system(index_bamfile)
        clean_bamfile = 'rm %s' % (bamfileoutput)
        os.system(clean_bamfile)
    p = Pool(nodes=int(${thread}))
    p.map(densebuilder_one, ftsize)

def avggene_multi(alignposition, fiveorthreeprime):
    #Create avggenes for all readlengths defined by ftsize
	#Run all avggene scripts in parallel using pathos
	
	alignpos = alignposition # '1' for start, '2' for stop
	assignment = fiveorthreeprime # '5' for fiveprime mapping, '3' for threeprime mapping
	norm = "uneq" # 'uneq' for no normalizaiton; 'eq' to give equal weight to all genes
	ribosome_shift = "0" 
	read_length_list = ftsize
	
	if norm == 'uneq':
		equalweight = '0'
	if norm == 'eq':
		equalweight = '1'
		
	filtermodule = str(0)
	exclusionmodule = str(0)
	thresh = str(0)  # rpkm threshold
	if alignpos == '1':
		regionlength5 = '50'
		regionlength3 = '150'
		cdslength = regionlength3
	if alignpos == '2':
		regionlength5 = '150'
		regionlength3 = '100'
		cdslength = regionlength5
		
	if norm_type == "raw":
		densitystring = "DensityUnnormalized"
	elif norm_type == "rpm":
		densitystring = "Density_rpm"
	
	outfolder = "${name}/avggene"+alignpos+"_ORF"+norm_type+"_"+ribosome_shift+"shift_"+assignment+norm+cdslength
	if not os.path.exists(outfolder):   os.makedirs(outfolder)
	
	def avggene_onesize(readin):
	    # Calculate average gene for one read size
		readsize = str(readin) # set read size
		trspdictfilestring = '%s/%s/density%sp_%sshift_%s/%s_%sf/%s_%sf_' %("${name}", densitystring, assignment, ribosome_shift, readsize, "${name}", readsize, "${name}", readsize)
		outfilebase = "%s/%s_%sf" % (outfolder, "${name}", readsize)
		commandstring = 'python ${params.makeavggene_main_path} --regionlength5 %s --regionlength3 %s --trspdictfilestring %s --UTRfilestring %s --filtermodule %s --exclusionmodule %s --threshold %s --alignpos %s --equalweight %s --outfilebase %s' %(regionlength5,regionlength3,trspdictfilestring,UTRfilestring,filtermodule,exclusionmodule,thresh,alignpos,equalweight,outfilebase)
		print commandstring
		os.system(commandstring)

	p = Pool(nodes=int(${thread}))
	p.map(avggene_onesize, read_length_list)

def get_riboshift(shiftsite, alignposition):
    alignpos = alignposition
    assignment = "5"  # must use 5' start mapping for this
    norm = "uneq"
    ribosome_shift = "0"
    shiftpos = shiftsite
    if shiftpos == "P":
        shift_adjust = 0
    if shiftpos == "A":
        shift_adjust = 3 # A site is shifted 3 nt's farther than P site from 5' end
    if shiftpos == "E":
        shift_adjust = -3
    cdslength = "150"
    avggene_paths = "${name}/avggene"+alignpos+"_ORF"+norm_type+"_"+ribosome_shift+"shift_"+assignment+norm+cdslength
    frame = pd.DataFrame()
    list_ = []
    for ft in ftsize:
		infile = "%s/%s_%sf_avg_1.csv" % (avggene_paths, "${name}", ft)
		print "infile:" + infile
		df = pd.read_csv(infile, index_col="position", header=0)
		df.columns = ["%d" % (ft)]
		list_.append(df)
    frame = pd.concat(list_, axis=1)
	# Get the max (or modal read number per position) of footprints relative to start codon
	# Adding a +1 to account of +/- 1nt from digestion heterogeneity
	# Only look between 11 to 16 nts upstream of start for 5' aligned reads
	# Define full length reads from 31-34 nt's in length
	# retrieve these using 31-15 == 16, and 34-15 == 19+1, thus 16:20
    offset_fl = 50+1+shift_adjust-frame.iloc[34:39,flmin-15:flmax+1-15].idxmax()
    print "offset_fl"
    print offset_fl
	# Define emptyA reads from 21-24 nt's in length
	# retrieve these using 21-15 == 6, and 24-15 == 9+1, thus 6:10
    offset_emptyA = 50+1+shift_adjust-frame.iloc[34:39,eAmin-15:eAmax+1-15].idxmax()
    print "offset_emptyA"
    print offset_emptyA
	# Define emptyE reads from 18-19 nt's in length - need to cut close to avoid empty A
	# retrieve these using 18-15 == 3, and 19-15 == 4, thus 3:5
	# also need to only look around the start codon for these reads since the 5' end is directly over the 'A' of the AUG
    offset_emptyE = 50+1+shift_adjust-frame.iloc[49:52,eEmin-15:eEmax+1-15].idxmax()
    print "offset_emptyE"
    print offset_emptyE
    outpath = "${name}/riboshift"
    if not os.path.exists(outpath):     os.makedirs(outpath)
    outfile_fl = "%s/full_length_%soffsets.csv" % (outpath, shiftpos)
    outfile_emptyA = "%s/emptyA_%soffsets.csv" % (outpath, shiftpos)
    outfile_emptyE = "%s/emptyE_%soffsets.csv" % (outpath, shiftpos)
    offset_fl.to_csv(outfile_fl, index = True, header = False)
    offset_emptyA.to_csv(outfile_emptyA, index = True, header = False)
    offset_emptyE.to_csv(outfile_emptyE, index = True, header = False)

def densebuilder_riboshift(ribosome_site):
    #Function to map densities for ranges of read lengths to single sites of the ribosome
    #Normally map to the A site of the ribosome
    assignment = "5" # should be five prime assignment to use shift dicts
    ribosome_shift = ribosome_site # 'A' or 'P'
    bamfileinput = "${bam}"
    bamfileoutpath = "${name}/${name}_star_default/"
    if not os.path.exists(bamfileoutpath): os.makedirs(bamfileoutpath)
    densitystring = "Density_rpm"
    totreads_countfile = "${name}/${name}_FPassigned_counts.txt"
    totreadcountf = open(totreads_countfile, "r")
    totreads = int(totreadcountf.read())
    totreadcountf.close()
    ## Set paths to riboshift dicts for each species of interest
    fl_infile = "${name}/riboshift/full_length_%soffsets.csv" % (ribosome_shift)
    eA_infile = "${name}/riboshift/emptyA_%soffsets.csv" % (ribosome_shift)
    eE_infile = "${name}/riboshift/emptyE_%soffsets.csv" % (ribosome_shift)
    fl_ribo_shift = pd.read_csv(fl_infile, index_col = 0, header=None)
    fl_ribo_shift.columns = ['riboshift']
    eA_ribo_shift = pd.read_csv(eA_infile, index_col = 0, header=None)
    eA_ribo_shift.columns = ['riboshift']
    eE_ribo_shift = pd.read_csv(eE_infile, index_col = 0, header=None)
    eE_ribo_shift.columns = ['riboshift']
    fl_sizes = fl_ribo_shift.index.tolist()
    eA_sizes = eA_ribo_shift.index.tolist()
    eE_sizes = eE_ribo_shift.index.tolist()
    fl_shiftdict = collections.OrderedDict()
    eA_shiftdict = collections.OrderedDict()
    eE_shiftdict = collections.OrderedDict()
    for i in fl_sizes:
        fl_shiftdict.update( {str(i):'['+str(fl_ribo_shift.loc[i,'riboshift'])+']'})
    for i in eA_sizes:
        eA_shiftdict.update( {str(i):'['+str(eA_ribo_shift.loc[i,'riboshift'])+']'})
    for i in eE_sizes:
        eE_shiftdict.update( {str(i):'['+str(eE_ribo_shift.loc[i,'riboshift'])+']'})
    allLenghts_shiftdict = collections.OrderedDict()
    allLenghts_shiftdict.update(eE_shiftdict)
    allLenghts_shiftdict.update(eA_shiftdict)
    allLenghts_shiftdict.update(fl_shiftdict)
    fp_assign_path = "${name}"
    
    def build_fl_dense(fldict):
	    readsize_dict = fldict
	    pop = "fl"
	    minlen = str(readsize_dict.keys()[0])
	    maxlen = str(readsize_dict.keys()[-1])
	    dict_conv = dict(readsize_dict)
	    # fl_shift:
	    riboshiftdict = str(dict_conv).replace(" ","")
	    # riboshiftdict = {'33':'17','32':'17','31':'16','34':'17'} ## manually set readlengths
	    bamfileoutput = '%s/%s_%s_match.bam' % (bamfileoutpath, "${name}", pop)
	    densityfilepath = '%s/%s/density%sp_%sshift_%sto%s/%s_%sto%sf' %("${name}",densitystring, assignment, ribosome_shift, minlen, maxlen, "${name}", minlen,maxlen)
	    densityfileout = '%s/%s_%sto%sf_' % (densityfilepath, "${name}", minlen, maxlen)
	    if not os.path.exists(densityfilepath): os.makedirs(densityfilepath)
	    commandstring = 'python ${params.densebuilder_main_path} --bamfileinput %s --GTFfile %s --twobitfile %s --assignment %s --riboshiftdict %s --threshold %s --totreads %s --outputdata %s --bamfileoutput %s' % (bamfileinput, "${params.single_transcript_gtf}", "${params.genome2bit}", assignment,riboshiftdict, threshold, totreads, densityfileout, bamfileoutput)
	    print commandstring
	    os.system(commandstring)
	    ## sort and index bamfile
	    sort_bamfile = 'samtools sort -o %s/%s_%s_match.sorted.bam %s' % (bamfileoutpath, "${name}", pop, bamfileoutput)
	    os.system(sort_bamfile)
	    index_bamfile = 'samtools index %s/%s_%s_match.sorted.bam' % (bamfileoutpath,"${name}", pop)
	    os.system(index_bamfile)
	    clean_bamfile = 'rm %s' % (bamfileoutput)
	    os.system(clean_bamfile)
    def build_eA_dense(eAdict):
        readsize_dict = eAdict
        pop = "eA"
        minlen = str(readsize_dict.keys()[0])
        maxlen = str(readsize_dict.keys()[-1])
        dict_conv = dict(readsize_dict)
        # eA_shift:
        riboshiftdict = str(dict_conv).replace(" ","")
        # riboshiftdict = {'33':'17','32':'17','31':'16','34':'17'}
        bamfileoutpath = '%s/%s_star_default' % ("${name}", "${name}")
        bamfileoutput = '%s/%s_%s_match.bam' % (bamfileoutpath, "${name}", pop)
        densityfilepath = '%s/%s/density%sp_%sshift_%sto%s/%s_%sto%sf' %("${name}",densitystring, assignment, ribosome_shift, minlen, maxlen, "${name}", minlen,maxlen)
        densityfileout = '%s/%s_%sto%sf_' % (densityfilepath, "${name}", minlen, maxlen)
        if not os.path.exists(densityfilepath): os.makedirs(densityfilepath)
        commandstring = 'python ${params.densebuilder_main_path} --bamfileinput %s --GTFfile %s --twobitfile %s --assignment %s --riboshiftdict %s --threshold %s --totreads %s --outputdata %s --bamfileoutput %s' % (bamfileinput, "${params.single_transcript_gtf}", "${params.genome2bit}", assignment,riboshiftdict, threshold, totreads, densityfileout, bamfileoutput)
        print commandstring
        os.system(commandstring)
        ## sort and index bamfile
        sort_bamfile = 'samtools sort -o %s/%s_%s_match.sorted.bam %s' % (bamfileoutpath, "${name}", pop, bamfileoutput)
        os.system(sort_bamfile)
        index_bamfile = 'samtools index %s/%s_%s_match.sorted.bam' % (bamfileoutpath,"${name}", pop)
        os.system(index_bamfile)
        clean_bamfile = 'rm %s' % (bamfileoutput)
        os.system(clean_bamfile)
    def build_eE_dense(eEdict):
        readsize_dict = eEdict
        pop = "eE"
        minlen = str(readsize_dict.keys()[0])
        maxlen = str(readsize_dict.keys()[-1])
        dict_conv = dict(readsize_dict)
        # eE_shift:
        riboshiftdict = str(dict_conv).replace(" ","")
        # riboshiftdict = {'33':'17','32':'17','31':'16','34':'17'}
        bamfileoutpath = '%s/%s_star_default' % (fp_assign_path, "${name}")
        bamfileoutput = '%s/%s_%s_match.bam' % (bamfileoutpath, "${name}", pop)
        densityfilepath = '%s/%s/density%sp_%sshift_%sto%s/%s_%sto%sf' %(fp_assign_path,densitystring, assignment, ribosome_shift, minlen, maxlen, "${name}", minlen,maxlen)
        densityfileout = '%s/%s_%sto%sf_' % (densityfilepath, "${name}", minlen, maxlen)
        if not os.path.exists(densityfilepath): os.makedirs(densityfilepath)
        commandstring = 'python ${params.densebuilder_main_path} --bamfileinput %s --GTFfile %s --twobitfile %s --assignment %s --riboshiftdict %s --threshold %s --totreads %s --outputdata %s --bamfileoutput %s' % (bamfileinput, "${params.single_transcript_gtf}", "${params.genome2bit}", assignment,riboshiftdict, threshold, totreads, densityfileout, bamfileoutput)
        print commandstring
        os.system(commandstring)
        ## sort and index bamfile
        sort_bamfile = 'samtools sort -o %s/%s_%s_match.sorted.bam %s' % (bamfileoutpath, "${name}", pop, bamfileoutput)
        os.system(sort_bamfile)
        index_bamfile = 'samtools index %s/%s_%s_match.sorted.bam' % (bamfileoutpath,"${name}", pop)
        os.system(index_bamfile)
        clean_bamfile = 'rm %s' % (bamfileoutput)
        os.system(clean_bamfile)
    def build_allLengths_dense(alldict):
        readsize_dict = alldict
        pop = "aL"
        minlen = str(readsize_dict.keys()[0])
        maxlen = str(readsize_dict.keys()[-1])
        dict_conv = dict(readsize_dict)
        # aL_shift:
        riboshiftdict = str(dict_conv).replace(" ","")
        # riboshiftdict = {'33':'17','32':'17','31':'16','34':'17'}
        bamfileoutpath = '%s/%s_star_default' % (fp_assign_path, "${name}")
        bamfileoutput = '%s/%s_%s_match.bam' % (bamfileoutpath, "${name}", pop)
        densityfilepath = '%s/%s/density%sp_%sshift_%s/%s_%sf' %(fp_assign_path,densitystring, assignment, ribosome_shift, pop, "${name}", pop)
        densityfileout = '%s/%s_%sf_' % (densityfilepath, "${name}", pop)
        if not os.path.exists(densityfilepath): os.makedirs(densityfilepath)
        commandstring = 'python ${params.densebuilder_main_path} --bamfileinput %s --GTFfile %s --twobitfile %s --assignment %s --riboshiftdict %s --threshold %s --totreads %s --outputdata %s --bamfileoutput %s' % (bamfileinput, "${params.single_transcript_gtf}", "${params.genome2bit}", assignment,riboshiftdict, threshold, totreads, densityfileout, bamfileoutput)
        print commandstring
        os.system(commandstring)
        ## sort and index bamfile
        sort_bamfile = 'samtools sort -o %s/%s_%s_match.sorted.bam %s' % (bamfileoutpath, "${name}", pop, bamfileoutput)
        os.system(sort_bamfile)
        index_bamfile = 'samtools index %s/%s_%s_match.sorted.bam' % (bamfileoutpath,"${name}", pop)
        os.system(index_bamfile)
        clean_bamfile = 'rm %s' % (bamfileoutput)
        os.system(clean_bamfile)
    p1 = multiprocessing.Process(target=build_fl_dense, args=(fl_shiftdict,))
    p1.start()
    p2 = multiprocessing.Process(target=build_eA_dense, args=(eA_shiftdict,))
    p2.start()
    p3 = multiprocessing.Process(target=build_eE_dense, args=(eE_shiftdict,))
    p3.start()
    p4 = multiprocessing.Process(target=build_allLengths_dense, args=(allLenghts_shiftdict,))
    p4.start()
    p1.join()
    p2.join()
    p3.join()
    p4.join()

def avggene_riboshift(alignposition, ribosome_site, normalization, thresh):
    #Create avggenes for all readlengths defined by ftsize
    #Run all avggene scripts in parallel using pathos
    alignpos = alignposition # '1' for start, '2' for stop
    assignment = "5" # should be 5' mapped at this point
    norm = normalization # 'uneq' for no normalizaiton; 'eq' to give equal weight to allgenes
    ribosome_shift = ribosome_site # 'A' or 'P' or maybe 'E' later
    
    if norm == 'uneq':
        equalweight = '0'
    if norm == 'eq':
        equalweight = '1' 	#equalweight = '0'    # 1 gives equal weight to all genes.
    filtermodule = str(0)
    exclusionmodule = str(0)
    if alignpos == '1':
        regionlength5 = '50'
        regionlength3 = '150'
        cdslength = regionlength3
    if alignpos == '2':
        regionlength5 = '150'
        regionlength3 = '100'
        cdslength = regionlength5
    if norm_type == "raw":
        densitystring = "DensityUnnormalized"
    elif norm_type == "rpm":
        densitystring = "Density_rpm"
    else:
        print "Normalization is not set for densebuilder!"
        sys.exit()
    fp_assign_path = "${name}"
    
    outfolder =  "%s/avggene%s_ORF%s_%sshift_%s%s%s" % (fp_assign_path, alignpos, norm_type, ribosome_shift, assignment, norm, cdslength) # norm should be 'uneq' fornow
    if not os.path.exists(outfolder):   os.makedirs(outfolder)
    ### define read sizes for fl, eA, and eE - need minlen and maxlen
    fl_infile = "%s/riboshift/full_length_%soffsets.csv" % (fp_assign_path,ribosome_shift)
    eA_infile = "%s/riboshift/emptyA_%soffsets.csv" % (fp_assign_path, ribosome_shift)
    eE_infile = "%s/riboshift/emptyE_%soffsets.csv" % (fp_assign_path, ribosome_shift)
    fl_ribo_shift = pd.read_csv(fl_infile, index_col = 0, header=None)
    fl_ribo_shift.columns = ['riboshift']
    eA_ribo_shift = pd.read_csv(eA_infile, index_col = 0, header=None)
    eA_ribo_shift.columns = ['riboshift']
    eE_ribo_shift = pd.read_csv(eE_infile, index_col = 0, header=None)
    eE_ribo_shift.columns = ['riboshift']
    # print eE_ribo_shift.index.tolist()
    fl_sizes = fl_ribo_shift.index.tolist()
    eA_sizes = eA_ribo_shift.index.tolist()
    eE_sizes = eE_ribo_shift.index.tolist()
    minlen_fl = fl_sizes[0]
    maxlen_fl = fl_sizes[-1]
    minlen_eA = eA_sizes[0]
    maxlen_eA = eA_sizes[-1]
    minlen_eE = eE_sizes[0]
    maxlen_eE = eE_sizes[-1]
    def avggene_fl():
        # Calculate average gene for one read size
        # readsize = str(readin)
        # change this to read in fl reads
        pop = "fl"
        minlen = minlen_fl
        maxlen = maxlen_fl
        trspdictfilestring = '%s/%s/density%sp_%sshift_%sto%s/%s_%sto%sf/%s_%sto%sf_' % (fp_assign_path, densitystring, assignment, ribosome_shift, minlen, maxlen, "${name}", minlen, maxlen, "${name}", minlen, maxlen)
        outfilebase = "%s/%s_%s_rpkmThresh%s_%sto%sf" % (outfolder, "${name}", pop,thresh, minlen, maxlen)
        print ""
        commandstring = 'python ${params.makeavggene_main_path} --regionlength5 %s --regionlength3 %s --trspdictfilestring %s --UTRfilestring %s --filtermodule %s --exclusionmodule %s --threshold %s --alignpos %s --equalweight %s --outfilebase %s' %(regionlength5,regionlength3,trspdictfilestring,UTRfilestring,filtermodule,exclusionmodule,thresh,alignpos,equalweight,outfilebase)
        print commandstring
        os.system(commandstring)
    
    def avggene_eA():
        # Calculate average gene for one read size
        # readsize = str(readin)
        # change this to read in fl reads
        pop = "eA"
        minlen = minlen_eA
        maxlen = maxlen_eA
        trspdictfilestring = '%s/%s/density%sp_%sshift_%sto%s/%s_%sto%sf/%s_%sto%sf_' % (fp_assign_path, densitystring, assignment, ribosome_shift, minlen, maxlen, "${name}", minlen, maxlen, "${name}", minlen, maxlen)
        outfilebase = "%s/%s_%s_rpkmThresh%s_%sto%sf" % (outfolder, "${name}", pop,thresh, minlen, maxlen)
        print ""
        commandstring = 'python ${params.makeavggene_main_path} --regionlength5 %s --regionlength3 %s --trspdictfilestring %s --UTRfilestring %s --filtermodule %s --exclusionmodule %s --threshold %s --alignpos %s --equalweight %s --outfilebase %s' %(regionlength5,regionlength3,trspdictfilestring,UTRfilestring,filtermodule,exclusionmodule,thresh,alignpos,equalweight,outfilebase)
        print commandstring
        os.system(commandstring)
    
    def avggene_eE():
        # Calculate average gene for one read size
        # readsize = str(readin)
        # change this to read in fl reads
        pop = "eE"
        minlen = minlen_eE
        maxlen = maxlen_eE
        trspdictfilestring = '%s/%s/density%sp_%sshift_%sto%s/%s_%sto%sf/%s_%sto%sf_' % (fp_assign_path, densitystring, assignment, ribosome_shift, minlen, maxlen, "${name}", minlen, maxlen, "${name}", minlen, maxlen)
        outfilebase = "%s/%s_%s_rpkmThresh%s_%sto%sf" % (outfolder, "${name}", pop,thresh, minlen, maxlen)
        commandstring = 'python ${params.makeavggene_main_path} --regionlength5 %s --regionlength3 %s --trspdictfilestring %s --UTRfilestring %s --filtermodule %s --exclusionmodule %s --threshold %s --alignpos %s --equalweight %s --outfilebase %s' %(regionlength5,regionlength3,trspdictfilestring,UTRfilestring,filtermodule,exclusionmodule,thresh,alignpos,equalweight,outfilebase)
        print commandstring
        os.system(commandstring)
    def avggene_aL():
        pop = "aL"
        # minlen = minlen_eE
        # maxlen = maxlen_eE
        trspdictfilestring = '%s/%s/density%sp_%sshift_%s/%s_%sf/%s_%sf_' % (fp_assign_path, densitystring, assignment, ribosome_shift, pop, "${name}", pop,"${name}", pop)
        outfilebase = "%s/%s_rpkmThresh%s_%s_%sf" % (outfolder, "${name}", thresh, pop,pop)
        commandstring = 'python ${params.makeavggene_main_path} --regionlength5 %s --regionlength3 %s --trspdictfilestring %s --UTRfilestring %s --filtermodule %s --exclusionmodule %s --threshold %s --alignpos %s --equalweight %s --outfilebase %s' %(regionlength5,regionlength3,trspdictfilestring,UTRfilestring, filtermodule,exclusionmodule,thresh,alignpos,equalweight,outfilebase)
        print commandstring
        os.system(commandstring)
    p1 = multiprocessing.Process(target=avggene_fl)
    p1.start()
    p2 = multiprocessing.Process(target=avggene_eA)
    p2.start()
    p3 = multiprocessing.Process(target=avggene_eE)
    p3.start()
    p4 = multiprocessing.Process(target=avggene_aL)
    p4.start()
    p1.join()
    p2.join()
    p3.join()
    p4.join()

def build_count_tables(alignposition="1", fiveorthreeprime="5", normalization="uneq", ribosome_site="A"):
    # wrapper to direct building of density tables
    # count tables are based on cds, utr5, and utr3 sizes defined by UTRfilesting
	# inset_choice contains adjustments to avoid start and stop codon peaks in data tables
	# output is a csv file with densities for each region and 3'UTR occupancies

    # mappingpos = "5"
	alignpos = alignposition # '1' for start, '2' for stop
	assignment = fiveorthreeprime # '5' for fiveprime mapping, '3' for threeprime mapping
	norm = normalization # 'uneq' for no normalizaiton; 'eq' to give equal weight to allgenes
	ribosome_shift = ribosome_site # 'A' or 'P' or maybe 'E' later
	flminlen = str(flmin)
	flmaxlen = str(flmax)
	eAminlen = str(eAmin)
	eAmaxlen = str(eAmax)
	eEminlen = str(eEmin)
	eEmaxlen = str(eEmax)
	fp_assign_path = "${name}"
	outfolder =  "%s/countTables" % (fp_assign_path) # make a new folder to store counttables
	if not os.path.exists(outfolder):   os.makedirs(outfolder)
	
	### require Density_rpm for densityfilestring
	densitystring = "Density_rpm"
	totreads_countfile = "%s/%s_FPassigned_counts.txt" % (fp_assign_path, "${name}")
	totreadcountf = open(totreads_countfile, "r")
	totreads = int(totreadcountf.read())
	totreadcountf.close()
	def fl_countTable():
		trspdictfilestring = '%s/%s/density%sp_%sshift_%sto%s/%s_%sto%sf/%s_%sto%sf_' %(
				fp_assign_path, densitystring, assignment, ribosome_shift, flminlen,flmaxlen, "${name}", flminlen, flmaxlen, "${name}", flminlen, flmaxlen)
		outfilestring = "%s/%s_fl_%s_%sto%s_countTable_rpkm" % (outfolder, "${name}",norm_type, flminlen, flmaxlen)
		print ""
		# need to input trspdictfilestring, UTRfilesting, cdsDenThresh,raw_dense_thresh, rpm_dense_thresh, norm_type, inset_choice
		commandstring = 'python ${params.riboseq_buildDenseTables_rpkm_path} --trspdictfilestring %s --UTRfilestring %s --cdsDenThresh %s --norm_type %s --raw_dense_thresh %s --rpm_dense_thresh %s --inset_choice %s --outfilestring %s --totreads %s' % (
			trspdictfilestring,UTRfilestring,cdsDenThresh,norm_type,raw_dense_thresh, rpm_dense_thresh,inset_choice,outfilestring, totreads)
		print commandstring
		os.system(commandstring)
	def eA_countTable():
		trspdictfilestring = '%s/%s/density%sp_%sshift_%sto%s/%s_%sto%sf/%s_%sto%sf_' %(
				fp_assign_path, densitystring, assignment, ribosome_shift, eAminlen,eAmaxlen, "${name}", eAminlen, eAmaxlen, "${name}", eAminlen, eAmaxlen)
		outfilestring = "%s/%s_eA_%s_%sto%s_countTable_rpkm" % (outfolder, "${name}",norm_type, eAminlen, eAmaxlen)
		print ""
		# need to input trspdictfilestring, UTRfilesting, cdsDenThresh,raw_dense_thresh, rpm_dense_thresh, norm_type, inset_choice
		commandstring = 'python ${params.riboseq_buildDenseTables_rpkm_path} --trspdictfilestring %s --UTRfilestring %s --cdsDenThresh %s --norm_type %s --raw_dense_thresh %s --rpm_dense_thresh %s --inset_choice %s --outfilestring %s --totreads %s' % (
			trspdictfilestring,UTRfilestring,cdsDenThresh,norm_type,raw_dense_thresh,rpm_dense_thresh,inset_choice,outfilestring, totreads)
		print commandstring
		os.system(commandstring)
	def eE_countTable():
		trspdictfilestring = '%s/%s/density%sp_%sshift_%sto%s/%s_%sto%sf/%s_%sto%sf_' %(
				fp_assign_path, densitystring, assignment, ribosome_shift, eEminlen,eEmaxlen, "${name}", eEminlen, eEmaxlen, "${name}", eEminlen, eEmaxlen)
		outfilestring = "%s/%s_eE_%s_%sto%s_countTable_rpkm" % (outfolder, "${name}",norm_type, eEminlen, eEmaxlen)
		print ""
		# need to input trspdictfilestring, UTRfilesting, cdsDenThresh,raw_dense_thresh, rpm_dense_thresh, norm_type, inset_choice
		commandstring = 'python ${params.riboseq_buildDenseTables_rpkm_path} --trspdictfilestring %s --UTRfilestring %s --cdsDenThresh %s --norm_type %s --raw_dense_thresh %s --rpm_dense_thresh %s --inset_choice %s --outfilestring %s --totreads %s' % (
			trspdictfilestring,UTRfilestring,cdsDenThresh,norm_type,raw_dense_thresh,rpm_dense_thresh,inset_choice,outfilestring, totreads)
		print commandstring
		os.system(commandstring)
	def aL_countTable():
		trspdictfilestring = '%s/%s/density%sp_%sshift_aL/%s_aLf/%s_aLf_' %(
				fp_assign_path, densitystring, assignment, ribosome_shift, "${name}","${name}",)
		outfilestring = "%s/%s_aL_%s_countTable_rpkm" % (outfolder, "${name}", norm_type)
		print ""
		# need to input trspdictfilestring, UTRfilesting, cdsDenThresh,raw_dense_thresh, rpm_dense_thresh, norm_type, inset_choice
		commandstring = 'python ${params.riboseq_buildDenseTables_rpkm_path} --trspdictfilestring %s --UTRfilestring %s --cdsDenThresh %s --norm_type %s --raw_dense_thresh %s --rpm_dense_thresh %s --inset_choice %s --outfilestring %s --totreads %s' % (
			trspdictfilestring,UTRfilestring,cdsDenThresh,norm_type,raw_dense_thresh,rpm_dense_thresh,inset_choice,outfilestring, totreads)
		print commandstring
		os.system(commandstring)
	p1 = multiprocessing.Process(target=fl_countTable)
	p1.start()
	p2 = multiprocessing.Process(target=eA_countTable)
	p2.start()
	p3 = multiprocessing.Process(target=eE_countTable)
	p3.start()
	p4 = multiprocessing.Process(target=aL_countTable)
	p4.start()
	p1.join()
	p2.join()
	p3.join()
	p4.join()

def build_count_tables_utr3adj(alignposition="1", fiveorthreeprime="5", normalization="uneq", ribosome_site="A"):
	# wrapper to direct building of density tables
	# count tables are based on cds, utr5, and utr3 sizes defined by UTRfilesting
	# inset_choice contains adjustments to avoid start and stop codon peaks in data tables
	# output is a csv file with densities for each region and 3'UTR occupancies
	# mappingpos = "5"
	stopcodons = "${params.single_transcript_stopcodons_csv}"
	alignpos = alignposition # '1' for start, '2' for stop
	assignment = fiveorthreeprime # '5' for fiveprime mapping, '3' for threeprime mapping
	norm = normalization # 'uneq' for no normalizaiton; 'eq' to give equal weight to allgenes
	ribosome_shift = ribosome_site # 'A' or 'P' or maybe 'E' later
	flminlen = str(flmin)
	flmaxlen = str(flmax)
	eAminlen = str(eAmin)
	eAmaxlen = str(eAmax)
	eEminlen = str(eEmin)
	eEmaxlen = str(eEmax)
	fp_assign_path = '${name}'
	outfolder =  "%s/countTables" % (fp_assign_path) # make a new folder to store counttables
	if not os.path.exists(outfolder):   os.makedirs(outfolder)
	
	### require Density_rpm for densityfilestring
	densitystring = "Density_rpm"
	totreads_countfile = "%s/%s_FPassigned_counts.txt" % (fp_assign_path, '${name}')
	totreadcountf = open(totreads_countfile, "r")
	totreads = int(totreadcountf.read())
	totreadcountf.close()
	def fl_countTable():
		trspdictfilestring = '%s/%s/density%sp_%sshift_%sto%s/%s_%sto%sf/%s_%sto%sf_' %(
				fp_assign_path, densitystring, assignment, ribosome_shift, flminlen,flmaxlen, '${name}', flminlen, flmaxlen, '${name}', flminlen, flmaxlen)
		outfilestring = "%s/%s_fl_%s_%sto%s_countTable_rpkm_utr3adj" % (outfolder, "${name}", norm_type, flminlen, flmaxlen)
		print ""
		commandstring = 'python ${params.riboseq_buildDenseTables_rpkm_utr3adj_path} --trspdictfilestring %s --UTRfilestring %s --cdsDenThresh %s --norm_type %s --raw_dense_thresh %s --rpm_dense_thresh %s --inset_choice %s --outfilestring %s --totreads %s --stopcodons %s' % (
			trspdictfilestring,UTRfilestring,cdsDenThresh,norm_type,raw_dense_thresh,rpm_dense_thresh,inset_choice,outfilestring, totreads, stopcodons)
		print commandstring
		os.system(commandstring)
	def eA_countTable():
		trspdictfilestring = '%s/%s/density%sp_%sshift_%sto%s/%s_%sto%sf/%s_%sto%sf_' %(fp_assign_path, densitystring, assignment, ribosome_shift, eAminlen,eAmaxlen, '${name}', eAminlen, eAmaxlen, '${name}', eAminlen, eAmaxlen)
		outfilestring = "%s/%s_eA_%s_%sto%s_countTable_rpkm_utr3adj" % (outfolder, "${name}", norm_type, eAminlen, eAmaxlen)
		print ""
		# need to input trspdictfilestring, UTRfilesting, cdsDenThresh,raw_dense_thresh, rpm_dense_thresh, norm_type, inset_choice
		commandstring = 'python ${params.riboseq_buildDenseTables_rpkm_utr3adj_path}  --trspdictfilestring %s --UTRfilestring %s --cdsDenThresh %s --norm_type %s --raw_dense_thresh %s --rpm_dense_thresh %s --inset_choice %s --outfilestring %s --totreads %s --stopcodons %s' % (trspdictfilestring,UTRfilestring,cdsDenThresh,norm_type,raw_dense_thresh,rpm_dense_thresh,inset_choice,outfilestring, totreads, stopcodons)
		print commandstring
		os.system(commandstring)
	def eE_countTable():
		trspdictfilestring = '%s/%s/density%sp_%sshift_%sto%s/%s_%sto%sf/%s_%sto%sf_' %(fp_assign_path, densitystring, assignment, ribosome_shift, eEminlen,eEmaxlen, '${name}', eEminlen, eEmaxlen, '${name}', eEminlen, eEmaxlen)
		outfilestring = "%s/%s_eE_%s_%sto%s_countTable_rpkm_utr3adj" % (outfolder, "${name}", norm_type, eEminlen, eEmaxlen)
		print ""
		# need to input trspdictfilestring, UTRfilesting, cdsDenThresh,raw_dense_thresh, rpm_dense_thresh, norm_type, inset_choice
		commandstring = 'python ${params.riboseq_buildDenseTables_rpkm_utr3adj_path} --trspdictfilestring %s --UTRfilestring %s --cdsDenThresh %s --norm_type %s --raw_dense_thresh %s --rpm_dense_thresh %s --inset_choice %s --outfilestring %s --totreads %s --stopcodons %s' % (trspdictfilestring,UTRfilestring,cdsDenThresh,norm_type,raw_dense_thresh,rpm_dense_thresh,inset_choice,outfilestring, totreads, stopcodons)
		print commandstring
		os.system(commandstring)
	def aL_countTable():
		trspdictfilestring = '%s/%s/density%sp_%sshift_aL/%s_aLf/%s_aLf_' %(fp_assign_path, densitystring, assignment, ribosome_shift, '${name}','${name}',)
		outfilestring = "%s/%s_aL_%s_countTable_rpkm_utr3adj" % (outfolder, '${name}',norm_type)
		print ""
		# need to input trspdictfilestring, UTRfilesting, cdsDenThresh,raw_dense_thresh, rpm_dense_thresh, norm_type, inset_choice
		commandstring = 'python ${params.riboseq_buildDenseTables_rpkm_utr3adj_path} --trspdictfilestring %s --UTRfilestring %s --cdsDenThresh %s --norm_type %s --raw_dense_thresh %s --rpm_dense_thresh %s --inset_choice %s --outfilestring %s --totreads %s --stopcodons %s' % (trspdictfilestring,UTRfilestring,cdsDenThresh,norm_type,raw_dense_thresh,rpm_dense_thresh,inset_choice,outfilestring, totreads, stopcodons)
		print commandstring
		os.system(commandstring)
	p1 = multiprocessing.Process(target=fl_countTable)
	p1.start()
	p2 = multiprocessing.Process(target=eA_countTable)
	p2.start()
	p3 = multiprocessing.Process(target=eE_countTable)
	p3.start()
	p4 = multiprocessing.Process(target=aL_countTable)
	p4.start()
	p1.join()
	p2.join()
	p3.join()
	p4.join()

def main():
    fiveorthreeprime = "5"
    densebuilder_multi_raw(fiveorthreeprime)
    densebuilder_multi_norm(fiveorthreeprime)
    avggene_multi("1", fiveorthreeprime)
    avggene_multi("2", fiveorthreeprime)
    get_riboshift("P", "1") # shiftsite, alignposition
    get_riboshift("A", "1")
    densebuilder_riboshift("A")
    avggene_riboshift("1","A", 'uneq', '0') 
    avggene_riboshift("2","A", 'uneq', '0')
    avggene_riboshift("1","A", 'eq', '10')
    avggene_riboshift("2","A", 'eq', '10')
    build_count_tables()
    build_count_tables_utr3adj()
    os.system("mv ${name} ${dirName}")

if __name__ == '__main__':
	# execute only if run as a script
	main()
"""
}

params.single_transcript_mRNA_csv =  ""  //* @input
params.single_transcript_UTR_csv =  ""  //* @input
params.genome2bit =  ""  //* @input
params.rphelper_path =  ""  //* @input

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 10
    $MEMORY = 30
}
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 240
    $CPU  = 10
    $MEMORY = 30
    $QUEUE = "short"
} 
//* platform
//* autofill

process riboseq_avggene_cdsNorm_start {

input:
 set val(name), file(dirName) from g_5_outputDir_g_9

output:
 set val(name), file("${newdirName}")  into g_9_outputDir_g_12

script:
newdirName = name + "riboseq_avggene_cdsNorm_start_stop"
"""
#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg') # set backend for matplotlib
import matplotlib.pyplot as plt 
plt.rcParams['pdf.fonttype'] = 42 # this keeps most text as actual text in PDFs, not outlines


import sys
sys.path.append("${params.rphelper_path}")
from Bio import SeqIO
import twobitreader
import GFF
import csv
import os
import struct
import pandas as pd
import rphelper as rph
import argparse
import importlib


threshold = 10
customSize = 30
pop = "fl"

inset_choice = 'default'
flmin = 28
flmax = 35 # set to exact readsize that is desired
eAmin = 21
eAmax = 24
eEmin = 18
eEmax = 19

assignment = "5"
ribosome_shift = "A"
norm_type = 'rpm'
densitystring = 'Density_rpm'

norm = 'eq'

if pop == "fl":
	minlen = str(flmin)
	maxlen = str(flmax)
elif pop == "eA":
	minlen = str(eAmin)
	maxlen = str(eAmax)
elif pop == "eE":
	minlen = str(eEmin)
	maxlen = str(eEmax)
elif pop == "custom":
	minlen = str(customSize)
	maxlen = str(customSize)
else:
	print "read lengths not set"



mRNAdf = pd.read_csv('${params.single_transcript_mRNA_csv}', sep=',', index_col=0)

def load_genomes(UTRfilestring, twobitfile):
	UTRdict= rph.readindict(open(UTRfilestring, "rU"))
	genome= twobitreader.TwoBitFile(twobitfile) 
	return UTRdict, genome

def build_avggene_first(UTRdict, region_choice, upstreamNTs, downstreamNTs):
    fp_assign_path = "${name}"
    samp = "${name}"
    trspfilestring = '%s/%s/density%sp_%sshift_%sto%s/%s_%sto%sf/%s_%sto%sf_' %(fp_assign_path, densitystring, assignment, ribosome_shift, minlen, maxlen, samp, minlen, maxlen, samp, minlen, maxlen)
    totreads_countfile = "${name}/${name}_FPassigned_counts.txt"
    totreadcountf = open(totreads_countfile, "r")
    totreads = int(totreadcountf.read())
    totreadcountf.close()
    ### This is where all the count files are loaded into a dictionary
    trspdict= rph.readcountsf(trspfilestring)
    ## create a list of 0's that is the length of the region of interest
    ## add 3 to account for the stop codon
    ## this will be added to for every transcript
    averagegene= [0 for num in range(0, 3+upstreamNTs+ downstreamNTs)] # add 3 for start or stop codon
    ## add counters and set to zero
    noUTRentry = 0 ### discard transcripts not in UTRdict
    zeroCdsdense = 0 ### discard transcripts with zero reads in CDS
    lowCdsdense = 0 ## Optional CDS density thresholding, in RPKM
    regionTooShort = 0 ### not enough 3'UTR in region of interest past first stop codon
    totalCountedTranscripts = 0 ## number included in final output
    
    ### calculate mRNA-region densities
    defaultInsets = { 'utr5Inset3' : 6, 'cdsInset5' : 18, 'cdsInset3' : 15, 'utr3Inset5' : 6 }
    zeroInsets    = { 'utr5Inset3' : 0, 'cdsInset5' : 0, 'cdsInset3' : 0, 'utr3Inset5' : 0 }
    customInsets  = { 'utr5Inset3' : 15, 'cdsInset5' : 24, 'cdsInset3' : 15, 'utr3Inset5' : 15 }
    
    if inset_choice == "default":
        insets = defaultInsets
    elif inset_choice == "zero":
        insets = zeroInsets
    elif inset_choice == "custom":
        insets = customInsets
    else:
        print "Insets were not set"
        sys.exit()
    
    ### Iterated through transcripts one at a time, retrieving counts in region of interest:
    for trsp in trspdict:
        ### Load in count file for the transcript here
        exonsplicedcounts = trspdict[trsp]
        
        if UTRdict.has_key(trsp)!=True: # check to make sure density file has an annotation in the UTR csv
            noUTRentry +=1
            continue
        
        mrnalen = int(UTRdict[trsp][3])
        cdslen = int(UTRdict[trsp][4])
        utr5len = int(UTRdict[trsp][5])
        utr3len = int(UTRdict[trsp][6])
        assert mrnalen == cdslen + utr5len + utr3len ## check that this is true
        
        ### define Coding sequence here
        cdsstart = utr5len
        cdsend = len(exonsplicedcounts) - utr3len # cdsend == first position of utr3
        if cdsstart == cdsend:
            print "Error, gene length is 0 for transcript "+ trsp
            sys.exit()
        
        ### Calculate Region Densities ###
        utr5lenMod = utr5len-insets['utr5Inset3']
        cdslenMod = cdslen-insets['cdsInset5']-insets['cdsInset3']
        utr3lenMod = utr3len-insets['utr3Inset5']
        mrnalenMod = utr5lenMod+cdslenMod+utr3lenMod
        
        utr5Counts = sum(exonsplicedcounts[:cdsstart-insets['utr5Inset3']])
        cdsCounts = sum(exonsplicedcounts[cdsstart+insets['cdsInset5']:cdsend-insets['cdsInset3']])
        utr3Counts = sum(exonsplicedcounts[cdsend+insets['utr3Inset5']:])
        mrnaCounts = utr5Counts+cdsCounts+utr3Counts
        
        ### RAW counts
        RAWutr5Counts = int(utr5Counts*(totreads/1E6))
        RAWutr3Counts = int(utr3Counts*(totreads/1E6))
        RAWcdsCounts = int(cdsCounts*(totreads/1E6))
        RAWmrnaCounts = int(mrnaCounts*(totreads/1E6))
        ### denisites
        # mrnaDensity = (mrnaCounts/mrnalenMod)
        cdsDensity = (cdsCounts/cdslenMod)
        # utr5Density = (utr5Counts/utr5lenMod)
        # utr3Density = (utr3Counts/utr3lenMod)
        #### RPKM densities
        # mrnaDensity_rpkm = (mrnaCounts/mrnalenMod) * 1000
        cdsDensity_rpkm = (cdsCounts/cdslenMod) * 1000
        # utr5Density_rpkm = (utr5Counts/utr5lenMod) * 1000
        # utr3Density_rpkm = (utr3Counts/utr3lenMod) * 1000
        ### throw out zero's
        if cdsDensity == 0:
            zeroCdsdense += 1
            continue
        if cdsDensity*float(1000)< int(threshold):	# Threshold on cds density: (thresholding on "rpkm")
            lowCdsdense += 1
            continue
        ### define vector in valid CDS region, normalize by cdsDensity
        cdsSplicedCounts = exonsplicedcounts[cdsstart+insets['cdsInset5']:cdsend-insets['cdsInset3']]
        cdsNormCounts = [rpf/cdsDensity for rpf in cdsSplicedCounts] ## just region of cds within insets for counts, sum/len == 1
        # print cdslenMod, "modified cds length"
        # print sum(cdsNormCounts), "total normalized counts, should equal length of cds"
        # print sum(cdsNormCounts)/len(cdsNormCounts) ## this should == 1
        ### define vector for whole transcript, normlaized by cdsDensity
        exonNormCounts = [rpf/cdsDensity for rpf in exonsplicedcounts] ## counts normalized by cds density, using only region within insets
        # print exonNormCounts
        # print sum(exonNormCounts)
        # print sum(exonNormCounts)/len(exonNormCounts) ## should typically be less than 1, unless greater density from utrs and start/stop codons
        ### for start codon metagenes
        if region_choice == 'start':
            ### check boundaries:
            # print mRNAdf.loc[trsp]['mRNAseqs'][cdsstart:cdsstart+3] ### start codon
            if len(exonNormCounts[cdsstart-upstreamNTs:cdsstart+3+downstreamNTs]) < upstreamNTs+3+downstreamNTs:
                regionTooShort +=1
                continue
            else:
                totalCountedTranscripts +=1
                avggene_counts = exonNormCounts[cdsstart-upstreamNTs:cdsstart+3+downstreamNTs]
                for i in range(len(avggene_counts)): ### add these counts to the running total
                    averagegene[i] += avggene_counts[i]
        
        if region_choice == 'stop':
            ### check boundaries
            # print mRNAdf.loc[trsp]['mRNAseqs'][cdsend-3:cdsend] ### stop codon
            if len(exonNormCounts[cdsend-3-upstreamNTs:cdsend+downstreamNTs]) < upstreamNTs+3+downstreamNTs:
                regionTooShort +=1
                continue
            else:
                totalCountedTranscripts +=1
                # print upstreamNTs+3+downstreamNTs
                # print len(exonNormCounts[cdsend-3-upstreamNTs:cdsend+downstreamNTs])
                avggene_counts = exonNormCounts[cdsend-3-upstreamNTs:cdsend+downstreamNTs]
                for i in range(len(avggene_counts)): ### add these counts to the running total
                    averagegene[i] += avggene_counts[i]
    
    averagegene_equal = [rpf/totalCountedTranscripts for rpf in averagegene] ### divide by total number of valid transcripts
    # print averagegene
    # print averagegene_equal
    positions = range(-upstreamNTs-1, downstreamNTs+2) # start or stop codon is [-1,0,1]
    df = pd.DataFrame({'position':positions, 'avg':averagegene_equal})
    df = df[['position', 'avg']]
    # print df['avg'].sum()/len(df)
    # print df 
    print "Avggene run compolete for sample %s" % samp
    print "Transcripts inculded %s" % totalCountedTranscripts
    print "Number of transcripts absent in UTRfile: %s" % noUTRentry
    print "Number of transcripts with zero CDS density: %s" % zeroCdsdense
    print "Number of transcripts below CDS density threshold: %s" % lowCdsdense 
    print "Number of transcripts too short avggene region: %s" % regionTooShort
    print "- - - - - - - -"
    print df
    print "- - - - - - - -"
    
    if region_choice == 'start':
        alignpos = "1"
    elif region_choice == 'stop':
        alignpos = "2"
    else:
        print "alignpos not set!!!"
    
    ### write csv file
    fp_assign_path = "${name}"
    avggene_csv_path = "%s/avggene%s_ORF%s_%sshift_%s%s150" % (fp_assign_path, alignpos, norm_type, ribosome_shift, assignment, norm)
    if not os.path.exists(avggene_csv_path):	os.makedirs(avggene_csv_path)
    # adding cdsNorm to indicate normalized densities
    csv_outfile = "%s/%s_%s_rpkmThresh%s_%sto%sf_avg_%s_cdsNorm.csv" %(avggene_csv_path, samp, pop, threshold, minlen, maxlen, alignpos)
    df.to_csv(csv_outfile, index=False)

def main():
    os.system("mv ${dirName} ${name}")
    UTRfilestring = "${params.single_transcript_UTR_csv}"
    twobitfile = "${params.genome2bit}"
    UTRdict, genomes = load_genomes(UTRfilestring, twobitfile)
    build_avggene_first(UTRdict, "start", 100, 150)
    build_avggene_first(UTRdict, "stop", 150, 100)
    os.system("mv ${name} ${newdirName}")

if __name__ == '__main__':
    main()
"""
}

params.motif_file_root_path =  ""  //* @input
params.codon_occupancy_main_path =  ""  //* @input
params.single_transcript_UTR_csv =  ""  //* @input

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 10
    $MEMORY = 30
}
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 240
    $CPU  = 10
    $MEMORY = 30
    $QUEUE = "short"
} 
//* platform
//* autofill

process codon_occupancies {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /${name}$/) "results/$filename"
}

input:
 set val(name), file(dirName) from g_9_outputDir_g_12

output:
 set val(name), file("${newdirName}")  into g_12_outputDir_g_6
 file "${name}"  into g_12_resultsdir

script:
newdirName = name + "codon_occupancies"
"""
#!/usr/bin/env python

import sys
import os
import numpy as np
import argparse
import importlib
from pathos.multiprocessing import ProcessingPool as Pool

motiffilerootpath = "${params.motif_file_root_path}"

motifs= " \\"['TTT','TTC','TTA','TTG','CTT','CTC','CTA','CTG','ATT','ATC','ATA','ATG','GTT','GTC','GTA','GTG','TCT','TCC','TCA','TCG','CCT','CCC','CCA','CCG','ACT','ACC','ACA','ACG','GCT','GCC','GCA','GCG','TAT','TAC','CAT','CAC','CAA','CAG','AAT','AAC','AAA','AAG','GAT','GAC','GAA','GAG','TGT','TGC','TGG','AGT','AGC','GGT','GGC','GGA','GGG','AGA','AGG','CGT','CGC','CGA','CGG','TAA','TAG','TGA']\\" "
cds5trim= '15'
cds3trim= '15'
seqwin= ' "[\\'0\\',\\'3\\']" ' # For density files shifted to the first nt of A sites, should be multiple of 3
exclusionfiles=  ' "[\\'0\\']" '
norm_type = "rpm"
flmin = 28
flmax = 35 # set to exact readsize that is desired
eAmin = 21
eAmax = 24
eEmin = 18
eEmax = 19
UTRfilestring = "${params.single_transcript_UTR_csv}"

## define handler function:

def get_codon_occ_fl(file):
	#function to get codon occupancies
	#footprint_pop: fl, eA, eE, or aL
	#ribosome_site: A, P, E, or 0 
	footprint_pop = 'fl'
	ribosome_site = 'A'
	fiveorthreeprime = '5'
	### set normalization type:
	if norm_type == "raw":
		densitystring = "DensityUnnormalized"
	elif norm_type == "rpm":
		densitystring = "Density_rpm"
	else:
		print "Normalization is not set!"
		sys.exit()

	### set footprint population:
	if footprint_pop == 'fl':
		readsizemin = str(flmin)
		readsizemax = str(flmax)
	elif footprint_pop == 'eA':
		readsizemin = str(eAmin)
		readsizemax = str(eAmax)
	elif footprint_pop == 'eE':
		readsizemin = str(eEmin)
		readsizemax = str(eEmax)
	elif footprint_pop == 'aL':
		readsizemin = 'aL'

	### density file paths:
	
	trspdictfilestring = '%s/%s/density%sp_%sshift_%sto%s/%s_%sto%sf/%s_%sto%sf_' % (
		file, densitystring, fiveorthreeprime, ribosome_site, readsizemin, 
		readsizemax, file, readsizemin, readsizemax,file,readsizemin,readsizemax)
	
	sample_name = file

	outlistfilepath = '%s/codon' % (file)
	if not os.path.exists(outlistfilepath): os.makedirs(outlistfilepath)
	outlistfile = '%s/%s_%sshift_%sto%s__%soccupancy_%scds5trim_%scds3trim_codonOccupancy.csv' % (
		outlistfilepath, file, ribosome_site, readsizemin, readsizemax, fiveorthreeprime, cds5trim, cds3trim)
	outfileparams = '%s/%s_%sshift_%sto%s__%soccupancy_%scds5trim_%scds3trim_codonOccupancy_params.txt' % (
		outlistfilepath, file, ribosome_site, readsizemin, readsizemax, fiveorthreeprime, cds5trim, cds3trim)

	commandstring= 'python ${params.codon_occupancy_main_path} \
						--motiffilerootpath %s \
						--motifs %s \
						--trspdictfilestring %s \
						--sample_name %s \
						--UTRfilestring %s \
						--cds5trim %s \
						--cds3trim %s \
						--seqwin %s \
						--exclusionfiles %s \
						--outfileparams %s \
						--outlistfile %s' % (

						motiffilerootpath, 
						motifs, 
						trspdictfilestring, 
						sample_name, 
						UTRfilestring, 
						cds5trim, 
						cds3trim, 
						seqwin, 
						exclusionfiles, 
						outfileparams, 
						outlistfile)
	print commandstring
	os.system(commandstring)


def get_codon_occ_eA(file):
	# function to get codon occupancies
	# footprint_pop: fl, eA, eE, or aL
	# ribosome_site: A, P, E, or 0 although zero is not very informative here
	print "starting sample %s" % (file)

	footprint_pop = 'eA'
	ribosome_site = 'A'
	fiveorthreeprime = '5'

	### set normalization type:
	if norm_type == "raw":
		densitystring = "DensityUnnormalized"
	elif norm_type == "rpm":
		densitystring = "Density_rpm"
	else:
		print "Normalization is not set!"
		sys.exit()

	### set footprint population:
	if footprint_pop == 'fl':
		readsizemin = str(flmin)
		readsizemax = str(flmax)
	elif footprint_pop == 'eA':
		readsizemin = str(eAmin)
		readsizemax = str(eAmax)
	elif footprint_pop == 'eE':
		readsizemin = str(eEmin)
		readsizemax = str(eEmax)
	elif footprint_pop == 'aL':
		readsizemin = 'aL'

	# what is codonocc reading these in as?
	### density file paths:
	
	trspdictfilestring = '%s/%s/density%sp_%sshift_%sto%s/%s_%sto%sf/%s_%sto%sf_' % (
		file, densitystring, fiveorthreeprime, ribosome_site, readsizemin, 
		readsizemax, file, readsizemin, readsizemax,file,readsizemin,readsizemax)
	
	sample_name = file

	outlistfilepath = '%s/codon' % (file)
	if not os.path.exists(outlistfilepath): os.makedirs(outlistfilepath)
	outlistfile = '%s/%s_%sshift_%sto%s__%soccupancy_%scds5trim_%scds3trim_codonOccupancy.csv' % (
		outlistfilepath, file, ribosome_site, readsizemin, readsizemax, fiveorthreeprime, cds5trim, cds3trim)
	outfileparams = '%s/%s_%sshift_%sto%s__%soccupancy_%scds5trim_%scds3trim_codonOccupancy_params.txt' % (
		outlistfilepath, file, ribosome_site, readsizemin, readsizemax, fiveorthreeprime, cds5trim, cds3trim)

	commandstring= 'python ${params.codon_occupancy_main_path} \
						--motiffilerootpath %s \
						--motifs %s \
						--trspdictfilestring %s \
						--sample_name %s \
						--UTRfilestring %s \
						--cds5trim %s \
						--cds3trim %s \
						--seqwin %s \
						--exclusionfiles %s \
						--outfileparams %s \
						--outlistfile %s' % (

						motiffilerootpath, 
						motifs, 
						trspdictfilestring, 
						sample_name, 
						UTRfilestring, 
						cds5trim, 
						cds3trim, 
						seqwin, 
						exclusionfiles, 
						outfileparams, 
						outlistfile)

	print commandstring
	os.system(commandstring)


def main():
    os.system("mv ${dirName} ${name}")
    samplelist = ["${name}"]
    p = Pool(nodes = ${thread})
    p.map(get_codon_occ_fl, samplelist)
    p2 = Pool(nodes = ${thread})
    p2.map(get_codon_occ_eA, samplelist)
    os.system("mv ${name} ${newdirName}")
    os.system("mkdir ${name}")
    os.system("rsync -avzu --exclude='Density_rpm' --exclude='DensityUnnormalized' ${newdirName}/ ${name}/")

if __name__ == '__main__':
	main()
"""
}

sample_order = params.aggr.sample_order
amino_acid_list = params.aggr.amino_acid_list
color_code_list = params.aggr.color_code_list
control_group_name = params.aggr.control_group_name
control_group = params.aggr.control_group
treatment_group_name = params.aggr.treatment_group_name
treatment_group = params.aggr.treatment_group
//* @style @array:{treatment_group_name, treatment_group} @multicolumn:{control_group_name,control_group},{treatment_group_name,treatment_group}

process aggr {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /${name}\/countTables\/.*.csv$/) "count_tables/$filename"
}

input:
 set val(name), file(dirName) from g_12_outputDir_g_6

output:
 file "${name}"  into g_6_outputDir_g_7, g_6_outputDir_g_10, g_6_outputDir_g_13, g_6_outputDir_g_15, g_6_outputDir_g_17, g_6_outputDir_g_52
 set val(name), file("${name}/countTables/*.csv")  into g_6_outputFileTab

"""
mv ${dirName} ${name}


"""
}


process figure4B {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*.pdf$/) "reports/$filename"
	else if (filename =~ /.*.sh$/) "figure_data/$filename"
	else if (filename =~ /.*.csv$/) "figure_data/$filename"
}

input:
 file nameAll from g_6_outputDir_g_10.collect()

output:
 file "*.pdf"  into g_10_outputFilePdf
 file "*.sh"  into g_10_script
 file "*.csv"  into g_10_csvout

script:
nameAll = nameAll.collect{ '"' + it + '"'}
"""
#!/usr/bin/env python
# script is based off of riboseq_readingFrame_utr3meta
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt 
plt.rcParams['pdf.fonttype'] = 42 # this keeps most text as actual text in PDFs, not outlines

## import dependencies
import sys, os
import math
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import scipy.stats as stats
pd.set_option('display.max_columns', 50)
import seaborn as sns
# from pylab import *
import argparse
import importlib
import matplotlib.cm as cm

samplelist = ${nameAll}
opt_sample_order = '${sample_order}'
new_sample_order = [x.strip() for x in opt_sample_order.split(',')]
if len(new_sample_order) > 1:
	samplelist = new_sample_order 


colorList = ['lightgray', 'gray', 'black']

def log_trans_b10(x):
	try:
		return math.log(x, 10)
	except:
		return float(-6.00)
#         return float("NaN")
def log_trans_b2(x):
	try:
		return math.log(x, 2)
	except:
		# return float("NaN")
		return float(-15.00) # set arbitrarily low value

def load_countTables():
	namelist = []
	dflist= []
	# dflu3 = []
	for samp in samplelist:
		dftemp = pd.read_csv('%s/countTables/%s_fl_rpm_28to35_countTable_rpkm_utr3adj.csv' % (samp, samp))
		dftemp['sampname'] = samp
		dftemp['cdsCountsLog2'] = dftemp['cdsCounts'].apply(log_trans_b2)
		dftemp['RAWcdsCountsLog2'] = dftemp['RAWcdsCounts'].apply(log_trans_b2)
		dftemp['utr3CountsLog2'] = dftemp['utr3Counts'].apply(log_trans_b2)
		dftemp['utr3AdjCountsLog2'] = dftemp['utr3AdjCounts'].apply(log_trans_b2)
		dftemp['utr3OccLog2'] = dftemp['utr3_occupancy'].apply(log_trans_b2)
		dftemp['RRTSlog2'] = dftemp['RRTS'].apply(log_trans_b2)

		dftemp = dftemp.loc[dftemp['RAWcdsCounts'] > 128]
		dflist.append(dftemp)
		namelist.append(samp)
		
	# ### combine into one master dataframe
	df = pd.concat(dflist, axis=0, ignore_index = True)
	return df, dflist, namelist


def plot_RRTS_boxplot(df, dflist, namelist):

	figout = "Fig4B.pdf"
	validsamps = namelist
	dfp = df.loc[df['sampname'].isin(validsamps)]
	dfp = dfp.sort_values('stopcodon', ascending=True)

	fig, ax = plt.subplots(figsize=(6, 6))
	ax = sns.boxplot(x=dfp['sampname'],  y=dfp['RRTS'], orient = 'v',
				 fliersize = 0, hue = dfp['stopcodon'], notch=True, order = validsamps, palette=colorList)

	ax = plt.gca()
	ax.set_ylim((0,1.1))
	plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
	plt.title('RRTS values of transcripts sorted by normal termination codon identity.')
	plt.xticks(rotation=60)
	plt.savefig(figout, format='pdf', bbox_inches = "tight")

def stat_table(df, dflist, namelist):
	columnsTitles = ['sampname', 'comparison', 'pval']
	dfres= pd.DataFrame(columns=columnsTitles)
	for samp in namelist:
		dfs = df.loc[df['sampname'].isin([samp])]
		dfs_TAA = dfs.loc[dfs['stopcodon'].isin(["TAA"])]['RRTS'].values
		dfs_TAG = dfs.loc[dfs['stopcodon'].isin(["TAG"])]['RRTS'].values
		dfs_TGA = dfs.loc[dfs['stopcodon'].isin(["TGA"])]['RRTS'].values
		zero_TAA = np.count_nonzero(dfs_TAA)
		zero_TAG = np.count_nonzero(dfs_TAG)
		zero_TGA = np.count_nonzero(dfs_TGA)
		pval_TAAvsTAG = "NA"
		pval_TAGvsTGA = "NA"
		pval_TAAvsTGA = "NA"
		
		# u_stat: The Mann-Whitney U statistic, equal to min(U for x, U for y).
        # pvalue: p-value assuming an asymptotic normal distribution.
		if ( zero_TAA == 0 and zero_TAG == 0 ) != True: u_stat_TAAvsTAG, pval_TAAvsTAG = stats.mannwhitneyu(dfs_TAA, dfs_TAG, alternative='two-sided')
		if ( zero_TAG == 0 and zero_TGA == 0 ) != True: u_stat_TAGvsTGA, pval_TAGvsTGA = stats.mannwhitneyu(dfs_TAG, dfs_TGA, alternative='two-sided')
		if ( zero_TAA == 0 and zero_TGA == 0 ) != True: u_stat_TAAvsTGA, pval_TAAvsTGA = stats.mannwhitneyu(dfs_TAA, dfs_TGA, alternative='two-sided')
		row1 = {'sampname':samp, 'comparison':'TAA vs TAG', 'pval':pval_TAAvsTAG}
		row2 = {'sampname':samp, 'comparison':'TAG vs TGA', 'pval':pval_TAGvsTGA}
		row3 = {'sampname':samp, 'comparison':'TAA vs TGA', 'pval':pval_TAAvsTGA}
		dfres = dfres.append(row1, ignore_index=True)
		dfres = dfres.append(row2, ignore_index=True)
		dfres = dfres.append(row3, ignore_index=True)
	dfres.to_csv("Figure4B_mannwhitney_twosided_stats.csv", index=False)

def main():
	df, dflist, namelist = load_countTables()
	print df.head()
	df.to_csv("Figure4B_data.csv", index=False)
	plot_RRTS_boxplot(df, dflist, namelist)
	stat_table(df, dflist, namelist)
	if os.path.exists(".command.sh"):
		os.system("cp .command.sh Fig4B.sh")

if __name__ == '__main__':
	main()

"""
}


process figures {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*.pdf$/) "reports/$filename"
	else if (filename =~ /.*.sh$/) "figure_data/$filename"
}

input:
 file nameAll from g_6_outputDir_g_7.collect()

output:
 file "*.pdf"  into g_7_outputFilePdf
 file "*.sh"  into g_7_script

script:
nameAll = nameAll.collect{ '"' + it + '"'}
"""
#!/usr/bin/env python
# This script is for taking input density files and plotting average genes around the stop codon
# can be used to overlay multiple samples on the same plot

## plotting
import matplotlib
matplotlib.use('Agg') ## set backend here
import matplotlib.pyplot as plt 
plt.rcParams['pdf.fonttype'] = 42 # this keeps most text as actual text in PDFs, not outlines

## import dependencies
import sys
from Bio import SeqIO
import twobitreader
import GFF
import csv
import os
import pandas as pd
import argparse
import importlib
import random
import matplotlib.cm as cm

### inputs:
alignposition = "2" 	# set 1 for start codon and 2 for stop codon
pop = "fl" 	# fl, eA, or eE for density files
customSize = 30
threshold = 10
ribosome_site = "A" 	# A, P, E, or 0
flmin = 28
flmax = 35
eAmin = 21
eAmax = 24
eEmin = 18
eEmax = 19
norm_type = "rpm"

samples = ${nameAll}
opt_sample_order = '${sample_order}'
new_sample_order = [x.strip() for x in opt_sample_order.split(',')]
if len(new_sample_order) > 1:
	samples = new_sample_order

print samples

sample_plot_names = samples
samples_plotted = '_vs_'.join(samples)

number_of_colors = len(samples)
col_list = []
for i in range(number_of_colors):
	col_list.append( cm.Spectral(i/10.,1))

# optional color list
color_code_list = '${color_code_list}'
new_code_list = [x.strip() for x in color_code_list.split(',')]
if len(new_code_list) > 1:
	col_list = new_code_list

if pop == "fl":
	minlen = str(flmin)
	maxlen = str(flmax)
elif pop == "eA":
	minlen = str(eAmin)
	maxlen = str(eAmax)
elif pop == "eE":
	minlen = str(eEmin)
	maxlen = str(eEmax)
elif pop == "custom":
	minlen = str(customSize)
	maxlen = str(customSize)
else:
	print "whoops, something went wrong here, horribly horribly wrong!"
	# pop = ftsize


### for already shifted samples
def avggene_riboshift_plot_overlay(alignposition, ribosome_site, fiveorthreeprime='5', normalization='uneq', threshold ="0"):
	# read_length_list = ftsize
	alignpos = alignposition # '1' for start, '2' for stop
	assignment = fiveorthreeprime # should be 5' mapped at this point
	norm = normalization # 'uneq' for no normalizaiton; 'eq' to give equal weight to all genes
	ribosome_shift = ribosome_site # 'A' or 'P' or maybe 'E' later

	sample_numb = len(samples)
	sample_index = range(0,len(samples))
	sample_dict = {}
	for i in sample_index:
		sample_dict[samples[i]] = sample_index[i]

	df_fl_list = []
	df_eA_list = []
	df_eE_list = []
	df_aL_list = []
	df_custom_list = []


	for file in samples:
		fp_assign_path = '%s' % (file)
		avggene_csv_path = "%s/avggene%s_ORF%s_%sshift_%s%s150" % (fp_assign_path, alignpos, norm_type, ribosome_shift, assignment, norm) # norm should be 'uneq' for now


		fl_avggene_csv = '%s/%s_fl_rpkmThresh%s_%sto%sf_avg_%s_cdsNorm.csv' % (avggene_csv_path, file, threshold, flmin, flmax, alignpos)
		# eA_avggene_csv = '%s/%s_eA_rpkmThresh%s_%sto%sf_avg_%s_cdsNorm.csv' % (avggene_csv_path, file, threshold, eAmin, eAmax, alignpos)
		# eE_avggene_csv = '%s/%s_eE_rpkmThresh%s_%sto%sf_avg_%s_cdsNorm.csv' % (avggene_csv_path, file, threshold, eEmin, eEmax, alignpos)
		# aL_avggene_csv = '%s/%s_rpkmThresh%s_aL_aLf_avg_%s_cdsNorm.csv' % (avggene_csv_path, file, threshold, alignpos)

		if pop == "custom":
			custom_avggene_csv = '%s/%s_%s_shiftCustom_rpkmThresh0_%sto%sf_avg_%s.csv' % (avggene_csv_path, file, customSize, minlen, maxlen, alignpos)
			print custom_avggene_csv
			custom_avggene_df = pd.read_csv(custom_avggene_csv, index_col = 0, header=0)
			df_custom_list.append(custom_avggene_df)

		fl_avggene_df = pd.read_csv(fl_avggene_csv, index_col = 0, header=0)
		# eA_avggene_df = pd.read_csv(eA_avggene_csv, index_col = 0, header=0)
		# eE_avggene_df = pd.read_csv(eE_avggene_csv, index_col = 0, header=0)
		# aL_avggene_df = pd.read_csv(aL_avggene_csv, index_col = 0, header = 0)

		df_fl_list.append(fl_avggene_df)
		# df_eA_list.append(eA_avggene_df)
		# df_eE_list.append(eE_avggene_df)
		# df_aL_list.append(aL_avggene_df)

	plot_outfile = "Fig2A.pdf"


	fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(6,4))


	counter = 0 ## for making offsets
	if pop == 'fl':
		for i in range(len(df_fl_list)):
			df_fl_list[i].plot.line(x=df_fl_list[i].index+counter, y=df_fl_list[i].columns.values[0], ax=ax, color = col_list[i], lw=1, use_index=True, label=sample_plot_names[i])
			# counter +=6
			# df_fl_list[i].plot.scatter(y=df_fl_list[i].columns.values[0], ax=ax, color = col_list[i], s=2, use_index=True)
	if pop == 'eA':
		for i in range(len(df_eA_list)):
			df_eA_list[i].plot.line(y=df_eA_list[i].columns.values[0], ax=ax, color = col_list[i], lw=1, use_index=True, label=sample_plot_names[i])
	if pop == 'custom':
		for i in range(len(df_fl_list)):
			df_custom_list[i].plot.line(x=df_custom_list[i].index+counter, y=df_custom_list[i].columns.values[0], ax=ax, color = col_list[i], lw=1, use_index=True, label=sample_plot_names[i])

	#plt.legend(loc=1, prop={'size': 6})
	plt.ylabel ('Normalized Reads')
	plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
	plt.annotate('Average gene (or metagene) analysis is performed.\\nAll transcripts at their annotated stop codons are aligned\\nand normalized ribosome densities are calculated in this window.', (0,0), (0, -40), xycoords='axes fraction', textcoords='offset points', va='top')
	plt.savefig(plot_outfile, format = 'pdf', bbox_inches = "tight")
	plt.close()

### for already shifted samples, for now
def avggene_riboshift_plot_overlay_zoom(alignposition, ribosome_site, fiveorthreeprime='5', normalization='uneq', threshold ="0"):
	# read_length_list = ftsize
	alignpos = alignposition # '1' for start, '2' for stop
	assignment = fiveorthreeprime # should be 5' mapped at this point
	norm = normalization # 'uneq' for no normalizaiton; 'eq' to give equal weight to all genes
	ribosome_shift = ribosome_site # 'A' or 'P' or maybe 'E' later

	sample_numb = len(samples)
	sample_index = range(0,len(samples))
	sample_dict = {}
	for i in sample_index:
		sample_dict[samples[i]] = sample_index[i]

	df_fl_list = []
	df_eA_list = []
	df_eE_list = []
	df_aL_list = []
	df_custom_list = []


	for file in samples:
		fp_assign_path = '%s' % (file)
		avggene_csv_path = "%s/avggene%s_ORF%s_%sshift_%s%s150" % (fp_assign_path, alignpos, norm_type, ribosome_shift, assignment, norm) # norm should be 'uneq' for now


		fl_avggene_csv = '%s/%s_fl_rpkmThresh%s_%sto%sf_avg_%s_cdsNorm.csv' % (avggene_csv_path, file, threshold, flmin, flmax, alignpos)
		# eA_avggene_csv = '%s/%s_eA_rpkmThresh%s_%sto%sf_avg_%s_cdsNorm.csv' % (avggene_csv_path, file, threshold, eAmin, eAmax, alignpos)
		# eE_avggene_csv = '%s/%s_eE_rpkmThresh%s_%sto%sf_avg_%s_cdsNorm.csv' % (avggene_csv_path, file, threshold, eEmin, eEmax, alignpos)
		# aL_avggene_csv = '%s/%s_rpkmThresh%s_aL_aLf_avg_%s_cdsNorm.csv' % (avggene_csv_path, file, threshold, alignpos)

		if pop == "custom":
			custom_avggene_csv = '%s/%s_%s_shiftCustom_rpkmThresh0_%sto%sf_avg_%s.csv' % (avggene_csv_path, file, customSize, minlen, maxlen, alignpos)
			print custom_avggene_csv
			custom_avggene_df = pd.read_csv(custom_avggene_csv, index_col = 0, header=0)
			df_custom_list.append(custom_avggene_df)

		fl_avggene_df = pd.read_csv(fl_avggene_csv, index_col = 0, header=0)
		# eA_avggene_df = pd.read_csv(eA_avggene_csv, index_col = 0, header=0)
		# eE_avggene_df = pd.read_csv(eE_avggene_csv, index_col = 0, header=0)
		# aL_avggene_df = pd.read_csv(aL_avggene_csv, index_col = 0, header = 0)

		df_fl_list.append(fl_avggene_df)
		# df_eA_list.append(eA_avggene_df)
		# df_eE_list.append(eE_avggene_df)
		# df_aL_list.append(aL_avggene_df)

	plot_outfile = "Fig2B.pdf"


	fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(4,4))


	ax.set_ylim(0, 0.35)
	ax.set_xlim(0,100)

	counter = 0 ## for making offsets
	if pop == 'fl':
		for i in range(len(df_fl_list)):
			df_fl_list[i].plot.line(x=df_fl_list[i].index+counter, y=df_fl_list[i].columns.values[0], ax=ax, color = col_list[i], lw=1, use_index=True, label=sample_plot_names[i])
			# counter +=6
			# df_fl_list[i].plot.scatter(y=df_fl_list[i].columns.values[0], ax=ax, color = col_list[i], s=2, use_index=True)
	if pop == 'eA':
		for i in range(len(df_eA_list)):
			df_eA_list[i].plot.line(y=df_eA_list[i].columns.values[0], ax=ax, color = col_list[i], lw=1, use_index=True, label=sample_plot_names[i])
	if pop == 'custom':
		for i in range(len(df_fl_list)):
			df_custom_list[i].plot.line(x=df_custom_list[i].index+counter, y=df_custom_list[i].columns.values[0], ax=ax, color = col_list[i], lw=1, use_index=True, label=sample_plot_names[i])

	plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
	plt.annotate('Zoomed version of Fig 2A.', (0,0), (0, -40), xycoords='axes fraction', textcoords='offset points', va='top')
	plt.ylabel ('Normalized Reads')
	#plt.legend(loc=1, prop={'size': 6})
	plt.savefig(plot_outfile, format = 'pdf', bbox_inches = "tight")
	plt.close()

def main():
	avggene_riboshift_plot_overlay(alignposition, ribosome_site, normalization='eq', threshold=threshold)
	avggene_riboshift_plot_overlay_zoom(alignposition, ribosome_site, normalization='eq', threshold=threshold)
	if os.path.exists(".command.sh"):
		os.system("cp .command.sh Fig2AB.sh")
	

if __name__ == '__main__':
	main()
"""
}


process figure2C {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*.pdf$/) "reports/$filename"
	else if (filename =~ /.*.sh$/) "figure_data/$filename"
}

input:
 file nameAll from g_6_outputDir_g_52.collect()

output:
 file "*.pdf"  into g_52_outputFilePdf
 file "*.sh"  into g_52_script

script:
treatment_group_str = treatment_group.collect{ '"' + it + '"'}
treatment_group_name_str = treatment_group_name.collect{ '"' + it + '"'}
"""
#!/usr/bin/env python
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt 
plt.rcParams['pdf.fonttype'] = 42 # this keeps most text as actual text in PDFs, not outlines

## import dependencies
import sys
from Bio import SeqIO
import twobitreader
import GFF
import csv
import os
import pandas as pd
pd.set_option('display.max_columns', 50)
pd.set_option('display.max_rows', 300)
import seaborn as sns
import argparse
import importlib
import matplotlib.cm as cm

### inputs:
alignposition = "2" 	# set 1 for start codon and 2 for stop codon
pop = 'fl'
ribosome_site = "A" 	# A, P, E, or 0
normalization = 'eq' # 'uneq' for rpm and 'eq' for rpkm
flmin = 28
flmax = 35
eAmin = 21
eAmax = 24
eEmin = 18
eEmax = 19

### Footprint assignment
norm_type = "rpm" # either 'raw' for no normalization or 'rpm' for normalization to all mapped reads in alignment BAM file
threshold = '0' 

## indexes
ctrl_sample = '${control_group}'
ctrl_sample_name = '${control_group_name}'
ctrl_sample_names = [ctrl_sample_name]
ctrl_samples = [x.strip() for x in ctrl_sample.split(',')]
treat_samples = ${treatment_group_str}
treat_samples_name = ${treatment_group_name_str}
samples = ctrl_samples + treat_samples
sample_names = ctrl_sample_names + treat_samples_name

#colorList
number_of_colors = len(samples)
colorList = []
for i in range(number_of_colors):
	colorList.append( cm.Spectral(i/10.,1))

# optional color list
color_code_list = '${color_code_list}'
new_code_list = [x.strip() for x in color_code_list.split(',')]
if len(new_code_list) > 1:
	colorList = new_code_list

# grouping samples for the plot
splittedSamples=[]
treatments=[]
colorVals=[]
for i in range(len(sample_names)):
	sampName = sample_names[i]
	colorVal = colorList[i]
	sampList = [x.strip() for x in samples[i].split(',')]
	for k in range(len(sampList)):
		splittedSamples.append(sampList[k])
		treatments.append(sampName)
		colorVals.append(colorVal)

print splittedSamples
print treatments
print colorVals
	

if pop == "fl":
	minlen = str(flmin)
	maxlen = str(flmax)
elif pop == "eA":
	minlen = str(eAmin)
	maxlen = str(eAmax)
elif pop == "eE":
	minlen = str(eEmin)
	maxlen = str(eEmax)
elif pop == "custom":
	minlen = customSize
	maxlen = customSize
else:
	print "whoops, something went wrong here, horribly horribly wrong!"


### for already shifted samples, for now
def avggene_riboshift_plot_overlay(alignposition, ribosome_site, normalization, fiveorthreeprime='5'):
	# read_length_list = ftsize
	alignpos = alignposition # '1' for start, '2' for stop
	assignment = fiveorthreeprime # should be 5' mapped at this point
	norm = normalization # 'uneq' for no normalization; 'eq' to give equal weight to all genes
	ribosome_shift = ribosome_site # 'A' or 'P' or maybe 'E' later

	df_fl_list = []
	df_eA_list = []
	df_eE_list = []
	df_aL_list = []

	df_custom_list = []
	
	for file in splittedSamples:
		fp_assign_path = file
		avggene_csv_path = "%s/avggene%s_ORF%s_%sshift_%s%s150" % (fp_assign_path, alignpos, norm_type, ribosome_shift, assignment, norm) # norm should be 'uneq' for now

	###

	## get paths to stored csv average gene files
		if pop == 'custom':
			if norm == 'uneq':
				custom_csv = '%s/%s_%s_shiftCustom_rpkmThresh0_%sto%sf_avg_%s.csv' % (avggene_csv_path, file, customSize, customSize, customSize, alignpos)
		if norm == 'uneq':
			fl_avggene_csv = '%s/%s_fl_rpkmThresh0_%sto%sf_avg_%s_cdsNorm.csv' % (avggene_csv_path, file, flmin, flmax, alignpos)
		elif norm == 'eq':
			fl_avggene_csv = '%s/%s_fl_rpkmThresh10_%sto%sf_avg_%s_cdsNorm.csv' % (avggene_csv_path, file, flmin, flmax, alignpos)

		else:
			print "no norm selected"

		fl_avggene_df = pd.read_csv(fl_avggene_csv, index_col = 0, header=0)
		df_fl_list.append(fl_avggene_df)

	#######

	counter = 0
	# for avdf in df_fl_list:
	relDenList = []

	if pop == 'fl':
		dfFrame = df_fl_list
	elif pop == 'eA':
		dfFrame = df_eA_list
	elif pop == 'eE':
		dfFrame = df_eE_list
	elif pop == 'custom':
		dfFrame = df_custom_list
	else:
		print "pop not set, no plot was made"

	for avdf in dfFrame:
		# print avdf

		# print avdf.loc[-148:-17] ## not zero based with loc
		# print len(avdf.loc[-148:-17])

		# print avdf.loc[5:100, 'avg'] ## not zero based with loc
		# print len(avdf.loc[5:100])

		### CDS - count regions from nt -148 to -17 
		### UTR3 - countr regiona from 5 to 100, inclusive


		rpfCounts = avdf['avg']

		cdsCounts = avdf.loc[-148:-17, 'avg']
		cdsDense = (sum(cdsCounts)/len(cdsCounts))

		utr3Counts = avdf.loc[5:100, 'avg'] ## do not include final position
		utr3Dense = (sum(utr3Counts)/len(utr3Counts))

		relDense = (utr3Dense/cdsDense)*100
		relDenList.append(relDense)

		counter +=1


	figout = "Fig2C.pdf"
	fig,ax = plt.subplots(figsize=(5,5))
	dfplt = pd.DataFrame(relDenList, index=splittedSamples, columns=['utrPer'])
	dfplt['treatments'] = treatments
	dfplt['colorVal'] = colorVals
	
	print relDenList
	# splittedSamples: ['c1bamtofq', 'e1bamtofq', 'e1bamtofq_dup']
	# relDenList:[0.7728843047744274, 4.729053958564587, 4.519081209823026]
	# return: av_vals
	# for each sample_names ["con", "exp1", "exp2"] take average of relDenList 
	av_vals=[]
	ind = 0;
	for i in range(len(sample_names)):
		sampList = [x.strip() for x in samples[i].split(',')]
		relDenLength = len(sampList)
		indStart = ind
		ind +=relDenLength
		renDenSplit = relDenList[indStart:ind]
		mean = sum(renDenSplit) / len(renDenSplit)
		av_vals.append(mean)
	print av_vals
	
	dfav = pd.DataFrame.from_dict({
			"tr":sample_names,
			'av':av_vals
		})


	### seaborn catplots:
	fig, ax = plt.subplots(figsize=(6,6))
	sns.boxplot(data=dfav, x='tr', y='av', showbox=False , width = 0.5, 
					 showcaps=False, color = "black")
	sns.swarmplot(data=dfplt, x='treatments', y='utrPer', size=8, ax=ax, palette=colorList)
	ax=plt.gca()

	for item in ax.get_xticklabels():
		item.set_rotation(90)
	ax.set_ylim(0, 25)
	plt.savefig(figout, format='pdf', bbox_inches = "tight")


##########
def main():
	avggene_riboshift_plot_overlay(alignposition, ribosome_site, normalization)
	if os.path.exists(".command.sh"):
		os.system("cp .command.sh Fig2C.sh")

if __name__ == '__main__':
    main()

"""
}

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 10
    $MEMORY = 30
}
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 240
    $CPU  = 10
    $MEMORY = 30
    $QUEUE = "short"
} 
//* platform
//* autofill

process figure3S1C {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*.pdf$/) "reports/$filename"
	else if (filename =~ /.*.sh$/) "figure_data/$filename"
}

input:
 file nameAll from g_6_outputDir_g_17.collect()

output:
 file "*.pdf"  into g_17_outputFilePdf
 file "*.sh"  into g_17_script

script:
nameAll = nameAll.collect{ '"' + it + '"'}
"""
#!/usr/bin/env python
# Script for plotting read size distributions for each region of an mRNA
# Can then look at 5'UTR vs CDS vs 3'UTR vs mRNA
# Relies on inputs from densebuilder: Unnormalized, 5'Mapped

import matplotlib
matplotlib.use('Agg') # set backend for matplotlib
import matplotlib.pyplot as plt 
plt.rcParams['pdf.fonttype'] = 42 # this keeps most text as actual text in PDFs, not outlines

import sys, os
sys.path.append("${params.rphelper_path}")
import pandas as pd
import pysam
import numpy as np
import csv
from datetime import datetime
import argparse
from pathos.multiprocessing import ProcessingPool as Pool
pd.set_option('display.max_columns', 40)
import argparse
import importlib
import rphelper as rph
import matplotlib.cm as cm

samplelist = ${nameAll}
opt_sample_order = '${sample_order}'
new_sample_order = [x.strip() for x in opt_sample_order.split(',')]
if len(new_sample_order) > 1:
	samplelist = new_sample_order 
ftsize = [15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40]


### set inputs here:
defaultInsets = { 'utr5Inset3' : 6, 'cdsInset5' : 18, 'cdsInset3' : 15, 'utr3Inset5' : 6 }
zeroInsets    = { 'utr5Inset3' : 0, 'cdsInset5' : 0, 'cdsInset3' : 0, 'utr3Inset5' : 0 }
# insets= defaultInsets
insets= zeroInsets
UTRfilestring = "${params.single_transcript_UTR_csv}"
UTRdict = rph.readindict(open(UTRfilestring, "rU"))

def region_size_dist_ftsize(readsize, sample):
	### build the counts for a given region 
	fp_assign_path = '%s' % (sample)
	totreads_countfile = "%s/%s_FPassigned_counts.txt" % (fp_assign_path, sample)
	totreadcountf = open(totreads_countfile, "r")
	totreads = int(totreadcountf.read())
	totreadcountf.close()
	print "total reads for sample %s = %s" % (sample, totreads)

	# for readsize in ftsize:
	readsize = str(readsize) # convert to string
	trspdictfilestring = '%s/DensityUnnormalized/density5p_0shift_%s/%s_%sf/%s_%sf_' %(
		fp_assign_path, readsize, sample, readsize, sample, readsize)

	bamfilepath_readsize = '%s/%s_star_default/%s_%s_match.sorted.bam' % (
		fp_assign_path, sample, sample, readsize)
	bamfile = pysam.AlignmentFile(bamfilepath_readsize, 'rb')
	read_count_bam = bamfile.count()
	print "total reads in bamfile for sample: %s, read length: %s, equals == %s" % (sample, readsize, read_count_bam)

	## build the trspdict now for a given readlength:
	trspdict = rph.readcountsf(trspdictfilestring)

	## add counters 
	noUTRentry = 0
	zeroUtrlen = 0
	zeroUtrlenInsets = 0
	zeroCdsdense = 0
	lowCdsdense = 0
	lowCdsCounts = 0 # adding cds raw read counter

	### counters for output
	totUtr5Counts = 0
	totCdsCounts = 0
	totUtr3Counts = 0
	totMrnaCounts = 0

	## iterate through every transcript in the gtf file
	for trsp in trspdict:
		if UTRdict.has_key(trsp)!=True: # check to make sure density file has an annotation in the UTR csv
			noUTRentry +=1
			continue

	# define base region sizes from UTRdict
		mrnalen = int(UTRdict[trsp][3])
		cdslen = int(UTRdict[trsp][4])
		utr5len = int(UTRdict[trsp][5])
		utr3len = int(UTRdict[trsp][6])

		### Not sure if I want to keep this here... see how many have lengths of zero first
		if utr5len == 0:
			zeroUtrlen +=1
			# print("transcript has zero utr5 len %s") % trsp
			# sys.exit() 
			continue
		if utr3len == 0:
			zeroUtrlen +=1
			continue

# get counts from density file
		exonsplicedcounts = trspdict[trsp]

		# set starts and ends 
		cdsstart = utr5len
		cdsend = len(exonsplicedcounts) - utr3len
		if cdsstart == cdsend:
			print "Error, gene length is 0 for transcript %s" % trsp
			sys.exit()


		# modify region lengths using insets:
		
		utr5len = utr5len-insets['utr5Inset3']
		cdslen = cdslen-insets['cdsInset5']-insets['cdsInset3']
		utr3len = utr3len-insets['utr3Inset5']
		mrnalen = utr5len+cdslen+utr3len

		if utr5len == 0:
			zeroUtrlenInsets +=1
			# print "transcript has zero utr5 len %s" % trsp
			# sys.exit() 
			continue
		if utr3len == 0:
			zeroUtrlenInsets +=1
			continue

		utr5Counts = sum(exonsplicedcounts[:cdsstart-insets['utr5Inset3']])
		utr3Counts = sum(exonsplicedcounts[cdsend+insets['utr3Inset5']:])
		cdsCounts = sum(exonsplicedcounts[cdsstart+insets['cdsInset5']:cdsend-insets['cdsInset3']])
		mrnaCounts = utr5Counts+cdsCounts+utr3Counts

		totUtr5Counts += utr5Counts
		totCdsCounts += cdsCounts
		totUtr3Counts += utr3Counts
		totMrnaCounts += mrnaCounts
	print "UTR5total = %s, CDStotal = %s, UTR3total = %s, mRNAtotal = %s " % (totUtr5Counts, totCdsCounts, totUtr3Counts, totMrnaCounts)
	return [int(readsize), totUtr5Counts, totCdsCounts, totUtr3Counts, totMrnaCounts]



def main():
	if os.path.exists(".command.sh"):
		os.system("cp .command.sh Fig3S1C.sh")
	counter = 0
	treatlist = samplelist
	for sample in samplelist:
		treat = treatlist[counter]
		sl = [sample] * len(ftsize)
		p = Pool(nodes=40) 
		results = p.map(region_size_dist_ftsize, ftsize, sl)

		dfcols = ['readLength', 'utr5', 'cds', 'utr3', 'mrna']
		df = pd.DataFrame.from_records(results, columns = dfcols)
		df = df.set_index('readLength')

		utr5reads = df['utr5'].sum()
		cdsreads = df['cds'].sum()
		utr3reads = df['utr3'].sum()
		mrnareads = df['mrna'].sum()

		df['utr5'] = df['utr5']/utr5reads
		df['cds'] = df['cds']/cdsreads
		df['utr3'] = df['utr3']/utr3reads
		df['mrna'] = df['mrna']/mrnareads

		df.drop('utr5', axis=1, inplace=True)
		df.drop('mrna', axis=1, inplace=True)
		outfig = 'Fig3S1C_%s.pdf' % (treat)

		fig, ax = plt.subplots()

		df.plot.line()
		plt.xticks([15,20,25,30,35,40])
		plt.title(sample)
		plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
		plt.ylabel ('Percent of Reads')
		plt.title ('Read size distribution')
		plt.savefig(outfig, format='pdf', bbox_inches="tight")
		plt.close()
		counter +=1

if __name__ == '__main__':
	main()
"""
}


process figure2S3A {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*.pdf$/) "reports/$filename"
	else if (filename =~ /.*.sh$/) "figure_data/$filename"
}

input:
 file nameAll from g_6_outputDir_g_15.collect()

output:
 file "*.pdf"  into g_15_outputFilePdf
 file "*.sh"  into g_15_script

script:
nameAll = nameAll.collect{ '"' + it + '"'}
"""
#!/usr/bin/env python
# This script is for taking input density files and plotting average genes around the start codon
# can be used to overlay multiple samples on the same plot

## plotting
import matplotlib
matplotlib.use('Agg') ## set backend here
import matplotlib.pyplot as plt
plt.rcParams['pdf.fonttype'] = 42 # this keeps most text as actual text in PDFs, not outlines

## import dependencies
import sys
from Bio import SeqIO
import twobitreader
import GFF
import csv
import os
import pandas as pd
import argparse
import importlib
import random
import matplotlib.cm as cm


### inputs:
alignposition = "1"     # set 1 for start codon and 2 for stop codon
pop = "fl"      # fl, eA, or eE for density files
customSize = 30
threshold = 10
ribosome_site = "A"     # A, P, E, or 0
flmin = 28
flmax = 35
eAmin = 21
eAmax = 24
eEmin = 18
eEmax = 19
norm_type = "rpm"

samples = ${nameAll}
opt_sample_order = '${sample_order}'
new_sample_order = [x.strip() for x in opt_sample_order.split(',')]
if len(new_sample_order) > 1:
	samples = new_sample_order  

sample_plot_names = samples
samples_plotted = '_vs_'.join(samples)
number_of_colors = len(samples)
col_list = []
for i in range(number_of_colors):
	col_list.append( cm.Spectral(i/10.,1))

# optional color list
color_code_list = '${color_code_list}'
new_code_list = [x.strip() for x in color_code_list.split(',')]
if len(new_code_list) > 1:
	col_list = new_code_list


if pop == "fl":
    minlen = str(flmin)
    maxlen = str(flmax)
elif pop == "eA":
    minlen = str(eAmin)
    maxlen = str(eAmax)
elif pop == "eE":
    minlen = str(eEmin)
    maxlen = str(eEmax)
elif pop == "custom":
    minlen = str(customSize)
    maxlen = str(customSize)
else:
    print "whoops, something went wrong here, horribly horribly wrong!"
    # pop = ftsize


### for already shifted samples
def avggene_riboshift_plot_overlay(alignposition, ribosome_site, fiveorthreeprime='5', normalization='uneq', threshold ="0"):
    # read_length_list = ftsize
    alignpos = alignposition # '1' for start, '2' for stop
    assignment = fiveorthreeprime # should be 5' mapped at this point
    norm = normalization # 'uneq' for no normalizaiton; 'eq' to give equal weight to all genes
    ribosome_shift = ribosome_site # 'A' or 'P' or maybe 'E' later

    sample_numb = len(samples)
    sample_index = range(0,len(samples))
    sample_dict = {}
    for i in sample_index:
        sample_dict[samples[i]] = sample_index[i]

    df_fl_list = []
    df_eA_list = []
    df_eE_list = []
    df_aL_list = []
    df_custom_list = []


    for file in samples:
        fp_assign_path = '%s' % (file)
        avggene_csv_path = "%s/avggene%s_ORF%s_%sshift_%s%s150" % (fp_assign_path, alignpos, norm_type, ribosome_shift, assignment, norm) # norm should be 'uneq' for now


        fl_avggene_csv = '%s/%s_fl_rpkmThresh%s_%sto%sf_avg_%s_cdsNorm.csv' % (avggene_csv_path, file, threshold, flmin, flmax, alignpos)
        # eA_avggene_csv = '%s/%s_eA_rpkmThresh%s_%sto%sf_avg_%s_cdsNorm.csv' % (avggene_csv_path, file, threshold, eAmin, eAmax, alignpos)
        # eE_avggene_csv = '%s/%s_eE_rpkmThresh%s_%sto%sf_avg_%s_cdsNorm.csv' % (avggene_csv_path, file, threshold, eEmin, eEmax, alignpos)
        # aL_avggene_csv = '%s/%s_rpkmThresh%s_aL_aLf_avg_%s_cdsNorm.csv' % (avggene_csv_path, file, threshold, alignpos)

        if pop == "custom":
            custom_avggene_csv = '%s/%s_%s_shiftCustom_rpkmThresh0_%sto%sf_avg_%s.csv' % (avggene_csv_path, file, customSize, minlen, maxlen, alignpos)
            print custom_avggene_csv
            custom_avggene_df = pd.read_csv(custom_avggene_csv, index_col = 0, header=0)
            df_custom_list.append(custom_avggene_df)

        fl_avggene_df = pd.read_csv(fl_avggene_csv, index_col = 0, header=0)


        df_fl_list.append(fl_avggene_df)


    plot_outfile = "Fig2S3A.pdf"
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(6,4))


    counter = 0 ## for making offsets
    if pop == 'fl':
        for i in range(len(df_fl_list)):
            df_fl_list[i].plot.line(x=df_fl_list[i].index+counter,
                                    y=df_fl_list[i].columns.values[0],
                                    ax=ax, color = col_list[i],
                                    lw=1,
                                    use_index=True,
                                    label=sample_plot_names[i])
            counter +=6
    ax.set_xlim(-50,100)
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.ylabel ('Normalized Reads')
    plt.annotate('Average gene (or metagene) analysis is performed.\\nAll transcripts at their annotated start codons are aligned\\nand normalized ribosome densities are calculated in this window.', (0,0), (0, -40), xycoords='axes fraction', textcoords='offset points', va='top')
    plt.savefig(plot_outfile, format = 'pdf', bbox_inches = "tight")
    plt.close()


def main():
    avggene_riboshift_plot_overlay(alignposition, ribosome_site, normalization='eq', threshold=threshold)
    if os.path.exists(".command.sh"):
    	os.system("cp .command.sh Fig2S3A.sh")

if __name__ == '__main__':
    main()

"""
}


process figure2S3B {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*.pdf$/) "reports/$filename"
	else if (filename =~ /.*.sh$/) "figure_data/$filename"
	else if (filename =~ /.*.csv$/) "figure_data/$filename"
}

input:
 file nameAll from g_6_outputDir_g_13.collect()

output:
 file "*.pdf"  into g_13_outputFilePdf
 file "*.sh"  into g_13_script
 file "*.csv"  into g_13_csvout

script:
nameAll = nameAll.collect{ '"' + it + '"'}
treatment_group_str = treatment_group.collect{ '"' + it + '"'}
treatment_group_name_str = treatment_group_name.collect{ '"' + it + '"'}
"""
#!/usr/bin/env python
# This script plots codon occupancies relative to untreated cells
import matplotlib
matplotlib.use('Agg') ## set backend here
import matplotlib.pyplot as plt 
import matplotlib.patches as mpatches
plt.rcParams['pdf.fonttype'] = 42 # this keeps most text as actual text in PDFs, not outlines

import sys, os
import pandas as pd
import pysam
import numpy as np
import math
import argparse
import importlib
from scipy import stats

readlength = "28to35"
shift = "A"

## indexes
ctrl_sample = '${control_group}'
ctrl_sample_name = '${control_group_name}'
ctrl_samples = [x.strip() for x in ctrl_sample.split(',')]
treat_samples = ${treatment_group_str}
treat_samples_name = ${treatment_group_name_str}
print ctrl_samples
print treat_samples

amino_acid_list = '${amino_acid_list}'
aa_list = [x.strip() for x in amino_acid_list.split(',')]


colorDict = {
	'A':'#000075',
	'C':'#808080',
	'D':'#e6194b',
	'E':'#3cb44b',
	'F':'#ffe119',
	'G':'#4363d8',
	'H':'#f58231',
	'I':'#911eb4',
	'K':'#46f0f0',
	'L':'#f032e6',
	'M':'#bcf60c',
	'N':'#fabebe',
	'P':'#008080',
	'Q':'#e6beff',
	'R':'#9a6324',
	'S':'#fffac8',
	'T':'#800000',
	'V':'#aaffc3',
	'W':'#808000',
	'Y':'#ffd8b1'
}

codonDict = {
	'A':['GCT', 'GCC', 'GCA', 'GCG'],
	'C':['TGT', 'TGC'],
	'D':['GAT', 'GAC'],
	'E':['GAA', 'GAG'],
	'F':['TTT', 'TTC'],
	'G':['GGT', 'GGC', 'GGA', 'GGG'],
	'H':['CAT', 'CAC'],
	'I':['ATT', 'ATC', 'ATA'],
	'K':['AAA', 'AAG'],
	'L':['CTT', 'CTC', 'CTA', 'CTG', 'TTA', 'TTG'],
	'M':['ATG'],
	'N':['AAT', 'AAC'],
	'P':['CCT', 'CCC', 'CCA', 'CCG'],
	'Q':['CAA', 'CAG'],
	'R':['AGA','AGG','CGT','CGC','CGA','CGG'],
	'S':['TCT', 'TCC', 'TCA', 'TCG'],
	'T':['ACT', 'ACC', 'ACA', 'ACG'],
	'V':['GTT', 'GTC', 'GTA', 'GTG'],
	'W':['TGG'],
	'Y':['TAT', 'TAC']
}


def log_trans_b10(x):
	try:
		return math.log(x, 10)
	except:
		return float(-6.00)
#         return float("NaN")
def log_trans_b2(x):
	try:
		return math.log(x, 2)
	except:
#         return float("NaN")
		return float(-15.00) # set arbitrarily low value

def corrfunc(x, y, **kws):
	r, _ = stats.pearsonr(x, y)
	rho, _ = stats.spearmanr(x, y)
	ax = plt.gca()
	ax.annotate("R = {:.3f}".format(r),
				xy=(.1, .9), xycoords=ax.transAxes)


def plot_codons():
	counter = 0
	for i in range(len(treat_samples)):
		comp_sample = treat_samples[i]
		comp_sample_name = treat_samples_name[i]
		comp_samples = [x.strip() for x in comp_sample.split(',')]
		dflist = []
		sublist = [] #keep all sample names that are going to be compared
		cnt_log2 = [] #all control sample names + "_log2"
		treat_log2 = [] #all treatment sample names + "_log2"
		for i in range(len(ctrl_samples)):
			ctrl_samp = ctrl_samples[i]
			df = pd.read_csv("%s/codon/%s_%sshift_%s__5occupancy_15cds5trim_15cds3trim_codonOccupancy.csv" % (ctrl_samp, ctrl_samp, shift, readlength), index_col=0)
			dflist.append(df)
			sublist.append(ctrl_samp)
			cnt_log2.append(ctrl_samp+"_log2")
			
		for k in range(len(comp_samples)):
			comp_samp = comp_samples[k]
			df = pd.read_csv("%s/codon/%s_%sshift_%s__5occupancy_15cds5trim_15cds3trim_codonOccupancy.csv" % (comp_samp, comp_samp, shift, readlength), index_col=0)
			dflist.append(df)
			sublist.append(comp_samp)
			treat_log2.append(comp_samp+"_log2")
				
		dfout = pd.concat(dflist, axis=1)
		dfout.drop(index=['TAA', 'TAG', 'TGA'],axis=0,inplace=True)
		dfout = dfout[sublist].apply(pd.to_numeric)
		
		for col in dfout.columns:
			dfout[col+'_log2'] = dfout[col].apply(log_trans_b2)


		dfout['ctrl_mean'] = dfout[cnt_log2].mean(axis=1)
		dfout['ctrl_std'] = dfout[cnt_log2].std(axis=1)
		dfout['ctrl_sem'] = dfout[cnt_log2].sem(axis=1)

		dfout['treat_mean'] = dfout[treat_log2].mean(axis=1)
		dfout['treat_std'] = dfout[treat_log2].std(axis=1)
		dfout['treat_sem'] = dfout[treat_log2].sem(axis=1)
		dfout.to_csv("Figure2S3B_data_"+ctrl_sample_name+"_vs_"+comp_sample_name+".csv")


		### plotting:
		fig, ax = plt.subplots(figsize=(6,8))

		for aa in codonDict:
			for cdn in codonDict[aa]:

				plt.scatter(x=dfout.loc[cdn,'ctrl_mean'], 
							y=dfout.loc[cdn,'treat_mean'], 
							color="black", 
							edgecolor="black", 
							linewidth = 0.5, 
							s=20
						   )
				plt.errorbar(x=dfout.loc[cdn,'ctrl_mean'], 
							y=dfout.loc[cdn,'treat_mean'], 
							xerr = dfout.loc[cdn,'ctrl_std'],
							yerr = dfout.loc[cdn,'treat_std'],
							color="black", 
							ecolor="black", 
							linewidth = 0.75
							)	

		for aa in aa_list:
			for cdn in codonDict[aa]:
	#           
				plt.scatter(x=dfout.loc[cdn,'ctrl_mean'], 
							y=dfout.loc[cdn,'treat_mean'], 
							color=colorDict[aa], 
							edgecolor="black", 
							linewidth = 0.5, 
							s=20
						   )

				plt.text(x=dfout.loc[cdn,'ctrl_mean']+0.00, 
						 y=dfout.loc[cdn,'treat_mean']+0.00, 
						 color=colorDict[aa], 
						 s=aa, 
						 fontsize = 16)
				plt.text(x=dfout.loc[cdn,'ctrl_mean']+0.00, 
						 y=dfout.loc[cdn,'treat_mean']-0.00, 
						 color=colorDict[aa], #"black", 
						 s=cdn.replace('T', 'U'), 
						 fontsize = 10, 
						 va="top", 
						 ha="left")

		plt.xlabel ("Codon Occupancy - " + ctrl_sample_name)
		plt.ylabel ("Codon Occupancy - " + comp_sample_name)
		plt.annotate('Codon occupancies are plotted on the graph. Numbers on the axis\\nare 2 to the power of the value (e.g. 0.5 means 2^0.5). The average\\noccupancy between two samples is plotted for all 61 sense codons.\\nPearson correlations are displayed at the top.', (0,0), (0, -50), xycoords='axes fraction', textcoords='offset points', va='top')
		
		lims = [
			np.min([-1.0, -1.0]),
			np.max([1.5, 1.5]),
		]

		ax.plot(lims, lims, 'k--', alpha=0.75, zorder=0)
		ax.set_aspect('equal')
		ax.set_xlim(lims)
		ax.set_ylim(lims)
		corrfunc(dfout['ctrl_mean'], dfout['treat_mean'])

		outfile = "Fig2S3B_%s_vs_%s.pdf" % (ctrl_sample_name, comp_sample_name)
		plt.savefig(outfile, format="pdf")

		counter +=1


def main():
	plot_codons()
	if os.path.exists(".command.sh"):
		os.system("cp .command.sh Fig2S3B.sh")

if __name__ == '__main__':
	main()

"""
}

mappingListQuoteSep = mapList.collect{ '"' + it + '"'}.join(",") 
rawIndexList = indexList.collect{ '"' + it + '"'}.join(",") 
process Sequential_Mapping_Module_Sequential_Mapping_Dedup_Bam_count {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*.counts.tsv$/) "sequential_mapping_counts/$filename"
}

input:
 file bam from g42_32_bam_file_g42_27.collect()
 file index from g42_32_bam_index_g42_27.collect()

output:
 file "*.counts.tsv"  into g42_27_outputFileTSV

shell:
'''
#!/usr/bin/env perl
use List::Util qw[min max];
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;

my @header;
my %all_files;

my @mappingList = (!{mappingListQuoteSep});
my @rawIndexList = (!{rawIndexList});
my %indexHash;
my $dedup = "";
@indexHash{@mappingList} = @rawIndexList;

chomp(my $contents = `ls *.bam`);
my @files = split(/[\\n]+/, $contents);
foreach my $file (@files){
        $file=~/(.*)@(.*)_sorted(.*)\\.bam/;
        my $mapper = $1; 
        my $name = $2; ##header
        print $3;
        if ($3 eq ".dedup"){
            $dedup = "dedup.";
        }
        push(@header, $name) unless grep{$_ eq $name} @header; #mapped element header
        push @{$all_files{$mapper}}, $file;
}


open OUT, ">header.tsv";
print OUT join ("\\t", "id","len",@header),"\\n";
close OUT;

foreach my $key (sort keys %all_files) {  
    my @array = @{ $all_files{$key} };  
        unless (-e "$indexHash{$key}.bed") {
        print "2: bed not found run makeBed\\n";
            if (-e "$indexHash{$key}.fa") { 
                makeBed("$indexHash{$key}.fa", $key, "$indexHash{$key}.bed");
            } elsif(-e "$indexHash{$key}.fasta"){
                makeBed("$indexHash{$key}.fasta", $key, "$indexHash{$key}.bed");
            }
        }
    
        my $bamFiles = join ' ', @array;
        print "bedtools multicov -bams $bamFiles -bed $indexHash{$key}.bed >$key.${dedup}counts.tmp\\n";
        `bedtools multicov -bams $bamFiles -bed $indexHash{$key}.bed >$key.${dedup}counts.tmp`;
        my $iniResColumn = int(countColumn("$indexHash{$key}.bed")) + 1;
        `awk -F \\"\\\\t\\" \\'{a=\\"\\";for (i=$iniResColumn;i<=NF;i++){a=a\\"\\\\t\\"\\$i;} print \\$4\\"\\\\t\\"(\\$3-\\$2)\\"\\"a}\\' $key.${dedup}counts.tmp> $key.${dedup}counts.tsv`;
        `sort -k3,3nr $key.${dedup}counts.tsv>$key.${dedup}sorted.tsv`;
        `cat header.tsv $key.${dedup}sorted.tsv> $key.${dedup}counts.tsv`;
}

sub countColumn {
    my ( \$file) = @_;
    open(IN, \$file);
    my $line=<IN>;
    chomp($line);
    my @cols = split('\\t', $line);
    my $n = @cols;
    close OUT;
    return $n;
}

sub makeBed {
    my ( \$fasta, \$type, \$bed) = @_;
    print "makeBed $fasta\\n";
    print "makeBed $bed\\n";
    open OUT, ">$bed";
    open(IN, \$fasta);
    my $name="";
    my $seq="";
    my $i=0;
    while(my $line=<IN>){
        chomp($line);
        if($line=~/^>(.*)/){
            $i++ if (length($seq)>0);
            print OUT "$name\\t1\\t".length($seq)."\\t$name\\t0\\t+\\n" if (length($seq)>0); 
            $name="$1";
            $seq="";
        } elsif($line=~/[ACGTNacgtn]+/){
            $seq.=$line;
        }
    }
    print OUT "$name\\t1\\t".length($seq)."\\t$name\\t0\\t+\\n" if (length($seq)>0); 
    close OUT;
}

'''


}

mappingListQuoteSep = mapList.collect{ '"' + it + '"'}.join(",") 
rawIndexList = indexList.collect{ '"' + it + '"'}.join(",") 
process Sequential_Mapping_Module_Sequential_Mapping_Bam_count {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*.counts.tsv$/) "sequential_mapping_counts/$filename"
}

input:
 file bam from g42_32_bam_file_g42_23.collect()
 file index from g42_32_bam_index_g42_23.collect()

output:
 file "*.counts.tsv"  into g42_23_outputFileTSV

shell:
'''
#!/usr/bin/env perl
use List::Util qw[min max];
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;

my @header;
my %all_files;

my @mappingList = (!{mappingListQuoteSep});
my @rawIndexList = (!{rawIndexList});
my %indexHash;
my $dedup = "";
@indexHash{@mappingList} = @rawIndexList;

chomp(my $contents = `ls *.bam`);
my @files = split(/[\\n]+/, $contents);
foreach my $file (@files){
        $file=~/(.*)@(.*)_sorted(.*)\\.bam/;
        my $mapper = $1; 
        my $name = $2; ##header
        print $3;
        if ($3 eq ".dedup"){
            $dedup = "dedup.";
        }
        push(@header, $name) unless grep{$_ eq $name} @header; #mapped element header
        push @{$all_files{$mapper}}, $file;
}


open OUT, ">header.tsv";
print OUT join ("\\t", "id","len",@header),"\\n";
close OUT;

foreach my $key (sort keys %all_files) {  
    my @array = @{ $all_files{$key} };  
        unless (-e "$indexHash{$key}.bed") {
        print "2: bed not found run makeBed\\n";
            if (-e "$indexHash{$key}.fa") { 
                makeBed("$indexHash{$key}.fa", $key, "$indexHash{$key}.bed");
            } elsif(-e "$indexHash{$key}.fasta"){
                makeBed("$indexHash{$key}.fasta", $key, "$indexHash{$key}.bed");
            }
        }
    
        my $bamFiles = join ' ', @array;
        print "bedtools multicov -bams $bamFiles -bed $indexHash{$key}.bed >$key.${dedup}counts.tmp\\n";
        `bedtools multicov -bams $bamFiles -bed $indexHash{$key}.bed >$key.${dedup}counts.tmp`;
        my $iniResColumn = int(countColumn("$indexHash{$key}.bed")) + 1;
        `awk -F \\"\\\\t\\" \\'{a=\\"\\";for (i=$iniResColumn;i<=NF;i++){a=a\\"\\\\t\\"\\$i;} print \\$4\\"\\\\t\\"(\\$3-\\$2)\\"\\"a}\\' $key.${dedup}counts.tmp> $key.${dedup}counts.tsv`;
        `sort -k3,3nr $key.${dedup}counts.tsv>$key.${dedup}sorted.tsv`;
        `cat header.tsv $key.${dedup}sorted.tsv> $key.${dedup}counts.tsv`;
}

sub countColumn {
    my ( \$file) = @_;
    open(IN, \$file);
    my $line=<IN>;
    chomp($line);
    my @cols = split('\\t', $line);
    my $n = @cols;
    close OUT;
    return $n;
}

sub makeBed {
    my ( \$fasta, \$type, \$bed) = @_;
    print "makeBed $fasta\\n";
    print "makeBed $bed\\n";
    open OUT, ">$bed";
    open(IN, \$fasta);
    my $name="";
    my $seq="";
    my $i=0;
    while(my $line=<IN>){
        chomp($line);
        if($line=~/^>(.*)/){
            $i++ if (length($seq)>0);
            print OUT "$name\\t1\\t".length($seq)."\\t$name\\t0\\t+\\n" if (length($seq)>0); 
            $name="$1";
            $seq="";
        } elsif($line=~/[ACGTNacgtn]+/){
            $seq.=$line;
        }
    }
    print OUT "$name\\t1\\t".length($seq)."\\t$name\\t0\\t+\\n" if (length($seq)>0); 
    close OUT;
}

'''


}


process Sequential_Mapping_Module_Deduplication_Summary {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /deduplication_summary.tsv$/) "sequential_mapping/$filename"
}

input:
 file flagstat from g42_32_log_file_g42_30.collect()
 val mate from g_1_mate_g42_30

output:
 file "deduplication_summary.tsv"  into g42_30_outputFileTSV

shell:
'''
#!/usr/bin/env perl
use List::Util qw[min max];
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;

my @header;
my %all_files;
my %tsv;
my %headerHash;
my %headerText;

my $i=0;
chomp(my $contents = `ls *_duplicates_stats.log`);
my @files = split(/[\\n]+/, $contents);
foreach my $file (@files){
    $i++;
    $file=~/(.*)@(.*)@(.*)_duplicates_stats\\.log/;
    my $mapOrder = int($1); 
    my $mapper = $2; #mapped element 
    my $name = $3; ##sample name
    push(@header, $mapper) unless grep{$_ eq $mapper} @header; 
        
    # my $duplicates;
    my $aligned;
    my $dedup; #aligned reads after dedup
    my $percent=0;
    if ("!{mate}" eq "pair" ){
        #first flagstat belongs to first bam file
        chomp($aligned = `cat $file | grep 'properly paired (' | sed -n 1p | awk '{sum+=\\$1+\\$3} END {print sum}'`);
        #second flagstat belongs to dedup bam file
        chomp($dedup = `cat $file | grep 'properly paired (' | sed -n 2p | awk '{sum+=\\$1+\\$3} END {print sum}'`);
    } else {
        chomp($aligned = `cat $file | grep 'mapped (' | sed -n 1p | awk '{sum+=\\$1+\\$3} END {print sum}'`);
        chomp($dedup = `cat $file | grep 'mapped (' | sed -n 2p | awk '{sum+=\\$1+\\$3} END {print sum}'`);
    }
    # chomp($duplicates = `cat $file | grep 'duplicates' | awk '{sum+=\\$1+\\$3} END {print sum}'`);
    # $dedup = int($aligned) - int($duplicates);
    if ("!{mate}" eq "pair" ){
       $dedup = int($dedup/2);
       $aligned = int($aligned/2);
    } 
    $percent = "0.00";
    if (int($aligned)  > 0 ){
       $percent = sprintf("%.2f", ($aligned-$dedup)/$aligned*100); 
    } 
    $tsv{$name}{$mapper}=[$aligned,$dedup,"$percent%"];
    $headerHash{$mapOrder}=$mapper;
    $headerText{$mapOrder}=["$mapper (Before Dedup)", "$mapper (After Dedup)", "$mapper (Duplication Ratio %)"];
}

my @mapOrderArray = ( keys %headerHash );
my @sortedOrderArray = sort { $a <=> $b } @mapOrderArray;

my $summary = "deduplication_summary.tsv";
open(OUT, ">$summary");
print OUT "Sample\\t";
my @headArr = ();
for my $mapOrder (@sortedOrderArray) {
    push (@headArr, @{$headerText{$mapOrder}});
}
my $headArrAll = join("\\t", @headArr);
print OUT "$headArrAll\\n";

foreach my $name (keys %tsv){
    my @rowArr = ();
    for my $mapOrder (@sortedOrderArray) {
        push (@rowArr, @{$tsv{$name}{$headerHash{$mapOrder}}});
    }
    my $rowArrAll = join("\\t", @rowArr);
    print OUT "$name\\t$rowArrAll\\n";
}
close(OUT);
'''
}

//* autofill
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 240
    $CPU  = 1
    $MEMORY = 8
    $QUEUE = "short"
}
//* platform
//* autofill

process Sequential_Mapping_Module_Sequential_Mapping_Summary {

input:
 set val(name), file(bowfile) from g42_32_bowfiles_g42_26
 val mate from g_1_mate_g42_26
 val filtersList from g42_32_filter_g42_26

output:
 file '*.tsv'  into g42_26_outputFileTSV_g42_13
 val "sequential_mapping_sum"  into g42_26_name_g42_13

shell:
'''
#!/usr/bin/env perl
open(my \$fh, '>', "!{name}.tsv");
print $fh "Sample\\tGroup\\tTotal Reads\\tReads After Sequential Mapping\\tUniquely Mapped\\tMultimapped\\tMapped\\n";
my @bowArray = split(' ', "!{bowfile}");
my $group= "\\t";
my @filterArray = (!{filtersList});
foreach my $bowitem(@bowArray) {
    # get mapping id
    my @bowAr = $bowitem.split("_");
    $bowCount = $bowAr[0] + -1;
    # if bowfiles ends with underscore (eg. bow_rRNA), parse rRNA as a group.
    my ($RDS_In, $RDS_After, $RDS_Uniq, $RDS_Multi, $ALGN_T, $a, $b, $aPer, $bPer)=(0, 0, 0, 0, 0, 0, 0, 0, 0);
    if ($bowitem =~ m/bow_([^\\.]+)$/){
        $group = "$1\\t";
        open(IN, $bowitem);
        my $i = 0;
        while(my $line=<IN>){
            chomp($line);
            $line=~s/^ +//;
            my @arr=split(/ /, $line);
            $RDS_In=$arr[0] if ($i=~/^1$/);
            # Reads After Filtering column depends on filtering type
            if ($i == 2){
                if ($filterArray[$bowCount] eq "Yes"){
                    $RDS_After=$arr[0];
                } else {
                    $RDS_After=$RDS_In;
                }
            }
            if ($i == 3){
                $a=$arr[0];
                $aPer=$arr[1];
                $aPer=~ s/([()])//g;
                $RDS_Uniq=$arr[0];
            }
            if ($i == 4){
                $b=$arr[0];
                $bPer=$arr[1];
                $bPer=~ s/([()])//g;
                $RDS_Multi=$arr[0];
            }
            $ALGN_T=($a+$b);
            $i++;
        }
        close(IN);
    } elsif ($bowitem =~ m/star_([^\\.]+)$/){
        $group = "$1\\t";
        open(IN2, $bowitem);
        my $multimapped;
		my $aligned;
		my $inputCount;
		chomp($inputCount = `cat $bowitem | grep 'Number of input reads' | awk '{sum+=\\$6} END {print sum}'`);
		chomp($uniqAligned = `cat $bowitem | grep 'Uniquely mapped reads number' | awk '{sum+=\\$6} END {print sum}'`);
		chomp($multimapped = `cat $bowitem | grep 'Number of reads mapped to multiple loci' | awk '{sum+=\\$9} END {print sum}'`);
		## Here we exclude "Number of reads mapped to too many loci" from multimapped reads since in bam file it called as unmapped.
		## Besides, these "too many loci" reads exported as unmapped reads from STAR.
		$RDS_In = int($inputCount);
		$RDS_Multi = int($multimapped);
        $RDS_Uniq = int($uniqAligned);
        $ALGN_T = $RDS_Uniq+$RDS_Multi;
		if ($filterArray[$bowCount] eq "Yes"){
            $RDS_After=$RDS_In-$ALGN_T;
        } else {
            $RDS_After=$RDS_In;
        }
    } elsif ($bowitem =~ m/bow1_([^\\.]+)$/){
        $group = "$1\\t";
        open(IN2, $bowitem);
        my $multimapped;
		my $aligned;
		my $inputCount;
		my $uniqAligned;
		chomp($inputCount = `cat $bowitem | grep '# reads processed:' | awk '{sum+=\\$4} END {print sum}'`);
		chomp($aligned = `cat $bowitem | grep '# reads with at least one reported alignment:' | awk '{sum+=\\$9} END {print sum}'`);
		chomp($uniqAligned = `cat $bowitem | grep '# unique mapped reads:' | awk '{sum+=\\$5} END {print sum}'`);
		## Here we exclude "Number of reads mapped to too many loci" from multimapped reads since in bam file it called as unmapped.
		## Besides, these "too many loci" reads exported as unmapped reads from STAR.
		$RDS_In = int($inputCount);
		$RDS_Multi = int($aligned) -int($uniqAligned);
		if ($RDS_Multi < 0 ){
		    $RDS_Multi = 0;
		}
        $RDS_Uniq = int($uniqAligned);
        $ALGN_T = int($aligned);
		if ($filterArray[$bowCount] eq "Yes"){
            $RDS_After=$RDS_In-$ALGN_T;
        } else {
            $RDS_After=$RDS_In;
        }
    }
    
    print $fh "!{name}\\t$group$RDS_In\\t$RDS_After\\t$RDS_Uniq\\t$RDS_Multi\\t$ALGN_T\\n";
}
close($fh);



'''

}


process Sequential_Mapping_Module_Merge_TSV_Files {

input:
 file tsv from g42_26_outputFileTSV_g42_13.collect()
 val outputFileName from g42_26_name_g42_13.collect()

output:
 file "${name}.tsv"  into g42_13_outputFileTSV_g42_14

script:
name = outputFileName[0]
"""    
awk 'FNR==1 && NR!=1 {  getline; } 1 {print} ' *.tsv > ${name}.tsv
"""
}


process Sequential_Mapping_Module_Sequential_Mapping_Short_Summary {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /sequential_mapping_short_sum.tsv$/) "sequential_mapping/$filename"
	else if (filename =~ /sequential_mapping_detailed_sum.tsv$/) "sequential_mapping/$filename"
}

input:
 file mainSum from g42_13_outputFileTSV_g42_14

output:
 file "sequential_mapping_short_sum.tsv"  into g42_14_outputFileTSV_g_26
 file "sequential_mapping_detailed_sum.tsv"  into g42_14_outputFile

shell:
'''
#!/usr/bin/env perl
use List::Util qw[min max];
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;

my @header;
my %all_rows;
my @seen_cols_short;
my @seen_cols_detailed;
my $ID_header;

chomp(my $contents = `ls *.tsv`);
my @files = split(/[\\n]+/, $contents);
foreach my $file (@files){
        open IN,"$file";
        my $line1 = <IN>;
        chomp($line1);
        ( $ID_header, my @h) = ( split("\\t", $line1) );
        my $totalHeader = $h[1];
        my $afterFilteringHeader = $h[2];
        my $uniqueHeader = $h[3];
        my $multiHeader = $h[4];
        my $mappedHeader = $h[5];
        push(@seen_cols_short, $totalHeader) unless grep{$_ eq $totalHeader} @seen_cols_short; #Total reads Header
        push(@seen_cols_detailed, $totalHeader) unless grep{$_ eq $totalHeader} @seen_cols_detailed; #Total reads Header

        my $n=0;
        while (my $line=<IN>) {
                
                chomp($line);
                my ( $ID, @fields ) = ( split("\\t", $line) ); 
                #SHORT
                push(@seen_cols_short, $fields[0]) unless grep{$_ eq $fields[0]} @seen_cols_short; #mapped element header
                $all_rows{$ID}{$fields[0]} = $fields[5];#Mapped Reads
                #Grep first line $fields[1] as total reads.
                if (!exists $all_rows{$ID}{$totalHeader}){    
                        $all_rows{$ID}{$totalHeader} = $fields[1];
                } 
                $all_rows{$ID}{$afterFilteringHeader} = $fields[2]; #only use last entry
                #DETAILED
                $uniqueHeadEach = "$fields[0] (${uniqueHeader})";
                $multiHeadEach = "$fields[0] (${multiHeader})";
                $mappedHeadEach = "$fields[0] (${mappedHeader})";
                push(@seen_cols_detailed, $mappedHeadEach) unless grep{$_ eq $mappedHeadEach} @seen_cols_detailed;
                push(@seen_cols_detailed, $uniqueHeadEach) unless grep{$_ eq $uniqueHeadEach} @seen_cols_detailed;
                push(@seen_cols_detailed, $multiHeadEach) unless grep{$_ eq $multiHeadEach} @seen_cols_detailed;
                $all_rows{$ID}{$mappedHeadEach} = $fields[5];
                $all_rows{$ID}{$uniqueHeadEach} = $fields[3];
                $all_rows{$ID}{$multiHeadEach} = $fields[4];
    }
    close IN;
    push(@seen_cols_short, $afterFilteringHeader) unless grep{$_ eq $afterFilteringHeader} @seen_cols_short; #After filtering Header
}


#print Dumper \\%all_rows;
#print Dumper \\%seen_cols_short;

printFiles("sequential_mapping_short_sum.tsv",@seen_cols_short,);
printFiles("sequential_mapping_detailed_sum.tsv",@seen_cols_detailed);


sub printFiles {
    my($summary, @cols_to_print) = @_;
    
    open OUT, ">$summary";
    print OUT join ("\\t", $ID_header,@cols_to_print),"\\n";
    foreach my $key ( keys %all_rows ) { 
        print OUT join ("\\t", $key, (map { $all_rows{$key}{$_} // '' } @cols_to_print)),"\\n";
        }
        close OUT;
}

'''


}

g_22_outputFileTSV_g_26= g_22_outputFileTSV_g_26.ifEmpty([""]) 
g42_14_outputFileTSV_g_26= g42_14_outputFileTSV_g_26.ifEmpty([""]) 
g_28_outputFileTSV_g_26= g_28_outputFileTSV_g_26.ifEmpty([""]) 
g0_11_outputFileTSV_g_26= g0_11_outputFileTSV_g_26.ifEmpty([""]) 
g0_21_outputFileTSV_g_26= g0_21_outputFileTSV_g_26.ifEmpty([""]) 
g0_16_outputFileTSV_g_26= g0_16_outputFileTSV_g_26.ifEmpty([""]) 

//* autofill
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 30
    $CPU  = 1
    $MEMORY = 10
    $QUEUE = "short"
}
//* platform
//* autofill

process Overall_Summary {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /overall_summary.tsv$/) "Summary/$filename"
}

input:
 file starSum from g_22_outputFileTSV_g_26
 file sequentialSum from g42_14_outputFileTSV_g_26
 file rsemSum from g_28_outputFileTSV_g_26
 file adapterSum from g0_11_outputFileTSV_g_26
 file trimmerSum from g0_21_outputFileTSV_g_26
 file qualitySum from g0_16_outputFileTSV_g_26

output:
 file "overall_summary.tsv"  into g_26_outputFileTSV

shell:
'''
#!/usr/bin/env perl
use List::Util qw[min max];
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;

my @header;
my %all_rows;
my @seen_cols;
my $ID_header;

chomp(my $contents = `ls *.tsv`);
my @rawFiles = split(/[\\n]+/, $contents);
my @files = ();
# order must be in this order for chipseq pipeline: bowtie->dedup
# rsem bam pipeline: dedup->rsem, star->dedup
# riboseq ncRNA_removal->star
my @order = ("adapter_removal","trimmer","quality","extractUMI","sequential_mapping","ncRNA_removal","bowtie","star","hisat2","tophat2", "dedup","rsem","kallisto");
for ( my $k = 0 ; $k <= $#order ; $k++ ) {
    for ( my $i = 0 ; $i <= $#rawFiles ; $i++ ) {
        if ( $rawFiles[$i] =~ /$order[$k]/ ) {
            push @files, $rawFiles[$i];
        }
    }
}

print Dumper \\@files;
##add rest of the files
for ( my $i = 0 ; $i <= $#rawFiles ; $i++ ) {
    push(@files, $rawFiles[$i]) unless grep{$_ == $rawFiles[$i]} @files;
}
print Dumper \\@files;

##Merge each file according to array order

foreach my $file (@files){
        open IN,"$file";
        my $line1 = <IN>;
        chomp($line1);
        ( $ID_header, my @header) = ( split("\\t", $line1) );
        push @seen_cols, @header;

        while (my $line=<IN>) {
        chomp($line);
        my ( $ID, @fields ) = ( split("\\t", $line) ); 
        my %this_row;
        @this_row{@header} = @fields;

        #print Dumper \\%this_row;

        foreach my $column (@header) {
            if (! exists $all_rows{$ID}{$column}) {
                $all_rows{$ID}{$column} = $this_row{$column}; 
            }
        }   
    }
    close IN;
}

#print for debugging
#print Dumper \\%all_rows;
#print Dumper \\%seen_cols;

#grab list of column headings we've seen, and order them. 
my @cols_to_print = uniq(@seen_cols);
my $summary = "overall_summary.tsv";
open OUT, ">$summary";
print OUT join ("\\t", $ID_header,@cols_to_print),"\\n";
foreach my $key ( keys %all_rows ) { 
    #map iterates all the columns, and gives the value or an empty string. if it's undefined. (prevents errors)
    print OUT join ("\\t", $key, (map { $all_rows{$key}{$_} // '' } @cols_to_print)),"\\n";
}
close OUT;

sub uniq {
    my %seen;
    grep ! $seen{$_}++, @_;
}

'''


}

params.run_FastQC =  "no"  //* @dropdown @options:"yes","no" @description:"FastQC provides quality control checks on raw sequence data."



process Adapter_Trimmer_Quality_Module_FastQC {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*.(html|zip)$/) "FastQC/$filename"
}

input:
 val mate from g_1_mate_g0_3
 set val(name), file(reads) from g_2_reads_g0_3

output:
 file '*.{html,zip}'  into g0_3_FastQCout

errorStrategy 'retry'
maxRetries 3

when:
(params.run_FastQC && (params.run_FastQC == "yes"))

script:
nameAll = reads.toString()
if (nameAll.contains('.gz')) {
    file =  nameAll - '.gz' - '.gz'
    runGzip = "ls *.gz | xargs -i echo gzip -df {} | sh"
} else {
    file =  nameAll 
    runGzip = ''
}
"""
${runGzip}
fastqc ${file} 
"""
}


workflow.onComplete {
println "##Pipeline execution summary##"
println "---------------------------"
println "##Completed at: $workflow.complete"
println "##Duration: ${workflow.duration}"
println "##Success: ${workflow.success ? 'OK' : 'failed' }"
println "##Exit status: ${workflow.exitStatus}"
}
