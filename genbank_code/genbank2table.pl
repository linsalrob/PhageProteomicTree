#!/usr/bin/perl -w

# get proteins out of a genbank file

use Bio::SeqIO;
use strict;

my $usage=<<EOF;
$0 [OPTIONS] <list of genbankfiles>

OPTIONS:
	-v   (for verbose output)
	-t   (to add the sequence accession to all the ids printed)
	-n   (to add the file name to all the ids printed)
	-d <filename>  file for DNA sequence of the genes
	-g <filename>  file for DNA sequence of the genome
	-p <filename>  file for the protein sequences
	-f <filename>  file for the functions
	-i <filename>  file for the sequence accessions and names

You must specify one of -d, -g, -p, -f or -i
EOF

die $usage unless ($ARGV[0]);

my @infiles;
my $verbose=0;

my $dnaF;
my $protF;
my $funcF;
my $genomeF;
my $adddesc;
my $addfilename;
my $idF;

while (@ARGV) {
	my $f = shift @ARGV;
	if ($f eq "-d") {$dnaF = shift @ARGV}
	elsif ($f eq "-p") {$protF = shift @ARGV}
	elsif ($f eq "-f") {$funcF = shift @ARGV}
	elsif ($f eq "-g") {$genomeF = shift @ARGV}
	elsif ($f eq "-i") {$idF = shift @ARGV}
	elsif ($f eq "-v") {$verbose = 1}
	elsif ($f eq "-t") {$adddesc = 1}
	elsif ($f eq "-n") {$addfilename = 1}
	else {push @infiles, $f}
}

unless ($idF || $protF || $funcF || $genomeF || $dnaF) {die "Please specify one of -f, -g, -p, -d.\n$usage\n"}

if ($protF && -e $protF) {die "$protF already exists. Not overwriting\n"}
if ($funcF && -e $funcF) {die "$funcF already exists. Not overwriting\n"}
if ($dnaF && -e $dnaF) {die "$dnaF already exists. Not overwriting\n"}
if ($genomeF && -e $genomeF) {die "$genomeF already exists. Not overwriting\n"}

if (defined $protF) {open(FA, ">$protF") || die "Can't open fasta $protF"}
if (defined $funcF) {open(PR, ">$funcF") || die "Can't open function $funcF"}
if (defined $dnaF) {open(SQ, ">$dnaF") || die "can't open sequences $dnaF"}
if (defined $genomeF) {open(GF, ">$genomeF") || die "can't open sequences $genomeF"}
if (defined $idF) {open(IDF, ">$idF") || die "Can't open ID file $idF"}

my $seedid=&seedid;
if ($seedid) {open(CO, ">contigs") || die "Can't open contigs"}

my $c;
foreach my $file (@infiles)
{
	my $desctag = "";
	my $sio=Bio::SeqIO->new(-file=>$file, -format=>'genbank');
	while (my $seq=$sio->next_seq) {
		if ($adddesc) {$desctag .= " [" . $seq->accession. "]"}
		my $shortfilename = $file;
		$shortfilename =~ s#^.*/##;
		$shortfilename =~ s#.gbk##;
		$shortfilename =~ s#.gb##;

		if ($addfilename) {
			$desctag .= " [" . $shortfilename . "]";
		}

		if (defined $genomeF) {print GF ">", $seq->display_name, "$desctag\n", $seq->seq, "\n"}
		if (defined $idF) {
			my $shortname = &shortname($seq->desc);
			print IDF join("\t", $seq->accession, $seq->desc, $shortname, $shortfilename), "\n";
		}
		my $seqname=$seq->display_name;
		my $source = "";
		#my $source = $seq->source;
		if ($verbose) {print STDERR "Parsing $seqname\n"}
		if ($seedid->{$seqname}) {print CO join("\t", $seedid->{$seqname}, $seqname)}
		my $organism; my $taxid;

		foreach my $feature ($seq->top_SeqFeatures()) {
			#my $ftype;
			#eval {$ftype=$feature->primary()};

			#next unless ($ftype eq "CDS");
			if ($feature->has_tag('organism')) {
				$organism = join(" ", $feature->each_tag_value('organism'));
				my $dbx = join(" ", $feature->each_tag_value('db_xref'));
				$dbx =~ m/taxon:(\d+)/;
				$taxid = $1;
			}

			$c++;
			my $id; # what we will call the sequence
			my ($trans, $gi, $geneid, $prod, $locus, $np);

			eval {$locus = join " ", $feature->each_tag_value("locus_tag")};
			eval {$prod  = join " ", $feature->each_tag_value("product")};
			eval {$trans = join " ", $feature->each_tag_value("translation")};
			unless ($trans) {$trans = $feature->seq->translate->seq()}
			eval {$np = join " ", $feature->each_tag_value("protein_id")};
			if (!$np) {
				my $fig = "";
				eval {$fig = join(" ", $feature->each_tag_value("db_xref"))};
				if ($fig) {
					$fig =~ m/SEED\:(fig\|\d+\.\d+\....\.\d+)/;
					$np = $1;
				}
			}
			if ($trans && !$np) {
				if ($locus) {
					if ($verbose) {print STDERR "Using LOCUS as NP\n"}
					$np = $locus;
				} else {
					if ($verbose) {print STDERR "No NP for $trans. Skipped\n"}
					next;
				}
			}
			elsif (!$trans && $np) {
				if ($verbose) {print STDERR "No translation for $np. Skipped\n"}
				next;
			}
			next unless ($trans && $np);

			eval {
				foreach my $xr ($feature->each_tag_value("db_xref")) 
				{
					($xr =~ /GI/) ? ($gi = $xr) : 1;
					($xr =~ /GeneID/) ? ($geneid = $xr) : 1;
				}
			};


			my $end = $feature->end;
			my $start = $feature->start;

			my $oids="";
			($locus)   && ($oids.="locus:$locus;");
			($geneid)  && ($oids.="$geneid;");
			($gi)      && ($oids.="$gi;");
			$oids =~ s/\;$//;

			unless ($prod)  {
				if ($verbose) {print STDERR "No product for $np\n"}
				$prod="hypothetical protein";
			}
			if (defined $dnaF) {print SQ ">$np$desctag\n", $feature->seq()->seq(), "\n"}
			if (defined $protF) {print FA ">$np$desctag\n$trans\n"}
			if (defined $funcF) {print PR "$np$desctag\t$start\t$end\t$prod\t$oids\t$seqname\t$organism\t$taxid\t", $seq->display_name, "\n"}
		}
	}
}


sub seedid {
	return unless (-e "../NCs.txt");
	open(IN, "../NCs.txt");
	my %seedid;
	while (<IN>) 
	{
		chomp;
		my ($id, $name, @acc)=split /\t/;
		map {$seedid{$_}=$id} @acc;
	}
	return \%seedid;
}


sub shortname {
	my $n=shift;
	$n .= " ";
	$n =~ s/, complete genome\.//i;
	$n =~ s/ complete genome\.//i;
	$n =~ s/\s+genome\s+/ /; 
	$n =~ s/^unclassified\s+//i;
	$n =~ s/unclassified\.$//i;
	$n =~ s/^\S+viridae\s+//;
	$n =~ s/^\S+virus\s+//;
	$n =~ s/^\S+virales\s+//;
	$n =~ s/^\S+virinae\s+//;
	$n =~ s/^\S+\-like\s+//;
	$n =~ s/^Viruses\s+//i; $n =~ s/\s+Viruses\.\s+//i;
	$n =~ s/\s+prophage\s+/ /i;
	$n =~ s/\s+bacteriophage\s+/ /i;
	$n =~ s/\s+phage\s+/ /i;
	$n =~ s/\s+temperate\s+/ /i;
	$n =~ s/\s+complex\s+/ /i;
	$n =~ s/\s+filamentous virus\s+/ fv /;
	$n =~ s/\s+spindle-shaped virus\s+/ ssv /;
	$n =~ s/^\s*phages\s+//;
	$n =~ s/\bcontig\d+\b//;
	$n =~ s/\bsequence\b//;
	$n =~ s/\bDNA\b//;
	$n =~ s/\s+phage\s+/ /;
	$n =~ s/\s+sensu lato\s+/ /;
	$n =~ s/\s+virus\s+/ /;
	$n =~ s/[\,\.\s]+/ /g;
	$n =~ s/\s+$//;
	$n =~ s/^(\S{3})\S*\s+/$1\. /;
	return $n;
}

