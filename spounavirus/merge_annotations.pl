#__perl__
#
# merge the annotations from Evelien's GFF files and the headers from the GenBank files that suck
#


use strict;
use Bio::SeqIO;

my $genbankF = shift || die "$0 <genbank file> <genbank from gff3 file> <output file>";
my $gffF = shift || die "$0 <genbank file> <genbank from gff3 file> <output file>";
my $newF = shift || die "$0 <genbank file> <genbank from gff3 file> <output file>";

my $sio=Bio::SeqIO->new(-file=>$genbankF, -format=>'genbank');
my $sout = Bio::SeqIO->new(-file=>">$newF", -format=>'genbank');
my @toremove;
my $seq=$sio->next_seq;
foreach my $feature ($seq->top_SeqFeatures()) {
	my $ftype;
	eval {$ftype=$feature->primary()};
	next unless ($ftype eq "CDS");
	push @toremove, $feature;
}
$seq->remove_SeqFeatures(@toremove);

my $ac = $seq->annotation;
my $comment = Bio::Annotation::Comment->new;
$comment->text('THIS GENOME WAS REANNOTATED USING PROKKA BY EVELIEN ADRIAENSSENS AND ROB EDWARDS');
$ac->add_Annotation('comment', $comment);
$seq->annotation($ac);

my $sin = Bio::SeqIO->new(-file=>$gffF, -format=>"genbank");
my $rseq = $sin->next_seq;
foreach my $feature ($rseq->top_SeqFeatures()) {
	my $ftype; 
	eval {$ftype=$feature->primary()};

	my $trans = $feature->seq->translate->seq;
	$feature->add_tag_value('translation', $trans);
	my $lt = join("_", $seq->display_id, $feature->get_tag_values('locus_tag'));

	$feature->add_tag_value('protein_id', $lt);

	$seq->add_SeqFeature($feature);
}
$sout->write_seq($seq);

