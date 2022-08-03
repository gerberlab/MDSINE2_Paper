#!/usr/bin/perl -w
use strict;

#################
#Replace RDP sequenceID with species name
################

if(!$ARGV[0] or !$ARGV[1]){
    print "Please type the tree XML file and reference species file:\nGut_Mouse_Initial-rd_v4_HITI_tog.xml\nGut_Mouse_16S_seq_info.csv\n";
    exit;
}
my $tree_infile=$ARGV[0];
my $taxa_infile=$ARGV[1];

my ($name)=($ARGV[0]=~/(.*)\.xml/);
my $outfile=$name."_speciesName.xml";
open(OUT,">$outfile") or die "$outfile:$!\n";

#Key:seqID
#value:speciesName
my %hash=();
get_taxa($taxa_infile,\%hash);

#Replace seqID with species name
replace_seqID($tree_infile,\%hash);
close(OUT);
exit;

#################
sub get_taxa{
    my ($file,$H)=@_;
    open(FH,"$file") or die "$file:$!\n";
    my $count=0;
    while(my $line=<FH>){
	$count++;
	chomp($line);
	if($count==1){
	    next;
	}
	my @array=split(',',$line);
	my $id=$array[0];
	my $name=$array[-1];
	$name=~s/ /_/g;
	$$H{$id}=$name;
    }
    close(FH);
    return;
}

sub replace_seqID{
    my ($file,$H)=@_;
    open(FH,"$file") or die "$file:$!\n";
    my $count=0;
    while(my $line=<FH>){
	$count++;
	if($count==1){
	    print OUT "$line";
	}else{
	    foreach my $id (sort keys %{$H}){
		$line=~s/$id/$$H{$id}/;
	    }
	    print OUT "$line";
	}
    }
    close(FH);
    return;
}
	

