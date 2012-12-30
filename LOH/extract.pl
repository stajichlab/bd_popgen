#!/usr/bin/perl -w
use strict;

use Statistics::R;
  
# Create a communication bridge with R and start R
my $R = Statistics::R->new();

my $out = $R->run(<<EOF
library(IRanges)
load(file="LOH_start_stops.RData")
ls()
EOF
    );

my $names = $R->run( "names<-names(contig.1)" );
$names= $R->get('names');
my @namelist = @$names;
for my $n ( 1..49 ) {
    for my $chr ( 1..10 ) {
	$R->set( "strain",$n );
	my $starts = $R->run( "start<-start(contig.$chr\[\[strain]])" );
	my $ends = $R->run( "end<-end(contig.$chr\[\[strain]])" );
	$starts = $R->get( "start" );
	$ends = $R->get( "end" );
	my @pairs;
	
	if( ref($starts) =~ /ARRAY/ ) {
	    push @pairs, @$starts;
	} elsif( $starts ne "integer(0)" ) {
	    push @pairs, $starts;
	}
	if( ref($ends) =~ /ARRAY/ ) {
	    my $i = 0;
	    for my $end ( @$ends ) {
		$pairs[$i] = [$pairs[$i], $end];
		$i++;
	    }
	} elsif( $ends ne "integer(0)" ) {
	    $pairs[0] = [$pairs[0], $ends];
	}
	for my $pairs ( @pairs ) {
	    print join("\t", $namelist[$n-1],
		       "Chr$chr",
		       @$pairs),"\n";
	    
	}
    }
}
