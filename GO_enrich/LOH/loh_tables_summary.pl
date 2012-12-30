#!/usr/bin/perl -w
use strict;
use Bio::Range;
use List::Util qw(min max sum);

#change this to change the block edges overlap that is allowed to merge two blocks
my $bigger = 10000;
#my $bigger = 55000;
my $dir = 'LOH_tables';

my %GPL = map { $_ => 1 } qw(MCT8 LBAbercrom NCRC106979 JEL423 JEL408 JEL429 JEL427 JEL310 EV001 CLFT024 JEL433 JEL275 JEL274 CLFT023 CLFT021 L2203 Au.Lcaerulae98 JEL271 JEL270 CJB7 CJB5.2 JEL627 MLA1 JEL267 CLFT026 CJB4 UKTvB MexMkt SRS812 JEL289 TSTS75 JEL238 JEL359 JEL261 IA042 PNP08489 X0711);

print join("\t", qw(CONTIG START END GPL_STRAINS STRAIN_COUNT)), "\n";
opendir(DIR,$dir) || die $!;
for my $file ( readdir(DIR) ) {
    next unless $file =~ /(\S+)\.txt/;
    my $stem = $1;
    open(my $fh => "$dir/$file") || die $!;
    my ($strain);
    my @data;
    while(<$fh>) {
	if(/^>\s+(\S+)/) {
	    $strain = $1;
#	    if( ! $GPL{$strain} ) {
#		warn($strain,"\n");
#	    }
	} elsif(/^start/) {
	    if( ! defined $strain ) {
		die("did not parse header properly\n");
	    }
	} else {
	    next unless $GPL{$strain};
	    my ($start,$end,$width) = split;
	    my $range = Bio::Range->new(-start => $start,
					-end   => $end);
	    my $skip = 0;
	    for my $r ( @data ) {
		if( $r->[0]->equals($range) ) {
#		    $r->[0] = $r->[0]->union($range);
		    $r->[1]->{$strain} = 1;
		    $skip = 1;
		    last;
		}
	    }
	    unless( $skip ) {
		push @data, [$range, { $strain => 1 }];
	    }
	}
    }
    # merge
    my $change = 1;
    while( $change ) {
	$change = 0;
	@data = grep { defined $_->[0] } @data;
	for( my $i =0; $i < scalar @data; $i++ ){
	    my $r_i = $data[$i];
	    next if ! defined $r_i->[0];
	   # warn("begin i ", $r_i->[0]->toString, " ", 
		# join(",",sort keys %{$r_i->[1]}), "\n");
	    
	    for( my $j = $i+1; $j < scalar @data; $j++ ) {
		next unless defined $data[$j];
		my $r_j = $data[$j];
		next if ! defined $r_j->[0];
	#	warn("begin j ", $r_j->[0]->toString, " ", 
	#	     join(",",sort keys %{$r_j->[1]}), "\n");
	    

		if( abs($r_i->[0]->start - $r_j->[0]->start) < $bigger &&
		    abs($r_i->[0]->end - $r_j->[0]->end) < $bigger ) {
		    # pretty close in space
		    warn("merging ", $r_i->[0]->toString, " ", 
			 join(",",sort keys %{$r_i->[1]}), " with ", 
			 $r_j->[0]->toString," ", 
			 join(",",sort keys %{$r_j->[1]}),"\n");
		    
		    $r_i->[0]->start( min( $r_i->[0]->start,
					   $r_j->[0]->start));

		    $r_i->[0]->end( max( $r_i->[0]->end,
					 $r_j->[0]->end));		    
		    for( keys %{$r_j->[1]} ) { $r_i->[1]->{$_} = 1; }
		    warn("after merging ", $r_i->[0]->toString, " ", 
			 join(",",sort keys %{$r_i->[1]}), " with ", 
			 $r_j->[0]->toString," ", 
			 join(",",sort keys %{$r_j->[1]}),"\n");
		    $data[$j] = [];
		    $data[$i] = $r_i;
		    $change = 1;
		}	       
	    }
#	    warn("\n");
	}
    }
    
    my %counts;
    @data = grep { defined $_->[0] } @data;
    for my $r ( sort { $a->[0]->start <=> $b->[0]->start ||
		       $a->[0]->end <=> $b->[0]->end 
		} @data ) {
	print join("\t", $stem, $r->[0]->start, $r->[0]->end, 
		   join(",",sort keys %{$r->[1]}),
		   scalar keys %{$r->[1]},
	    ),"\n";
    }
#    last;
}
