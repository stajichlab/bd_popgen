use strict;
use warnings;
use Getopt::Long;
use List::Util qw(sum min max);

my %d;
my %strains;
my $Ovlrange = 10_000;
my $minsize  = 1000;

GetOptions(
    'r|range:i' => \$Ovlrange,
    'm|min:i'   => \$minsize,
    );

while(<>){
    my ($strain,$chrom,$start,$end) = split;
    if( abs($end-$start) < $minsize ) {
	warn("$strain $chrom $start $end -- too short\n");
	next;
    }
    $strains{$strain}++;
    push @{$d{$chrom}}, [$strain,$start,$end];    
}
my %blocks;
while( my ($chrom,$data) = each %d ) {
    for my $bl ( &cluster($data) ) {
	$blocks{join("-",$chrom,$bl->[0],$bl->[1])} = $bl->[2];
    }
}
my @blnames =  sort keys %blocks;
my @strains = sort keys %strains;
print "#NEXUS\n";
print "[DATE ".localtime,"]\n";
print "[SITES ",join("\t", qw(STRAIN), @blnames), "]\n";
print "begin data;\n";
printf "dimensions ntax=%d nchar=%d;\n",
    scalar @strains, scalar @blnames;
print "format datatype=standard gap=- interleave=no;\n";
print "matrix\n";
for my $strain ( @strains ) {
    print $strain, "\t", join("",map { $blocks{$_}->{$strain} || 0 } @blnames), "\n";
}
print ";\n";
print "end;\n";
sub cluster {
    my $d = shift;
    my @clusters;
    for my $row ( sort { $a->[1] <=> $b->[1] } @$d ) {
	my ($str,$start,$end) = @$row;
	if( ! @clusters ) {
	    push @clusters, [$start, $end, { $str => 1}];
	} else {
	    my $cl;
	    for my $c ( @clusters ) {
		if( &overlaps($start,$end, $c->[0],$c->[1],$Ovlrange) ) {
		    $c->[2]->{$str}++;
		    $c->[0] = min($start,$c->[0]);
		    $c->[1] = max($start,$c->[1]);
		    $cl = 1;
		    last;
		}
	    }
	    unless($cl) {
		push @clusters, [$start, $end, { $str => 1}];
	    }
	}
    }
    my $continue = 1;
    while( $continue ) {
	$continue = 0;
	for(my $cl1=0;$cl1 < scalar @clusters; $cl1++ ) {	    
	    for(my $cl2 = $cl1+1; $cl2 < scalar @clusters; $cl2++ ) {
		if( &overlaps($clusters[$cl1]->[0],
			      $clusters[$cl1]->[1],
			      $clusters[$cl2]->[0],
			      $clusters[$cl2]->[1],
			      $Ovlrange) ) {
		    for my $str ( keys %{$clusters[$cl1]->[2]} ) {
			warn("str is $str\n");
			$clusters[$cl1]->[2]->{$str}++;
		    }
		    $clusters[$cl1]->[0] = min($clusters[$cl1]->[0],
					       $clusters[$cl2]->[0]);

		    $clusters[$cl1]->[1] = max($clusters[$cl1]->[1],
					       $clusters[$cl2]->[2]);
		    $continue = 1;
		}			      
	    }
	}
    }
    @clusters;
}

sub overlaps {
    my ($left_b,$left_e,
	$right_b,$right_e,
	$wiggle) = @_;
    my $rc = 0;
    if( abs($left_b - $right_b) < $wiggle &&
	abs($left_b - $right_b) < $wiggle ) {
	$rc = 1;
    }
    $rc;
}

sub mean {
    my @vals = @_;
    sum(@vals) / scalar @vals;
}
