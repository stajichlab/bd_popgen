use strict;
use warnings;
my %seen;
my %trans2gene;
open(my $fh => 'BD_trans_to_gene.tab') || die("cannot open lookup file");
while (<$fh>) {
  my ($trans,$gene) = split;
  $trans2gene{$trans} = $gene;
}
while (<>) {
  chomp;
  my ($trans,@rest) = split(/\t/,$_);
  if ( ! $trans2gene{$trans} )  {
    die("cannot find $trans\n");
  }
  while ($rest[12] = /(GO:\d+)/g) {
    next if( $seen{$trans."-".$1}++);
    print join("\t",$1, 'IEA',$trans2gene{$trans}),"\n";
  }
}
