while(<>) {
    my @row = split;
    print if $row[-1] >= 25;
}
