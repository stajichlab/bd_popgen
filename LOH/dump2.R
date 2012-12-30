library(IRanges)
load(file="LOH_start_stops.RData")
ls()
##### check out contig 2 - this is in a "list" format
## note that each item in the list is an IRanges object
strains <- names(contig.1)
##check out the starts/ends of LOHs for a strain
# note the double brackets pulls an item out of a list, and the output is no longer a list

chroms <- array(1:10,c(1,10))
chroms[1] <- contig.1
chroms
for (strain in 1:3) {
  chr = 1
  while( chr < 11 ) {
  chrom <- chroms[chr]
  numbreaks <- length(start(chrom[[strain]]))
  if( numbreaks > 0 ) {
    for ( brk in 1:numbreaks ) {
      print(c(strains[strain],chr,brk,start(chrom[[strain]][brk]),
              end(chrom[[strain]][brk])
              ))
    }
  chr = chr + 1	
  }
}
}
