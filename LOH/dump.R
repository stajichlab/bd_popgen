library(IRanges)
load(file="LOH_start_stops.RData")
ls()
##### check out contig 2 - this is in a "list" format
## note that each item in the list is an IRanges object
strains <- names(contig.1)
##check out the starts/ends of LOHs for a strain
# note the double brackets pulls an item out of a list, and the output is no longer a list

contig.5[[47]]

chroms <- c(contig.1,contig.2,contig.3,contig.4,contig.5,contig.6,contig.7,
            contig.8,contig.9,contig.10)
for (strain in 1:4) {
  chr = 1;
  chrom <- contig.1
  numbreaks <- length(start(chrom[[strain]]))
  if( numbreaks > 0 ) {
    for ( brk in 1:numbreaks ) {
      write(c(strains[strain],chr,brk,start(chrom[[strain]][brk]),
		   end(chrom[[strain]][brk])
              ), sep="\t")
    }
  }
  chr = chr + 1
  chrom <- contig.2
  numbreaks <- length(start(chrom[[strain]]))
  for ( brk in 1:numbreaks ) {
    print(c(strains[strain],chr,brk,start(chrom[[strain]][brk]),
            end(chrom[[strain]][brk])
            ))
  }

  chr = chr + 1
  chrom <- contig.3
  numbreaks <- length(start(chrom[[strain]]))
  for ( brk in 1:numbreaks ) {
    print(c(strains[strain],chr,brk,start(chrom[[strain]][brk]),
            end(chrom[[strain]][brk])
            ))
  }

  chr = chr + 1
  chrom <- contig.4
  numbreaks <- length(start(chrom[[strain]]))
  for ( brk in 1:numbreaks ) {
    print(c(strains[strain],chr,brk,start(chrom[[strain]][brk]),
            end(chrom[[strain]][brk])
            ))
  }

  chr = chr + 1
  chrom <- contig.5
  numbreaks <- length(start(chrom[[strain]]))
  for ( brk in 1:numbreaks ) {
    print(c(strains[strain],chr,brk,start(chrom[[strain]][brk]),
            end(chrom[[strain]][brk])
            ))
  }

  chr = chr + 1
  chrom <- contig.6
  numbreaks <- length(start(chrom[[strain]]))
  for ( brk in 1:numbreaks ) {
    print(c(strains[strain],chr,brk,start(chrom[[strain]][brk]),
            end(chrom[[strain]][brk])
            ))
  }

  chr = chr + 1
  chrom <- contig.7
  numbreaks <- length(start(chrom[[strain]]))
  for ( brk in 1:numbreaks ) {
    print(c(strains[strain],chr,brk,start(chrom[[strain]][brk]),
            end(chrom[[strain]][brk])
            ))
  }

  chr = chr + 1
  chrom <- contig.8
  numbreaks <- length(start(chrom[[strain]]))
  for ( brk in 1:numbreaks ) {
    print(c(strains[strain],chr,brk,start(chrom[[strain]][brk]),
            end(chrom[[strain]][brk])
            ))
  }

  chr = chr + 1
  chrom <- contig.9
  numbreaks <- length(start(chrom[[strain]]))
  if( numbreaks > 0 ) {
    for ( brk in 1:numbreaks ) {
      print(c(strains[strain],chr,brk,start(chrom[[strain]][brk]),
              end(chrom[[strain]][brk])
              ))
    }
  }


  chr = chr + 1
  chrom <- contig.10
  print(names(chrom))
  numbreaks <- length(start(chrom[[strain]]))
  if( numbreaks > 0 ) {
    for ( brk in 1:numbreaks ) {
      print(c(strains[strain],chr,brk,start(chrom[[strain]][brk]),
              end(chrom[[strain]][brk])
              ))
    }
  }
}
