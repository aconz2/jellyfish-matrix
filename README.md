# jellyfish-matrix #
The purpose of this tool is to create a binary presence/absence matrix across many samples.
The input is a list of .jf jellyfish files which must all have the same k-mers and in the same order. See the aconz2/jellyfish-mirror tool for this purpose.
In a single sample, presence is indicated by some lower bound, defaulting to 1. Typically you may set this higher so that a k-mer is only considered present if it appears 3 times for example.
In the output, a k-mer is kept if it is considered present in at least n samples, but no more than m samples. This is to filter k-mers which are uniformative either because they are not present in enough or too many of the samples.
The output will use n bits in the counter field for n input files.
The counter field for a particular k-mer is the binary concatenation of each sample. So 101 = 5 would indicate the k-mer is present in the 1st and 3rd samples. Ordering is maintained left to right as per the input order.
