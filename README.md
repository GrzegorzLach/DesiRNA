# DesiRNA
RNA sequence design software

Requirements: python, pycosat, RNAlib (from the ViennaRNA package)

usage: python DesiRNA -f [input filename] [options]

  -s SET, --set SET                                         Name of a file containing set of the preferred motifs.
  
  -m NUMBER_OF_MUTATIONS, --mutations NUMBER_OF_MUTATIONS   Number of mutation cycles.
  
  -v {0,1}, --verbose {0,1}                                 Displaying solution from SATsolver.
  
  -p PROB, --probability PROB                               Probability of mutating the whole fragment from library of motifs.
  
  -c {0,1}, --complex {0,1}                                 Include RNA-RNA complex formation in  DesiRNA scoring function.
  
  -n NUMBER_SEQ, --number_of_sequences NUMBER_SEQ           Number of best sequences to return.


with the example.fas structured as:

\> structure 1

[RNA secondary structure in dot-bracket notation]

\> structure 2

[RNA secondary structure in dot-bracket notation]

...

[RNA squence constraints in IUPAC notation]

prevented [1st prevented RNA subsequence in IUPAC notation]

prevented [2nd prevented RNA subsequence in IUPAC notation]

...
