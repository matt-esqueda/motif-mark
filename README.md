# motif-mark

Motif Mark with OOP

Inputs: fasta_file -f, motif_file -m, oneline_fasta_output -o
Outputs: PNG file 

This python script uses OOP code to create an image of motifs on gene sequences. 
The script takes a FASTA file and motifs file and outputs a single image per FASTA file with the same prefix.
It is capable of handling ambiguous nucleotides and all features (motifs, introns, exons) are to scale.

Run in this environment: conda activate my_pycairo
