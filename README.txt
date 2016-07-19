# StructureFold
Multiple RNA secondary structures prediction from high throughput RNA structure profiling data

Authors: Yin Tang
         David C. Tack

StructureFold is a series of software packages that automates the process of predicting RNA secondary structure for a transcript or an entire transcriptome, with or without the inclusion of constraints on the structure(s) provided by wet bench experimentation. The process consists of mapping the raw reads of RNA structural data on every transcript in the dataset to the transcriptome, getting RT stop counts on each nucleotide, calculating structural reactivities on the nucleotides, and predicting the RNA structures.


Please cite:
Tang Y, Bouvier E, Kwok CK, Ding Y, Nekrutenko A, Bevilacqua PC, Assmann SM. StructureFold: genome-wide RNA secondary structure mapping and reconstruction in vivo. Bioinformatics. 2015;31:2668-75.


#Package contents
Main Scripts
        The main programs of the package

Accessory Scripts
        Fast ways to organize your data or convert it into a useful format

        generate_PPV_file.py
                Requires RNAStructure package. Iterates through two directories of .ct files and computes all PPV values.

        generate_coverages.py
                Creates lists of coverage overlap between all coverage files in a directory.
                Overlap lists are useful for other tools.
                
