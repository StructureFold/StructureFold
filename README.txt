# StructureFold
Multiple RNA secondary structures prediction from high throughput RNA structure profiling data

Authors: Yin Tang
         David C. Tack

StructureFold is a series of software packages that automates the process of predicting RNA secondary structure for a transcript or an entire transcriptome, with or without the inclusion of constraints on the structure(s) provided by wet bench experimentation. The process consists of mapping the raw reads of RNA structural data on every transcript in the dataset to the transcriptome, getting RT stop counts on each nucleotide, calculating structural reactivities on the nucleotides, and predicting the RNA structures.


Please cite:
Tang Y, Bouvier E, Kwok CK, Ding Y, Nekrutenko A, Bevilacqua PC, Assmann SM. StructureFold: genome-wide RNA secondary structure mapping and reconstruction in vivo. Bioinformatics. 2015;31:2668-75.


#Package contents
         All scripts contain a help function. Run script with -h for all options.

Main Scripts
        The main programs of the package
        
        reactivitiy_calculator.py
               Derives structural reactivity from RT stop count files.
               Requires RT stop files from - and + DMS treated samples
               IN: <rtsc>, transcript.fasta, transcript list
               OUT: <.react>, normalization_scales
               
         batch_folder.py
               Uses a restraint file, a list of transcripts, and their sequences to predict secondary structure.
               Requires either RNAStructure or ViennaPackage
               IN: <.react>, transcript.fasta, transcript list
               OUT: <.ct>, <.ps>
               

Accessory Scripts
        Fast ways to organize your data or convert it into a useful format.

        generate_PPV_file.py
                Requires RNAStructure package. 
                Iterates through two directories of .ct files and computes all PPV values.
                Only transcripts appearing in both directories are compared.
                IN: <.ct>
                OUT: <.csv>
         
         reactivity_stats.py
                Uses all .rtsc files in a directory to generate basic statsistics.
                Requires an overlap list as a key file; only transcrips on this list will be processed
                Options to remove last n bp of transcripts, enforce minimum length.
                IN: directory of <.rtsc>, overlap
                OUT: <.csv>

         generate_coverages.py
                Creates lists of coverage overlap between all coverage files in a directory.
                Overlap lists are useful for other tools.
                IN: directory of coverages
                OUT: directory of overlaps of coverage
                
         get_coverage.py
                Calculates transcript coverage using .react files and the .fasta of the transcripts
                IN: <.rtsc>, transcript .fastq
                OUT: coverages
                
         get_specificity.py
                Calculates the stop specificity using .react files and the .fasta of the transcripts
                IN: <.rtsc>, transcript .fastq
                OUT: specificities
         
         get_replicate_correlation.py
                 Creates a csv to compare the correlation of RT stop counts.
                 IN: <.rtsc>
                 OUT: <.csv>
         
         get_abundance.py
                  Creates a csvs of abundances.
                  
         reset_n_bp.py
                  sets the last n bases of all .react files in the directory to NA.
                  used when the molecular technique cannot resolve the reactivity at the last n bp.
                  
         process_sams.py
                  Removes reads where there are more than 3 mismatches, a mismatch on first base,
                  or are in the wrong orientation. Requires SamTools. Workson a directory.
                  IN: <.sam>
                  OUT: <.sam>
         
         combine_RT_stop.py
                  Combines two .rtsc files. 
                  IN: <.rtsc>
                  OUT: <.rtsc>
