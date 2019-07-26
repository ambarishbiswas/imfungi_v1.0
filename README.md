#### imFungi version:V1.0 
identify metagenomic Fungi (imFungi): A novel fungi prediction tool for (meta)genomes.


#### imFungi syntax:

     docker run --rm -it --name imfungi ambarishbiswas/imfungi:v1 -h

#### imFungi examples:

     docker run -v $(pwd):/tmp --rm -it --name imfungi ambarishbiswas/imfungi:v1 -test fast

     docker run -v $(pwd):/tmp --rm -it --name imfungi ambarishbiswas/imfungi:v1 -test full

     docker run -v $(pwd):/tmp --rm -it --name imfungi ambarishbiswas/imfungi:v1 -sample_id ABCD -f test.fastq -o EFGH

     docker run -v /tmp:/tmp --rm -it --name imfungi ambarishbiswas/imfungi:v1 -sample_id ABCD -1 test_R1.fastq -2 test_R2.fastq -o EFGH

     docker run -v /tmp:/tmp --rm -it --name imfungi ambarishbiswas/imfungi:v1 -sample_id SRR390728 -sra SRR390728 -o SRR390728_output

     docker run -v /tmp/ABCD:/tmp --rm -it --name imfungi ambarishbiswas/imfungi:v1 -sample_id XYZ -ena ERR2105523 -o output_of_ERR2105523


#### imFungi commandline parmeters:

  Options must be provided:
  
	-sample_id                            XYZ1                      [Every imFungi run must be supplied with a sample_id.]  
 	-o/-out_dir                           a_folder_name             [A folder with the provided name will be created in the current directory; No '/' or '\' is accepted]

  Options for input sequence(s) in FASTA format:
  
 	-i/-f                                 input_fasta_sequence_file [either gzip compressed (i.e. with extension .gz) or uncompressed FASTA sequence file. Supported extensions are .fa, .fna, 
                                      	                          	.fasta ]



  Options for input sequence(s) in FASTQ format:
  
 	-s/-r                                 input_fastq_sequence_file [either gzip compressed (i.e. with extension .gz) or uncompressed FASTQ file. Supported extensions are .fq, .fastq ]
 	-1                                    forward_fastq_file        [Specify the forward reads FASTQ file of paired-end libraries. Supported extensions are .fq, .fastq with/without .gz]
 	-2                                    reverse_fastq_file        [Specify the reverse reads FASTQ file of paired-end libraries. Supported extensions are .fq, .fastq with/without .gz]
 	--12                                  interleaved_fastq_file    [Specify a FASTQ file with forward and reverse reads interleaved. Supported extensions are .fq, .fastq with/without .gz]


  Options for publicly available (meta)genomic NGS reads:
  
 	-sra                                  SRA_Accession             [Accession of (meta)genomic NGS reads from NCBI SRA database]
 	-ena                                  ENA_Accession             [Accession of (meta)genomic NGS reads from ENA database]


  Options for reads specific operations:
  
 	-read_error_correction                0/1                       [Default is set to 1; imFungi uses BBmap suite of tools for error correction] 
 	-min_read_length                      20                        [Default value set to 20]
 	-subsample_reads_to                   N                         [A positive integer 'N' can be provided to create a subsample_reads file using randomly selected reads; This is useful where 
                                      	                          	multiple metagenomic libraries with different reads depth are being analyzed; ]

 	-subsample_reads_min_length           N                         [Default value set to 75; Any positive integer is supported; Remaining reads after read_error_correction Reads will be checked 
                                      	                          	for length >= N before including them in the subsample_reads file ]



  Options for reads assembly:
 	-use_reference_sequence_from_assembly 0/1                       [Default is set to 1, which means assembly will be done on the inputted reads/sequences in the inputted prior constructing 
                                      	                          	a sample specific reference library of marker sequence, which will then be used for classifying archaea, bacteria and viral 
                                      	                          	reads ]

 	-assembler                            Spades/Megahit            [Default set to SPAdes; to use Megahit assembler specify '-assembler megahit']
 	-assembled_contigs_file               FASTA_seq_file            [specify an existing (multi)FASTA sequence file. Taxonomy prediction of the contigs will be carried out prior constructing 
                                      	                          	a sample specific reference library of marker sequence. Either gzip compressed (i.e. with extension .gz) or uncompressed 
                                      	                          	FASTA sequence file is expected ]

  Options for blast based intron search:
 	-blastn_word_size                     28                        [Default is set to 28; Although a smaller word length (e.g. 21) can identify more fungal sequences, but it also increases the risk of false predictions;]


  Options for Taxonomy search:
 	-taxa_assignment_of_contigs           0/1                       [Default is set to 1, i.e. the taxonomy will be predicted on assembled contigs;]
 	-taxa_assignment_of_reads             0/1                       [Default is set to 1, i.e. the taxonomy will be predicted on the inputted reads (or reads from the subsampled reads file if 
                                      	                          	opted); ]


  Options for additional reporting:
 	-kraken_style_output                  0/1                       [Default is set to 0; Use '1' to create a output file similar to what kraken, CLARKS and kaiju provides]


  Other options:
 	-T/-threads                           N                         [By default the program uses 4 threads; Any positive integer is accepted]
 	-q/-quiet                             0/1                       [Default is set to 0, which shows program step-by-step logs; Use '1' to turn of the logging; Note, a file log.log will still 
                                      	                          	be created in the output_folder ]

 	-h/--h/-help/--help                                             [Shows this help]
 	-v/-version                                                     [Shows imFungi Version ]
 	-rm                                                             [removes the output folder specified by -o or -out_dir (if exists);]
 	-clean                                                          [removes all intermediate temporary files once a process is finished;]
 	-test                                                           [supported options are: full, fast or custom; If -test is provided, imFungi uses pre-supplied fastq files for testing the pipeline; 
                                      	                          	the option 'full' sets parameters for most-comprehensive imFungi run; The option 'fast' does the least-comprehensive imFungi 
                                      	                          	run where reads are classified directly using imFungi RLMS; The 'custom' option takes all other user parameters except input 
                                      	                          	fasta/fastq related parameters; ]



General information:
  For all queries, bugs and updates (and even suggestions how to make imFungi better) please email Ambarish Biswas [ambarishbiswas@gmail.com]
  Latest version and other informationcan can be found at: https://github.com/ambarishbiswas/imfungi_v1.0

  The imFungi V1.0 tool is developed as a proof-of-concept for the upcoming/associated paper and is free for academic use. If you use imFungi for your reasearch, please cite the associated paper. 
