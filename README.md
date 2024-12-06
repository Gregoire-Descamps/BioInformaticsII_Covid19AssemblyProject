This is a projec given in 1st year Data management in Biosiences, Bioinformatics course.


## Instructions
For the assembly, only the R1 fastq file should be used (not the R2)

 - Step 1: Indexing k-mers found in the reads of the R1 fastq file in a dictionary. The keys in the dictionary are the k-mers and the values are the number of occurrences of the k-mers in the reads. E.g.: {"ATGGATTA": 3, "GATGATAG": 26}. Your script should be tested with at least three k values: [20, 25, 30]. Make sure you index only once a k-mer and its reverse complement, for instance by choosing randomly one of them to include in the dictionary.

 - Step 2: Construct the DBG by looking for perfect overlaps of k-1 length between the k-mers of the index you have constructed in step 1. For that, start by a random k-mer from your index and append it to the right and then to the left, using the k-1 overlap strategy. Be aware of using each k-mer once in your graph.

    E.g.:
        You start with "ATGGATTA"
        You append it to the right

      ATGGATTA
       |||||||
       TGGATTAT
        |||||||
        GGATTATC
         |||||||
         GATTATCG

        Once there is no more possibility to append to the right, you append to the left

     TATGGATT
      |||||||
      ATGGATTA
       TGGATTAT
        GGATTATC
         GATTATCG

    For each overlap check, do not forget to take into consideration the reverse complement of the k-mers.

 - Step 3: Once there is no more possibility to append to both right and left, generate the consensus sequence from your overlapping k-mers as follows and save it as a contig in an output fasta file.

     TATGGATT
      ATGGATTA
       TGGATTAT
        GGATTATC
         GATTATCG
    ___________________
       TATGGATTATCG --> consensus sequence

 - Repeat Step 2 and Step 3 until all the k-mers in your dictionary have been checked.

### Remarks:

 - Rare k-mers in Illumina sequencing could reflect a sequencing error. This is why, do not integrate a k-mer in your graph if its number of occurrences is below a threshold (you can test 2 and 20 as thresholds).
 - There is a possibility that a k-mer in your graph has more than one overlapping k-mer. In this case, you will have forks (more than one edge outgoing from the same vertex) in your graph as follows:

    TATGGATT --> ATGGATTC
             --> ATGGATTG

In this case, append both branches and propose a strategy to construct different contigs with these branches or to select one of them to construct one contig. You can eliminate too short branches (tips).

### Output:

 - As an output of your script, you should get a fasta file containing one or multiple contig(s): this will be your SARS-CoV-2 genomic assembly.
 -  You should assess your assembly(ies) quality by computing N50 and L50 and by comparing it/them to the provided SARS-CoV-2 reference assembly using a dotplot. Online tools are available for this comparison.


## Links to data URL

Genomic assembly: https://www.ebi.ac.uk/ena/browser/view/GCA_011545535.1

Raw reads: https://www.ebi.ac.uk/ena/browser/view/SRR21719088
