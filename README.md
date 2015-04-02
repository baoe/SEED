##### Contents
[Overview] (#overview)  
[Copy right] (#copyright)  
[How to cite SEED?] (#cite)  
[Short Manual] (#manual)  

<a name="overview"/>
### Overview
SEED is a software for clustering large sets of Next Generation Sequences (NGS) with hundreds of millions of reads in a time and memory efficient manner. Its algorithm joins highly similar sequences into clusters that can differ by up to three mismatches and three overhanging residues.

<a name="copyright"/>
###Copy right
SEED is under the [Artistic License 2.0](http://opensource.org/licenses/Artistic-2.0).

<a name="cite"/>
### How to cite SEED?
If you use SEED, please cite the following paper:  
Bao E, Jiang T, Kaloshian I, Girke T (2011) SEED: Efficient Clustering of Next Generation Sequences. Bioinformatics: [epub](http://www.hubmed.org/display.cgi?uids=21810899).

<a name="manual"/>
### Short manual
1. System requirements

   SEED is suitable for 32-bit or 64-bit machines with Windows, OS X or Linux operating systems. At least 4GB of system memory is recommended for clustering larger data sets.

2. Installation

   The downloaded .cpp file can be compiled as follows:  
   * On Mac/UNIX/Linux systems, execute on the command line: `g++ -o SEED SEED.cpp`
   * On Windows systems, the code can be compiled under the Visual C++ environment.

3. Input

   Only FASTQ format is supported in the current version. The sequence length should be between 21 bp and 100 bp with the max variation of 5 bp.

4. Using SEED

   ```
   SEED --input input.fastq --output output.txt [--mismatch M] [--shift S] [--QV1 L] [--QV2 U] [--fast/short] [--reverse] [--input2 input2.fastq]
   ```

   --mismatch is the maximum number of mismatches allowed from the center sequence in each cluster (0 - 3, default 3).  
   --shift is the maximum number of shifts allowed from the center sequence in each cluster (0 - 6, default 3).  
   --QV1 is the threshold for the base call quality values (QV) that are provided in the FASTQ files as Phred scores. SEED ignores those mismatches where the sum of the Phred scores of the mismatching bases is lower than the specified QV1 threshold value (0 - 2 * 93). The default value for QV1 is 0.  
   --QV2 is another QV threshold. It prevents co-clustering of sequences where the sum of all mismatched positions is higher than the threshold value (0 - 6 * 93). The default value for QV2 is 6 * 93.  
   --fast uses a bigger spaced seed weight to save running time. It is only applicable for sequences longer than 58 bp and may need more memory.  
   --short is to use a smaller spaced seeds weight for sequences as short as 21 bp. This setting often results in longer compute times.  
   -- reverse is to co-cluster sequences in sense and anti-sense orientation (reverse and complement).  
   --input2 specifies the paired sequences so that paired-end library can be clustered. In current implementation, no shift is allowed for this option, and if --reverse option is specified minimum sequence lengths of both pairs should be the same.

5. Output

   SEED outputs two files: a SEED file and a FASTQ file. The outputted FASTQ file has the same format as the input FASTQ file, but it contains only the center sequences and their quality scores for each cluster with one or more members. In other words, it is the filtered version of the input FASTQ file where the redundant sequences have been removed. The SEED file has a tabular format that is explained in the following table. The third column in this table is only available if the --reverse argument has been specified.

   |Cluster ID                   |Sequence ID                  | Is Reversed                 |
   |:----------------------------|:----------------------------|:----------------------------|
   |Center sequence for cluster 0|                             |                             |
   |0                            |Sequence  ID from input file |1                            |
   |0                            |Sequence  ID from input file |0                            |
   |Center sequence for cluster 1|                             |                             |
   |1                            |Sequence  ID from input file |1                            |
   |1                            |Sequence  ID from input file |0                            |
