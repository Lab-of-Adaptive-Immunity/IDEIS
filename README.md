
-[//]--[//]--[//]--[//]--[//]--[//]--[//]--[//]--[//]--[//]--[//]--[//]--[//]--[//]--[//]--[//]-  
/]--[//]--[//]--[//]--[//]--[//]--[//]--[//]--[//]--[//]--[//]--[//]--[//]--[//]--[//]--[//]--[/  

                       _______  _____    _______  _______  _____
                      [__   __]|  __ \  |  _____][__   __]/  ___\
                         | |   | |  \ \ | |         | |   | |
                         | |   | |   | \| |___      | |   \ \__
                         | |   | |   | ||  ___]     | |    \__ \
                         | |   | |   | /| |         | |       \ \
                       __| |__ | |__/ / | |_____  __| |__ ____| |
                      [_______]|_____/  |_______][_______]\_____/


-[//]--[//]--[//]--[//]--[//]--[//]--[//]--[//]--[//]--[//]--[//]--[//]--[//]--[//]--[//]--[//]-  
/]--[//]--[//]--[//]--[//]--[//]--[//]--[//]--[//]--[//]--[//]--[//]--[//]--[//]--[//]--[//]--[/  



This script runs a pipeline that identifies CD45 or PTPRC isoforms in 10X single cell data processed by Cell Ranger beforehand. The software uses a list of exons, a list of  rules and the length of reads used for mapping to build a list of customized and optimized transcripts, to which then reads from CD45/PTPRC locus are mapped using salmon alevin. The quantification then tells us which isoforms are present in which cells. Main purpose is to identify CD45RA and CD45RO isoforms in Homo Sapiens from 5' sequencing, but it also works in Mus Musculus and for exons B and C. However, due to the location of key exons in CD45/PTPRC the usage is extremely limited in 3' sequencing.

 /]]]]]]]]]]]]\  
| Requirements  
/

To run this, you need to have following software installed:

* python3.7 or higher (tested on version 3.8.10)
  * python packages os, sys, glob, gzip, pathlib, argparse, inspect, shutil and subprocess
* salmon alevin 1.9
* samtools
* R with Rscript
* bamtofastq-1.4.1 from 10X

You can download bamtofastq-1.4.1 here: <https://github.com/10XGenomics/bamtofastq/releases/tag/v1.4.1>. Once you dowload it, put in main directory (where 'IDEIS_main.py' is located).

If you don't have these, you won't be able to use this script. If you have them installed locally, you need to precise *absolute* paths to them in the script config.txt.

Moreover, it assumes that the data for which you want to use it were already analyzed by cellranger count, as it requires .bam generated by it. If not, you need to run cellranger count on the entirety of the dataset. This is because the alignment is needed to filter the reads that align with CD45/Ptprc gene.


 /]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]\  
| Installation and preparation before first run  
/

This is a simple python script, so no installation is needed. However, you need to get bamtofastq-1.4.1 (download it at https://github.com/10XGenomics/bamtofastq/releases/tag/v1.4.1 and put it in directory with IDEIS_main.py in it). The first step is the creation of salmon indices or references that this script needs, so you have to have files that allow you to create such references. Those are .fasta files containing exons of interest and files with rulesets to build them. Both files are provided for Homo Sapiens (version GRCh38), Mouse C57BL6/j (version GRCm38) and BALB/cJ (version GRCm38) - for others you need to prepare your own references. You only need references for exons A,B,C and those coming before and after (the length is arbitrary, but should be enough to cover the length of read). For custom references for your own desired species, the exons may be obtained on site of Ensembl:

https://www.ensembl.org/index.html

Where you can search for desired species and get latest exons for it. The rule set may be copied from other species. You'll need also the genome range for prior sub-filtering; Ensembl site should help you with finding the range where Ptprc is located. More information below.


 /]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]\  
| Custom file with exons and rules  
/

This software builds artificial "transcriptome" from provided exons and rule set. Rule set provides instructions how to chain exons to generate transcripts, while exons are the sequences themselves.

For custom analyis you'll need yoyur own exon file. In case of Ptprc, this means exons R,A,B,C and O. Exons R and O can technically be (and for provided files are) multiple exons chained accordingly to sufficiently cover the length of read. The headers should be copied from either Mus Musculus or Homo Sapiens (the only difference is the capitalisation). Under the headers you copy the sequence from Ensembl accordingly, and use rule-set for species you copied headers from. Finally, search range of Ptprc gene for species of interest (approximation is sufficient as long as it fiably covers the whole gene). Now you should be ready to build reference using the optional parameters

Example of exon R for Mus Musculus:

\>Ptprc-R R
TTTGTTCTTAGGGTAAGAGAGTAGGAAACTTGCTCCCCATCTGATAAGACAGAGTGCAAA
GGAGACCCTATTTCTTAGGGGCACAGCTGATCTCCAGATATGACCATGGGTTTGTGGCTC
AAACTTCTGGCCTTTGGATTTGCCCTTCTGGACACAGAAGTCTTTGTCACAG
GGCAAACACCTACACCCAGTGATG

If you desire to customize headers, you need to also cusomize rulesets afterwards. A rule is a string that defines how transcriptome should be built and looks like this:

A:R,A,B,C,O:Ptprc-RABCO-1:Ptprc-RA

The important fields are separated by ':':

* First field indicates the exon of interest for this transcript. This will be the central part of newly generated transcript.
* Second field shows the chaining of exons. The exon from the first field is the part that interests us, while the chaining will be used to build flanks before and after exon of interest of convenient length determined from the length of provided reads and desired offset. For this reason it is important that the flanks are long enough.
* Third field shows the name of transcript.
* Last field shows the name of gene, here the name of isoform of interest, in this case RA.

In some cases it is not the exon that interests us but junction. In such case the first field is left empty and the junction is instead indicated as a dash '-', and the rest of fields stays the same. Example:

:R,-,O:PTPRC-RO-13:PTPRC-RO


 /]]]]]]]]]]]]]]]]\  
| Running software  
/

Run python IDEIS_main.py -h from its directory to see all options or just run python IDEIS_main.py it without any option.

If you don't have created references, just point out the directory where you want to create them. First run will create them for you, and any subsetquent runs will re-use these references if the same path for directory is used and read + cutoff parameters are the same; the references are adjusted to the length of read and desired cutoff. You specify the path to reference by 'reference' (R) parameter.

The .bam file you need is located in outs/ of your directory created by cellranger, along with .bam.bai file. Usually they are named as 'possorted_genom_bam.bam' and 'possorted_genom_bam.bam.bai'. These are the alignment files you need to point to this script to 'bam_path' (B) parameter. 

The parameter 'output' (O) describes where output will be created. Complete output has following structure:

* Fastqs(.fastq generated from .bam)  
* Results(results containing following subdirectories):
* * alevin (alevin output)  
* * Iso_counts (raw matrix in)  
* * Datasets (if --data-set was provided; contains data set wirth extra assay)  
* * whitelist.csv (if --data-set was provided; whitelist generated from data set)  

You are required to use one of whitelisting options:
- using parameter --force-cells; this will force the number of cells for analysis. Least reliable, so use only in emergencies.
- using parameter --whitelist, a path to whitelist, which will be used to limit the cells of the interest; recommended option.
- using parameter --data-set, a path to .rds with SeuratObject, from which whitelist will be generated from the column names of the object; recommended option.

WARNING! If the script is rerun with same options, the extraction of fastqs will be skipped (as their presence is auto-detected) and the script will throw a warning, but it will re-do the rest of analysis while overwriting results. You have been warned.

 /]]]]]]]]]]]]]]]]]]]]\  
| Example  
/

The most basic command is to use just pre-defined options. For example, when using sequences from C57BL/6 mouse:

* python IDEIS_main.py PTPRC_mouse_alevin possorted_genome_bam.bam Output

For human, you need to specify 'Homo-Sapiens' for option -g:

* python IDEIS_main.py PTPRC_human_alevin possorted_genome_bam.bam Output -g Homo-Sapiens

This would employ default --force-cells parameter, which is probably the least fiable option. Instead, you should use either whitelist or sata set containing SeuratObject to determine the cells of interest:

* python IDEIS_main.py PTPRC_human_alevin possorted_genome_bam.bam Output --whitelist whitelist.csv -g Homo-Sapiens

* python IDEIS_main.py PTPRC_human_alevin possorted_genome_bam.bam Output --data-set SeuratObject.rds -g Homo-Sapiens

If you want to use custom reference, you need to provide your fasta of exons (option --fasta, path to file), list of rules (--rule-set, path to file), and gene range (--gene-range, str in format chrX:Y-Z, similar to string used for samtools filtering).

* python IDEIS_main.py PTPRC_custom_alevin possorted_genome_bam.bam Output -g Custom --gene-range chrX:Y-Z --fasta Custom_fasta_of_exons --rule-set Custom-rule-set--whitelist whitelist.csv

To see all options, execute python IDEIS_main.py -h.

 /]]]]]]]]]]]]]]]]]]]]\  
| Technical Limits  
/

- The method does not work very well on 3' sequenced data. This is because such methods do not support sequencing (much) end of transcript where CD45 key exons are located. You might get some data from particularly deep sequenced data, but still do not expect much.
- The software does not make use of R1 even if sequenced beyond recommended settings (for example cases where R1 has 150bp). Though tests so far show no great impact of this caveat, we're searching for imporvement of algorithm so it will also use R1 data.

 /]]]]]]]]]]]]]]]\  
| Troubleshooting  
/ 

Problem: My data have too low mapping rate.
Solution: There may be multiple causes to this:

* Using software on 3' data. Here the usual mapping rate is about 0.5%. This is because the key exons are located from 5' end of the sequence. As said above, don't expect much from this software when using it on 3' data.
* Using bad gene range for custom genomes. Please be sure you delimit gene for range where CD45/PTPRC is truly located, otherwise the reads will not correspond to the gene.
* Using .bam mapped to incompatible reference when mapping on one of three base species. The references used here are GRCh38 for Homo Sapiens and GRCm38 for both mice. If you use, for example, GRCm39, you need to provide custom reference as shown above (just copy moiuse exon list and replace it with appropriate exons; or you can just precise gene range).

 /]]]]]]]]]]]\  
| License  
/

IDEIS (c) by Lab of Adaptive Immunity from Institute of Molecular Genetics of the Czech Academy of Sciences

IDEIS and all its parts are licensed under a Creative Commons Attribution 4.0 International License.

You should have received a copy of the license along with this work. If not, see <https://creativecommons.org/licenses/by/4.0/>. 

 /]]]]]]]]]]]\  
| Acknowledgements  
/

This project was supported by the National Institute of Virology and Bacteriology (Programme EXCELES, LX22NPO5103 to Ondrej Stepanek) - funded by the European Union - Next Generation EU. 

 /]]]]]]]\
| Contact
/

Feel free to send a mail to juraj.michalik@img.cas.cz if you encounter any bugs.
