# BIO-410-FINAL-PROJECT
## Background
The data consist of 6 samples from the organism *Zaire ebolavirus*.The organism is a species of viruses in the family *Filoviridae*. (https://eol.org/pages/741151)
The purpose of this project was to create a phylogenetic tree from 6 samples of Zaire ebolavirus in order to determine the evolutionary relationships between the samples.
## Methodes
**1.The Anaconda/Miniconda program was installed to facilitate the easy installation of bioinformatics tools, specifically MEGAHIT.**

**2.We install the MEGAHIT assembler from the Bioconda channel, which assembles paired-end reads into contigs (pieces of the genome).**

 `conda install -c bioconda megahit` 
-	To run any commands, we must open WSL and navigate to the folder containing our FASTQ files.
-	Then we will be given a command `megahit -1 sim_t1_1.fq -2 sim_t1_2.fq -o t1_out.` 
-	This command will be run 6 times (t1 through t6).
`-1 sim_t1_1.fq`
This is the forward reads file for sample t1.
`-2 sim_t1_2.fq`
This is the reverse reads file for sample t1.
-	Paired-end sequencing will produce two files: _1.fq (forward reads) and _2.fq (reverse reads).
-	Then we will get -o t1_out. This will tell MEGAHIT where to put the assembly results. It will create a folder name `t1_out`.
`megahit -1 sim_t1_1.fq -2 sim_t1_2.fq -o t1_out`
`megahit -1 sim_t2_1.fq -2 sim_t2_2.fq -o t2_out`
`megahit -1 sim_t3_1.fq -2 sim_t3_2.fq -o t3_out`
`megahit -1 sim_t4_1.fq -2 sim_t4_2.fq -o t4_out`
`megahit -1 sim_t5_1.fq -2 sim_t5_2.fq -o t5_out`
`megahit -1 sim_t6_1.fq -2 sim_t6_2.fq -o t6_out`

-	Each sample produces `t1_out/final.contigs.fa`, which is the fasta file containing assembled genome fragments. 

**3. In R program, the following code was executed:**
**Biostrings allows R to read DNA sequences from FASTA files, while DECIPHER performs alignment and builds phylogenetic trees. To install these packages:**
BiocManager::install("Biostrings")
BiocManager::install("DECIPHER")

**Load required libraries, so we can use them.**
library(Biostrings)
library(DECIPHER)

**Read in megahit assembly results for each sample.**
- `allcontigs <- c()`  Creates an empty object allcontigs, loops from sample 1 to 6, reads each sample’s `final.contigs.fa`,  and appends all contigs into one big list.
- `allcontigs <- c()`
- `for(i in 1:6){`
- `contigs <- readDNAStringSet(paste0('C:/Users/sabye/Downloads/t', i, '_out/final.contigs.fa'))`
- `allcontigs <- c(allcontigs, contigs)}`
  
**This loop outputs a list, so we flatten the list. The loop creates a list of DNAStringSets.**
`allcontigs <- do.call(c, allcontigs)`

**There are some small fragments of genomes in the assembly**
- We are only going to align the biggest parts (aka the ones that are bigger than 5 kbp; the genome itself is around 18 kbp)using the command: 
- `toalign <- allcontigs[which(nchar(allcontigs) > 5000)]`
- `names(toalign) <- 1:length(toalign)`
-	The viral genome is ~18 kb. Small contigs (<5 kb) are junk or assembly fragments, so it keep only the large pieces.
-	`names(toalign)` gives each sequence a simple numeric name.

**We command `alignment <- AlignSeqs(toalign)`**
-	This will make DECIPHER align all sequences using progressive alignment. It will compute k-mer distances, cluster similar sequences, align them, and refine the alignment. 
-	This produces a multiple sequence alignment.

**View the alignment.**
-	 We command `BrowseSeqs(alignment, htmlFile = 'final_project.html')`
-	This will create an HTML file showing the alignment visually.
  
**Build a Phylogenetic Tree**
-	We command `tree <- Treeline(alignment, method = 'ML', showPlot = T)`
-	The program builds a maximum likelihood (ML) phylogenetic tree, tests multiple substitution models (JC69, K80, HKY85, GTR, etc.), picks the best model (T92), optimizes tree likelihood, and plots the final tree.
-	This gives us the evolutionary relationship between our assembled genomes.

**Understanding the model selection output**
- The selected model was: T92
Frequency(A) = Frequency(T) = 0.293
Frequency(C) = Frequency(G) = 0.207
Transition rates = 0.973
Transversion rates = 1
The T92 model considers unequal base frequencies and varying transition/transversion rates. Its parameters encompass nucleotide frequencies, mutation rates, and the optimal evolutionary model for our sequences.
## Results
![ALT TEXT](https://github.com/renny6018-bit/BIO-410-FINAL-PROJECT/blob/f24f511bf6a8749949a41a9633aad1d27af4bbeb/Rplot02%20final%20project%20tree.png)
### Explanation
The phylogenetic tree reveals that the closest connections are found between species 2 and 3, as well as between species 1 and 4. This indicates a strong relationship among these species. The connecting lines, which are the shortest, illustrate that these groups share a common node. On the other hand, sample 5 represents the closest secondary relationship to samples 2 and 3, while sample 6 serves as the secondary connection to samples 1 and 4.
Based on the number of terminal nodes at the end of each branch of this phylogenetic tree, it can be infere that the samples represent six distinct individuals.
