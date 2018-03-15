
##Introduction to Genome Assembly

### Summary

Welcome to the denovo assembly tutorial! This is a multi-part course of interactive tutorials to help you get to grips with genome assembly. We will guide you through each of the necessary steps to putting together a full genomic sequence, from quality checking your reads, to constructing contiguous sequences (contigs) and scaffolds. Along the way, we will use multiple computational tools to show you the differences your tools can have on your final assembly. We will also guide you through common problems you may face, show you examples of good and bad quality reads and assemblies, and give you a chance to test what you have learnt. By the end of this course, you will have the necessary understanding and ability to assemble your own genomes. This tutorial has been designed with the total beginner to genomics in mind, so we hope you will find it easy to follow, no matter your starting ability. Let's get started!

 Why assemble?

Modern sequencing technologies have come a long way since Sanger and colleagues first introduced their method of DNA sequencing in 1977. We can now easily produce huge amounts of genetic data in a short space of time and at relatively little cost thanks to high-throughput sequencing techniques. However, no sequencing technologies produce a full string of genomic code; instead, they generate a number of fragments containing a short nucleotide sequence. These fragments need to be fitted together like a puzzle by computational tools in order to generate the full sequence. There are different levels of sequence assembly, starting with reads straight from the sequencer, contigs pieced together from the reads, scaffolds that join together the contigs, and the full genome assembly containing no missing sequence. The tests we want to do on the sequence later on determines how assembled we need our genome, but many applications will require an assembly to scaffold level or above. 

Challenges and problems faced 

Assembling a genome comes with numerous problems, from the initial sequencing all the way to full genome assembly. There is ultimately a compromise in high-throughput sequencing between accuracy and time, so while our current technologies will generate a lot of data relatively quickly, Sanger sequencing remains the most accurate. A genome is sequenced in fragments, which can be long or short depending on the sequencer; short fragments have a high base pair accuracy whereas long fragments have a low base pair accuracy. Errors occur when the sequencer misidentifies a nucleotide, skips some sequence or inserts some sequence. We can try to work out where the errors lie computationally using sequence coverage. There are two main reasons why we have to sequence multiple times the actual size of the genome:

1)	Coverage: A sequencer doesn't just sequence the genome in one long string; the genomic material is broken down into multiple copies of random fragments in the sequencing preparation steps. There are too many fragments to be read by the sequencer, so the sequencer samples a random selection of the fragments from the fragment pool presented to it. Because of this, in order to cover the whole genome, the coverage level needs to be set so that every bit of genome is likely to be read. Imagine throwing a 6-sided die and trying to see all sides at least once - you would need to throw the die way more than 6 times to cover all 6 sides, because some sides you will see repeatedly. The same principle applies to genome sequencing.

2)	Errors: We're assuming all of the observations made are perfect, but in reality there are parts of the genome that look like each other. This leads to errors, and errors detract from coverage so all parts of the genome need to be sampled enough to distinguish between what's an error and what's a true observation. Imagine your die had two 2's but no 5; by throwing the die you could tell me very quickly that there was an error, but by observation alone you would find it hard to tell which 2 is the true observation and which is an error. If you roll the 2 enough times, the sticker making it look like a 2 peels off to reveal a 5. You gave your die enough coverage to tell what the true observation in that position should be.

In this first tutorial, we will guide you through some of the basic principles and methods we use to begin assembling a genome. This will give you an understanding of the concepts behind some of the tools we use. We will start with a simple artificial genome (a random string of nucleotides) to introduce our methods, before increasing the complexity and adding to your understanding. 

Within this first tutorial, our genomes are homozygous and contain 4 chromosomes with 50x coverage provided as a set of paired 150bp reads (called paired-end reads). They have been randomly synthesised to represent increasing complexity when it comes to assembly, which we will introduce step-by-step. The first genome we will be working with, genome_1, has a 1% random error rate and a fragment size of 350bp but no complexities. We will introduce duplications to you in genome_2, which also contains a random error rate of 1% and fragment size of 350bp.

Contents of this tutorial:

1.	Looking at k-mer coverage 
2.	A simple assembly
3.	De Bruijn Graph, coverage and N50
4.	Re-assembling
5.	Comparing the PE reads to assembly
6.	Checking insert size and distribution
7.	Introducing some complexity - Duplications
8.	Changing parameters
9.	Changing the k-mer
10.	Combining cutoff and k-mer alterations

This tutorial assumes a very basic understanding of command line; if you have no previous experience in this, you may find this http://rik.smith-unna.com/command_line_bootcamp/?id=ga33gnprbdt useful.

1. Looking at k-mer coverage

A k-mer is any sequential combination of nucleotides in our sequence of k length. By this we mean that our sequence is split into multiple k-length sections, each containing a shift of one nucleotide along the sequence so that the second nucleotide becomes the first, and the last nucleotide becomes the second-to-last, with a new nucleotide filling the last position. Eg:

If our sequence is this:
```
GGTCAGTCAGTCTTCGATTG
```

And we have k-mers of 7 nucleotides in length (7-mers), our first k-mer would be:

```
GGTCAGT
And our second:
    GTCAGTC
And our third:
       TCAGTCA
And so on…
          CAGTCAG
             AGTCAGT
                GTCAGTC
                   TCAGTCT
```

Each k-mer appears only once and we would hope to find that most k-mers are a unique combination of nucleotides. This we achieve by changing the length of the k-mer; a longer k-mer parameter will more likely give unique k-mers but at the cost of coverage, so what we aim for is something in the middle. Estimating the k-mer coverage gives us an idea of the k-mer variance and errors in our reads. Change directory to genome_1 and type in the following to produce a k-mer histogram:

```
kat hist -o scer_pe_hist -h 80 -t 8 -m 27 -H 100000000 pe_reads_R1.fastq pe_reads_R2.fastq
```

In this command, "-h" specifies the high count value for the histogram and "-H" specifies the hash size, needed for k-mer counting. Our k-mer length (-m) is 27, which is suitable for the dataset we are using. You can see a full list of flags by typing kat hist --help. We can see from the output how many times each of the k-mers appear in our reads, with a symmetrical distribution around the 60 value. This shows us that our estimated k-mer coverage is about the 60 mark and that we have a reasonably small variance, which is good. We will want to use the k-mers in this distribution for our assembly. The sharp peak at ~5 are k-mers from erroneous reads, which we will always find in a real dataset due to errors in our sequence; here, we know there are 1% errors in our synthesised set, which is what we see.

2. A simple assembly

Now we can move onto constructing a simple assembly of our simple reads. VELVET is an assembler package containing algorithms for the assembly of short reads. It works by creating De Bruijn graphs, which, in basic terms, splits the sequences down into k-mers and then constructs a graph of k-mers that connect by overlapping. This makes VELVET a good tool for taking a quick look at how our data can assemble. Velveth is the first command we use, which breaks our sequence into individual k-mers ready for the next command. Velvetg is the second command, which uses the k-mers produced by velveth to construct a De Bruijn graph. It performs simplification and corrects errors on the graph, before finally constructing and outputting the contigs. 

To run the assembler we use the paired-end read command:

```
/velvet/velveth assm_dir 31 -fastq -shortPaired -separate pe_reads_R1.fastq pe_reads_R2.fastq
```

Followed by: 

```
/velvet/velvetg assm_dir
```

Great! Now we've assembled our reads. You can view the new files in the assm_dir directory! The files we will be using later are contigs.fa and LastGraph. contigs.fa contains the information about the contig assembly we just made with VELVET, which we can use for quality checking. LastGraph contains information for visualisation of the assembly, which we will look at in the next step. 

 3. De Bruijn Graph, coverage and N50

Remember that VELVET works by producing a Dr Bruijn graph? You can have a look at the De Bruijn graph in BANDAGE by going: 

```
file -> load graph -> select "Last Graph" in the assm_dir directory -> hit "Draw graph." 
```

Notice how many fragments there are with alternative paths, each one of these is called a node. Without specifying a cutoff point for the node coverage, the algorithm can't decide which path is the true path, so includes a web of possibilities in the output.

Take some time to assess your sequence. You can zoom into sections and click on the different coloured nodes to view their length and coverage (depth) on the right hand side. You can also drag nodes around to see more clearly where they interact with other nodes. See if you can tell where the true path through the graph nodes lies by looking at their coverages. If you set the "Scope" to "Min: 10" on the left hand side before you hit "Draw graph," you'll notice that all the low coverage nodes have been removed, leaving only the sequence with a coverage of >10. 

If you click "More info" on the left hand side of BANDAGE, you'll get a pop-up window of some statistics on the assembly, including the N50 value. N50 is a commonly used statistic to assess the contiguity of an assembly. An N50 of 1kb means that 50% of the assembly is contained in sequences of 1kb or larger. A high N50 indicates long contigs and a low N50 indicates a fragmented assembly. 

N50 is often used to imply quality of the set of contigs; in the ideal assembly, we would have one contig representing each chromosome, giving us a high N50 value. Therefore, high N50 is often taken as being solely indicative of good quality sequence and low N50, to be of poor quality sequence. However, N50 is not so reliable and genome assembly not so straight forward as to rely on one statistic for quality representation. N50 is not a measure of the accuracy of the assembled sequence and can also be dramatically increased by removing shorter contigs. This leaves the statistic vulnerable to manipulation and can easily become misleading when no other information on the assembly is considered. In short, a sequence with high N50 does not necessarily represent the best possible assembly from a given set of reads, and a sequence with a lower N50 may, in fact, provide the more comprehensive assembly. What we can mostly tell from the N50 is how easily we would be able to work with our data later on, assuming the assembly was assembled correctly; longer sequences are more easy to work with because they contain many full genes

In our assembly, we have an N50 value of 55bp, which we can see indicates a fragmented assembly in the De Bruijn graph. Here, our N50 is indicative of a poor assembly, but does not conclude that the assembly is poor. With other evidence such as the De Bruijn graph, we can tell that this is not ideal; because N50 represents how easily we would possibly work with our assembly, a small N50 like this means that genes will be fragmented. To improve this, we would like to remove the low-coverage errors and keep only the true sequence.



4. Re-assembling

Nevermind! This is a learning experience and we will just have to alter our parameters. The first step of our assembly is the same:

```
/velvet/velveth new_assm_dir 31 -fastq -shortPaired -separate pe_reads_R1.fastq pe_reads_R2.fastq
```

But the second step includes a few adjustments:

```
/velvet/velvetg new_assm_dir -cov_cutoff auto
```

We are using ```-cov_cutoff``` to cut out the k-mers of lower coverage - the low coverage nodes in our De Bruijn graph - these are most likely erroneous and may be causing problems in our assembly. Take another look at the k-mer spectra histogram we produced in part 1 - remember those erroneous reads we identified at ~5? Those are very low coverage, and are what we need to cut out. We'll let the computer decide where the cutoff should be for now by specifying "auto." Check the De Bruijn graph in BANDAGE once again. Notice the difference? We now have a clean graph showing the four chromosomes from our data, with no ambiguities.

5. Comparing the PE reads to assembly

A k-mer spectra plot is a useful tool for being able to compare the k-mer content between the assembly and the original reads it was assembled from. This allows us to check for any missing content or misassembly. We can use both sets of contigs we previously generated for comparison. 

Type the following:

```
kat comp -o scer_pe_v2_ctgs -t 8 -m 27 -H 100000000 -I 100000000 "pe_reads_R1.fastq pe_reads_R2.fastq" assm_dir/contigs.fa
```

```
kat plot spectra-cn -o scer_pe_v2_ctgs.png -x 100 scer_pe_v2_ctgs.mx
```

This has generated a spectra-cn plot, showing us the distribution of content in our reads and how it appears in our assembly. Take a look the the .png produced. Because we first used default parameters for our assembler, we have caused some problems with our assembly. Notice how there is a lot of black in the centre of our main distribution - this is content from the PE reads that does not appear in our assembly. We can also see a lot of duplicated data (orange, green) on the edge of our main distribution. If we continued with this set of contigs, our assembly would be poor quality and would not represent all the data we obtained from the PE reads. Adjusting your parameters to your specific needs are important! This is a good example of why we shouldn't just use the default parameters on a tool.

Now try running the spectra-cn again:

```
kat comp -o scer_pe_v2_ctgs_new -t 8 -m 27 -H 100000000 -I 100000000 'pe_reads_R1.fastq pe_reads_R2.fastq' new_assm_dir/contigs.fa
```

```
kat plot spectra-cn -o scer_pe_v2_ctgs_new.png -x 100 scer_pe_v2_ctgs_new.mx
```

This looks better! You can see we have a k-mer distribution peaking at around 60, with most of our content appearing once and very little of our content missing. This means that most of the k-mer content within the PE reads is also in our assembly, and only appears once - excellent! The peak at <5 shows erroneous k-mers, and the large separation between the erroneous k-mers and the distribution peak tells us that this is a good dataset.

6. Checking insert size and distribution

The final step in our simplest assembly is to check the insert size of the paired-end reads and distribution. This will tell us how many of our reads map back to our assembly and what the original DNA fragment sizes were. Hopefully, we will find the majority of them will as they are from the same origin, so there should be no errors. 

First, bwa needs to make an index of our assembly sequence:

```
bwa index -p assembly_index -a bwtsw new_assm_dir/contigs.fa
```

Here, "-p" specifies the prefix of the output and "-a bwtsw" is the algorithm for constructing the index. Next, we can align our reads to the reference, using the index we just created:

```
bwa mem -SP -t 8 assembly_index pe_reads_R1.fastq pe_reads_R2.fastq > pe2assembly.sam
```

Here, "-SP" specifies paired-end mode for our reads and "-t" the number of threads. We specify that we want our aligned read to be outputted to the .sam format with the last part of the command "> pe2assembly.sam".

Hopefully we will find that a reasonable percentage of our reads map to the assembly. We can find out how many reads are mapping by using samtools:

```
cat pe2assembly.sam | grep -v '^@' | cut -f 1 | sort | uniq | wc -l
```

This output gives us the number of mapped reads.

We can also use the SAM output to create the raw data we need to draw an insert size histogram: 


```
grep -v '@SQ' pe2assembly.sam | grep -v '@PG' | awk -v binsize=20 '{if ($5==60) {if ($9>0) {print int($9/binsize)}else{print int($9/binsize*-1)}}}' | sort -n | uniq -c | awk -v binsize=20 '{print $2*binsize","$1}' > pe2assembly.is
```

Then plot this with a plotting tool we created to view the histogram:

```
python plotting_toolPE.py <path to is file>
```

You can see from this graph that the insert sizes are very symmetrically distributed at ~350bp, which is what we expect. The distribution is fairly narrow because the insert sizes of the PE reads don't vary very much from the average, so should produce a decent assembly. 

7. Introducing some complexity - Duplications

Great! You've learnt the basics to quality check and assemble the contigs of a simple genome! Generally, assemblies are not as straight-forward as that one, however, so we will introduce some of the complexities you might find in a real dataset. 

In this next dataset we have introduced duplications. Duplications can occur biologically in a genome, whereby a section of the DNA is copied and reinserted, creating duplicate regions. These show up on our spectra-cn as an extra bump beside the main distribution. Let's take a look at how this affects our data. 

You'll need to change directory to genome_2, then take a look at the k-mer coverage of these reads:

```
kat hist -o 350bp_fragment/scer_pe_hist -h 80 -t 8 -m 27 -H 100000000 350bp_fragment/pe_reads_R1.fastq 350bp_fragment/pe_reads_R2.fastq
```

Again, we can see a fairly symmetrical distribution centred at ~60 with only a little variance. Let's assemble it:

```
/velvet/velveth 350bp_fragment/assm_dir 31 -fastq -shortPaired -separate 350bp_fragment/pe_reads_R1.fastq 350bp_fragment/pe_reads_R2.fastq
/velvet/velvetg 350bp_fragment/assm_dir
```

Take a look at our LastGraph in BANDAGE. Notice how much messier this graph is compared to our simple dataset? This graph isn't distinguishing between likely errors and the true path, offering multiple connections between nodes. The duplicates we are seeing in this dataset make this more confusing for the assembler, as there are multiple possible connections between repeat fragments.

Run the spectra-cn to compare the PE reads to the assembly:

```
kat comp -o 350bp_fragment/scer_pe_v2_ctgs_31noparam -t 8 -m 27 -H 100000000 -I 100000000 '350bp_fragment/pe_reads_R1.fastq 350bp_fragment/pe_reads_R2.fastq' 350bp_fragment/assm_dir/contigs.fa
kat plot spectra-cn -o scer_pe_v2_ctgs_31noparam.png -x 100 scer_pe_v2_ctgs_31noparam.mx
```

Whoa, that's a bad plot! Look at the amount of black missing sequence there is, and how large those duplications around the distribution edges are. We can see the extra duplications around ~120 but it's not very clear.

8. Changing parameters

We are experiencing greater complexity in this dataset due to the duplications, so we will need to change our assembly approach. One thing we can alter are the parameters in our assembly, adjusting the coverage cutoff value to include only the high-coverage fragments. This will prevent the assembler from misaligning the fragments and will have a noticeable effect on our graphs.

Let's try adjusting the parameters to see if we can make a better assembly. Look at the k-mer spectra histogram we created in step 7; here you can see the erroneous data at <5 with low coverage. You can also take a look at the LastGraph in BANDAGE and click on some of the nodes within webs of connections. You'll notice on the De Bruijn graph that most of the alternate nodes have a coverage of <2. So, we will adjust the coverage cutoff to 2, meaning any nodes with a coverage of less than 2 will be discounted. 

```
/velvet/velvetg 350bp_fragment/assm_dir -cov_cutoff 2
```

```
kat comp -o 350bp_fragment/scer_pe_v2_ctgs_31cutoff2 -t 8 -m 27 -H 100000000 -I 100000000 '350bp_fragment/pe_reads_R1.fastq 350bp_fragment/pe_reads_R2.fastq' 350bp_fragment/assm_dir/contigs.fa
kat plot spectra-cn -o output.png -x 100 input.mx
```

That's an improvement, the duplications from the edges of the distribution and black errors in the centre of the distribution have reduced. 

Let's check the De Bruijn graph. This looks quite fragmented; we don't see the same tangled mess as we did before, but our nodes are not connected up into the four chromosomes we would expect. The coverage cutoff has resolved some of the ambiguities, preventing the misidentified paths we saw before. However, our k-mer size is too small, so the fragments are not overlapping enough for the algorithm to place them in sequential order.

9. Changing the k-mer

Our k-mer value in the last run was 31, which was too small; too many of the nodes look like repeats with a small k-mer size. A low k-mer value gives greater coverage, but at the expense of the k-mer's ability to overlap and form longer strings (ie chromosomes). Let's try changing this to 61 and reassembling. We will only adjust the one parameter at a time so that we can easily tell the effect it is having on our data.

```
/velvet/velveth 350bp_fragment/assm_dir 61 -fastq -shortPaired -separate 350bp_fragment/pe_reads_R1.fastq 350bp_fragment/pe_reads_R2.fastq
/velvet/velvetg 350bp_fragment/assm_dir
```

```
kat comp -o 350bp_fragment/scer_pe_v2_ctgs_61noparams -t 8 -m 27 -H 100000000 -I 100000000 '350bp_fragment/pe_reads_R1.fastq 350bp_fragment/pe_reads_R2.fastq' 350bp_fragment/assm_dir/contigs.fa
kat plot spectra-cn -o output.png -x 100 input.mx
```


The De Bruijn graph looks much better. We can see our four distinct chromosomes, but it can still be improved by adding a coverage cutoff. We're also back to having a high error rate in the spectra-cn plot. With such a large black section, there are many k-mers in our PE reads that have not been included in the assembly.

10. Combining cutoff and k-mer alterations

So if we combine the k-mer value at 61 with a coverage cutoff of 2…

```
/velvet/velvetg 350bp_fragment/assm_dir -cov_cutoff 2
```

```
kat comp -o 350bp_fragment/scer_pe_v2_ctgs_61cutoff2 -t 8 -m 27 -H 100000000 -I 100000000 '350bp_fragment/pe_reads_R1.fastq 350bp_fragment/pe_reads_R2.fastq' 350bp_fragment/assm_dir/contigs.fa
kat plot spectra-cn -o output.png -x 100 input.mx
```

That's much better! There is very little missing data here, and our plot shows no duplicated k-mers around the distribution edges. What we can now clearly see is that we have some duplicated data centred at ~120. This is normal in a real dataset due to biological duplications within the genome. We can also see a tiny bump of heterozygosity at ~25. If you check the De Bruijn graph, you can see some heterozygous content forming bubbles in the chromosome. Click on the individual nodes and it will tell you the node coverage, or depth, on the right hand side. Some of them are more likely to be true heterozygosity than others - try and identify which nodes may be errors (low coverage). 

Glossary:

```
High-throughput Sequencing (HTS) - the process by which genomic material can be translated into a string of sequential nucleotides. This is carried out in a variety of methods by sequencing machines, or "sequencers". 

Nucleotide - the four bases that make up a strand of DNA, adenine (A), guanine (G), cytosine (C) and thymine (T).

Fragments - small portions of DNA sequenced by HTS. HTS systems cannot sequence an entire genome in one string, so smaller fragments are made. These can be short-read (~100-600 base pairs) or long-read (~10 - 15 kilo-base pairs).

Inserts - a section of DNA that is inserted into a vector, such as E. coli or yeast, for multiplication or manipulation. 

Library (Insert Library) - a collection of inserts that make up an entire genome.

Reads - a set of reads is produced from an insert library. Reads are a computational sequence of base pairs corresponding to the inserts. 

Paired-end (PE) reads - reads that are read once forwards, and once backwards. This method is used to increase accuracy of the mapping reads, especially over repetitive regions. 

Long Mate Pair (LMP) reads - similar to PE reads, but a longer sequence is covered. Used in conjunction with PE reads, it can provide greater genome coverage. 

Raw reads - the genomic sequence reads straight from the sequencer. No modifications have been made and the reads inevitably will have errors. 

Genome Assembly - the computational process of overlapping, and subsequently, piecing together of all the inserts of a library to make up an entire genome.

Contigs (Contiguous Sequences) - portions of a partially-assembled genome produced from overlapping inserts. The creation of contigs forms the first level in genome assembly.

Scaffolds - the joining together of contigs, separated by gaps of known length. Scaffolding forms the second level of genome assembly.

Computational tools (or just 'tools') - programmes on a computer that perform a specific task needed by the user. These involve the use of algorithms to perform tasks such as calculating statistics, aligning sequences or constructing assemblies.

Pipeline - a set of computational tools that work in tandem to fulfil an overreaching goal through a step-by-step approach. This can be scripted to run automatically or can be run manually by the user. 

Bias - the over- or under-representation of a measurement in statistics. 

Coverage (depth) - the number of times any given nucleotide in any given fragment of DNA is read (or covered). The higher the coverage, the more reliable the base call can be considered. 

Insertions - the addition of nucleotides into a sequence. This can occur biologically or as an error in sequencing. 

Deletions - the loss of nucleotides in a sequence. This can occur biologically or as an error in sequencing. 

Indels - Insertion/deletion events.

K-mer - a string of 'k' length from a sequence. 

K-mer spectra - a plot that provides information to analyze how much and what type of k-mer content from reads is present in an assembly. It decomposes the k-mer spectrum of a read data set by the frequency in which the k-mers are encountered in the assembly.

Variance - the spread of values around the average point in a set of data. 

Error - can refer to miscalled nucleotides or stretches of nucleotides in a sequence, or the amount by which an observation deviates from its expected values.

```