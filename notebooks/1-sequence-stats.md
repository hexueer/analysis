# Some code
That will calculate the mean and standard deviation of the lengths and gc content of your coronavirus genomes

##### Get genomes
```julia
using BioinformaticsBISC195

genomes = parse_fasta("../data/selectedSARSgenomes.fasta")[2]
```

##### Calculate the mean sequence length and standard deviation of the genomes

```julia
using Statistics

seqLengths = map(length, genomes)

mean_cov2_length = mean(seqLengths)
std_cov2_length = std(seqLengths)
```

##### Calculate the mean and standard deviation of GC content of the genomes

```julia
seqGC = map(gc_content, genomes)

mean_cov2_gc = mean(seqGC)
std_cov2_gc = std(seqGC)
```