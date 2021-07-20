# Some code 
That will calculate the mean and standard deviation of the lengths and gc content of the coronavirus genomes found in .../data/selectedSARSgenomes.fasta.

#### Get genomes
First, we will make use of our BioinformaticsBISC195 package which has a function to parse through fasta files. This function call will give us vector to hold the headers as well as a vector to hold sequences.

```julia
using BioinformaticsBISC195

genomeData = parse_fasta("../data/selectedSARSgenomes.fasta")
headerVec = genomeData[1]
sequenceVec = genomeData[2]
```

#### Length Analysis
Next, we will calculate the minimum, maximum, mean, and standard deviation of the sequence length of the genomes. This will give us an idea of what the sequence data looks like in terms of length. To literally visualize the data, I also plotted a histogram of the lengths :D

```julia
using Statistics
using Plots

seqLengths = map(length, sequenceVec)

min__cov_length = minimum(seqLengths) # 23327
max__cov_length = maximum(seqLengths) # 30484
mean_cov_length = mean(seqLengths) # 29837.373197326768
std_cov_length = std(seqLengths) # 168.65155310071873
histogram(seqLengths, legend = false, xlabel = "Sequence Length", ylabel = "Number of Genomes", title = "Coronavirus Genomes By Sequence Length") # no outliers under 1000, attributes modified
```

#### Filter out data with sequence length <25K
We have two methods we can attempt:

1. Sort the vectors and cut off right above the 25K mark
```julia
lenPermVec = sortperm(seqLengths)
orderedSeqLengths = seqLengths[lenPermVec] # check ordered!
headerOrderedByLen = headerVec[lenPermVec]
seqOrderedByLen = sequenceVec[lenPermVec]
len, index = 0
while len < 25000
    index += 1
    len = orderedSeqLengths[index]
end
headerVec = headerOrderedByLen[index:end]
sequenceVec = seqOrderedByLen[index:end]
```

2. Findall sequences less than 25K and delete them
```julia
lengthLessThan25KInd = findall(x -> x<=25000, seqLengths)
count = 0
for ind in lengthLessThan25KInd
    deleteat!(headerVec, ind-count)
    deleteat!(sequenceVec, ind-count)
    count += 1
end
@assert length(headerVec) == length(sequenceVec) # check!
```

Let's re-plot our data now that we've trimmed off some outliers!
```julia
histogram(seqLengths, legend = false, xlabel = "Sequence Length (above 25K)", ylabel = "Number of Genomes", title = "Coronavirus Genomes By Sequence Length")
```

#### GC Content Analysis
The last thing we will do in this notebook is calculate the mean and standard deviation of GC content of the genomes.

```julia
seqGC = map(gc_content, genomes)

mean_cov_gc = mean(seqGC)
std_cov_gc = std(seqGC)
```