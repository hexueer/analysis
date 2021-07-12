# Some code
That will calculate the mean and standard deviation of the lengths and gc content of your coronavirus genomes

##### Get genomes
```julia
using BioinformaticsBISC195

genomeData = parse_fasta("../data/selectedSARSgenomes.fasta")
headerVec = genomeData[1]
sequenceVec = genomeData[2]
```

##### Calculate the minimum, maximum, mean, and standard deviation of the sequence length of the genomes
##### I also plotted a histogram of the lengths :D

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

##### Filter out data with sequence length <25K
```julia
# sort and cut off at above 25K mark
# lenPermVec = sortperm(seqLengths)
# orderedSeqLengths = seqLengths[lenPermVec] # check ordered!
# headerOrderedByLen = headerVec[lenPermVec]
# seqOrderedByLen = sequenceVec[lenPermVec]
# len, index = 0
# while len < 25000
#     index += 1
#     len = orderedSeqLengths[index]
# end
# headerVec = headerOrderedByLen[index:end]
# sequenceVec = seqOrderedByLen[index:end]

# findall less than 25K and delete
lengthLessThan25KInd = findall(x -> x<=25000, seqLengths)
count = 0
for ind in lengthLessThan25KInd
    deleteat!(headerVec, ind-count)
    deleteat!(sequenceVec, ind-count)
    count += 1
end
@assert length(headerVec) == length(sequenceVec) # check!

# plot again!
histogram(seqLengths, legend = false, xlabel = "Sequence Length (above 25K)", ylabel = "Number of Genomes", title = "Coronavirus Genomes By Sequence Length")
```

##### Calculate the mean and standard deviation of GC content of the genomes

```julia
seqGC = map(gc_content, genomes)

mean_cov_gc = mean(seqGC)
std_cov_gc = std(seqGC)
```