# GC Content Analysis
In this notebook, we will be analyzing the genomes found in *.../data/selectedSARSgenomes.fasta* through the lens of gc content across time.

First, let's grab our data! In the Kmer Analysis, we parsed and saved all the data in files for easy loading.

```julia
using JLD2, FileIO

headerVec = load("../data/genomeData.jld2", "headerVec")
sequenceVec = load("../data/genomeData.jld2", "sequenceVec")
```

Next, let's create a corresponding vector that contains the gc content of each sequence.

```julia
using BioinformaticsBISC195

gcContentVec = map(gc_content, sequenceVec)
```
It would also be nice to have vectors with the accession number and collection date;

```julia
accessionVec, collectDateVec = [], []
for header in headerVec
    headerInfo = getHeaderAttrib(header, "|", [1, 5])
    push!(accessionVec, headerInfo[1])
    push!(collectDateVec, headerInfo[2])
end
```

Now, let's make a dataframe! Perhaps something that looks like this?
|Accession|GC Content|Collection Date|
|---|---|---|
|id 1|gc_content 1|date 1|
|id 2|gc_content 2|date 2|
|id 3|gc_content 3|date 3|

```julia
using DataFrames

gcContentDf = DataFrame(id = accessionVec, gc_content = gcContentVec, collection_date = collectDateVec)
```

Haha, as you can see, I am no where near done! Instead you can look at my cats:
![Cat #1](/assets/cat#1.jpg)
![Cat #2](/assets/cat#2.jpg)

Now that we have a table to organize our data, let's plot it for easy viewing!

```julia
using StatsPlots
```