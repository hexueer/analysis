# Kmer Analysis
In this notebook, we will be analyzing the genomes found in *.../data/selectedSARSgenomes.fasta* through the lens of kmers across geograpghical region.

Let's start by parsing through the fasta file to retrieve the sequence data!

```julia
using BioinformaticsBISC195

genomeData = parse_fasta("../data/selectedSARSgenomes.fasta")
```

For my own convenience and the sanctity of my laptop, let's write all this data into files!

```julia
using JLD2, FileIO

save("../data/genomeData.jld2", Dict("headerVec" => genomeData[1], "sequenceVec" => genomeData[2]))
```

Next, we will create a new vector that contains the corresponding kmer data to each sequence. Each vector item will be a dictionary that contains the counts of every unique kmer of length *k* that is contained within that specific sequence. Then, we will save that data in a file!

```julia
# kmers of length 3, 4, 5, 6, 7
k3merVec, k4merVec, k5merVec, k6merVec, k7merVec = [], [], [], [], []
for sequence in sequenceVec
    push!(k3merVec, getKmerCount(sequence, 3))
    push!(k4merVec, getKmerCount(sequence, 4))
    push!(k5merVec, getKmerCount(sequence, 5))
    push!(k6merVec, getKmerCount(sequence, 6))
    push!(k7merVec, getKmerCount(sequence, 7))
end

save("../data/kmerData.jld2", Dict("3mer" => k3merVec, "4mer" => k4merVec, "5" => k5merVec, "6mer" => k6merVec, "7mer" => k7merVec))
```
Now that we have collected all the data into files, we can just load it whenever we need to continue working.

```julia
headerVec = load("../data/genomeData.jld2", "headerVec")
sequenceVec = load("../data/genomeData.jld2", "sequenceVec")
k3merVec = load("../data/kmerData.jld2", "3mer")
k4merVec = load("../data/kmerData.jld2", "4mer")
k5merVec = load("../data/kmerData.jld2", "5")
k6merVec = load("../data/kmerData.jld2", "6mer")
k7merVec = load("../data/kmerData.jld2", "7mer")
```

Since we have all the data readily available now, let's build a table to organize it! We probably want something like this:
|Unique Kmer|*k* Length|Count|Accession|Location|
|------|------|------|------|------|
|kmer 1|length 1|count 1|id 1|location 1|
|kmer 2|length 2|count 2|id 2|location 2|
|kmer 3|length 3|count 3|id 3|location 3|

We'll use DataFrames.jl to do this! The easiest way to set up a dataframe is to pass a vector for each column, so first, we have to build our vectors.

```julia
kmer, k, count, id, location = [], [], [], [], []
for num in 3:6
    num == 3 && (kmerVec = k3merVec)
    num == 4 && (kmerVec = k4merVec)
    num == 5 && (kmerVec = k5merVec)
    num == 6 && (kmerVec = k6merVec)
    for index in 1:length(sequenceVec)
        for kmerKey in keys(kmerVec[index])
            push!(kmer, kmerKey)
            push!(k, num)
            push!(count, kmerVec[index][kmerKey])
            headerInfo = getHeaderAttrib(headerVec[index], "|", [1,4])
            push!(id, headerInfo[1])
            push!(location, headerInfo[2])
        end
    end
end

save("../data/kmerDf.jld2", Dict("kmer" => kmer, "k" => k, "count" => count, "id" => id, "location" => location)) # in case we want to skip all the previous steps
```

Now we have our vectors! You might have noticed that I excluded the 7mers; that's because my laptop was struggling! So I gave up on that dream. Moving on, let's set up our dataframe and see what it looks like. So exciting!!

```julia
using DataFrames

kmerDf = DataFrame(kmer = kmer, k = k, count = count, id = id, location = location)
```
![kmerDf](/assets/kmerDf.png)

As you can tell, the location column is taken straight from the header so it's not consistent. Some of them contain city names while some only contain country names. Let's take a look at all the unique locations by turning the vector into a set.

```julia
Set(kmerDf.location)
```
![kmerLocationSet](/assets/kmerLocationSet.png)

Now that we have an idea of what the locations are, and how they are formatted (country name before first colon, if any), let's add a new column to generalize by country.

```julia
country = map(getCountryFromLocation, kmerDf.location)
kmerDf.country = country
Set(kmerDf.country) # see unique countries
```
![kmerDfCountry](/assets/kmerDfCountry.png)

We have narrowed down from 240 locations to 60 countries! Now, let's add just one more column to organize them by continent.

```julia
using Countries
# blah blah code
```

Now that we have a lovely table, let's plot our data!

```julia
using Plots

```