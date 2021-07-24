# Kmer Analysis
In this notebook, we will be analyzing the genomes found in *.../data/selectedSARSgenomes.fasta* through the lens of kmers across geograpghical region.

Let's start by parsing through the fasta file to retrieve the sequence data!

```julia
julia> using BioinformaticsBISC195

julia> genomeData = parse_fasta("../data/selectedSARSgenomes.fasta")
```

For my own convenience and the sanctity of my laptop, let's write all this data into files!

```julia
julia> using JLD2, FileIO

julia> save("../data/genomeData.jld2", Dict("headerVec" => genomeData[1], "sequenceVec" => genomeData[2]))
julia> headerVec = load("../data/genomeData.jld2", "headerVec")
5774-element Vector{Any}:
 "NC_045512.2 |Severe acute respiratory syndrome-related coronavirus||China|2019-12"
 "MZ573077.1 |Severe acute respiratory syndrome-related coronavirus||Bahrain|2021-07-05"
 "MZ573079.1 |Severe acute respiratory syndrome-related coronavirus||Bahrain|2021-07-05"
 "MZ573080.1 |Severe acute respiratory syndrome-related coronavirus||Bahrain|2021-07-05"
 "MZ572200.1 |Severe acute respiratory syndrome-related coronavirus||India|2021-06-01"
 "MZ572201.1 |Severe acute respiratory syndrome-related coronavirus||India|2021-06-01"
 "MZ572203.1 |Severe acute respiratory syndrome-related coronavirus||India|2021-06-01"
 "MZ572204.1 |Severe acute respiratory syndrome-related coronavirus||India|2021-06-04"
 "MZ572206.1 |Severe acute respiratory syndrome-related coronavirus||India|2021-06-09"
 "MZ572207.1 |Severe acute respiratory syndrome-related coronavirus||India|2021-06-14"
 "MZ571142.1 |Severe acute respiratory syndrome-related coronavirus|oronasopharynx|Morocco|2021-04-22"
 "MZ562707.1 |Severe acute respiratory syndrome-related coronavirus|oronasopharynx|Pakistan|2021-04-30"
 "MZ562746.1 |Severe acute respiratory syndrome-related coronavirus||India: Madhya Pradesh|2021-03-07"
 "MZ562747.1 |Severe acute respiratory syndrome-related coronavirus||India: Madhya Pradesh|2021-03-07"
 "MZ562748.1 |Severe acute respiratory syndrome-related coronavirus||India: Madhya Pradesh|2021-03-08"
 "MZ562749.1 |Severe acute respiratory syndrome-related coronavirus||India: Madhya Pradesh|2021-03-08"
 "MZ562750.1 |Severe acute respiratory syndrome-related coronavirus||India: Madhya Pradesh|2021-03-08"
 "MZ562751.1 |Severe acute respiratory syndrome-related coronavirus||India: Madhya Pradesh|2021-03-07"
 "MZ562752.1 |Severe acute respiratory syndrome-related coronavirus||India: Madhya Pradesh|2021-02-08"
 "MZ562753.1 |Severe acute respiratory syndrome-related coronavirus||India: Madhya Pradesh|2021-02-08"
 ⋮
 "AY502925.1 |Severe acute respiratory syndrome-related coronavirus||Taiwan|"
 "AY502926.1 |Severe acute respiratory syndrome-related coronavirus||Taiwan|"
 "AY502927.1 |Severe acute respiratory syndrome-related coronavirus||Taiwan|"
 "AY502928.1 |Severe acute respiratory syndrome-related coronavirus||Taiwan|"
 "AY502929.1 |Severe acute respiratory syndrome-related coronavirus||Taiwan|"
 "AY502930.1 |Severe acute respiratory syndrome-related coronavirus||Taiwan|"
 "AY502931.1 |Severe acute respiratory syndrome-related coronavirus||Taiwan|"
 "AY502932.1 |Severe acute respiratory syndrome-related coronavirus||Taiwan|"
 "AY313906.1 |Severe acute respiratory syndrome-related coronavirus||China: Jiangmen, Guangdong|"
 "AY304486.1 |Severe acute respiratory syndrome-related coronavirus||Hong Kong|"
 "AY304488.1 |Severe acute respiratory syndrome-related coronavirus||Hong Kong|"
 "AY304495.1 |Severe acute respiratory syndrome-related coronavirus||Hong Kong|"
 "AY348314.1 |Severe acute respiratory syndrome-related coronavirus||Taiwan|"
 "AY291451.1 |Severe acute respiratory syndrome-related coronavirus||Taiwan|"
 "AY283794.1 |Severe acute respiratory syndrome-related coronavirus||Singapore|"
 "AY283795.1 |Severe acute respiratory syndrome-related coronavirus||Singapore|"
 "AY283796.1 |Severe acute respiratory syndrome-related coronavirus||Singapore|"
 "AY283797.1 |Severe acute respiratory syndrome-related coronavirus||Singapore|"
 "AY283798.2 |Severe acute respiratory syndrome-related coronavirus||Singapore|"
 "AY278554.2 |Severe acute respiratory syndrome-related coronavirus||China: Hong Kong, Prince of Wales Hospital|"

julia> sequenceVec = load("../data/genomeData.jld2", "sequenceVec")
# output not shown here because it is too long and honestly not helpful to look at
```

Next, we will create a new vector that contains the corresponding kmer data to each sequence. Each vector item will be a dictionary that contains the counts of every unique kmer of length *k* that is contained within that specific sequence. Then, we will save that data in a file!

```julia
# kmers of length 3, 4, 5, 6, 7
julia> k3merVec, k4merVec, k5merVec, k6merVec, k7merVec = [], [], [], [], []
julia> for sequence in sequenceVec
    push!(k3merVec, getKmerCount(sequence, 3))
    push!(k4merVec, getKmerCount(sequence, 4))
    push!(k5merVec, getKmerCount(sequence, 5))
    push!(k6merVec, getKmerCount(sequence, 6))
    push!(k7merVec, getKmerCount(sequence, 7))
end

# save to file
julia> save("../data/kmerData.jld2", Dict("3mer" => k3merVec, "4mer" => k4merVec, "5" => k5merVec, "6mer" => k6merVec, "7mer" => k7merVec))

# load when you want
julia> k3merVec = load("../data/kmerData.jld2", "3mer")
julia> k4merVec = load("../data/kmerData.jld2", "4mer")
julia> k5merVec = load("../data/kmerData.jld2", "5")
julia> k6merVec = load("../data/kmerData.jld2", "6mer")
julia> k7merVec = load("../data/kmerData.jld2", "7mer")
```

Right now, our k#merVecs look like this:
```julia
julia> k3merVec
5686-element Vector{Any}:
 Dict{Any, Any}("TGT" => 858, "GAC" => 340, "GAA" => 535, "TTC" => 518, "ACA" => 809, "TAG" => 427, "GTA" => 469, "GTG" => 552, "CCT" => 344, "GCT" => 521…)
 Dict{Any, Any}("TGT" => 859, "GAC" => 340, "GAA" => 534, "TTC" => 516, "ACA" => 810, "TAG" => 426, "GTA" => 470, "GTG" => 552, "CCT" => 343, "GCT" => 521…)
 Dict{Any, Any}("TGT" => 859, "GAC" => 340, "GAA" => 535, "TTC" => 518, "ACA" => 807, "TAG" => 426, "GTA" => 469, "GTG" => 553, "CCT" => 342, "GCT" => 520…)
 Dict{Any, Any}("TGT" => 859, "GAC" => 340, "GAA" => 535, "TTC" => 516, "ACA" => 810, "TAG" => 426, "GTA" => 469, "GTG" => 553, "CCT" => 343, "GCT" => 520…)
 Dict{Any, Any}("TGT" => 859, "GAC" => 340, "GAA" => 535, "TTC" => 516, "ACA" => 808, "TAG" => 426, "GTA" => 469, "GTG" => 553, "CCT" => 343, "GCT" => 521…)
 Dict{Any, Any}("TGT" => 859, "GAC" => 340, "GAA" => 533, "TTC" => 517, "ACA" => 808, "TAG" => 427, "GTA" => 469, "GTG" => 553, "CCT" => 343, "GCT" => 520…)
 Dict{Any, Any}("TGT" => 857, "GAC" => 338, "GAA" => 534, "TTC" => 513, "ACA" => 804, "TAG" => 428, "GTA" => 468, "GTG" => 552, "GCT" => 519, "CCT" => 344…)
 Dict{Any, Any}("TGT" => 857, "GAC" => 338, "GAA" => 534, "TTC" => 513, "ACA" => 804, "TAG" => 428, "GTA" => 468, "GTG" => 552, "GCT" => 519, "CCT" => 344…)
 Dict{Any, Any}("TGT" => 857, "GAC" => 338, "GAA" => 534, "TTC" => 513, "ACA" => 804, "TAG" => 428, "GTA" => 468, "GTG" => 552, "GCT" => 519, "CCT" => 344…)
 Dict{Any, Any}("TGT" => 857, "GAC" => 338, "GAA" => 534, "TTC" => 513, "ACA" => 804, "TAG" => 428, "GTA" => 468, "GTG" => 552, "GCT" => 519, "CCT" => 344…)
 Dict{Any, Any}("TGT" => 857, "GAC" => 338, "GAA" => 534, "TTC" => 513, "ACA" => 804, "TAG" => 428, "GTA" => 468, "GTG" => 552, "GCT" => 519, "CCT" => 344…)
 Dict{Any, Any}("TGT" => 857, "GAC" => 338, "GAA" => 534, "TTC" => 513, "ACA" => 804, "TAG" => 428, "GTA" => 468, "GTG" => 552, "GCT" => 519, "CCT" => 344…)
 Dict{Any, Any}("TGT" => 857, "GAC" => 338, "GAA" => 534, "TTC" => 513, "ACA" => 804, "TAG" => 428, "GTA" => 468, "GTG" => 552, "GCT" => 519, "CCT" => 344…)
 Dict{Any, Any}("TGT" => 857, "GAC" => 338, "GAA" => 534, "TTC" => 513, "ACA" => 804, "TAG" => 428, "GTA" => 468, "GTG" => 552, "GCT" => 519, "CCT" => 344…)
 Dict{Any, Any}("TGT" => 858, "GAC" => 339, "GAA" => 533, "TTC" => 512, "ACA" => 806, "TAG" => 423, "GTA" => 465, "GTG" => 553, "GCT" => 520, "CCT" => 342…)
 Dict{Any, Any}("TGT" => 856, "GAC" => 338, "GAA" => 534, "TTC" => 510, "ACA" => 804, "TAG" => 423, "GTA" => 467, "GTG" => 554, "GCT" => 520, "CCT" => 343…)
 Dict{Any, Any}("TGT" => 855, "GAC" => 340, "GAA" => 534, "TTC" => 509, "ACA" => 806, "TAG" => 422, "GTA" => 466, "GTG" => 553, "GCT" => 518, "CCT" => 342…)
 Dict{Any, Any}("TGT" => 856, "GAC" => 338, "GAA" => 534, "TTC" => 510, "ACA" => 805, "TAG" => 423, "GTA" => 467, "GTG" => 554, "GCT" => 520, "CCT" => 343…)
 Dict{Any, Any}("TGT" => 856, "GAC" => 338, "GAA" => 534, "TTC" => 510, "ACA" => 804, "TAG" => 423, "GTA" => 467, "GTG" => 554, "GCT" => 520, "CCT" => 343…)
 Dict{Any, Any}("TGT" => 856, "GAC" => 338, "GAA" => 534, "TTC" => 510, "ACA" => 804, "TAG" => 423, "GTA" => 467, "GTG" => 554, "GCT" => 520, "CCT" => 343…)
 ⋮
 Dict{Any, Any}("TGT" => 795, "GAC" => 367, "GAA" => 454, "TTC" => 563, "ACA" => 780, "TAG" => 351, "GTA" => 438, "GTG" => 550, "CCT" => 346, "GCT" => 622…)
 Dict{Any, Any}("TGT" => 795, "GAC" => 366, "GAA" => 454, "TTC" => 563, "ACA" => 780, "TAG" => 351, "TCN" => 1, "GTA" => 438, "NGC" => 1, "GTG" => 550…)
 Dict{Any, Any}("TGT" => 795, "GAC" => 366, "GAA" => 454, "TTC" => 563, "ACA" => 780, "TAG" => 351, "GTA" => 438, "GTG" => 550, "CCT" => 346, "GCT" => 621…)
 Dict{Any, Any}("TGT" => 795, "GAC" => 367, "GAA" => 454, "TTC" => 562, "ACA" => 780, "TAG" => 351, "GTA" => 438, "GTG" => 550, "CCT" => 346, "GCT" => 622…)
 Dict{Any, Any}("TGT" => 794, "GAC" => 367, "GAA" => 454, "TTC" => 562, "ACA" => 780, "TAG" => 351, "GTA" => 438, "GTG" => 550, "CCT" => 346, "GCT" => 622…)
 Dict{Any, Any}("TGT" => 795, "GAC" => 367, "GAA" => 454, "TTC" => 562, "ACA" => 780, "TAG" => 351, "GTA" => 438, "GTG" => 551, "CCT" => 346, "GCT" => 622…)
 Dict{Any, Any}("TGT" => 794, "GAC" => 367, "GAA" => 454, "TTC" => 562, "ACA" => 780, "TAG" => 351, "TCN" => 1, "GTA" => 438, "GTG" => 550, "CCT" => 345…)
 Dict{Any, Any}("TGT" => 794, "GAC" => 367, "GAA" => 454, "TTC" => 562, "ACA" => 780, "TAG" => 351, "GTA" => 438, "GTG" => 550, "CCT" => 346, "GCT" => 622…)
 Dict{Any, Any}("TGT" => 795, "GAC" => 368, "GAA" => 454, "TTC" => 560, "ACA" => 780, "TAG" => 351, "GTA" => 440, "GTG" => 551, "CCT" => 346, "GCT" => 622…)
 Dict{Any, Any}("TGT" => 789, "GAC" => 366, "GAA" => 461, "TTC" => 563, "ACA" => 777, "TAG" => 350, "GTA" => 436, "GTG" => 549, "CCT" => 349, "GCT" => 622…)
 Dict{Any, Any}("TGT" => 789, "GAC" => 365, "GAA" => 459, "TTC" => 563, "ACA" => 775, "TAG" => 349, "GTA" => 436, "GTG" => 549, "CCT" => 349, "GCT" => 622…)
 Dict{Any, Any}("TGT" => 795, "GAC" => 367, "GAA" => 455, "TTC" => 563, "ACA" => 781, "TAG" => 351, "GTA" => 438, "GTG" => 550, "CCT" => 345, "GCT" => 622…)
 Dict{Any, Any}("TGT" => 789, "GAC" => 366, "GAA" => 449, "TTC" => 560, "ACA" => 779, "TAG" => 345, "GTA" => 435, "GTG" => 547, "GCT" => 621, "CCT" => 342…)
 Dict{Any, Any}("TGT" => 795, "GAC" => 367, "GAA" => 454, "TTC" => 563, "ACA" => 780, "TAG" => 351, "GTA" => 438, "GTG" => 550, "CCT" => 346, "GCT" => 622…)
 Dict{Any, Any}("TGT" => 795, "GAC" => 367, "GAA" => 455, "TTC" => 563, "ACA" => 778, "TAG" => 350, "GTA" => 438, "GTG" => 550, "CCT" => 345, "GCT" => 622…)
 Dict{Any, Any}("TGT" => 795, "GAC" => 367, "GAA" => 453, "TTC" => 562, "ACA" => 778, "TAG" => 350, "GTA" => 438, "GTG" => 551, "CCT" => 345, "GCT" => 623…)
 Dict{Any, Any}("TGT" => 796, "GAC" => 367, "GAA" => 454, "TTC" => 563, "ACA" => 780, "TAG" => 350, "GTA" => 438, "GTG" => 550, "CCT" => 344, "GCT" => 622…)
 Dict{Any, Any}("TGT" => 795, "GAC" => 367, "NAA" => 1, "GAA" => 453, "TTC" => 562, "ACA" => 778, "TAG" => 350, "GTA" => 438, "GTG" => 550, "CCT" => 345…)
 Dict{Any, Any}("TGT" => 795, "GAC" => 367, "GAA" => 454, "TTC" => 563, "ACA" => 778, "TAG" => 351, "GTA" => 438, "GTG" => 550, "CCT" => 345, "GCT" => 622…)
 Dict{Any, Any}("TGT" => 794, "GAC" => 368, "GAA" => 454, "TTC" => 563, "ACA" => 781, "TAG" => 350, "GTA" => 437, "GTG" => 550, "CCT" => 346, "GCT" => 622…)
```

Since we have all the data readily available now, let's build a table to organize it! We probably want something like this:
|Unique Kmer|*k* Length|Count|Accession|Location|
|------|------|------|------|------|
|kmer 1|length 1|count 1|id 1|location 1|
|kmer 2|length 2|count 2|id 2|location 2|
|kmer 3|length 3|count 3|id 3|location 3|

We'll use DataFrames.jl to do this! The easiest way to set up a dataframe is to pass a vector for each column, so first, we have to build our vectors.

```julia
julia> kmer, k, count, id, location = [], [], [], [], []
julia> for num in 3:6
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

kmer, k, count, id, location = [], [], [], [], []
for index in 1:length(k3merVec)
    for kmerKey in keys(k3merVec[index])
        push!(kmer, kmerKey)
        push!(k, 3)
        push!(count, k3merVec[index][kmerKey])
        headerInfo = getHeaderAttrib(headerVec[index], "|", [1,4])
        push!(id, headerInfo[1])
        push!(location, headerInfo[2])
    end
end

julia> save("../data/kmerDf.jld2", Dict("kmer" => kmer, "k" => k, "count" => count, "id" => id, "location" => location)) # in case we want to skip all the previous steps

# load!
julia> kmer = load("../data/kmerDf.jld2", "kmer")
julia> k = load("../data/kmerDf.jld2", "k")
julia> count = load("../data/kmerDf.jld2", "count")
julia> id = load("../data/kmerDf.jld2", "id")
julia> location = load("../data/kmerDf.jld2", "location")
```

Now we have our vectors! You might have noticed that I excluded the 7mers; that's because my laptop was struggling! So I gave up on that dream. Moving on, let's set up our dataframe and see what it looks like. So exciting!!

```julia
julia> using DataFrames

julia> kmerDf = DataFrame(kmer = kmer, k = k, count = count, id = id, location = location)

julia> kmerDf
29619676×5 DataFrame
      Row │ kmer    k    count  id           location                          
          │ Any     Any  Any    Any          Any                               
──────────┼────────────────────────────────────────────────────────────────────
        1 │ TGT     3    858    NC_045512.2  China
        2 │ GAC     3    340    NC_045512.2  China
        3 │ GAA     3    535    NC_045512.2  China
        4 │ TTC     3    518    NC_045512.2  China
        5 │ ACA     3    809    NC_045512.2  China
        6 │ TAG     3    427    NC_045512.2  China
        7 │ GTA     3    469    NC_045512.2  China
        8 │ GTG     3    552    NC_045512.2  China
        9 │ CCT     3    344    NC_045512.2  China
       10 │ GCT     3    521    NC_045512.2  China
       11 │ GGC     3    223    NC_045512.2  China
       12 │ AGG     3    329    NC_045512.2  China
       13 │ CGG     3    76     NC_045512.2  China
       14 │ AGC     3    301    NC_045512.2  China
       15 │ ATT     3    773    NC_045512.2  China
       16 │ GAT     3    440    NC_045512.2  China
       17 │ CAG     3    438    NC_045512.2  China
       18 │ ATG     3    725    NC_045512.2  China
    ⋮     │   ⋮      ⋮     ⋮         ⋮                       ⋮
 29619659 │ TCAAGT  6    8      AY278554.2   China: Hong Kong, Prince of Wale…
 29619660 │ GAACAG  6    7      AY278554.2   China: Hong Kong, Prince of Wale…
 29619661 │ TACTGT  6    14     AY278554.2   China: Hong Kong, Prince of Wale…
 29619662 │ ATATTA  6    7      AY278554.2   China: Hong Kong, Prince of Wale…
 29619663 │ GCTACA  6    18     AY278554.2   China: Hong Kong, Prince of Wale…
 29619664 │ GCCAAT  6    6      AY278554.2   China: Hong Kong, Prince of Wale…
 29619665 │ GGATCA  6    3      AY278554.2   China: Hong Kong, Prince of Wale…
 29619666 │ CCCAGC  6    3      AY278554.2   China: Hong Kong, Prince of Wale…
 29619667 │ TTTAGT  6    13     AY278554.2   China: Hong Kong, Prince of Wale…
 29619668 │ CTGCGA  6    1      AY278554.2   China: Hong Kong, Prince of Wale…
 29619669 │ TCGTGC  6    6      AY278554.2   China: Hong Kong, Prince of Wale…
 29619670 │ TTAACG  6    3      AY278554.2   China: Hong Kong, Prince of Wale…
 29619671 │ ATTGTA  6    12     AY278554.2   China: Hong Kong, Prince of Wale…
 29619672 │ AAGTCG  6    2      AY278554.2   China: Hong Kong, Prince of Wale…
 29619673 │ TAAATC  6    5      AY278554.2   China: Hong Kong, Prince of Wale…
 29619674 │ GCTTGT  6    18     AY278554.2   China: Hong Kong, Prince of Wale…
 29619675 │ TTGTCC  6    5      AY278554.2   China: Hong Kong, Prince of Wale…
 29619676 │ CAGGTT  6    14     AY278554.2   China: Hong Kong, Prince of Wale…
                                                          29619640 rows omitted
```

As you can tell, the location column is taken straight from the header so it's not consistent. Some of them contain city names while some only contain country names. There also seems to be a lot slightly different ones, which would make it very difficult to plot the data as it is now. Let's take a look at all the unique locations by turning the vector into a set.

```julia
julia> uniqueLocations = Set(kmerDf.location)
Set{Any} with 240 elements:
  "Uruguay"
  "India: Dahod"
  "China: Anhui, Fuyang"
  "Philippines: NCR, Cavite City"
  "India: Junagadh"
  "Brazil: Amazonas, Manaus"
  "India: Khambhaliya"
  "India: Delhi"
  "India: Dhandhuka"
  "India: Jamnagar"
  "India: Gujarat, Mehsana"
  "India: Gujarat, Kodinar"
  "Bangladesh: Khulna"
  "India: Gujarat, Sutrapada"
  "Colombia: Antioquia"
  "Thailand"
  "Cambodia"
  "Philippines: Caloocan City"
  "India: Vadodara"
  "Japan: Osaka"
  "India: Dhansura"
  "India: Kadi"
  "Chile"
  "Pakistan: Gilgit Baltistan"
  "India: Bhuj"
  "Saudi Arabia: Jeddah"
  "Thailand: Bangkok"
  "Iraq"
  "India: Savli"
  "India: Botad"
  "Argentina"
  "China: Southern China"
  "India: Kapadvanj"
  "Tunisia: Bizerte"
  "India: Gondal"
  "Jordan: Amman"
  "Colombia"
  "Togo"
  "India: Anand"
  "Japan: Kanagawa, Sagamihara"
  ⋮ 
```

Now that we have an idea of what the locations are, and how they are formatted (country name before first colon, if any), let's add a new column to generalize by country. Of course, since we are automating this process, we are reproducing how the NCBI has labeled the data, even though there may be errors or we might not agree (e.g. Hong Kong labeled under China). Because of this, our analysis might miss some of the important nuances to the data.

```julia
julia> kmerDf.country = map(x -> split(x, ":")[1], kmerDf.location)

julia> uniqueCountries = Set(kmerDf.country) # see set of unique countries
Set{SubString{String}} with 60 elements:
  "Uruguay"
  "Thailand"
  "Cambodia"
  "Chile"
  "Iraq"
  "Argentina"
  "Colombia"
  "Togo"
  "Zambia"
  "India"
  "South Korea"
  "Malaysia"
  "Iran"
  "United Arab Emirates"
  "Israel"
  "Turkey"
  "Ethiopia"
  "Sierra Leone"
  "Venezuela"
  "Brazil"
  "Uzbekistan"
  "Ghana"
  "China"
  "Somalia"
  "Uganda"
  "Libya"
  "Qatar"
  "Egypt"
  "Tunisia"
  "Jordan"
  "Mali"
  "Pakistan"
  "Gambia"
  "Philippines"
  "Oman"
  "Kenya"
  "Singapore"
  "Peru"
  "Myanmar"
  "Ecuador"
  ⋮ 
```

We have narrowed down from 240 locations to 60 countries! That's great, but perhaps we should narrow it down a bit more... I found a JSON file *countries.json* populated with country data that was created by a GitHub user named lukes. From that mapping, I'll make a new JSON that fits my needs more closely.

```julia
julia> using JSON

# read from lukes JSON
julia> countryInfo = Dict()
julia> open("../data/countries.json", "r") do f
    global countryInfo
    info = read(f, String)
    countryInfo = JSON.parse(info)
end

# make my own dicts with just country, region, and sub-region info
julia> country_to_region = Dict()
julia> country_to_subregion = Dict()
julia> for country in countryInfo
    country_to_region[country["name"]] = country["region"]
    country_to_subregion[country["name"]] = country["sub-region"]
end

# DONT RUN THIS -- IT WILL NOT WORK ANYMOREEE, this writes to file but I also edited manually
julia> open("../data/country_to_region.json", "w") do f
        write(f, JSON.json(country_to_region))
end
julia> open("../data/country_to_subregion.json", "w") do f
        write(f, JSON.json(country_to_subregion))
end
```

The dictionary is obviously not perfect and does not perfectly match some of the country names used. I made edits manually work for my set of sequence data. Of course, feel free to edit as you need for your set of data.

```julia
# just run this to load
julia> country_to_region = Dict()
julia> country_to_subregion = Dict()

julia> open("../data/country_to_region.json", "r") do f
    global country_to_region
    info = read(f, String)
    country_to_region = JSON.parse(info)
end
julia> open("../data/country_to_subregion.json", "r") do f
    global country_to_subregion
    info = read(f, String)
    country_to_subregion = JSON.parse(info)
end
```

We've narrowed down the geographical categories by a lot! Here are the unique regions and sub-regions that came out of our data!

```julia
julia> Set(map(x -> country_to_region[x], collect(uniqueCountries)))
Set{Any} with 3 elements:
  "Africa"
  "Asia"
  "Americas"

julia> Set(map(x -> country_to_subregion[x], collect(uniqueCountries)))
Set{Any} with 8 elements:
  "Sub-Saharan Africa"
  "Northern Africa"
  "South-eastern Asia"
  "Latin America and the Caribbean"
  "Southern Asia"
  "Western Asia"
  "Central Asia"
  "Eastern Asia"
```

Now, let's add just one more column to the dataframe to organize the data by region and subregion. Also, it'd be nice if the types were more specific.

```julia
julia> kmerDf.region = map(x -> country_to_region[x], kmerDf.country)
julia> kmerDf.subregion = map(x -> country_to_subregion[x], kmerDf.country)

julia> kmerDf[!, :kmer] = convert.(String, kmerDf[:, :kmer])
julia> kmerDf[!, :k] = convert.(Int, kmerDf[:, :k])
julia> kmerDf[!, :count] = convert.(Int, kmerDf[:, :count])
julia> kmerDf
29619676×8 DataFrame
      Row │ kmer    k    count  id           location           country    region  subregion    
          │ String     Int64  Int64    Any          Any                SubStrin…  String  String       
──────────┼───────────────────────────────────────────────────────────────────────────────────────
        1 │ TGT     3    858    NC_045512.2  China              China      Asia    Eastern Asia
        2 │ GAC     3    340    NC_045512.2  China              China      Asia    Eastern Asia
        3 │ GAA     3    535    NC_045512.2  China              China      Asia    Eastern Asia
        4 │ TTC     3    518    NC_045512.2  China              China      Asia    Eastern Asia
        5 │ ACA     3    809    NC_045512.2  China              China      Asia    Eastern Asia
        6 │ TAG     3    427    NC_045512.2  China              China      Asia    Eastern Asia
        7 │ GTA     3    469    NC_045512.2  China              China      Asia    Eastern Asia
        8 │ GTG     3    552    NC_045512.2  China              China      Asia    Eastern Asia
        9 │ CCT     3    344    NC_045512.2  China              China      Asia    Eastern Asia
       10 │ GCT     3    521    NC_045512.2  China              China      Asia    Eastern Asia
       11 │ GGC     3    223    NC_045512.2  China              China      Asia    Eastern Asia
       12 │ AGG     3    329    NC_045512.2  China              China      Asia    Eastern Asia
       13 │ CGG     3    76     NC_045512.2  China              China      Asia    Eastern Asia
       14 │ AGC     3    301    NC_045512.2  China              China      Asia    Eastern Asia
       15 │ ATT     3    773    NC_045512.2  China              China      Asia    Eastern Asia
       16 │ GAT     3    440    NC_045512.2  China              China      Asia    Eastern Asia
       17 │ CAG     3    438    NC_045512.2  China              China      Asia    Eastern Asia
       18 │ ATG     3    725    NC_045512.2  China              China      Asia    Eastern Asia
    ⋮     │   ⋮      ⋮     ⋮         ⋮               ⋮               ⋮         ⋮          ⋮
 29619659 │ TCAAGT  6    8      AY278554.2   China: Hong Kong…  China      Asia    Eastern Asia
 29619660 │ GAACAG  6    7      AY278554.2   China: Hong Kong…  China      Asia    Eastern Asia
 29619661 │ TACTGT  6    14     AY278554.2   China: Hong Kong…  China      Asia    Eastern Asia
 29619662 │ ATATTA  6    7      AY278554.2   China: Hong Kong…  China      Asia    Eastern Asia
 29619663 │ GCTACA  6    18     AY278554.2   China: Hong Kong…  China      Asia    Eastern Asia
 29619664 │ GCCAAT  6    6      AY278554.2   China: Hong Kong…  China      Asia    Eastern Asia
 29619665 │ GGATCA  6    3      AY278554.2   China: Hong Kong…  China      Asia    Eastern Asia
 29619666 │ CCCAGC  6    3      AY278554.2   China: Hong Kong…  China      Asia    Eastern Asia
 29619667 │ TTTAGT  6    13     AY278554.2   China: Hong Kong…  China      Asia    Eastern Asia
 29619668 │ CTGCGA  6    1      AY278554.2   China: Hong Kong…  China      Asia    Eastern Asia
 29619669 │ TCGTGC  6    6      AY278554.2   China: Hong Kong…  China      Asia    Eastern Asia
 29619670 │ TTAACG  6    3      AY278554.2   China: Hong Kong…  China      Asia    Eastern Asia
 29619671 │ ATTGTA  6    12     AY278554.2   China: Hong Kong…  China      Asia    Eastern Asia
 29619672 │ AAGTCG  6    2      AY278554.2   China: Hong Kong…  China      Asia    Eastern Asia
 29619673 │ TAAATC  6    5      AY278554.2   China: Hong Kong…  China      Asia    Eastern Asia
 29619674 │ GCTTGT  6    18     AY278554.2   China: Hong Kong…  China      Asia    Eastern Asia
 29619675 │ TTGTCC  6    5      AY278554.2   China: Hong Kong…  China      Asia    Eastern Asia
 29619676 │ CAGGTT  6    14     AY278554.2   China: Hong Kong…  China      Asia    Eastern Asia
                                                                            29619640 rows omitted
```

The biggest issue now is that we have multiple entries for the same kmer, since added each kmer for each sequence into the dataframe. That's okay, but now we want to be selective about what we will plot, since we can't plot everything! Let's pick the top 50 most frequent ones!

```julia
kmerTotalCounts = Dict()
for kStr in Set(kmerDf.kmer)
    kmerTotalCounts[kStr] = sum(kmerDf[occursin(kStr).(kmerDf.kmer), :].count)
end

top50kmers = sort(collect(kmerTotalCounts), by = x -> x[2])[length(kmerTotalCounts)-49:end]
```

Wow, we have a lot of data! Now that we have a lovely table, let's plot our data!

```julia
julia> using StatsPlots
top50kmersDf = kmerDf[in(first.(top50kmers)).(kmerDf.kmer), :]
top5kmers = first.(top50kmers)[1:5]
top5kmersDf = kmerDf[in(top5kmers).(kmerDf.kmer), :]

groupedbar(kmerDf.kmer, kmerDf.count, bar_position = :stack, group = kmerDf.subregion, xlabel = "Most Common 3-mers", ylabel = "Frequency", title = "50 Most Common 3-mers by Frequency", xrotation = 90, color_palette = palette(:Paired_8), linecolor = nothing, legend = :topleft, foreground_color_legend = nothing, size = (2000, 1300), left_margin = 20mm, bottom_margin = 10mm)

@df top5kmersDf groupedhist(:count, group = :subregion, bar_position = :dodge, color_palette = palette(:Paired_8), linecolor = nothing)
```

### Almost done!! On to the next analysis!