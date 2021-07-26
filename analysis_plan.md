# Ivy's Analysis Plan!!

This is an analysis repository to analyze **Sars-CoV-related genomes** using functions from the *BioinformaticsBISC195* package.

We have calculated and plotted sequence lengths and GC content, and calculated kmer composition. Here, we will discuss *two* additional sequence features to analyze:

## 1. We will analyze the kmer compositions of these genomes across geographical regions.

### What we have already done:
* We have already built a function `getKmers` to determine the unique kmer compositions of a sequence according to a given *k* length.

### Next steps:
0. First, we will use the `parse_fasta` function to parse the fasta file containing our coronavirus genomes and retrieve a vector of header attributes and a vector of the corresponding sequences.

1. We must record the number of occurrences of each kmer in each sequence(I have settled on only collecting the 3-mers, 4-mers, 5-mers, and 6-mers since there is way too much data otherwise).

    We will do this by creating a modified version of `getKmers` called `getKmerCount` which will take a sequence and a *k* length and return a dictionary of all the unique kmers of *k* length in that sequence, mapped as key-value pairs to their frequency counts.

2. We will organize each set of unique kmers by geographic region and order them by frequency in DataFrame. Perhaps something that looks like this:

    |Unique Kmer|*k* Length|Count|Accession|Location|
    |------|------|------|------|------|
    |kmer 1|length 1|count 1|id 1|location 1|
    |kmer 2|length 2|count 2|id 2|location 2|
    |kmer 3|length 3|count 3|id 3|location 3|

    To do this, we will define a new function `getHeaderAttrib` that takes a single header string, the desired delimiter (in our case, it will always be "|"), and an array of indices and returns a vector of the header attribute values at those indices. This will allow us to only loop through every header once and take the data we need, as in the accession number and location.

    We will also need to transform the location into a region somehow. To do this, we will grab only the country name from each location string. Then, we will build upon a JSON file created by the Github user *lukes* which maps countries to their identifying attributes such as region, sub-region, and country code. Using that mapping, we will match each kmer in the DataFrame to the associated region and subregion.

3. Using StatsPlots, we will then plot multiple bar charts for each *k* length by kmer and frequency.
    
    
## 2. We will analyze the gc content of these genomes over time.

### What we have already done:
* We have built a function that calculates the gc content of a sequence. 

### Next steps:
0. We first will use the `parse_fasta` function to parse the fasta file containing our coronavirus genomes and retrieve a vector of header attributes and a vector of the corresponding sequences. (While doing the first analysis, I ended up saving all the data to files to make the process quicker for myself, so I just loaded the files at the beginning of this analysis)

1. We will prepare our data for analysis. Some of the sequences do not have the collection date information available. In that case, we will want to exclude those points from our data as to avoid any assumptions.

    To do this, we will define a general function `filterByEmptyAttrib` which takes the header vector, sequence vector, specific delimiter, and the index of the header attribute to filter by. This function will check which headers in the header vector do not contain attribute information at the given index (for our data, it will be index 5 for collection date) and return an updated version of the header vector with those entries removed as well as an uodated version of the sequence vector with the corresponding entries removed.

2. We must take the gc content of each sequence in our data by mapping our sequence vector with the `gc_content` function to generate a new vector that contains the corresponding gc content of each sequence.

3. Next, we will want to get associated accession number and collection date to each sequence. We will do this using the `getHeaderAttrib` function that we defined in the last analysis.

4. To be able to order the sequences and their gc contents chronologically, it would be best for the time data to be held in date objects, not strings. We will map the collection date vector to transform each string into a date object with format *yyyy-mm-dd*.

5. At this point, all the data we need will be prepared so we can organize it in a DataFrame like this:

    |Accession|GC Content|Collection Date|
    |---|---|---|
    |id 1|gc_content 1|date 1|
    |id 2|gc_content 2|date 2|
    |id 3|gc_content 3|date 3|

2. Using StatsPlots, we will plot each sequence by their gc content and collection date on a scatter plot. Depending on how the data looks, we might want to refine it or zoom in as to focus on some important points in time (such as the recent pandemic).

### Packages to be used and their purpose:
|Package|Purpose|
|--------|-------|
|BioinformaticsBISC195|interact with genomes/data held in fasta files|
|DataFrames|tables|
|StatsPlots|plots|
|Statistics|general statistics functions|
|JLD2, FileIO|save and load large data|
|JSON|parse and write JSON files|
|Dates|turn strings into Date objects|

To finish off this plan, here are my cats!
![Cat#1](notebooks/assets/cat1.jpg)
![Cat#2](notebooks/assets/cat2.jpg)
