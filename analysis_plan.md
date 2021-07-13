# Ivy's Analysis Plan!!

This is an analysis repository to analyze **Sars-CoV-related genomes** using functions from the *BioinformaticsBISC195* package.

We have calculated and plotted sequence lengths and GC content, and calculated kmer composition. Here, we will discuss *two* additional sequence features to analyze:

## 1. We will analyze the kmer compositions of these genomes across geographical regions.

### What we have already done:
* We have already built a function to determine the unique kmer compositions of a sequence according to a given *k* length.

### Next steps:
1. We must make a small edit to record the number of occurrences of each one.
2. We will organize each set of unique kmers by geographic region and order them by frequency in table.
3. We will then take the most common ones (number depending on the outcome of the data) to plot on a bar graph by *k* length and frequency.
    
*Note:* I also wanted to ask: what is a de Bruijn graph? I read articles but I still don't understand.
    
## 2. We will analyze the gc content of these genomes over time.

### What we have already done:
* We have built a function that calculates the gc content of a sequence. 

### Next steps:
1. We must take the gc content of each sequence in our data, organize them in an ordered table with their associated sequence and collection date.
2. We will plot each sequence by their gc content and collection date on a scatter plot.
3. We will estimate a linear regression to model the data calculate whether the trend is statistically significant.


### Packages to be used and their purpose:
|Package|Purpose|
|--------|-------|
|DataFrames.jl|tables|
|Plots.jl|plots|
|Gadfly|attempt at special visualizations|
|Statistics.jl|general statistics functions|
|SciKitLearn.jl|modeling data and significance tests|
