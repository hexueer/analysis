# Ivy's Analysis Plan!!

This is an analysis repository to analyze **Sars-CoV-related genomes** using functions from the *BioinformaticsBISC195* package.

We have calculated and plotted sequence lengths and GC content, and calculated kmer composition. Here, we will discuss *two* additional sequence features to analyze:

1. We will analyze the kmer compositions of these genomes across different geographical regions. We have already built a function to determine the unique kmer compositions of a sequence according to a given *k* length. We must make a small edit to record the number of occurrences of each one. Next, we will organize each set of unique kmers by geographic region and order them by frequency in table. We will then take the most common ones (number depending on the outcome of the data) to plot on a bar graph by *k* length and frequency. I also wanted to ask: what is a de Bruijn graph? I read articles but I still don't understand.

2. We will analyze the gc content of these genomes over time. We have built a function that calculates the gc content of a sequence. Next, we must take the gc content of each sequence in our data, organize them in an ordered table with their associated sequence and collection date. Then, we will plot each sequence by their gc content and collection date on a scatter plot. This will bring us to a good place to estimate a linear regression and conduct a significance test.

The tables and plots will be made using DataFrames.jl and Plots.jl (I will also attempt more intricate visualisations using Gadfly).

Any statistical analysis will be done using Statistics.jl and SciKitLearn.jl.