# Genomidi

R package for listening genomic annotations. 

# Dependences

+ R
+ GenomicRanges
+ Midicomp
+ Timidity

# Quick start

Basic steps to reveal genomic symphony.

1. Build a database with `GenomicRanges` package
2. Define `Track` and `Symphony` objects
3. Compile and listen your `mid` file



# Docker container

1. Clone this repo to setup your working directory.

```
git clone https://github.com/fgajardoe/genomidi.git
cd genomidi
```

2. Setup your `data` folder

```
cd data
wget ... # all annotations you wish
```

3. Get the docker image

```
docker pull fgajardoe/genomidi:latest
```

Alternatively you can build it your own with the `Dockerfile` included in the repo

4. Write your genomic symphony.

An example file is included in `symphony.R`.
Make sure your files in `data/` are the same declared in your script.

5. Compile and get your `mid` file

```
docker run -it -v `pwd`:/home/genomidi/ fgajardoe/genomidi:lastest Rscript symphony.R
```

Where `pwd` is your working directory (that is, the root of this repository).


# Genomidi objects

## Track

A `Track` object represents a set of genomic annotations of the same type, for example, the information you get from a `bed` file for genes, or transposable elements, or CHIPseq peaks, and so on.

```
Track.genes=Track(
 DB=db.gr,
 Aes=list(Name="names", Note="note", Vol="fixed_vol"),
 Name="Gene", 
 Channel=1,
 Instrument=1 # Piano
)
```

There is a few things to have in mind while setting up a `Track` object:

_In construction_

## Symphony

```
S=Symphony(
 Label="The-repetitive-order",
 Chr="chr19",
 Start=19490627,
 End=20123884,
 Tempo=100000,
 TrackList=list(T.ATAC,T.GENE,T.LINE,T.SINE,T.DNA,T.LTR)
)
```

# Useful links

+ [Note names, MIDI numbers and frequencies](https://newt.phys.unsw.edu.au/jw/notes.html)
+ [Midi instruments](http://fmslogo.sourceforge.net/manual/midi-instrument.html)


# Licence

GPLv3
