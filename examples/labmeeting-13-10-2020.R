source("genomidi.R")


# Load DBs
TEs.bed=read.table("data/danRer10.fa.align.recalculated.bed.order")
TEs.gr=GRanges(seqnames=TEs.bed$V1, ranges=IRanges(start=TEs.bed$V2,end=TEs.bed$V3),strand=TEs.bed$V6,names=TEs.bed$V4,kimura=TEs.bed$V7,fixed_vol=rep(100,times=length(TEs.bed$V7)))

GENEs.bed=read.table("data/refGene.bed.without-gene-name")
GENEs.gr=GRanges(seqnames=GENEs.bed$V1, ranges=IRanges(start=GENEs.bed$V2,end=GENEs.bed$V3),strand=GENEs.bed$V6,names=GENEs.bed$V4,note=rep(50,times=length(GENEs.bed$V7)),fixed_vol=rep(127,times=length(GENEs.bed$V7)))   #vol and notes range: 0-127


ATACSEQ.bed=read.table("data/WT_atac.peaks.col1-6.liftover-danRer10.asterix-strand.bed")
ATACSEQ.gr=GRanges(seqnames=ATACSEQ.bed$V1, ranges=IRanges(start=ATACSEQ.bed$V2,end=ATACSEQ.bed$V3),strand=ATACSEQ.bed$V6,names=rep("ATAC",times=length(ATACSEQ.bed$V1)),note=rep(57,times=length(ATACSEQ.bed$V1)),fixed_vol=rep(127,times=length(ATACSEQ.bed$V1)))   #vol and notes range: 0-127



# Define all tracks we want to listen

T.ATAC=Track(
 DB=ATACSEQ.gr,
 Aes=list(Name="names", Note="note", Vol="fixed_vol"),
 Name="ATAC", 
 Channel=6,
 Instrument=53 # Voice Ooh
# Instrument=19 # Organo de iglesia
)

T.GENE=Track(
 DB=GENEs.gr,
 Aes=list(Name="names", Note="note", Vol="fixed_vol"),
 Name="GENE", 
 Channel=5,
 Instrument=53 # Voice Ooh
)

T.LINE=Track(
 DB=TEs.gr,
 Aes=list(Name="names", Note="kimura", Vol="fixed_vol"),
 Name="LINE", 
 Channel=1,
 Instrument=43 # contrabajo
)

T.SINE=Track(
 DB=TEs.gr,
 Aes=list(Name="names", Note="kimura", Vol="fixed_vol"),
 Name="SINE", 
 Channel=2,
 Instrument=40 # violin
)

T.DNA=Track(
 DB=TEs.gr,
 Aes=list(Name="names", Note="kimura", Vol="fixed_vol"),
 Name="DNA", 
 Channel=3,
 Instrument=42 #cello
)

T.LTR=Track(
 DB=TEs.gr,
 Aes=list(Name="names", Note="kimura", Vol="fixed_vol"),
 Name="LTR", 
 Channel=4,
 Instrument=46 # orchestal harp
)

S=Symphony(
 Label="hoxc13a-at-chr-3",
 Chr="chr3",
 Start=23200000,
 End=25000000,
 Tempo=100000,
 TrackList=list(T.ATAC,T.GENE,T.LINE,T.SINE,T.DNA,T.LTR)
)

# And compose them by bringing all elements comprending the selected region
composeSymphony(S)

# Generate output files and reproduce
playSymphony(S)

# Save workspace
save.image()


