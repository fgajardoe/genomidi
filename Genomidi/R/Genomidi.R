# Dependencies
library(GenomicRanges)

# Prior definitions

Track=setRefClass("Track",fields=list(DB="GRanges",Aes="list",Name="character",Channel="numeric",Instrument="numeric",Elements="GRanges"))

Symphony=setRefClass("Symphony",
		     fields=list(Label="character",Chr="character",Start="numeric",End="numeric",Tempo="numeric",TrackList="list"))


# Functions

composeSymphony=function(symphony){

	# Get GRanges object for the region
	region.gr=GRanges(seqnames=symphony$Chr, ranges=IRanges(start=symphony$Start,end=symphony$End))

	TrackNames=c()
	# load al instruments

	# Pull elements corresponding to each track in the region
	for(track in symphony$TrackList){

		# Find all elements overlapping region
		region.elements=subsetByOverlaps(track$DB,region.gr)

		start(region.elements)=start(region.elements)-symphony$Start+1
		end(region.elements)=end(region.elements)-symphony$Start

		# Exclude elements not overlapping completely
		region.elements=region.elements[(start(region.elements) > 0) & (end(region.elements) <= symphony$End),]

		# Generate new GenomicRanges objects for tracks
		e=paste("region.elements[region.elements$",track$Aes$Name," == track$Name,]",sep="")
		track$Elements=eval(parse(text=e), environment())
		
		#track$Elements$n=FUN.note(paste("track$Elements$",track$Aes$Note,sep=""))
		#track$Elements$v=FUN.vol(paste("track$Elements$",track$Aes$Vol,sep=""))

	}
}


playSymphony=function(symphony,play=FALSE){
	write.table(printSymphony(symphony),paste(symphony$Label,".asc",sep=""),quote=F,sep=" ",row.names=F,col.names=F)
	if(play==FALSE){
		system(paste("midicomp -c ",symphony$Label,".asc ",symphony$Label,".mid", sep="")) # only generate mid file
	}
	else{
		system(paste("midicomp -c ",symphony$Label,".asc ",symphony$Label,".mid && timidity ",symphony$Label,".mid",sep="")) # generate mid file and play it
	}
}


printSymphony=function(symphony){
	t=t(data.frame(c("MFile","1",as.character(length(symphony$TrackList)+1),"960",""),c("MTrk","","","",""),c("0","Tempo",as.character(format(symphony$Tempo,scientific=F)),"",""),c("TrkEnd","","","","")))
	colnames(t)=c("time","action","param1","param2","param3")
	for(track in symphony$TrackList){
		t=rbind(t,printTrack(track))
	}
	rownames(t)=NULL
	return(t)
}

printTrack=function(track){

	# printing track start and instrument (header)
	ch=paste("ch=",as.character(track$Channel),sep="")
	p=paste("p=",as.character(track$Instrument),sep="")

	track.name=paste("\"",track$Name,"\"",sep="")

	track.header=t(data.frame(c("MTrk","","","",""),c("0","PrCh",ch,p,""),c("0","Meta","TrkName",track.name,"")))
	colnames(track.header)=c("time","action","param1","param2","param3") 

	# printing notes
	n=eval(parse(text=paste("track$Elements$",track$Aes$Note,sep="")),environment())
	v=eval(parse(text=paste("track$Elements$",track$Aes$Vol,sep="")),environment())
	n=paste("n=",as.integer(n),sep="")
	v=paste("v=",as.integer(v),sep="")

	# pressed (On)
	s=start(track$Elements)
	on.df=data.frame(s,"On",ch,n,v)
	colnames(on.df)=c("time","action","param1","param2","param3")

	# lifted (Off)
	e=end(track$Elements)
	off.df=data.frame(e,"Off",ch,n,"v=0")
	colnames(off.df)=c("time","action","param1","param2","param3")
	
	# concatenate header and notes
	df=rbind(on.df,off.df)
	df=df[order(df$time),]
	df$time=as.character(df$time)
	df=rbind(track.header,df)

	# concatenate end of the track
	track.end=t(data.frame(c("TrkEnd","","","","")))
	colnames(track.end)=c("time","action","param1","param2","param3") 
	df=rbind(df,track.end)

	# reset rownames and return
	rownames(df)=NULL
	return(df)
}




