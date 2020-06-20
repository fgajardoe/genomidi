chr=$1
start=$2
end=$3
feature=$4


tail -n +2 $1.midibed.raw-values |awk -F'\t' '{start='${start}'; end='${end}'; if($3>start && $4<end){print $12"\n"$13 }}' |sort -gk1,1 > $1.region-${start}-${end}.pivot
tail -n +2 $1.midibed.raw-values|perl -ne'my $f="'$feature'"; if(m/$f/g){ print $_;} ' |awk -F'\t' '{start='${start}'; end='${end}'; if($3>start && $4<end){print $12" O\n"$13" O" }}' |grep -Ff - $1.midibed.raw-values.4midi > $1.4midi.subset
tail -n +2 $1.midibed.raw-values|perl -ne'my $f="'$feature'"; if(m/$f/g){ print $_;} ' |awk -F'\t' '{start='${start}'; end='${end}'; if($3>start && $4<end){print $0}}'  > $1.region-${start}-${end}.${feature}.raw-values

Rscript correct-offset.R $1.4midi.subset $1.region-${start}-${end}.pivot
#Rscript correct-offset.R $1.4midi.subset
tail -n +2 $1.4midi.subset.time0 | cut -f2- -d' '|awk '{print $6" "$2" "$3" "$4" "$5}'|sort -u|sort -gk1,1 > $1.region-${start}-${end}.${feature}.4midi
rm $1.4midi.subset $1.4midi.subset.time0
sh generate-midi.sh $1.region-${start}-${end}.${feature}.4midi
timidity $1.region-${start}-${end}.${feature}.4midi.mid
