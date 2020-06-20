awk 'BEGIN{while(getline < "te-genome-contrib.only-tes.tab"){gcontrib[$1]=$2}}{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"gcontrib[$4]"\t"$7}' danRer10.fa.out.bed |grep -v ")n"> danRer10.fa.out.bed.midibed
awk '{print $0 > $1".midibed"}' danRer10.fa.out.bed.midibed
rm chrUn_KN1*
