cat header.txt $1 | awk '{print $0}END{print "TrkEnd"}' |midicomp -d -c $1.mid
