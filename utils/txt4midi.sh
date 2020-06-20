tail -n+2 $1 |awk '{print $12" On ch=1 n="$10" v="$11"\n"$13" Off ch=1 n="$10" v=0"}'  |sort -gk1,1 > $1.4midi
