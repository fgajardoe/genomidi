#ffmpeg -i $1 -acodec libmp3lame -ab 64k $1.mp3 # -> algo raro, pero bueno resulta de ese

timidity $1 -Ow -o - |ffmpeg -i - -acodec libmp3lame -ab 64k $1.mp3
