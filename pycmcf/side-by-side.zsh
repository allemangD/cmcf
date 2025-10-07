#ffmpeg -frameraet 5 -i out/stage%03d
mkdir -p side-by-side
for n in $(seq -f %03.0f 1 25); do
  montage out/stage$n-edges{,-rev}.png -geometry +0 -gravity center -crop 85%x100% -tile 2x1 side-by-side/stage$n-edges.png;
done

ffmpeg -framerate 5 -i side-by-side/stage%03d-edges.png -crf 19 stages-sbs-edges.webm -y

#f
#for n in {0..025}; do
#  echo $n;
#done

#for n in 000 005 010 015 020 025; do
#    montage export/stage$n{,-rev}.png -geometry +0 -gravity center -crop 85%x100% -tile 2x1 side-by-side-$n.png;
#done
