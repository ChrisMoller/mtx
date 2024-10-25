
set terminal jpeg enhanced size 512,512 font NimbusRoman 10

set output "genr.jpg"

set title "cov"

plot  [0:120][0:120] "genr.data" using 1:2 with points pt 7 ps .75 t 'data', \
   "genre1.data" using 1:2 with lines t "1", \
   "genre2.data" using 1:2 with lines t "2"

#plot  "genr.data" using 1:2 with points pt 7 ps .75 t 'data', \
#      "genraxes.data" using 1 with lines t "1", \
#      "genraxes.data" using 2 with lines t "2"

#plot  "genraxes.data" using 1:2 with lines t "1", \

#   xa - ec[0,0]/2   ya - ec[0,1]/2
#   xa + ec[0,0]/2   ya + ec[0,1]/2

#   xa - ec[1,0]/2   ya - ec[1,1]/2
#   xa + ec[1,0]/2   ya + ec[1,1]/2
