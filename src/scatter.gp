
set terminal jpeg enhanced size 512,480 font NimbusRoman 10

set output "s5000.jpg"

set title "horizontal sdev 5.0, vertical sdev 2.0, 5000 samples"

plot [-20:20] [-20:20]  "s5000.data" with points pt 7 ps .15 not

set output "s500.jpg"

set title "horizontal sdev 5.0, vertical sdev 2.0, 500 samples"

plot [-20:20] [-20:20]  "s500.data" u 1:2 with points pt 7 ps 1 not
