m←c pca e;x;y;c;d;el;ec;xb;yb;xa;ya;co;av;s
⍝c←100
x←4×grand c⍴1
y←4×grand c⍴1
xb←x+⍳c
yb←y+e×⍳c
xa←(+/xb)÷⍴xb
ya←(+/yb)÷⍴yb
av←xa,ya
d←2 c⍴xb,yb
(⍉2 c⍴d) print 'pca.data'
co←covm d
el←eval co
ec←20×evec co

'eigenvectors' print 'pcaeigensystem.txt'
(ec÷20) print '>pcaeigensystem.txt'
'eigenvalues' print '>pcaeigensystem.txt'
el print '>pcaeigensystem.txt'

xa← 100×(1↑el)÷+/el
xb← 100×(¯1↑el)÷+/el

s←⍕'set label "λ0 = ', (1↑el), ' (', xa, '%)" at graph .2,.8'
s print 'pcalabels.gp'
s←⍕'set label "λ1 = ', (¯1↑el), ' (', xb, '%)" at graph .2,.75'
s print '>pcalabels.gp'

(2 2⍴(av-(ec[⎕io;])÷2),av+(ec[⎕io;])÷2) print 'pcae1.data'
(2 2⍴(av-(ec[⎕io+1;])÷2),av+(ec[⎕io+1;])÷2) print 'pcae2.data'



