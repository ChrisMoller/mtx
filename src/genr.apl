m←c genr e;x;y;c;d;el;ec;xb;yb;xa;ya;co;av;s
⍝c←100
x←4×grand c⍴1
y←4×grand c⍴1
xb←x+⍳c
yb←y+e×⍳c
xa←(+/xb)÷⍴xb
ya←(+/yb)÷⍴yb
av←xa,ya
d←2 c⍴xb,yb
(⍉2 c⍴d) print 'genr.data'
co←covm d
el←eval co
ec←20×evec co

xa← 100×(1↑el)÷+/el
xb← 100×(¯1↑el)÷+/el

s←⍕'set label "λ0 = ', (1↑el), ' (', xa, '%)" at graph .2,.8'
s print 'labels.gp'
s←⍕'set label "λ1 = ', (¯1↑el), ' (', xb, '%)" at graph .2,.75'
s print '>labels.gp'

(2 2⍴(av-(ec[⎕io;])÷2),av+(ec[⎕io;])÷2) print 'genre1.data'
(2 2⍴(av-(ec[⎕io+1;])÷2),av+(ec[⎕io+1;])÷2) print 'genre2.data'

