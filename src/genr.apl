m←genr;x;y;c;d;el;ec
c←100
x←4×grand c⍴1
y←4×grand c⍴1
xb←x+⍳c
yb←y+.5×⍳c
xa←(+/xb)÷⍴xb
ya←(+/yb)÷⍴yb
av←xa,ya
d←2 c⍴xb,yb
(⍉2 c⍴d) print 'genr.data'
co←covm d
el←eval co
ec←20×evec co
(2 2⍴(av-(ec[⎕io;])÷2),av+(ec[⎕io;])÷2) print 'genre1.data'
(2 2⍴(av-(ec[⎕io+1;])÷2),av+(ec[⎕io+1;])÷2) print 'genre2.data'
