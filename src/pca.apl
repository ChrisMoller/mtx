m←c pca e;x;y;c;d;el;ec;xb;yb;xa;ya;co;av;s;⎕io;m0;m1;h
⎕io←0
x←4×grand c⍴1
y←4×grand c⍴1
xb←x+⍳c
yb←y+e×⍳c
xa←(+/xb)÷⍴xb
ya←(+/yb)÷⍴yb
av←xa,ya
d←2 c⍴xb,yb
x←(⍉2 c⍴d) print 'pca.data'
co←covm d
el←eval co
ec←20×evec co

x←'eigenvectors' print 'pcaeigensystem.txt'
x←(ec÷20) print '>pcaeigensystem.txt'
x←'eigenvalues' print '>pcaeigensystem.txt'
x←el print '>pcaeigensystem.txt'

xa← 100×(1↑el)÷+/el
xb← 100×(¯1↑el)÷+/el

h←'w' ⎕fio[3] 'pcalabels.gp'   ⍝ open

s←⍕'set label "λ0 = ', (1↑el), ' (', xa, '%)" at graph .2,.8'
x←'%s' s ⎕fio[22] h  ⍝ fwrite
x←⎕fio[16] h ⍝ flush
x←10⎕fio[42] h     ⍝ newline
⍝x←⎕fio[16] h ⍝ flush

s←⍕'set label "λ1 = ', (¯1↑el), ' (', xb, '%)" at graph .2,.75'
x←'%s' s ⎕fio[22] h  ⍝ fwrite
x←⎕fio[16] h ⍝ flush
x←10⎕fio[42] h     ⍝ newline
⍝x←⎕fio[16] h ⍝ flush

m0←atand ec[0;1]÷ec[0;0]
m1←atand ec[1;1]÷ec[1;0]

s←⍕'set label "angle 0 = ', m0, ' degress" at graph .2,.7'
x←'%s' s ⎕fio[22] h  ⍝ fwrite
x←⎕fio[16] h ⍝ flush
x←10⎕fio[42] h     ⍝ newline

s←⍕'set label "angle 1 = ', m1, ' degrees" at graph .2,.65'
x←'%s' s ⎕fio[22] h  ⍝ fwrite
x←⎕fio[16] h ⍝ flush
x←10⎕fio[42] h     ⍝ newline

x←⎕fio[4] h


x←(2 2⍴(av-(ec[0;])÷2),av+(ec[0;])÷2) print 'pcae1.data'
x←(2 2⍴(av-(ec[1;])÷2),av+(ec[1;])÷2) print 'pcae2.data'

