m←c pca a;x;y;c;d;el;e;ec;xb;yb;xa;ya;co;av;s;⎕io;m0;m1;h
⎕io←0
e←tand a
x←4×grand c⍴1
y←4×grand c⍴1
xb←x+⍳c
yb←y+e×⍳c
xa←(+/xb)÷⍴xb
ya←(+/yb)÷⍴yb
av←xa,ya
d←2 c⍴xb,yb
⊣(⍉2 c⍴d) print 'pca.data'
co←covm d
el←eval co
ec←20×evec co

⊣'eigenvectors' print 'pcaeigensystem.txt'
⊣(ec÷20) print '>pcaeigensystem.txt'
⊣'eigenvalues' print '>pcaeigensystem.txt'
⊣el print '>pcaeigensystem.txt'

xa← 100×(1↑el)÷+/el
xb← 100×(¯1↑el)÷+/el

m0←atand ÷/⌽1 2↑ec
m1←atand ÷/⌽¯1 2↑ec

h←'w' ⎕fio[3] 'pcalabels.gp'   ⍝ open

s←⍕"set label \"eigenvalue λ0 = ", (1↑el), " (", xa, "%)\" at graph .1,.8\n"
⊣'%s' s ⎕fio[22] h  ⍝ fwrite

s←⍕"set label \"eigenvalue λ1 = ", (¯1↑el), " (", xb, "%)\" at graph .1,.75\n"
⊣'%s' s ⎕fio[22] h  ⍝ fwrite

s←⍕"set label \"slope λ0 = ", m0, " degress (nominally ", a, ")\" at graph .1,.7\n"
⊣'%s' s ⎕fio[22] h  ⍝ fwrite

s←⍕"set label \"slope λ1 = ", m1, " degrees (nominally ", (a-90), ")\" at graph .1,.65\n"
⊣'%s' s ⎕fio[22] h  ⍝ fwrite

⊣⎕fio[4] h


⊣(2 2⍴(av-(,1 2↑ec)÷2),av+(,1 2↑ec)÷2) print 'pcae1.data'
⊣(2 2⍴(av-(,¯1 2↑ec)÷2),av+(,¯1 2↑ec)÷2) print 'pcae2.data'

