% 2022-12-22
[tau1, w1, yint1] = DataLinearModel (input1, epsilon0);
yynew1 = yint1 -tau1(2)*xx1';
S = yynew1;
[mode1, modefreq1, freqs1, Ss1] = imodeR([inf(S), sup(S)]);
outer1 = infsup( min(inf(S)), max(sup(S)) )
inner1 = infsup( max(inf(S)), min(sup(S)) )
wid(outer1)
%
[tau2, w2, yint2] = DataLinearModel (input2, epsilon0);
yynew2 = yint2 -tau2(2)*xx1';
S = yynew2;
[mode2, modefreq2, freqs2, Ss2] = imodeR([inf(S), sup(S)]);
outer2 = infsup( min(inf(S)), max(sup(S)) )
inner2 = infsup( max(inf(S)), min(sup(S)) )
wid(outer2)
R21outer =  outer2 / outer1
wid(R21outer)
R21inner =  inner2 / inner1

inner2/inner1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S = yint1;
Outer1 = infsup( min(inf(S)), max(sup(S)) )
Inner1 = infsup( max(inf(S)), min(sup(S)) )
wid(Inner1)
wid(Outer1)
S = yint2;
Outer2 = infsup( min(inf(S)), max(sup(S)) )
Inner2 = infsup( max(inf(S)), min(sup(S)) )
wid(Inner2)
wid(Outer2)

R21Outer =  Outer2 / Outer1
wid(R21Outer)
