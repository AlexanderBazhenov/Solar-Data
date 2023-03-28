% Created: 2022-09-14
% get 2 files
function [out1, out2, out3, out4, out5] = getSolar2 (FN1, FN2)
% FN1 =  FN1, FN2 = FN2
% FName
FN1str=strrep(FN1,'_','-');
FN2str=strrep(FN2,'_','-');
%
blankpos = findstr(FN1,'_');
cpos = findstr(FN1,'c');
mpos = findstr(FN1,'m');
FNstr = FN1(blankpos+1:cpos-2);
FNstr = strrep(FNstr,'_','-');
%
Lambdastr = FN1(blankpos+1:mpos-2);
Threadstr = FN1(mpos+2:cpos-2);
% remove mm
mpos = findstr(Threadstr,'m');
if (mpos>1)
   Threadstr = Threadstr(1:mpos-1); 
end
% Data
x1 = csvread (FN1);
x1(1,:)=[];
x2 = csvread (FN2);
x2(1,:)=[];
% OUTs
out1 = x1;
out2 = x2;
out3 = FNstr;
out4 = Lambdastr;
out5 = Threadstr;
endfunction
