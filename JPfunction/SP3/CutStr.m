function Out=CutStr(Str,Sym,Order)

%Made by Sangho Kim 2006-08-02
%Cutting String with Symbol by folloed order number

len=length(Str);

Sym_Order=findstr(Str,Sym);

Size_Order=size(Sym_Order,2);

New_Mat=ones(1,Size_Order+2);

New_Mat(1)=0;

New_Mat(1,2:Size_Order+1)=Sym_Order;

New_Mat(Size_Order+2)=len+1;

Out=Str(New_Mat(Order)+1:New_Mat(Order+1)-1);