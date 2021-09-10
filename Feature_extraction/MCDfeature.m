function [set] =MCDfeature( A)
wei=21;
A=cell2mat(A);
set=[];
for j=1:21
B=A(:,j);
chen=[];
for i=1:label
    C=B(i);
str=strrep(C,'A','01');
str=strrep(str,'C','02');
str=strrep(str,'D','03');
str=strrep(str,'E','04');
str=strrep(str,'F','05');
str=strrep(str,'G','06');
str=strrep(str,'H','07');
str=strrep(str,'I','08');
str=strrep(str,'K','09');
str=strrep(str,'L','10');
str=strrep(str,'M','11');
str=strrep(str,'N','12');
str=strrep(str,'P','13');
str=strrep(str,'Q','14');
str=strrep(str,'R','15');
str=strrep(str,'S','16');
str=strrep(str,'T','17');
str=strrep(str,'V','18');
str=strrep(str,'W','19');
str=strrep(str,'Y','20');
str=strrep(str,'X','21');
chen=[chen;str];
str=[];
end
shuju=str2num(chen);
set=[set,shuju];
end
% s=zeros(label,(wei-1));
% for j=1:(wei-1)
%     for i=1:label
%    s(i,j)= wei*(set(i,j)-1)+set(i,j+1);
%     end
% end
end

    

