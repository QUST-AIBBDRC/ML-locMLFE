clear all
clc
%%%%%
fid=fopen('quweidian22.txt');
string=fscanf(fid,'%s'); %文件输入
%匹配的字符串
firstmatches=findstr(string,'>')+7;%开始位置
endmatches=findstr(string,'>')-1;
firstnum=length(firstmatches); %firstnum=endnum序列的条数
endnum=length(endmatches);
  for k=1:firstnum-1
    j=1;
    lensec(k)=endmatches(k+1)-firstmatches(k)+1;%每条序列的长度
   for mm=firstmatches(k):endmatches(k+1)
        sequence(k,j)=string(mm); %字符序列
        j=j+1;
   end
   
  end
  %上面是提取每条序列，下面调用选用sequence(1,:)
  Auto=[];
for i=1:firstnum-1

% M1(i,:)=Auto1(sequence(i,1:lensec(i)),OriginData,lag);
% M2(i,:)=Auto2(sequence(i,1:lensec(i)),OriginData,lag);
% M3(i,:)=Auto3(sequence(i,1:lensec(i)),OriginData,lag);

M(i,:)=MCDZD(sequence(i,1:lensec(i)));
Auto=[Auto;M(i,:)];
end
save quweidian22.data_Auto.mat Auto


