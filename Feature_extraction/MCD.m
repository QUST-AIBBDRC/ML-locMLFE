clear all
clc
%%%%%
fid=fopen('quweidian22.txt');
string=fscanf(fid,'%s'); %�ļ�����
%ƥ����ַ���
firstmatches=findstr(string,'>')+7;%��ʼλ��
endmatches=findstr(string,'>')-1;
firstnum=length(firstmatches); %firstnum=endnum���е�����
endnum=length(endmatches);
  for k=1:firstnum-1
    j=1;
    lensec(k)=endmatches(k+1)-firstmatches(k)+1;%ÿ�����еĳ���
   for mm=firstmatches(k):endmatches(k+1)
        sequence(k,j)=string(mm); %�ַ�����
        j=j+1;
   end
   
  end
  %��������ȡÿ�����У��������ѡ��sequence(1,:)
  Auto=[];
for i=1:firstnum-1

% M1(i,:)=Auto1(sequence(i,1:lensec(i)),OriginData,lag);
% M2(i,:)=Auto2(sequence(i,1:lensec(i)),OriginData,lag);
% M3(i,:)=Auto3(sequence(i,1:lensec(i)),OriginData,lag);

M(i,:)=MCDZD(sequence(i,1:lensec(i)));
Auto=[Auto;M(i,:)];
end
save quweidian22.data_Auto.mat Auto


