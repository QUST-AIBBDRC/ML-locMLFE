clear all
clc
L=5;
%%%%%
fid=fopen('testlablezx.txt');
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
for i=1:firstnum-1
eb1(i,:)= ebgw1(sequence(i,1:lensec(i)),L);
eb2(i,:)= ebgw2(sequence(i,1:lensec(i)),L);
eb3(i,:)= ebgw3(sequence(i,1:lensec(i)),L);
end
save testlablezx.EBGW_5.mat eb1 eb2 eb3

