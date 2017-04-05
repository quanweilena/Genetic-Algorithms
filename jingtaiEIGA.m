clear all;
c=[80 82 85 70 72 70 66 50 55 25 50 55 40 48 50 32 22 60 30 32 40 38 35 32 25 28 30 22 50 30 45 30 60 50 20 65 20 25 30 10 20 25 15 10 10 10 4 4 2 1 84 81 85 77 72 70 66 50 55 25 55 55 40 48 43 32 25 60 30 32 40 37 35 12 25 28 30 52 50 30 35 30 60 50 20 75 20 25 30 10 20 25 15 10 10 10 4 4 2 3 ];
value=[220 208 198 192 180 180 165 162 160 158 155 130 125 122 120 118 115 110 105 101 100 100 98 96 95 90 88 82 80 77 75 73 72 70 69 66 65 63 60 58 56 50 30 20 15 10 8 5 3 1 220 208 198 192 180 180 165 162 160 158 155 130 125 122 120 118 115 110 105 101 100 100 98 96 95 90 88 82 80 77 75 73 72 70 69 66 65 63 60 58 56 50 30 20 15 10 8 5 3 1];
g=value./c;%��ֵ������
[h,N]=size(c);L=N;
m1=1000;%����������
itt=1;xunhuan=30; pc = 0.5; pm = 0.01;%�Ŵ�����λ200��������Ϊ0.5��������Ϊ0.01
sol1=1;  it = 1;updatef=1;d=[];ave=[];xd=[];theshold=[];dead=[];xave=[];E=[];
GER=300;
ger =100;%����
ri=0.1;%����������
R=ri*N;%�����������

disp(sprintf('��Ӣ�����Ŵ��㷨 Elitism-based Immigrant GA(EIGA)'));
disp(sprintf('���������� Number of generations: %d',GER));
disp(sprintf('��Ⱥ��ģ Population size: %d',N));
disp(sprintf('������� Crossover probability: %.3f',pc));
disp(sprintf('������� Mutation probability: %.3f',pm));
disp(sprintf('��Ӣ������� Ratio of elite immigrants: %.3f',ri));



while itt<=xunhuan
    c=[80 82 85 70 72 70 66 50 55 25 50 55 40 48 50 32 22 60 30 32 40 38 35 32 25 28 30 22 50 30 45 30 60 50 20 65 20 25 30 10 20 25 15 10 10 10 4 4 2 1 84 81 85 77 72 70 66 50 55 25 55 55 40 48 43 32 25 60 30 32 40 37 35 12 25 28 30 52 50 30 35 30 60 50 20 75 20 25 30 10 20 25 15 10 10 10 4 4 2 3 ];
value=[220 208 198 192 180 180 165 162 160 158 155 130 125 122 120 118 115 110 105 101 100 100 98 96 95 90 88 82 80 77 75 73 72 70 69 66 65 63 60 58 56 50 30 20 15 10 8 5 3 1 220 208 198 192 180 180 165 162 160 158 155 130 125 122 120 118 115 110 105 101 100 100 98 96 95 90 88 82 80 77 75 73 72 70 69 66 65 63 60 58 56 50 30 20 15 10 8 5 3 1];
v = round(rand(N,L));%����01�������
sol1=1;  it = 1;updatef=1;

 u=0;
v=greedy(v,c,g,m1);%̰���㷨�޸���
fit =v*value';%��Ӧ�Ⱥ�������01���в������ܼ�ֵ
t0 = clock;%t0���ص�ǰ��ʱ��

while it <= GER

%��Ӣ����  �滻R�����ĸ���
fit =v*value';%��Ӧ�Ⱥ�������01���в������ܼ�ֵ
sortfit=sortrows(fit);%����Ӧ������
threshold=sortfit(R);
dead=[];
E=[];
%ѡ��Ӣ
for i=1:(N-1)
    if fit(i)<= fit(i+1)
    E=v((i+1),:);    
    end
end
for i=1:N
   if threshold >= fit(i) 
      dead=[dead;[i]]; 
   end
end
numdead=length(dead);%��̭�ĸ�����Ŀ

for i=1:numdead
    ii=dead(i);
    %v(ii,:)=round(rand(1,L));
    % Mutation����
    M=rand(1,L)<=pm;%������N��L��ά��01����Ϊ1��λ�ý��б���
    ee=E-2.*(E.*M)+M;%
    v(ii,:)=ee;
end
v=greedy(v,c,g,m1);%̰���㷨�޸���
  
    
    
% ѡ����,���ö�����ѡ��ѡ��ʵ��
totalfit=sum(fit);
fitvalue=fit/totalfit;
fitvalue=cumsum(fitvalue);
for i=1:N
   p=rand(1); sindex=1;
   while p > fitvalue(sindex) %���ؿ��巽������
      sindex=sindex+1; 
   end
   newv(i,:)=v(sindex,:); 
end
for i=1:N
   v(i,:)=newv(i,:);
end
% Crossver����
for i=1:N
   cindex(i)=i;
end
for i=1:N %����Ҫ��Եĸ�������ţ�����N��˳���������ԭ��˳����ң�ʹ��������������Ϊ����ĸ���
   point=unidrnd(N-i+1);
   temp=cindex(i);
   cindex(i)=cindex(i+point-1);
   cindex(i+point-1)=temp;
end

for i=1:2:N
   p=rand(1);
   if(p<pc)
      point=unidrnd(L-1)+1;%1<point<L ���������
      for j=point:(L-1) 
         ch=v(cindex(i),j);
         v(cindex(i),j)=v(cindex(i+1),j); %cindex�����ڵ�����Ϊ�������������
         v(cindex(i+1),j)=ch;
      end
   end
end
% Mutation����
M=rand(N,L)<=pm;%������N��L��ά��01����Ϊ1��λ�ý��б���
v=v-2.*(v.*M)+M;%
v=greedy(v,c,g,m1);%̰���㷨�޸���



 fit = v*value';
 ave(itt,it)=(sum(fit))/N;
[sol1,indb1] = max(fit); %

if updatef>=sol1
   sol1=updatef;
end
updatef=sol1;
d(itt,it)=sol1;
	it = it + 1;

end;
itt=itt+1;
end
xave=sum(ave);
xd=sum(d);
xd=xd/xunhuan;
xave=xave/xunhuan;
figure(1);
subplot(211)
hold on;
plot(xd,'b','linewidth',1.5),xlabel('��������'),ylabel('�����Ӧ��');
legend('SGA','RIGA','EIGA');
subplot(212)
hold on;
plot(xave,'b--','linewidth',1.5),xlabel('��������'),ylabel('ƽ����Ӧ��');
legend('SGA','RIGA','EIGA');
disp(sprintf('ȫ�����Ž�: %d',updatef));
T = etime(clock,t0); 
disp(fprintf('the total time is: %2.4f',T));



    
