
c=[80 82 85 70 72 70 66 50 55 25 50 55 40 48 50 32 22 60 30 32 40 38 35 32 25 28 30 22 50 30 45 30 60 50 20 65 20 25 30 10 20 25 15 10 10 10 4 4 2 1 84 81 85 77 72 70 66 50 55 25 55 55 40 48 43 32 25 60 30 32 40 37 35 12 25 28 30 52 50 30 35 30 60 50 20 75 20 25 30 10 20 25 15 10 10 10 4 4 2 3 ];
value=[220 208 198 192 180 180 165 162 160 158 155 130 125 122 120 118 115 110 105 101 100 100 98 96 95 90 88 82 80 77 75 73 72 70 69 66 65 63 60 58 56 50 30 20 15 10 8 5 3 1 220 208 198 192 180 180 165 162 160 158 155 130 125 122 120 118 115 110 105 101 100 100 98 96 95 90 88 82 80 77 75 73 72 70 69 66 65 63 60 58 56 50 30 20 15 10 8 5 3 1];
g=value./c;%价值重量比
[h,N]=size(c);L=N;
m1=1000;%背包承重量
itt=1;xunhuan=30; pc = 0.5; pm = 0.01;%遗传代数位200，交叉率为0.5，变异率为0.01
p1=0.3;p3=10;
sol1=1;  it = 1;updatef=1;d=[];sss=floor(1/p1);ave=[];
GER=800;
ger =200;%代数

disp(sprintf('Number of generations: %d',GER));
disp(sprintf('Population size: %d',N));
disp(sprintf('Crossover probability: %.3f',pc));
disp(sprintf('Mutation probability: %.3f',pm));
while itt<=xunhuan
v = round(rand(N,L));%产生01随机矩阵
sol1=1;  it = 1;updatef=1;

 u=0;
v=greedy(v,c,g,m1);%贪婪算法修复解
fit =v*value';%适应度函数采用01序列产生的总价值
t0 = clock;%t0返回当前的时间

while it <= GER
    ss=mod(it,ger);%数据更新
    if ss==0
        for i=1:p1*N
 uu=ceil(rand()*sss);
 u=u+uu;
 value(u)=value(u)+(2*rand()-1)*p3; 
        end
        g=value./c;  
     updatef=1;
     u=0;
      
    end
% 选择复制,采用赌轮盘选择法选择实现
totalfit=sum(fit);
fitvalue=fit/totalfit;
fitvalue=cumsum(fitvalue);
for i=1:N
   p=rand(1); sindex=1;
   while p > fitvalue(sindex) %蒙特卡洛方法抽样
      sindex=sindex+1; 
   end
   newv(i,:)=v(sindex,:); 
end
for i=1:N
   v(i,:)=newv(i,:);
end
% Crossver交叉
for i=1:N
   cindex(i)=i;
end
for i=1:N %产生要配对的父代的序号；经过N次顺序调换，将原有顺序打乱，使相邻两个个体作为交叉的父代
   point=unidrnd(N-i+1);
   temp=cindex(i);
   cindex(i)=cindex(i+point-1);
   cindex(i+point-1)=temp;
end

for i=1:2:N
   p=rand(1);
   if(p<pc)
      point=unidrnd(L-1)+1;%1<point<L 产生交叉点
      for j=point:(L-1) 
         ch=v(cindex(i),j);
         v(cindex(i),j)=v(cindex(i+1),j); %cindex中相邻的两个为两个父代的序号
         v(cindex(i+1),j)=ch;
      end
   end
end
% Mutation变异
M=rand(N,L)<=pm;%产生（N，L）维的01矩阵，为1的位置进行变异
v=v-2.*(v.*M)+M;%
v=greedy(v,c,g,m1);%贪婪算法修复解
 fit = v*value';
 ave(itt,it)=(sum(fit))/N;
[sol1,indb1] = max(fit); %

%if updatef>=sol1
%   sol1=updatef;
%end
%updatef=sol1;
d(itt,it)=sol1;
	it = it + 1;

end;
itt=itt+1;
end
ave=sum(ave);
d=sum(d);
d=d/xunhuan;
ave=ave/xunhuan;
figure(1);
hold on;
plot(d,'k','linewidth',1.5),xlabel('进化代数'),ylabel('适应度');
hold on;
plot(ave,'k--','linewidth',1.5);
legend('最佳适应度','平均适应度');
T = etime(clock,t0); 
disp(fprintf('the total time is: %2.4f',T));



    
