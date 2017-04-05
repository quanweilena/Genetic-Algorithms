function sample=greedy(sample,c,g,m1)
[m n]=size(sample);

for i=1:m
    label=g;
    label=label.*sample(i,:);
    if sum(sample(i,:).*c)>m1
        %ÐÞ¸´
        
       
        while sum(sample(i,:).*c)>m1
           min=10000;
           for j=1:n
                if label(j)~=0&&label(j)<min
                    min=label(j);
                    modify=j;
                end
           end
            sample(i,modify)=0;
            label(modify)=0;
        end
    end
end
weight=sample*c';

        