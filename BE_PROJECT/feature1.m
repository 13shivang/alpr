function [ p ] = feature1( A )
% A=im2bw(A);
A=imresize(A,[48 24]);
add=0;
flag=0;
k=1;

for i=1:24
    flag=flag+1;
    for j=1:48
        add=add+A(j,i);%sum of i th column
    end;
    if(mod(i,2)==0)
        ver(k) = add;
        k = k+1;
        add=0;
    end;
    if(flag==2)
        i=i+2;
        flag=0;
    end;
end;
        
add=0;
flag=0;
k=1;

for i=1:48
    flag=flag+1;
    for j=1:24
        add=add+A(i,j);%sum of i th row
    end;
    if(mod(i,2)==0)
        hor(k) = add;
        k = k+1;
        add=0;
    end;
    if(flag==2)
        i=i+2;
        flag=0;
    end;
end;
        
zone = mat2cell(A, [12 12 12 12], [6 6 6 6]);% dividing into zones of size 12*6
sum = 0;
for i = 1:4
    for j = 1:4     %for each zone
        temp = zone{i,j};
        for k = 1:12
            for l = 1:6
                sum = sum + temp(k, l);
            end;
        end;
        zone_add(i,j) = sum; %assign sum of each pixel of a particular zone
        sum = 0;
    end;
end;

p = [ver hor reshape(zone_add,[1,16])];

end

