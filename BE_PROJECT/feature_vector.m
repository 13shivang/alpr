function [a] = feature_vector( I )
z=(feature1(I));
I=imresize(I,[42 24]);
holes=hole(I);
I_ori=I;
I = padarray(I,[2 2],'both');
I=imfill(I,'holes');
se = strel('disk',3);
I = imerode(I,se);
% figure();imshow(I);
%%%%%%%%%%%%%%%%%%%%%%HORIZONTAL LINE DETECTION%%%%%%%%%%%%%%%%%%%%%%%%
[~,ang]=imgradient(I);
ang=round(ang);
%Blakens the pixels with angles other than -90 or 90 as horizontal lines
%are at 90 degree
for i=1:46
    for j=1:28
        if(ang(i,j)==90||ang(i,j)==-90)
            I_hori(i,j)=255;
        else
            I_hori(i,j)=0;
        end
    end
end

I_hori=bwareaopen(I_hori,15); %Removes small unwanted pixels
%figure();imshow(I_hori);
se = strel('disk',3);
I_hori = imclose(I_hori,se);
%figure();imshow(I_hori);
[L,hori] = bwlabel(I_hori); %counts the horizontal lines


%%%%%%%%%%%%%%%%%%%%%%%VERTICAL LINE DETECTION%%%%%%%%%%%%%%%%%%%%%%%%
I_verti=rot90(I,-1); %Rotates image by 90 in anti-clock wise sense and
                     %then detecting horizontal lines
[~,ang1]=imgradient(I_verti);
ang1=round(ang1);
%Blakens the pixels with angles other than -90 or 90 as horizontal lines
%are at 90 degree
for i=1:28
    for j=1:46
        if(ang1(i,j)==90||ang1(i,j)==-90)
            I_verti(i,j)=255;
        else
            I_verti(i,j)=0;
        end
    end
end
I_verti=rot90(I_verti); %Rotating the image by 90 in clockwise direction
se = strel('disk',1);
I_verti = imclose(I_verti,se);
%figure();imshow(I_verti);
I_verti=bwareaopen(I_verti,24); %Removes small unwanted pixels
%figure();imshow(I_verti);
se = strel('disk',5);
I_verti = imclose(I_verti,se);
%figure();imshow(I_verti);
[L,verti] = bwlabel(I_verti);
I=I_ori;

approx = imresize(I, [7 4]);
approx = reshape(approx, [1 28]);

[~,ang]=imgradient(I);
ang=round(ang);
%r=/
%l=\
temp_l1=0;
temp_l2=0;

temp_r1=0;
temp_r2=0;

l=0;
r=0;

for i=1:42
    for j=1:24
        if((ang(i,j)>20 && ang(i,j)<=45)|(ang(i,j)>-160 && ang(i,j)<=-135))
            temp_l1=temp_l1+1;
        elseif((ang(i,j)>45 && ang(i,j)<=70)|(ang(i,j)>-135 && ang(i,j)<=-110))
            temp_l2=temp_l2+1;
        end   
        if((ang(i,j)<-20 && ang(i,j)>=-45)|(ang(i,j)<160 && ang(i,j)>=135))
            temp_r1=temp_r1+1;
        elseif((ang(i,j)<-45 && ang(i,j)>=-70)|(ang(i,j)<135 && ang(i,j)>=110))
            temp_r2=temp_r2+1;
        end
    end
end

% temp_l1
% temp_l2
% 
% temp_r1
% temp_r2

if(temp_l1>45|temp_l2>45|temp_r1>45|temp_r2>45)
    if(abs(temp_l1-temp_r1)<20)
        l=1;
        r=1;
    elseif((temp_l1>temp_r1)&&temp_l1>40)
        l=1;
    elseif((temp_l1<temp_r1)&&temp_r1>40)
        r=1;
    end
end
a =transpose( [ holes verti hori r l z ])

end

