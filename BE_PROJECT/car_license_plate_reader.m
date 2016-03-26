clc;
close all;
clear all;
I=imread('carplate2.jpg');
I=imresize(I,[500 NaN]); 
I=rgb2gray(I);
I1=I;
 figure;
 imshow(I);

[rows cols] = size(I);

%--------- Dilate and Erode Image in order to remove noise-----------
tempI=I;
stre=strel('line',1,3);   
I=imdilate(I,stre);
%-------------------- PROCESS EDGES IN HORIZONTAL DIRECTION--------------
difference = 0;
total_sum = 0;
max_horz = 0;
maximum = 0;
for j = 2:cols
    sum = 0;
    for i = 2:rows
        difference = abs(uint32(I(i, j) - I(i-1,j)));
        if(difference > 50)
            sum = sum + difference;
        end
    end
    horz1(j) = sum;
   %----------------------------- Find Peak horizontal Value---------------
    if(sum > maximum)
        max_horz = j;
        maximum = sum;
    end
    total_sum = total_sum + sum;
end
average = total_sum / cols;


%------- Smoothen the Horizontal Histogram by applying Low Pass Filter
horz = horz1;
for i = 21:(cols-21)
sum = 0;
    for j = (i-20):(i+20)
        sum = sum + horz1(j);
    end
    horz(i) = sum / 40;
end

for j = 1:cols
    if(horz(j) < .35*average)
        horz(j) = 0;
        for i = 1:rows
            I(i, j) = 0;
        end
    end
end
%--------------- PROCESS EDGES IN VERTICAL DIRECTION
difference = 0;
total_sum = 0;
maximum = 0;
max_vert = 0;
for i = 2:rows
    sum = 0;
    for j = 2:cols
        difference = abs(uint32(I(i, j) - I(i, j-1)));
        if(difference > 20)
            sum = sum + difference;
        end
    end
    vert1(i) = sum;
    %-------------------------- Find Peak in Vertical Histogram
    if(sum > maximum)
        max_vert = i;
        maximum = sum;
    end
    total_sum = total_sum + sum;
end
average = total_sum / rows;
%-------Smoothen the Vertical Histogram by applying Low Pass Filter
sum = 0;
vert = vert1;
temp_vert=0;
 if(average<0.3*maximum)%To check if the original image is already isolated or not
     temp_vert=1;
    for i = 21:(rows-21)
        sum = 0;
                    for j = (i-20):(i+20)
                        sum = sum + vert1(j);
                    end
         vert(i) = sum / 40;
    end

    for i = 1:rows
        if(vert(i) < average)
            vert(i) = 0;
            for j = 1:cols
                I(i, j) = 0;
            end
        end
    end
 else
     I=I1;
 end

if(temp_vert==1)
    %--------------------- Find Probable candidates for Number Plate
    %--------Saving the co-ordinates having non-zero horizontal histogram
    j = 1;
    for i = 2:cols-2
        if(horz(i) ~= 0 && horz(i-1) == 0 && horz(i+1) == 0)
            column(j) = i;
            column(j+1) = i;
            j = j + 2;
        elseif((horz(i) ~= 0 && horz(i-1) == 0) || (horz(i) ~= 0 && horz(i+1) == 0))
            column(j) = i;
            j = j+1;
        end
    end

    %--------------Saving the co-ordinates having non-zero vertical histogram
    j = 1;
    for i = 2:rows-2
        if(vert(i) ~= 0 && vert(i-1) == 0 && vert(i+1) == 0)
            row(j) = i;
            row(j+1) = i;
            j = j + 2;
        elseif((vert(i) ~= 0 && vert(i-1) == 0) || (vert(i) ~= 0 && vert(i+1) == 0))
            row(j) = i;
            j = j+1;
        end
    end

    %-------------------converting odd size of the hori co-ordinate array to even
    [temp column_size] = size (column);
    if(mod(column_size, 2))
        column(column_size+1) = cols;
    end

    %----------------converting odd size of the verti co-ordinate array to even
    [temp row_size] = size (row);
    if(mod(row_size, 2))
        row(row_size+1) = rows;
    end
    %------------------------------Region of Interest Extraction
    %---------------------------Check each probable candidate
    for i = 1:2:row_size
        for j = 1:2:column_size
    %------------------------If it is not the most probable region remove it from image
            
            if(~((max_horz >= column(j) && max_horz <= column(j+1)) && (max_vert >=row(i) && max_vert <= row(i+1))))
            for m = row(i):row(i+1)
                for n = column(j):column(j+1)
                    I(m, n) = 0;
                end
            end
          end
        end
    end


     
        str=strel('disk', 1);
        I=imerode(I,str);
        % Get all rows and columns where the image is nonzero
        [nonZeroRows nonZeroColumns] = find(I);
        % Get the cropping parameters
        topRow = min(nonZeroRows(:));
        bottomRow = max(nonZeroRows(:));
        leftColumn = min(nonZeroColumns(:));
        rightColumn = max(nonZeroColumns(:));
        % Extract a cropped image from the original.
        croppedImage = I(topRow-5:bottomRow+5, leftColumn-5:rightColumn+5);
        figure, imshow(croppedImage);
        I=croppedImage;
  
end;
 
%----------------------------------------------------------------------------% 
f=I;
f=imresize(f,[400 NaN]);                   %%image loading unit
%no of rows =400 
imshow(f);
%g=rgb2gray(f);
g=medfilt2(f,[3 3]);
%2d median filter
%reduce salt n pepper noise
%the median value in the 3-by-3 neighborhood

%**********************************
stre=strel('disk',1);    %strucruring element
%structuring element of disk shape(matrix of 1) with radius 1

gi=imdilate(g,stre);
%fill gaps

ge=imerode(g,stre);   %%%% morphological image processing
%shrink or eliminate irrelevant details

gdiff=imsubtract(gi,ge);
gdiff=mat2gray(gdiff);
%amin and amax as min and mx value of image

gdiff=conv2(gdiff,[1 1;1 1]); %convolution of matrix   %brighten

gdiff=imadjust(gdiff,[0.5 0.7],[0 1],.1);
%maps intensity from .5-.7 to 0-1 with curve of gamma .1

B=logical(gdiff);%so that size works
[a1 b1]=size(B);
figure(2)
imshow(B)

er=imerode(B,strel('line',100,0));  
%img , structuring element
%removing usless horizontal line
%figure(3)
%imshow(er)
out1=imsubtract(B,er);

F=imfill(out1,'holes');      %filling the object
H=bwmorph(F,'thin',1); 
%morphological operation on binary image
%thin it once  %to remove useless thin lines
H=imerode(H,strel('line',3,90));
%structuring element of lenf=gth 3 and degree 90
%figure(4)
%imshow(H)

final=bwareaopen(H,floor((a1/15)*(b1/15)));  
%remove objects that hav fewer pixels than second argument
%floor rounds the value

final(1:floor(.9*a1),1:2)=1;

final(a1:-1:(a1-20),b1:-1:(b1-2))=1;
%vertically a1 to a1-20 in steps of -1
%white rectangle at bottom left
yyy=template(2);
figure(5)
imshow(final)

Iprops=regionprops(final,'BoundingBox','Image');
% drawing boundary box around 
% each connected component in the image
hold on                            %add to the exsisting graph
for n=1:size(Iprops,1)               %sixe of iprops in first dimension
    rectangle('Position',Iprops(n).BoundingBox,'EdgeColor','g','LineWidth',2); 
end
hold off                            %reset axis properties

NR=cat(1,Iprops.BoundingBox);   %%Data storage section
%string of all boundary boxes, 
%the properties in the string are x coordinate,y coordinate, x width, y width
%white box of size bounding box
%NR= number of regions*4

[r ttb]=connn(NR);
%ttb is array of useful NR


if ~isempty(r)

    xlow=floor(min(reshape(ttb(:,1),1,[])));
    xhigh=ceil(max(reshape(ttb(:,1),1,[])));
    xadd=ceil(ttb(size(ttb,1),3));
    ylow=floor(min(reshape(ttb(:,2),1,[])));    %%%%%area selection
    yadd=ceil(max(reshape(ttb(:,4),1,[])));
    final1=H(ylow:(ylow+yadd+(floor(max(reshape(ttb(:,2),1,[])))-ylow)),xlow:(xhigh+xadd));
    [a2 b2]=size(final1);
    final1=bwareaopen(final1,floor((a2/20)*(b2/20)));
    figure(6)
    imshow(final1)

    
   
    Iprops1=regionprops(final1,'BoundingBox','Image');
    NR3=cat(1,Iprops1.BoundingBox);
  
    
    ydim=(reshape(NR3(:,2),1,[]));
    siz=size(ydim);
    
    
%     j=0;
%     for i=1:siz(2)
%         if(abs(ydim(i)-t)<10)
%             j=j+1;
%         end;
%     end;
%     if(j==1)
%         for i=1:siz(2)
%             if(ydim(i)==t)
%                 ydim(i)=401;
%                 t=min(ydim);
%             end;
%         end;
% %     end;
%     j=1;


%     NR3
%     while(j<siz(2))
%         t=min(ydim);        
%         for i=1:siz(2)
%             if(abs(ydim(i)-t) <10)
%                 sort1(j)=i;
%                 ydim(i)=401;
%                 j=j+1;
%             end;
%         end;
%     end;
%     sort1
%     temp=Iprops1;
%     for i=1:siz(2)
%         Iprops1(i,:)=temp(ceil(sort1(i)),:);
%     end;
   
 
    I1={Iprops1.Image};
    %%
  
    carnum=[];
    if (size(NR3,1)>size(ttb,1))
        [r2 to]=connn2(NR3);
        
        for i=1:size(Iprops1,1)
            
            ff=find(i==r2);
            if ~isempty(ff)
                N1=I1{1,i};
                letter=readLetter(N1,2);
            else
                N1=I1{1,i};
                letter=readLetter(N1,1);
            end
            if ~isempty(letter)
                carnum=[carnum letter];
            end
        end
    else
        for i=1:size(Iprops1,1)
            N1=I1{1,i};
            letter=readLetter(N1,1);
            carnum=[carnum letter];
        end
    end
    
    %%
    
    fid1 = fopen('carnum.txt', 'wt');
    fprintf(fid1,'%s',carnum);
    fclose(fid1);
    winopen('carnum.txt')
   


else
    fprintf('license plate recognition failure\n');
    fprintf('Characters are not clear \n');
end