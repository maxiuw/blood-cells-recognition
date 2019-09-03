clear all; close all; clc;
NumMec=79186;

listaF=dir('../svpi2017_TP2_img_*.png'); % read the list of images

txt=fopen('tp2_79186.txt','wt'); % open txt

for v=1:max(size(listaF)) %loop to read the images form the list
    
    [X,map]=imread(strcat('../',listaF(v).name)); %read img and put them in the order
    img=im2double(imread(strcat('../',listaF(v).name)));
    temp=listaF(v).name(22);
    
    if      temp=='0'
        NumImg=listaF(v).name(23); % img nr
        NumSec=listaF(v).name(18:20); % seq nr
    else
        NumImg=listaF(v).name(22:23); % img nr
        NumSec=listaF(v).name(18:20); % seq nr
    end
    
    
    [m,n,k] = size(img);
    % calling the red channel
    b_red = img(:,:,1);
    % calling the Green Channel
    b_green = img(:,:,2);
    % calling the blue channel
    b_blue = img(:,:,3);
    % creating negative for red channel
    % Creating a null array
    b_red_neg = zeros(m,n);
    for i = 1:m
        for j = 1:n
            b_red_neg(i,j) = 1 - b_red(i,j);
        end
    end
    % creating negative for green channel
    % Creating a null array
    b_green_neg = zeros(m,n);
    for i = 1:m
        for j = 1:n
            b_green_neg(i,j) = 1 - b_green(i,j);
        end
    end
    % creating negative for blue channel
    % Creating a null array
    b_blue_neg = zeros(m,n);
    for i = 1:m
        for j = 1:n
            b_blue_neg(i,j) = 1 - b_blue(i,j);
        end
    end
    % combining all the "Negative" channel to form the Negative of a given RGB
    % or Color Image
    % Creating a Null array
    d = zeros(m,n,k);
    d(:,:,1) = b_red_neg;
    d(:,:,2) = b_green_neg;
    d(:,:,3) = b_blue_neg;
    %figure(); imshow(d);
    figure();
    subplot(2,2,1); imhist(b_red_neg);
    subplot(2,2,2); imhist(b_green_neg);
    subplot(2,2,3); imhist(b_blue_neg);
    %  HSV=rgb2hsv(img);
    %  A=rgb2hsv(img);
    %  B=imbinarize(A,0.2);
    % B = im2bw(C, 0.8);  %%binarization
    % C=imfill(d, 'holes');
    B = im2bw(d, 0.7);  %%binarization
    C=imcomplement(B);
    F1=[0 1 1 0; 1 1 1 1; 1 1 1 1;0 1 1 0];
    %F2=[0 0 0; 0 -1 1;0 1 0];
    filter_1=filter2(F1,C);
    %filter_2=filter2(F2,filter_1);
    %filter_3=filter2(F2,filter_2);
    E=imfill(filter_1, 'holes');
    [L n]=bwlabel(E);
    s=regionprops(L, 'All'); 
    T=[s.Area]; %pole
    P=[s.Perimeter]; %obwod
    ff=4*pi*T./P./P; %form factor
    idxA = find(T > 200); %find the regions with more than 100px
    regions= ismember(L,idxA); %regions with more than 100px
   % F=imclose(regions,ones(1)); %erase singlepixels
%     F_1=imfill(F, 'holes');
    filter_2=filter2(F1,regions);
    F=imfill(filter_2, 'holes');
   % F_1=imclose(F,ones(4)); %erase singlepixels
    %  N=[0 1 0; 1 1 1; 0 1 0];
    %  D=imopen(C,N);
    %  L= bwlabel(C);
    %  S=regionprops(L);
    
      figure(); imshow(filter_2);
    
    
end








% - binarize the picture
% - detect what is the values of the each shape like ff itp %1?
% - show all the objects with the values
% -
