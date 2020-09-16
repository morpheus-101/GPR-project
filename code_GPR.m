%code for generating Bscans
clear all;
close all; 
M = 76;
XT=ones(1601,M); 
for i = 1:M
filename= sprintf('%d.CSV',i); 
Ascan = csvread(filename,3,0); 
Mag = sqrt(Ascan(:,2).^2 + Ascan(:,3).^2); 
if (Ascan(:,2)~=0) 
phase = tanh(Ascan(:,3)./ Ascan(:,2)); 
else
phase = +1;
end
rx_sig = Mag.*exp(1i.*phase); 
XT(:,i) = Mag; 
end
time = Ascan(:,1);
%--------------------------------
%setting the scale for the depth axis
[rowsf, columnsf] = size(time);
depth = zeros(rowsf,1);
for i= 1 : rowsf
    depth(i,1) = ((3*10^8)/sqrt(5))*(time(i,1)/2)*(10^3);
end
depthMax = max(depth);
%--------------------------------
%setting scale for the column axis
[rowsx,columnsx] = size(XT);
columnsb = zeros;
for p = 1 : columnsx
    columnsb(1,p) = p;
end
%--------------------------------
Mod_XT = zeros(1601,M);
for i = 1 : rowsx
    for j = 1 : columnsx
        if(i < 1435 && j < columnsx)
           Mod_XT(i,j) = XT(i,j);
        end
    end
end

figure(500);
imagesc(columnsb, depth ,XT);
colorbar
xlabel('Column number')
ylabel('Depth(mm)')
title('Bscan')

figure(111);
%Mod_XT = imnoise(Mod_XT,'gaussian',0,0.0000010);
imagesc(columnsb, depth ,Mod_XT);
colorbar
xlabel('Column number')
ylabel('Depth(mm)')
title('Bscan with noise')
%--------------------------------------------------------------------------------------------------------------------------------------------------
%average subtraction code
imwrite (mat2gray(Mod_XT), 'NEW IMAGE.png');
[XT2,map] = imread ('NEW IMAGE.png');
[rows,columns] = size(XT2);
avg=0;
XTnew = zeros;
for i = 1: rows
    sum=0;
    for j = 1: columns
        for k= -2:2
             z=i+k;
            if(z>0 && z<=rows)    
               sum = sum + XT2(i,j);
            end
        end
       XTnew(i,j)=XT2(i,j)-(1/5).*sum;
    end
end
colorbar;
figure (2);
imagesc(columnsb ,depth ,XTnew); 
xlabel('Column number')
ylabel('Depth(mm)')
title('Bscan clutter reduced')
maxclutter = max(XTnew);
%-------------------------------------------------------------------------------------------------------------------------
% Mean filter
MF=ones(7,7)/49;
MF_scan = imfilter(XTnew,MF);
figure(3);
imagesc(columnsb,depth ,MF_scan); 
colorbar
xlabel('Column number')
ylabel('Depth(mm)')
title('Bscan after mean filtering')
%----------------------------------------------------------------------------------------------------------------------------
% Gaussian filter
q2=conv2(XTnew,MF);
GF = fspecial('gaussian', [3 3], 0.5);
GF_scan = imfilter(q2,GF);
axis tight;
figure(4);
imagesc(columnsb,depth,GF_scan);
colorbar
xlabel('Column number')
ylabel('Depth(mm)')
title('Bscan after applying Gaussian filter')
%--------------------------------------------------------------------------------------------------------------------------------
% Erosion technique 
SE = strel('disk',15,0);
XTerode = imerode(GF_scan,SE);
figure(5);
imagesc(columnsb ,depth ,XTerode);
xlabel('Column number')
ylabel('Depth(mm)')
title('Bscan with erosion')
%------------------------------------------------------------------------------------------------------------------------------
%Thresholding operation
XTint= zeros ( rowsx, columnsx);
for i = 1: rowsx
    for j = 1: columnsx
        if XTerode(i,j)>= 110
            XTint(i,j) = XTerode(i,j);
        end
    end
end
figure(6);
imagesc(columnsb ,depth, XTint);
xlabel('Column number')
ylabel('Depth(mm)')
title('Bscan after thresholding')

%-----------------------------------------------------------------------------------------------------------------------------------
%Energy projection plots
figure(7)
plot(columnsb, XTint);
xlabel('Column number')
ylabel('Intensity')
figure(8)
plot(depth, XTint);
xlabel('Depth(mm)')
ylabel('Intensity')
%---------------------------------------------------------------------------------------------------------------------------------------------------
%Setting thresholds for estimating the number of objects
%taking the transpose of XTint and generating the pos array
XTintt= XTint.';
[rowst , columnst] = size(XTintt);
pos = zeros(1,columnst);
for i= 1:columnst
    for j = 1: rowst
        if XTintt (j,i) ~= 0
            pos (1,i) =  i;
        end
    end
end
%--------------------------------------------------------------------------
%counting the number pixels which represent the mine from the pos array
ngroups = zeros(1, columnst);
k=1;
count=0;
for i= 1 : columnst
    while k < columnst
        if pos(1,k) ~= 0
            count = count+1;
            k=k+1;
            
        else
            k = k+1;
            ngroups(1,i) = count;
            count=0;
            break;
        end
    end
end
%--------------------------------------------------------------------------
%generating posw array similar to pos array
posw = zeros(1,columnsx + 1);
for i= 1:columnsx
    for j = 1: rowsx
        if XTint (j,i) ~= 0
            posw (1,i) =  i;
        end
    end
end
%--------------------------------------------------------------------------
%generating an array ngw which is similar to ngroups array 
ngw = zeros(1, columnsx);
kw=1;
countw=0;
for i= 1:columnsx
    while kw < columnsx + 2
        if posw(1,kw) ~= 0
            countw = countw+1;
            kw=kw+1;
           
        else
            kw = kw+1;
            ngw(1,i) = countw;
            countw=0;
            break;
        end
    end
end
%--------------------------------------------------------------------------
%considering the X-Z(column number v/s intensity plot)
tolerance1 = 8;
nobj = 0;
for i = 1:columnsx
    if ngw(1,i) >= tolerance1 
        nobj = nobj + 1;
    end
end
%----------------------------------------------------------------------------------------------------------------------------------------
% extracting individual objects

if nobj == 0
    fprintf('No target is present in the scanned area\n');
%--------------------------------------------------------------------------    
elseif nobj == 1
    fprintf('1 target is present in the scanned area\n');
    maxvaltar = zeros (rowsx , 1);
for i = 1:rowsx
    for j = 1: columnsx
        if maxvaltar(i,1)< XTint(i,j)
            maxvaltar(i,1) = XTint(i,j);
        end
    end
end
    
    [~,tarposrow]= max(maxvaltar);
    [~,tarposcol]= max(max(XTint));
    
    coordinateX= ((3*10^8)/sqrt(5))*(time(tarposrow,1)/2)*(10^2); 
    coordinateY = tarposcol *0.5;
    fprintf("Co-ordinates of the target = [%f,%f]\n",coordinateX,coordinateY);
    
    depthtarg= ((3*10^8)/sqrt(5))*(time((tarposrow-45),1)/2)*(10^3);
    fprintf('Depth of the target =%f\n', depthtarg);
    
    ET = zeros(rows,columnsx);
    for i= (tarposrow - 90):(tarposrow + 90) 
        %for j = (tarposcol - (tarposcol-3)) : (tarposcol +(tarposcol-3))
        for j = 3 : 68
            ET (i,j) = XTint (i,j);
        end
    end
    figure(9);
    imagesc(columnsb, depth, ET);
    xlabel('Column No.');
    ylabel('Depth(mm)');
    title ('Extracted signature of target');
%--------------------------------------------------------------------------    
elseif nobj == 2
     fprintf('2 targets are present in the scanned area\n');
     maxvaltar = zeros (rowsx , 1);
for i = 1:rowsx
    for j = 1: columnsx
        if maxvaltar(i,1)< XTint(i,j)
            maxvaltar(i,1) = XTint(i,j);
        end
    end
end
    
    [~,tarposrow]= max(maxvaltar);
    [~,tarposcol]= max(max(XTint));
    
    ET = zeros(rows,columnsx);
    for i= (tarposrow - 150):(tarposrow + 125) 
        for j = (tarposcol - 30) : (tarposcol + 8)
            ET (i,j) = XTint (i,j);
        end
    end 
    
   ET2 = XTint - ET;
   
   maxval1st = zeros(rowsx,1);
 for i = 1:rowsx
    for j = 1: columnsx
        if maxval1st(i,1)< ET(i,j)
            maxval1st(i,1) = ET(i,j);
        end
    end
end
    
    [inttar1st,tarposrow1st]= max(maxval1st);
    [~,tarposcol1st] = max(max(ET));
    
    coordinateX1= ((3*10^8)/sqrt(5))*(time(tarposrow1st,1)/2)*(10^2);
    coordinateY1 = tarposcol1st *0.5;
    fprintf("Co-ordinates of the target = [%f,%f]\n",coordinateX1,coordinateY1);
    
    depth1st = ((3*10^8)/sqrt(5))*(time((tarposrow1st-45),1)/2)*(10^3);
    fprintf('Depth of the target 1 =%fmm\n', depth1st);
    
    figure(10);
    imagesc(columnsb, depth, ET);
    xlabel('Column No.');
    ylabel('Depth(mm)');
    title ('Extracted signature of target 1');
    
   maxval2nd = zeros(rowsx,1);
 for i = 1:rowsx
    for j = 1: columnsx
        if maxval2nd(i,1)< ET2(i,j)
            maxval2nd(i,1) = ET2(i,j);
        end
    end
end
    
    [inttar2nd,tarposrow2nd]= max(maxval2nd);
    [~,tarposcol2nd] = max(max(ET2));
    
    coordinateX2= ((3*10^8)/sqrt(5))*(time(tarposrow2nd,1)/2)*(10^2); 
    coordinateY2 = tarposcol2nd *0.5;
    fprintf("Co-ordinates of the target = [%f,%f]\n",coordinateX2,coordinateY2);
    
    depth2nd = ((3*10^8)/sqrt(5))*((time((tarposrow2nd-45),1))/2)*(10^3);
    fprintf('Depth of the target 2 =%f mm\n', depth2nd);
    
    figure(11);
    imagesc(columnsb, depth, ET2);
    xlabel('Column No.');
    ylabel('Depth(mm)');
    title ('Extracted signature of target 2');
    separation = ((tarposcol2nd-10)*0.5)- ((tarposcol1st+10)*0.5);
    fprintf("Separation between the two targets is %f cm\n",separation);
end
%-----------------------------------------------------------------------------------------------------------------
% finding out the height
height= zeros(nobj,1);
k=1;
for i = 1: columnsx
    if ngw(1, i) ~=0
        height(k,1) = ngw(1,i)* 0.5;
        k = k+1;
    end
end
for j = 1: nobj
    fprintf('Width of the target %d =%f cm\n',j,height(j,1));
end
%--------------------------------------------------------------------------------------------------------------------------------------------------------
% post processing
  if nobj == 0
         
  elseif nobj == 1
        GF = imgaussfilt(ET);
        figure(12);
        imagesc(columnsb, depth, GF);
        xlabel('Column No.');
        ylabel('Depth(mm)');
        title ('Gaussian filtered data');
        
        summ = 0;
        maxcol = max (GF);
        [rowsg, columnsg] = size(GF);
        for n = 1: columnsg
            summ = maxcol(1,n) + summ;
        end
        maxintensity = summ / columnsg;
        threshold50 = (3*maxintensity)/4;
        TM = zeros(rowsg,columnsg);
        for m = 1: rowsg
            for j = 1:columnsg
                if GF(m,j) >= threshold50
                    TM(m,j) = GF(m,j);
                end
            end
        end
        figure(13);
        imagesc(columnsb, depth, TM);
        xlabel('Column No.');
        ylabel('Depth(mm)');
        title ('Thresholded data');
%-------------------------------------------------------------------------------------------------------------         
  elseif nobj == 2
        GF1 = imgaussfilt(ET);
        figure(14);
        imagesc(columnsb, depth, GF1);
        xlabel('Column No.');
        ylabel('Depth(mm)');
        title ('Gaussian filtered data for target 1');
        
        summ = 0;
        maxcol = max (GF1);
        [rowsg, columnsg] = size(GF1);
        for n = 1: columnsg
            summ = maxcol(1,n) + summ;
        end
        maxintensity = summ / columnsg;
        threshold50 = (3*maxintensity)/4;
        TM = zeros(rowsg, columnsg);
        for m = 1: rowsg
            for j = 1:columnsg
                if GF1(m,j) >= threshold50
                    TM(m,j) = GF1(m,j);
                end
            end
        end
        figure(15);
        imagesc(columnsb, depth, TM);
        xlabel('Column No.');
        ylabel('Depth(mm)');
        title ('Thresholded data for target 1');
        
        
        GF2 = imgaussfilt(ET2);
        figure(16);
        imagesc(columnsb, depth, GF2);
        xlabel('Column No.');
        ylabel('Depth(mm)');
        title ('Gaussian filtered data for target 2');
        
        maxcol2 = max (GF2);
        summ1=0;
        [rowsg2, columnsg2] = size(GF2);
        for n = 1: columnsg2
            summ1 = maxcol2(1,n) + summ1;
        end
        maxintensity2 = summ1 / columnsg2;
        threshold502 = (3*maxintensity2)/4;
        TM2 = zeros(rowsg2,columnsg2);
        for m = 1: rowsg2
            for j = 1:columnsg2
                if GF2(m,j) >= threshold502
                    TM2(m,j) = GF2(m,j);
                end
            end
        end
        figure(17);
        imagesc(columnsb, depth, TM2);
        xlabel('Column No.');
        ylabel('Depth(mm)');
        title ('Thresholded data for target 2');
  end
%--------------------------------------------------------------------------------------------------------------------------  
%PCA
if(nobj == 1)
    [m, n]=size (GF); 
    mn = mean(GF,2);
    X=GF - repmat(mn,1,n);
    Z = 1/sqrt(n-1) * X';
    covZ = Z' * Z; 
    [eigvec, eigval] = eig(covZ);
    [U,S,V]=svd(covZ);
    variances=diag (S).*diag(S); 
    figure(18);
    bar (variances (1:50)); 
    title('principle components of the target')
elseif(nobj == 2)
    [m, n]=size (GF1); 
    mn = mean(GF1,2);
    X=GF1 - repmat(mn,1,n);
    Z = 1/sqrt(n-1) * X';
    covZ = Z' * Z; 
    [eigvec, eigval] = eig(covZ);
    [U,S,V]=svd(covZ);
    variances=diag (S).*diag(S); 
    figure(19);
    bar (variances (1:50)); 
    title('principle components of target 1')
    
    [m, n]=size (GF2); 
    mn = mean(GF2,2);
    X=GF2 - repmat(mn,1,n);
    Z = 1/sqrt(n-1) * X';
    covZ = Z' * Z; 
    [eigvec, eigval] = eig(covZ);
    [U,S,V]=svd(covZ);
    variances=diag (S).*diag(S); 
    figure(20);
    bar (variances (1:50)); 
    title('priniple components of target 2')
end

%-----------------------------------------------------------------------------------------------------------------------------------------
%generating mean of specific ROI 
total=0;
denominator=0;
for i = (tarposrow-10) : (tarposrow+10)
    for j = (tarposcol-5) : (tarposcol+5)
        if(XTint(i,j)>0)
            total = total + XTint(i,j);
            denominator = denominator + 1;
        end
    end
end

ROImean = total/denominator;


