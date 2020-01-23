function NumMec= tp2_79168

%SVPI_processTP2_2018('tp2_79186.m')
clear all; close all; clc;
NumMec=79186;

listaF=dir('../svpi2018_TP2_img_*.png'); % read the list of images

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
    
    
    ObjBord=0;% objects which touching the border of the image
    ObjOK=0; %number of the object which are ok
    acanthocyte=0;
    degmacyte=0;
    echinocyte=0;
    eliptocyte=0;
    mega =0;
    normo=0;
    schistocyte=0;
    sickle=0;
    sphere=0;
    stomato=0;
    target=0;
    teardrop=0;
    
    
    % making better contrast and emphasising edges,better for separation
    I1=rgb2gray(img);
    hy = fspecial('prewitt'); %or prewitt/laplacian/sobel
    hx = hy';
    Iy = imfilter(I1, hy, 'symmetric');
    Ix = imfilter(I1, hx, 'symmetric');%replicate
    gradmag = sqrt(Ix.^2 + Iy.^2);
    %gradmag=imbinarize(gradmag,'adaptive');
    
    [E51, tresh] = edge(gradmag,'log', 0.001);
    E=imclose(E51,ones(3)); 
    E_1=imfill(E,'holes');
    E3 = bwareaopen(E_1, 300); % obj bigger than 350px will disappear
    E4=imclearborder(E3); %clearing boarders
    
    seD = strel('diamond',1); %finishing, making edges and shape sharper
    Efinal = imerode(E4,seD);
    Efinal = imerode(Efinal,seD);
    Efinal=bwareaopen(Efinal, 250);
         figure();imshow(Efinal);
    
    %obtaining more properties like saturation ect
    HSV=rgb2hsv(img);
    H=HSV(:,:,1);
    S=HSV(:,:,2);
    V=HSV(:,:,3);
    [L,N]=bwlabel(Efinal);
    % s=regionprops(L);
    
    for k=1:N
        M= L ==k;
        %  MS=S*0;
        %     MS(M)=S(M);
        % figure;imshow((H(M)))
        %inv1=mean(invmoments(H(M)));
        hue=mean(H(M));
        sat=mean(S(M));
        intens=mean(V(M));
        
        
        K(k,1)=hue;
        K(k,2)=sat;
        K(k,3)=intens;%vals de S saturation of alll the pixels, H,V
        %K(k,4)=inv1/1.0e+09 ;
        %         text(s(k).Centroid(1),s(k).Centroid(2),...
        %              sprintf(' %s \n %s \n %s',hue,sat,intens), 'Color',[0 0 1], 'FontSize', 10)
    end
    
    %    borders objects
    se = strel('disk',4);
    b1=im2bw((gradmag),0.13);
    b2=imclose(b1,se);
    b3=imfill(b2,'holes');
    b4=imclearborder(b3);
    bw=b3-b4;
    bw_a = padarray(bw,[1 1],1,'pre');
    bw_a_filled = imfill(bw_a,'holes');
    bw_a_filled = bw_a_filled(2:end,2:end);    
    bw_b = padarray(padarray(bw,[1 0],1,'pre'),[0 1],1,'post');
    bw_b_filled = imfill(bw_b,'holes');
    bw_b_filled = bw_b_filled(2:end,1:end-1);    
    bw_c = padarray(bw,[1 1],1,'post');
    bw_c_filled = imfill(bw_c,'holes');
    bw_c_filled = bw_c_filled(1:end-1,1:end-1);    
    bw_d = padarray(padarray(bw,[1 0],1,'post'),[0 1],1,'pre');
    bw_d_filled = imfill(bw_d,'holes');
    bw_d_filled = bw_d_filled(1:end-1,2:end);
    bw_filled = bw_a_filled | bw_b_filled | bw_c_filled | bw_d_filled;
    
    b = regionprops(bw_filled, 'All');
    ObjBord=numel(b);
    
    
    %     figure();subplot(2,2,1); imshow(I1);title('grayimg');
    %     subplot(2,2,2);imshow(gradmag,[]), title('Gradient magnitude (gradmag)');
    %     subplot(2,2,3); imshow(Efinal);title('final');
    %     subplot(2,2,4);imshow(bw_filled); title('Final filled result');
    
    
    s = regionprops(Efinal, 'All');
    %Obj = (Efinal == 1);   % 1 is the label number of the first object.
    %     figure, imshow(Obj); title(sprintf('%s %s', NumSec,NumImg));
    
    # 
    for k = 1:numel(s)
        if verLessThan('matlab','8.3') %test Matlab(8.3->2014a)
            per = s(k).PerimeterOld;
        else %newer versions of Matlab
            per = s(k).Perimeter;
        end
        centr = s(k).Centroid;
        solid= s(k).Solidity;  %solidity vector
        area=s(k).Area;     %vector area
        per=s(k).Perimeter; %vector perimeter
        ecc=s(k).Eccentricity;
        s(k).FormFactor=4*pi*area./per./per;  %form factor
        ff=s(k).FormFactor;
        s(k).inv=invmoments(s(k).Image);
        text(centr(1), centr(2), sprintf('%d ',...
            k), ...
            'HorizontalAlignment', 'right', ...
            'VerticalAlignment', 'top','color', '[0.6 0.8 1]');
        
    end
    
    % creating db of the properties of the blood cells
    % during this task, we could not use the outside database 
    % it did not have the access to the template of the bloodcells
     %solidity  %ecc    %ff %3 columns od H,S,V  %abs log(hu moment)/10
    T(:,:,5)=[
        0.8314    0.4134    0.6008    0.4742    0.1875    0.6801    0.1785
        0.7208    0.8074    0.4848    0.5855    0.2020    0.5586    0.1529
        0.9215    0.5170    0.8600    0.3474    0.1963    0.7117    0.1815
        0.9768    0.6673    0.9708    0.3211    0.1757    0.7237    0.1788
        0.9866    0.3725    1.0224    0.0404    0.1869    0.7268    0.1835
        0.8677    0.5919    0.6945    0.6371    0.3580    0.7645    0.1764
        0.8940    0.5829    0.7881    0.0708    0.1986    0.7536    0.1796
        0.8917    0.5762    0.7741    0.0541    0.2007    0.7522    0.1791
        0.8695    0.4878    0.7212    0.0504    0.2077    0.7556    0.1783
        0.9758    0.8913    0.8309    0.4523    0.2059    0.7463    0.1544
        0.9798    0.8532    0.8866    0.4333    0.2097    0.7211    0.1636
        0.9573    0.9253    0.7366    0.7935    0.1673    0.6933    0.1416
        0.9889    0.5774    1.0143    0.7522    0.1857    0.6902    0.1817
        0.9854    0.4102    1.0188    0.6136    0.1768    0.7099    0.1833
        0.9896    0.6414    0.9971    0.8568    0.1792    0.6984    0.1802
        0.9850    0.2455    1.0375    0.0381    0.1978    0.7173    0.1837
        0.9825    0.3557    1.0323    0.0427    0.2010    0.7196    0.1835
        0.9891    0.3016    1.0473    0.0345    0.2058    0.6930    0.1836
        0.7975    0.7800    0.6145    0.6151    0.2109    0.6970    0.1573
        0.9332    0.7002    0.8115    0.7972    0.2194    0.6184    0.1721
        0.8019    0.6936    0.6087    0.7789    0.2097    0.6083    0.1618
        0.9834    0.5452    0.9989    0.0655    0.1953    0.7659    0.1818
        0.5317    0.9702    0.2795    0.6982    0.2653    0.5731    0.0418
        0.8024    0.9808    0.3749    0.7441    0.2755    0.6513    0.0702
        0.9820    0.3212    1.0470    0.8129    0.2538    0.6111    0.1835
        0.9842    0.4363    1.0523    0.7584    0.2480    0.5835    0.1831
        0.9775    0.3311    1.0188    0.8169    0.2585    0.5942    0.1834
        0.9882    0.2455    1.0405    0.5845    0.2324    0.7106    0.1837
        0.9862    0.5783    1.0155    0.5570    0.2307    0.7402    0.1816
        0.9871    0.3077    1.0363    0.5104    0.2261    0.7365    0.1836
        0.9099    0.9589    0.4831    0.0411    0.1770    0.6744    0.1102
        0.9487    0.2057    0.9265    0.7977    0.5480    0.6659    0.1824
        0.9844    0.2970    1.0237    0.0454    0.1769    0.6952    0.1836
        0.9623    0.8498    0.8007    0.6931    0.1957    0.6983    0.1617
        0.9612    0.7096    0.8527    0.7589    0.2809    0.7628    0.1735
        0.9753    0.6470    0.9811    0.7993    0.2290    0.5880    0.1781];
    T(:,:,4)=[
        0.8535    0.3871    0.6161    0.4237    0.1675    0.6075    0.1790
        0.7573    0.7908    0.5269    0.4943    0.1705    0.4716    0.1564
        0.9211    0.5076    0.8559    0.3160    0.1786    0.6473    0.1816
        0.9780    0.6551    0.9618    0.2967    0.1624    0.6687    0.1793
        0.9874    0.3657    1.0232    0.0380    0.1758    0.6836    0.1835
        0.8839    0.5688    0.7111    0.5754    0.3233    0.6904    0.1770
        0.9032    0.5722    0.7974    0.0644    0.1807    0.6855    0.1799
        0.9047    0.5649    0.7907    0.0492    0.1823    0.6835    0.1796
        0.8781    0.4821    0.7317    0.0463    0.1908    0.6939    0.1786
        0.9778    0.8821    0.8466    0.4136    0.1882    0.6824    0.1570
        0.9814    0.8473    0.8957    0.3986    0.1930    0.6633    0.1646
        0.9623    0.9189    0.7566    0.7204    0.1519    0.6294    0.1445
        0.9895    0.5689    1.0154    0.7107    0.1754    0.6521    0.1818
        0.9862    0.4030    1.0198    0.5777    0.1665    0.6684    0.1834
        0.9892    0.6369    0.9957    0.8159    0.1706    0.6651    0.1804
        0.9861    0.2380    1.0372    0.0353    0.1830    0.6635    0.1837
        0.9838    0.3483    1.0333    0.0395    0.1856    0.6646    0.1836
        0.9900    0.2931    1.0483    0.0317    0.1889    0.6361    0.1836
        0.8147    0.7648    0.6348    0.5494    0.1884    0.6226    0.1589
        0.9385    0.6973    0.8230    0.7223    0.1988    0.5603    0.1724
        0.8280    0.6764    0.6345    0.6992    0.1882    0.5461    0.1638
        0.9843    0.5386    1.0000    0.0617    0.1839    0.7210    0.1819
        0.5810    0.9677    0.3180    0.5764    0.2191    0.4731    0.0536
        0.8344    0.9781    0.4129    0.6325    0.2342    0.5536    0.0788
        0.9836    0.3113    1.0471    0.7385    0.2305    0.5552    0.1836
        0.9857    0.4242    1.0516    0.6816    0.2229    0.5245    0.1832
        0.9794    0.3330    1.0179    0.7474    0.2365    0.5436    0.1834
        0.9891    0.2428    1.0390    0.5394    0.2145    0.6558    0.1837
        0.9871    0.5689    1.0160    0.5207    0.2157    0.6920    0.1818
        0.9880    0.3008    1.0363    0.4739    0.2099    0.6839    0.1836
        0.9232    0.9549    0.5183    0.0356    0.1534    0.5844    0.1145
        0.9522    0.1880    0.9333    0.7461    0.5126    0.6229    0.1826
        0.9854    0.2915    1.0240    0.0425    0.1654    0.6498    0.1836
        0.9603    0.8445    0.8093    0.6378    0.1801    0.6425    0.1626
        0.9636    0.7081    0.8623    0.7104    0.2629    0.7140    0.1737
        0.9781    0.6339    0.9895    0.7058    0.2023    0.5192    0.1786];
    T(:,:,3)=[
        0.8079    0.4599    0.5680    0.5399    0.2115    2.3044    0.1785
        0.6865    0.8128    0.4553    0.7140    0.2436    2.0201    0.1524
        0.9020    0.5115    0.8295    0.3836    0.2154    2.3464    0.1815
        0.9687    0.6697    0.9601    0.3479    0.1881    2.3256    0.1785
        0.9437    0.2502    0.9265    0.8474    0.5822    2.1229    0.1824
        0.8653    0.6203    0.7084    0.6926    0.3904    2.5041    0.1755
        0.8968    0.5802    0.8001    0.0783    0.2186    2.4914    0.1798
        0.8995    0.5761    0.8007    0.0568    0.2188    2.4620    0.1792
        0.8510    0.4991    0.7028    0.0555    0.2263    2.4697    0.1779
        0.9761    0.8907    0.8318    0.4842    0.2199    2.3932    0.1548
        0.9739    0.8676    0.8673    0.4677    0.2281    2.3533    0.1606
        0.9520    0.9347    0.7085    0.8831    0.1859    2.3164    0.1363
        0.9859    0.5802    1.0106    0.7975    0.1967    2.1932    0.1817
        0.9857    0.3977    1.0242    0.6449    0.1862    2.2431    0.1834
        0.9853    0.6440    0.9894    0.8943    0.1870    2.1869    0.1801
        0.9776    0.2418    1.0236    0.0410    0.2127    2.3139    0.1837
        0.9811    0.3452    1.0310    0.0462    0.2171    2.3326    0.1835
        0.9813    0.2768    1.0384    0.0374    0.2230    2.2545    0.1836
        0.7705    0.8017    0.5833    0.6929    0.2361    2.3408    0.1561
        0.9489    0.7207    0.8032    0.8846    0.2439    2.0572    0.1710
        0.7927    0.7063    0.5938    0.8662    0.2329    2.0233    0.1606
        0.9088    0.9618    0.4464    0.0477    0.2053    2.3500    0.1073
        0.4615    0.9737    0.2250    0.8882    0.3381    2.1892    0.0231
        0.7565    0.9842    0.3232    0.9086    0.3363    2.3863    0.0572
        0.9870    0.3353    1.0548    0.8926    0.2791    2.0173    0.1835
        0.9764    0.4240    1.0452    0.8395    0.2746    1.9378    0.1832
        0.9780    0.3465    1.0090    0.8925    0.2833    1.9552    0.1831
        0.9844    0.2358    1.0390    0.6355    0.2529    2.3197    0.1837
        0.9786    0.5861    1.0027    0.5943    0.2467    2.3748    0.1815
        0.9838    0.2870    1.0313    0.5486    0.2437    2.3826    0.1836
        0.9860    0.5509    1.0003    0.0693    0.2076    2.4432    0.1815
        0.9869    0.3611    1.0273    0.0428    0.1979    2.3097    0.1835
        0.9828    0.2883    1.0237    0.0464    0.1884    2.2216    0.1836
        0.9666    0.8550    0.7956    0.7516    0.2127    2.2782    0.1606
        0.9605    0.7169    0.8513    0.8107    0.3002    2.4435    0.1734
        0.9625    0.6461    0.9617    0.9058    0.2597    1.9999    0.1780];
    
    T(:,:,2)=[
        0.8270    0.4637    0.6050    0.5325    0.2104    1.5266    0.1786
        0.6775    0.8383    0.4486    0.7089    0.2440    1.3458    0.1491
        0.8893    0.5106    0.7940    0.3850    0.2162    1.5673    0.1811
        0.9664    0.6732    0.9514    0.3448    0.1879    1.5491    0.1784
        0.9446    0.2165    0.9227    0.8508    0.5847    1.4205    0.1823
        0.8644    0.6296    0.7035    0.7013    0.3947    1.6905    0.1754
        0.8849    0.5837    0.7785    0.0756    0.2182    1.6571    0.1796
        0.8906    0.5847    0.7714    0.0584    0.2189    1.6427    0.1790
        0.8582    0.4846    0.7094    0.0538    0.2267    1.6490    0.1782
        0.9717    0.8987    0.8144    0.4949    0.2242    1.6261    0.1522
        0.9686    0.8653    0.8622    0.4725    0.2287    1.5738    0.1611
        0.9431    0.9326    0.7173    0.8746    0.1843    1.5282    0.1375
        0.9850    0.5788    1.0073    0.7933    0.1962    1.4595    0.1817
        0.9799    0.4114    1.0145    0.6486    0.1863    1.4956    0.1833
        0.9862    0.6481    0.9887    0.8972    0.1875    1.4614    0.1800
        0.9845    0.2920    1.0329    0.0410    0.2131    1.5453    0.1836
        0.9824    0.3450    1.0372    0.0461    0.2171    1.5555    0.1835
        0.9788    0.2865    1.0334    0.0376    0.2246    1.5138    0.1835
        0.7748    0.7996    0.5837    0.6920    0.2361    1.5623    0.1562
        0.9327    0.7108    0.7932    0.8813    0.2429    1.3668    0.1714
        0.7845    0.7085    0.5863    0.8664    0.2331    1.3526    0.1599
        0.9150    0.9588    0.4602    0.0459    0.1974    1.5060    0.1100
        0.4677    0.9718    0.2208    0.8790    0.3353    1.4445    0.0304
        0.7662    0.9844    0.3257    0.9009    0.3339    1.5775    0.0569
        0.9861    0.3353    1.0541    0.8882    0.2774    1.3352    0.1836
        0.9775    0.4139    1.0475    0.8460    0.2756    1.2971    0.1832
        0.9777    0.3711    1.0232    0.8993    0.2849    1.3093    0.1833
        0.9733    0.2648    1.0185    0.6324    0.2526    1.5459    0.1836
        0.9810    0.5796    1.0048    0.5933    0.2465    1.5819    0.1816
        0.9815    0.2844    1.0322    0.5545    0.2441    1.5901    0.1836
        0.9835    0.5427    1.0007    0.0698    0.2068    1.6224    0.1819
        0.9826    0.3694    1.0204    0.0429    0.1983    1.5428    0.1834
        0.9843    0.3014    1.0279    0.0484    0.1885    1.4816    0.1836
        0.9584    0.8486    0.7909    0.7460    0.2116    1.5099    0.1617
        0.9551    0.7102    0.8435    0.8090    0.2993    1.6265    0.1733
        0.9756    0.6506    0.9829    0.8954    0.2563    1.3171    0.1779];
   
    T(:,:,1) =[
        0.8066    0.4494    0.5754    0.5389    0.2123    0.7704    0.1783
        0.6865    0.8239    0.4500    0.7125    0.2445    0.6761    0.1502
        0.9121    0.5228    0.8421    0.3831    0.2161    0.7838    0.1814
        0.9748    0.6792    0.9644    0.3480    0.1906    0.7853    0.1782
        0.9468    0.2258    0.9286    0.8554    0.5876    0.7141    0.1823
        0.8455    0.6187    0.6679    0.7109    0.3994    0.8530    0.1759
        0.8859    0.5941    0.7850    0.0783    0.2196    0.8334    0.1793
        0.8856    0.5888    0.7657    0.0598    0.2215    0.8310    0.1787
        0.8584    0.4877    0.7069    0.0552    0.2275    0.8275    0.1781
        0.9734    0.9005    0.8126    0.4977    0.2263    0.8207    0.1515
        0.9793    0.8591    0.8793    0.4733    0.2291    0.7876    0.1625
        0.9558    0.9316    0.7202    0.8790    0.1853    0.7681    0.1383
        0.9882    0.5861    1.0127    0.7975    0.1969    0.7318    0.1815
        0.9845    0.4178    1.0152    0.6533    0.1883    0.7558    0.1832
        0.9890    0.6455    0.9951    0.9009    0.1884    0.7344    0.1801
        0.9837    0.2536    1.0356    0.0413    0.2144    0.7776    0.1836
        0.9830    0.3632    1.0362    0.0465    0.2185    0.7822    0.1835
        0.9905    0.3106    1.0515    0.0377    0.2252    0.7585    0.1836
        0.7706    0.7994    0.5834    0.6958    0.2383    0.7875    0.1560
        0.9310    0.7029    0.7970    0.8855    0.2437    0.6869    0.1719
        0.7797    0.7123    0.5853    0.8754    0.2358    0.6836    0.1599
        0.8970    0.9628    0.4435    0.0481    0.2072    0.7899    0.1057
        0.4740    0.9727    0.2372    0.8739    0.3324    0.7172    0.0282
        0.7814    0.9834    0.3338    0.8947    0.3313    0.7833    0.0607
        0.9801    0.3318    1.0456    0.9003    0.2808    0.6763    0.1834
        0.9823    0.4491    1.0514    0.8502    0.2780    0.6542    0.1829
        0.9753    0.3279    1.0173    0.8974    0.2839    0.6528    0.1834
        0.9872    0.2482    1.0412    0.6360    0.2529    0.7732    0.1837
        0.9852    0.5879    1.0132    0.5975    0.2476    0.7941    0.1814
        0.9861    0.3150    1.0355    0.5516    0.2443    0.7959    0.1835
        0.9823    0.5520    0.9952    0.0698    0.2080    0.8155    0.1817
        0.9869    0.3795    1.0254    0.0431    0.1991    0.7746    0.1834
        0.9832    0.3028    1.0216    0.0487    0.1899    0.7460    0.1835
        0.9590    0.8555    0.7833    0.7563    0.2135    0.7619    0.1606
        0.9589    0.7104    0.8452    0.8130    0.3009    0.8171    0.1734
        0.9719    0.6600    0.9754    0.9133    0.2617    0.6722    0.1775];
        
    # comparing each object to the possible solution and adding 1 to the total amount
    # of particular objects if x object was found
    for k=1:numel(s)
        s(k).prop=[s(k).Solidity  s(k).Eccentricity   s(k).FormFactor...
            K(k,1) K(k,2) K(k,3) abs(log(s(k).inv(1,1)))/10];
        for i=1:36
            diff(i,1,k)=norm((T(i,:,1)-s(k).prop),'fro');% use to norm and put invmatrix to T
            diff(i,2,k)=norm((T(i,:,2)-s(k).prop),'fro');
            diff(i,3,k)=norm((T(i,:,3)-s(k).prop),'fro');
            diff(i,4,k)=norm((T(i,:,4)-s(k).prop),'fro');
            diff(i,5,k)=norm((T(i,:,5)-s(k).prop),'fro');
            solution1(i,k)=min(diff(i,:,k));
            
        end
        solutionk(k,1)= find(((solution1(:,k)-min(solution1(:,k)))==0)==1);
        
        if solutionk(k,1)==1 || solutionk(k,1)==2  || solutionk(k,1)==3
            
            solution2nd=find(diff(:,1,k)==min(setdiff((diff(:,1,k)),min((diff(:,1,k))))));
            if  solution2nd==9
                echinocyte=echinocyte+1;
            elseif solution2nd~=7 || solution2nd~=8 || solution2nd~=9
                sorted(:,1)=sort(diff(:,1,k));
                Sr=((diff(:,1,k)-sorted(3,1)==0)==1);
                Sr1=find(Sr==1);
                if Sr1==7 || Sr1==8 || Sr1==8
                    echinocyte=echinocyte+1;
                else
                    acanthocyte=acanthocyte+1;
                end
            end
            
        elseif solutionk(k,1)==4 || solutionk(k,1)==5 ||solutionk(k,1)==6
            degmacyte=degmacyte+1;
        elseif solutionk(k,1)==7 || solutionk(k,1)==8 ||solutionk(k,1)==9
            echinocyte=echinocyte+1;
        elseif solutionk(k,1)==10 || solutionk(k,1)==11 ||solutionk(k,1)==12
            eliptocyte=eliptocyte+1;
        elseif solutionk(k,1)==13 || solutionk(k,1)==14 ||solutionk(k,1)==15
            mega=mega+1;
        elseif solutionk(k,1)==16 || solutionk(k,1)==17 ||solutionk(k,1)==18
            normo=normo+1;
        elseif solutionk(k,1)==19 || solutionk(k,1)==20 ||solutionk(k,1)==21
            schistocyte=schistocyte+1;
        elseif solutionk(k,1)==22 || solutionk(k,1)==23 ||solutionk(k,1)==24
            sickle=sickle+1;
        elseif solutionk(k,1)==25 || solutionk(k,1)==26 ||solutionk(k,1)==27
            sphere=sphere+1;
        elseif solutionk(k,1)==28 || solutionk(k,1)==29 ||solutionk(k,1)==30
            stomato=stomato+1;
        elseif solutionk(k,1)==31 || solutionk(k,1)==32 ||solutionk(k,1)==33
            target=target+1;
        elseif solutionk(k,1)==34 || solutionk(k,1)==35 ||solutionk(k,1)==36
            teardrop=teardrop+1;
        end
        
        
        
    end
    
    ObjOK=acanthocyte+degmacyte+echinocyte+eliptocyte+mega+normo+schistocyte+sickle+sphere+stomato+target+teardrop;
    Dif=1;
    %d
    %  %   fprintf(txt,'%d,%s,%s,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%s \n',NumMec,NumSec,NumImg,TotNM,TotCB,TotQR,R0,R90,R180,R270,ReflCB,BadCB,TotDigCB,CBL,CBR,CBG,StringCB);
    %
    fprintf(txt,'%d,%d,%s,%s,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d \n',...
        NumMec,Dif,NumSec,NumImg,ObjBord,ObjOK,acanthocyte,degmacyte,echinocyte,...
        eliptocyte,mega,normo,schistocyte,sickle,sphere,stomato,target,teardrop);
    
    
    %figure();imshow(E3);% title('regions');
    %distancias mahal
    %get the color of gray and black of the image
end

end %of main function

function phi =  invmoments(F)
%INVMOMENTS Compute invariant moments of image.
%   PHI = INVMOMENTS(F) computes the moment invariants of the image
%   F. PHI is a seven-element row vector containing the moment
%   invariants as defined in equations (11.3-17) through (11.3-23) of
%   Gonzalez and Woods, Digital Image Processing, 2nd Ed.
%
%   F must be a 2-D, real, nonsparse, numeric or logical matrix.

%   Copyright 2002-2004 R. C. Gonzalez, R. E. Woods, & S. L. Eddins
%   Digital Image Processing Using MATLAB, Prentice-Hall, 2004
%   $Revision: 1.5 $  $Date: 2003/11/21 14:39:19 $

if (ndims(F) ~= 2) | issparse(F) | ~isreal(F) | ~(isnumeric(F) | ...
        islogical(F))
    error(['F must be a 2-D, real, nonsparse, numeric or logical ' ...
        'matrix.']);
end

F = double(F);

phi = compute_phi(compute_eta(compute_m(F)));
end

%-------------------------------------------------------------------%
function m = compute_m(F)

[M, N] = size(F);
[x, y] = meshgrid(1:N, 1:M);

% Turn x, y, and F into column vectors to make the summations a bit
% easier to compute in the following.
x = x(:);
y = y(:);
F = F(:);

% DIP equation (11.3-12)
m.m00 = sum(F);
% Protect against divide-by-zero warnings.
if (m.m00 == 0)
    m.m00 = eps;
end
% The other central moments:
m.m10 = sum(x .* F);
m.m01 = sum(y .* F);
m.m11 = sum(x .* y .* F);
m.m20 = sum(x.^2 .* F);
m.m02 = sum(y.^2 .* F);
m.m30 = sum(x.^3 .* F);
m.m03 = sum(y.^3 .* F);
m.m12 = sum(x .* y.^2 .* F);
m.m21 = sum(x.^2 .* y .* F);
end
%-------------------------------------------------------------------%
function e = compute_eta(m)

% DIP equations (11.3-14) through (11.3-16).

xbar = m.m10 / m.m00;
ybar = m.m01 / m.m00;

e.eta11 = (m.m11 - ybar*m.m10) / m.m00^2;
e.eta20 = (m.m20 - xbar*m.m10) / m.m00^2;
e.eta02 = (m.m02 - ybar*m.m01) / m.m00^2;
e.eta30 = (m.m30 - 3 * xbar * m.m20 + 2 * xbar^2 * m.m10) / m.m00^2.5;
e.eta03 = (m.m03 - 3 * ybar * m.m02 + 2 * ybar^2 * m.m01) / m.m00^2.5;
e.eta21 = (m.m21 - 2 * xbar * m.m11 - ybar * m.m20 + ...
    2 * xbar^2 * m.m01) / m.m00^2.5;
e.eta12 = (m.m12 - 2 * ybar * m.m11 - xbar * m.m02 + ...
    2 * ybar^2 * m.m10) / m.m00^2.5;
end
%-------------------------------------------------------------------%
function phi = compute_phi(e)

% DIP equations (11.3-17) through (11.3-23).

phi(1) = e.eta20 + e.eta02;
phi(2) = (e.eta20 - e.eta02)^2 + 4*e.eta11^2;
phi(3) = (e.eta30 - 3*e.eta12)^2 + (3*e.eta21 - e.eta03)^2;
% phi(4) = (e.eta30 + e.eta12)^2 + (e.eta21 + e.eta03)^2;
% phi(5) = (e.eta30 - 3*e.eta12) * (e.eta30 + e.eta12) * ...
%          ( (e.eta30 + e.eta12)^2 - 3*(e.eta21 + e.eta03)^2 ) + ...
%          (3*e.eta21 - e.eta03) * (e.eta21 + e.eta03) * ...
%          ( 3*(e.eta30 + e.eta12)^2 - (e.eta21 + e.eta03)^2 );
% phi(6) = (e.eta20 - e.eta02) * ( (e.eta30 + e.eta12)^2 - ...
%                                  (e.eta21 + e.eta03)^2 ) + ...
%          4 * e.eta11 * (e.eta30 + e.eta12) * (e.eta21 + e.eta03);
% phi(7) = (3*e.eta21 - e.eta03) * (e.eta30 + e.eta12) * ...
%          ( (e.eta30 + e.eta12)^2 - 3*(e.eta21 + e.eta03)^2 ) + ...
%          (3*e.eta12 - e.eta30) * (e.eta21 + e.eta03) * ...
%          ( 3*(e.eta30 + e.eta12)^2 - (e.eta21 + e.eta03)^2 );

end
