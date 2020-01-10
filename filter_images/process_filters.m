%img = imread('C:\Users\heidi\Downloads\20190220_103858.jpg');imshow(img)
%hEllipse = imellipse(gca,[510 754 1680 1680]);

%img = imread('C:\work\LTER\filter_images\20190220_103549.jpg'); imshow(img)
%hEllipse = imellipse(gca,[535.8 757.7 1680 1680]);
%disp('move the cirlce as needed')
%pause 
imgpath = '\\sosiknas1\Lab_data\LTER\20190201_EN627\GFF_images\Paul_camera\';
files = dir([imgpath '*.jpg']);
load([imgpath 'SummaryTable'])
defaultPos = [2066 -4 4200 4200];

for filenum = 30:length(files)
    img = imread([imgpath files(filenum).name]);
    [~,rownum] = intersect(files(filenum).name, SummaryTable.File);
    figure(1), clf, imshow(img)
    if isempty(rownum)
        hEllipse = imellipse(gca,defaultPos);
        disp('move the circle as needed')
        pause
        newEntry.File = files(filenum).name;
        newEntry.CirclePosition = hEllipse.getPosition;
    else
        hEllipse = imellipse(gca,SummaryTable.CirclePosition(rownum));
    end
    
    %img = imread('\\sosiknas1\Lab_data\LTER\20190201_EN627\GFF_images\Paul_camera\EN627_C22_N09_a_IMG_6517.jpg'); imshow(img)
    %hEllipse = imellipse(gca,[2066 -4 4200 4200]);
    
    
    %hEllipse = imellipse(gca,[500 750 1680 1680]);
    %hEllipse = imellipse(gca,[510 754 1680 1680]);
    %hpos = hEllipse.getPosition;
    
    img2 = img;
    mask = imcomplement(hEllipse.createMask());
    
    mask3 = cat(3, mask, mask, mask);img2= img;img2(mask3) = 255; imshow(img2)
    
    %[BW, masked_img] = createMask_test2(img2);
    [BW, masked_img] = createMask4(img2);
    %[L, num] = bwlabel(BW, 4);
    %rgb = label2rgb(L); imshow(rgb)
    
    r = regionprops(BW, 'Area', 'FilledArea', 'PixelIdxList', 'MajorAxisLength', 'MinorAxisLength');
    B = bwboundaries(BW,'noholes');
    L = bwlabel(BW);
    %imshow(img2)
    %hold on
    %for k = 1:length(B)
    %   boundary = B{k};
    %   plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', .5)
    %end
    
    t = find([r.FilledArea]>600 & [r.FilledArea]<3000 );
    imshow(img2)
    hold on
    for k = 1:length(t)
        boundary = B{t(k)};
        plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', .5)
    end
    
    BW2 = BW;
    rdelete = r; rdelete(t) = [];
    for rnum = 1:length(rdelete)
        BW2(rdelete(rnum).PixelIdxList) = 0;
    end
    
    disp('Click blobs to remove')
    g = ginput;
    g = round(g);
    for rnum = 1:size(g,2);
        t = L(g(rnum,2),g(rnum,1));
        BW2(r(t).PixelIdxList) = 0;
    end
    
    
end

