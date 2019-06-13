

%% CREATE A NICE FIGURE OF ALL ROIs.
I=imread('/Users/jannisborn/Desktop/HIFO/Allen_CCF_outline.png');
BW = rgb2gray(I);

level=graythresh(BW);
BW = imbinarize( BW , 0.02 );
BW=~BW;
img = 255 * repmat(uint8(BW), 1, 1, 3);
figure
hold on;
set(gcf,'color','w');

bonds = [2, 4, 8, 10, 15];
colors = {[237, 144, 177], [144, 187, 237], [150, 237, 144], ...
    [200, 200, 200], [253, 255, 120]};
names = {rois(:).name};

bond_ind = 1;
for k = 1:length(rois)
    if k>bonds(bond_ind)
        bond_ind=bond_ind+1;
    end
    mask_w = imresize(rois(k).mask, size(BW));
    c1 = ones(size(img,1), size(img,2))*colors{bond_ind}(1);
    c2 = ones(size(img,1), size(img,2))*colors{bond_ind}(2);
    c3 = ones(size(img,1), size(img,2))*colors{bond_ind}(3);
    fill = cat(3,c1,c2,c3).*mask_w;

    for x = 1:size(img,1)
        for y = 1:size(img,2)
            if ~isequal(squeeze(fill(x,y,:)), [0;0;0])
                img(x,y,:) = fill(x,y,:);
            end
        end
    end

    
   
    
    img = im2uint8(img);
    for j = 1:k
        pos = rois(j).center .* size(BW)./size(rois(j).mask);
        if j==8 || j==10
            pos(2) = pos(2)-10;
        end
        img = insertText(img, [pos(2), pos(1)], names(j), 'BoxColor', ...
            'white', 'BoxOpacity',0, 'TextColor', 'black', 'FontSize', 30, ...
            'Font', 'Trebuchet MS Bold', 'AnchorPoint', 'Center');
    end 
    
    imshow(img, 'InitialMagnification',200);
    waitforbuttonpress
    
    
end



