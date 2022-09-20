function [sumtable, dist_emonet] = imagefeatures_table(directory, ext)
% extracts luminance, entropy, and spatial frequency of visual stimuli and
% saves them as a table, also computes the distance matrix between the
% stimuli based on those features
% directory = directory in which to look for images
% ext = string that specifies the file extension, e.g., 'bmp', 'jpeg'
% requires Image Processing Toolbox (for entropy), SpatialFrequency function, Emonet toolbox
% Spatial Frequ works only with RGB, so this should maybe be optional
% by Leonie Koban, 2022
 

%% Input files
fs = dir([directory, '/*.', ext]);

for fi = 1:numel(fs)
    f(fi) = cellstr(fs(fi).name);
end

f = f';
%% Get and load the EmoNet model

model_filepath=which('netTransfer_20cat.mat');
if isempty(model_filepath)
    fprintf('Please download EmoNet from https://sites.google.com/colorado.edu/emonet \n');
end
load(model_filepath);
% display the network layers
netTransfer.Layers

[probs_all, fc8_all] = deal(NaN * zeros(length(f), 20));
Category = cell(length(f), 1);
[CategoryNum, MaxProb] = deal(zeros(length(f), 1));

%% 

for img = 1:numel(f)
    
    I = imread(f{img});
    if size(I,3)==3
        Luminance(img,1) = calcluminance(I);
    end
    
    entro(img,1) = entropy(I);
    SpatFrequ(img,1) = SpatialFrequency(f{img});
    
    I = readAndPreprocessImage(I);
    probs = netTransfer.predict(I);
    output_table = table(netTransfer.Layers(23).Classes, probs','VariableNames',{'EmotionCategory','Probability'});
    probs_all(img, :) = double(probs);
    
    [~, wh] = max(probs);
    Category{img} = output_table.EmotionCategory(wh);
    CategoryNum(img) = wh;
    MaxProb(img) = probs(wh);

    fc8_activation = netTransfer.activations(I,'fc');
    fc8_activation = squeeze(double(fc8_activation));
    
    fc8_all(img, :) = fc8_activation';

    close all
    
end

%%
imgnum = (1:img)';

sumtable = table(f, imgnum, Luminance, entro, SpatFrequ, Category, CategoryNum, MaxProb, fc8_all, 'VariableNames', {'Image', 'ImageNum', 'Luminance', 'Entropy', 'SpatFrequ', 'Category', 'CatNumber', 'Prob', 'Layer8_Activ'});

disp(sumtable)

dist_emonet = (1-corr(fc8_all'))./2;
figure; imagesc(dist_emonet)


%%
function L = calcluminance(II)
L = .2126 * mean(mean(II(:,:,1))) + 0.7152 * mean(mean(II(:,:,2))) + 0.0722 *mean(mean(II(:,:,3)));

