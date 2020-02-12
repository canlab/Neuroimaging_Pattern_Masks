painpathways = load_atlas('painpathways');

pbn = select_atlas_subset(painpathways, {'pbn'});

amy = select_atlas_subset(painpathways, {'Amy'});

s1 = select_atlas_subset(painpathways_fine, {'Ctx_1' 'Ctx_3b'});
s2 = select_atlas_subset(painpathways_fine, {'Ctx_OP1'});

% Load amygdala subregions - CEa
painpathways_fine = load_atlas('painpathways_finegrained');
cea = select_atlas_subset(painpathways_fine, {'Amygdala_CM'});

vplm = select_atlas_subset(painpathways, {'Thal_VPLM_R'});

%dpins = select_atlas_subset(painpathways, {'dpIns_R'}); 
% TW: I don't love this - kind of messy for rendering... 

dpins = select_atlas_subset(painpathways_fine, {'Ctx_Ig_R'});

% Create a 'dummy' region for spinal input, connected to PBN
[hpatch, cl] = cluster_image_sphere([0 -42 -70], 'color', [1 0 0], 'radius', 3);
delete(hpatch);                         % We don't really need the sphere, so delete it.
spinal_r = cluster2region(cl);

% Below are two ways of rendering a pathway connecting these two atlas
% regions (2 right, 2 left-hem for each) on a surface.

% NOTE: We should create an object method for brainpathways objects that
% uses this. sigmat could use the .connections_apriori or .connectivity
% properties. 

%% Method 1: Use cluster_nmdsfig_glassbrain

% Create a brain surface on which to render pathways
% -------------------------------------------------------------------------

create_figure('brain'); 
han_surf = addbrain('hires left');
han_surf = [han_surf addbrain('brainstem')];
han_surf = [han_surf addbrain('hires right')];
set(han_surf(1), 'FaceAlpha', 1);
set(han_surf(2), 'FaceAlpha', .1, 'FaceColor', [.5 .5 .5]);
set(han_surf(end), 'FaceAlpha', .05);
view(60, 10); 

anat = fmri_data(which('keuken_2014_enhanced_for_underlay.img'), 'noverbose');
[phan, ~, isosurf, isocap] = isosurface(anat, 'thresh', 140, 'nosmooth', 'ylim', [-Inf -27], 'zlim', [0 Inf]);

lightFollowView; lightRestoreSingle;

% han_pbn = isosurface(pbn);
% han_amy = isosurface(amy);

% Create a region object from the selected atlas regions
% -------------------------------------------------------------------------

r = [atlas2region(pbn) atlas2region(amy)];

r = [r spinal_r];

% Render regions and connecting lines on the surface
% -------------------------------------------------------------------------

classes = ones(length(r), 1);   % which color/size group each region belongs to
colors = {[1 0 0]};             %  colors is cell array of colors for each class (text or rgb vector)
radiusvals = [2 2 4 4 2];         % radius of sphere for each region

% sigmat is an important input. It specifies which regions to connect to
% which other regions with lines. it should be a n_regions x n_regions
% matrix. It is signed, so that values of 1 are plotted in one color, and
% -1 in another color, to depict positive and negative associations.

sigmat = zeros(length(r));      % matrix of which regions to connect with lines
sigmat(1, 3) = 1;               % sigmat has signed values for significant relationships among clusters
sigmat(2, 4) = 1; 
sigmat(1, 5) = 1;               % spinal input
sigmat(2, 5) = 1;               % spinal input
sigmat2 = [];                   % sigmat2 is just like sigmat, but uses lighter line colors.
                                % this is another way of customizing the plot.
                                
[mov, linehandles, linehandles2] = cluster_nmdsfig_glassbrain(r,classes,colors,sigmat,sigmat2, 'existingfig', 'nobrain', 'radius', radiusvals);
 
set(linehandles, 'Color', [1 .3 0]);

drawnow, snapnow;

%   var args:
%   'samefig' or 'existingfig' : do not create new figure
%   'blobs' to image blobs
%   'spheres' to image spheres (default)
%   'nobrain' to avoid creating any new brain surface objects; just use
%   existing
%   'radius', followed by radius of sphere, or vector of radii
%   'movie', make movie
%   'straight', no bend
 
%% spinothalamic
% -------------------------------------------------------------------------

r = [spinal_r atlas2region(vplm) atlas2region(dpins) atlas2region(s1) atlas2region(s2)];
% 1 1 1 4 2

sigmat =   [0 0 0 0 0 0 0 0 0;
            0 0 0 0 0 0 0 0 0;
            0 0 0 0 0 0 0 0 0;
            0 0 0 0 0 0 0 0 0;
            0 0 0 0 0 0 0 0 0;                        
            0 0 0 0 0 0 0 0 0;
            0 0 0 0 0 0 0 0 0;
            0 0 0 0 0 0 0 0 0;            
            0 0 0 0 0 0 0 0 0];
        
sigmat =   [0 1 0 0 0 0 0 0 0;
            0 0 1 1 1 1 1 1 1;
            0 0 0 0 0 0 0 0 0;
            0 0 0 0 0 0 0 0 0;
            0 0 0 0 0 0 0 0 0;                        
            0 0 0 0 0 0 0 0 0;
            0 0 0 0 0 0 0 0 0;
            0 0 0 0 0 0 0 0 0;            
            0 0 0 0 0 0 0 0 0];  
        
            % Render regions and connecting lines on the surface
% -------------------------------------------------------------------------

classes = ones(length(r), 1);   % which color/size group each region belongs to
colors = {[.2 .2 1]};             %  colors is cell array of colors for each class (text or rgb vector)
radiusvals = [2 2 3 3 3 3 3 3 3 3];         % radius of sphere for each region

% sigmat is an important input. It specifies which regions to connect to
% which other regions with lines. it should be a n_regions x n_regions
% matrix. It is signed, so that values of 1 are plotted in one color, and
% -1 in another color, to depict positive and negative associations.

sigmat2 = [];                   % sigmat2 is just like sigmat, but uses lighter line colors.
                                % this is another way of customizing the plot.
                                
[mov, linehandles, linehandles2] = cluster_nmdsfig_glassbrain(r,classes,colors,sigmat,sigmat2, 'existingfig', 'nobrain', 'radius', radiusvals);
 
set(linehandles, 'Color', [0 .3 1]);

drawnow, snapnow;

%% A variant:
% Render regions and connecting lines using region isosurfaces
% -------------------------------------------------------------------------

create_figure('brain'); 
han_surf = addbrain('hires left');
han_surf = [han_surf addbrain('brainstem')];
han_surf = [han_surf addbrain('hires right')];
set(han_surf(1), 'FaceAlpha', 1);
set(han_surf(2), 'FaceAlpha', .1, 'FaceColor', [.5 .5 .5]);
set(han_surf(end), 'FaceAlpha', .05);
view(60, 10); 
lightFollowView; lightRestoreSingle;

[mov, linehandles, linehandles2] = cluster_nmdsfig_glassbrain(r,classes,colors,sigmat,sigmat2, 'existingfig', 'nobrain', 'blobs');
 
set(linehandles, 'Color', [1 .3 0]);

drawnow, snapnow;

%%
%mov = movie_tools('lines',startcoords,endcoords,mov,color,bendval,movlength,startt,endt,handles,targetaz,targetel);
x = r(1).mm_center;
y = r(3).mm_center;
han = [han_surf];
mov = movie_tools('lines',x,y,mov,'b',[0 -.1 0],2,0,.8,han,90,10);


%% Method 2: Rendering using object methods, and draw lines using nmdsfig_tools

% Create a brain surface on which to render pathways
% -------------------------------------------------------------------------

create_figure('brain'); 
han_surf = addbrain('hires left');
han_surf = [han_surf addbrain('brainstem')];
han_surf = [han_surf addbrain('hires right')];
set(han_surf(1), 'FaceAlpha', 1);
set(han_surf(2), 'FaceAlpha', .1, 'FaceColor', [.5 .5 .5]);
set(han_surf(end), 'FaceAlpha', .05);
view(60, 10); 
lightFollowView; lightRestoreSingle;

han_pbn = isosurface(pbn);
han_amy = isosurface(amy);

color = [.9 .3 .1];
thickness = 2;
bendpercent = .1;

for rindx = 1:length(r)
    for rindx2 = rindx + 1:length(r)
        
        if sigmat(rindx, rindx2)
            
            x = [; r(rindx2).mm_center];
            
            out = nmdsfig_tools('connect3d',r(rindx),r(rindx2), 'color', color, 'thickness', thickness, 'bendpercent', bendpercent);
            
        end
        
    end
    
end

drawnow, snapnow;

 
 

  