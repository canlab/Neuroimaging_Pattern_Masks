painpathways = load_atlas('painpathways');

% select right-hemisphere pathways for illustration?
pbn = select_atlas_subset(painpathways, {'pbn'});

amy = select_atlas_subset(painpathways, {'Amy'});

% Load amygdala subregions - CEa
painpathways_fine = load_atlas('painpathways_finegrained');
cea = select_atlas_subset(painpathways_fine, {'Amygdala_CM'});

vplm = select_atlas_subset(painpathways, {'Thal_VPLM_R'});

%dpins = select_atlas_subset(painpathways, {'dpIns_R'}); 
% TW: I don't love this - kind of messy for rendering... 

dpins = select_atlas_subset(painpathways_fine, {'Ctx_Ig_R'});

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
lightFollowView; lightRestoreSingle;

% Create a region object from the selected atlas regions
% -------------------------------------------------------------------------

r = [atlas2region(pbn) atlas2region(amy)];

% Render regions and connecting lines on the surface
% -------------------------------------------------------------------------

classes = ones(length(r), 1);   % which color/size group each region belongs to
colors = {[1 0 0]};             %  colors is cell array of colors for each class (text or rgb vector)
radiusvals = [2 2 4 4];         % radius of sphere for each region

% sigmat is an important input. It specifies which regions to connect to
% which other regions with lines. it should be a n_regions x n_regions
% matrix. It is signed, so that values of 1 are plotted in one color, and
% -1 in another color, to depict positive and negative associations.

sigmat = zeros(length(r));      % matrix of which regions to connect with lines
sigmat(1, 3) = 1;               % sigmat has signed values for significant relationships among clusters
sigmat(2, 4) = 1; 
sigmat2 = [];                   % sigmat2 is just like sigmat, but uses lighter line colors.
                                % this is another way of customizing the plot.
                                
[mov, linehandles, linehandles2] = cluster_nmdsfig_glassbrain(r,classes,colors,sigmat,sigmat2, 'existingfig', 'nobrain', 'radius', radiusvals);
 
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
 
% A variant:
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
 
drawnow, snapnow;

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


% Render regions on surface
han_pbn = isosurface(pbn);
han_amy = isosurface(amy);
set(han_amy, 'FaceAlpha', .3, 'FaceColor', [.8 .2 .1]);
han_cea = isosurface(cea);

% Create region object and matrix of connections
r = [atlas2region(pbn) atlas2region(cea)];

sigmat = zeros(length(r));      % matrix of which regions to connect with lines
sigmat(1, 3) = 1;               % sigmat has signed values for significant relationships among clusters
sigmat(2, 4) = 1; 

% Create a 'dummy' region for spinal input, connected to PBN
[hpatch, cl] = cluster_image_sphere([0 -42 -70], 'color', [1 0 0], 'radius', 3);
delete(hpatch);                         % We don't really need the sphere, so delete it.
spinal_r = cluster2region(cl);
r(end+1) = spinal_r;  
sigmat(:, end + 1) = [1 1 0 0]';
sigmat(end + 1, :) = [1 1 0 0 0];

cea = select_atlas_subset(painpathways_fine, {'Amygdala_CM'});

% Draw streamlines connecting regions

color = [.9 .4 .1];         % cannot exceed .9 for any value, because of random addition below.
thickness = .1;
bendpercent = .1;
nstreamlines = 30;

for rindx = 1:length(r)
    for rindx2 = rindx + 1:length(r)
        
        if sigmat(rindx, rindx2)
            
            for i = 1:nstreamlines % number of streamlines
                randval = unifrnd(-.1, .1);
                out = nmdsfig_tools('connect3d',r(rindx).mm_center + 20*randval,r(rindx2).mm_center + 20*randval, color + randval, thickness, bendpercent+randval);
            end
            
        end
        
    end
    
end

camzoom(1.5)
campan(0, -.5)
camzoom(1.5)
campan(0, -.5)
drawnow, snapnow;

%saveas(gcf, 'painpathways_surface_rendering1_pbn_amy.png');

%% Add spinothalamic pathway

r_st = [spinal_r atlas2region(vplm) atlas2region(dpins) ];
 
sigmat = zeros(length(r_st));      % matrix of which regions to connect with lines
sigmat(1, 2) = 1;               % sigmat has signed values for significant relationships among clusters
sigmat(2, 3) = 1; 

color = [.15 .4 .85];         % cannot exceed .9 for any value, because of random addition below.
thickness = .1;
bendpercent = .1;
nstreamlines = 30;

for rindx = 1:length(r_st)
    for rindx2 = rindx + 1:length(r_st)
        
        if sigmat(rindx, rindx2)
            
            for i = 1:nstreamlines % number of streamlines
                randval = unifrnd(-.15, .15);
                out = nmdsfig_tools('connect3d',r_st(rindx).mm_center + 20*randval, r_st(rindx2).mm_center + 20*randval, color + randval, thickness, bendpercent+randval);
            end
            
        end
        
    end
    
end
 
vplm_han = isosurface(vplm);
set(vplm_han, 'FaceColor', [.1 .3 1]);

dpins_han = isosurface(dpins);
set(dpins_han, 'FaceColor', [.1 .3 1]);

camzoom(.5);
  view(182, 12);
  lightRestoreSingle;
  
  saveas(gcf, 'painpathways_surface_rendering1_pbn_amy_st.png');

%% Alternate: Cutaways, ST pathway alone

pp = addbrain('cutaway', 'ycut_mm', -30);

set(pp(1), 'Visible', 'off');
set(pp([2 3 8:12]), 'FaceAlpha', 0);
set(pp([4 7]), 'FaceAlpha', .4);
set(pp([7]), 'FaceAlpha', .7);

color = [.15 .4 .85];         % cannot exceed .9 for any value, because of random addition below.
thickness = .1;
bendpercent = .1;
nstreamlines = 30;

for rindx = 1:length(r_st)
    for rindx2 = rindx + 1:length(r_st)
        
        if sigmat(rindx, rindx2)
            
            for i = 1:nstreamlines % number of streamlines
                randval = unifrnd(-.15, .15);
                out = nmdsfig_tools('connect3d',r_st(rindx).mm_center + 20*randval, r_st(rindx2).mm_center + 20*randval, color + randval, thickness, bendpercent+randval);
            end
            
        end
        
    end
    
end
 
vplm_han = isosurface(vplm);
set(vplm_han, 'FaceColor', [.1 .3 1]);

dpins_han = isosurface(dpins);
set(dpins_han, 'FaceColor', [.1 .3 1]);

camzoom(1.5);
view(169, -3);
lightRestoreSingle;

saveas(gcf, 'painpathways_surface_rendering1_st.png');
