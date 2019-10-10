painpathways = load_atlas('painpathways');

pbn = select_atlas_subset(painpathways, {'pbn'});

amy = select_atlas_subset(painpathways, {'Amy'});

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

% han_pbn = isosurface(pbn);
% han_amy = isosurface(amy);

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

han_pbn = isosurface(pbn);
han_amy = isosurface(amy);

color = [.1 .1 .1];
thickness = 2;
bendpercent = .1;

for rindx = 1:length(r)
    for rindx2 = rindx + 1:length(r)
        
        if sigmat(rindx, rindx2)
            
            x = [; r(rindx2).mm_center];
            
            out = nmdsfig_tools('connect3d',r(rindx).mm_center,r(rindx2).mm_center, color, thickness, bendpercent);
            
        end
        
    end
    
end

drawnow, snapnow;

 
 

  