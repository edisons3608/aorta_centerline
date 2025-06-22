%% Automated Aortic Centerline Tracking & Diameter Measurement
% Edison Sun
% 6/22/2025

%% Loading specific nifti file. Extracting metadata.
nii      = niftiinfo('path');
vol      = niftiread(nii);
mask     = vol > 0;                       % binary aorta mask
voxelSz  = nii.PixelDimensions;           % [dx, dy, dz] in mm
T        = nii.Transform.T;               % 4×4 affine → world

%% Develop point cloud and linear indices from matrix.
maskLinIdx = find(mask);                  % linear indices of all mask voxels
[idxX,idxY,idxZ] = ind2sub(size(mask), maskLinIdx);
voxels    = [idxX,idxY,idxZ,ones(numel(idxX),1)]';
worldPts  = (T * voxels)';                
worldPts  = worldPts(:,1:3);              % N×3 list of points in mm

worldPts = get_pointcloud(vol);

%% Get centroids of two lowest slices (axially)
uniqueZ = unique(idxZ);
z1 = uniqueZ(1); z2 = uniqueZ(2);
pts1 = worldPts(idxZ==z1,:);
pts2 = worldPts(idxZ==z2,:);
c1   = mean(pts1,1);
c2   = mean(pts2,1);

%% Params
maxSteps    = 100;
stepSize    = 2;
wobbleAng   = 40;
angleInc    = 2;
slabHalfThk = stepSize/2;

centerline = nan(maxSteps,3);
diameters  = nan(maxSteps,1);
centerline(1:2,:) = [c1; c2];

nbrOffsets = [1 0 0; -1 0 0; 0 1 0; 0 -1 0; 0 0 1; 0 0 -1];

for i = 3:maxSteps
    % propagate
    dirVec = centerline(i-1,:) - centerline(i-2,:);
    dirVec = dirVec / norm(dirVec);
    candPt = centerline(i-1,:) + dirVec*stepSize;
    
    % build local frame
    if abs(dot(dirVec,[1 0 0])) < 0.9
        arb = [1 0 0];
    else
        arb = [0 1 0];
    end
    u = cross(dirVec,arb); u = u/norm(u);
    v = cross(dirVec,u);
    
    bestArea = inf;
    bestXY   = [];
    bestCompPts = [];
    
    % wobble search
    for alpha = -wobbleAng:angleInc:wobbleAng
      Ra = axang2rotm([u,deg2rad(alpha)]);
      for beta = -wobbleAng:angleInc:wobbleAng
        Rb = axang2rotm([v,deg2rad(beta)]);
        n  = (Rb*(Ra*dirVec'))'; n = n/norm(n);
        
        % select slab voxels
        d = (worldPts - candPt)*n';
        slabIdx = abs(d) <= slabHalfThk;
        if ~any(slabIdx), continue; end
        
        % flood‐fill in 3D to isolate connected comp around candPt
        slabLin   = maskLinIdx(slabIdx);
        slabSubs  = [idxX(slabIdx), idxY(slabIdx), idxZ(slabIdx)];
        % seed: nearest slab voxel to candPt
        distToSeed = vecnorm(worldPts(slabIdx,:) - candPt,2,2);
        [~,seedLoc] = min(distToSeed);
        seedIdx     = slabLin(seedLoc);
        
        % build a quick lookup for slab membership
        slabMask3D = false(size(mask));
        slabMask3D(slabLin) = true;
        
        % BFS flood‐fill
        visited = false(size(mask));
        queue   = seedIdx;
        visited(seedIdx) = true;
        compLin = seedIdx;
        while ~isempty(queue)
          cur = queue(1); queue(1) = [];
          [x,y,z] = ind2sub(size(mask), cur);
          for k = 1:6
            nb = [x y z] + nbrOffsets(k,:);
            if any(nb<1) || nb(1)>size(mask,1) || nb(2)>size(mask,2) || nb(3)>size(mask,3)
              continue;
            end
            nbLin = sub2ind(size(mask), nb(1),nb(2),nb(3));
            if slabMask3D(nbLin) && ~visited(nbLin)
              visited(nbLin) = true;
              queue(end+1)   = nbLin;    
              compLin(end+1) = nbLin;     
            end
          end
        end
        
        % get the world‐pts of that component
        compMask = ismember(maskLinIdx, compLin);
        compPts  = worldPts(compMask,:);
        if size(compPts,1) < 3, continue; end
        
        % project to (u,v) and area
        rel = compPts - candPt;
        xy  = [rel*u', rel*v'];
        try
          k    = convhull(xy(:,1), xy(:,2));
          area = polyarea(xy(k,1), xy(k,2));
        catch
          continue;
        end
        
        if area < bestArea
          bestArea = area;
          bestXY   = xy;
          bestCompPts = compPts;
        end
      end
    end
    
    % check
    if isempty(bestXY)
      warning('Step %d: no valid slice → stopping.', i);
      break;
    end
    
    % record
    centerline(i,:) = mean(bestCompPts,1);
    diameters(i)    = max(pdist(bestXY));
    
    % f) bounds
    if any(centerline(i,:) < min(worldPts)) || any(centerline(i,:) > max(worldPts))
      warning('Step %d exited bounds → stopping.', i);
      break;
    end
end
function coords = get_pointcloud(img)
    linIdx = find(img == 1);            % linear indices of all voxels == 1
    [i, j, k] = ind2sub(size(img), linIdx);
    coords = [i, j, k];


end

%% Trim & visualize
%last = find(~isnan(diameters),1,'last');
last = 80;
ctr  = centerline(1:last,:);
dmd  = diameters(1:last);
dist = (0:last-1)'*stepSize;

figure('Position',[100 100 1000 400]);
subplot(1,2,1);
scatter3(worldPts(:,1),worldPts(:,2),worldPts(:,3),1,'.k'); hold on;
plot3(ctr(1:last,1),ctr(1:last,2),ctr(1:last,3),'r-','LineWidth',2);
%plot3(ctr(1:75,1),ctr(1:75,2),ctr(1:75,3),'r-','LineWidth',2);
%plot3(ctr(1:54,1),ctr(1:54,2),ctr(1:54,3),'r-','LineWidth',2);
axis equal; grid on;
xlabel('X'); ylabel('Y'); zlabel('Z');
title('Aortic Centerline');

subplot(1,2,2);
plot(1:last,dmd,'-o','LineWidth',1.5);
grid on;
xlabel('Steps');
ylabel('Max diameter');
title('Diameter Profile');
