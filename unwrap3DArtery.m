function new_proj = unwrap3DArtery( input_arg, method , ExternalInternal)
% ---------------------------------------------------------------------------- %
% Inputs:
%   - input_arg : either a 3D image of an artery or the path to this image
%   - method    : 1 or 2 (the method for the separation between the two
%                 laminae).
%   - ExternalInternal: choice of which laminae to extract, Default is 2 for
%                 the internal, with 1 extract external
%
% Outputs:
%   - proj      : the projection of the artery on a cylinder
% ---------------------------------------------------------------------------- %

    image_thickness = 2.05; % in pixels

    if isa(input_arg, 'char')
        dir0 = input_arg;
        dir1 = dir(strcat(input_arg,'/*.tif'));
        numSubDir = size(dir1,1);

        firstSample = imread(strcat(dir0,'/',(dir1(1).name)));
        [rows,cols,levs] = size(firstSample);

        disp(strcat(num2str(numSubDir),...
                    ' images in folder, with dimensions: ', num2str(rows),...
                    ' x  ', num2str(cols),...
                    ' x  ', num2str(levs)));
        
        artery(rows,cols,numSubDir) = 0;
        % %counterDir=31;
        for counterDir=1:numSubDir
            disp(counterDir) 
            tempDir = strcat(dir0,'/',dir1(counterDir).name);
            dataIn = imread((tempDir));
            artery(:,:,counterDir) = double(dataIn);
        end
    else
        artery = input_arg;
    end

    [length, width, depth] = size(artery);

    % default values in case these were not passed as arguments
    if ~exist('method','var'); method =2; end
    if ~exist('ExternalInternal','var'); ExternalInternal =2; end
    
    
    
% ============================================================================ %
%                         ROTATION OF THE ARTERY
% ============================================================================ %
    disp(' 1 / 8 - ROTATION OF THE ARTERY');

    step = 10;
    
    % --------------- FIRST ROTATION ---------------
    
    % projection on the x,y plan
    proj = mean(artery,3);
    
    [X,Y,W_diag] = deal([]); 
    for i = 1:step:(length - step)
        for j = 1:step:(width - step)
            X(end+1,:) = [1,i];
            Y(end+1,:) = j;
            W_diag(end+1,:) = proj(i,j); 
        end
    end
    W = spdiags(W_diag,0,size(X,1),size(X,1));
    
    clear W_diag;
    
    % weighted linear regression
    Theta = (X'*W*X)\(X'*W*Y);
    clear X Y W;
    
    % Y = X'*Theta <=> (i,j) = u + v*t with u and v two vectors ant t in R
    % u = K * (0, Theta_1)
    % v = K * (1, Theta_2)
    % with K an constant
    v = [1, Theta(2)]';
    v = v / sqrt(v'*v);
    
    % the angle around the z axis
    alpha = - atan(([0,1] * v) / ([1,0] * v));
    
    % rotation
    for k = 1:depth
        artery(:,:,k) = imrotate(artery(:,:,k),...
                                 alpha/pi*180, 'bicubic', 'crop');
    end
    
    % --------------- SECOND ROTATION ---------------
    
    % projection on the x,k plan
    proj = reshape(mean(artery,2), length, depth);
    
    [X,Y,W_diag] = deal([]); 
    for i = 1:step:(length - step)
        for k = 1:step:(depth - step)
            X(end+1,:) = [1,i];
            Y(end+1,:) = k;
            W_diag(end+1,:) = proj(i,k)^4; 
        end
    end
    W = spdiags(W_diag,0,size(X,1),size(X,1));
    
    clear W_diag;
    
    % weighted linear regression
    Theta = (X'*W*X)\(X'*W*Y);
    
    clear X Y W;
    
    % Y = X'*Theta <=> (i,k) = u + v*t with u and v two vectors ant t in R
    % u = K * (0, Theta_1)
    % v = K * (1, Theta_2)
    % with K an constant
    v = [1, Theta(2)]';
    v = v / sqrt(v'*v);
    
    % the angle around the y axis
    alpha = -atan(([0,1] * v) / ([1,0] * v));
    
    % rotation
    for j = 1:width
        artery(:,j,:) = imrotate(reshape(artery(:,j,:), length, depth),... 
                                 alpha/pi*180, 'bicubic', 'crop');
    end
    
% ============================================================================ %
%                        PROJECTION ON THE Y,Z PLAN
% ============================================================================ %
    disp(' 2 / 8 - PROJECTION ON THE Y,Z PLAN');
    
    proj = reshape(sum(artery, 1), width, depth)./...
           reshape(sum((artery>0)*1, 1), width, depth);
    
% ============================================================================ %
%                            POLYNOMIAL REGRESSION
% ============================================================================ %
    disp(' 3 / 8 - POLYNOMIAL REGRESSION');
    
    % finding the limits of the points which will be used for the regression
    j_min = find(proj(1:end/2, depth) > 0.1*max(proj(1:end/2, depth)),...
                 1, 'last');
    j_max = find(proj(end/2:end, depth) > 0.1*max(proj(end/2:end, depth)), ...
                 1, 'first') + width/2;
             
    % creation of the variables for the weighted linear regression
    [X,Y,W_diag] = deal([]);
    for j = j_min:step:j_max
        for k = 1:depth
            x = zeros(1,3);
            for d = 0:2
                x(d+1) = j^(2 - d);
            end
            X(end+1,:) = x;
            Y(end+1,:) = k*image_thickness;
            W_diag(end+1,:) = proj(j,k);
        end
    end
    
    %clear j_min j_max proj;
    
    W = spdiags(W_diag,0,size(W_diag,1),size(W_diag,1));
    
    clear W_diag;
    
    % weighted polynomial regression
    P = (X'*W*X)\(X'*W*Y);
    
    clear X Y W;
    
% ============================================================================ %
%                   FINDING THE CYLINDER FOR THE PROJECTION
% ============================================================================ %
    disp(' 4 / 8 - FINDING THE CYLINDER FOR THE PROJECTION');
    
    % P(x) = a*(x-a)^2 + c, the curvature radius = 1 / f''(x)
    % Then the calculation are easy
    
    % coordinates of the circle's center
    C = [0, ... 
         - P(2) / (2 * P(1)), ...
         1 / (2 * P(1)) + P(3) - P(2)^2 / (4 * P(1))];
    
    % circle's radius
    R = 1 / (2 * P(1));
    
    clear P;

% ============================================================================ %
%                            SEPARATOR BUILDING
% ============================================================================ %
    disp(' 5 / 8 - SEPARATOR BUILDING');
    
    switch method
        
        case 1
            % computiong the projection on the cylinder
            proj_2 = zeros(2*width, floor(2*R)); 

            for j = 1:width
                for k = 1:depth
                    theta = -atan( (C(2) - j) / (C(3) - k*image_thickness) );
                    p_j = round(theta * R + (C(2)-length/2) + length );

                    p_k = round(1+sqrt((C(2)-j)^2+(C(3)-k*image_thickness)^2));

                    p_j = min(2*width, p_j);

                    proj_2(p_j, p_k) = proj(j,k);  
                end
            end

            hist = mean(proj_2,1);
            hist = hist ./ max(hist);
            hist = conv(hist,1/10*ones(10,1), 'same');
            
            [pks,locs,w,p] = findpeaks(hist);
            
            switch ExternalInternal
                case 1
                    for i = size(pks,2):-1:1
                        if pks(i) > 0.1
                            mu = locs(i);
                            break;
                        end
                    end
                    
                case 2
                    for i = 1:size(pks,2)
                        if pks(i) > 0.1
                            mu = locs(i);
                            break;
                        end
                    end
            end

        case 2
            % Loop over a variation of canny scales, then select the first scale 
            % at which there are only 2 edges (presumably both will be close to 
            % the external and internal lamina on the outside of them.
            % Under this condition, the first edge to be found is the
            % External (1) and the second is the Internal (2)
            
            % In some special cases where the artery is close to the edge,
            % it may be that instead of 2 edges there will be 3 as the
            % external is broken into 2 sections.
            numEdges = zeros(21,1);
            for k=1:21
                edgesArtery = edge(proj, 'canny', [], k);
                [~, numEdges(k)] = bwlabel(edgesArtery);
            end
            
            if numEdges(end)==2
                scaleEdges = find(numEdges == 2, 1, 'first');               
                % calculate the edges again, and Label to detect the smallest
                % (internal) and the largest (external).
                try
                    edgesArtery = edge(proj, 'canny', [], scaleEdges);
                catch
                    q=1;
                end
                
                [edgesL] = bwlabel(edgesArtery);
                
                % the boundary between internal and external can be calculated by a
                % dilation of the internal lamina.
                validRegion = imdilate(edgesL==ExternalInternal,ones(55,35));
            elseif numEdges(end)==3
                 scaleEdges = find(numEdges == 3, 1, 'first');               
                % calculate the edges again, and Label to detect the smallest
                % (internal) and the largest (external).
                try
                    edgesArtery = edge(proj, 'canny', [], scaleEdges);
                 catch
                     q=1;
                 end
                 
                [edgesL] = bwlabel(edgesArtery);
                
                % the boundary between internal and external can be calculated by a
                % dilation of the internal lamina.
                if ExternalInternal ==2
                    validRegion = imdilate(edgesL==3,ones(55,35));
                else
                    validRegion = imdilate((edgesL==1)|(edgesL==2),ones(55,35));
                    validRegion = imclose(validRegion,ones(120,1));
                end
            end
            
    end

% ============================================================================ %
%                               PROJECTION
% ============================================================================ %
    disp(' 6 / 8 - PROJECTION');
    
    % computiong the projection on the cylinder
    proj = zeros(length, 2*width); 
    for j = 1:width
        for k = 1:depth
            theta = -atan( (C(2) - j) / (C(3) - k*image_thickness) );
            p_j = round( theta * R + (C(2)-width/2) + width );
            p_k = round(sqrt((C(2) - j)^2 + (C(3) - k*image_thickness)^2));
            for i = 1:length 
                z = i;
                p_i = z ;
                
                switch method
                    case 1
                        if abs(p_k - mu) < 15
                            proj(p_i, p_j) = proj(p_i, p_j) + artery(i,j,k); 
                        end
                    case 2
                        if validRegion(j,k) == 1
                            proj(p_i, p_j) = proj(p_i, p_j) + artery(i,j,k); 
                        end 
                end
            end
        end
    end

% ============================================================================ %
%                      SELECTION OF THE CENTRAL REGION
% ============================================================================ %
    disp(' 7 / 8 - SELECTION OF THE CENTRAL REGION');
    
    hist = sum(proj, 1)>0;
    center = round(sum((1:size(hist,2)) .* hist)/ sum(hist)); 
    
    vert_margin = 0; % the margins in pixels for the final selection 
    hori_margin = 0;

    proj = proj(1+vert_margin:...
                length - vert_margin,...
                1+center - width/2 + hori_margin:...
                center + width/2 - hori_margin);
              
% ============================================================================ %
%                         REMOVAL OF THE ARTEFACTS
% ============================================================================ %
    disp(' 8 / 8 - REMOVAL OF THE ARTEFACTS');

    new_proj = proj;
    smooth_proj = imgaussfilt(proj, [1,5]);
    
    range = 10;
    for j = 2:size(proj,2)
        for i = 1:size(proj,1)
            a = max(1, i-range);
            b = min(size(proj,1), i+range); 
            u = smooth_proj(a:b,j);
            v = proj(a:b,j);
            if v'*v ~= 0
                lambda = sqrt((u'*u) / (v'*v));
                new_proj(i,j) = proj(i,j)*lambda;
            end
        end
    end
    
    % normalisation
    mu = mean(new_proj(:));
    new_proj = new_proj-mu;
    sigma = std(new_proj(:));
    new_proj = new_proj/sigma;
    % Low pass filtering
    new_proj = imfilter (new_proj, [ 0.0625 0.1250 0.0625;  0.1250 0.2500 0.1250; 0.0625 0.1250 0.0625]);
end



