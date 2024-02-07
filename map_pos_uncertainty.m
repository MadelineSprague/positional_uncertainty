% Madeline Sprague (2-7-2024) - mjsprague@uri.edu
%
% This script takes inputs defining a given region of interest surrounding
% any number of sound source locations, and generates a map of positional
% uncertainty for an acoustic receiver (e.g. a RAFOS float) at every grid
% cell within the region. This work is motivated by the recognition that
% each arrival measured from a sound source has an uncertainty that is a
% sum of several different effects and that when you combine uncertainties
% between arrivals from different sources, this produces a final positional
% uncertainty that is dependent on the geometry of the receiver relative to
% the sources. 
%
% This script makes use of the 'intersect' function from the Antenna
% toolbox and 'latlon2local' from the Automated Driving toolbox. If these
% aren't installed, MATLAB will prompt an install when the functions are
% applied. 
%
% Inputs: 
%
%     lon, lat: double (Nx1)
%         Arrays of sound source positions
% 
%     mindepth: double 
%         Subsurface acoustic propagations often do not transmit well in
%         some shallow environments, thus mindepth allows you to skip
%         computing values for grid cells with a depth shallower than, for
%         example, 100 m (mindepth = 100). 
% 
%     maxrange: double
%         The max range that you want to model from each source. The bounds
%         of the map will be determined by the outermost limits of
%         <maxrange> from the outermost sources. 
% 
%     gridspace: double 
%         The distance between grid cells in kilometers. The main
%         determining factor of runtime along with maxrange. 
% 
%     varargin: lon_grid and lat_grid, double arrays
%         The user may enter lon_grid, lat_grid as final arguments and
%         override the automated grid that is constructed using lat/lon,
%         maxrange, and gridspace. Allows for easier comparisons within a
%         given static region if source positions are moved. 
%
% Outputs: 
%
%     results: struct 
%         d_pos_radius: double array
%             The plotted variable describing the positional uncertainty in
%             kilometers. 
%         lon_grid, lat_grid: double array
%             Grids of longitude/latitude determined by the function or
%             inputted via varargin. 
%         mask: logical array
%             A mask for grid cells that are shallower than the user's
%             chosen mindepth. 
%         bathy: double array
%             The bathymetry depths (positive downwards) of the region. 
        
function [results] = map_pos_uncertainty(lon, lat, mindepth, maxrange, gridspace, varargin)

tic;            % monitor runtime 
c = 1.48;       % sound speed (km/s)
dc = 0.003 * c; % sound speed uncertainty 
warning('off', 'MATLAB:polyshape:repairedBySimplify'); 

if isempty(varargin)

    % define lat/lon bounds
    
        min_lat = min(lat - maxrange./(110.574)); 
        max_lat = max(lat + maxrange./(110.574)); 
    
        min_lon = min(lon - maxrange./(110.574 .* cosd(lat))); 
        max_lon = max(lon + maxrange./(110.574 .* cosd(lat))); 
    
    % define the lat/lon grid of receivers
    
        lon_dist = (max_lon - min_lon) * 110.574 * cosd(min_lat); 
        lat_dist = (max_lat - min_lat) * 110.574; 
    
        lon_vec  = linspace(min_lon, max_lon, round(lon_dist/gridspace)); 
        lat_vec  = linspace(min_lat, max_lat, round(lat_dist/gridspace))'; 
    
        [lon_grid, lat_grid] = meshgrid(lon_vec, lat_vec); 

elseif ~isempty(varargin)

        lon_grid = varargin{1}; 
        lat_grid = varargin{2}; 

end

% determine depths at each grid cell and create a mask of grid cells that
% are below the minimum depth threshold 

    [bathy, ~, ~] = gebco_bathy_3d(lon_grid, lat_grid); 
    mask = bathy > -mindepth; 

% determine the ranges, from each source, at each grid cell. At the same
% time, we determine the arrival uncertainty, from each source, at each
% grid cell as a function of range. 

    nsource = length(lon); 
    range   = NaN([size(lon_grid) nsource]); % distance from source

    for k = 1:nsource % 1 range array for each source

        range_temp = NaN(size(lon_grid)); 

         for i = 1:numel(lon_grid)

            range_temp(i) = m_lldist([lon(k) lon_grid(i)], [lat(k) lat_grid(i)]); 

         end

         range(:,:,k) = range_temp; 

    end
    
    D       = 14*24*60*60;                                       % deployment time of both tags and sources
    d_arr   = (20/1e9 + 1/1e6)*D + (1.8)*0.3 + 0.0058 .* range;  % formulation of arrival uncertainty (described in readme)
    d_range = d_arr * (c + dc);                                  % the range uncertainty associated with each grid cell for each source 

    min_range = range - d_range; 
    max_range = range + d_range; 
    
% now we use a cartesian grid to find the overlap zone from the uncertainty
% rings surrounding each source. We only need to find the cartesian
% coordinates for the sources, so we redefine the lat/lon grid into
% cartesian points then identify the source locations within them, using
% the latlon2local function. 

    [x_grid,y_grid,~] = latlon2local(lat_grid, lon_grid, ...
        zeros(size(lon_grid)), [lat(1) lon(1) 0]);             % origin located at 1st inputted source position
    x_grid = x_grid/1000;                                      % convert to kilometers
    y_grid = y_grid/1000;                                      % convert to kilometers
    sx = interp2(lon_grid, lat_grid, x_grid, lon, lat);        % find source x values
    sy = interp2(lon_grid, lat_grid, y_grid, lon, lat);        % find source y values

% run the algorithm 

    areas    = NaN(size(lon_grid)); 
    ntheta   = 300; 
    theta    = linspace(0, 2*pi, ntheta); 
    costheta = cos(theta); 
    sintheta = sin(theta); 
    numcells = sum(double(mask == 0), 'all'); 

    for i = 1:height(lon_grid)

        if i == round(height(lon_grid)/20) % runtime indicator 

            runtime_prelim = toc; 
            numcells_completed = sum(mask(1:i,:) == 0, 'all'); 
            percent_completion = numcells_completed / numcells * 100; 
            runtime_expected = runtime_prelim * 100/percent_completion; 

            disp([char(string(percent_completion)) '% completed;']); 

            % expected run time display 

                if runtime_expected >= 60*60 
                    runtime_expected_h = runtime_expected/(60*60); % runtime in hours 
                    disp(['Expected runtime: ' char(string(runtime_expected_h)) ' hours'])
                elseif runtime_expected >= 60 & runtime_expected < 60*60
                    runtime_expected_m = runtime_expected/60;      % runtime in minutes 
                    disp(['Expected runtime: ' char(string(runtime_expected_m)) ' minutes'])
                elseif runtime_expected < 60
                    disp(['Expected runtime: ' char(string(runtime_expected)) ' seconds'])
                end

            % expected time of completion 

                etd = now + runtime_expected / 86400; 
                disp(['Expected time of completion: ' char(datestr(etd))]); 

        end
        
        for j = 1:width(lon_grid)

            if mask(i,j) == 1 % omit points shallower than minimum depth 
                areas(i,j) = NaN; 
            end
        
            for k = 1:length(lon) % generate polyshape objects for each source ring 
                
                rings(k) = polyshape(sx(k) + [min_range(i,j,k) .* costheta, max_range(i,j,k) .* costheta], ...
                                     sy(k) + [min_range(i,j,k) .* sintheta, max_range(i,j,k) .* sintheta]); 

            end
        
            intersect_shape = intersect(rings);      % shape of intersection
            areas(i,j)      = area(intersect_shape); % area of intersection (km^2)
        
        end
        
    end
    
% convert area into equivalent radius 

    results              = struct; 
    d_pos_radius         = sqrt(areas/pi); 
    results.d_pos_radius = d_pos_radius; % equivalent radius (km)
    
% save other variables

    results.lon_grid = lon_grid; 
    results.lat_grid = lat_grid; 
    results.mask     = mask; 
    results.bathy    = bathy; 

% map the results 

    figure; hold on 

        easy_m_proj('mercator', lon_grid, lat_grid); 
        m_pcolor(lon_grid, lat_grid, d_pos_radius); 
        m_scatter(lon, lat, 'r', 'fill', 'MarkerEdgeColor', 'k'); 
        cbar = colorbar; cbar.Label.String = 'Geopositional Uncertainty (km)'; 
        m_gshhs_h('patch', [.75 .75 .75]); 
        labelformat([14 16])
        m_grid('box', 'fancy'); 

% runtime display 

    runtime = toc; 

    if runtime >= 60*60 
        runtime_h = runtime/(60*60); % runtime in hours 
        disp(['Runtime: ' char(string(runtime_h)) ' hours'])
    elseif runtime >= 60 & runtime < 60*60
        runtime_m = runtime/60;      % runtime in minutes 
        disp(['Runtime: ' char(string(runtime_m)) ' minutes'])
    elseif runtime < 60
        disp(['Runtime: ' char(string(runtime)) ' seconds'])
    end

end