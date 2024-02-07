% I'm tired of writing m_proj(...) and having to calculate the limits 
% myself. This function lets you just input: 
%
%   easy_m_proj(projection, lon, lat) 
%
% where lon, lat is whatever grid or vector of lon/lat values you are
% using, the function simply does the following: 
% 
%   m_proj(projection, 'long', [min(lon) max(lon)], 'lat', [min(lat)
%   max(lat)]) 
%
% which is often kind of annoying because you're working with arrays and
% you get this long ugly expression. This is merely laziness. 
%
% If you wish to project with some margin along the edge positions of your
% data, enter a fourth argument with the number of degrees of margin you
% want (one value for both latitude and longitude). 

function easy_m_proj(projection, lon, lat, varargin) 

switch nargin
    case 3

        m_proj(projection, 'long', [min(lon, [], 'all') max(lon, [], 'all')], ...
                           'lat',  [min(lat, [], 'all') max(lat, [], 'all')]); 

    case 4 
        
        margin = varargin{1}; 
        m_proj(projection, 'long', [min(lon, [], 'all')-margin max(lon, [], 'all')+margin], ...
                           'lat',  [min(lat, [], 'all')-margin max(lat, [], 'all')+margin]); 

end