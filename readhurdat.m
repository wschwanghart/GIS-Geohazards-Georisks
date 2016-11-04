function S = readhurdat(file,geom)
% read best track hurricane database (hurdat2) into mapping structure
% 
% Syntax 
%
%     S = readhurdat(geom,file)
%
% Description
%
%     This function reads the best track hurricane database (hurdat2) into 
%     a mapping structure array. The data can be downloaded here:
%     http://www.nhc.noaa.gov/data/ (see section Best Track Data (HURDAT2))
%     The mapping structure can be exported as shapefile and used in a GIS
%     environment. shapewrite requires the Mapping Toolbox.
%     
% Input arguments
%
%     file      path to file
%     geom      type of the mapping structure.
%               'line2'  creates a line for each track. The resulting 
%               mapping structure contains several fields including maximum
%               windspeed, minimum pressure, etc.
%               'line' creates a line for each track segment. The
%               resulting structure contains for each segment windspeed,
%               pressure, date, hour, as well as wind radii maximum extents
%               for each quadrant. Note that all units are changed from
%               knots to km/h and from miles to km. This is the default.
%               
% Output arguments
%
%     S         mapping structure
%
% Example
%
%     S = readhurdat('hurdat2-1851-2015-070616.txt','line2');
%     geoshow(S)
%
% See also: Understanding Vector Geodata in the MATLAB documentation
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 23. February, 2013

narginchk(1,2);

if nargin == 1;
    geom = 'line';
else
    geom = validatestring(geom,{'point','line','line2'});
end

% get file identifier
fid = fopen(file);

switch geom
    case 'point'
        dd = 0;
    case 'line'
        dd = 1;
end

try
    
% read first line
tline = fgetl(fid);

% counter
c = 0;
% second counter
c2 = 0;

% go through the file as long there is one
while ischar(tline)
    % increase counter by one
    c = c+1;
    % split the string using commas into own strings
    HeaderLine = regexp(tline,',','split');
    
    % there may be extra lines at the bottom of the file
    if numel(HeaderLine) <2
        return
    end
    
    switch geom
        case 'line2'
            
            S(c).Geometry = 'Line';
            S(c).ID = c;
            S(c).Basin = HeaderLine{1}(1:2);
            S(c).ATCFn = HeaderLine{1}(3:4);
            S(c).Year  = HeaderLine{1}(5:8);
            S(c).Name  = strtrim(HeaderLine{2});
            
            S(c).nrtrackentries = str2double(HeaderLine{3});
            
            % read data into a cell array
            C  = cell(S(c).nrtrackentries,21);
            for r = 1:S(c).nrtrackentries;
                tline = fgetl(fid);
                C(r,:) = regexp(tline,',','split');
            end
            
            S(c).StaDate  = C{1,1};
            S(c).StaHours = C{1,2};
            S(c).EndDate  = C{end,1};
            S(c).EndHours = C{end,1};
            
            S(c).Lat = getlat(C(:,5));
            S(c).Lon = getlon(C(:,6));
            
            % Maximum sustained wind
            S(c).MaxWS = knots2kmh(max(str2double(C(:,7))));
            
            % Minimum pressure
            pr = str2double(C(:,8));
            pr(pr==-999) = nan;
            if all(isnan(pr))
                
                S(c).MinPr = -999;
            else
                S(c).MinPr = min(pr);
            end
            
            
            
        case {'point','line'}
            
            ID = c;
            Basin = HeaderLine{1}(1:2);
            ATCFn = HeaderLine{1}(3:4);
            Year  = HeaderLine{1}(5:8);
            Name  = strtrim(HeaderLine{2});
            nrtrackentries = str2double(HeaderLine{3});
            
            C  = cell(nrtrackentries,21);
            for r = 1:nrtrackentries;
                tline = fgetl(fid);
                C(r,:) = regexp(tline,',','split');
            end
            
            for r = 1:(nrtrackentries-dd);
                c2 = c2+1;
                S(c2).Geometry = geom;
                S(c2).ID = ID;
                S(c2).Date  = C{r,1};
                S(c2).Hours = C{r,2};
                S(c2).Basin = Basin;
                S(c2).ATCFn = ATCFn;
                S(c2).Year = Year;
                S(c2).Name = Name;
                
                switch geom
                    case 'point'
                        S(c2).Lat = getlat(C(r,5));
                        S(c2).Lon = getlon(C(r,6));
                    case 'line'
                        S(c2).Lat = getlat(C([r r+1],5));
                        S(c2).Lon = getlon(C([r r+1],6));
                end
                
                S(c2).WS = knots2kmh(str2double(C(r,7)));
                
                % Minimum pressure
                pr = str2double(C(r,8));
                pr(pr==-999) = nan;
               
                S(c2).Pr = pr;
                S(c2).RecId = strtrim(C{r,3});
                S(c2).status = strtrim(C{r,4});
                
                S(c2).WR34_NE = knots2kmh(str2double(C(r,9)));
                S(c2).WR34_SE = knots2kmh(str2double(C(r,10)));
                S(c2).WR34_SW = knots2kmh(str2double(C(r,11)));
                S(c2).WR34_NW = knots2kmh(str2double(C(r,12)));
                
                S(c2).WR50_NE = knots2kmh(str2double(C(r,13)));
                S(c2).WR50_SE = knots2kmh(str2double(C(r,14)));
                S(c2).WR50_SW = knots2kmh(str2double(C(r,15)));
                S(c2).WR50_NW = knots2kmh(str2double(C(r,16)));
                
                S(c2).WR64_NE = knots2kmh(str2double(C(r,17)));
                S(c2).WR64_SE = knots2kmh(str2double(C(r,18)));
                S(c2).WR64_SW = knots2kmh(str2double(C(r,19)));
                S(c2).WR64_NW = knots2kmh(str2double(C(r,20)));
            
            end

    
    end
            
    tline = fgetl(fid);
end

catch ME
    fclose(fid);
    rethrow(ME);
    
end

end

function lat = getlat(str)
    if iscell(str)
        lat = zeros(size(str));
        for rr = 1:numel(str);
            lat(rr) = getlat(str{rr});
        end
        return
    end
    if strcmp(str(end),'N')
        sgn = 1;
    else
        sgn = -1;
    end
lat = sgn*str2double(str(1:end-1));
end

function lon = getlon(str)
    if iscell(str)
        lon = zeros(size(str));
        for rr = 1:numel(str);
            lon(rr) = getlon(str{rr});
        end
        return
    end
    if strcmp(str(end),'W')
        lon = -(str2double(str(1:end-1)));
        sgn = -1;
    else
        lon = str2double(str(1:end-1));
        sgn = 1;
    end
end

function a = knots2kmh(b)
% check for missing values
I = b==-99 | b==-999;
a = 1.852*b;
a(I) = nan;
end
