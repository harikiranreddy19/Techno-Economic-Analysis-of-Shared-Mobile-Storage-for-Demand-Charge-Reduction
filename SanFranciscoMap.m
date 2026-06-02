%%
clc; clear all; close all;
% url = 'https://raw.githubusercontent.com/OpenDataDE/State-zip-code-GeoJSON/master/ca_california_zip_codes_geo.min.json';
% geojson_str = webread(url);
% geojson_data = loadjson(geojson_str);
% features = geojson_data.features;
% save('CaliforniaFeatures')
load('CaliforniaFeatures.mat')
%%
AllZIP = [94102,94103,94104,94105,94107,94108,94109,94110,94111,94112,94114,94115,...,
       94116,94117,94118,94121,94122,94123,94124,94127,94131,94132,94133,...,
       94134,94129,94130,94158];
for i = 1:numel(features)
    [check, ind] = ismember(double(string(features{i}.properties.ZCTA5CE10)), AllZIP);
    if check
        SFGeometry{ind,1} = features{i}.geometry.coordinates;
    end
end
save('SFMap','AllZIP','SFGeometry')
