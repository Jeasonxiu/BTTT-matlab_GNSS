fid_out = fopen('daum_map','w');

% load 1215
% gd = NEV1(:,1:2);
% gd = NEV1(:,2:3);
n = length(gd);

for t = 1:n
    if gd(t,1) ~= -1
        fprintf(fid_out,'{latlng: new daum.maps.LatLng( %11.8f , %12.8f )},\n', gd(t,1), gd(t,2));
    end
end