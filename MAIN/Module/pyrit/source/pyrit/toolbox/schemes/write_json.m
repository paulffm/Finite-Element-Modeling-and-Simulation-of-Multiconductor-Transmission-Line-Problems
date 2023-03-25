% Script that iterates over the possible integration orders of the dunavant
% integration schemes and writes a json file for each integration order. 
% They can be combined into a single json file with the python script
% 'combine_json_files.py'. 

for k = 1:20
    [xyz, stats] = set_dunavant_standard(k);

    jsonencode(xyz);
    a.ID = 2;
    a.xyz = xyz;
    
    text = jsonencode(a);

    fileID = fopen(['dunavant_' num2str(k) '.json'],'w');
    fprintf(fileID,text);
    fclose(fileID);
end
