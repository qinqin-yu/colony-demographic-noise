folders = dir2('/Volumes/PORTKeio3/');
for i = 2:length(folders)
    folders(i).name
    get_growth_layer_depth(['/Volumes/PORTKeio3/' folders(i).name]);
end