function get_colony_size(path_folder)
    %Colony size is in units of micron^2
    
    %path_folder = uigetdir('/Users/qinqinyu/Documents/hallatschek_lab/data_changed/','Main folder'); 
    %path_folder = '/Volumes/PORTKeioDat/20190910_1024_KQY3a';

    colony_mat_folder = [path_folder filesep 'colony_mat'];
    files = dir2(colony_mat_folder);
    nfiles = length(files);

    %Make folder to store results if it doesn't already exist
    if ~isdir([path_folder filesep 'summary_data'])
        mkdir([path_folder filesep 'summary_data']);
    end

    A = nan(96,1);

    %Conversion factor
    conv = 0.1240; 
    %0.1240 pixels = 1 um for 8x; 

    %Loop through all files
    for i=1:nfiles%1:nfiles
        %Figure out the well position
        filename = files(i).name;
        k = strfind(filename, '_');
        strain_name = filename(1:k-1);
        full_strain_name = filename(1:k-1);
        k2 = well2num(strain_name);

        %Load in the front points and fitted circle
        load([colony_mat_folder filesep filename]);

        %Create a new image that has the area inside the boundary filled in
        I = zeros(1040,1388);
        for j = 1:length(x)
            I(y(j)*conv,x(j)*conv) = 1;
        end
        I_fill = imfill(I,'holes');
        A(str2num(strain_name)+1) = sum(sum(I_fill))/conv^2;

    end
    csvwrite([path_folder filesep 'summary_data/colony_area.csv'], A)
end
