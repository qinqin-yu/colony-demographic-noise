deltar = 1000;%50;   %microns

% folders = dir2('/Volumes/PORTKeio3/');
% for i = 2:length(folders)
%     folders(i).name
%     get_msd_extrap_front(['/Volumes/PORTKeio3/' folders(i).name], deltar);
% end

%local_front_roughness_fit('/Volumes/PORTKeio3/20200213_2030_mreBc', 1)   %starts at file 23
get_msd_extrap_front('/Users/qinqinyu/Documents/hallatschek_lab/data_changed/20180921_ecoli_beads_timelapse_second_try/20180920_bOH1_bOH2_red_1um_16x_1hr_2h-23h/bf/', deltar);



%get_msd_extrap_front('/Volumes/PORTKeio3/20200214_0900_DE5b', deltar);
%get_msd_extrap_front('/Volumes/PORTKeioDat/20190912_2100_KQY2b', deltar)
%get_msd_extrap_front('/Volumes/PORTKeioDat/20190913_0930_KQY2c', deltar)
%get_msd_extrap_front('/Volumes/PORTKeioDat/20190913_2130_KQY1c', deltar)
