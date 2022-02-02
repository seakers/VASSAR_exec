%%
clear
close
clc

%% Read contents of text file
file_loc = 'C:\\SEAK Lab\\SEAK Lab Github\\VASSAR\\VASSAR_exec_heur\\results\\full instrument duty cycle satisfaction architectures.txt';
data = regexp(fileread(file_loc),'\r?\n','split');

%% Analyze architectures
instruments_list = ["ACE_ORCA","ACE_POL","ACE_LID","CLAR_ERB","ACE_CPR","DESD_SAR","DESD_LID","GACM_VIS","GACM_SWIR","HYSP_TIR","POSTEPS_IRS","CNES_KaRIN"];
orbits_list = ["LEO-600-polar-NA","SSO-600-SSO-AM","SSO-600-SSO-DD","SSO-800-SSO-PM","SSO-800-SSO-DD"];

n_instr_all = size(instruments_list,2);
n_orbs_all = size(orbits_list,2);

for i = 2:size(data,2)-1
    current_arch = data{1,i};
    disp(strcat('Architecture: ',current_arch))
    n_instr = 0;
    for j = 1:n_instr_all
        for k = 1:n_orbs_all
            if strcmp(current_arch(n_orbs_all*(j-1)+k),'1')
                disp(strcat('Instrument: ',instruments_list(j),' in orbit: ',orbits_list(k)))
                n_instr = n_instr + 1;
            end
        end
    end
    disp(strcat('total number of instruments: ',num2str(n_instr)))
    disp('----------')
end
