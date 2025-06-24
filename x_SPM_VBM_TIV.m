clear;clc

path = '/share/inspurStorage/home1/ISTBI_data/IMAGEN/Smoke/T1_Followup2_VBM8';
sub = dir(fullfile(path,'sub-*'));

seg = cell(length(sub),1);
for i = 1:length(sub)
    disp(i)
    clear seg8
    seg8 = dir(fullfile(path,sub(i).name,'psub*seg8.txt'));
    seg{i} = fullfile(seg8.folder,seg8.name);
end
matlabbatch{1}.spm.tools.vbm8.tools.calcvol.data = seg;
matlabbatch{1}.spm.tools.vbm8.tools.calcvol.calcvol_name = 'Followup2_volumes.txt';

spm_jobman('run',matlabbatch);

Data = importdata('/share/inspurStorage/home1/ISTBI_data/IMAGEN/Smoke/T1_Followup2_VBM8/Followup2_volumes.txt');
T = table;
T.SubID = cellfun(@(x) replace(x,'sub-',''),{sub.name}','un',0);
T = [T,array2table(Data,'VariableNames',{'FU2_GM','FU2_WM','FU2_CSF','FU2_TIV'})];
writetable(T,'/share/inspurStorage/home1/ISTBI_data/IMAGEN/Smoke/T1_Pairwise_Longitudinal_Registration/FU2_TIV.csv');