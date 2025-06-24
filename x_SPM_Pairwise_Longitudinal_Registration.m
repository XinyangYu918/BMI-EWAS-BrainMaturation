clear;clc

addpath(genpath('/share/inspurStorage/home1/ISTBI_data/toolbox/spm12'));
path = '/share/inspurStorage/home1/ISTBI_data/IMAGEN/Smoke/T1_Pairwise_Longitudinal_Registration';
sub = dir(fullfile(path,'sub-*'));

finish = zeros(length(sub),1);
parfor i = 1:length(sub)
    disp(i)
    nii = dir(fullfile(path,sub(i).name,'avg*'));
    if length(nii) == 1
        finish(i) = i;
    end
end
sub(finish ~= 0) = [];

sub_job = cell(length(sub),1);
for i = 1:length(sub)
    disp(sub(i).name)
    BL_T1w = dir(fullfile(path,sub(i).name,'*baseline*T1w.nii'));
    FU_T1w = dir(fullfile(path,sub(i).name,'*followup*T1w.nii'));
    clear matlabbatch
    matlabbatch{1}.spm.tools.longit.pairwise.vols1 = {[fullfile(BL_T1w.folder,BL_T1w.name),',1']};
    matlabbatch{1}.spm.tools.longit.pairwise.vols2 = {[fullfile(FU_T1w.folder,FU_T1w.name),',1']};
    matlabbatch{1}.spm.tools.longit.pairwise.tdif = 5;
    matlabbatch{1}.spm.tools.longit.pairwise.noise = NaN;
    matlabbatch{1}.spm.tools.longit.pairwise.wparam = [0 0 100 25 100];
    matlabbatch{1}.spm.tools.longit.pairwise.bparam = 1000000;
    matlabbatch{1}.spm.tools.longit.pairwise.write_avg = 1;
    matlabbatch{1}.spm.tools.longit.pairwise.write_jac = 1;
    matlabbatch{1}.spm.tools.longit.pairwise.write_div = 1;
    matlabbatch{1}.spm.tools.longit.pairwise.write_def = 0;
    sub_job{i} = matlabbatch;
end

error = zeros(1,length(sub_job));
parfor i = 1:length(sub_job)
    try
        spm_jobman('run',sub_job{i});
    catch
        error(i) = i;
    end
end
error(error == 0) = [];