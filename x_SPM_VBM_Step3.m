clear;clc

rmpath(genpath('/share/inspurStorage/home1/ISTBI_data/toolbox/spm12'));
addpath(genpath('/share/inspurStorage/home1/ISTBI_data/toolbox/spm8'));
path = '/share/inspurStorage/home1/ISTBI_data/IMAGEN/Smoke/T1_Pairwise_Longitudinal_Registration';
sub = dir(fullfile(path,'sub-*'));

sub_job = cell(length(sub),1);
for i = 1:length(sub)
    disp(sub(i).name)
    Avg_T1w = dir(fullfile(path,sub(i).name,'avg*T1w.nii'));
    clear matlabbatch
    matlabbatch{1}.spm.tools.vbm8.write.data = {[fullfile(Avg_T1w.folder,Avg_T1w.name),',1']};
    matlabbatch{1}.spm.tools.vbm8.write.extopts.dartelwarp.normhigh.darteltpm = {'/share/inspurStorage/home1/ISTBI_data/IMAGEN/Smoke/T1_Pairwise_Longitudinal_Registration/Template_1.nii'};
    matlabbatch{1}.spm.tools.vbm8.write.extopts.sanlm = 2;
    matlabbatch{1}.spm.tools.vbm8.write.extopts.mrf = 0.15;
    matlabbatch{1}.spm.tools.vbm8.write.extopts.cleanup = 1;
    matlabbatch{1}.spm.tools.vbm8.write.extopts.print = 0;
    matlabbatch{1}.spm.tools.vbm8.write.output.GM.native = 1;
    matlabbatch{1}.spm.tools.vbm8.write.output.GM.warped = 1;
    matlabbatch{1}.spm.tools.vbm8.write.output.GM.modulated = 2;
    matlabbatch{1}.spm.tools.vbm8.write.output.GM.dartel = 0;
    matlabbatch{1}.spm.tools.vbm8.write.output.WM.native = 1;
    matlabbatch{1}.spm.tools.vbm8.write.output.WM.warped = 1;
    matlabbatch{1}.spm.tools.vbm8.write.output.WM.modulated = 2;
    matlabbatch{1}.spm.tools.vbm8.write.output.WM.dartel = 0;
    matlabbatch{1}.spm.tools.vbm8.write.output.CSF.native = 1;
    matlabbatch{1}.spm.tools.vbm8.write.output.CSF.warped = 1;
    matlabbatch{1}.spm.tools.vbm8.write.output.CSF.modulated = 2;
    matlabbatch{1}.spm.tools.vbm8.write.output.CSF.dartel = 0;
    matlabbatch{1}.spm.tools.vbm8.write.output.bias.native = 0;
    matlabbatch{1}.spm.tools.vbm8.write.output.bias.warped = 1;
    matlabbatch{1}.spm.tools.vbm8.write.output.bias.affine = 0;
    matlabbatch{1}.spm.tools.vbm8.write.output.label.native = 0;
    matlabbatch{1}.spm.tools.vbm8.write.output.label.warped = 0;
    matlabbatch{1}.spm.tools.vbm8.write.output.label.dartel = 0;
    matlabbatch{1}.spm.tools.vbm8.write.output.jacobian.warped = 1;
    matlabbatch{1}.spm.tools.vbm8.write.output.warps = [1 0];
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
