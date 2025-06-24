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
    matlabbatch{1}.spm.tools.vbm8.estwrite.data = {[fullfile(Avg_T1w.folder,Avg_T1w.name),',1']};
    matlabbatch{1}.spm.tools.vbm8.estwrite.opts.tpm = {'/share/inspurStorage/home1/ISTBI_data/toolbox/spm8/toolbox/Seg/TPM.nii'};
    matlabbatch{1}.spm.tools.vbm8.estwrite.opts.ngaus = [2 2 2 3 4 2];
    matlabbatch{1}.spm.tools.vbm8.estwrite.opts.biasreg = 0.0001;
    matlabbatch{1}.spm.tools.vbm8.estwrite.opts.biasfwhm = 60;
    matlabbatch{1}.spm.tools.vbm8.estwrite.opts.affreg = 'mni';
    matlabbatch{1}.spm.tools.vbm8.estwrite.opts.warpreg = 4;
    matlabbatch{1}.spm.tools.vbm8.estwrite.opts.samp = 3;
    matlabbatch{1}.spm.tools.vbm8.estwrite.extopts.dartelwarp.normhigh.darteltpm = {'/share/inspurStorage/home1/ISTBI_data/toolbox/spm8/toolbox/vbm8/Template_1_IXI550_MNI152.nii'};
    matlabbatch{1}.spm.tools.vbm8.estwrite.extopts.sanlm = 2;
    matlabbatch{1}.spm.tools.vbm8.estwrite.extopts.mrf = 0.15;
    matlabbatch{1}.spm.tools.vbm8.estwrite.extopts.cleanup = 1;
    matlabbatch{1}.spm.tools.vbm8.estwrite.extopts.print = 0;
    matlabbatch{1}.spm.tools.vbm8.estwrite.output.GM.native = 0;
    matlabbatch{1}.spm.tools.vbm8.estwrite.output.GM.warped = 0;
    matlabbatch{1}.spm.tools.vbm8.estwrite.output.GM.modulated = 0;
    matlabbatch{1}.spm.tools.vbm8.estwrite.output.GM.dartel = 2;
    matlabbatch{1}.spm.tools.vbm8.estwrite.output.WM.native = 0;
    matlabbatch{1}.spm.tools.vbm8.estwrite.output.WM.warped = 0;
    matlabbatch{1}.spm.tools.vbm8.estwrite.output.WM.modulated = 0;
    matlabbatch{1}.spm.tools.vbm8.estwrite.output.WM.dartel = 2;
    matlabbatch{1}.spm.tools.vbm8.estwrite.output.CSF.native = 0;
    matlabbatch{1}.spm.tools.vbm8.estwrite.output.CSF.warped = 0;
    matlabbatch{1}.spm.tools.vbm8.estwrite.output.CSF.modulated = 0;
    matlabbatch{1}.spm.tools.vbm8.estwrite.output.CSF.dartel = 0;
    matlabbatch{1}.spm.tools.vbm8.estwrite.output.bias.native = 0;
    matlabbatch{1}.spm.tools.vbm8.estwrite.output.bias.warped = 0;
    matlabbatch{1}.spm.tools.vbm8.estwrite.output.bias.affine = 0;
    matlabbatch{1}.spm.tools.vbm8.estwrite.output.label.native = 0;
    matlabbatch{1}.spm.tools.vbm8.estwrite.output.label.warped = 0;
    matlabbatch{1}.spm.tools.vbm8.estwrite.output.label.dartel = 0;
    matlabbatch{1}.spm.tools.vbm8.estwrite.output.jacobian.warped = 0;
    matlabbatch{1}.spm.tools.vbm8.estwrite.output.warps = [0 0];
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

