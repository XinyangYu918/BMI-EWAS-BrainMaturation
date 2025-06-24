clear;clc

path = '/share/inspurStorage/home1/ISTBI_data/IMAGEN/Smoke/T1_Pairwise_Longitudinal_Registration';
sub = dir(fullfile(path,'sub-*'));

parfor i = 1:length(sub)
    disp(i)
    jd = dir(fullfile(path,sub(i).name,'jd*'));
    gray = dir(fullfile(path,sub(i).name,'p1*'));
    Gray_JD = fullfile(path,sub(i).name,['Gray_JD_',sub(i).name,'_T1w.nii.gz']);
    unix(['fslmaths ',fullfile(jd.folder,jd.name),' -mul ',fullfile(gray.folder,gray.name),' ',Gray_JD]);
    gunzip(Gray_JD);
    unix(['rm ',Gray_JD]);
end

clear;clc

path = '/share/inspurStorage/home1/ISTBI_data/IMAGEN/Smoke/T1_Pairwise_Longitudinal_Registration';
sub = dir(fullfile(path,'sub-*'));
addpath(genpath('/share/inspurStorage/home1/ISTBI_data/toolbox/spm8'));

sub_job = cell(length(sub),1);
for i = 1:length(sub)
    disp(sub(i).name)
    Flowfields = dir(fullfile(path,sub(i).name,'u_*.nii'));
    Gray_JD = dir(fullfile(path,sub(i).name,'Gray_JD*.nii'));
    clear matlabbatch
    matlabbatch{1}.spm.tools.dartel.mni_norm.template = {'/share/inspurStorage/home1/ISTBI_data/IMAGEN/Smoke/T1_Pairwise_Longitudinal_Registration/Template_6.nii'};
    matlabbatch{1}.spm.tools.dartel.mni_norm.data.subjs.flowfields = {fullfile(Flowfields.folder,Flowfields.name)};
    matlabbatch{1}.spm.tools.dartel.mni_norm.data.subjs.images = {{fullfile(Gray_JD.folder,Gray_JD.name)}};
    matlabbatch{1}.spm.tools.dartel.mni_norm.vox = [NaN NaN NaN];
    matlabbatch{1}.spm.tools.dartel.mni_norm.bb = [NaN NaN NaN
                                                NaN NaN NaN];
    matlabbatch{1}.spm.tools.dartel.mni_norm.preserve = 0;
    matlabbatch{1}.spm.tools.dartel.mni_norm.fwhm = [8 8 8];
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