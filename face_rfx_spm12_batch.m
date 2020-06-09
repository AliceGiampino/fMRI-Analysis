% This batch script analyses the Face Group dataset available from the SPM
% website:
%   http://www.fil.ion.ucl.ac.uk/spm/data/face_rfx/
% as described in the SPM manual:
%   http://www.fil.ion.ucl.ac.uk/spm/doc/spm12_manual.pdf#Chap:data:faces_group
%__________________________________________________________________________
% Copyright (C) 2014 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: face_rfx_spm12_batch.m 16 2014-09-30 12:53:58Z guillaume $


% Directory containing Face Group data
%--------------------------------------------------------------------------
data_path = fileparts(mfilename('E:\PhD\Statistical Inference II\HW5\'));
if isempty(data_path), data_path = pwd; end
fprintf('%-40s:', 'Downloading Face dataset...');
urlwrite('http://www.fil.ion.ucl.ac.uk/spm/download/data/face_rfx/face_rfx.zip','face_rfx.zip');
unzip(fullfile(data_path,'face_rfx.zip'));
data_path = fullfile(data_path,'face_rfx');
fprintf(' %30s\n', '...done');

% Initialise SPM 
%--------------------------------------------------------------------------
spm('Defaults','fMRI');
spm_jobman('initcfg');
%spm_get_defaults('cmdline',true);


% Canonical HRF
%==========================================================================
f = spm_select('FPList', fullfile(data_path,'cons_can'), '^con_.*\.img$') ;

clear matlabbatch

matlabbatch{1}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.parent = cellstr(data_path);
matlabbatch{1}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.name = 'canonical';

matlabbatch{2}.spm.stats.factorial_design.dir = cellstr(fullfile(data_path,'canonical'));
matlabbatch{2}.spm.stats.factorial_design.des.t1.scans = cellstr(f);

matlabbatch{3}.spm.stats.fmri_est.spmmat = cellstr(fullfile(data_path,'canonical','SPM.mat'));

matlabbatch{4}.spm.stats.con.spmmat = cellstr(fullfile(data_path,'canonical','SPM.mat'));
matlabbatch{4}.spm.stats.con.consess{1}.fcon.name = 'Faces vs Baseline: Canonical HRF';
matlabbatch{4}.spm.stats.con.consess{1}.fcon.weights = [1];

matlabbatch{5}.spm.stats.results.spmmat = cellstr(fullfile(data_path,'canonical','SPM.mat'));
matlabbatch{5}.spm.stats.results.conspec.contrasts = Inf;
matlabbatch{5}.spm.stats.results.conspec.threshdesc = 'FWE';

spm_jobman('run',matlabbatch);


% Informed Basis Set
%==========================================================================
f = spm_select('FPList', fullfile(data_path,'cons_informed'), '^con_.*\.img$') ;

clear matlabbatch

matlabbatch{1}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.parent = cellstr(data_path);
matlabbatch{1}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.name = 'informed';

matlabbatch{2}.spm.stats.factorial_design.dir = cellstr(fullfile(data_path,'informed'));
matlabbatch{2}.spm.stats.factorial_design.des.fd.fact.name = 'Basis';
matlabbatch{2}.spm.stats.factorial_design.des.fd.fact.levels = 3;
matlabbatch{2}.spm.stats.factorial_design.des.fd.fact.dept = 1;
matlabbatch{2}.spm.stats.factorial_design.des.fd.icell(1).levels = 1;
matlabbatch{2}.spm.stats.factorial_design.des.fd.icell(1).scans = cellstr(f(1:12,:));
matlabbatch{2}.spm.stats.factorial_design.des.fd.icell(2).levels = 2;
matlabbatch{2}.spm.stats.factorial_design.des.fd.icell(2).scans = cellstr(f(13:24,:));
matlabbatch{2}.spm.stats.factorial_design.des.fd.icell(3).levels = 3;
matlabbatch{2}.spm.stats.factorial_design.des.fd.icell(3).scans = cellstr(f(25:36,:));

matlabbatch{3}.spm.stats.fmri_est.spmmat = cellstr(fullfile(data_path,'informed','SPM.mat'));

matlabbatch{4}.spm.stats.con.spmmat = cellstr(fullfile(data_path,'informed','SPM.mat'));
matlabbatch{4}.spm.stats.con.consess{1}.fcon.name = 'Faces vs Baseline: Informed';
matlabbatch{4}.spm.stats.con.consess{1}.fcon.weights = eye(3);
matlabbatch{4}.spm.stats.con.consess{2}.tcon.name = 'Informed: t canonical contrast';
matlabbatch{4}.spm.stats.con.consess{2}.tcon.weights = 1;
matlabbatch{4}.spm.stats.con.consess{3}.fcon.name = 'Ftemp';
matlabbatch{4}.spm.stats.con.consess{3}.fcon.weights = [0 1 0];
matlabbatch{4}.spm.stats.con.consess{4}.fcon.name = 'Fdisp';
matlabbatch{4}.spm.stats.con.consess{4}.fcon.weights = [0 0 1];

matlabbatch{5}.spm.stats.results.spmmat = cellstr(fullfile(data_path,'informed','SPM.mat'));
matlabbatch{5}.spm.stats.results.conspec(1).contrasts = Inf;
matlabbatch{5}.spm.stats.results.conspec(1).threshdesc = 'FWE';

spm_jobman('run',matlabbatch);


% FIR Basis Set
%==========================================================================
clear matlabbatch

matlabbatch{1}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.parent = cellstr(data_path);
matlabbatch{1}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.name = 'FIR';

matlabbatch{2}.spm.stats.factorial_design.dir = cellstr(fullfile(data_path,'FIR'));
matlabbatch{2}.spm.stats.factorial_design.des.fd.fact.name = 'TimeBin';
matlabbatch{2}.spm.stats.factorial_design.des.fd.fact.levels = 12;
matlabbatch{2}.spm.stats.factorial_design.des.fd.fact.dept = 1;
for i=1:12
	matlabbatch{2}.spm.stats.factorial_design.des.fd.icell(i).levels = i;
	f = spm_select('FPList', fullfile(data_path,'cons_fir'), sprintf('^con_fir_bin%02d.*\\.img$',i));
	matlabbatch{2}.spm.stats.factorial_design.des.fd.icell(i).scans = cellstr(f);
end

matlabbatch{3}.spm.stats.fmri_est.spmmat = cellstr(fullfile(data_path,'FIR','SPM.mat'));

matlabbatch{4}.spm.stats.con.spmmat = cellstr(fullfile(data_path,'FIR','SPM.mat'));
matlabbatch{4}.spm.stats.con.consess{1}.fcon.name = 'Faces vs Baseline: FIR';
matlabbatch{4}.spm.stats.con.consess{1}.fcon.weights = eye(12);

xBF.dt = 1;
xBF.name = 'hrf (with time and dispersion derivatives)';
xBF.length = 32;
xBF.order = 1;
xBF = spm_get_bf(xBF);
all = xBF.bf(2:2:24,:)';
can = all(1,:);
tem = all(2,:);
dis = all(3,:);
disp(corrcoef(all')');
nullcan = eye(12) - pinv(can)*can;
nullall = eye(12) - pinv(all)*all;

matlabbatch{4}.spm.stats.con.consess{2}.tcon.name = 'Can-weighted FIR';
matlabbatch{4}.spm.stats.con.consess{2}.tcon.convec = can;
matlabbatch{4}.spm.stats.con.consess{3}.fcon.name = 'Null space of canonical HRF';
matlabbatch{4}.spm.stats.con.consess{3}.fcon.weights = nullcan;
matlabbatch{4}.spm.stats.con.consess{4}.fcon.name = 'Null space of informed basis';
matlabbatch{4}.spm.stats.con.consess{4}.fcon.weights = nullall;

matlabbatch{5}.spm.stats.results.spmmat = cellstr(fullfile(data_path,'FIR','SPM.mat'));
matlabbatch{5}.spm.stats.results.conspec.contrasts = Inf;
matlabbatch{5}.spm.stats.results.conspec.threshdesc = 'FWE';

spm_jobman('run',matlabbatch);
