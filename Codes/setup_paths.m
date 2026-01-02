function setup_paths()
%SETUP_PATHS Add required paths for HVSR joint inversion
%
% Call this function at the start of any script to ensure all required
% MATLAB functions and CPS binaries are accessible.
%
% USAGE:
%   setup_paths();

    base_path = '/Users/birotimi/Library/CloudStorage/OneDrive-SyracuseUniversity/Desktop/PHD/HVSR-project/HVSR-Joint-Inversion/Codes';
    cps_bin = '/Users/birotimi/Library/CloudStorage/OneDrive-SyracuseUniversity/Desktop/PHD/HVSR-project/HVSR-Joint-Inversion/AGU_models/bin_v3.30';
    hvf_bin = '/Users/birotimi/Library/CloudStorage/OneDrive-SyracuseUniversity/Desktop/PHD/HVSR-project/HVSR-Joint-Inversion/AGU_models/bin_HV_DFA';
    func_path = '/Users/birotimi/Library/CloudStorage/OneDrive-SyracuseUniversity/Desktop/PHD/HVSR-project/HVSR-Joint-Inversion/AGU_models/functions';
    
    addpath(genpath(base_path));
    addpath(func_path);
    
    current_path = getenv('PATH');
    if ~contains(current_path, hvf_bin)
        setenv('PATH', [hvf_bin, ':', cps_bin, ':', current_path]);
    end
end