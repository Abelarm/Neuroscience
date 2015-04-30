function NameFiles = GV_getNameFiles_Salerno(Subject,varargin)
% Multivariate decoding of single button press
% project on hemodynamic modeling for Multi-Voxel Pattern Analysis (MVPA)
% Giancarlo Valente, Maastricht, 2014
% giancarlo.valente@maastrichtuniversity.nl
if isempty(varargin)
    MaskPostFix = '_Pre&PostCGsmall';
else
    MaskPostFix = varargin{1};
end

switch lower(Subject)
    case 'gv'
        dirdata = '/Users/Luigi/Desktop/progetto_neuroscienze/Data/GV/';
        NameFiles.NameMask = fullfile(dirdata,['GV',MaskPostFix,'.msk']);
        
                NameFiles.NameVtc = {fullfile(dirdata,'GV_20130417_03_MOTORMVPA_Slow_v2_SCCAI2_3DMCTS_THPGLMF9c_TAL.vtc'),...
                    fullfile(dirdata,'GV_20130417_06_MOTORMVPA_Slow_v1_SCCAI2_3DMCTS_THPGLMF9c_TAL.vtc')};
                NameFiles.NamePrt = {fullfile(dirdata,'GV_MotorMVPA_03_SlowER_v2_vols.prt'),...
                    fullfile(dirdata,'GV_MotorMVPA_06_SlowER_v1_vols.prt')};
                NameFiles.NameSave = fullfile(dirdata,['GV',MaskPostFix,'_SLOW']);
       
        
        
        
        
    case  'az'
        dirdata = '/Users/Luigi/Desktop/progetto_neuroscienze/Data/AZ/';
        NameFiles.NameMask = fullfile(dirdata,['AZ',MaskPostFix,'.msk']);
       
                NameFiles.NameVtc = {fullfile(dirdata,'AZ_20130315_04_MOTORMVPA_SE_v2_SCCAI2_3DMCTS_THPGLMF9c_TAL.vtc'),...
                    fullfile(dirdata,'AZ_20130315_05_MOTORMVPA_SE_v1_SCCAI2_3DMCTS_THPGLMF9c_TAL.vtc')};
                NameFiles.NamePrt = {fullfile(dirdata,'AZ_20130315_MOTORMVPA_05_SlowER_v2_vols.prt'),...
                    fullfile(dirdata,'AZ_20130315_MOTORMVPA_06_SlowER_v1_vols.prt')};
                NameFiles.NameSave = fullfile(dirdata,['AZ',MaskPostFix,'_SLOW']);
          
        
    case 'je'
        
        dirdata = '/Users/Luigi/Desktop/progetto_neuroscienze/Data/JE/';
        NameFiles.NameMask = fullfile(dirdata,['JE',MaskPostFix,'.msk']);
       
                NameFiles.NameVtc = {fullfile(dirdata,'JE_20130313_MOTORMVPA_02_SE_v1_SCCAI2_3DMCTS_THPGLMF9c_TAL.vtc'),...
                    fullfile(dirdata,'JE_20130313_MOTORMVPA_07_SE_v2_SCCAI2_3DMCTS_THPGLMF9c_TAL.vtc')};
                NameFiles.NamePrt = {fullfile(dirdata,'JE_20130313_MotorMVPA_02_SlowER_v1_vols.prt'),...
                    fullfile(dirdata,'JE_20130313_MotorMVPA_06_SlowER_v2_vols.prt')};
                NameFiles.NameSave = fullfile(dirdata,['JE',MaskPostFix,'_SLOW']);
            

        
    otherwise
            error('Wrong subject name')
        
        
end