function [cbf]=DSC_mri_cbf(conc,aif,mask,options,app)

% Calculates parametric maps of Cerebral Blood Flow (CBF) for a subject
%
% Input parameters - conc (4D Matrix), contains the DSC concentration trends of all voxels.
%                  - aif, concentration trend in the site chosen as arterial
%                  - mask (3D Matrix), contains the matrix used to mask the brain volume to be analyzed
%
% Options is the struct that contains the method options, the significant ones are:
%
% The calculation of the CBF parameter requires a deconvolution operation, 
% in the toolbox there are some methodologies suitable for this purpose: SVD and cSVD (Block-circulant version).
%
% options.deconv.method - Must be a cell array containing the names of the algorithms 
% intended to be used for the analysis. Ex. SVD, cSVD, or SS
%
%                       For each method, a struct with the specific parameters of the method 
%                       itself must be inserted in the options, as in this example:
%
%                       options.deconv.<method name>.<par_1>
%                       options.deconv.<method name>.<...>
%                       options.deconv.<method name>.<par_n>
%
% options.time - Represents the time vector of the DSC exam (each sample represents the acquisition of an entire brain volume)
%
% options.display - level 1 Shows the processing progress, level 2 Also provides information on the parameters set for the algorithms used
%
% Output parameters: 
% cbf - different sub-structs, one for each method chosen to be used, each sub-struct contains a map field that distinguishes the calculated cbv map, for example:
%      cbf.<method name>.map
%
%      residual, a 4D matrix that contains the residuals (it's an optional field that must be requested in the options: parameter options.deconv.<method name>.residual) for example:
%      cbf.<method name>.residual



method={'SVD';'cSVD';'oSVD'};


for alg=1:size(options.deconv.method,1)
    switch options.deconv.method{alg,:}
        case method{1,:} %SVD
            cbf.svd=DSC_mri_SVD(conc,aif,mask.data,options,app);
        case method{2,:} %cSVD
            cbf.csvd=DSC_mri_cSVD(conc,aif,mask.data,options,app);
        case method{3,:} %oSVD          
            cbf.osvd=DSC_mri_oSVD(conc,aif,mask.data,options,app);
            
        %AGGIORNARE SE SI INTENDE AGGIUNGERE UN NUOVO METODO
        %case method{n,:} 
        %cbf.<nome metodo>=DSC_mri_<nome_metodo>(...);
        
        otherwise %metodo non riconosciuto
            str_errore='';
            for m=1:size(method,1)
                %DA SISTEMARE
                str_errore=strcat(str_errore,method{m,:},' or ');
            end
            error(['Deconvolution method not recognized. Use: ' str_errore])
    end
end