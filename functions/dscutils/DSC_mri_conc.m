function [conc,S0map,bolus]=DSC_mri_conc(volumes,mask,options)

% Computes concentrations and S0 maps in DSC-MRI examinations.
%
% Input parameters:
% - volumes (4D matrix) containing the time course of the DSC signal for all voxels.
% - mask (3D matrix) containing the matrix to mask non-interest voxels for the study
%
% Options is the struct containing method options, the significant ones are:
%
% par.kvoi - Proportionality constant for calculating tracer concentration in VOI, by default it is
% considered unknown and set to 1.
%
% S0 - parameter series for determining the S0 calculation threshold
% S0.nSamplesMin - number of samples definitely acquired
% before injection
% S0.nSamplesMax - number of samples after which I stop anyway
% S0.thresh; - Add a sample if its difference from the mean is
% less than the threshold
%
% Output parameters:
% - conc: 4D matrix of concentrations
% - S0map: 3D matrix of S0s

%#ok<*AGROW>

[S0map,bolus] = DSC_mri_S0(volumes, mask, options);

conc = zeros(size(volumes));

ind=find(mask);
k=options.nR*options.nC*options.nS;

for t=1:options.nT
    step1=volumes(ind+k*(t-1))./S0map(ind);
    conc(ind+k*(t-1))=-(options.par.kvoi/options.te).*log(step1);
end


end


function [S0map,bolus]=DSC_mri_S0(volumes,mask,options)

% The function calculates the injection bolus timing and S0 from the data.
% 1) On the average time course, it computes the bolus injection timing: it calculates
% the mean of the first n samples and adds the n+1th sample if its percentage difference
% from the mean is below a given threshold.
% 2) It calculates S0 as the mean of the first n samples for all voxels.

nSamplesMin=options.S0.nSamplesMin;
nSamplesMax=options.S0.nSamplesMax;
thresh=options.S0.thresh;
mean_signal=zeros(size(options.time));

for s=1:options.nS
    for t=1:options.nT
        indMask=find(mask(:,:,s));
        mean_signal(s,t)=mean(volumes(indMask+(options.nR*options.nC*options.nS)*(t-1)+(options.nS*(s-1))));
        clear('temp','indMask');
    end
end

for s=1:options.nS

    ciclo=true;
    pos=nSamplesMin;
    
    while ciclo
    
        mean_val  =mean(mean_signal(s,1:pos));
        if abs((mean_val-mean_signal(s,pos+1))/mean_val)<thresh
            pos = pos+1;
        else
            ciclo = false;
            pos = pos-1;
        end
        if pos == nSamplesMax
            ciclo = false;
            pos = pos-1;
        end
    
    end

    S0map(:,:,s)=mask(:,:,s).*mean(volumes(:,:,s,1:pos),4);

    bolus(s) = pos;

end

end

