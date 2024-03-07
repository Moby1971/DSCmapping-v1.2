function [res_svd] = DSC_mri_SVD(conc,aif,mask,options,app)

% Calculates parametric maps of Cerebral Blood Flow (CBF) for a subject
% and uses the SINGULAR VALUE DECOMPOSITION method with truncation for deconvolution
%
% Input parameters - conc (4D Matrix), contains the DSC concentration trends of all voxels.
%                  - aif, concentration trend in the site chosen as arterial
%                  - mask (3D Matrix), contains the matrix used to mask the brain volume to be analyzed
% Options is the struct that contains the method options, the significant ones are:
%
% options.deconv.svd.threshold - truncation threshold as a percentage of
%                                the maximum eigenvalue, in Ostergaard and
%                                Calamante et al. it is fixed at 20%
%
% options.deconv.SVD.residual - if set to 1, it also outputs the 4D matrix
%                               of residuals, otherwise, they are not calculated

aifVett=zeros(options.nT,1);
aifVett(1)=aif(1);
aifVett(options.nT)=aif(options.nT);

for k=2:(options.nT-1)
    aifVett(k)=(aif(k-1)+4*aif(k)+aif(k+1))/6;
end

G = toeplitz(aifVett,[aifVett(1) zeros(1,options.nT-1)]);

[U,S,V] = svd(G);

eigenV=diag(S);
threshold=options.deconv.SVD.threshold*max(eigenV);   
newEigen=zeros(size(eigenV));
for k=1:length(eigenV)
    if eigenV(k)>=threshold
        newEigen(k)=1/eigenV(k);
    end
end

Ginv=V*diag(newEigen)*(U');

res_svd.map=zeros(options.nR,options.nC,options.nS);
if options.deconv.SVD.residual
	res_svd.residual=zeros(options.nR,options.nC,options.nS,options.nT);
end


% MAIN LOOP
for s=1:options.nS

    for r=1:options.nR

        tic;

        for c=1:options.nC

            if mask(r,c,s)

                vettConc=reshape(conc(r,c,s,:),options.nT,1);
                vettRes=(1/options.tr)*Ginv*vettConc;

                res_svd.map(r,c,s)=max(abs(vettRes));
                if options.deconv.SVD.residual
                    res_svd.residual(r,c,s,:)=vettRes;
                end

            end

        end

        app.elapsedTime = app.elapsedTime + toc;
        app.counter = app.counter + options.step_cbf;
        app.d.updateProgressBar(app);

    end

end


end