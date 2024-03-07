function [res_csvd]=DSC_mri_cSVD(conc,aif,mask,options,app)
% ultima modifica: Denis Peruzzo 07/06/2010

%Funzione del pacchetto DSC_mri - DSC_mri_cbf
%Autore: Castellaro Marco - Universit� di Padova - DEI
%
%Calcola le mappe parametriche di Cerebral Blood Flow (CBF) per un soggetto
%e per la deconvoluzione utilizza il metodo della SINGULAR VALUE
%DECOMPOSITION versione BLOCK-CIRCULANT con troncamento
%
%Parametri in ingresso - conc (Matrice 4D), contiene gli andamenti delle
%                        concentrazioni DSC di tutti i voxel.
%                      - aif, andamento delle concentrazioni nel sito
%                        scelto come arteriale
%                      - mask (Matrice 3D), contiene la matrice utilizzata
%                        per mascherare il volume cerebrale da analizzare
%Options � la sruct che contiene le opzioni del metodo, quelli
%significativi sono:
%
%options.deconv.svd.threshold - soglia di troncamento percentuale riferita
%                               all'autovalore massimo, in Wu et al.
%                               Calamante � fissato al 10%
%
%options.deconv.SVD.residual - se a 1 produce in uscita anche la matrice 4D
%                              dei residui, altrimenti non vengono
%                              calcolati


nTpad=2*options.nT;
columnG=zeros(nTpad,1);
columnG(1)=aif(1);
columnG(options.nT)=(aif(options.nT-1)+4*aif(options.nT))/6;
columnG(options.nT+1)=aif(options.nT)/6;

for k=2:(options.nT-1)
    columnG(k)=(aif(k-1)+4*aif(k)+aif(k+1))/6;
end

rowG=zeros(1,nTpad);
rowG(1)=columnG(1);

for k=2:nTpad
    rowG(k)=columnG(nTpad+2-k);
end

G=toeplitz(columnG,rowG);

[U,S,V]=svd(G);

eigenV=diag(S);
threshold=options.deconv.cSVD.threshold*max(eigenV);
newEigen=zeros(size(eigenV));

for k=1:length(eigenV)
    if eigenV(k)>=threshold
        newEigen(k)=1/eigenV(k);
    end
end

Ginv=V*diag(newEigen)*(U');

res_csvd.map=zeros(options.nR,options.nC,options.nS);
if options.deconv.SVD.residual
    res_csvd.residual=zeros(options.nR,options.nC,options.nS,nTpad);
end


% MAIN LOOP
for s=1:options.nS

    for r=1:options.nR

        tic;

        for c=1:options.nC

            if mask(r,c,s)

                vettConc=zeros(nTpad,1);
                vettConc(1:options.nT)=reshape(conc(r,c,s,:),options.nT,1);
                vettRes=(1/options.tr)*Ginv*vettConc;

                res_csvd.map(r,c,s)=max(abs(vettRes));

                if options.deconv.cSVD.residual
                    res_csvd.residual(r,c,s,:)=vettRes;
                end

            end

        end

        app.elapsedTime = app.elapsedTime + toc;
        app.counter = app.counter + options.step_cbf;
        app.d.updateProgressBar(app);

    end

end



end