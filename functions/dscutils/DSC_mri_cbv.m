function [cbv]=DSC_mri_cbv(conc,aif,mask,options,app)
%Funzione del pacchetto DSC_mri - DSC_mri_cbv
%Autore: Castellaro Marco - Università di Padova - DEI
%
%Calcola le mappe parametriche di Cerebral Blood Volume per un soggetto
%
%Parametri in ingresso: conc (Matrice 4D) che contiene gli andamenti delle
%concentrazioni DSC di tutti i voxel.
%Options è la sruct che contiene le opzioni del metodo, quelli
%significativi sono:
%
%options.time - Rappresenta il vettore dei tempi dell'esame DSC (ogni
%               campione rappresenta l'acquisizione di un intero volume cerebrale
%
%options.par.kh - Parametro che rappresenta la dipendenza dall'ematocrito
%                 del Cerebral Blood Volume (CBV), per default è settato ad
%                 uno, in questo caso si ottengono stime RELATIVE del
%                 parametro CBV.
%
%options.par.rho - Parametro che rappresenta la dipendenza dalla densità
%                  del sangue del Cerebral Blood Volume (CBV), per default
%                  è settato ad uno, in questo caso si ottengono stime
%                  RELATIVE del parametro CBV.
%
%Parametri in uscita:
%cbv - (matrice 3D) che contiene al suo interno la mappa parametrica
%      calcolata



cbv=zeros(options.nR,options.nC,options.nS);


for s=1:options.nS

    tic;

    cbv(:,:,s)=(options.par.kh/options.par.rho)*mask.data(:,:,s).* ...
        (trapz(options.time,conc(:,:,s,:),4)./trapz(options.time,aif));

    app.elapsedTime = app.elapsedTime + toc;
    app.counter = app.counter + options.step_cvb;
    app.d.updateProgressBar(app);

end

end
