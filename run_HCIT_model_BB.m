close all;
clear all;
clc;

generate_matrices = 1;   %%%%%% flag to compute control matrices or not
genInfCubeFlag = 0;

%mainProgPath = '/Users/Pikachu/Dropbox/Kasdin_Lab/simulations';
%largeFilePath = '~/Documents/MATLAB/'; % Store G matrix NOT in Dropbox--usually >1GB file size
mainProgPath = '~/Workspace/PHCIL/FPWC/'; % Neil's workstation
largeFilePath = '~/Data/FPWC/';           % Neil's workstation

%Nitr = 30; % Number of control iterations
Nitr = 3;
controller = 'linesearch'; % Don't change for now: stroke min control method
c_range = [-8 0]; % log scale PSF plotting range
plotflag = 1; % flag to plot PSF correction in real time

SPfile = './SPs/SP_AFTA_loqo_hN1k_erNo_c8_r3_4WA13_60deg.fits';
model_params.SP0 = fitsread(SPfile);
BigN = 1000;
model_params.sampling = 3; % pixels per lambda0/D in focal plane
model_params.IP_OWA = 9; % OWA in focal plane, in lambda0/D
lambda0 = 550e-9; % nominal wavelength
model_params.lambda0 = lambda0*1.00; % wavelength used
model_params.Dpup = 48e-3; % meters
model_params.Ddm = model_params.Dpup; % meters
model_params.z_dm1_dm2 = 1; % meters
model_params.fl_M1 = 1.5;  % focal lengths of OAPs
model_params.fl_M2 = model_params.fl_M1/2;
model_params.fl_M3 = 0.774;
model_params.Nact = 16;  % Number of actuators across the DMs
% VtoH1 = 5e-9*ones(Nact); % 5 nm/V in surface change
% VtoH2 = 5e-9*ones(Nact); % 5 nm/V in surface change
model_params.VtoH1 = ones(model_params.Nact); % 
model_params.VtoH2 = ones(model_params.Nact); % 
model_params.Ein = ones(BigN); % Input field at DM1
model_params.abFlag = 1;  % Flag to include aberrations on each optic
PSD_DM1 = (100*1e-9)*fitsread('./errormaps/psd_DM1_5nmRMS_N1000.fits');
PSD_DM2 = (100*1e-9)*fitsread('./errormaps/psd_DM2_5nmRMS_N2000.fits');
PSD_OAP1 = (100*1e-9)*fitsread('./errormaps/psd_OAP1_5nmRMS_N2000.fits');
PSD_OAP2 = (100*1e-9)*fitsread('./errormaps/psd_OAP2_5nmRMS_N2000.fits');
PSD_SP = (100*1e-9)*fitsread('./errormaps/psd_SP_10nmRMS_N1000.fits');
model_params.errmaps = containers.Map({'DM1','DM2','OAP1','OAP2','SP'},...
                                      {PSD_DM1,PSD_DM2,PSD_OAP1,PSD_OAP2,PSD_SP});
fracBW = 0.1; % fractional bandwidth of correction
Nlambda = 3; % number of wavelength samples across correction band
%Nlambda = 1;
if Nlambda > 1
    lambda_vec = linspace(model_params.lambda0 - model_params.lambda0*fracBW/2,...
                          model_params.lambda0 + model_params.lambda0*fracBW/2,...
                          Nlambda) % vector of corrected wavelengths
else
    lambda_vec = [model_params.lambda0];
end
li_ref = round(Nlambda/2);

if(genInfCubeFlag==0)
    cd(largeFilePath);
    load infCubeData
    model_params.infCube = infCube;
    cd(mainProgPath)
else
    model_params.infCube = 1; % place holder value
end
genGmat_False = 0; genGmat_True = 1;
aberr_False = 0; aberr_True = 1;
genInfCube_False = 0; genInfCube_True = 1;

DM_config.num_dms=2; % 1 or 2, number of DMs to use.
DM_config.which_dm=1;  % used if num_dms==1. DM1 is at pupil; DM2 is after pupil
DM_config.DM1V = zeros(model_params.Nact);
DM_config.DM2V = zeros(model_params.Nact);
if (DM_config.num_dms == 2)
    IPsideCor = 'LR'; %'L', 'R', or 'LR'  % which side of image plane to correct
    IPsideScore = 'LR'; %'L', 'R', or 'LR'
elseif (DM_config.num_dms == 1)
    IPsideCor = 'R'; %'L', 'R', or 'LR'
    IPsideScore = 'R'; %'L', 'R', or 'LR'
end

G_mat_fname = sprintf('G_stroke_%dDM_%dpcntBW_at%dnm_Nlambda%02d.mat', DM_config.num_dms,...
                      round(fracBW*100), round(model_params.lambda0*1e9), Nlambda);

% Evaluate aberrated image at central wavelength
I00_unity = 1.;
[E_foc_ab_cent, Lam0D] = HCIT_model(model_params, DM_config, I00_unity, model_params.lambda0, aberr_True, genGmat_False, genInfCubeFlag);
I00_cent = max(max(abs(E_foc_ab_cent).^2));
Im_cent = abs(E_foc_ab_cent).^2/I00_cent;
E_foc_ab_cent = E_foc_ab_cent/sqrt(I00_cent);

E_foc_ab = repmat(E_foc_ab_cent, 1, 1, Nlambda)*0;
Im = E_foc_ab*0;
I00 = zeros(1,Nlambda);
% Evaluate aberrated image at all wavelenth samples in passband
for li=1:Nlambda,
    [E_foc_ab(:,:,li), ~] = HCIT_model(model_params, DM_config, I00_unity, lambda_vec(li), aberr_True, genGmat_False, genInfCube_False);
    I00(li) = max(max(abs(E_foc_ab(:,:,li)).^2));
    Im(:,:,li) = abs(E_foc_ab(:,:,li)).^2/I00(li);
    E_foc_ab(:,:,li) = E_foc_ab(:,:,li)/sqrt(I00(li));
end
Im_bandavg = mean(Im, 3); % wavelength average of normalized images

figure; imagesc(Lam0D,Lam0D,log10(Im_bandavg),[-8 0]); axis square; colorbar;
 title('Uncorrected, band-averaged PSF','Fontsize',24,'Interpreter','LaTeX');
xlabel('x ($\lambda_0$/D)','FontSize',16,'Interpreter','LaTeX'); 
ylabel('y ($\lambda_0$/D)','FontSize',16,'Interpreter','LaTeX');
axis equal; axis tight; axis xy;
% side_lim = IP_OWA;%28.7687;
% xlim([-side_lim side_lim]); ylim([-side_lim side_lim])
set(gca,'FontSize',18,'FontName','Times','FontWeight','Normal')

%figure; imagesc(Lam0D,Lam0D,log10(Im_BB(:,:,1)),[-8 0]); axis square; colorbar;
%figure; imagesc(Lam0D,Lam0D,log10(Im_BB(:,:,3)),[-8 0]); axis square; colorbar;

switch controller
    case 'linesearch'
        mu0 = 1; mu_it = 40; muFac=1.05; target_frac = 0.60;
    case 'fminbnd'
        muLow=1e-6; muHigh=1e6; target_frac = 0.10;
end

%%

Nimg = model_params.sampling*model_params.IP_OWA; % there are (2*Nimg+1) points across the focal plane
flatFlag = 0;           % 1 for trapezoid, 0 for doughnut segment, 3 for rounded claw shape, 4 for rounded doughnut segment
flatFlagScore = 0; 
Mwa_cor = [4 7 0];    % [IWA, OWA, XI_min] for correction area
Mwa_score = [4.0 7 0];      % [IWA, OWA, XI_min] for contrast measurement area
MWAangle=(56)*pi/180;   % angle corrected over, for CorMask
MWAfac = 1;   % amount of CorMask angle to use for ScoreMask 

CorMask = AngleMask_v3(Nimg,Nimg,max(Lam0D),max(Lam0D),Mwa_cor,1*MWAangle,flatFlag,IPsideCor);
ScoreMask = AngleMask_v3(Nimg,Nimg,max(Lam0D),max(Lam0D),Mwa_score,MWAfac*MWAangle,flatFlagScore,IPsideScore);
ScoreMask_Left = AngleMask_v3(Nimg,Nimg,max(Lam0D),max(Lam0D),Mwa_score,MWAfac*MWAangle,flatFlagScore,'L');
ScoreMask_Right = AngleMask_v3(Nimg,Nimg,max(Lam0D),max(Lam0D),Mwa_score,MWAfac*MWAangle,flatFlagScore,'R');

area = sum(sum(ScoreMask));
area_left = sum(sum(ScoreMask_Left));
area_right = sum(sum(ScoreMask_Right));

Maskline = ScoreMask(:).'; %reshape(ScoreMask,1,(2*Nimg+1)^2);   
cor_ele = find(CorMask~=0);  
score_ele = find(ScoreMask~=0);  

% figure; imagesc(ScoreMask); axis square; colormap gray;

% sum(sum(I_foc2.*ScoreMask))/area

%%
% contrast=zeros(1,Nitr+1);
% contrast_left=zeros(1,Nitr+1);
% contrast_right=zeros(1,Nitr+1);
% 
% contrast(1) = sum(sum(Im.*ScoreMask))/area;
% contrast_left(1) = sum(sum(Im.*ScoreMask_Left))/area_left;
% contrast_right(1)=sum(sum(Im.*ScoreMask_Right))/area_right;

contrast = zeros(Nitr+1, Nlambda);
contrast_left = zeros(Nitr+1, Nlambda);
contrast_right = zeros(Nitr+1, Nlambda);
contrast_bandavg = zeros(Nitr+1);
contrast_left_bandavg = zeros(Nitr+1);
contrast_right_bandavg = zeros(Nitr+1);

for li=1:Nlambda,
    contrast(1,li) = sum(sum( Im(:,:,li).*ScoreMask ))/area;
    contrast_left(1,li) = sum(sum( Im(:,:,li).*ScoreMask_Left ))/area_left;
    contrast_right(1,li) = sum(sum( Im(:,:,li).*ScoreMask_Right ))/area_right;
end

contrast_bandavg(1) = mean(contrast(1,:));
contrast_left_bandavg(1) = mean(contrast_left(1,:));
contrast_right_bandavg(1) = mean(contrast_right(1,:));

% CTarget = target_frac*contrast_right(1);
%CTarget = target_frac*contrast(1);
CTarget = target_frac * contrast_bandavg(1);

%Gstar1_cent = zeros(Nact^2,length(cor_ele));
%Gstar2_cent = zeros(Nact^2,length(cor_ele));
Gstar1 = zeros(model_params.Nact^2,length(cor_ele),Nlambda);
Gstar2 = zeros(model_params.Nact^2,length(cor_ele),Nlambda);
% Eim = zeros(2*Nimg+1,2*Nimg+1,Nitr); 
%FieldActual = zeros(length(cor_ele),Nitr);
%Field = zeros(length(cor_ele),Nitr);
FieldActual = zeros(length(cor_ele),Nlambda,Nitr);
Field = zeros(length(cor_ele),Nlambda,Nitr);

% Eim(:,:,1)=E0;
for li=1:Nlambda,
    E_foc_ab_mono = E_foc_ab(:,:,li);
    FieldActual(:,li,1) = E_foc_ab_mono(cor_ele);
end

dDM1V = zeros(model_params.Nact,model_params.Nact);  % delta voltage on DM1
dDM2V = zeros(model_params.Nact,model_params.Nact);  % delta voltage on DM2
DM1Vcor_array = zeros(model_params.Nact,model_params.Nact,Nitr+1);
DM2Vcor_array = zeros(model_params.Nact,model_params.Nact,Nitr+1);
Im_seq = repmat(Im,1,1,1,Nitr+1);
Im_bandavg_seq = repmat(Im_bandavg,1,1,Nitr+1);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Begin the Control Iterations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%function [Eout,Lam0D] = HCIT_model(Ein,I00,SP0,DM1V,DM2V,VtoH1,VtoH2,Ddm,Nact,...
%    sampling,lambda0,lambda,z_dm1_dm2,fl_M1,fl_M2,fl_M3,Dpup,IP_OWA,abFlag,errmaps,...
%    generate_matrices,num_dms,which_dm,genInfCubeFlag,infCube)

for Itr=1:Nitr
    fprintf(['Iteration: ' num2str(Itr) ', Avg contrast @ ref wav. / band-avg''d: %.8e / %.8e\n'], contrast(Itr,li_ref), contrast_bandavg(Itr));
if (Itr==1) && (generate_matrices==1)
fprintf('Creating Influence Matrices ... '); tic

genGmat_DM_config.num_dms = 1;
abFlag_False = 0;
for li=1:Nlambda,
    fprintf('\nChannel %d/%d (%.2f nm):\n', li, Nlambda, lambda_vec(li)*1e9);
    model_params.lambda = lambda_vec(li);
    model_params.I00 = I00(li);
    if (DM_config.num_dms==2) || (DM_config.which_dm==1)        % DM1, compute Jacobian
        genGmat_DM_config.which_dm = 1;
        fprintf(' DM1 ...');
        for q = 1:model_params.Nact^2,
            DMSweep = zeros(model_params.Nact);
            DMSweep(q) = 1;
            genGmat_DM_config.DM1V = DMSweep;
            genGmat_DM_config.DM2V = zeros(model_params.Nact);
            [Etemp1,~] = HCIT_model(model_params, genGmat_DM_config, I00(li), lambda_vec(li), aberr_False, genGmat_True, genInfCube_False);
            Gstar1(q,:,li) = conj(Etemp1(cor_ele)); % Re-order into a vector for the Jacobian matrix
        end   
    end
    if (DM_config.num_dms==2) || (DM_config.which_dm==2)        % DM2 (after pupil), compute Jacobian
        fprintf(' DM2 ...');
        genGmat_DM_config.which_dm = 2;
        for q = 1:model_params.Nact^2,   
            DMSweep = zeros(model_params.Nact);
            DMSweep(q) = 1;
            genGmat_DM_config.DM1V = zeros(model_params.Nact);
            genGmat_DM_config.DM2V = DMSweep;
            [Etemp2,~] = HCIT_model(model_params, genGmat_DM_config, I00(li), lambda_vec(li), aberr_False, genGmat_True, genInfCube_False);
            Gstar2(q,:,li) = conj(Etemp2(cor_ele)); % Re-order into a vector for the Jacobian matrix
        end
    end
end

cd(largeFilePath);
save(G_mat_fname, 'Gstar1', 'Gstar2', 'cor_ele', 'score_ele', 'I00')
cd(mainProgPath)

fprintf(' done. Time: %.3f\n',toc);

elseif (Itr==1) && (generate_matrices==0)
    cd(largeFilePath);
    load(G_mat_fname)
    cd(mainProgPath)
end

if(Itr==1)
    G1 = zeros(length(cor_ele), model_params.Nact*model_params.Nact, Nlambda);
    G2 = zeros(length(cor_ele), model_params.Nact*model_params.Nact, Nlambda);
    if DM_config.num_dms==1
        M = zeros(model_params.Nact*model_params.Nact, model_params.Nact*model_params.Nact, Nlambda);
        EyeM = repmat(eye(model_params.Nact*model_params.Nact,model_params.Nact*model_params.Nact), 1, 1, Nlambda);
    elseif DM_config.num_dms==2
        M = zeros(2*model_params.Nact*model_params.Nact, 2*model_params.Nact*model_params.Nact, Nlambda);
        EyeM = repmat(eye(2*model_params.Nact*model_params.Nact,2*model_params.Nact*model_params.Nact), 1, 1, Nlambda);
    end
    for li=1:Nlambda,
        if DM_config.num_dms==1 && DM_config.which_dm==1
            G1(:,:,li) = ( squeeze(Gstar1(:,:,li)).*repmat(Maskline(cor_ele),model_params.Nact*model_params.Nact,1) )';
            M(:,:,li) = real( squeeze(G1(:,:,li))'*squeeze(G1(:,:,li)) );%/I00;
        elseif DM_config.num_dms==1 && DM_config.which_dm==2
            G2(:,:,li) = ( squeeze(Gstar2(:,:,li)).*repmat(Maskline(cor_ele),model_params.Nact*model_params.Nact,1) )';  
            M(:,:,li) = real( sueeze(G2(:,:,li))'*squeeze(G2(:,:,li)) );%/I00;
        elseif DM_config.num_dms==2
            G1(:,:,li) = ( squeeze(Gstar1(:,:,li)).*repmat(Maskline(cor_ele),model_params.Nact*model_params.Nact,1) )';
            G2(:,:,li) = ( squeeze(Gstar2(:,:,li)).*repmat(Maskline(cor_ele),model_params.Nact*model_params.Nact,1) )';
            MatrixInfluence11 = real( squeeze(G1(:,:,li))'*squeeze(G1(:,:,li)) );%/I00;
            MatrixInfluence12 = real( squeeze(G1(:,:,li))'*squeeze(G2(:,:,li)) );%/I00;
            MatrixInfluence22 = real( squeeze(G2(:,:,li))'*squeeze(G2(:,:,li)) );%/I00;
            M(:,:,li) = [[MatrixInfluence11 MatrixInfluence12];[MatrixInfluence12.' MatrixInfluence22]];
        end
        EyeM(:,:,li) = EyeM(:,:,li) * max(max(diag(squeeze(M(:,:,li)))));
    end
end


% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimation 
Field(:,:,Itr) = FieldActual(:,:,Itr);

% Iinco2D = zeros(2*Nimg+1);
% Iinco2D(cor_ele) = abs(Iinco);
% % Fest = zeros(2*Nimg+1);
% % Fest(cor_ele) = Field(:,Itr);
% % Fact = zeros(2*Nimg+1);
% % Fact(cor_ele) = FieldActual(:,Itr);
% % figure; imagesc(abs(Fest)); axis square; colorbar; title('abs(E est)')
% % figure; imagesc(abs(Fact)); axis square; colorbar; title('abs(E actual)')
% % figure; imagesc(angle(Fest)); axis square; colorbar; title('angle(E est)')
% % figure; imagesc(angle(Fact)); axis square; colorbar; title('angle(E actual)')
% % pause(2);

% % %EimProj is Im{b0} in matrix notation.
% if num_dms==1 && which_dm==1
%     EimProj = imag(hole1*Field(:,Itr));
% elseif num_dms==1 && which_dm==2
%     EimProj = imag(hole2*Field(:,Itr));
% elseif num_dms==2
%     EimProj1 = imag(hole1*Field(:,Itr));
%     EimProj2 = imag(hole2*Field(:,Itr));
%     EimProj = [EimProj1; EimProj2];
% end

% %EimProj is Im{b0} in matrix notation.

if DM_config.num_dms==1
    RealGstarEab = zeros(model_params.Nact*model_params.Nact, Nlambda);
%    I00_mat = repmat(I00, Nact*Nact, 1);
%    I00_cube = repmat( reshape(I00,1,1,Nlambda), Nact*Nact, Nact*Nact, 1 );
elseif DM_config.num_dms==2
    RealGstarEab = zeros(2*model_params.Nact*model_params.Nact, Nlambda);
%    I00_mat = repmat(I00, 2*Nact*Nact, 1);
%    I00_cube = repmat( reshape(I00,1,1,Nlambda), 2*Nact*Nact, 2*Nact*Nact, 1 );
end

for li=1:Nlambda,
    if DM_config.num_dms==1 && DM_config.which_dm==1
        RealGstarEab(:,li) = real(squeeze(G1(:,:,li))'*Field(:,li,Itr));
    elseif DM_config.num_dms==1 && DM_config.which_dm==2
        RealGstarEab(:,li) = real(squeeze(G2(:,:,li))'*Field(:,li,Itr));
    elseif DM_config.num_dms==2
        RealGstarEab1 = real(squeeze(G1(:,:,li))'*Field(:,li,Itr));
        RealGstarEab2 = real(squeeze(G2(:,:,li))'*Field(:,li,Itr));
        RealGstarEab(:,li) = [RealGstarEab1; RealGstarEab2];
    end
end

RealGstarEab_bandavg = mean(RealGstarEab, 2);
M_bandavg = mean(M, 3);
EyeM_bandavg = mean(EyeM, 3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Stroke Minimization Control Algorithm - Takes Estimate and Determines
%Necessary Control to achieve the targeted contrast value
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    fprintf('Control beginning ... '); tic

    switch controller
        case 'fminbnd'%         %  %%%%%%%%%%   Use fminbnd instead of the line search
%         fprintf(['Contrast target: ' num2str(CTarget) '\n'])
%         mu = fminbnd( @(mu) StrokeMinFminBnd_1DMpup(mu, M,EyeM, EimProj, Nact,...
%             infdx_dm2, DM2Vcor, SP0crop, allAb, FTpre, FTpost, I00, ...
%             ScoreMask,area,CTarget) , muLow, muHigh,'iter');
%         fprintf('Optimal Lagrange multiplier = %.4e \n',mu);
%          DMcalc = -(EyeM/mu + M)\EimProj;

        case 'Emin'   %%%%%%%%%%%%%%
%         muLow=1e-2;
%         muHigh=1e6;
%         mu = fminbnd( @(mu) EnergyMinFminBnd_2DM(mu, M, area, EimProj, Nact,Nact,...
%             infdx_dm1,infdx_dm2, DM1Vcor, DM2Vcor, SP0,EDM1nom,num_dms,which_dm, ...
%             FTpre,FTpost, I00, ScoreMask,Ddm,lambda,z1to2,...
%             Npup,Npup_dm,EyeM), muLow, muHigh,'iter');
%         fprintf('Optimal Lagrange multiplier = %.4e \n',mu);
%         DMcalc = -(EyeM/mu + M)\EimProj;
        
        case 'linesearch'    %%%%%%%%%%%%%% begin line search method  
        converge_flag = 0;
        while converge_flag == 0
            fprintf(['Contrast target: ' num2str(CTarget) '\n'])

            mu=mu0;
            k = 1;
            muf = zeros(1,mu_it);
        while (k<4)||((((muf(k-1)-muf(k-2))>0)&&((muf(k-2)-muf(k-3))>0))&&k<=mu_it)||(((muf(k-1)-muf(k-2))<0)&&k<=mu_it);
            %DMcalc = -((  EyeM/mu + M)\RealGstarEab);
            %quad = (DMcalc.'*M*DMcalc)/area;
            %lin = (DMcalc.'*RealGstarEab)/area;
            %Tot = quad + 2*lin + contrast(Itr);
            DMcalc = -(EyeM_bandavg/mu + M_bandavg)\RealGstarEab_bandavg;
            quad = (DMcalc.'*M_bandavg*DMcalc)/area;
            lin = (DMcalc.'*RealGstarEab_bandavg)/area;
            Tot = quad + 2*lin + contrast_bandavg(Itr);
            if Tot < CTarget
                mu = mu/muFac;
            else
                mu = mu*muFac;
            end
            Cfin(k) = Tot;
            muf(k) = mu;
            k = k+1;
        end
        k=k-1;
        
        if muf(end)==mu0*1.05^mu_it
            fprintf('Stroke Minimization Did Not Converge \n')
            CTarget = CTarget*1.15;
        elseif k<1
            CTarget = CTarget*.85;
            fprintf('Too easy for me ... increase target \n')
        else
            converge_flag = 1;
            fprintf(['Stroke Minimization Has Converged in ' num2str(k) ' iterations. mu = ' num2str(mu) '\n'])
        end
        end
        mu0 = mu;  % for next iteration
    %%%%% end of line search
    end
    fprintf(' done. Time: %.3f\n',toc);

    % Simulate the physics
    if (DM_config.num_dms==1 && DM_config.which_dm==1)
        dDM1V = reshape(DMcalc.',model_params.Nact,model_params.Nact);
        DM_config.DM1V = DM_config.DM1V + dDM1V;
    elseif (DM_config.num_dms==1 && DM_config.which_dm==2)
        dDM2V = reshape(DMcalc.',model_params.Nact,model_params.Nact);
        DM_config.DM2V = DM_config.DM2V + dDM2V; % units of phase
    elseif DM_config.num_dms==2
        dDM1V = reshape(DMcalc(1:model_params.Nact^2).',model_params.Nact,model_params.Nact);
        dDM2V = reshape(DMcalc(model_params.Nact^2+1:end).',model_params.Nact,model_params.Nact);
        DM_config.DM1V = DM_config.DM1V + dDM1V; % in radians
        DM_config.DM2V = DM_config.DM2V + dDM2V; % in radians
    end
    
    for li=1:Nlambda,
        [Eout,~] = HCIT_model(model_params, DM_config, I00(li), lambda_vec(li), aberr_True, genGmat_False, genInfCube_False);
        Im(:,:,li) = abs(Eout).^2;
        FieldActual(:,li,Itr+1) = Eout(cor_ele);   % Actual field in esimtation area
    end

    Im_seq(:,:,:,Itr+1) = Im;
    Im_bandavg_seq(:,:,Itr+1) = mean(Im,3);
    
%     if Itr==1 
%         figure; imagesc(Lam0D,Lam0D,log10(I_array(:,:,1)),c_range); ch=colorbar; 
%         title('Aberrated PSF before Correction','FontSize',24,'Interpreter','LaTeX');
%         xlabel('x ($\lambda$/D)','FontSize',16,'Interpreter','LaTeX'); 
%         ylabel('y ($\lambda$/D)','FontSize',16,'Interpreter','LaTeX');
%         xlim([-IP_OWA IP_OWA]); ylim([-IP_OWA IP_OWA])
%         set(gca,'FontSize',18,'FontName','Times','FontWeight','Normal')
%         axis equal; axis tight;
%     end
    
   if(plotflag)     % Real-time image plane plot

    figure(6); imagesc(Lam0D,Lam0D,log10(Im_bandavg_seq(:,:,Itr+1)),c_range); axis square; ch=colorbar; 
    title('Corrected, band-averaged PSF','FontSize',24,'Interpreter','LaTeX');
    xlabel('x ($\lambda$/D)','FontSize',16,'Interpreter','LaTeX'); 
    ylabel('y ($\lambda$/D)','FontSize',16,'Interpreter','LaTeX');
    axis equal; axis tight; axis xy;
%    xlim([-IP_OWA IP_OWA]); ylim([-IP_OWA IP_OWA])
    set(gca,'FontSize',18,'FontName','Times','FontWeight','Normal')
    
%     figure(7); imagesc(Lam0D,Lam0D,log10(Iinco2D),c_range); ch=colorbar; 
%     title('Incoherent Light Estimate','FontSize',24,'Interpreter','LaTeX');
%     xlabel('x ($\lambda$/D)','FontSize',16,'Interpreter','LaTeX'); 
%     ylabel('y ($\lambda$/D)','FontSize',16,'Interpreter','LaTeX');
%     xlim([-IP_OWA IP_OWA]); ylim([-IP_OWA IP_OWA])
%     set(gca,'FontSize',18,'FontName','Times','FontWeight','Normal')
%     axis equal; axis tight;
        pause(2); % Pause to let the plot update
   end

    DM1Vcor_array(:,:,Itr) = DM_config.DM1V;
    DM2Vcor_array(:,:,Itr) = DM_config.DM2V;
    %CorScore_Right = sum(sum(Im.*CorMask_Right))/area;
    for li=1:Nlambda,
        contrast(Itr+1,li) = sum(sum( Im(:,:,li).*ScoreMask ))/area;
        contrast_left(Itr+1,li) = sum(sum( Im(:,:,li).*ScoreMask_Left ))/area_left;
        contrast_right(Itr+1,li) = sum(sum( Im(:,:,li).*ScoreMask_Right ))/area_right;
    end
    contrast_bandavg(Itr+1) = mean(contrast(Itr+1,:));
    contrast_left_bandavg(Itr+1) = mean(contrast_left(Itr+1,:));
    contrast_right_bandavg(Itr+1) = mean(contrast_right(Itr+1,:));

    CTarget = target_frac*contrast(Itr+1); %_right(Itr+1);
    fprintf('Contrast: %.3e Left Contrast: %.3e Right Contrast: %.3e \n \n',...
            contrast(Itr+1),contrast_left(Itr+1),contrast_right(Itr+1));

end

figure(8); semilogy(0:Nitr,contrast(:,li_ref),0:Nitr,contrast_left(:,li_ref),0:Nitr,contrast_right(:,li_ref),...
                    0:Nitr,contrast_bandavg,0:Nitr,contrast_left_bandavg, 0:Nitr,contrast_right_bandavg,...
                    'MarkerSize',19,'LineWidth',1.5);
title('Contrast in Dark Hole','FontSize',24,'Interpreter','LaTeX');
xlabel('Iteration','FontSize',16,'Interpreter','LaTeX'); 
ylabel('Contrast','FontSize',16,'Interpreter','LaTeX');
legend('Overall (lambda_{ref})','Left (lambda_{ref})', 'Right (lambda_{ref})',...
       'Overall (band-avg''d)','Left (band-avg''d)','Right (band-avg''d)',...
       'Location','best');
% xlim([1 Nitr]);
% ylim([ 0.8*contrast_both_des  1.5*contrast(1)])
% ylim([ 0.8*contrast(end)  1.2*contrast(1)])
set(gca,'FontSize',18,'FontName','Times','FontWeight','Normal');
% print -depsc 'contrast_curves.eps'
