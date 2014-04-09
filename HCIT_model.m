% Flags:
%  -Aberrations on/off
% INPUTS:
% -Ein
% -
% -DM commands
% -DM gains
% -SP file name
% -sampling in focal plane
% -lambda0  (central wavelength)
% -wavelength
% -z_dm1_dm2
% -fl_M1, fl_M2, fl_M3
% -Dpup
% -Nact
% -Ddm  (later)
% -IP_OWA

% OUTPUTS:
% -Eout
% -LamD  (focal plane grid coordinates)



% SPfile = './SPs/SP_AFTA_loqo_hN1k_erNo_c8_r3_4WA13_60deg.fits';
% sampling = 3; % pixels per lambda0/D
% IP_OWA = 15; % lambda0/D
% lambda0 = 550e-9;
% lambda = lambda0*1.10;
% Dpup = 48e-3; % meters
% Ddm = Dpup; % meters
% z_dm1_dm2 = 1; % meters
% fl_M1 = 1.5;
% fl_M2 = fl_M1/2;
% fl_M3 = 0.774;
% Nact = 24;
% % DM1V = 
% VtoH1 = 5e-9; % 5 nm/V in surface change
% VtoH2 = 5e-9; % 5 nm/V in surface change

% Npup=2*1000/2;
% % % % % % % function P=ABB2D(N,N1Waves,N2Waves,Nt,RMS,Kolmogorov)
% phaseAb1 = ABB2D(2*Npup,2,50,10,0.05,1);
% % phaseAb2 = ABB2D(2*Npup,1,50,3,0.09,1);
% % figure; imagesc(phaseAb1); axis square; axis xy; colorbar;
% std(phaseAb1(:))
% fitswrite(phaseAb1,'./errormaps/psd_OAP2_5nmRMS_N2000.fits');

function [Eout,Lam0D] = HCIT_model(params, DM_config, I00, lambda, abFlag, generate_matrices, genInfCubeFlag)

Ein = params.Ein;
SP0 = params.SP0;
VtoH1 = params.VtoH1;
VtoH2 = params.VtoH2;
Ddm = params.Ddm;
Dpup = params.Dpup;
Nact = params.Nact;
sampling = params.sampling;
lambda0 = params.lambda0;
z_dm1_dm2 = params.z_dm1_dm2;
fl_M1 = params.fl_M1;
fl_M2 = params.fl_M2;
fl_M3 = params.fl_M3;
IP_OWA = params.IP_OWA;
errmaps = params.errmaps;
infCube = params.infCube;
num_dms = DM_config.num_dms;
which_dm = DM_config.which_dm;
DM1V = DM_config.DM1V;
DM2V = DM_config.DM2V;

%function [Eout,Lam0D] = HCIT_model(Ein,I00,SP0,DM1V,DM2V,VtoH1,VtoH2,Ddm,Nact,...
%    sampling,lambda0,lambda,z_dm1_dm2,fl_M1,fl_M2,fl_M3,Dpup,IP_OWA,abFlag,errmaps,...
%    generate_matrices,num_dms,which_dm,genInfCubeFlag,infCube)

if(size(Ein,1)~=size(SP0,1))
    disp('Error: Input field and SP dimensions do not match.');
    return
end

Npup = length(SP0)/2;
% SP0pad = padarray(SP0,[2*Npup,2*Npup,0]); % Make 3x as wide

if(abFlag)
% %    X  % Extra factor of 2 because in reflection
     EDM1error = exp(1i*errmaps('DM1')*(2*pi/lambda));
     EDM2error = exp(1i*errmaps('DM2')*(2*pi/lambda));
     EOAP1error = exp(1i*errmaps('OAP1')*(2*pi/lambda));
     EOAP2error = exp(1i*errmaps('OAP2')*(2*pi/lambda));
     ESPerror = exp(1i*errmaps('SP')*(2*pi/lambda));
else
     EDM1error = ones(size(errmaps('DM1')));
     EDM2error = ones(size(errmaps('DM2')));
     EOAP1error = ones(size(errmaps('OAP1')));
     EOAP2error = ones(size(errmaps('OAP2')));
     ESPerror = ones(size(errmaps('SP')));
end


samplingFlag=0;
pupFac = fl_M2/fl_M1;
z_dm2_M1 = fl_M1-z_dm1_dm2;
Dsp = fl_M1/fl_M2;


dx_dm = Dpup/2/Npup;
xs_dm = (-Npup:Npup-1)'*dx_dm + dx_dm/2;
Nimg = sampling*IP_OWA;

% DM1U = zeros(Nact);
% DM2U = zeros(Nact);

% DM1U = DM1V*VtoH1*2*(2*pi/lambda);  % DM1 command in radians
% DM2U = DM2V*VtoH2*2*(2*pi/lambda);  % DM2 command in radians
% DM1H = DM1V*VtoH1*2;  % DM1 command in meters, x2 because its a mirror
% DM2H = DM2V*VtoH2*2;  % DM2 command in meters



% Influence Functions based on a Gaussian
InfFuncSigma = (Ddm/Nact*1.4);

infdx = zeros(Nact,length(xs_dm));
for q = 1:Nact
    x_cent = q-Nact/2-1/2;
    infdx(q,:) = exp(-4*log(2)*((xs_dm-Ddm*x_cent/(Nact-1)).^2)./(InfFuncSigma)^2);
end

% padarray(u0,[P*(pad_fac-1)/2 Q*(pad_fac-1)/2,0])
% figure; imagesc(SP0pad); colormap gray; axis square;

temp = fitsread('./pupils/pupil_D1Kpix.fits'); sideBuffer = 12;
AFTApupil = temp(sideBuffer+1:1000+sideBuffer,sideBuffer+1:1000+sideBuffer);
DM1stop = flipud(AFTApupil);
% DM1stop = ones(2*Npup);


%  DM1stop=zeros(2*Npup);
% for ii=1:2*Npup
%     for jj=1:2*Npup
%         if ( ((jj-Npup-0.5)^2 + (ii-Npup-0.5)^2) <= (Npup+1)^2) 
%             DM1stop(ii,jj)=1;
%         end
%     end
% end
% figure; imagesc(DM1stop); colormap gray; axis square; axis xy;

EinStop = Ein.*DM1stop;
company = 'custom';
pitch_other = Dpup/Nact;

NptsBuffer = round(2*Npup/10);
NptsZeroPad = Npup-NptsBuffer;
EinStopPad = padarray(EinStop,[Npup,Npup,0],0);
EDM1errorPad = padarray(EDM1error,[Npup,Npup,0],0);

% Don't forget factor of 2 from being in reflection


if(generate_matrices)  % When making control matrices
    if(num_dms==1 && which_dm==1)
        DM1surfU = makeDMsurf_v1(Nact,2*Npup,NptsBuffer,DM1V,VtoH1,company,pitch_other,genInfCubeFlag,infCube); % in radians
        DM1surfUpad = padarray(DM1surfU,[NptsZeroPad,NptsZeroPad,0],0);
        E_dm1 = (1i*DM1surfUpad).*EinStopPad;
        E_dm2 = prop_PTP(E_dm1,2*Dpup,lambda,z_dm1_dm2,samplingFlag);
    elseif(num_dms==1 && which_dm==2)
        E_dm1 =  padarray( EinStop,[Npup,Npup,0]);
        DM2surfU = makeDMsurf_v1(Nact,2*Npup,NptsBuffer,DM2V,VtoH2,company,pitch_other,genInfCubeFlag,infCube);
        DM2surfUpad = padarray(DM2surfU,[NptsZeroPad,NptsZeroPad,0],0);
        E_dm2 = (1i*DM2surfUpad).*prop_PTP(E_dm1,2*Dpup,lambda,z_dm1_dm2,samplingFlag);  
    end
else
    DM1surfU = makeDMsurf_v1(Nact,2*Npup,NptsBuffer,DM1V,VtoH1,company,pitch_other,genInfCubeFlag,infCube);
    DM1surfUpad = padarray(DM1surfU,[NptsZeroPad,NptsZeroPad,0],0);
    E_dm1 = exp(1i*DM1surfUpad).*EinStopPad.*EDM1errorPad;
    
    DM2surfU = makeDMsurf_v1(Nact,2*Npup,NptsBuffer,DM2V,VtoH2,company,pitch_other,genInfCubeFlag,infCube);
    DM2surfUpad = padarray(DM2surfU,[NptsZeroPad,NptsZeroPad,0],0);
    E_dm2 = exp(1i*DM2surfUpad).*EDM2error.*prop_PTP(E_dm1,2*Dpup,lambda,z_dm1_dm2,samplingFlag);
end



% if(generate_matrices)  % When making control matrices
%     if(num_dms==1 && which_dm==1)
%         E_dm1 =  1i*padarray( Ein.*(DM1stop.*(infdx.'*DM1U*infdx)),[Npup,Npup,0]); % Make twice as big as pupil    
%         E_dm2 = prop_PTP(E_dm1,2*Dpup,lambda,z_dm1_dm2,samplingFlag);
%     elseif(num_dms==1 && which_dm==2)
%         E_dm1 =  padarray( Ein.*DM1stop,[Npup,Npup,0]);
%         
%         DM2part = 1i*padarray(infdx.'*DM2U*infdx,[Npup,Npup,0],0); % Only zero pad
%         E_dm2 = DM2part.*prop_PTP(E_dm1,2*Dpup,lambda,z_dm1_dm2,samplingFlag);  
%     end
% else
%     E_dm1 =  padarray( EDM1error.*Ein.*(DM1stop.*(exp(1i*infdx.'*DM1U*infdx))),[Npup,Npup,0]); % Make twice as big as pupil
% 
%     DM2part = padarray( exp(1i*infdx.'*DM2U*infdx),[pad1,pad1,0],1); % Make mirror larger
%     DM2part = padarray( DM2part,[pad2,pad2,0],0);
%     E_dm2 = EDM2error.*DM2part.*prop_PTP(E_dm1,2*Dpup,lambda,z_dm1_dm2,samplingFlag);
% end

E_M1 = EOAP1error.*prop_PTP(E_dm2,2*Dpup,lambda,z_dm2_M1,samplingFlag);
% figure; imagesc(abs(E_dm1));axis square; colorbar;

E_M2 = EOAP2error.*prop_SWS(E_M1,2*Dpup,lambda,fl_M1,fl_M2);
% figure; imagesc(abs(E_M2)); axis square; colorbar; title('E M2');

E_pup2 =  prop_PTP(E_M2,pupFac*2*Dpup,lambda,fl_M2,samplingFlag);

E_pup2_crop_SP = ESPerror.*SP0.*E_pup2(Npup+2:3*Npup+1,Npup+2:3*Npup+1);

% figure; imagesc(abs(E_pup2_crop_SP)); axis xy; colormap gray; axis square;

% % % E_pup2_crop_SP = abs(E_pup2_crop_SP); % Goal of phase retrieval
% E_exitPup_phase = angle(ESPerror.*E_pup2(Npup+1:3*Npup,Npup+1:3*Npup));
% % E_exitPup_phase(E_exitPup_phase<0) = E_exitPup_phase(E_exitPup_phase<0)+pi;
% figure; imagesc(angle(E_exitPup_phase)); axis square; axis xy; colorbar;
%     title('Phase at exit pupil');

%Image Plane Parameters
dx = Dsp/2/Npup;
xs = (-Npup:Npup-1)'*dx + dx/2;
dxi = IP_OWA*lambda0*fl_M3/Dsp/Nimg;  % pixel scaling fixed in terms of lambda0
% DZrow = 2*Nimg+1;
% DZcol= 2*Nimg+1;
xis = (-Nimg:Nimg)'*dxi;
Lam0D = xis/lambda0*Dsp/fl_M3;

%-----% FTs For Each Wavelength %-----% FT(X) = FTpre*X*FTpost
x_xi = xs*xis';
FTpre = dx*(exp(-2*pi*1i*x_xi/(lambda*fl_M3))).';
FTpost = dx/(1i*lambda*fl_M3)*exp(-2*pi*1i*x_xi /(lambda*fl_M3));

E_foc2 = FTpre*E_pup2_crop_SP*FTpost;
% I_foc2 = abs(E_foc2/max(max(abs(E_foc2)))).^2;
% 
% figure; imagesc(log10(I_foc2),[-8 0]); axis square; colorbar;
Eout = E_foc2/sqrt(I00);

