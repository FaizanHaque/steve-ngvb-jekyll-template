%Radius of the earth
a = 6371; % in km

%dzero returns a sparse matrix (?)
d0 = @(x) spdiags(x(:),0,length(x(:)),length(x(:)));




% setup the multiplying factors for the operators
zt = [0;5;10;15;20;25;30;35;40;45;50;55;60;65;70;75;80;85;90;95;100;125;...
    150;175;200;225;250;275;300;325;350;375;400;425;450;475;500;550;...
    600;650;700;750;800;850;900;950;1000;1050;1100;1150;1200;1250;1300;...
    1350;1400;1450;1500;1550;1600;1650;1700;1750;1800;1850;1900;1950;...
    2000;2100;2200;2300;2400;2500;2600;2700;2800;2900;3000;3100;3200;...
    3300;3400;3500;3600;3700;3800;3900;4000;4100;4200;4300;4400;4500;...
    4600;4700;4800;4900;5000;5100;5200;5300;5400;5500];
%    zt = zt(1:nz);
nz = length(zt);
lambda_beg = 0.5;
lambda_end = 359.5;

phi_beg = -89.5;
phi_end = 89.5;
lambdat = [lambda_beg:1:lambda_end]';
phit = [phi_beg:1:phi_end]';
nlambda = length(lambdat);
nphi = length(phit);

elambda = ones(nlambda,1);
ephi = ones(nphi,1);
ez = ones(nz,1);

dphit = phit([2:end,1]) - phit; 
dphi = zeros( nphi,1,1); 
dphi(:,1,1)= dphit; 
dphi = dphi(:,elambda,ez);


dlambdat = lambdat([2:end,1]) - lambdat; 
dlambda = zeros(1, nlambda,1); 
dlambda(1,:,1)= dlambdat; 
dlambda = dlambda(ephi,:,ez);


dzt = zt([2:end,1]) - zt; 
dz = zeros(1,1, nz); 
dz(1,1,:) =dzt; 
dz = dz(ephi,elambda,:);


Difactor = d0(1./ (a.*dphi));
Djfactor =  d0(1./(a.*(cos(dphit).*dlambda))); %is this the correct dphit?
Dkfactor = d0(1./dz);


%Initialize nmax
%ny = 5;
%nx =4 ;
%nz = 102;

%Ii = zeros(ny,1,1);
%Ii = 1:ny;
%Ii (:,ones (nx,1),ones(nz,1));

%[Ij,Ii,Ik] = meshgrid(1:nx,1:ny,1:nz


%II = zeros (ny,nx,nz);
%II(:) = 1:ny*nx*nz;

II = zeros (nphi,nlambda,nz);
II(:) = 1:nphi*nlambda*nz;


%The lower case are the indices
iE = II (:, [ 2:end,1],:);
iW = II (:, [ end, 1:end-1],:);
iN = II ( [ 2:end,1],:,:);
iS = II ( [ end, 1:end-1],:,:);
iD = II (:,:, [ 2:end,1]);
iU = II (:,:, [ end, 1:end-1]);

%Need an sparse identity matrix (only stores the nonzero values
%I = speye(ny*nx*nz);
I = speye(nphi*nlambda*nz);
%
IE= I(iE(:),:);
IW= I(iW(:),:);
IN= I(iN(:),:);
IS= I(iS(:),:);
IU= I(iU(:),:);
ID= I(iD(:),:);


%
%a= 6371 e3;
%Latitude
%phi = -89.5:1:89.5;
%Longitude

%Derivatives with respect to the indices
FDj = IE - I;
FDi = IN - I;
FDk = ID - I;


%Derivatives with respect to the indices
BDj = abs(IW - I);
BDi = abs(IS - I);
BDk = abs(IU - I);


%Derivatives

BddPhi = BDi.*Difactor;
BddLambda = BDj.*Difactor;
BddZ = BDk.*Dkfactor;

FddPhi = FDi.*Difactor;
FddLambda = FDj.*Difactor;
FddZ = BDk.*Dkfactor;

%Central Derivatives:

CDPhi = 0.5*(BddPhi+FddPhi);
CDLambda = 0.5*(BddPhi+FddLambda);
CDZ = 0.5*(BddZ + FddZ);


%Secondorder Derivatives
D2Phi = FddPhi*BddPhi;
D2Lambda = FddPhi*BddLambda;
D2Z = FddZ*BddZ;

D2PhiBF = BddPhi*FddPhi;
D2LambdaBF = BddPhi*FddLambda;
D2Z = BddZ*FddZ;


%nabla-operators
Fgradh = [FddPhi, FddLambda]';
Bgradh = [BddPhi, BddLambda]';
Cgradh = [CDPhi, CDLambda]';



FcurlH = [-FddPhi, FddLambda];
BcurlH = [-BddPhi, BddLambda];
CcurlH = [-CDPhi, CDLambda];



LaplacianH = D2Phi + D2Lambda;
Laplacian3D = D2Phi + D2Lambda + D2Z;







%   Fdz(1,1,:) = [2.5;zt(2:end)-zt(1:end-1)]




%[M,X,Y,Z] readwoa13(16);


