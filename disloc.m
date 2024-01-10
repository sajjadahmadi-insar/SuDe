function varargout=disloc(varargin)
% Rectangular Dislocation
if nargin < 9
	error('Not enough input arguments.')
end

if nargin > 10
	error('Too many input arguments.')
end

if any(~cellfun(@isnumeric,varargin(1:9)))
	error('Input arguments E,N,DEPTH,STRIKE,DIP,LENGTH,WIDTH,RAKE,SLIP must be numeric.')
end

if ~isscalar(varargin{5})
	error('DIP argument must be scalar.')
end

% Default values for optional input arguments
nu = 0.25;	% isotropic Poisson's ratio

% Assigns input arguments
e = varargin{1};
n = varargin{2};
depth = varargin{3};
strike = varargin{4}*pi/180;
dip = varargin{5}*pi/180;
L = varargin{6};
W = varargin{7};
rake = varargin{8}*pi/180;
slip = varargin{9};

% Defines dislocation in the fault plane system
U1 = cos(rake).*slip;
U2 = sin(rake).*slip;

% Converts fault coordinates (E,N,DEPTH) relative to centroid
% into Okada's reference system (X,Y,D)
d = depth + sin(dip).*W/2;
ec = e + cos(strike).*cos(dip).*W/2;
nc = n - sin(strike).*cos(dip).*W/2;
x = cos(strike).*nc + sin(strike).*ec + L/2;
y = sin(strike).*nc - cos(strike).*ec + cos(dip).*W;

% Variable substitution (independent from xi and eta)
p = y.*cos(dip) + d.*sin(dip);
q = y.*sin(dip) - d.*cos(dip);

% Displacements
if any(nargout==[3, 5, 7, 9])
	ux = -U1/(2*pi) .* chinnery(@ux_ss,x,p,L,W,q,dip,nu) ... % strike-slip
		- U2/(2*pi) .* chinnery(@ux_ds,x,p,L,W,q,dip,nu); ... % dip-slip

	uy = -U1/(2*pi) .* chinnery(@uy_ss,x,p,L,W,q,dip,nu) ... % strike-slip
		- U2/(2*pi) .* chinnery(@uy_ds,x,p,L,W,q,dip,nu); ... % dip-slip

	uz = -U1/(2*pi) .* chinnery(@uz_ss,x,p,L,W,q,dip,nu) ... % strike-slip
		- U2/(2*pi) .* chinnery(@uz_ds,x,p,L,W,q,dip,nu); ... % dip-slip

	% Rotation from Okada's axes to geographic
	ue = sin(strike).*ux - cos(strike).*uy;
	un = cos(strike).*ux + sin(strike).*uy;
end
 
% Assigns output arguments
switch nargout
	case 3
		varargout = {ue, un, uz};
	otherwise
		disp('Unvalid number of output arguments.')
end

% =================================================================
% Chinnery's notation

% -----------------------------------------------------------------
function u=chinnery(f,x,p,L,W,q,dip,nu)
u = feval(f,x,p,q,dip,nu) ...
	- feval(f,x,p-W,q,dip,nu) ...
	- feval(f,x-L,p,q,dip,nu) ...
	+ feval(f,x-L,p-W,q,dip,nu);

% =================================================================
% Displacement subfunctions

% strike-slip displacement subfunctions

% -----------------------------------------------------------------
function u=ux_ss(xi,eta,q,dip,nu)
R = sqrt(xi.^2 + eta.^2 + q.^2);
u = xi.*q./(R.*(R + eta)) ...
	+ I1(xi,eta,q,dip,nu,R).*sin(dip);
k = find(q~=0);
u(k) = u(k) + atan(xi(k).*eta(k)./(q(k).*R(k)));

% -----------------------------------------------------------------
function u=uy_ss(xi,eta,q,dip,nu)
R = sqrt(xi.^2 + eta.^2 + q.^2);
u = (eta.*cos(dip) + q.*sin(dip)).*q./(R.*(R + eta)) ...
	+ q.*cos(dip)./(R + eta) ...
	+ I2(eta,q,dip,nu,R).*sin(dip);

% -----------------------------------------------------------------
function u=uz_ss(xi,eta,q,dip,nu)
R = sqrt(xi.^2 + eta.^2 + q.^2);
db = eta.*sin(dip) - q.*cos(dip);
u = (eta.*sin(dip) - q.*cos(dip)).*q./(R.*(R + eta)) ...
	+ q.*sin(dip)./(R + eta) ...
	+ I4(db,eta,q,dip,nu,R).*sin(dip);

% dip-slip displacement subfunctions
% -----------------------------------------------------------------
function u=ux_ds(xi,eta,q,dip,nu)
R = sqrt(xi.^2 + eta.^2 + q.^2);
u = q./R ...
	- I3(eta,q,dip,nu,R).*sin(dip).*cos(dip);

% -----------------------------------------------------------------
function u=uy_ds(xi,eta,q,dip,nu)
R = sqrt(xi.^2 + eta.^2 + q.^2);
u = (eta.*cos(dip) + q.*sin(dip)).*q./(R.*(R + xi)) ...
	- I1(xi,eta,q,dip,nu,R).*sin(dip).*cos(dip);
k = find(q~=0);
u(k) = u(k) + cos(dip).*atan(xi(k).*eta(k)./(q(k).*R(k)));

% -----------------------------------------------------------------
function u=uz_ds(xi,eta,q,dip,nu)
R = sqrt(xi.^2 + eta.^2 + q.^2);
db = eta.*sin(dip) - q.*cos(dip);
u = db.*q./(R.*(R + xi)) ...
	- I5(xi,eta,q,dip,nu,R,db).*sin(dip).*cos(dip);
k = find(q~=0);
u(k) = u(k) + sin(dip).*atan(xi(k).*eta(k)./(q(k).*R(k)));

% I... displacement subfunctions
% -----------------------------------------------------------------
function I=I1(xi,eta,q,dip,nu,R)
db = eta.*sin(dip) - q.*cos(dip);
if cos(dip) > eps
	I = (1 - 2*nu) * (-xi./(cos(dip).*(R + db))) ...
		- sin(dip)./cos(dip).*I5(xi,eta,q,dip,nu,R,db);
else
	I = -(1 - 2*nu)/2 * xi.*q./(R + db).^2;
end

% -----------------------------------------------------------------
function I=I2(eta,q,dip,nu,R)
I = (1 - 2*nu) * (-log(R + eta)) - I3(eta,q,dip,nu,R);

% -----------------------------------------------------------------
function I=I3(eta,q,dip,nu,R)
yb = eta.*cos(dip) + q.*sin(dip);
db = eta.*sin(dip) - q.*cos(dip);
if cos(dip) > eps
	I = (1 - 2*nu) * (yb./(cos(dip)*(R + db)) - log(R + eta)) ...
		+ sin(dip)./cos(dip) * I4(db,eta,q,dip,nu,R);
else
	I = (1 - 2*nu)/2 * (eta./(R + db) + yb.*q./(R + db).^2 - log(R + eta));
end

% -----------------------------------------------------------------
function I=I4(db,eta,q,dip,nu,R)
if cos(dip) > eps
	I = (1 - 2*nu) * 1./cos(dip) * (log(R + db) - sin(dip).*log(R + eta));
else
	I = -(1 - 2*nu) * q./(R + db);
end

% -----------------------------------------------------------------
function I=I5(xi,eta,q,dip,nu,R,db)
X = sqrt(xi.^2 + q.^2);
if cos(dip) > eps
	I = (1 - 2*nu) * 2./cos(dip) ...
		.* atan((eta.*(X + q.*cos(dip)) + X.*(R + X).*sin(dip)) ...
			./(xi.*(R + X).*cos(dip)));
	I(xi==0) = 0;
else
	I = -(1 - 2*nu) * xi.*sin(dip)./(R + db);
end
