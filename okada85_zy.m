function varargout=okada85_zy(varargin)
%OKADA85 Surface deformation due to a finite rectangular source.
%	[uE,uN,uZ,uZE,uZN,uNN,uNE,uEN,uEE] = OKADA85(...
%	   E,N,DEPTH,STRIKE,DIP,LENGTH,WIDTH,RAKE,SLIP,OPEN)经度，纬度，震源深度、走向角、倾角、断层长度、宽度、滑动角、合成滑动量
%	computes displacements, tilts and strains at the surface of an elastic
%	half-space, due to a dislocation defined by RAKE, SLIP, and OPEN on a 
%	rectangular fault defined by orientation STRIKE and DIP, and size LENGTH and
%	WIDTH. The fault centroid is located (0,0,-DEPTH).
%
%	   E,N    : coordinates of observation points in a geographic referential 
%	            (East,North,Up) relative to fault centroid (units are described below)
%	   DEPTH  : depth of the fault centroid (DEPTH > 0)
%	   STRIKE : fault trace direction (0 to 360?relative to North), defined so 
%	            that the fault dips to the right side of the trace
%	   DIP    : angle between the fault and a horizontal plane (0 to 90?
%	   LENGTH : fault length in the STRIKE direction (LENGTH > 0)
%	   WIDTH  : fault width in the DIP direction (WIDTH > 0)
%	   RAKE   : direction the hanging wall moves during rupture, measured relative
%	            to the fault STRIKE (-180 to 180?.
%	   SLIP   : dislocation in RAKE direction (length unit)
%	   OPEN   : dislocation in tensile component (same unit as SLIP)
%
%	returns the following variables (same matrix size as E and N):
%	   uN,uE,uZ        : displacements (unit of SLIP and OPEN)
%	   uZE,uZN         : tilts (in rad * FACTOR)
%	   uNN,uNE,uEN,uEE : horizontal strains POSITIVE = COMPRESSION (unit of FACTOR)
%
%	Length unit consistency: E, N, DEPTH, LENGTH, and WIDTH must have the same 
%	unit (e.g. km) which can be different from that of SLIP and OPEN (e.g. m) but
%	with a possible FACTOR on tilt and strain results (in this case, an 
%	amplification of km/m = 1000). To have FACTOR = 1 (tilt in radians and 
%	correct strain unit), use the same length unit for all aforesaid variables.
%
%	[...] = OKADA85(...,NU) specifies Poisson's ratio NU (default is 0.25 for
%	an isotropic medium).
%
%	Formulas and notations from Okada [1985] solution excepted for strain 
%	convention (here positive strain means compression), and for the fault 
%	parameters after Aki & Richards [1980], e.g.:
%	      DIP=90, RAKE=0   : left lateral (senestral) strike slip
%	      DIP=90, RAKE=180 : right lateral (dextral) strike slip
%	      DIP=70, RAKE=90  : reverse fault
%	      DIP=70, RAKE=-90 : normal fault
%
%	It is also possible to produce partial outputs, with following syntax:
%	   [uE,uN,uZ] = OKADA85(...) for displacements only;
%	   [uE,uN,uZ,uZE,uZN] = OKADA85(...) for displacements and tilts;
%	   [uE,uN,uZ,uNN,uNE,uEN,uEE] = OKADA85(...) for displacements and strains;
%	   [uZE,uZN] = OKADA85(...) for tilts only;
%	   [uZE,uZN,uNN,uNE,uEN,uEE] = OKADA85(...) for tilts and strains;
%	   [uNN,uNE,uEN,uEE] = OKADA85(...) for strains only.
%
%	Note that vertical strain components can be obtained with following equations:
%	   uNZ = -uZN;
%	   uEZ = -uZE;
%	   uZZ = -(uEE + uNN)*NU/(1-NU);
%
%	[...] = OKADA85(...,'plot') or OKADA85(...) without output argument 
%	produces a 3-D figure with fault geometry and dislocation at scale (if
%	all of the fault parameters are scalar).
%
%	Equations are all vectorized excepted for argument DIP which must be
%	a scalar (beacause of a singularity in Okada's equations); all other
%	arguments can be scalar or matrix of the same size.
%
%	Example:
%
%	   [E,N] = meshgrid(linspace(-10,10,50));
%	   [uE,uN,uZ] = okada85(E,N,2,30,70,5,3,-45,1,1,'plot');
%	   figure, surf(E,N,uN)用surf函数画出来的图 叫surf 三维着色表面图、三维表面图、表面图
%
%	considers a 5x3 fault at depth 2, N30?strike, 70?dip, and unit dislocation
%	in all directions (reverse, senestral and open). Displacements are computed
%	on a regular grid from -10 to 10, and North displacements are plotted as a
%	surface.

%


nu = 0.25;	% isotropic Poisson's ratio

% Assigns input arguments
e = varargin{1};
n = varargin{2};
depth = varargin{3};
strike = varargin{4}*pi/180;	% converting STRIKE in radian
dip = varargin{5}*pi/180;	% converting DIP in radian ('delta' in Okada's equations)
L = varargin{6};
W = varargin{7};
rake = varargin{8}*pi/180;	% converting RAKE in radian
slip = varargin{9};
U3 = varargin{10};


% Defines dislocation in the fault plane system
U1 = cos(rake).*slip;
U2 = sin(rake).*slip;

% Converts fault coordinates (E,N,DEPTH) relative to centroid
% into Okada's reference system (X,Y,D)
d = depth + sin(dip).*W/2;	% d is fault's top edge
ec = e + cos(strike).*cos(dip).*W/2;
nc = n - sin(strike).*cos(dip).*W/2;
x = cos(strike).*nc + sin(strike).*ec + L/2;
y = sin(strike).*nc - cos(strike).*ec + cos(dip).*W;
%将经纬度转化为x,y坐标
% Variable substitution (independent from xi and eta)
p = y.*cos(dip) + d.*sin(dip);
q = y.*sin(dip) - d.*cos(dip);

% Displacements
if any(nargout==[3, 5, 7, 9])
	ux = -U1/(2*pi) .* chinnery(@ux_ss,x,p,L,W,q,dip,nu) ... % strike-slip
		- U2/(2*pi) .* chinnery(@ux_ds,x,p,L,W,q,dip,nu) ... % dip-slip
		+ U3/(2*pi) .* chinnery(@ux_tf,x,p,L,W,q,dip,nu); ... % tensile fault

	uy = -U1/(2*pi) .* chinnery(@uy_ss,x,p,L,W,q,dip,nu) ... % strike-slip
		- U2/(2*pi) .* chinnery(@uy_ds,x,p,L,W,q,dip,nu) ... % dip-slip
		+ U3/(2*pi) .* chinnery(@uy_tf,x,p,L,W,q,dip,nu); ... % tensile fault

	uz = -U1/(2*pi) .* chinnery(@uz_ss,x,p,L,W,q,dip,nu) ... % strike-slip
		- U2/(2*pi) .* chinnery(@uz_ds,x,p,L,W,q,dip,nu) ... % dip-slip
		+ U3/(2*pi) .* chinnery(@uz_tf,x,p,L,W,q,dip,nu); ... % tensile fault

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Notes for I... and K... subfunctions:
%
%	1. original formulas use Lame's parameters as mu/(mu+lambda) which
%	   depends only on the Poisson's ratio = 1 - 2*nu
%	2. tests for cos(dip) == 0 are made with "cos(dip) > eps" 
%	   because cos(90*pi/180) is not zero but = 6.1232e-17 (!)
%	   NOTE: don't use cosd and sind because of incompatibility
%	   with Matlab v6 and earlier...


% =================================================================
% Chinnery's notation [equation (24) p. 1143]

% -----------------------------------------------------------------
function u=chinnery(f,x,p,L,W,q,dip,nu)
u = feval(f,x,p,q,dip,nu) ...
	- feval(f,x,p-W,q,dip,nu) ...
	- feval(f,x-L,p,q,dip,nu) ...
	+ feval(f,x-L,p-W,q,dip,nu);


% =================================================================
% Displacement subfunctions

% strike-slip displacement subfunctions [equation (25) p. 1144]

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

% dip-slip displacement subfunctions [equation (26) p. 1144]
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

% tensile fault displacement subfunctions [equation (27) p. 1144]
% -----------------------------------------------------------------
function u=ux_tf(xi,eta,q,dip,nu)
R = sqrt(xi.^2 + eta.^2 + q.^2);
u = q.^2 ./(R.*(R + eta)) ...
	- I3(eta,q,dip,nu,R).*sin(dip).^2;

% -----------------------------------------------------------------
function u=uy_tf(xi,eta,q,dip,nu)
R = sqrt(xi.^2 + eta.^2 + q.^2);
u = -(eta.*sin(dip) - q.*cos(dip)).*q./(R.*(R + xi)) ...
	- sin(dip).*xi.*q./(R.*(R + eta)) ...
	- I1(xi,eta,q,dip,nu,R).*sin(dip).^2;
k = find(q~=0);
u(k) = u(k) + sin(dip).*atan(xi(k).*eta(k)./(q(k).*R(k)));

% -----------------------------------------------------------------
function u=uz_tf(xi,eta,q,dip,nu)
R = sqrt(xi.^2 + eta.^2 + q.^2);
db = eta.*sin(dip) - q.*cos(dip);
u = (eta.*cos(dip) + q.*sin(dip)).*q./(R.*(R + xi)) ...
	+ cos(dip).*xi.*q./(R.*(R + eta)) ...
	- I5(xi,eta,q,dip,nu,R,db).*sin(dip).^2;
k = find(q~=0);
u(k) = u(k) - cos(dip).*atan(xi(k).*eta(k)./(q(k).*R(k)));


% I... displacement subfunctions [equations (28) (29) p. 1144-1145]
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


