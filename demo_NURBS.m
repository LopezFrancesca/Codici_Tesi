
function demo_NURBS(a,sigma,XYW,XYWR)

%--------------------------------------------------------------------------
% OBJECT
%--------------------------------------------------------------------------
% Numerical experiment in "Filtered hyperinterpolation on NURBS-shaped planar domains".
% Region: NURBS.
%--------------------------------------------------------------------------
% Usage:
% >> demo_NURBS;
%--------------------------------------------------------------------------
% Note:
% The routine uses 'binornd' that requires Statistics and Machine Learning
% Toolbox.
%--------------------------------------------------------------------------
% Dates:
% Written on June 8, 2023: A. Sommariva and F. Lopez.
%--------------------------------------------------------------------------
%

clear; clf;
warning off

%--------------------------------------------------------------------------
% Degrees of precision in numerical experiments: can be a vector.
%--------------------------------------------------------------------------
LV=10;     % Hyperinterpolant tot degree.
NV=2*LV;   % Degree of the rule used in hyperinterpolation.
NR=18;     % Degree of precision of the reference rule (estimate L2 error).


%--------------------------------------------------------------------------
% Noise.
%--------------------------------------------------------------------------

if nargin <3, XYW=[]; end
if nargin <4, XYWR=[]; end
noise=0;

if noise
    if nargin <1,a=1.5; end     % defining impulse noise (in experiment 2)
    if nargin <2,sigma=0.5; end % defining gaussian noise (in experiment 2)
    fprintf('\n a=%.2f \n sigma=%.2f \n',a, sigma)
else
    a=0; sigma=0; % no noise.
end


% * Function to approximate:
% 1. degree L poly., 2. degree floor(L/2)-1 poly. 3. test functions
% (see line 265 approx).
funct_example=3 ;

% No table or stats.
display_stats=1;

% Domain
domain_example=20 ;

%--------------------------------------------------------------------------
% Special settings.
%--------------------------------------------------------------------------

% Plot domain and nodes: do_plot=1 (yes), do_plot=0 (no).
do_plot=1;

% Number of tests for reconstructing functions of a type on this domain.
ntests=30;




% ........................ Main code below ................................


% ......................... Define domain .................................

[geometry_NURBS,domain_str]=define_NURBS_domain(domain_example);


% ........ Numerical approximation, varying the degree in "nV" ............

AEinfMV=[]; AE2MV=[]; % vectors used for statistics
JzMV=[]; HzMV=[];

for k=1:length(NV)
    N=NV(k); % Quadrature points.
    L=LV(k); % Hyperinterpolant degree.


    % Define quadrature rule for hyperinterpolation at L, with  ADE=N.
    if isempty(XYW),
        XYW = cub_NURBS(N,geometry_NURBS);
    end
    X=XYW(:,1); Y=XYW(:,2); W=XYW(:,3);

    % Test points
    if isempty(XYWR),
        XYWR= cub_NURBS(NR,geometry_NURBS);
    end
    XR=XYWR(:,1); YR=XYWR(:,2); WR=XYWR(:,3);

    % Compute bounding box
    structure_RS=geometry_NURBS;
    [xyw0, ~, ~, ~, ~,bbox] = cub_NURBS(0,structure_RS);

    plot_NURBS(geometry_NURBS); hold on; plot(X,Y,'ro');


    % rectangle is [bbox(1), bbox(2)] x [bbox(3), bbox(4)];
    xLimit=[bbox(1), bbox(2)]; yLimit=[bbox(3), bbox(4)];
    dbox=[xLimit' yLimit'];

    % Compute orthonormal basis matrix at nodes.
    jvec=1:(L+1)*(L+2)/2;
    [U,~,Q,R,~,degs] = dORTHVAND(L,[X Y],W,jvec,[],dbox);

    % .. testing AE_L2err hyperinterpolation error for each "f" at "deg" ..
    poly_coeffs=[];

    for j=1:ntests

        % ... define function to approximate ...
        g=define_function(funct_example,L);

        % ... evaluate function to approximate ...
        gXY=feval(g,X,Y);

        % ... add noise (if present) ...

        % a) add impulse noise
        pert_impulse=0;
        if a > 0
            pert_impulse=a*(1-2*rand(length(gXY),1))*binornd(1,0.5);
            while norm(pert_impulse) == 0
                pert_impulse=a*(1-2*rand(length(gXY),1))*binornd(1,0.5);
            end
        end

        % b) add gaussian noise
        pert_gauss=0;
        if sigma > 0
            var=sigma^2;
            pert_gauss=sqrt(var)*randn(size(gXY));
            while norm(pert_gauss) == 0
                pert_gauss=sqrt(var)*randn(size(gXY));
            end
        end

        % add gaussian + impulse noise
        pert=pert_impulse+pert_gauss;

        % perturbed values
        gXY_pert=gXY+pert;

        % ... determine polynomial hyperinterpolant ...
        % compute hyperinterpolant coefficients
        coeff0=Q'*(sqrt(W).*gXY_pert);

        % ... test hyperinterpolant with or withour filters ...

        for ktest=[1 2]
            switch ktest
                case 1
                    hypermode='filtered';
                    coeff=hyperfilter(hypermode,coeff0,degs);
                case 2
                    hypermode='hyperinterpolation';
                    coeff=coeff0;
            end

            % evaluate polynomial at reference points.
            gXYR=feval(g,XR,YR);
            VR=chebvand(L,[XR YR],dbox);
            pXYR = (VR(:,jvec)/R)*coeff;

            % errors
            AEinfV(ktest,j)=norm(gXYR-pXYR,inf); % absolute error (inf norm)
            AE2V(ktest,j)=sqrt(WR'*((gXYR-pXYR).^2)); % absolute error (2 norm)
            beta0V(ktest,j)=sum(abs(coeff) > 0);
            % evaluate J(coeff) and H(coeff), that are error relevant
            % parameters, as observed in Thm 5.1.

        end


    end

    % averages of the errors (vectors 5 x 1)
    AEinfM=mean(AEinfV,2);
    AE2M=mean(AE2V,2);

    if display_stats

        AEinf=AEinfM([1 2],:); AE2=AE2M([1 2],:); 

        fprintf('\n       ........ table at degree: %2.0f ........ \n \n ',N);
        HypType=categorical({'filtered';  'hyperint.'});
        T = table(HypType,AEinf,AE2); 
        disp(T)





    end

    AEinfMV=[AEinfMV AEinfM]; AE2MV=[AE2MV AE2M]; 

end

Wg_norm2=(norm(sqrt(W).*gXY,2))^2;








function g=define_function(funct_example,L)

% function to test

switch funct_example

    case 1 % test exactness hyperinterpolation
        nexp=L;
        c0=0.4474; 
        %c0=rand(1);
        c1=0.8275;
        %c1=rand(1); 
        c2=0.9318;
        %c2=rand(1);
        g=@(x,y) (c0+c1*x+c2*y).^nexp;

    case 2 % test exactness filt. hyperinterpolation
        nexp=max(floor(L/2),0);
        c0=0.6623; 
        %c0=rand(1); 
        c1=0.4254;
        %c1=rand(1);
        c2=0.4195;
        %c2=rand(1);
        g=@(x,y) (c0+c1*x+c2*y).^nexp;

    case 3 % function of that type

        funct_example_sub=2;

        switch funct_example_sub
            case 1
                g=@(x,y) exp(-(x.^2+y.^2));
            case 2
                g=@(x,y) sin(-(x.^2+y.^2));
            case 3
                g=@(x,y) 1+0*x+0*y;
            case 4
                g=@(x,y) sqrt((x-0.5).^2+(y-0.5).^2);
            case 5
                g=@(x,y) (0.2*x+0.5*y).^19;
            case 6
                g=@(x,y) exp((x-0.5).^2+(y-0.5).^2);
            case 7
                g=@(x,y) exp(-100*((x-0.5).^2+(y-0.5).^2));
            case 8
                g=@(x,y) cos(30*(x+y));
            case 9
                g=@(x,y) cos(5*(x+y));
            case 10
                g=@(x,y) exp((x-0.5).^1+(y-0.5).^1);
            case 11
                g=@(x,y) exp((x-0.5).^3+(y-0.5).^3);
            case 12
                g=@(x,y) (0.2*x+0.5*y).^15;
            case 13
                g=@(x,y) 1./(x.^2+y.^2);
            case 14
                % g=@(x,y) (x+y).^ade;
                g=@(x,y) (1+x+0.5*y).^ade;
            case 16
                x0=0.5; y0=0.5;
                g=@(x,y) exp(-((x-x0).^2+(y-y0).^2));
            case 17
                x0=0.5; y0=0.5;
                g=@(x,y) ((x-x0).^2 + (y-y0).^2).^(3/2);
        end


end









