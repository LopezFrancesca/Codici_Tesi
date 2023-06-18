

function [geometry_NURBS,domain_str]=define_NURBS_domain(example)

switch example

    case 1 % M-shaped domain

        domain_str='S(3): M shaped domain';
        order=3;
        P=[-1 1; -1 0; -1 -1; -0.6 -1; -0.6 0; -0.6 0.6; 0 0; 0.6 0.6; ...
            0.6 0; 0.6 -1; 1 -1; 1 0; 1 1; 0.6 1; 0 0.4; -0.6 1; -1 1];
        kk=size(P,1)-order;
        knots_mid=linspace(0.1,0.9,kk);
        knots=[0 0 0 knots_mid 1 1 1]; % 12-3
        w=1000*ones(1,size(P,1));

        geometry_NURBS=makeNURBSarc('free','P',...
            P,'knots',knots,'weights',w,'order',order);

    case 2

        domain_str='S(4): domain defined by a disk, an ellipse, a segment';
        geometry_NURBS(1)=makeNURBSarc('disk_arc',...
            'center',[0 0],'angles',[0 pi/2],'radius',1);

        % add arc of an ellipse
        geometry_NURBS(2)=makeNURBSarc('elliptical_arc',...
            'center',[0 0],'angles',[pi/2 3*pi/2],'ell_axis',[0.5 1],...
            'tilt_angle',0);

        % add a final segment
        Pend=lastpointNURBSPL(geometry_NURBS(2));
        Pinit=firstpointNURBSPL(geometry_NURBS(1));
        geometry_NURBS(3)=makeNURBSarc('segment','vertices',[Pend; Pinit]);

        % join pieces
        geometry_NURBS=joinNURBSPLarcs(geometry_NURBS);

    case 3 % lune-like
        domain_str='S(5): lune-like';
        P=[1 0; 0.45 0.37; 1 1; 0 1; -0.8 0.9; -1 0; 0.3 0.5; 0 -1; ...
            1 -1;  1 0];
        knots=[0 0 0 .15 .25 .25 .5 .5 .75 .75 1 1 1];
        c=1/sqrt(2); w=[1 c c 1 c 1 c 1 c 1];
        order=3;

        geometry_NURBS=makeNURBSarc('free','P',...
            P,'knots',knots,'weights',w,'order',order);


    case 4 % square/disk difference
        domain_str='domain 4';

        % add arc of a disk
        geometry_NURBS(1)=makeNURBSarc('disk_arc',...
            'center',[0 0],'angles',[pi/2 0],'radius',0.25);

        % compute first point of the piecewise NURBS domain
        Pinit=firstpointNURBSPL(geometry_NURBS(1));

        % compute last point of the so made NURBS
        Pend=lastpointNURBSPL(geometry_NURBS(1));


        % add arc of an ellipse
        v=[0.5 0;0.5 0.5;0 0.5];
        vertices=[Pend; v; Pinit];
        geometry_NURBS(2)=makeNURBSarc('polygonal_arc',...
            'vertices',vertices);

        % join piecewise NURBS
        geometry_NURBS=joinNURBSPLarcs(geometry_NURBS);





    case 5 % polygon/disk union
        domain_str='domain 5';

        % add arc of a disk
        geometry_NURBS(1)=makeNURBSarc('disk_arc',...
            'center',[0.25 0.25],'angles',[pi pi+pi/2],'radius',0.25);

        % compute first point of the piecewise NURBS domain
        Pinit=firstpointNURBSPL(geometry_NURBS(1));

        % compute last point of the so made NURBS
        Pend=lastpointNURBSPL(geometry_NURBS(1));

        % add arc of an ellipse
        v=[0.4 0.05; 0.5 0.25; 0.45 0.45; 0.3 0.5;0.1 0.45]; % polyg. vert.
        vertices=[Pend; v; Pinit];
        geometry_NURBS(2)=makeNURBSarc('polygonal_arc',...
            'vertices',vertices);

        % join piecewise NURBS
        geometry_NURBS=joinNURBSPLarcs(geometry_NURBS);





    case 6 % polygon/disk union
        domain_str='domain 6';

        % add arc of a disk
        geometry_NURBS(1)=makeNURBSarc('disk_arc',...
            'center',[0 0],'angles',[pi/2 0],'radius',0.25);

        % compute first point of the piecewise NURBS domain
        Pinit=firstpointNURBSPL(geometry_NURBS(1));

        % compute last point of the so made NURBS
        Pend=lastpointNURBSPL(geometry_NURBS(1));

        % add arc of an ellipse
        v=[0.25 0.2; 0.2 0.25];  % polyg. vert.
        vertices=[Pend; v; Pinit];
        geometry_NURBS(2)=makeNURBSarc('polygonal_arc',...
            'vertices',vertices);

        % join piecewise NURBS
        geometry_NURBS=joinNURBSPLarcs(geometry_NURBS);



    case 7 % polygon/disk union
        domain_str='domain 7';

        % add arc of a disk
        geometry_NURBS(1)=makeNURBSarc('disk_arc',...
            'center',[0 0],'angles',[pi/2 0],'radius',0.25);

        % compute first point of the piecewise NURBS domain
        Pinit=firstpointNURBSPL(geometry_NURBS(1));

        % compute last point of the so made NURBS
        Pend=lastpointNURBSPL(geometry_NURBS(1));

        % add arc of an ellipse
        v=[0.4 0.05; 0.5 0.25; 0.45 0.45; 0.3 0.5;0.1 0.45]; % polyg. vert.
        vertices=[Pend; v; Pinit];
        geometry_NURBS(2)=makeNURBSarc('polygonal_arc',...
            'vertices',vertices);

        % join piecewise NURBS
        geometry_NURBS=joinNURBSPLarcs(geometry_NURBS);


    case 8 % polygon/disk union
        domain_str='domain 8: square';

        v=[0 0; 1 0; 1 1; 0 1; 0 0];
        geometry_NURBS(1)=makeNURBSarc('polygonal_arc',...
            'vertices',v);


    case 9 % polygon/disk union
        domain_str='domain 9: union of rectangles';

        % add arc of an ellipse
        v=[0 0; 1 0; 1 1; 2 1; 2 2; 1 2; 1 1.5; 0 1.5; 0 1.25; -0.5 1.25; ...
            -0.5 0; 0 0];
        geometry_NURBS(1)=makeNURBSarc('polygonal_arc',...
            'vertices',v);


    case 10
        domain_str='domain 10: square';
        P=[0 0; 1 0; 1 1; 0 1; 0 0];
        geometry_NURBS=makeNURBSarc('polygonal_arc','vertices',P);

    case 11 % circle
        domain_str='domain 11: circle';
        geometry_NURBS=makeNURBSarc('disk_arc',...
            'center',[0 0],'angles',[0 2*pi],'radius',1);

    case 12 % lune-like
        domain_str='domain 12: lune-like';
        P=[1 0; 0.45 0.37; 1 1; 0 1; -0.8 0.9; -1 0; 0.3 0.5; 0 -1; ...
            1 -1;  1 0];
        order=3;
        L=size(P,1); tt=linspace(0,1,L-order+2);
        ttinit=zeros(1,order); ttmid=tt(2:end-1);  ttfin=ones(1,order);
        knots=[ttinit ttmid ttfin];
        c=1/sqrt(2); w=[1 c c 1 c 1 c 1 c 1];

        geometry_NURBS=makeNURBSarc('free',...
            'P',P,'knots',knots,'weights',w,'order',order);


    case 13 % tau-like-symbol
        domain_str='domain 13: tau-like-symbol';
        P=[1 0; 0.45 0.37; 1 1; 0 1; -0.8 0.9; -1 0; 0.3 0.5;  ...
            0.2 -0.45; -0.4 -0.4; 0 -1; 1 -1;  1 0];
        order=3;
        L=size(P,1); tt=linspace(0,1,L-order+2);
        ttinit=zeros(1,order); ttmid=tt(2:end-1);  ttfin=ones(1,order);
        knots=[ttinit ttmid ttfin];
        c=1/sqrt(2); w=[1 c c 1 c 1 c 1 c 1 c 1];

        geometry_NURBS=makeNURBSarc('free','P',...
            P,'knots',knots,'weights',w,'order',order);


    case 14 % cubic-domain
        domain_str='domain 14: cubic-domain: leaf';
        P=[1 0; 1 1; 0 1; -1 1; -1 0; -1 -1; 0 -1; 1 -1;  1 0];
        order=4;
        L=size(P,1); tt=linspace(0,1,L-order+2);
        ttinit=zeros(1,order); ttmid=tt(2:end-1);  ttfin=ones(1,order);
        knots=[ttinit ttmid ttfin];
        c=1/sqrt(2); w=[1 c 1 c 1 c 1 c 1];

        geometry_NURBS=makeNURBSarc('free','P',...
            P,'knots',knots,'weights',w,'order',order);



    case 15 % cubic-domain: nut
        domain_str='domain 15: cubic-domain: slanted nut';
        P=[1 0; 1 1; 0 1; -1 0.6; -1 0; -1 -1; 0 -1; 1 -1;  1 0];
        order=4;
        L=size(P,1); tt=linspace(0,1,L-order+2);
        ttinit=zeros(1,order); ttmid=tt(2:end-1);  ttfin=ones(1,order);
        knots=[ttinit ttmid ttfin];
        c=1/sqrt(2); w=[1 c 1 c 1 c 1 c 1];

        geometry_NURBS=makeNURBSarc('free','P',...
            P,'knots',knots,'weights',w,'order',order);


    case 16 % cubic-domain: slanted boomerang
        domain_str='domain 16: cubic-domain: golf club';
        P=[1 0; 1 1; 0 1;  -1 0;  -1 -1; -0.4 0.6; 0 0.5; 1 -1;  1 0];
        order=4;
        L=size(P,1); tt=linspace(0,1,L-order+2);
        ttinit=zeros(1,order); ttmid=tt(2:end-1);  ttfin=ones(1,order);
        knots=[ttinit ttmid ttfin];
        c=1/sqrt(2); w=[1 c/2 1 c/4 1 c/4 1 c 1];

        geometry_NURBS=makeNURBSarc('free','P',...
            P,'knots',knots,'weights',w,'order',order);

    case 17 % cubic-domain: nut
        domain_str='domain 17: cubic-domain: dino head';
        P=[1 0; 1 1; 0 1; -0.4 0.6; -1 0; -1 -1; 0 -1; 0 0.5; 1 -1;  1 0];
        order=4;
        L=size(P,1); tt=linspace(0,1,L-order+2);
        ttinit=zeros(1,order); ttmid=tt(2:end-1);  ttfin=ones(1,order);
        knots=[ttinit ttmid ttfin];
        c=1/sqrt(2); w=[1 c/2 1 c/4 1 c/4 1 1 c 1];


        geometry_NURBS=makeNURBSarc('free','P',...
            P,'knots',knots,'weights',w,'order',order);


    case 18 % cubic-domain: slanted boomerang
        domain_str='domain 18: cubic-domain: tadpole';
        P=[1 0; 1 1; 0 1;  -1 0;  -1 -1; -0.4 0.6; 0 -0.5; 1 -1;  1 0];
        order=4;
        L=size(P,1); tt=linspace(0,1,L-order+2);
        ttinit=zeros(1,order); ttmid=tt(2:end-1);  ttfin=ones(1,order);
        knots=[ttinit ttmid ttfin];
        c=1/sqrt(2); w=[1 1 1 1 1 1 1 1 1];

        geometry_NURBS=makeNURBSarc('free','P',...
            P,'knots',knots,'weights',w,'order',order);

    case 19 % cubic-domain: nut
        domain_str='domain 19: cubic-domain: nose';
        P=[1 0; 1 1; 0 1; -0.4 0.6; -1 0; -1 -1; 0 0.5; 1 -1;  1 0];
        order=4;
        L=size(P,1); tt=linspace(0,1,L-order+2);
        ttinit=zeros(1,order); ttmid=tt(2:end-1);  ttfin=ones(1,order);
        knots=[ttinit ttmid ttfin];
        c=1/sqrt(2); w=[1 5 5 5 5 5 5 5 1];


        geometry_NURBS=makeNURBSarc('free','P',...
            P,'knots',knots,'weights',w,'order',order);

    case 20 % rough-ball
        domain_str='domain 20: rough ball';
        order=3;
        sidesL=20;
        t=linspace(0,2*pi,sidesL); t=t';
        P=[cos(t) sin(t)]; P=[P(1:end-1,:); P(1,:)];
        knots_mid=linspace(0.1,0.9,sidesL-order);
        knots=[zeros(1,order) knots_mid ones(1,order)];
        M=size(P,1);
        % w=rand(1,M);
        w=[ 8.258169774895474e-01
            5.383424352600571e-01
            9.961347166268855e-01
            7.817552875318368e-02
            4.426782697754463e-01
            1.066527701805844e-01
            9.618980808550537e-01
            4.634224134067444e-03
            7.749104647115024e-01
            8.173032206534330e-01
            8.686947053635097e-01
            8.443584551091032e-02
            3.997826490988965e-01
            2.598704028506542e-01
            8.000684802243075e-01
            4.314138274635446e-01
            9.106475944295229e-01
            1.818470283028525e-01
            2.638029165219901e-01
            1.455389803847170e-01]';

        geometry_NURBS=makeNURBSarc('free','P',...
            P,'knots',knots,'weights',w,'order',order);

    case 21 % L-shaped domain

        domain_str='domain 21: L shaped domain';
        P=[-1 1; -1 -1; 1 -1; 1 -0.6; -0.6 -0.6; -0.6 1; -1 1]; % 7 points
        order=3;
        L=size(P,1); tt=linspace(0,1,L-order+2);
        ttinit=zeros(1,order); ttmid=tt(2:end-1);  ttfin=ones(1,order);
        knots=[ttinit ttmid ttfin];
        w=[0.2 5 10 10 100 5 0.2];

        geometry_NURBS=makeNURBSarc('free','P',...
            P,'knots',knots,'weights',w,'order',order);

    case 22 % variable order
        domain_str='domain 22, defined by a disk, an ellipse and a segment';

        % add arc of a disk
        geometry_NURBS(1)=makeNURBSarc('disk_arc',...
            'center',[0 0],'angles',[0 pi/2],'radius',1);

        % compute first point of the piecewise NURBS domain
        Pinit=firstpointNURBSPL(geometry_NURBS(1));

        % compute last point of the so made NURBS
        Pend=lastpointNURBSPL(geometry_NURBS(1));

        % add arc of an ellipse
        geometry_NURBS(2)=makeNURBSarc('elliptical_arc',...
            'center',Pend-[1 0],'angles',[0 pi+pi/4],...
            'ell_axis',[1 2],'tilt_angle',0);

        % compute last point of the so made NURBS
        Pend=lastpointNURBSPL(geometry_NURBS(2));

        % "close" the boundary with a segment
        geometry_NURBS(3)=makeNURBSarc('segment','extrema',...
            [Pend; Pinit]);

        % join piecewise NURBS
        geometry_NURBS=joinNURBSPLarcs(geometry_NURBS);



    case 23 % variable order

        domain_str1='domain 23, defined by a disk, an ellipse, a segment';
        domain_str2=' and a free NURBS (variable order)';
        domain_str=strcat(domain_str1,domain_str2);

        % add arc of a disk
        geometry_NURBS(1)=makeNURBSarc('disk_arc',...
            'center',[0 0],'angles',[0 pi/2],'radius',1);

        % compute first point of the piecewise NURBS domain
        Pinit=firstpointNURBSPL(geometry_NURBS(1));

        % compute last point of the so made NURBS
        Pend=lastpointNURBSPL(geometry_NURBS(end));

        % add arc of an ellipse
        geometry_NURBS(2)=makeNURBSarc('elliptical_arc',...
            'center',Pend-[1 0],'angles',[0 pi+pi/4],'ell_axis',[1 2],...
            'tilt_angle',0);

        % compute last point of the so made NURBS
        Pend=lastpointNURBSPL(geometry_NURBS(end));

        % add segment
        geometry_NURBS(3)=makeNURBSarc('segment','extrema',[Pend; 0 Pend(2)]);

        % compute first point of the piecewise NURBS domain
        Pinit=firstpointNURBSPL(geometry_NURBS);

        % "close" the boundary with a "free" NURBS.
        geometry_NURBS(4)=makeNURBSarc('free',...
            'P',[0 Pend(2); -1.9 0.3; -1.8 0.5; Pinit],...
            'knots',[0 0 0 0 1 1 1 1],'weights',[1 1 2 1],'order',4);

        % join piecewise NURBS
        geometry_NURBS=joinNURBSPLarcs(geometry_NURBS);

    case 24 % variable order

        domain_str='domain 24, defined by a disk, an ellipse, a segment';
        geometry_NURBS(1)=makeNURBSarc('disk_arc',...
            'center',[0 0],'angles',[0 pi/2],'radius',1);

        % add arc of an ellipse
        geometry_NURBS(2)=makeNURBSarc('elliptical_arc',...
            'center',[0 0],'angles',[pi/2 3*pi/2],'ell_axis',[0.5 1],...
            'tilt_angle',0);

        % add a final segment
        Pend=lastpointNURBSPL(geometry_NURBS(2));
        Pinit=firstpointNURBSPL(geometry_NURBS(1));
        geometry_NURBS(3)=makeNURBSarc('segment','vertices',[Pend; Pinit]);

        % join pieces
        geometry_NURBS=joinNURBSPLarcs(geometry_NURBS);


    case 25 % cubic-domain: nut
        domain_str='domain 25: 15-inverse cubic-domain: slanted nut';
        P=[1 0; 1 1; 0 1; -1 0.6; -1 0; -1 -1; 0 -1; 1 -1;  1 0];
        P=P(:,[2 1]);
        order=4;
        L=size(P,1); tt=linspace(0,1,L-order+2);
        ttinit=zeros(1,order); ttmid=tt(2:end-1);  ttfin=ones(1,order);
        knots=[ttinit ttmid ttfin];
        c=1/sqrt(2); w=[1 c 1 c 1 c 1 c 1];

        geometry_NURBS=makeNURBSarc('free','P',...
            P,'knots',knots,'weights',w,'order',order);

    case 26 % circular segment
        domain_str='circular segment';

        C=[0 0]; alpha=pi/4; beta=pi; r=1;

        % add arc of a disk
        geometry_NURBS(1)=makeNURBSarc('disk_arc',...
            'center',C,'angles',[alpha beta],'radius',r);

        % compute first point of the piecewise NURBS domain
        Pinit=firstpointNURBSPL(geometry_NURBS(1));

        % compute last point of the so made NURBS
        Pend=lastpointNURBSPL(geometry_NURBS(1));

        % "close" the boundary with a segment
        geometry_NURBS(2)=makeNURBSarc('segment','extrema',...
            [Pend; Pinit]);

        % join pieces
        geometry_NURBS=joinNURBSPLarcs(geometry_NURBS);




end









