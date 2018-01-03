
function [ int , err , nr_points ] = quadcc ( f , a , b , tol )
%QUADCC  evaluates an integral using adaptive quadrature. The
%   algorithm uses Clenshaw-Curtis quadrature rules of increasing
%   degree in each interval and bisects the interval if either the
%   function does not appear to be smooth or a rule of maximum
%   degree has been reached. The error estimate is computed from the
%   L2-norm of the difference between two successive interpolations
%   of the integrand over the nodes of the respective quadrature rules.
%
%   INT = QUADCC ( F , A , B , TOL ) approximates the integral
%   of F in the interval [A,B] up to the relative tolerance TOL. The
%   integrand F should accept a vector argument and return a vector
%   result containing the integrand evaluated at each element of the
%   argument.
%
%   [INT,ERR,NPOINTS] = QUADCC ( F , A , B , TOL ) returns
%   ERR, an estimate of the absolute integration error as well as
%   NPOINTS, the number of function values for which the integrand
%   was evaluated. The value of ERR may be larger than the requested
%   tolerance, indicating that the integration may have failed.
%
%   QUADCC halts with a warning if the integral is or appears
%   to be divergent.
%
%   Reference: "Increasing the Reliability of Adaptive Quadrature
%       Using Explicit Interpolants", P. Gonnet, ACM Transactions on
%       Mathematical Software, 37 (3), art. no. 26, 2010.

% Copyright (C) 2012 Pedro Gonnet

    % declare persistent variables
    persistent n xi_1 xi_2 xi_3 xi_4 b_1_def b_2_def b_3_def b_4_def ...
        V_1 V_2 V_3 V_4 V_1_inv V_2_inv V_3_inv V_4_inv xx Vcond T_left T_right w U
    
    % have the persistent variables been declared already?
    if ~exist('U') || isempty(U)
    
        % the nodes and newton polynomials
        n = [4,8,16,32];
        xi_1 = -cos([0:n(1)]/n(1)*pi)';
        xi_2 = -cos([0:n(2)]/n(2)*pi)';
        xi_3 = -cos([0:n(3)]/n(3)*pi)';
        xi_4 = -cos([0:n(4)]/n(4)*pi)';
        b_1_def = [0., .233284737407921723637836578544e-1, 0., -.831479419283098085685277496071e-1, 0.,0.0541462136776153483932540272848 ]';
        b_2_def = [0., .883654308363339862264532494396e-4, 0., .238811521522368331303214066075e-3, 0., .135365534194038370983135068211e-2, 0., -.520710690438660595086839959882e-2, 0.,0.00341659266223572272892690737979 ]';
        b_3_def = [0., .379785635776247894184454273159e-7, 0., .655473977795402040043020497901e-7, 0., .103479954638984842226816620692e-6, 0., .173700624961660596894381303819e-6, 0., .337719613424065357737699682062e-6, 0., .877423283550614343733565759649e-6, 0., .515657204371051131603503471028e-5, 0.,-.203244736027387801432055290742e-4, 0.,0.0000134265158311651777460545854542 ]';
        b_4_def = [0., .703046511513775683031092069125e-13, 0., .110617117381148770138566741591e-12, 0., .146334657087392356024202217074e-12, 0., .184948444492681259461791759092e-12, 0., .231429962470609662207589449428e-12, 0., .291520062115989014852816512412e-12, 0., .373653379768759953435020783965e-12, 0., .491840460397998449461993671859e-12, 0., .671514395653454630785723660045e-12, 0., .963162916392726862525650710866e-12, 0., .147853378943890691325323722031e-11, 0., .250420145651013003355273649380e-11, 0., .495516257435784759806147914867e-11, 0., .130927034711760871648068641267e-10, 0., .779528640561654357271364483150e-10, 0., -.309866395328070487426324760520e-9, 0., .205572320292667201732878773151e-9]';

        % compute the coefficients
        V_1 = [ ones(size(xi_1)) xi_1 ];
        V_2 = [ ones(size(xi_2)) xi_2 ];
        V_3 = [ ones(size(xi_3)) xi_3 ];
        V_4 = [ ones(size(xi_4)) xi_4 ];
        xil = (xi_4 - 1) / 2; xir = (xi_4 + 1) / 2;
        Vl = [ ones(size(xil)) xil ]; Vr = [ ones(size(xir)) xir ];
        xx = linspace(-1,1,500)'; Vx = [ ones(size(xx)) xx ];
        for i=3:n(1)+1
            V_1 = [ V_1 , ((2*i-3) / (i-1) * xi_1 .* V_1(:,i-1) - (i-2) / (i-1) * V_1(:,i-2)) ];
        end;
        for i=3:n(2)+1
            V_2 = [ V_2 , ((2*i-3) / (i-1) * xi_2 .* V_2(:,i-1) - (i-2) / (i-1) * V_2(:,i-2)) ];
        end;
        for i=3:n(3)+1
            V_3 = [ V_3 , ((2*i-3) / (i-1) * xi_3 .* V_3(:,i-1) - (i-2) / (i-1) * V_3(:,i-2)) ];
        end;
        for i=3:n(4)+1
            V_4 = [ V_4 , ((2*i-3) / (i-1) * xi_4 .* V_4(:,i-1) - (i-2) / (i-1) * V_4(:,i-2)) ];
            Vx = [ Vx , ((2*i-3) / (i-1) * xx .* Vx(:,i-1) - (i-2) / (i-1) * Vx(:,i-2)) ];
            Vl = [ Vl , ((2*i-3) / (i-1) * xil .* Vl(:,i-1) - (i-2) / (i-1) * Vl(:,i-2)) ];
            Vr = [ Vr , ((2*i-3) / (i-1) * xir .* Vr(:,i-1) - (i-2) / (i-1) * Vr(:,i-2)) ];
        end;
        for i=1:n(1)+1, V_1(:,i) = V_1(:,i) * sqrt(4*i-2)/2; end;
        for i=1:n(2)+1, V_2(:,i) = V_2(:,i) * sqrt(4*i-2)/2; end;
        for i=1:n(3)+1, V_3(:,i) = V_3(:,i) * sqrt(4*i-2)/2; end;
        for i=1:n(4)+1
            V_4(:,i) = V_4(:,i) * sqrt(4*i-2)/2;
            Vx(:,i) = Vx(:,i) * sqrt(4*i-2)/2;
            Vl(:,i) = Vl(:,i) * sqrt(4*i-2)/2;
            Vr(:,i) = Vr(:,i) * sqrt(4*i-2)/2;
        end;
        V_1_inv = inv(V_1); V_2_inv = inv(V_2); V_3_inv = inv(V_3); V_4_inv = inv(V_4);
        Vcond = [ cond(V_1) , cond(V_2) , cond(V_3) , cond(V_4) ];

        % shift matrix
        T_left = V_4_inv * Vl; 
        T_right = V_4_inv * Vr;

        % compute the integral
        w = [ sqrt(2) , zeros(1,n(4)) ] / 2; % legendre
        
        % set-up the downdate matrix
        k = [0:n(4)]';
        U = diag(sqrt((k+1).^2 ./ (2*k+1) ./ (2*k+3))) + diag(sqrt(k(3:end).^2 ./ (4*k(3:end).^2-1)),2);
        
    end; % if exist('n')
        
    % create the original datatype
    ivals = struct( ...
        'a', [], 'b',[], ...
        'c', [], 'c_old', [], ...
        'fx', [], ...
        'int', [], ...
        'err', [], ...
        'tol', [], ...
        'depth', [], ...
        'rdepth', [], ...
        'ndiv', [] );
    
    % compute the first interval
    points = (a+b)/2 + (b-a) * xi_4 / 2;
    fx = f(points); nans = [];
    for i=1:length(fx), if ~isfinite(fx(i)), nans = [ i , nans ]; fx(i) = 0.0; end; end;
    ivals(1).c = zeros(n(4)+1,4);
    ivals(1).c(1:n(4)+1,4) = V_4_inv * fx;
    ivals(1).c(1:n(3)+1,3) = V_3_inv * fx([1:2:n(4)+1]);
    for i=nans, fx(i) = NaN; end;
    ivals(1).fx = fx;
    ivals(1).c_old = zeros(size(fx));
    ivals(1).a = a; ivals(1).b = b;
    ivals(1).int = (b-a) * w * ivals(1).c(:,4);
    c_diff = norm(ivals(1).c(:,4) - ivals(1).c(:,3));
    ivals(1).err = (b-a) * c_diff;
    if c_diff / norm(ivals(1).c(:,4)) > 0.1
        ivals(1).err = max( ivals(1).err , (b-a) * norm(ivals(1).c(:,4)) );
    end;
    ivals(1).tol = tol;
    ivals(1).depth = 4;
    ivals(1).ndiv = 0;
    ivals(1).rdepth = 1;
    
    % init some globals
    int = ivals(1).int; err = ivals(1).err; nr_ivals = 1;
    int_final = 0; err_final = 0; err_excess = 0;
    i_max = 1; nr_points = n(4)+1;
    ndiv_max = 20;
    
    % do we even need to go this way?
    if err < int * tol, return; end;
    
    % main loop
    while true
    
        % get some global data
        a = ivals(i_max).a; b = ivals(i_max).b;
        depth = ivals(i_max).depth; 
        split = 0;
        
        % depth of the old interval
        if depth == 1
            points = (a+b)/2 + (b-a)*xi_2/2;
            fx(1:2:n(2)+1) = ivals(i_max).fx;
            fx(2:2:n(2)) = f(points(2:2:n(2)));
            fx = fx(1:n(2)+1);
            nans = [];
            for i=1:length(fx), if ~isfinite(fx(i)), fx(i) = 0.0; nans = [ i , nans ]; end; end;
            c_new = V_2_inv * fx;
            if length(nans) > 0
                b_new = b_2_def; n_new = n(2);
                for i=nans
                    b_new(1:end-1) = (U(1:n(2)+1,1:n(2)+1) - diag(ones(n(2),1)*xi_2(i),1)) \ b_new(2:end);
                    b_new(end) = 0; c_new = c_new - c_new(n_new+1) / b_new(n_new+1) * b_new(1:end-1);
                    n_new = n_new - 1; fx(i) = NaN;
                end;
            end;
            ivals(i_max).fx = fx;
            nc = norm(c_new);
            nr_points = nr_points + n(2)-n(1);
            ivals(i_max).c(1:n(2)+1,2) = c_new;
            c_diff = norm(ivals(i_max).c(:,1) - ivals(i_max).c(:,2));
            ivals(i_max).err = (b-a) * c_diff;
            int_old = ivals(i_max).int; ivals(i_max).int = (b-a) * w(1) * c_new(1);
            if nc > 0 && c_diff / nc > 0.1
                split = 1;
            else
                ivals(i_max).depth = 2;
            end;
        elseif depth == 2
            points = (a+b)/2 + (b-a)*xi_3/2;
            fx(1:2:n(3)+1) = ivals(i_max).fx;
            fx(2:2:n(3)) = f(points(2:2:n(3)));
            fx = fx(1:n(3)+1);
            nans = [];
            for i=1:length(fx), if ~isfinite(fx(i)), fx(i) = 0.0; nans = [ i , nans ]; end; end;
            c_new = V_3_inv * fx;
            if length(nans) > 0
                b_new = b_3_def; n_new = n(3);
                for i=nans
                    b_new(1:end-1) = (U(1:n(3)+1,1:n(3)+1) - diag(ones(n(3),1)*xi_3(i),1)) \ b_new(2:end);
                    b_new(end) = 0; c_new = c_new - c_new(n_new+1) / b_new(n_new+1) * b_new(1:end-1);
                    n_new = n_new - 1; fx(i) = NaN;
                end;
            end;
            ivals(i_max).fx = fx;
            nc = norm(c_new);
            nr_points = nr_points + n(3)-n(2);
            ivals(i_max).c(1:n(3)+1,3) = c_new;
            c_diff = norm(ivals(i_max).c(:,2) - ivals(i_max).c(:,3));
            ivals(i_max).err = (b-a) * c_diff;
            int_old = ivals(i_max).int; ivals(i_max).int = (b-a) * w(1) * c_new(1);
            if nc > 0 && c_diff / nc > 0.1
                split = 1;
            else
                ivals(i_max).depth = 3;
            end;
        elseif depth == 3
            points = (a+b)/2 + (b-a)*xi_4/2;
            fx(1:2:n(4)+1) = ivals(i_max).fx;
            fx(2:2:n(4)) = f(points(2:2:n(4)));
            fx = fx(1:n(4)+1);
            nans = [];
            for i=1:length(fx), if ~isfinite(fx(i)), fx(i) = 0.0; nans = [ i , nans ]; end; end;
            c_new = V_4_inv * fx;
            if length(nans) > 0
                b_new = b_4_def; n_new = n(4);
                for i=nans
                    b_new(1:end-1) = (U(1:n(4)+1,1:n(4)+1) - diag(ones(n(4),1)*xi_4(i),1)) \ b_new(2:end);
                    b_new(end) = 0; c_new = c_new - c_new(n_new+1) / b_new(n_new+1) * b_new(1:end-1);
                    n_new = n_new - 1; fx(i) = NaN;
                end;
            end;
            ivals(i_max).fx = fx;
            nc = norm(c_new);
            nr_points = nr_points + n(4)-n(3);
            ivals(i_max).c(:,4) = c_new;
            c_diff = norm(ivals(i_max).c(:,3) - ivals(i_max).c(:,4));
            ivals(i_max).err = (b-a) * c_diff;
            int_old = ivals(i_max).int; ivals(i_max).int = (b-a) * w(1) * c_new(1);
            if nc > 0 && c_diff / nc > 0.1
                split = 1;
            else
                ivals(i_max).depth = 4;
            end;
        else
            split = 1;
        end;
        
        % can we safely ignore this interval?
        if points(2) <= points(1) || points(end) <= points(end-1) || ...
            ivals(i_max).err < abs(ivals(i_max).int) * eps * Vcond(ivals(i_max).depth)
            err_final = err_final + ivals(i_max).err;
            int_final = int_final + ivals(i_max).int;
            ivals(i_max) = ivals(nr_ivals);
            nr_ivals = nr_ivals - 1;
        elseif split
            m = (a + b) / 2;
            nr_ivals = nr_ivals + 1;
            ivals(nr_ivals).a = a; ivals(nr_ivals).b = m;
            ivals(nr_ivals).tol = ivals(i_max).tol / sqrt(2);
            ivals(nr_ivals).depth = 1; ivals(nr_ivals).rdepth = ivals(i_max).rdepth + 1;
            ivals(nr_ivals).c = zeros(n(4)+1,4);
            fx = [ ivals(i_max).fx(1) ; f((a+m)/2+(m-a)*xi_1(2:end-1)/2) ; ivals(i_max).fx((end+1)/2) ];
            nr_points = nr_points + n(1)-1; nans = [];
            for i=1:length(fx), if ~isfinite(fx(i)), fx(i) = 0.0; nans = [ i , nans ]; end; end;
            c_new = V_1_inv * fx;
            if length(nans) > 0
                b_new = b_1_def; n_new = n(1);
                for i=nans
                    b_new(1:end-1) = (U(1:n(1)+1,1:n(1)+1) - diag(ones(n(1),1)*xi_1(i),1)) \ b_new(2:end);
                    b_new(end) = 0; c_new = c_new - c_new(n_new+1) / b_new(n_new+1) * b_new(1:end-1);
                    n_new = n_new - 1; fx(i) = NaN;
                end;
            end;
            ivals(nr_ivals).fx = fx;
            ivals(nr_ivals).c(1:n(1)+1,1) = c_new; nc = norm(c_new);
            ivals(nr_ivals).c_old = T_left * ivals(i_max).c(:,depth);
            c_diff = norm(ivals(nr_ivals).c(:,1) - ivals(nr_ivals).c_old);
            ivals(nr_ivals).err = (m - a) * c_diff;
            ivals(nr_ivals).int = (m - a) * ivals(nr_ivals).c(1,1) * w(1);
            ivals(nr_ivals).ndiv = ivals(i_max).ndiv + (abs(ivals(i_max).c(1,1)) > 0 && ivals(nr_ivals).c(1,1) / ivals(i_max).c(1,1) > 2);
            if ivals(nr_ivals).ndiv > ndiv_max && 2*ivals(nr_ivals).ndiv > ivals(nr_ivals).rdepth, warning(sprintf('possibly divergent integral in the interval [%e,%e]! (h=%e)',a,m,m-a)); int = sign(int) * Inf; return; end;
            nr_ivals = nr_ivals + 1;
            ivals(nr_ivals).a = m; ivals(nr_ivals).b = b;
            ivals(nr_ivals).tol = ivals(i_max).tol / sqrt(2);
            ivals(nr_ivals).depth = 1; ivals(nr_ivals).rdepth = ivals(i_max).rdepth + 1;
            ivals(nr_ivals).c = zeros(n(4)+1,4);
            fx = [ ivals(i_max).fx((end+1)/2) ; f((m+b)/2+(b-m)*xi_1(2:end-1)/2) ; ivals(i_max).fx(end) ];
            nr_points = nr_points + n(1)-1; nans = [];
            for i=1:length(fx), if ~isfinite(fx(i)), fx(i) = 0.0; nans = [ i , nans ]; end; end;
            c_new = V_1_inv * fx;
            if length(nans) > 0
                b_new = b_1_def; n_new = n(1);
                for i=nans
                    b_new(1:end-1) = (U(1:n(1)+1,1:n(1)+1) - diag(ones(n(1),1)*xi_1(i),1)) \ b_new(2:end);
                    b_new(end) = 0; c_new = c_new - c_new(n_new+1) / b_new(n_new+1) * b_new(1:end-1);
                    n_new = n_new - 1; fx(i) = NaN;
                end;
            end;
            ivals(nr_ivals).fx = fx;
            ivals(nr_ivals).c(1:n(1)+1,1) = c_new; nc = norm(c_new);
            ivals(nr_ivals).c_old = T_right * ivals(i_max).c(:,depth);
            c_diff = norm(ivals(nr_ivals).c(:,1) - ivals(nr_ivals).c_old);
            ivals(nr_ivals).err = (b - m) * c_diff;
            ivals(nr_ivals).int = (b - m) * ivals(nr_ivals).c(1,1) * w(1);
            ivals(nr_ivals).ndiv = ivals(i_max).ndiv + (abs(ivals(i_max).c(1,1)) > 0 && ivals(nr_ivals).c(1,1) / ivals(i_max).c(1,1) > 2);
            if ivals(nr_ivals).ndiv > ndiv_max && 2*ivals(nr_ivals).ndiv > ivals(nr_ivals).rdepth, warning(sprintf('possibly divergent integral in the interval [%e,%e]! (h=%e)',m,b,b-m)); int = sign(int) * Inf; return; end;
            ivals(i_max) = ivals(nr_ivals);
            nr_ivals = nr_ivals - 1;
        end;
        
        % compute the running err and new max
        i_max = 1; i_min = 1; err = err_final; int = int_final;
        for i=1:nr_ivals
            if ivals(i).err > ivals(i_max).err, i_max = i;
            elseif ivals(i).err < ivals(i_min).err, i_min = i; end; 
            err = err + ivals(i).err;
            int = int + ivals(i).int;
        end;
        
        % nuke smallest element if stack is larger than 200
        if nr_ivals > 200
            err_final = err_final + ivals(i_min).err;
            int_final = int_final + ivals(i_min).int;
            ivals(i_min) = ivals(nr_ivals);
            if i_max == nr_ivals, i_max = i_min; end;
            nr_ivals = nr_ivals - 1;
        end;
        
        % get up and leave?
        if err == 0 || err < abs(int) * tol || (err_final > abs(int) * tol && err - err_final < abs(int) * tol) || nr_ivals == 0
            break;
        end;
        
    end; % main loop
    
    % clean-up and tabulate
    err = err + err_excess;
    
    
end


