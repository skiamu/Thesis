
function testsuite ( nr_runs )
%TESTSUITE tests the integrators QUAD, QUADL, QUADGK and INTEGRAL against
%   QUADCC and produces data for the "justdata.tex" file.
%
%   TESTSUITE(NR) performs NR runs. If NR is not specified, 1000 runs are
%   evaluated. Note that this may take some time.

    %% Lyness-Kaganove tests

    % wrappers and helpers
    v = version();
    if v(1) == '7'
        fflush = inline('dos(''sync'');','x');
    end;

    % nr of runs
    if nargin < 1
        nr_runs = 100;
    end
    
    % Set the tolerances we will be using
    tols = 10.^(-[3,6,9,12]);

    % tolerance
    % for tol = 1.0e-12
    for tol=tols

        % file descriptor for output
        fd = fopen( sprintf( 'testsuite_%.0e.tex' , tol ) ,'w+t');
        fdt = fopen('testsuite_times_%.0e.tex','w+t');

        % init rand
        rng(6178,'twister');

        %% Singularity
        hits = zeros(7,1); digits = zeros(7,1); fails = zeros(7,1);
        warns = zeros(7,2); nevals = zeros(7,1); time = zeros(7,1);
        for runs=1:nr_runs

            % set up integrand
            alpha = -rand(1)*0.5;
            lambda = rand(1);
            f = @(x) abs(x-lambda).^alpha;
            xx = linspace(0,1,500);
            % plot(xx,f(xx)); replot;

            % compute the exact solution
            int_exact = (lambda^(alpha+1) + (1-lambda)^(alpha+1)) / (alpha+1);
            rtol = abs(int_exact * tol);

            % call quad
            tic; [int_quad,neval_quad] = quad(f,0,1,rtol); time(1) = time(1) + toc;
            if length(lastwarn) > 10, warn_quad = 1; lastwarn(''); else warn_quad = 0; end;
            if abs((int_quad - int_exact)/int_exact) < tol
                if warn_quad, warns(1,1) = warns(1,1) + 1; end;
                hits(1) = hits(1) + 1;
            else
                if warn_quad, warns(1,2) = warns(1,2) + 1; end;
                fails(1) = fails(1) + 1;
            end;
            digits(1) = digits(1) + log10(max(abs((int_quad - int_exact)/int_exact),eps));
            nevals(1) = nevals(1) + neval_quad;

            % call quadl
            tic; [int_quadl,neval_quadl] = quadl(f,0,1,rtol); time(2) = time(2) + toc;
            if length(lastwarn) > 10, warn_quadl = 1; lastwarn(''); else warn_quadl = 0; end;
            if abs((int_quadl - int_exact)/int_exact) < tol
                if warn_quadl, warns(2,1) = warns(2,1) + 1; end;
                hits(2) = hits(2) + 1;
            else
                if warn_quadl, warns(2,2) = warns(2,2) + 1; end;
                fails(2) = fails(2) + 1;
            end;
            digits(2) = digits(2) + log10(max(abs((int_quadl - int_exact)/int_exact),eps));
            nevals(2) = nevals(2) + neval_quadl;

            % call quadgk
            tic; [int_quadgk,err_quadgk,neval_quadgk] = myquadgk(f,0,1,'AbsTol',rtol,'RelTol',tol); time(3) = time(3) + toc;
            if length(lastwarn) > 10, warn_quadgk = 1; lastwarn(''); else warn_quadgk = 0; end;
            if abs((int_quadgk - int_exact)/int_exact) < tol
                if warn_quadgk || err_quadgk > rtol, warns(3,1) = warns(3,1) + 1; end;
                hits(3) = hits(3) + 1;
            else
                if warn_quadgk || err_quadgk > rtol, warns(3,2) = warns(3,2) + 1; end;
                fails(3) = fails(3) + 1;
            end;
            digits(3) = digits(3) + log10(max(abs((int_quadgk - int_exact)/int_exact),eps));
            nevals(3) = nevals(3) + neval_quadgk;

            % call integralAls Actionheld kämpfte Arnold Schwarzenegger gegen das Böse. Heute will er mit seiner neuen Organisation von der Schweiz aus den Planeten retten. «Umweltschutz ist sexy \u2013 und lohnt sich», lautet seine 
            tic; [int_integral,neval_integral] = myintegral(f,0,1,'AbsTol',rtol,'RelTol',tol); time(4) = time(4) + toc;
            if length(lastwarn) > 10, warn_integral = 1; lastwarn(''); else warn_integral = 0; end;
            if abs((int_integral - int_exact)/int_exact) < tol
                if warn_integral, warns(4,1) = warns(4,1) + 1; end;
                hits(4) = hits(4) + 1;
            else
                if warn_integral, warns(4,2) = warns(4,2) + 1; end;
                fails(4) = fails(4) + 1;
            end;
            digits(4) = digits(4) + log10(max(abs((int_integral - int_exact)/int_exact),eps));
            nevals(4) = nevals(4) + neval_integral;

            % call quadcc
            tic; [int_quadcc,err_quadcc,neval_quadcc] = quadcc(f,0,1,tol); time(5) = time(5) + toc;
            if length(lastwarn) > 10, warn_quadcc = 1; fprintf(fd,'d) div. integral detected (alpha = %e, lambda = %e)...\n',alpha,lambda); lastwarn(''); else warn_quadcc = 0; end;
            if abs((int_quadcc - int_exact)/int_exact) < tol
                if warn_quadcc || err_quadcc > rtol, warns(5,1) = warns(5,1) + 1; end;
                hits(5) = hits(5) + 1;
            else
                if warn_quadcc || err_quadcc > rtol, warns(5,2) = warns(5,2) + 1; end;
                fails(5) = fails(5) + 1;
                % disp([alpha,lambda,int_exact,int_quadcc,err_quadcc]); return;
            end;
            digits(5) = digits(5) + log10(max(abs((int_quadcc - int_exact)/int_exact),eps));
            nevals(5) = nevals(5) + neval_quadcc;

            % verbosity
            if mod(runs,10) == 0
                disp(sprintf('testsuite: done %i runs...',runs));
            end;

        end;

        % dump the results
        fprintf(fd,'Eqn (\\ref{eqn:val_sing}) & $%d\\,(%d)$ & $%d\\,(%d)$ & $%.2f$ & $%d\\,(%d)$ & $%d\\,(%d)$ & $%.2f$ & $%d\\,(%d)$ & $%d\\,(%d)$ & $%.2f$ & $%d\\,(%d)$ & $%d\\,(%d)$ & $%.2f$ & $%d\\,(%d)$ & $%d\\,(%d)$ & $%.2f$ \\\\ \\hline\n', ...
            hits(1), warns(1,1), fails(1), warns(1,2), nevals(1)/nr_runs, ...
            hits(2), warns(2,1), fails(2), warns(2,2), nevals(2)/nr_runs, ...
            hits(3), warns(3,1), fails(3), warns(3,2), nevals(3)/nr_runs, ...
            hits(4), warns(4,1), fails(4), warns(4,2), nevals(4)/nr_runs, ...
            hits(5), warns(5,1), fails(5), warns(5,2), nevals(5)/nr_runs);
        fflush(fd);
        fprintf(fdt,'Eqn (\\ref{eqn:val_sing}) & $%.2f\\,(%.2f)$ & $%.2f\\,(%.2f)$ & $%.2f\\,(%.2f)$ & $%.2f\\,(%.2f)$ & $%.2f\\,(%.2f)$ & $%.2f\\,(%.2f)$ & $%.2f\\,(%.2f)$ \\\\ \\hline\n', ...
            time(1)/nr_runs*1000, time(1)/nevals(1)*1000000, ...
            time(2)/nr_runs*1000, time(2)/nevals(2)*1000000, ...
            time(3)/nr_runs*1000, time(3)/nevals(3)*1000000, ...
            time(4)/nr_runs*1000, time(4)/nevals(4)*1000000, ...
            time(5)/nr_runs*1000, time(5)/nevals(5)*1000000);
        fflush(fdt);
        
        % Store the data for the scatter plots
        data_quad = [ hits(1) , fails(1) , nevals(1)/nr_runs ];
        data_quadl = [ hits(2) , fails(2) , nevals(2)/nr_runs ];
        data_quadgk = [ hits(3) , fails(3) , nevals(3)/nr_runs ];
        data_integral = [ hits(4) , fails(4) , nevals(4)/nr_runs ];
        data_quadcc = [ hits(5) , fails(5) , nevals(5)/nr_runs ];


        % init rand
        rng(6178,'twister');

        %% Discontinuity
        hits = zeros(7,1); digits = zeros(7,1); fails = zeros(7,1);
        warns = zeros(7,2); nevals = zeros(7,1); time = zeros(7,1);
        for runs=1:nr_runs

            % set up integrand
            alpha = rand(1);
            lambda = rand(1);
            % f = inline(sprintf('(x>%22.16e).*exp(%22.16e*x)',lambda,alpha),'x');
            f = @(x) (x > lambda) .* exp( alpha * x );
            xx = linspace(0,1,500);
            % plot(xx,f(xx)); replot;

            % compute the exact solution
            if alpha < 0.01
                int_exact = 1-lambda+(-1/2*lambda^2+1/2+(1/6-1/6*lambda^3+(-1/24*lambda^4+1/24+(-1/120*lambda^5+1/120+(-1/720*lambda^6+1/720+(1/5040-1/5040*lambda^7+(-1/40320*lambda^8+1/40320+(-1/362880*lambda^9+1/362880)*alpha)*alpha)*alpha)*alpha)*alpha)*alpha)*alpha)*alpha;
            else
                int_exact = - (exp(alpha*lambda) - exp(alpha)) / alpha;
            end;
            rtol = abs(int_exact * tol);

            % call quad
            tic; [int_quad,neval_quad] = quad(f,0,1,rtol); time(1) = time(1) + toc;
            if length(lastwarn) > 10, warn_quad = 1; lastwarn(''); else warn_quad = 0; end;
            if abs((int_quad - int_exact)/int_exact) < tol
                if warn_quad, warns(1,1) = warns(1,1) + 1; end;
                hits(1) = hits(1) + 1;
            else
                if warn_quad, warns(1,2) = warns(1,2) + 1; end;
                fails(1) = fails(1) + 1;
            end;
            digits(1) = digits(1) + log10(max(abs((int_quad - int_exact)/int_exact),eps));
            nevals(1) = nevals(1) + neval_quad;

            % call quadl
            tic; [int_quadl,neval_quadl] = quadl(f,0,1,rtol); time(2) = time(2) + toc;
            if length(lastwarn) > 10, warn_quadl = 1; lastwarn(''); else warn_quadl = 0; end;
            if abs((int_quadl - int_exact)/int_exact) < tol
                if warn_quadl, warns(2,1) = warns(2,1) + 1; end;
                hits(2) = hits(2) + 1;
            else
                if warn_quadl, warns(2,2) = warns(2,2) + 1; end;
                fails(2) = fails(2) + 1;
            end;
            digits(2) = digits(2) + log10(max(abs((int_quadl - int_exact)/int_exact),eps));
            nevals(2) = nevals(2) + neval_quadl;

            % call quadgk
            tic; [int_quadgk,err_quadgk,neval_quadgk] = myquadgk(f,0,1,'AbsTol',rtol,'RelTol',tol); time(3) = time(3) + toc;
            if length(lastwarn) > 10, warn_quadgk = 1; lastwarn(''); else warn_quadgk = 0; end;
            if abs((int_quadgk - int_exact)/int_exact) < tol
                if warn_quadgk || err_quadgk > rtol, warns(3,1) = warns(3,1) + 1; end;
                hits(3) = hits(3) + 1;
            else
                if warn_quadgk || err_quadgk > rtol, warns(3,2) = warns(3,2) + 1; end;
                fails(3) = fails(3) + 1;
            end;
            digits(3) = digits(3) + log10(max(abs((int_quadgk - int_exact)/int_exact),eps));
            nevals(3) = nevals(3) + neval_quadgk;

            % call integral
            tic; [int_integral,neval_integral] = myintegral(f,0,1,'AbsTol',rtol,'RelTol',tol); time(4) = time(4) + toc;
            if length(lastwarn) > 10, warn_integral = 1; lastwarn(''); else warn_integral = 0; end;
            if abs((int_integral - int_exact)/int_exact) < tol
                if warn_integral, warns(4,1) = warns(4,1) + 1; end;
                hits(4) = hits(4) + 1;
            else
                if warn_integral, warns(4,2) = warns(4,2) + 1; end;
                fails(4) = fails(4) + 1;
            end;
            digits(4) = digits(4) + log10(max(abs((int_integral - int_exact)/int_exact),eps));
            nevals(4) = nevals(4) + neval_integral;

            % call quadcc
            tic; [int_quadcc,err_quadcc,neval_quadcc] = quadcc(f,0,1,tol); time(5) = time(5) + toc;
            if length(lastwarn) > 10, warn_quadcc = 1; fprintf(fd,'d) div. integral detected (alpha = %e, lambda = %e)...\n',alpha,lambda); lastwarn(''); else warn_quadcc = 0; end;
            if abs((int_quadcc - int_exact)/int_exact) < tol
                if warn_quadcc || err_quadcc > rtol, warns(5,1) = warns(5,1) + 1; end;
                hits(5) = hits(5) + 1;
            else
                if warn_quadcc || err_quadcc > rtol, warns(5,2) = warns(5,2) + 1; end;
                fails(5) = fails(5) + 1;
                % disp([alpha,lambda,int_exact,int_quadcc,err_quadcc]); return;
            end;
            digits(5) = digits(5) + log10(max(abs((int_quadcc - int_exact)/int_exact),eps));
            nevals(5) = nevals(5) + neval_quadcc;

            % verbosity
            if mod(runs,10) == 0
                disp(sprintf('testsuite: done %i runs...',runs));
            end;

        end;

        % dump the results
        fprintf(fd,'Eqn (\\ref{eqn:val_disc}) & $%d\\,(%d)$ & $%d\\,(%d)$ & $%.2f$ & $%d\\,(%d)$ & $%d\\,(%d)$ & $%.2f$ & $%d\\,(%d)$ & $%d\\,(%d)$ & $%.2f$ & $%d\\,(%d)$ & $%d\\,(%d)$ & $%.2f$ & $%d\\,(%d)$ & $%d\\,(%d)$ & $%.2f$ \\\\ \\hline\n', ...
            hits(1), warns(1,1), fails(1), warns(1,2), nevals(1)/nr_runs, ...
            hits(2), warns(2,1), fails(2), warns(2,2), nevals(2)/nr_runs, ...
            hits(3), warns(3,1), fails(3), warns(3,2), nevals(3)/nr_runs, ...
            hits(4), warns(4,1), fails(4), warns(4,2), nevals(4)/nr_runs, ...
            hits(5), warns(5,1), fails(5), warns(5,2), nevals(5)/nr_runs);
        fflush(fd);
        fprintf(fdt,'Eqn (\\ref{eqn:val_disc}) & $%.2f\\,(%.2f)$ & $%.2f\\,(%.2f)$ & $%.2f\\,(%.2f)$ & $%.2f\\,(%.2f)$ & $%.2f\\,(%.2f)$ & $%.2f\\,(%.2f)$ & $%.2f\\,(%.2f)$ \\\\ \\hline\n', ...
            time(1)/nr_runs*1000, time(1)/nevals(1)*1000000, ...
            time(2)/nr_runs*1000, time(2)/nevals(2)*1000000, ...
            time(3)/nr_runs*1000, time(3)/nevals(3)*1000000, ...
            time(4)/nr_runs*1000, time(4)/nevals(4)*1000000, ...
            time(5)/nr_runs*1000, time(5)/nevals(5)*1000000);
        fflush(fdt);
        
        % Store the data for the scatter plots
        data_quad = [ data_quad ; hits(1) , fails(1) , nevals(1)/nr_runs ];
        data_quadl = [ data_quadl ; hits(2) , fails(2) , nevals(2)/nr_runs ];
        data_quadgk = [ data_quadgk ; hits(3) , fails(3) , nevals(3)/nr_runs ];
        data_integral = [ data_integral ; hits(4) , fails(4) , nevals(4)/nr_runs ];
        data_quadcc = [ data_quadcc ; hits(5) , fails(5) , nevals(5)/nr_runs ];


        % init rand
        rng(6178,'twister');

        %% C_0 function
        hits = zeros(7,1); digits = zeros(7,1); fails = zeros(7,1);
        warns = zeros(7,2); nevals = zeros(7,1); time = zeros(7,1);
        for runs=1:nr_runs

            % set up integrand
            alpha = 4*rand(1);
            lambda = rand(1);
            % f = inline(sprintf('exp(-%22.16e*abs(x-%22.16e))',alpha,lambda),'x');
            f = @(x) exp( -alpha * abs(x - lambda) );
            xx = linspace(0,1,500);
            % plot(xx,f(xx)); replot;

            % compute the exact solution
            int_exact = -(-2+exp(-alpha*lambda)+exp(alpha*(lambda-1)))/alpha;
            rtol = abs(int_exact * tol);

            % call quad
            tic; [int_quad,neval_quad] = quad(f,0,1,rtol); time(1) = time(1) + toc;
            if length(lastwarn) > 10, warn_quad = 1; lastwarn(''); else warn_quad = 0; end;
            if abs((int_quad - int_exact)/int_exact) < tol
                if warn_quad, warns(1,1) = warns(1,1) + 1; end;
                hits(1) = hits(1) + 1;
            else
                if warn_quad, warns(1,2) = warns(1,2) + 1; end;
                fails(1) = fails(1) + 1;
            end;
            digits(1) = digits(1) + log10(max(abs((int_quad - int_exact)/int_exact),eps));
            nevals(1) = nevals(1) + neval_quad;

            % call quadl
            tic; [int_quadl,neval_quadl] = quadl(f,0,1,rtol); time(2) = time(2) + toc;
            if length(lastwarn) > 10, warn_quadl = 1; lastwarn(''); else warn_quadl = 0; end;
            if abs((int_quadl - int_exact)/int_exact) < tol
                if warn_quadl, warns(2,1) = warns(2,1) + 1; end;
                hits(2) = hits(2) + 1;
            else
                if warn_quadl, warns(2,2) = warns(2,2) + 1; end;
                fails(2) = fails(2) + 1;
            end;
            digits(2) = digits(2) + log10(max(abs((int_quadl - int_exact)/int_exact),eps));
            nevals(2) = nevals(2) + neval_quadl;

            % call quadgk
            tic; [int_quadgk,err_quadgk,neval_quadgk] = myquadgk(f,0,1,'AbsTol',rtol,'RelTol',tol); time(3) = time(3) + toc;
            if length(lastwarn) > 10, warn_quadgk = 1; lastwarn(''); else warn_quadgk = 0; end;
            if abs((int_quadgk - int_exact)/int_exact) < tol
                if warn_quadgk || err_quadgk > rtol, warns(3,1) = warns(3,1) + 1; end;
                hits(3) = hits(3) + 1;
            else
                if warn_quadgk || err_quadgk > rtol, warns(3,2) = warns(3,2) + 1; end;
                fails(3) = fails(3) + 1;
            end;
            digits(3) = digits(3) + log10(max(abs((int_quadgk - int_exact)/int_exact),eps));
            nevals(3) = nevals(3) + neval_quadgk;

            % call integral
            tic; [int_integral,neval_integral] = myintegral(f,0,1,'AbsTol',rtol,'RelTol',tol); time(4) = time(4) + toc;
            if length(lastwarn) > 10, warn_integral = 1; lastwarn(''); else warn_integral = 0; end;
            if abs((int_integral - int_exact)/int_exact) < tol
                if warn_integral, warns(4,1) = warns(4,1) + 1; end;
                hits(4) = hits(4) + 1;
            else
                if warn_integral, warns(4,2) = warns(4,2) + 1; end;
                fails(4) = fails(4) + 1;
            end;
            digits(4) = digits(4) + log10(max(abs((int_integral - int_exact)/int_exact),eps));
            nevals(4) = nevals(4) + neval_integral;

            % call quadcc
            tic; [int_quadcc,err_quadcc,neval_quadcc] = quadcc(f,0,1,tol); time(5) = time(5) + toc;
            if length(lastwarn) > 10, warn_quadcc = 1; fprintf(fd,'d) div. integral detected (alpha = %e, lambda = %e)...\n',alpha,lambda); lastwarn(''); else warn_quadcc = 0; end;
            if abs((int_quadcc - int_exact)/int_exact) < tol
                if warn_quadcc || err_quadcc > rtol, warns(5,1) = warns(5,1) + 1; end;
                hits(5) = hits(5) + 1;
            else
                if warn_quadcc || err_quadcc > rtol, warns(5,2) = warns(5,2) + 1; end;
                fails(5) = fails(5) + 1;
                % disp([alpha,lambda,int_exact,int_quadcc,err_quadcc]); return;
            end;
            digits(5) = digits(5) + log10(max(abs((int_quadcc - int_exact)/int_exact),eps));
            nevals(5) = nevals(5) + neval_quadcc;

            % verbosityAls Actionheld kämpfte Arnold Schwarzenegger gegen das Böse. Heute will er mit seiner neuen Organisation von der Schweiz aus den Planeten retten. «Umweltschutz ist sexy \u2013 und lohnt sich», lautet seine 
            if mod(runs,10) == 0
                disp(sprintf('testsuite: done %i runs...',runs));
            end;

        end;

        % dump the results
        fprintf(fd,'Eqn (\\ref{eqn:val_c0}) & $%d\\,(%d)$ & $%d\\,(%d)$ & $%.2f$ & $%d\\,(%d)$ & $%d\\,(%d)$ & $%.2f$ & $%d\\,(%d)$ & $%d\\,(%d)$ & $%.2f$ & $%d\\,(%d)$ & $%d\\,(%d)$ & $%.2f$ & $%d\\,(%d)$ & $%d\\,(%d)$ & $%.2f$ \\\\ \\hline\n', ...
            hits(1), warns(1,1), fails(1), warns(1,2), nevals(1)/nr_runs, ...
            hits(2), warns(2,1), fails(2), warns(2,2), nevals(2)/nr_runs, ...
            hits(3), warns(3,1), fails(3), warns(3,2), nevals(3)/nr_runs, ...
            hits(4), warns(4,1), fails(4), warns(4,2), nevals(4)/nr_runs, ...
            hits(5), warns(5,1), fails(5), warns(5,2), nevals(5)/nr_runs);
        fflush(fd);
        fprintf(fdt,'Eqn (\\ref{eqn:val_c0}) & $%.2f\\,(%.2f)$ & $%.2f\\,(%.2f)$ & $%.2f\\,(%.2f)$ & $%.2f\\,(%.2f)$ & $%.2f\\,(%.2f)$ & $%.2f\\,(%.2f)$ & $%.2f\\,(%.2f)$ \\\\ \\hline\n', ...
            time(1)/nr_runs*1000, time(1)/nevals(1)*1000000, ...
            time(2)/nr_runs*1000, time(2)/nevals(2)*1000000, ...
            time(3)/nr_runs*1000, time(3)/nevals(3)*1000000, ...
            time(4)/nr_runs*1000, time(4)/nevals(4)*1000000, ...
            time(5)/nr_runs*1000, time(5)/nevals(5)*1000000);
        fflush(fdt);

        % Store the data for the scatter plots
        data_quad = [ data_quad ; hits(1) , fails(1) , nevals(1)/nr_runs ];
        data_quadl = [ data_quadl ; hits(2) , fails(2) , nevals(2)/nr_runs ];
        data_quadgk = [ data_quadgk ; hits(3) , fails(3) , nevals(3)/nr_runs ];
        data_integral = [ data_integral ; hits(4) , fails(4) , nevals(4)/nr_runs ];
        data_quadcc = [ data_quadcc ; hits(5) , fails(5) , nevals(5)/nr_runs ];


        % init rand
        rng(6178,'twister');

        %% single peak
        hits = zeros(7,1); digits = zeros(7,1); fails = zeros(7,1);
        warns = zeros(7,2); nevals = zeros(7,1); time = zeros(7,1);
        for runs=1:nr_runs

            % set up integrand
            alpha = -3*rand(1) - 3;
            lambda = rand(1)+1;
            % f = inline(sprintf('10^%22.16e./((x-%22.16e).^2+10^(2*%22.16e))',alpha,lambda,alpha),'x');
            f = @(x) 10^alpha ./ ( (x-lambda).^2 + 10^(2*alpha) );
            xx = linspace(1,2,500);
            % plot(xx,f(xx)); replot;

            % compute the exact solution
            int_exact = 10^alpha*(atan((lambda-1)/sqrt(100^alpha))-atan((-2+lambda)/sqrt(100^alpha)))/sqrt(100^alpha);
            rtol = abs(int_exact * tol);

            % call quad
            tic; [int_quad,neval_quad] = quad(f,1,2,rtol); time(1) = time(1) + toc;
            if length(lastwarn) > 10, warn_quad = 1; lastwarn(''); else warn_quad = 0; end;
            if abs((int_quad - int_exact)/int_exact) < tol
                if warn_quad, warns(1,1) = warns(1,1) + 1; end;
                hits(1) = hits(1) + 1;
            else
                if warn_quad, warns(1,2) = warns(1,2) + 1; end;
                fails(1) = fails(1) + 1;
            end;
            digits(1) = digits(1) + log10(max(abs((int_quad - int_exact)/int_exact),eps));
            nevals(1) = nevals(1) + neval_quad;

            % call quadl
            tic; [int_quadl,neval_quadl] = quadl(f,1,2,rtol); time(2) = time(2) + toc;
            if length(lastwarn) > 10, warn_quadl = 1; lastwarn(''); else warn_quadl = 0; end;
            if abs((int_quadl - int_exact)/int_exact) < tol
                if warn_quadl, warns(2,1) = warns(2,1) + 1; end;
                hits(2) = hits(2) + 1;
            else
                if warn_quadl, warns(2,2) = warns(2,2) + 1; end;
                fails(2) = fails(2) + 1;
            end;
            digits(2) = digits(2) + log10(max(abs((int_quadl - int_exact)/int_exact),eps));
            nevals(2) = nevals(2) + neval_quadl;

            % call quadgk
            tic; [int_quadgk,err_quadgk,neval_quadgk] = myquadgk(f,1,2,'AbsTol',rtol,'RelTol',tol); time(3) = time(3) + toc;
            if length(lastwarn) > 10, warn_quadgk = 1; lastwarn(''); else warn_quadgk = 0; end;
            if abs((int_quadgk - int_exact)/int_exact) < tol
                if warn_quadgk || err_quadgk > rtol, warns(3,1) = warns(3,1) + 1; end;
                hits(3) = hits(3) + 1;
            else
                if warn_quadgk || err_quadgk > rtol, warns(3,2) = warns(3,2) + 1; end;
                fails(3) = fails(3) + 1;
            end;
            digits(3) = digits(3) + log10(max(abs((int_quadgk - int_exact)/int_exact),eps));
            nevals(3) = nevals(3) + neval_quadgk;

            % call integral
            tic; [int_integral,neval_integral] = myintegral(f,1,2,'AbsTol',rtol,'RelTol',tol); time(4) = time(4) + toc;
            if length(lastwarn) > 10, warn_integral = 1; lastwarn(''); else warn_integral = 0; end;
            if abs((int_integral - int_exact)/int_exact) < tol
                if warn_integral, warns(4,1) = warns(4,1) + 1; end;
                hits(4) = hits(4) + 1;
            else
                if warn_integral, warns(4,2) = warns(4,2) + 1; end;
                fails(4) = fails(4) + 1;
            end;
            digits(4) = digits(4) + log10(max(abs((int_integral - int_exact)/int_exact),eps));
            nevals(4) = nevals(4) + neval_integral;

            % call quadcc
            tic; [int_quadcc,err_quadcc,neval_quadcc] = quadcc(f,1,2,tol); time(5) = time(5) + toc;
            if length(lastwarn) > 10, warn_quadcc = 1; fprintf(fd,'d) div. integral detected (alpha = %e, lambda = %e)...\n',alpha,lambda); lastwarn(''); else warn_quadcc = 0; end;
            if abs((int_quadcc - int_exact)/int_exact) < tol
                if warn_quadcc || err_quadcc > rtol, warns(5,1) = warns(5,1) + 1; end;
                hits(5) = hits(5) + 1;
            else
                if warn_quadcc || err_quadcc > rtol, warns(5,2) = warns(5,2) + 1; end;
                fails(5) = fails(5) + 1;
                % disp([alpha,lambda,int_exact,int_quadcc,err_quadcc]); return;
            end;
            digits(5) = digits(5) + log10(max(abs((int_quadcc - int_exact)/int_exact),eps));
            nevals(5) = nevals(5) + neval_quadcc;

            % verbosity
            if mod(runs,10) == 0
                disp(sprintf('testsuite: done %i runs...',runs));
            end;

        end;

        % dump the results
        fprintf(fd,'Eqn (\\ref{eqn:val_peak}) & $%d\\,(%d)$ & $%d\\,(%d)$ & $%.2f$ & $%d\\,(%d)$ & $%d\\,(%d)$ & $%.2f$ & $%d\\,(%d)$ & $%d\\,(%d)$ & $%.2f$ & $%d\\,(%d)$ & $%d\\,(%d)$ & $%.2f$ & $%d\\,(%d)$ & $%d\\,(%d)$ & $%.2f$ \\\\ \\hline\n', ...
            hits(1), warns(1,1), fails(1), warns(1,2), nevals(1)/nr_runs, ...
            hits(2), warns(2,1), fails(2), warns(2,2), nevals(2)/nr_runs, ...
            hits(3), warns(3,1), fails(3), warns(3,2), nevals(3)/nr_runs, ...
            hits(4), warns(4,1), fails(4), warns(4,2), nevals(4)/nr_runs, ...
            hits(5), warns(5,1), fails(5), warns(5,2), nevals(5)/nr_runs);
        fflush(fd);
        fprintf(fdt,'Eqn (\\ref{eqn:val_peak}) & $%.2f\\,(%.2f)$ & $%.2f\\,(%.2f)$ & $%.2f\\,(%.2f)$ & $%.2f\\,(%.2f)$ & $%.2f\\,(%.2f)$ & $%.2f\\,(%.2f)$ & $%.2f\\,(%.2f)$ \\\\ \\hline\n', ...
            time(1)/nr_runs*1000, time(1)/nevals(1)*1000000, ...
            time(2)/nr_runs*1000, time(2)/nevals(2)*1000000, ...
            time(3)/nr_runs*1000, time(3)/nevals(3)*1000000, ...
            time(4)/nr_runs*1000, time(4)/nevals(4)*1000000, ...
            time(5)/nr_runs*1000, time(5)/nevals(5)*1000000);
        fflush(fdt);

        % Store the data for the scatter plots
        data_quad = [ data_quad ; hits(1) , fails(1) , nevals(1)/nr_runs ];
        data_quadl = [ data_quadl ; hits(2) , fails(2) , nevals(2)/nr_runs ];
        data_quadgk = [ data_quadgk ; hits(3) , fails(3) , nevals(3)/nr_runs ];
        data_integral = [ data_integral ; hits(4) , fails(4) , nevals(4)/nr_runs ];
        data_quadcc = [ data_quadcc ; hits(5) , fails(5) , nevals(5)/nr_runs ];


        % init rand
        rng(6178,'twister');

        %% Four peaks
        hits = zeros(7,1); digits = zeros(7,1); fails = zeros(7,1);
        warns = zeros(7,2); nevals = zeros(7,1); time = zeros(7,1);
        for runs=1:nr_runs

            % set up integrand
            alpha = -2*rand(1) - 3;
            lambda = rand(4,1)+1;
            % f = inline(sprintf('10^%22.16e * ( 1./((x-%22.16e).^2+10^(2*%22.16e)) + 1./((x-%22.16e).^2+10^(2*%22.16e)) + 1./((x-%22.16e).^2+10^(2*%22.16e)) + 1./((x-%22.16e).^2+10^(2*%22.16e)))',alpha,lambda(1),alpha,lambda(2),alpha,lambda(3),alpha,lambda(4),alpha),'x');
            f = @(x) 10^alpha * ( 1./((x-lambda(1)).^2+10^(2*alpha)) + 1./((x-lambda(2)).^2+10^(2*alpha)) + 1./((x-lambda(3)).^2+10^(2*alpha)) + 1./((x-lambda(4)).^2+10^(2*alpha)));
            xx = linspace(1,2,500);
            % plot(xx,f(xx)); replot;

            % compute the exact solution
            int_exact = 10^alpha*(atan((lambda(1)-1)/sqrt(100^alpha))-atan((-2+lambda(1))/sqrt(100^alpha)))/sqrt(100^alpha);
            int_exact = int_exact + 10^alpha*(atan((lambda(2)-1)/sqrt(100^alpha))-atan((-2+lambda(2))/sqrt(100^alpha)))/sqrt(100^alpha);
            int_exact = int_exact + 10^alpha*(atan((lambda(3)-1)/sqrt(100^alpha))-atan((-2+lambda(3))/sqrt(100^alpha)))/sqrt(100^alpha);
            int_exact = int_exact + 10^alpha*(atan((lambda(4)-1)/sqrt(100^alpha))-atan((-2+lambda(4))/sqrt(100^alpha)))/sqrt(100^alpha);
            rtol = abs(int_exact * tol);

            % call quad
            tic; [int_quad,neval_quad] = quad(f,0,1,rtol); time(1) = time(1) + toc;
            if length(lastwarn) > 10, warn_quad = 1; lastwarn(''); else warn_quad = 0; end;
            if abs((int_quad - int_exact)/int_exact) < tol
                if warn_quad, warns(1,1) = warns(1,1) + 1; end;
                hits(1) = hits(1) + 1;
            else
                if warn_quad, warns(1,2) = warns(1,2) + 1; end;
                fails(1) = fails(1) + 1;
            end;
            digits(1) = digits(1) + log10(max(abs((int_quad - int_exact)/int_exact),eps));
            nevals(1) = nevals(1) + neval_quad;

            % call quadl
            tic; [int_quadl,neval_quadl] = quadl(f,1,2,rtol); time(2) = time(2) + toc;
            if length(lastwarn) > 10, warn_quadl = 1; lastwarn(''); else warn_quadl = 0; end;
            if abs((int_quadl - int_exact)/int_exact) < tol
                if warn_quadl, warns(2,1) = warns(2,1) + 1; end;
                hits(2) = hits(2) + 1;
            else
                if warn_quadl, warns(2,2) = warns(2,2) + 1; end;
                fails(2) = fails(2) + 1;
            end;
            digits(2) = digits(2) + log10(max(abs((int_quadl - int_exact)/int_exact),eps));
            nevals(2) = nevals(2) + neval_quadl;

            % call quadgk
            tic; [int_quadgk,err_quadgk,neval_quadgk] = myquadgk(f,1,2,'AbsTol',rtol,'RelTol',tol); time(3) = time(3) + toc;
            if length(lastwarn) > 10, warn_quadgk = 1; lastwarn(''); else warn_quadgk = 0; end;
            if abs((int_quadgk - int_exact)/int_exact) < tol
                if warn_quadgk || err_quadgk > rtol, warns(3,1) = warns(3,1) + 1; end;
                hits(3) = hits(3) + 1;
            else
                if warn_quadgk || err_quadgk > rtol, warns(3,2) = warns(3,2) + 1; end;
                fails(3) = fails(3) + 1;
            end;
            digits(3) = digits(3) + log10(max(abs((int_quadgk - int_exact)/int_exact),eps));
            nevals(3) = nevals(3) + neval_quadgk;

            % call integral
            tic; [int_integral,neval_integral] = myintegral(f,1,2,'AbsTol',rtol,'RelTol',tol); time(4) = time(4) + toc;
            if length(lastwarn) > 10, warn_integral = 1; lastwarn(''); else warn_integral = 0; end;
            if abs((int_integral - int_exact)/int_exact) < tol
                if warn_integral, warns(4,1) = warns(4,1) + 1; end;
                hits(4) = hits(4) + 1;
            else
                if warn_integral, warns(4,2) = warns(4,2) + 1; end;
                fails(4) = fails(4) + 1;
            end;
            digits(4) = digits(4) + log10(max(abs((int_integral - int_exact)/int_exact),eps));
            nevals(4) = nevals(4) + neval_integral;

            % call quadcc
            tic; [int_quadcc,err_quadcc,neval_quadcc] = quadcc(f,1,2,tol); time(5) = time(5) + toc;
            if length(lastwarn) > 10, warn_quadcc = 1; fprintf(fd,'d) div. integral detected (alpha = %e, lambda = %e)...\n',alpha,lambda); lastwarn(''); else warn_quadcc = 0; end;
            if abs((int_quadcc - int_exact)/int_exact) < tol
                if warn_quadcc || err_quadcc > rtol, warns(5,1) = warns(5,1) + 1; end;
                hits(5) = hits(5) + 1;
            else
                if warn_quadcc || err_quadcc > rtol, warns(5,2) = warns(5,2) + 1; end;
                fails(5) = fails(5) + 1;
                % disp([alpha,lambda,int_exact,int_quadcc,err_quadcc]); return;
            end;
            digits(5) = digits(5) + log10(max(abs((int_quadcc - int_exact)/int_exact),eps));
            nevals(5) = nevals(5) + neval_quadcc;

            % verbosity
            if mod(runs,10) == 0
                disp(sprintf('testsuite: done %i runs...',runs));
            end;

        end;

        % dump the results
        fprintf(fd,'Eqn (\\ref{eqn:val_4peak}) & $%d\\,(%d)$ & $%d\\,(%d)$ & $%.2f$ & $%d\\,(%d)$ & $%d\\,(%d)$ & $%.2f$ & $%d\\,(%d)$ & $%d\\,(%d)$ & $%.2f$ & $%d\\,(%d)$ & $%d\\,(%d)$ & $%.2f$ & $%d\\,(%d)$ & $%d\\,(%d)$ & $%.2f$ \\\\ \\hline\n', ...
            hits(1), warns(1,1), fails(1), warns(1,2), nevals(1)/nr_runs, ...
            hits(2), warns(2,1), fails(2), warns(2,2), nevals(2)/nr_runs, ...
            hits(3), warns(3,1), fails(3), warns(3,2), nevals(3)/nr_runs, ...
            hits(4), warns(4,1), fails(4), warns(4,2), nevals(4)/nr_runs, ...
            hits(5), warns(5,1), fails(5), warns(5,2), nevals(5)/nr_runs);
        fflush(fd);
        fprintf(fdt,'Eqn (\\ref{eqn:val_4peak}) & $%.2f\\,(%.2f)$ & $%.2f\\,(%.2f)$ & $%.2f\\,(%.2f)$ & $%.2f\\,(%.2f)$ & $%.2f\\,(%.2f)$ & $%.2f\\,(%.2f)$ & $%.2f\\,(%.2f)$ \\\\ \\hline\n', ...
            time(1)/nr_runs*1000, time(1)/nevals(1)*1000000, ...
            time(2)/nr_runs*1000, time(2)/nevals(2)*1000000, ...
            time(3)/nr_runs*1000, time(3)/nevals(3)*1000000, ...
            time(4)/nr_runs*1000, time(4)/nevals(4)*1000000, ...
            time(5)/nr_runs*1000, time(5)/nevals(5)*1000000);
        fflush(fdt);

        % Store the data for the scatter plots
        data_quad = [ data_quad ; hits(1) , fails(1) , nevals(1)/nr_runs ];
        data_quadl = [ data_quadl ; hits(2) , fails(2) , nevals(2)/nr_runs ];
        data_quadgk = [ data_quadgk ; hits(3) , fails(3) , nevals(3)/nr_runs ];
        data_integral = [ data_integral ; hits(4) , fails(4) , nevals(4)/nr_runs ];
        data_quadcc = [ data_quadcc ; hits(5) , fails(5) , nevals(5)/nr_runs ];


        % init rand
        rng(6178,'twister');

        %% Oscillatory
        hits = zeros(7,1); digits = zeros(7,1); fails = zeros(7,1);
        warns = zeros(7,2); nevals = zeros(7,1); time = zeros(7,1);
        for runs=1:nr_runs

            % set up integrand
            alpha = 1.8 + rand(1)*0.4;
            lambda = rand(1);
            beta = 10^alpha / max(lambda^2,(1-lambda)^2);
            % f = inline(sprintf('2*%22.16e*(x-%22.16e).*cos(%22.16e*(x-%22.16e).^2)',beta,lambda,beta,lambda),'x');
            f = @(x) 2*beta*(x-lambda).*cos(beta*(x-lambda).^2);
            xx = linspace(0,1,500);
            % plot(xx,f(xx)); replot;

            % compute the exact solution
            int_exact = -sin(beta*lambda^2)+sin(beta*(lambda-1)^2);
            rtol = abs(int_exact * tol);

            % call quad
            tic; [int_quad,neval_quad] = quad(f,0,1,rtol); time(1) = time(1) + toc;
            if length(lastwarn) > 10, warn_quad = 1; lastwarn(''); else warn_quad = 0; end;
            if abs((int_quad - int_exact)/int_exact) < tol
                if warn_quad, warns(1,1) = warns(1,1) + 1; end;
                hits(1) = hits(1) + 1;
            else
                if warn_quad, warns(1,2) = warns(1,2) + 1; end;
                fails(1) = fails(1) + 1;
            end;
            digits(1) = digits(1) + log10(max(abs((int_quad - int_exact)/int_exact),eps));
            nevals(1) = nevals(1) + neval_quad;

            % call quadl
            tic; [int_quadl,neval_quadl] = quadl(f,0,1,rtol); time(2) = time(2) + toc;
            if length(lastwarn) > 10, warn_quadl = 1; lastwarn(''); else warn_quadl = 0; end;
            if abs((int_quadl - int_exact)/int_exact) < tol
                if warn_quadl, warns(2,1) = warns(2,1) + 1; end;
                hits(2) = hits(2) + 1;
            else
                if warn_quadl, warns(2,2) = warns(2,2) + 1; end;
                fails(2) = fails(2) + 1;
            end;
            digits(2) = digits(2) + log10(max(abs((int_quadl - int_exact)/int_exact),eps));
            nevals(2) = nevals(2) + neval_quadl;

            % call quadgk
            tic; [int_quadgk,err_quadgk,neval_quadgk] = myquadgk(f,0,1,'AbsTol',rtol,'RelTol',tol); time(3) = time(3) + toc;
            if length(lastwarn) > 10, warn_quadgk = 1; lastwarn(''); else warn_quadgk = 0; end;
            if abs((int_quadgk - int_exact)/int_exact) < tol
                if warn_quadgk || err_quadgk > rtol, warns(3,1) = warns(3,1) + 1; end;
                hits(3) = hits(3) + 1;
            else
                if warn_quadgk || err_quadgk > rtol, warns(3,2) = warns(3,2) + 1; end;
                fails(3) = fails(3) + 1;
            end;
            digits(3) = digits(3) + log10(max(abs((int_quadgk - int_exact)/int_exact),eps));
            nevals(3) = nevals(3) + neval_quadgk;

            % call integral
            tic; [int_integral,neval_integral] = myintegral(f,0,1,'AbsTol',rtol,'RelTol',tol); time(4) = time(4) + toc;
            if length(lastwarn) > 10, warn_integral = 1; lastwarn(''); else warn_integral = 0; end;
            if abs((int_integral - int_exact)/int_exact) < tol
                if warn_integral, warns(4,1) = warns(4,1) + 1; end;
                hits(4) = hits(4) + 1;
            else
                if warn_integral, warns(4,2) = warns(4,2) + 1; end;
                fails(4) = fails(4) + 1;
            end;
            digits(4) = digits(4) + log10(max(abs((int_integral - int_exact)/int_exact),eps));
            nevals(4) = nevals(4) + neval_integral;

            % call quadcc
            tic; [int_quadcc,err_quadcc,neval_quadcc] = quadcc(f,0,1,tol); time(5) = time(5) + toc;
            if length(lastwarn) > 10, warn_quadcc = 1; fprintf(fd,'d) div. integral detected (alpha = %e, lambda = %e)...\n',alpha,lambda); lastwarn(''); else warn_quadcc = 0; end;
            if abs((int_quadcc - int_exact)/int_exact) < tol
                if warn_quadcc || err_quadcc > rtol, warns(5,1) = warns(5,1) + 1; end;
                hits(5) = hits(5) + 1;
            else
                if warn_quadcc || err_quadcc > rtol, warns(5,2) = warns(5,2) + 1; end;
                fails(5) = fails(5) + 1;
                % disp([alpha,lambda,int_exact,int_quadcc,err_quadcc]); return;
            end;
            digits(5) = digits(5) + log10(max(abs((int_quadcc - int_exact)/int_exact),eps));
            nevals(5) = nevals(5) + neval_quadcc;

            % verbosity
            if mod(runs,10) == 0
                disp(sprintf('testsuite: done %i runs...',runs));
            end;

        end;

        % dump the results
        fprintf(fd,'Eqn (\\ref{eqn:val_oscill}) & $%d\\,(%d)$ & $%d\\,(%d)$ & $%.2f$ & $%d\\,(%d)$ & $%d\\,(%d)$ & $%.2f$ & $%d\\,(%d)$ & $%d\\,(%d)$ & $%.2f$ & $%d\\,(%d)$ & $%d\\,(%d)$ & $%.2f$ & $%d\\,(%d)$ & $%d\\,(%d)$ & $%.2f$ \\\\ \\hline\n', ...
            hits(1), warns(1,1), fails(1), warns(1,2), nevals(1)/nr_runs, ...
            hits(2), warns(2,1), fails(2), warns(2,2), nevals(2)/nr_runs, ...
            hits(3), warns(3,1), fails(3), warns(3,2), nevals(3)/nr_runs, ...
            hits(4), warns(4,1), fails(4), warns(4,2), nevals(4)/nr_runs, ...
            hits(5), warns(5,1), fails(5), warns(5,2), nevals(5)/nr_runs);
        fflush(fd);
        fprintf(fdt,'Eqn (\\ref{eqn:val_oscill}) & $%.2f\\,(%.2f)$ & $%.2f\\,(%.2f)$ & $%.2f\\,(%.2f)$ & $%.2f\\,(%.2f)$ & $%.2f\\,(%.2f)$ & $%.2f\\,(%.2f)$ & $%.2f\\,(%.2f)$ \\\\ \\hline\n', ...
            time(1)/nr_runs*1000, time(1)/nevals(1)*1000000, ...
            time(2)/nr_runs*1000, time(2)/nevals(2)*1000000, ...
            time(3)/nr_runs*1000, time(3)/nevals(3)*1000000, ...
            time(4)/nr_runs*1000, time(4)/nevals(4)*1000000, ...
            time(5)/nr_runs*1000, time(5)/nevals(5)*1000000);
        fflush(fdt);

        %% Close the files again
        fclose(fd);
        fclose(fdt);

        % Store the data for the scatter plots
        data_quad = [ data_quad ; hits(1) , fails(1) , nevals(1)/nr_runs ];
        data_quadl = [ data_quadl ; hits(2) , fails(2) , nevals(2)/nr_runs ];
        data_quadgk = [ data_quadgk ; hits(3) , fails(3) , nevals(3)/nr_runs ];
        data_integral = [ data_integral ; hits(4) , fails(4) , nevals(4)/nr_runs ];
        data_quadcc = [ data_quadcc ; hits(5) , fails(5) , nevals(5)/nr_runs ];
        
        %% Construct the scatter plot for this tolerance
        % avg nr of evals
        avgs = data_quad + data_quadl + data_quadgk + data_integral + data_quadcc;
        avgs = avgs / 5;

        % plot
        figure(1); loglog(data_quad(:,3)./avgs(:,3),data_quad(:,1)./avgs(:,1),'xk');
        hold on; loglog(data_quadl(:,3)./avgs(:,3),data_quadl(:,1)./avgs(:,1),'ok');
        hold on; loglog(data_quadgk(:,3)./avgs(:,3),data_quadgk(:,1)./avgs(:,1),'*k');
        hold on; loglog(data_integral(:,3)./avgs(:,3),data_integral(:,1)./avgs(:,1),'+k');
        hold on; loglog(data_quadcc(:,3)./avgs(:,3),data_quadcc(:,1)./avgs(:,1),'or');
        axis square; axis(axis .* [1/1.1,1.1,1/1.1,1.1]);
        legend('quad','quadl','quadgk','integral','quadcc','Location','Best');
        xlabel('nr. eval / avg. nr. eval'); ylabel('nr. correct / avg. nr. correct');
        title( sprintf( '\\tau = %.0e' , tol ) );
        set( gcf , 'PaperSize' , 12*[ 1 1 ] );
        set( gcf , 'PaperPosition' , 12*[ 0 0 1 1 ] );
        print( '-dpdf' , sprintf( 'scatter_%.0e.pdf' , tol ) );
        hold off; 

    end % loop over tolerances
    
    
    %% Battery tests

    % init the functions
    functions = [ 'f1 ' ; 'f2 ' ; 'f3 ' ; 'f4 ' ; 'f5 ' ; 'f6 ' ; 'f7 ' ; 'f8 ' ; ...
        'f9 ' ; 'f10' ; 'f11' ; 'f12' ; 'f13' ; 'f14' ; 'f15' ; 'f16' ; ...
        'f17' ; 'f18' ; 'f19' ; 'f20' ; 'f21' ; 'f22' ; 'f23' ; 'f24' ; ...
        'f25' ];
    ranges = [ 0 1 ; 0 1 ; 0 1 ; -1 1 ; -1 1 ; 0 1 ; 0 1 ; 0 1 ; 0 1 ; ...
        0 1 ; 0 1 ; 0 1 ; 0 1 ; 0 10 ; 0 10 ; 0 10 ; 0 1 ; 0 pi ; ...
        0 1 ; -1 1 ; 0 1 ; 0 1 ; 0 1 ; 0 3 ; 0 5 ];
    f_exact = [ 1.7182818284590452354 , 0.7 , 2/3 , 0.4794282266888016674 , ...
        1.5822329637296729331 , 0.4 , 2 , 0.86697298733991103757 , ...
        1.1547005383792515290 , 0.69314718055994530942 , 0.3798854930417224753 , ...
        0.77750463411224827640 , 0.49898680869304550249 , ...
        0.5 , 1 , 0.13263071079267703209e+08 , 0.49898680869304550249 , ...
        0.83867634269442961454 , -1 , 1.5643964440690497731 , ...
        0.16349494301863722618 , -0.63466518254339257343 , ...
        0.013492485649467772692 , 17.664383539246514971 , 7.5 ];
    
    % version-specific stuff    
    v = version();
    if v(1) == '7'
        fflush = inline('dos(''sync'');','x');
        replot = inline('drawnow;');
    end;
    
    % init result data
    neval = zeros(length(functions),length(tols),5);
    fails = zeros(length(functions),length(tols),5);
    warns = zeros(length(functions),length(tols),5);
    lastwarn('');
    
    % loop over all functions
    for fid=1:length(functions)
    
        % init some stuff, plot the function
        a = ranges(fid,1); b = ranges(fid,2);
        f = eval(sprintf('@%s',functions(fid,:)));
        xx = linspace(a,b,500); plot(xx,f(xx),'-g'); replot;

        % loop over all tolerances
        for tid=1:length(tols)

            % verbosity
            disp(sprintf('calling %s for tols=%e...',functions(fid,:),tols(tid)));
        
            % required absolute tolerance
            rtol = abs(f_exact(fid)*tols(tid));

            % call quad
            [int_quad,neval_quad] = quad(f,a,b,rtol);
            if length(lastwarn) > 10, warns(fid,tid,1) = 1; warning(''); end;
            if ~isfinite(int_quad) || abs(int_quad - f_exact(fid)) > rtol, fails(fid,tid,1) = 1; end;
            neval(fid,tid,1) = neval_quad;
                    
            % call quadl
            [int_quadl,neval_quadl] = quadl(f,a,b,rtol);
            if length(lastwarn) > 10, warns(fid,tid,1) = 1; warning(''); end;
            if ~isfinite(int_quadl) || abs(int_quadl - f_exact(fid)) > rtol, fails(fid,tid,2) = 1; end;
            neval(fid,tid,2) = neval_quadl;
                    
            % call quadgk
            [int_quadgk,err_quadgk,neval_quadgk] = myquadgk(f,a,b,'AbsTol',rtol,'RelTol',tols(tid));
            if length(lastwarn) > 10 || err_quadgk > rtol, warns(fid,tid,5) = 1; warning(''); end;
            if ~isfinite(int_quadgk) || abs(int_quadgk - f_exact(fid)) > rtol, fails(fid,tid,3) = 1; end;
            neval(fid,tid,3) = neval_quadgk;
                    
            % call integral
            [int_integral,neval_integral] = myintegral(f,a,b,'AbsTol',rtol,'RelTol',tols(tid));
            if length(lastwarn) > 10, warns(fid,tid,5) = 1; warning(''); end;
            if ~isfinite(int_integral) || abs(int_integral - f_exact(fid)) > rtol, fails(fid,tid,4) = 1; end;
            neval(fid,tid,4) = neval_integral;
                    
            % call quadcc
            [int_quadcc,err_quadcc,neval_quadcc] = quadcc(f,a,b,tols(tid));
            if length(lastwarn) > 10 || err_quadcc > rtol, warns(fid,tid,4) = 1; warning(''); end;
            if ~isfinite(int_quadcc) || abs(int_quadcc - f_exact(fid)) > rtol, fails(fid,tid,5) = 1; end;
            neval(fid,tid,5) = neval_quadcc;
                    
        end % loop over all tolerances
    
    end % loop over all functions
    
    % Open a file for writing
    fd = fopen( 'battery.tex' , 'w+t' );
    
    % print some nice output
    for fid = 1:length(functions)
        fprintf(fd,'%s ',functions(fid,:));
        for tid = 1:length(tols)
            for i=1:5
                if fails(fid,tid,i)
                    fprintf(fd,' & \\st{%d}',neval(fid,tid,i));
                elseif neval(fid,tid,i)/(1-fails(fid,tid,i)) == min(neval(fid,tid,:)./(1-fails(fid,tid,:)))
                    fprintf(fd,' & {\\bf %d}',neval(fid,tid,i));
                else
                    fprintf(fd,' & %d',neval(fid,tid,i));
                end;
            end;
        end; % for tid
        fprintf( fd , '\\\\\n' );
    end; % for fid
    
    % close the file
    fclose(fd);
    
end % testsuite
    
    
function [ q , errbnd , nr_points ] = myquadgk ( f , varargin )

    % Wrapper function to integrate
    function out = wrapper ( x )
        nr_points = nr_points + prod( size( x ) );
        out = f( x );
    end
    
    % Init the counter
    nr_points = 0;
    
    % Call the quadrature routine with the wrapper
    [ q , errbnd ] = quadgk( @wrapper , varargin{:} );
    
end % myquadgk


function [ q , nr_points ] = myintegral ( f , varargin )

    % Wrapper function to integrate
    function out = wrapper ( x )
        nr_points = nr_points + prod( size( x ) );
        out = f( x );
    end
    
    % Init the counter
    nr_points = 0;
    
    % Call the quadrature routine with the wrapper
    q = integral( @wrapper , varargin{:} );
    
end % myintegral


%-------------------------------------------------------------------------------
%  define the test functions for the battery test

function res = f1 ( x )
    res = exp(x);
end

    
function res = f2 ( x )
    res = double(x >= 0.3);
end


function res = f3 ( x )
    res = sqrt(x);
end


function res = f4 ( x )
    res = (23/25) * cosh(x) - cos(x);
end


function res = f5 ( x )
    res = 1 ./ (x.^4 + x.^2 + 0.9);
end


function res = f6 ( x )
    res = sqrt( x.^3 );
end


function res = f7 ( x )
    res = 1 ./ sqrt(x);
end


function res = f8 ( x )
    res = 1 ./ (1 + x.^4);
end


function res = f9 ( x )
    res = 2 ./ (2 + sin(10*pi*x));
end


function res = f10 ( x )
    res = 1 ./ (1 + x);
end


function res = f11 ( x )
    res = 1 ./ (1 + exp(x));
end


function res = f12 ( x )
    res = x ./ (exp(x) - 1);
    % res = zeros(size(x));
    % for k=1:length(x)
    %     if (x(k) > 1.0e-3)
    %         res(k) = x(k) / (exp(x(k)) - 1);
    %     else
    %         res(k) = 1.0 - 1/2*x(k) + 1/12*x(k)^2 - 1/720*x(k)^4 + 1/30240*x(k)^6 - 1/1209600*x(k)^8;
    %     end
    % end
end


function res = f13 ( x )
    res = sin(100 * pi * x) ./ (pi * x);
end


function res = f14 ( x )
    res = sqrt(50) * exp(-50*pi*x.^2);
end


function res = f15 ( x )
    res = 25 * exp(-25*x);
end


function res = f16 ( x )
    res = 50 / pi * (2500 * x.^2 + 1);
end


function res = f17 ( x )
    res = 50 * (sin(50*pi*x) ./ (50*pi*x)).^2;
end


function res = f18 ( x )
    res = cos( cos(x) + 3*sin(x) + 2*cos(2*x) + 3*sin(2*x) + 3*cos(3*x) );
end


function res = f19 ( x )
    res = log(x);
end


function res = f20 ( x )
    res = 1 ./ (x.^2 + 1.005);
end


function res = f21 ( x )
    res = 1 ./ cosh( 10 * (x - 0.2) * 2 ) + ...
          1 ./ cosh( 100 * (x - 0.4) * 4 ) + ...
          1 ./ cosh( 1000 * (x - 0.6) * 8 );
end


function res = f22 ( x )
    res = 4 * pi^2 * x .* sin(20*pi*x) .* cos(2*pi*x);
end


function res = f23 ( x )
    res = 1 ./ (1 + (230*x - 30).^2);
end


function res = f24 ( x )
    res = floor(exp(x));
end


function res = f25 ( x )
    res = (x < 1) .* (x + 1) + ...
          (1 <= x & x <= 3) .* (3 - x) + ...
          (x > 3) * 2;
end

