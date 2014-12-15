function [ot, ou, oh, process] = trajopt_midpoint(plant, N, Tmin, Tmax, xDim, uDim, jcost, startpos, endpos)
%% Initialization of parameters
umax = 10;
alpha = 0.1; betasucc = 1.2; betafail = 0.5; rhoinit = 1;
lambda = 200;
tlambda = 10;
lambda_inc = 1.5;
Kp = 60;
Kc = 200;
Kt = 50;
rhothres = 0.0001;
xtol = 0.0001;
ftol = 0.0001;
ctol = 0.000001;
phis = [];
phihats = [];
etas = [];
etahats = [];
tcosts = [];
jcosts = [];
ucosts = [];
rhos = [];
deltahats = [];
deltas = [];
converge_cnt = 0;
%% Set the inital trajectory.
ot = [startpos rand(xDim,N-2)*2*pi-pi endpos];
ou = rand(uDim,N)*10-5;
oT = 4;
oh = ones(1,N-1)*oT/(N-1);
%initial value of phi.
rho = rhoinit;
%% Convex sequential optimization
cvx_quiet(true);
tic
for epoch = 1:Kp
    rho = rhoinit;
    for j = 1:Kc
        eta = zeros(xDim, N);
        A = cell(N-1,1);
        B = cell(N-1,1);
        S = cell(N-1,1);
        F = cell(N-1,1);
        for t = 1:N-1
            x0 = ot(:,t); u0 = ou(:,t);
            x1 = ot(:,t+1); u1 = ou(:,t+1);
            [f,df] = plant.dynamics(0,(x0 + x1)/2,(u0 + u1)/2);
            F{t} = f;
            A{t} = full(df(:,2:1+xDim));
            B{t} = full(df(:,2+xDim:end));
            eta(:,t) = x1 - x0 - oh(t)*f;
        end
        oldphi = tlambda*sum(oh)+(jcost(ot) + sum_square(ou(:))) + lambda*sum(abs(eta(:)));
        
        for i = 1:Kt
            cvx_begin
            variable nt(xDim, N);
            variable u(uDim, N);
            variable h(1,N-1);
            etahat = cvx(zeros(xDim, N-1));
            % Initial and final conditions.
            nt(:,1) == startpos;
            nt(:,N) == endpos;
            for t = 1:N-1
                x0 = ot(:,t); u0 = ou(:,t);
                x1 = ot(:,t+1); u1 = ou(:,t+1);
                etahat(:,t) = nt(:,t+1) - nt(:,t) - F{t}*(h(t)) - A{t}*(nt(:,t) +nt(:,t+1)-x0-x1)/2*oh(t) - B{t}*(u(:,t)+u(:,t+1)-u0-u1)/2*oh(t);
                %etahat(:,t) == 0;
            end
            h >= Tmin/(N-1);
            h <= Tmax/(N-1);
            % Trust region constraints.
            abs(nt(:) - ot(:)) <= rho;
            abs(u(:) - ou(:)) <= rho;
            abs(h(:) - oh(:)) <= rho;
            % Torque limit.
            abs(u(:)) <= umax;
            minimize(tlambda*sum(h)+(jcost(nt) + sum_square(u(:)))+ lambda*sum(abs(etahat(:))));
            cvx_end
            
            % Calculate actual costs.
            eta = zeros(xDim, N-1);
            for t = 1:N-1
                % Penalty for ynamics constraints
                x0 = nt(:,t); u0 = u(:,t);
                x1 = nt(:,t+1); u1 = u(:,t+1);
                [f] = plant.dynamics(0,(x0+x1)/2,(u0+u1)/2);
                eta(:,t) = nt(:,t+1) - nt(:,t) - h(t)*f;
            end
            cons = sum(abs(eta(:)));
            phi = tlambda*sum(h)+(jcost(nt) + sum_square(u(:))/4)+ lambda*cons;
            delta = oldphi - phi;
            phihat = cvx_optval;
            deltahat = oldphi - phihat;
            phis = [phis phi];
            phihats = [phihats phihat];
            etas = [etas sum(abs(eta(:)))];
            etahats = [etahats sum(abs(etahat(:)))];
            jcosts = [jcosts jcost(ot)];
            ucosts = [ucosts sum_square(ou(:))];
            deltas = [deltas delta];
            deltahats = [deltahats deltahat];
            tcosts = [tcosts sum(oh)];
            deltax = sum(abs(nt(:) - ot(:)));
            disp(cvx_status);
            disp('   epoch       j          i     etahat   eta      rho')
            disp([epoch j i sum(abs(etahat(:))) cons  rho])
            disp('  phihat       phi  deltahat     delta       deltax      rho')
            disp([phihat phi deltahat delta deltax rho])
            
            rhos = [rhos rho];
            if isnan(phihat) || isinf(phihat) || delta <= alpha*deltahat
                rho = betafail*rho;
            else
                rho = betasucc*rho;
                % Accept the solution
                ot = nt; ou = u; oh = h; oT = sum(h);
                %oldphi = phi;
                break;
            end
            if(rho < rhothres)
                break;
            end
        end
        if (abs(deltax) < xtol || (abs(delta) < ftol && abs(deltahat) < ftol) || rho < rhothres)
            break;
        end
    end
    if (cons <= ctol)
        break;
    else
        lambda = lambda*lambda_inc;
    end
end
tot = toc
process.phis = phis;
process.jcosts = jcosts;
process.ucosts = ucosts;
process.tcosts = tcosts;
process.rhos = rhos;
process.etas = etas;
process.phihats = phihats;
process.etahats = etahats;
process.deltas = deltas;
process.deltahats = deltahats;