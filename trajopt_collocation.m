function [ot, ou, oh, process] = trajopt_collocation(plant, N, Tmin, Tmax, xDim, uDim, jcost, startpos, endpos)
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

xdot = cell(N,1);
dxdot = cell(N,1);
for t = 1:N
    [f,df] = plant.dynamics(t,ot(:,t),ou(:,t));
    xdot{t} = f;
    dxdot{t} = df;
end
eta = zeros(xDim, N);
F = cell(N-1,1);
dF = cell(N-1,1);
for t = 1:N-1
    x0 = ot(:,t); u0 = ou(:,t);
    x1 = ot(:,t+1); u1 = ou(:,t+1);
    [f,df] = col_constraint_fun(plant,oh(t),x0,x1,u0,u1,xdot{t},dxdot{t},xdot{t+1},dxdot{t+1});
    F{t} = f;
    dF{t} = df;
    eta(:,t) = f;
end
oldphi = tlambda*sum(oh)+(jcost(ot) + sum_square(ou(:))) + lambda*sum(abs(eta(:)));
converge_cnt = 0;
%% Convex sequential optimization
tic
for epoch = 1:Kp
    %rho = rhoinit;
    for j = 1:Kc
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
                etahat(:,t) = F{t} + dF{t}*[h(t)-oh(t);nt(:,t)-ot(:,t);nt(:,t+1)-ot(:,t+1);u(:,t)-ou(:,t);u(:,t+1)-ou(:,t+1)];
            end
            h >= Tmin/(N-1);
            h <= Tmax/(N-1);
            % Trust region constraints.
            abs(nt(:) - ot(:)) <= rho;
            abs(u(:) - ou(:)) <= rho;
            abs(h(:) - oh(:)) <= rho;
            % Torque limit.
            abs(u(:)) <= umax-0.1;
            minimize(tlambda*sum(h)+(jcost(nt) + sum_square(u(:)))+ lambda*sum(abs(etahat(:))));
            cvx_end
            if isnan(cvx_optval) || isinf(cvx_optval)
                rho = betafail*rho;
                continue;
            end
            % Calculate actual costs.
            xdot = cell(N,1);
            dxdot = cell(N,1);
            for t = 1:N
                [f,df] = plant.dynamics(0,nt(:,t),u(:,t));
                xdot{t} = f;
                dxdot{t} = df;
            end
            eta = zeros(xDim, N);
            F2 = cell(N-1,1);
            dF2 = cell(N-1,1);
            for t = 1:N-1
                x0 = nt(:,t); u0 = u(:,t);
                x1 = nt(:,t+1); u1 = u(:,t+1);
                [f,df] = col_constraint_fun(plant,h(t),x0,x1,u0,u1,xdot{t},dxdot{t},xdot{t+1},dxdot{t+1});
                F2{t} = f;
                dF2{t} = df;
                eta(:,t) = f;
            end
            
            cons = sum(abs(eta(:)));
            phi = tlambda*sum(h)+(jcost(nt) + sum_square(u(:)))+ lambda*cons;
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
            if  delta < 0 || delta <= alpha*deltahat
                rho = betafail*rho;
            else
                disp('accept')
                rho = betasucc*rho;
                % Accept the solution
                ot = nt; ou = u; oh = h; oT = sum(h);
                oldphi = phi;
                F = F2; dF = dF2;
                break;
            end
            if(rho < rhothres)
                break;
            end
        end
        if (rho < rhothres)
            break;
        end
        if (abs(deltax) < xtol || (abs(delta) < ftol && abs(deltahat) < ftol))
            converge_cnt = converge_cnt + 1;
            if converge_cnt >= 5
                break;
            end
        end
    end
    if (cons <= ctol)
        sco_status = 'success';
        break;
    else
        sco_status = 'fail';
        if (cons <= ctol*10)
            sco_status = 'partial success';
        end
        lambda = lambda*lambda_inc;
        oldphi = tlambda*sum(oh)+(jcost(ot) + sum_square(ou(:))) + lambda*sum(abs(eta(:)));
    end
end
toc
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