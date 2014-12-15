%% Get the plant and visualizer
NL = 3;
plant = PlanarNLink(NL);
v = plant.constructVisualizer();
v.axis=[-1.2 1.2 -1.2 1.2]+(NL-1)*[-1 1 -1.2 1.2];
xG = [pi;zeros(2*NL-1,1)];

%% Initialization of parameters
N=41;
% squared displacement matrix
Dis=eye(N)*2;
Dis(1)=1;
Dis(N,N)=1;
Dis(2:N,1:N-1)=Dis(2:N,1:N-1)+eye(N-1)*-1;
Dis(1:N-1,2:N)=Dis(1:N-1,2:N)+eye(N-1)*-1;
jcost = @(nt) 0;
if NL == 2
    jcost = @(nt) nt(1,:)*Dis*nt(1,:)' + nt(2,:)*Dis*nt(2,:)';
end
if NL == 3
    jcost = @(nt) nt(1,:)*Dis*nt(1,:)' + nt(2,:)*Dis*nt(2,:)' +nt(3,:)*Dis*nt(3,:)';
end
if NL == 4
    jcost = @(nt) nt(1,:)*Dis*nt(1,:)' + nt(2,:)*Dis*nt(2,:)' +nt(3,:)*Dis*nt(3,:)' +nt(4,:)*Dis*nt(4,:)';
end
Tmin = 2; Tmax = 8; xDim = NL*2; uDim = NL-1;
startpos = zeros(xDim,1); endpos = double(xG);

%% call trajopt_forward/backward/midpoint/collocation
%[ot, ou, oh, process] = trajopt_forward(plant, N, Tmin, Tmax, xDim, uDim, jcost, startpos, endpos);
%[ot, ou, oh, process] = trajopt_midpoint(plant, N, Tmin, Tmax, xDim, uDim, jcost, startpos, endpos);
[ot, ou, oh, process] = trajopt_collocation(plant, N, Tmin, Tmax, xDim, uDim, jcost, startpos, endpos);

%% show results
phis = process.phis;
jcosts = process.jcosts;
ucosts = process.ucosts;
tcosts = process.tcosts;
rhos = process.rhos;
etas = process.etas;
phihats = process.phihats;
etahats = process.etahats;
deltas = process.deltas;
deltahats = process.deltahats;

%% Play back the trajectory
T = sum(oh);

x = ot;
xtraj = PPTrajectory(foh([0 cumsum(oh)],x));
xtraj = xtraj.setOutputFrame(plant.getStateFrame);
v.playback(xtraj)

% Save swf file
%v.playbackSWF(xtraj,'test.swf');
%% Play back the trajectory simulated with inputs
utraj = PPTrajectory(foh([0 cumsum(oh)],ou));
utraj = utraj.setOutputFrame(plant.getInputFrame);
sys = cascade(utraj,plant);
uxtraj=simulate(sys,utraj.tspan,zeros(xDim,1));
v.playback(uxtraj)

%% Time-varying LQR
Q = diag([ones(NL,1)*10;ones(NL,1)]);
R = diag(ones(NL-1,1));
c = plant.tvlqr(xtraj,utraj,Q,R,diag([ones(NL*2,1)*10]));
sys_cl = plant.feedback(c);

x0_test = zeros(NL*2,1)%+.02*(rand(4,1) - 1);
%x0_test(3) = x0(3)-0.1;
traj_test_1 = sys_cl.simulate(xtraj.tspan,x0_test);
v.playback(traj_test_1);
%v.playbackSWF(traj_test_1,'acrobot_mid_lqrtest2.swf');

%% plot produce graphs.
% phi.
colors = [ 255 0 0; 0 0 0; 0 0 255;  245 168 36; 159 174 29; 0 255 0; 167 0 130; 255 242 17;]/255;
Kmax = size(phis,2);
figure(1); cla;
plot(phis, 'Color', colors(1,:), 'LineWidth', 1);hold on;
plot(phihats,'Color', colors(3,:), 'LineWidth', 1)
axis auto; a = axis; a(1) = 1; a(2) = Kmax; axis(a);
xlabel('k','FontSize',18); ylabel('phi','FontSize',18);
legend('phi','phihat');
set(gca, 'FontSize', 18);
print -depsc phis.eps

% Jcost.
figure(2);
plot(jcosts, 'Color', colors(1,:), 'LineWidth', 1)
axis auto; a = axis; a(1) = 1; a(2) = Kmax; axis(a);
ylabel('Jx','FontSize',18); xlabel('k','FontSize',18);
set(gca, 'FontSize', 18);
print -depsc Jcosts.eps

% Torque cost.
figure(3);
semilogy(ucosts, 'Color', colors(1,:), 'LineWidth', 1)
axis auto; a = axis; a(1) = 1; a(2) = Kmax; axis(a);
xlabel('k','FontSize',18); ylabel('Ju','FontSize',18);
set(gca, 'FontSize', 18);
print -depsc ucosts.eps

% rho (in degrees).
figure(4);
semilogy(rhos, 'Color', colors(1,:), 'LineWidth', 1);
axis auto; a = axis; a(1) = 1; a(2) = Kmax; axis(a);
xlabel('k','FontSize',18); ylabel('rho','FontSize',18);
set(gca, 'FontSize', 18);
print -depsc rhos.eps

% eta.
figure(5);
plot(etas, 'Color', colors(1,:), 'LineWidth', 1);hold on;
plot(etahats,'Color', colors(3,:), 'LineWidth', 1)
axis auto; a = axis; a(1) = 1; a(2) = Kmax; axis(a);
xlabel('k','FontSize',18); ylabel('eta','FontSize',18);
legend('eta','etahat')
set(gca, 'FontSize', 18);
print -depsc etas.eps
% eta.
figure(6);
plot(deltas, 'Color', colors(1,:), 'LineWidth', 1);hold on;
plot(deltahats,'Color', colors(3,:), 'LineWidth', 1)
axis auto; a = axis; a(1) = 1; a(2) = Kmax; axis(a);
xlabel('k','FontSize',18); ylabel('delta','FontSize',18);
legend('delta','deltahat')
set(gca, 'FontSize', 18);
print -depsc deltas.eps
% Tcost.
figure(7);
plot(tcosts, 'Color', colors(1,:), 'LineWidth', 1)
axis auto; a = axis; a(1) = 1; a(2) = Kmax; axis(a);
ylabel('T','FontSize',18); xlabel('k','FontSize',18);
set(gca, 'FontSize', 18);
print -depsc Tcosts.eps
