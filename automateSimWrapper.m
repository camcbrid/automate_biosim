clc

%repression cascade with ribosomes multi-equilibria
tic
A1 = [0 -1;
    0 0;
    -1 0];
u1 = 90;
pfun1 = @paramsRepCasc;
disp('fmincon')
out1fmincon = runfminconbio(A1,pfun1,u1);
disp('ODE_static')
out1static = runODE_static(A1,pfun1,u1);
disp('ODE_sweep')
out1sweep = runODE_sweep(A1,pfun1,[86,94]);
disp('fsolve')
out1fsolve = runfsolvebio(A1,pfun1,u1);
disp('ODE_mulambda')
out1mu = runODE_mulambda(A1,pfun1,u1);
toc
pause;

%activation cascade with protease multi-equilibria
tic
A2 = [0 1;
    0 0;
    1 0];
u2 = 0.07;
pfun2 = @paramsActCasc;
disp('fmincon')
out2fmincon = runfminconbio(A2,pfun2,u2);
disp('ODE_static')
out2static = runODE_static(A2,pfun2,u2);
disp('ODE_sweep')
out2sweep = runODE_sweep(A2,pfun2,[0.015,0.1]);
disp('fsolve')
out2fsolve = runfsolvebio(A2,pfun2,u2);
disp('ODE_mulambda')
out2mu = runODE_mulambda(A2,pfun2,u2);
toc
pause;

%oscillations with hasty oscillator
tic
A3 = [1 1;
    -1 -1];
pfun3 = @paramsARosc;
u3 = 0;
disp('fmincon')
out3fmincon = runfminconbio(A3,pfun3,u3);
disp('ODE_static')
out3static = runODE_static(A3,pfun3,u3);
disp('ODE_sweep')
out3sweep = runODE_sweep(A3,pfun3,[0,1]);
disp('fsolve')
out3fsolve = runfsolvebio(A3,pfun3,u3);
disp('ODE_mulambda')
out3mu = runODE_mulambda(A3,pfun3,u3);
toc
pause;

