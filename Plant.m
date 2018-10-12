function [yret]=Plant(i,Y,T,U)
    u=U;
    bias=[-0.00001 0.00004 0.03];
    tspan=[T(i) T(i+1)]; %setting up time span
    [t,y] = ode23s(@(t,y) fcc_fn_to_solve_odemodel(t,y,u), tspan, Y); %Using ode23s to solve for state values for the given input
    yret=y(size(y,1),:); %Isolating the last computation of the numerical solver.
end
