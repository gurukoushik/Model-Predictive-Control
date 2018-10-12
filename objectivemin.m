function[u]=objectivemin(A,B,C,Hc,Hp,X0,yref,weights)
currdir = [pwd filesep];
filename = [currdir , 'objmin.m'];
%Calculating sizes of the inputs
m=size(A,1);
n=size(A,2);
p=size(B,1);
q=size(B,2);
r=size(C,1);
s=size(C,2);
%Initializing required variables
Xnew=zeros(m,1);
Y=zeros(r,1);
Xold=X0;
obj=0;
U = sym('u', [Hc*q, 1]);
%Calculating the objective fuction by summing them over each step for the prediction horizon. 
for i=1:Hp
if i<=Hc
Xnew=A*Xold+B*U(i*q-q+1:i*q,:);
Y=C*Xnew;
Xold=Xnew;
phi=10*i*diag(weights);
obj=obj+(transpose(Y-yref))*phi*((Y-yref));
end
%The inputs are the same after the control horizon.
if i<=Hp&&i>=Hc
Xnew=A*Xold+B*U(Hc*q-q+1:Hc*q,:);
Y=C*Xnew;
Xold=Xnew;
phi=10*i*diag(weights);
obj=obj+(transpose(Y-yref))*phi*((Y-yref));
end
end
Constraint(); %Calling the constraint function to calculate the constraint matrices to pass to fmincon.
matlabFunction(obj,'file',filename,'vars',{U});
start=zeros(Hc*q,1); %initial values for inputs
%Optimising the objective with fmincon 
%uval= fmincon(@objmin,start,-ConstraintA,ConstraintB,ConstraintA,TerminalB,zeros(Hc*q,1),1000*ones(Hc*q,1));
uval= fmincon(@objmin,start,[],[],[],[],zeros(Hc*q,1),1000*ones(Hc*q,1));
u=uval(1:q,:);
end