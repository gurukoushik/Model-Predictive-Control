m=size(A,1);
n=size(A,2);
p=size(B,1);
q=size(B,2);
r=size(C,1);
s=size(C,2);
%Constraint matrix A to set up the constraint in the form AX<=B
ConstraintA=zeros(r*Hp,q*Hc);
for l=1:Hc
    for k=1:Hc   
    ConstraintA(k*r-r+1:k*r,l*q-q+1:l*q)=C*(A^(k-l-1))*B;
    end
end
for l=1:Hc-1
    for k=Hc:Hp   
    ConstraintA(k*r-r+1:k*r,l*q-q+1:l*q)=C*(A^(k-l-1))*B;
    end
end
sum=C*B;
for i=Hc+1:Hp
    sum=sum+C*A^(i-Hc-1)*B;
    ConstraintA(i*r-r+1:i*r,Hc*q-q+1:Hc*q)=sum;
end

for l=1:Hc
    for k=1:Hp
    if l>k 
    ConstraintA(k*r-r+1:k*r,l*q-q+1:l*q)=zeros(r,q);    
    end
    end
end
%Constraint matrix B to set up the constraint in the form AX<=B
ConstraintB=zeros(r*Hp,1);
for w=1:Hp
    ConstraintB(w*r-r+1:r*w,:)=C*(A^(w))*X0;
end
%Constraint matrix B to set up the constraint in the form AX<=B
TermB=zeros(r*Hp,1);
TermB(r*Hp-r+1:r*Hp,:)=yref;
TerminalB=TermB-ConstraintB;