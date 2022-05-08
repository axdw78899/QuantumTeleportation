% Quantum Teleportation 
clc;
clear all;
Id = [1 0; 0 1];
ket0 = [1; 0];
ket1 = [0; 1];
CNOT = [1 0 0 0; 0 1 0 0; 0 0 0 1; 0 0 1 0];
Hadamard = (1/sqrt(2))*[1 1; 1 -1];
%Qubit To Teleport
r=rand(1,2);
r=r/norm(r);
alfa=r(1,1);
beta=r(1,2);
qBitToTeleport = alfa * ket0 + beta * ket1;
qBitToTeleportDensityMatrix = qBitToTeleport * qBitToTeleport';

disp('Input');
disp(qBitToTeleportDensityMatrix);

p=0.01;
M1 = sqrt(1-p)*[1 0; 0 1];
M2 = sqrt(p)*[0 1; 1 0];
%[0 1; 1 0]
%[1 0; 0 -1]
%[0 -1i; 1i 0]
maxEntangledState = (kron(ket0,ket0)+kron(ket1,ket1))/sqrt(2);

% Density Matrix
system = kron(qBitToTeleport,maxEntangledState);
systemDensityMatrix = system *system';

% Applying noise
K=0;
DensityMatrix=0;
[out] = bin_listq1(3);
for i=1:length(out)
K=getError(out(i,:),M1,M2);
DensityMatrix = DensityMatrix + K * systemDensityMatrix * K';
end

% Applying CNOT Gate
CNOT = kron(CNOT, Id);
DensityMatrix = CNOT * DensityMatrix * CNOT';

% Applying noise
K=0;
DensityMatrix2=0;
[out] = bin_listq1(3);
for i=1:length(out)
K=getError(out(i,:),M1,M2);
DensityMatrix2 = DensityMatrix2 + K * DensityMatrix * K';
end

H = kron(kron(Hadamard,Id),Id);
DensityMatrix2 = H * DensityMatrix2 * H';

% Applying noise
K=0;
DensityMatrix3=0;
[out] = bin_listq1(3);
for i=1:length(out)
K=getError(out(i,:),M1,M2);
DensityMatrix3 = DensityMatrix3 + K * DensityMatrix2 * K';
end
%measure
[OutcomeDensityMatrix1,prob1,outcome1] = measure(DensityMatrix3, [1 0 0],ket0);
[OutcomeDensityMatrix2,prob2,outcome2] = measure(OutcomeDensityMatrix1, [0 1 0],ket0);

OutcomeDensityMatrix2(find(OutcomeDensityMatrix2==0))=[];
Output = reshape(OutcomeDensityMatrix2,[2,2]);
outcome=[outcome1 outcome2];

disp('outcome');
disp(outcome);

if(outcome==[1 1])
    Output=[0 -1i; 1i 0]*Output*[0 -1i; 1i 0];
elseif(outcome==[0 1])
    Output=[0 1; 1 0]*Output*[0 1; 1 0];
elseif(outcome==[1 0])
    Output=[1 0; 0 -1]*Output*[1 0; 0 -1];
else
    Output=Output;
end
Fidelity=trace(Output*qBitToTeleportDensityMatrix);
disp('Output');
disp(Output);
disp('Fidelity');
disp(Fidelity);
