clear;

tolerance = 10^(-6);

bus = importdata('bus.dat');
branch = importdata('branch.dat');

N = length(bus);
L = length(branch);

fromBus = branch(:,2);
toBus = branch(:,3);
refBus = find(bus(:,2) == 3);

A = zeros(L,N);
for l = 1:L
    A(l,fromBus(l)) = -1;
    A(l,toBus(l)) = 1;
end
A(:,refBus) = [];

bmva = 100;
deg2rad = pi/180;

Bhat = diag(bmva*deg2rad ./branch(:,6));
B = -A' * Bhat * A;
LSF = Bhat * A * inv(B);
for i = 1:12
    [i, length(find(abs(LSF) < 10^(-i)))]
end
LSF(find(abs(LSF) < tolerance)) = 0;
LSF = [LSF(:,1:refBus-1),zeros(L,1),LSF(:,refBus:end)];

dlmwrite('shift_factor.dat', LSF, 'delimiter', '\t', 'precision', 6);


