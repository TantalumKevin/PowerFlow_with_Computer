% 平衡、PV、PQ
    P = [ 0  0.7 -2.8];
    Q = [ 0    0 -1.2];
    U = [ 1 1.05    1];
delta = [ 0    0    0];
Y = -[20j -10j -10j;
    -10j 20j -10j;
    -10j -10j 20j];
G = real(Y);
B = imag(Y);


for it = 1:10
    deltaN = zeros(3);
    for i = 1:3
        for j = 1:3
            deltaN(i,j) = delta(i) - delta(j);
        end
    end
    deltaP = zeros(3,1);
    deltaQ = zeros(3,1);
    for i = 1:3
        deltaP(i) = P(i) - U(i) * sum(U .* ((G(i,:) .* cos(deltaN(i,:)) + (B(i,:) .* sin(deltaN(i,:))))));
        deltaQ(i) = Q(i) - U(i) * sum(U .* ((G(i,:) .* sin(deltaN(i,:)) - (B(i,:) .* cos(deltaN(i,:))))));
    end
    if max(abs([deltaP deltaQ])) < 1 
        break
    end
    H = zeros(3);
    N = zeros(3);
    M = zeros(3);
    L = zeros(3);

    for i = 1:3
        for j = 1:3
            if i ~= j
                H(i,j) = U(i)*U(j)*(G(i,j) * sin(deltaN(i,j)) - B(i,j) * cos(deltaN(i,j)));
                N(i,j) = U(i)*U(j)*(G(i,j) * cos(deltaN(i,j)) + B(i,j) * sin(deltaN(i,j)));
                M(i,j) = -N(i,j);
                L(i,j) =  H(i,j);
            else
                H(i,j) = -U(i)*U(j)*B(i,j)-Q(i);
                N(i,j) = +U(i)*U(j)*G(i,j)+P(i);
                M(i,j) = -U(i)*U(j)*G(i,j)+P(i);
                L(i,j) = -U(i)*U(j)*B(i,j)+Q(i);
            end
        end
    end
    temp1 = [H N];
    temp2 = [M L];
    J = [temp1(2:3,2:3) temp1(2:3,6);temp2(3,2:3) temp2(3,6)];
    correctionPQ = [P(1,2:3)';Q(3)];
    correctionU = J \ correctionPQ;
    U = [1 1.05 U(3)*(1+correctionU(3))];
    delta = [0 delta(2)+correctionU(1) delta(3)+correctionU(2)];
end