function [Gaa, Gav, Gva, Gvv] = d2Sbus_dV2(Ybus, V, lam)

n = length(V);
Ibus    = Ybus * V;
diaglam = sparse(1:n, 1:n, lam, n, n);
diagV   = sparse(1:n, 1:n, V, n, n);

A = sparse(1:n, 1:n, lam .* V, n, n);
B = Ybus * diagV;
C = A * conj(B);
D = Ybus' * diagV;
E = conj(diagV) * (D * diaglam - sparse(1:n, 1:n, D*lam, n, n));
F = C - A * sparse(1:n, 1:n, conj(Ibus), n, n);
G = sparse(1:n, 1:n, ones(n, 1)./abs(V), n, n);

Gaa = E + F;
Gva = 1j * G * (E - F);
Gav = Gva.';
Gvv = G * (C + C.') * G;
