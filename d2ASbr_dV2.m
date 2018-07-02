function [Haa, Hav, Hva, Hvv] = d2ASbr_dV2(dSbr_dVa, dSbr_dVm, Sbr, Cbr, Ybr, V, lam)

nl1 = length(lam);
nb = length(V);

diaglam = sparse(1:nl1, 1:nl1, lam, nl1, nl1);
diagSbr_conj = sparse(1:nl1, 1:nl1, conj(Sbr), nl1, nl1);

lam = diagSbr_conj*lam;
nl = length(lam);

diaglam1 = sparse(1:nl, 1:nl, lam, nl, nl);
diagV   = sparse(1:nb, 1:nb, V, nb, nb);


A = Ybr' * diaglam1 * Cbr;
B = conj(diagV) * A * diagV;
D = sparse(1:nb, 1:nb, (A*V) .* conj(V), nb, nb);
E = sparse(1:nb, 1:nb, (A.'*conj(V)) .* V, nb, nb);
F = B + B.';
G = sparse(1:nb, 1:nb, ones(nb, 1)./abs(V), nb, nb);

Saa = F - D - E;
Sva = 1j * G * (B - B.' - D + E);
Sav = Sva.';
Svv = G * F * G;

Haa = 2 * real( Saa + dSbr_dVa.' * diaglam * conj(dSbr_dVa) );
Hva = 2 * real( Sva + dSbr_dVm.' * diaglam * conj(dSbr_dVa) );
Hav = 2 * real( Sav + dSbr_dVa.' * diaglam * conj(dSbr_dVm) );
Hvv = 2 * real( Svv + dSbr_dVm.' * diaglam * conj(dSbr_dVm) );
