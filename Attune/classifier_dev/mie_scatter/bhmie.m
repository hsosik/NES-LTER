function [S1, S2, Qext, Qsca] = bhmie(x, m, theta)
% function [S1, S2, Qext, Qsca] = bhmie(x, m, theta)
% Mie scattering amplitudes for a homogeneous sphere (Bohren & Huffman 1983,
% BHMIE), vectorized over scattering angles.
% x     - size parameter, pi d n_medium / lambda
% m     - complex refractive index relative to the medium (n + i k, k >= 0)
% theta - scattering angles (radians), any shape
% S1,S2 - complex amplitudes (perpendicular, parallel), size of theta;
%         unpolarized intensity ~ (|S1|^2 + |S2|^2) / 2
%
% Nov 2026, classifier_dev/mie_scatter

mu = cos(theta(:))';
nstop = ceil(x + 4 * x^(1/3) + 2);
y = m * x;
nmx = ceil(max(nstop, abs(y)) + 15);
% logarithmic derivative D_n(mx), downward recurrence
D = zeros(nmx, 1);
for n = nmx - 1:-1:1
    D(n) = (n + 1) / y - 1 / (D(n + 1) + (n + 1) / y);
end
% Riccati-Bessel functions psi_n(x), chi_n(x); start n = -1, 0
psi0 = cos(x); psi1 = sin(x);
chi0 = -sin(x); chi1 = cos(x);
xi1 = complex(psi1, -chi1);
pim1 = zeros(size(mu)); pin = ones(size(mu)); % pi_0, pi_1
S1 = zeros(size(mu)); S2 = zeros(size(mu));
Qext = 0; Qsca = 0;
for n = 1:nstop
    fn = (2 * n + 1) / (n * (n + 1));
    psi = (2 * n - 1) * psi1 / x - psi0;
    chi = (2 * n - 1) * chi1 / x - chi0;
    xi = complex(psi, -chi);
    a = ((D(n) / m + n / x) * psi - psi1) / ((D(n) / m + n / x) * xi - xi1);
    b = ((D(n) * m + n / x) * psi - psi1) / ((D(n) * m + n / x) * xi - xi1);
    Qsca = Qsca + (2 * n + 1) * (abs(a)^2 + abs(b)^2);
    Qext = Qext + (2 * n + 1) * real(a + b);
    taun = n * mu .* pin - (n + 1) * pim1;
    S1 = S1 + fn * (a * pin + b * taun);
    S2 = S2 + fn * (a * taun + b * pin);
    pnext = ((2 * n + 1) / n) * mu .* pin - ((n + 1) / n) * pim1;
    pim1 = pin; pin = pnext;
    psi0 = psi1; psi1 = psi; chi0 = chi1; chi1 = chi; xi1 = complex(psi1, -chi1);
end
Qsca = 2 / x^2 * Qsca; Qext = 2 / x^2 * Qext;
S1 = reshape(S1, size(theta)); S2 = reshape(S2, size(theta));
end
