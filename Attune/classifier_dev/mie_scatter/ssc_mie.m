function I = ssc_mie(d, m, lambda, nmed, halfangle)
% function I = ssc_mie(d, m, lambda, nmed, halfangle)
% Side scatter of spheres collected over a cone around 90 deg (axis
% perpendicular to the beam), unpolarized: integral over the cone of
% (|S1|^2 + |S2|^2) / 2 d(solid angle). Arbitrary units, same for all
% particles at one wavelength (divide by a bead to get bead units).
% d         - diameters (um), vector
% m         - relative complex refractive index (scalar, or one per d)
% lambda    - vacuum wavelength (um)
% nmed      - medium refractive index
% halfangle - cone half-angle (deg)
%
% Nov 2026, classifier_dev/mie_scatter

th = linspace(0, pi, 3601); % 0.05 deg grid for |S|^2 (large spheres oscillate)
% cone directions v = (cos b, sin b cos p, sin b sin p) around x; beam along z
b = linspace(0, deg2rad(halfangle), 60); p = linspace(0, 2 * pi, 121); p(end) = [];
[B, P] = meshgrid(b, p);
w = sin(B) * (b(2) - b(1)) * (p(2) - p(1));
thc = acos(sin(B) .* sin(P)); % scattering angle of each cone direction
if isscalar(m), m = repmat(m, size(d)); end
I = zeros(size(d));
for i = 1:numel(d)
    x = pi * d(i) * nmed / lambda;
    [S1, S2] = bhmie(x, m(i), th);
    i0 = (abs(S1).^2 + abs(S2).^2) / 2;
    I(i) = sum(w(:) .* interp1(th, i0, thc(:)));
end
end
