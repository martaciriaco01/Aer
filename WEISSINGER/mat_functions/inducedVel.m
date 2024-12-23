function u = inducedVel(induced, inductor, finite)
%  INDUCEDVEL - Ballielo, Ciriaco
%  Restituisce la Velocità Indotta dal Vortice a Ferro
%  di Cavallo di un Pannello sul Punto di Controllo di
%  un Altro Dato Pannello (Per Unità di Circolazione)
%  Applicando la Legge di Biot-Savart.
%  Si Tiene Conto del Vortice Finito e dei Due Semi-Infiniti,
%  a meno che sia Specificato di Non Considerare quello Finito.
%  
%  Syntax
%    u = INDUCEDVEL(control, inductor)
%    u = INDUCEDVEL(induced, inductor, false)
%
%  Input Arguments
%    induced - Struttura del Pannello Indotto
%    inductor - Struttura del Pannello Induttore
%    finite - Booleano per Considerare il Vortice Finito
%             true | false (Default: true)
%
%  Output Arguments
%    u - Vettore Velocità Indotta sul Punto di Controllo

if(nargin == 2)
    finite = true;
end

infinite = [1e+4; 0; 0]; % Approssimazione Lunghezza Vortice Semi-Infinito

u = [0; 0; 0];

% Contributo Primo Vortice Semi-Infinito
r0 = -infinite;
r1 = induced.xc - (inductor.xv1 + infinite);
r2 = induced.xc - inductor.xv1;
r1xr2 = cross(r1, r2);
u = u + 1 / (4 * pi) * dot(r0, r1 / norm(r1) - r2 / norm(r2)) * r1xr2 / norm(r1xr2)^2;

% Contributo Secondo Vortice Semi-Infinito
r0 = infinite;
r1 = induced.xc - inductor.xv2;
r2 = induced.xc - (inductor.xv2 + infinite);
r1xr2 = cross(r1, r2);
u = u + 1 / (4 * pi) * dot(r0, r1 / norm(r1) - r2 / norm(r2)) * r1xr2 / norm(r1xr2)^2;

% Contributo Vortice Portante
if(finite)
    r0 = inductor.xv2 - inductor.xv1;
    r1 = induced.xc - inductor.xv1;
    r2 = induced.xc - inductor.xv2;
    r1xr2 = cross(r1, r2);
    u = u + 1 / (4 * pi) * dot(r0, r1 / norm(r1) - r2 / norm(r2)) * r1xr2 / norm(r1xr2)^2;
end

end