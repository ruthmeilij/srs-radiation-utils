function EspacioFase()

% Tipo de particula (7: fotones, 3: electrones)
tipo_part = 7;
% Cantidad de particulas en cada haz (en cada circulo)
ParticleCount = 1000;
% Tamaño de campo (radio de los circulos en cm)
size_campo = 0.4;
% Divergencia (radio en porcentaje)
BeamTargetScaleFactor = 2;
% Desplazamiento del PTV en cm
displacement_vec = [-2.33,1.14,-10.44];

%% Asignacion de energia cinetica inicial 

% Asignar un tipo de energia 
metodoEnergia = 2;

% Caso 1: Beam monoenergetico

if metodoEnergia == 1
    monoE = 6;                                                              
end

% Caso 2: Beam con espectro de energias segun file externo

if metodoEnergia == 2
    % leer file externo (2da col: energia, 3ra col: prob)
    spect = dlmread('e_linac.dat');

    E = spect(:,2);
    %cambio de unidades para que sea en MeV (en caso de ser necesario)
    %E = E * 1E-6;
    
    % Probabilidad
    %P = spect(:,3);
    
    % Probabilidad normalizada en caso de que no lo este
    P = spect(:,3) / sum(spect(:,3));
end


%%  Aca se arman los arcos

% Nota sobre las especificaciones: 

% 'radius' es el radio del arco (en cm)
% 'angulo_inicial' es el angulo inicial del arco 
% 'angulo_final' es el angulo final del arco 
% 'equisp_beam' es el equiespaciamiento de los haces 
% 'archangle' es la posicion angular del arco en el espacio
% 'peso_por_arco' es el peso estadistico del arco entero
% 'peso_por_haz' es el peso estadistico de cada haz dentro del arco
% La cantidad de lineas será la cantidad de arcos totales


especificaciones = [
    struct('radius', 50, 'angulo_inicial', 230, 'angulo_final', 250, 'equisp_beam', 10, 'archangle', 0, 'peso_por_arco', 2, 'peso_por_haz', [])
    struct('radius', 50, 'angulo_inicial', 50, 'angulo_final', 70, 'equisp_beam', 10, 'archangle', 0, 'peso_por_arco', 2, 'peso_por_haz', [])
    struct('radius', 50, 'angulo_inicial', 90, 'angulo_final', 240, 'equisp_beam', 10, 'archangle', 30, 'peso_por_arco', 1, 'peso_por_haz', [])
    struct('radius', 50, 'angulo_inicial', 110, 'angulo_final', 200, 'equisp_beam', 10, 'archangle', 50, 'peso_por_arco', 1, 'peso_por_haz', [])
];

% Cantidad de arcos
ArchCount = max(size(especificaciones));

for i = 1 : ArchCount
    arco = especificaciones(i);
    numero_de_haces = (arco.angulo_final - arco.angulo_inicial) / arco.equisp_beam + 1;
    dim_peso_por_haz = size(arco.peso_por_haz);
    assert(isempty(arco.peso_por_haz) || (dim_peso_por_haz(2) == numero_de_haces), 'Faltan especificar pesos %d estadisticos para cada haz en el arco %d', numero_de_haces - dim_peso_por_haz(2), i);
end

BeamCount = zeros(ArchCount,1);
for i = 1 : ArchCount
    arco = especificaciones(i);
    % Cantidad de haces a lo largo de cada arco
    BeamCount(i) = BeamCount(i) + ((arco.angulo_final - arco.angulo_inicial)/ arco.equisp_beam) + 1;
end

% Formo los arcos definiendo los centros de los haces
BeamCount_size = size(BeamCount);
BeamPositionsList = zeros(BeamCount_size(1),BeamCount_size(2), 3);
for i = 1 : ArchCount
    arco = especificaciones(i);
    for CurrentBeamIndex=1: BeamCount(i)
        paso = (deg2rad(arco.equisp_beam) * (CurrentBeamIndex-1));
        
        % el arco se forma en el plano XY
        x = arco.radius * cos(deg2rad(arco.angulo_inicial) + paso);% * sin(pi/2);
        y = arco.radius * sin(deg2rad(arco.angulo_inicial) + paso);% * sin(pi/2);
        z = 0; %r * cos(pi/2);
              
        BeamPositionsList(i, CurrentBeamIndex, :) = [x, y, z];
    end
end

%% Armo la matriz completa 

% Row definitions
ParticleType = 1;
BeamEnergy = 2;
PesoEstadistico = 9;

OutputMatrix = zeros(sum(BeamCount) * ParticleCount, 9);

% PRIMERA COLUMNA DEL ESPACIO DE FASES -> Tipo de particula 
OutputMatrix(:,ParticleType) = tipo_part;

OutputWriteIndex = 1;

for CurrentArch = 1 : ArchCount
    arco = especificaciones(CurrentArch);
    
    for CurrentBeamIndex = 1 : BeamCount(CurrentArch)
        BeamPosition = squeeze(BeamPositionsList(CurrentArch,CurrentBeamIndex, :));
        BeamPositionNormalized = BeamPosition/arco.radius;
        %BeamPositionNormalized = BeamPosition/norm(BeamPosition);
        
        % Genero puntos random acotados a un circulo
        RandomCircle = generateRandomCircle(size_campo, BeamPositionNormalized, ParticleCount);
        % Defino el haz (Beam) y su respectiva seccion eficaz (BeamTarget)
        BeamTarget = RandomCircle * BeamTargetScaleFactor;
        Beam = RandomCircle + BeamPosition';      
        
        rot = eul2rotm([0 -deg2rad(arco.archangle) 0], 'XYZ');
        for rowI = 1 : size(Beam, 1)
            Beam(rowI, :) = Beam(rowI, :) * rot;
            BeamTarget(rowI, :) = BeamTarget(rowI, :) * rot;
        end

        if isempty(arco.peso_por_haz) 
            peso_por_haz = 1;
        else
            peso_por_haz = arco.peso_por_haz(CurrentBeamIndex);
        end
        
        for ParticleIndex = 1: ParticleCount
                      
            switch metodoEnergia
                % SEGUNDA COLUMNA DEL ESPACIO DE FASES -> Energia elegida
                case 1
                    OutputMatrix(OutputWriteIndex,BeamEnergy) = monoE;
                case 2
                    OutputMatrix(OutputWriteIndex,BeamEnergy) = rowEnergy(E,P);
            end
            
            % Defino las posiciones iniciales de las particulas
            ParticlePosition = Beam(ParticleIndex,:) + displacement_vec;
            % TERCERA, CUARTA Y QUINTA COLUMNA DEL ESPACIO DE FASES
            OutputMatrix(OutputWriteIndex,3:5) = ParticlePosition;
            
            % Defino los cosenos directores / direcciones iniciales
            TargetVector = ParticlePosition - (BeamTarget(ParticleIndex, :) + displacement_vec);
            DCosine = TargetVector / norm(TargetVector);
            % SEXTA, SEPTIMA Y OCTAVA COLUMNA DEL ESPACIO DE FASES
            OutputMatrix(OutputWriteIndex, 6:8) = -DCosine;
            
            OutputMatrix(OutputWriteIndex, PesoEstadistico) = arco.peso_por_arco * peso_por_haz;
                      
            OutputWriteIndex = OutputWriteIndex + 1;
        end
    end
end

%% Guardo la matriz
dlmwrite('esp_fase.dat', OutputMatrix, '\t')


%% VISUALIZAR

hold on

%posiciones iniciales

x0 = OutputMatrix(:,3);
y0 = OutputMatrix(:,4);
z0 = OutputMatrix(:,5);
% 
plot3(x0, y0, z0, '.', 0,0,0,'rx');
grid
xlabel('x [cm]','FontSize',10);
ylabel('y [cm]','FontSize',10);
zlabel('z [cm]','FontSize',10);
%xlim([-15 25]);
%ylim([-5 25]);
%zlim([-50 50]);

% direcciones iniciales

% dirx0 = OutputMatrix(:,6);
% diry0 = OutputMatrix(:,7);
% dirz0 = OutputMatrix(:,8);
% 
% quiver3(x0, y0, z0, dirx0, diry0, dirz0, 20)
% xlabel('x','FontSize',20);
% ylabel('y','FontSize',20);
% zlabel('z','FontSize',20);


end

function r = generateRandomCircle(radius, pointTo, N)

  r = radius * sqrt(rand(1, N));
  th = rand(1, N) * 2 * pi;
  x = r .* cos(th);
  y = r .* sin(th);
  z = zeros(N, 1);

  % If the spherical cap is centered around the north pole, we're done.
  if all(pointTo(:) == [0;0;1])
      r = [x', y', z];
      return;
  end

  % Find the rotation axis `u` and rotation angle `rot` [1]
  u = normc(cross([0;0;1], normc(pointTo)));
  rot = acos(dot(normc(pointTo), [0;0;1]));

  % Convert rotation axis and angle to 3x3 rotation matrix [2]
  % [2] See https://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle
  crossMatrix = @(x,y,z) [0 -z y; z 0 -x; -y x 0];
  R = cos(rot) * eye(3) + sin(rot) * crossMatrix(u(1), u(2), u(3)) + (1-cos(rot))*(u * u');

  % Rotate [x; y; z] from north pole to `pointTo`.
  r = R * [x; y; z'];
  r = r';
end
function y = normc(x)
y = bsxfun(@rdivide, x, sqrt(sum(x.^2)));
end

function e = rowEnergy(E, P)
    x = rand;
    acc = P(1);
    i = 1;

    while x > acc
        i = i + 1;
        acc = acc + P(i);
    end
            
    e = E(i);
end