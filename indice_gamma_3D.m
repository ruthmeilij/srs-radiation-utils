clear all;
close all;
clc;


%% Lectura de archivos 

% Cortes Axiales
% cantidad total = 10

% file_list = [
%     ["Axial 9.9.dcm", "slice88.dat"];
%     ["Axial 10.03.dcm", "slice87.dat"];
%     ["Axial 10.15.dcm", "slice86.dat"];
%     ["Axial 10.28.dcm", "slice85.dat"];
%     ["Axial 10.40.dcm", "slice84.dat"];
%     ["Axial 10.53.dcm", "slice83.dat"];
%     ["Axial 10.65.dcm", "slice82.dat"];
%     ["Axial 10.78.dcm", "slice81.dat"];
%     ["Axial 10.90.dcm", "slice80.dat"];
%     ["Axial 11.03.dcm", "slice79.dat"];
% ];

% Cortes Frontales
% cantidad total = 16
% 
% file_list = [
%     ["Frontal -0.45.dcm", "frontal slice260.dat"];
%     ["Frontal -0.55.dcm", "frontal slice261.dat"];
%     ["Frontal -0.64.dcm", "frontal slice262.dat"];
%     ["Frontal -0.74.dcm", "frontal slice263.dat"];
%     ["Frontal -0.84.dcm", "frontal slice264.dat"];
%     ["Frontal -0.94.dcm", "frontal slice265.dat"];
%     ["Frontal -1.03.dcm", "frontal slice266.dat"];
%     ["Frontal -1.13.dcm", "frontal slice267.dat"];
%     ["Frontal -1.23.dcm", "frontal slice268.dat"];
%     ["Frontal -1.33.dcm", "frontal slice269.dat"];
%     ["Frontal -1.42.dcm", "frontal slice270.dat"];
%     ["Frontal -1.52.dcm", "frontal slice271.dat"];
%     ["Frontal -1.62.dcm", "frontal slice272.dat"];
%     ["Frontal -1.72.dcm", "frontal slice273.dat"];
%     ["Frontal -1.82.dcm", "frontal slice274.dat"];
%     ["Frontal -1.91.dcm", "frontal slice275.dat"];
%  ];

% Cortes Sagitales
% cantidad total = 14

file_list = [
    ["Sagital -3.dcm", "sagital slice225.dat"];
    ["Sagital -2.97.dcm", "sagital slice226.dat"];
    ["Sagital -2.88.dcm", "sagital slice227.dat"];
    ["Sagital -2.78.dcm", "sagital slice228.dat"];
    ["Sagital -2.68.dcm", "sagital slice229.dat"];
    ["Sagital -2.58.dcm", "sagital slice230.dat"];
    ["Sagital -2.48.dcm", "sagital slice231.dat"];
    ["Sagital -2.39.dcm", "sagital slice232.dat"];
    ["Sagital -2.29.dcm", "sagital slice233.dat"];
    ["Sagital -2.19.dcm", "sagital slice234.dat"];
    ["Sagital -2.09.dcm", "sagital slice235.dat"];
    ["Sagital -1.9.dcm", "sagital slice236.dat"];
    ["Sagital -1.8.dcm", "sagital slice237.dat"];
    ["Sagital -1.7.dcm", "sagital slice238.dat"];
];

tensor_ref_norm = [];
tensor_test_norm = [];

for file_index = 1:size(file_list, 1)
    
    files = file_list(file_index, :);
    DOSIS_1 = dicomread(files(1));
    INPUT_2 = dlmread(files(2));
    %DOSIS_2 = dicomread(files(2));

    %% Procesamiento de datos (matrices cuadradas y misma dimension)
    
%     Dim = size(INPUT_2); 
%     Lado = sqrt(max(Dim));
%     for i=1:Lado
%         for j=1:Lado
%             indice = Lado*(i-1) + j;
%             DOSIS_2(Lado - i + 1,j) = INPUT_2(indice,3);
%             %      
%         end
%     end

    
    %% Procesamiento de datos (matrices no cuadradas)
    
    DOSIS_2 = zeros(199,512);
    ErDOSIS_2 = zeros(199,512);
    for i=173:-1:1
        for j=1:512
            indice = 512*(i-1) + j;
            DOSIS_2(i,j) = INPUT_2(indice,3);
            ErDOSIS_2(i,j) = INPUT_2(indice,4).*INPUT_2(indice,3)/100;     
        end
    end
    
    
    dim1 = size(DOSIS_1);
    dim2 = size(DOSIS_2);

    DOSIS_2 = flip(DOSIS_2, 1);
    DOSIS_2 = flip(DOSIS_2, 2);
    
    ErDOSIS_2 = flip(DOSIS_2, 1);
    ErDOSIS_2 = flip(DOSIS_2, 2);

    % Normalizo al maximo
    % Establezco datos de referencia y datos de prueba

    tensor_ref_norm(:, :, file_index) = 100*DOSIS_1/(max(max(DOSIS_1)));
    tensor_test_norm(:, :, file_index) = 100*DOSIS_2/(max(max(DOSIS_2)));
    
    %tensor_test_norm_er(:, :, file_index) = 100*ErDOSIS_2/(max(max(ErDOSIS_2)));
    tensor_test_er(:, :, file_index) = ErDOSIS_2;

end

%% Lectura de tolerancias 

TolDosis = 2; 
%TolDistancia = 2.062; % 2 mm / si 0.97 mm = 1 pix => 2 mm = 2.062 pix
TolDistancia = 1.8; % con relación de aspecto 1:1 

window = 2;

%% INDICE GAMMA EN TRES DIMENSIONES 
%NOTA: Considero incertidumbre de los datos de prueba

dim_ref_tensor = size(tensor_ref_norm);
dim_test_tensor = size(tensor_test_norm);

IGamma3 = zeros(dim_ref_tensor) + realmax;

centro = dim_ref_tensor / 2;
radio_umbral = 10;

PuntosDescartadosDentroDeLaEsfera = 0;

for xx = 1:dim_ref_tensor(1,1)
    for xy = 1:dim_ref_tensor(1,2)
        for xz = 1:dim_ref_tensor(1,3)
            
            if (tensor_test_norm(xx,xy,xz) <= 5)
                if (norm(centro - [xx, xy, xz]) <= radio_umbral)
                    PuntosDescartadosDentroDeLaEsfera = PuntosDescartadosDentroDeLaEsfera + 1;
                end
                IGamma3(xx,xy,xz) = -1;
                continue
            end
            
%             if (norm(centro - [xx, xy, xz]) > radio_umbral)
%                 IGamma3(xx,xy,xz) = -1;
%                 continue
%             end

            for yx = max(1,xx-window):min(xx+window,dim_test_tensor(1,1))
                for yy = max(1,xy-window):min(xy+window,dim_test_tensor(1,2))
                    for yz = max(1,xz-window):min(xz+window,dim_test_tensor(1,3))

                        %(DT(y) − DR(x))
                        DD = tensor_test_norm(yx,yy,yz) - tensor_ref_norm(xx,xy,xz);
                        %if DD == 0
                        if (abs(DD) <= tensor_test_er(yx,yy,yz))
                            IGamma3(xx,xy,xz) = 0;
                            break
                        end  
                        %(y - x)
                        H = (xx-yx)^2+(xy-yy)^2+(xz-yz)^2;
                        %Gamma al cuadrado
                        MaybeGamma3 = H/(TolDistancia)^2 + ((DD)/(TolDosis))^2;
                        if(MaybeGamma3 < IGamma3(xx,xy,xz))
                            IGamma3(xx,xy,xz) = MaybeGamma3;
                        end
                    end
                end
            end
        end
    end
end

%Termino de calcular gamma haciendo la raiz
IGamma3 = sqrt(IGamma3);

%% Visualizacion

% -----------------------------------------
% FIGURA 1: INDICE GAMMA 3D

IGamma3_t = IGamma3(IGamma3 ~= sqrt(-1));
total = numel(IGamma3_t);
percent_accept = sum(IGamma3_t(:,:,:) <= 1, 'all');
percent = (100 * percent_accept) / total;

PuntosQuePasan = [];
PuntosQueNoPasan = [];
IGammaPuntos = [];
Puntos = [];
dimGamma3 = size(IGamma3);
for j=1:dimGamma3(1,3)
    [r, c] = find(IGamma3(:, :, j) <= 1 & IGamma3(:, :, j) ~= sqrt(-1));
    PuntosQuePasan = [PuntosQuePasan;[r , c, ones(length(r),1) * j]];
    [r, c] = find(IGamma3(:, :, j) > 1 & IGamma3(:, :, j) ~= sqrt(-1));
    PuntosQueNoPasan = [PuntosQueNoPasan; [r , c, ones(length(r),1) * j]];
    [r, c] = find(IGamma3(:, :, j) > 0 & IGamma3(:, :, j) ~= sqrt(-1));
    Puntos = [Puntos; [r , c, ones(length(r),1) * j]];
end
figure;
ax = pcshow(Puntos, IGamma3(IGamma3 > 0  & IGamma3 ~= sqrt(-1)));
title('Cálculo de Índice Gamma en 3D');
xlim([0 199]);
ylim([0 512]);
zlim([0 10]);
xlabel('z [pix]','FontSize',13);
ylabel('y [pix]','FontSize',13);
zlabel('cantidad de slice','FontSize',13);
daspect([20 20 1]);


% -----------------------------------------
%FIGURA 2: DISTRIBUCION DE DOSIS 3D

PuntosMatRef = zeros(2000000, 3);
DosisPuntosRef = zeros(2000000, 1);
npm = 1;
dimMat = size(tensor_ref_norm);
for j=1:dimMat(1,3)
    for k=1:dimMat(1,1)
        for l=1:dimMat(1,2)           
            if(tensor_ref_norm(k, l, j) == 0)
                continue;
            end
            PuntosMatRef(npm, :) = [k l j];
            DosisPuntosRef(npm) = tensor_ref_norm(k, l, j);
            npm = npm + 1;
        end
    end
end
figure;
ax = pcshow(PuntosMatRef(1:npm, :), DosisPuntosRef(1:npm, :));
xlim([0 199]);
ylim([0 512]);
zlim([0 10]);
title('Distribución de dosis de referencia');
xlabel('z [pix]','FontSize',13);
ylabel('y [pix]','FontSize',13);
zlabel('cantidad de slice','FontSize',13);
daspect([20 20 1]);

PuntosMatTest = zeros(2000000, 3);
DosisPuntosTest = zeros(2000000, 1);
npm = 1;
dimMat = size(tensor_test_norm);
for j=1:dimMat(1,3)
    for k=1:dimMat(1,1)
        for l=1:dimMat(1,2)           
            %if(matrix_test_norm(k, l, j) <= matrix_test_norm_er(k, l, j))
             if(tensor_test_norm(k, l, j) <= 2)
                continue;
            end
            PuntosMatTest(npm, :) = [k l j];
            DosisPuntosTest(npm) = tensor_test_norm(k, l, j);
            npm = npm + 1;
        end
    end
end
figure;
ax = pcshow(PuntosMatTest(1:npm, :), DosisPuntosTest(1:npm, :));
xlim([0 199]);
ylim([0 512]);
zlim([0 10]);
title('Distribución de dosis de prueba');
xlabel('z [pix]','FontSize',13);
ylabel('y [pix]','FontSize',13);
zlabel('cantidad de slice','FontSize',13);
daspect([20 20 1]);


% -----------------------------------------
% FIGURA 3: MAPA DE PUNTOS QUE PASAN Y NO PASAN EL CRITERIO
f = figure;
colormap(f, autumn)
ax = pcshow(PuntosQueNoPasan);
xlim([0 199]);
ylim([0 512]);
zlim([0 10]);
title('Mapa de puntos que no pasan el criterio');
xlabel('z [pix]','FontSize',13);
ylabel('y [pix]','FontSize',13);
zlabel('cantidad de slice','FontSize',13);
daspect([20 20 1]);

f = figure;
colormap(f, winter)
ax = pcshow(PuntosQuePasan);
xlim([0 199]);
ylim([0 512]);
zlim([0 10]);
annotation('textbox', [0.02, .06, 1, 0], 'string', sprintf('Porcentaje de aceptacion: %.3f%% [%d/%d]', percent, percent_accept, total), 'FitBoxToText','on', 'EdgeColor', [1 1 1], 'Color', [1 1 1])
title('Mapa de puntos que pasan el criterio');
xlabel('z [pix]','FontSize',13);
ylabel('y [pix]','FontSize',13);
zlabel('cantidad de slice','FontSize',13);
daspect([20 20 1]);




