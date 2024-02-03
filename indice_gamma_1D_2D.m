clear all;
close all;
clc;

%% Lectura de archivos 

INPUT_1 = dicomread('dosis_1er_arco.dat');
INPUT_2 = dlmread('dosis_2do_arco.dat');


%% Procesamiento de datos (matrices cuadradas y misma dimension)
% 
Dim = size(INPUT_2); 
Lado = sqrt(max(Dim));
for i=1:Lado
    for j=1:Lado
        indice = Lado*(i-1) + j;
        DOSIS_1(Lado - i+1,j) = INPUT_1(indice,3);
        ErDOSIS_1(Lado - i+1,j) = INPUT_1(indice,4).*INPUT_1(indice,3)/100; 
        %
        DOSIS_2(Lado - i + 1,j) = INPUT_2(indice,3);
        ErDOSIS_2(Lado - i + 1,j) = INPUT_2(indice,4).*INPUT_2(indice,3)/100;
        %      
    end
end

%% Procesamiento de datos (matrices no cuadradas)

% si mi vector columna tiene 89262 elementos
% 
% DOSIS_2 = zeros(173,512);
% ErDOSIS_2 = zeros(173,512);
% for i=173:-1:1
%     for j=1:512
%         indice = 512*(i-1) + j;
%         DOSIS_2(i,j) = INPUT_2(indice,3);
%         ErDOSIS_2(i,j) = INPUT_2(indice,4).*INPUT_2(indice,3)/100;
%         %      
%     end
% end


%% Extraer perfiles (1D) (valor medio)

dim1 = size(DOSIS_1);
dim2 = size(DOSIS_2);

PERFIL_HOR_1 = DOSIS_1(:,round(dim1(1,1)/2));
%PERFIL_HOR_1_er = ErDOSIS_1(:,round(dim1(1,1)/2));

PERFIL_HOR_2 = DOSIS_2(:,round(dim2(1,1)/2));
PERFIL_HOR_2_er = ErDOSIS_2(:,round(dim2(1,1)/2));


%% Normalizo al maximo
% Ademas, establezco datos de referencia y datos de prueba

vec_ref_norm = 100*PERFIL_HOR_1/(max(max(PERFIL_HOR_1)));
vec_test_norm = 100*PERFIL_HOR_2/(max(max(PERFIL_HOR_2)));

test_error = PERFIL_HOR_2_er;

matrix_ref_norm = 100*DOSIS_1/(max(max(DOSIS_1)));
matrix_test_norm = 100*DOSIS_2/(max(max(DOSIS_2)));

m_test_error = ErDOSIS_2;


%% Lectura de tolerancias 

TolDosis = 0.02; 
TolDistancia = 2.062; % 2 mm / si 1 pix = 0.97 mm => 2 mm = 2.062 pix
%TolDistancia = 1.8; % con relación de aspecto 1:1 

window = 2;

%% UNA DIMENSIoN
% NOTA: considero la incertidumbre para DD pero no para DTA

IGamma = zeros(1,max(size(vec_ref_norm))) + realmax;


for x = 1:max(size(vec_ref_norm))    
    for y = max(1,x-window):min(x+window,max(size(vec_ref_norm)))
        %(DT(y) − DR(x))
        DD = vec_test_norm(y) - vec_ref_norm(x);

        %(y - x)
        DTA = y - x;
        %Gamma al cuadrado
        MaybeGamma = (DTA/TolDistancia)^2 + (DD/TolDosis)^2;
        if(MaybeGamma < IGamma(x))
            IGamma(x) = MaybeGamma;
        end
    end
    
end

%Termino de calcular gama haciendo la raiz
IGamma = sqrt(IGamma);

%% DOS DIMENSIONES
% NOTA: considero la incertidumbre para DD pero no para DTA

dim_ref = size(matrix_ref_norm);
dim_test = size(matrix_test_norm);

IGamma2 = zeros(dim_ref) + realmax;

for xx = 1:dim_ref(1,1)
    for xy = 1:dim_ref(1,2)
        for yx = max(1,xx-window):min(xx+window,dim_test(1,1))
            for yy = max(1,xy-window):min(xy+window,dim_test(1,2))
                %(DT(y) − DR(x))
                DD = matrix_test_norm(yx,yy) - matrix_ref_norm(xx,xy);
                if(abs(DD) <= m_test_error(yx,yy))
                    IGamma2(xx,xy) = 0;
                    break
                end
          
                %(y - x)
                H = (xx-yx)^2+(xy-yy)^2;
                %Gamma al cuadrado
                MaybeGamma2 = H/(TolDistancia^2) + (DD/TolDosis)^2;
                if(MaybeGamma2 < IGamma2(xx,xy))
                    IGamma2(xx,xy) = MaybeGamma2;
                end
            end
        end

    end

end

%Termino de calcular gamma haciendo la raiz
IGamma2 = sqrt(IGamma2);



