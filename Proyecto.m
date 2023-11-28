clear all; close all;
%Seleccionar Imagen
[Nom_I,PathName1] = uigetfile('C:\Users\Joan\Documents\Proyecto_Vision\*.*','Seleccione el archivo imagen');
%Leer Imagen
I=imread(Nom_I);
%información de la imagen
InfoI=imfinfo(Nom_I); 
W_I=InfoI.Width ;
H_I=InfoI.Height ;
Form=InfoI.Format;
%conversión de RGB a GrayScale.
I= rgb2gray(I);

%ploteo de la información
figure(1);
subplot(2,2,1);imshow(I);title(['Image, W=' num2str(W_I) 'H=' num2str(H_I)]);
subplot(2,2,2);imhist(I);title('General Histogram');
SH = sum(I);
SV = sum(I');
subplot(2,2,3);plot(0:W_I-1,SH);xlabel('Horizontal Axe');ylabel('Cumulative sum of grays');title('Horizomtal profile');grid on;
subplot(2,2,4);plot(0:H_I-1,SV);xlabel('Vertical Axe');ylabel('Cumulative sum of grays');title('Vertical profile');grid on;

%filtros paso bajo

Mask = 131; %131
order = 2;
diffV = zeros(1,length(SV));
aux = SV;
for i=1:order
    init = aux(1);
    finit = aux(length(aux));
    for j=1:length(diffV)
        for k=(0:Mask-1)-(Mask-1)/2
            if j + k<=0
                diffV(j) = diffV(j) + init;
            else
                if j+k>length(aux)
                    diffV(j) = diffV(j) + finit;
                else
                    diffV(j) = diffV(j) + aux(j+k);
                end
            end
        end
        diffV(j) = diffV(j)/Mask;
    end
    aux = diffV;
end

Mask = 21; %131
order = 2;
aux = diffV;
for i=1:order
    init = aux(1);
    finit = aux(length(aux));
    for j=1:length(diffV)
        for k=(0:Mask-1)-(Mask-1)/2
            if j + k<=0
                diffV(j) = diffV(j) + init;
            else
                if j+k>length(aux)
                    diffV(j) = diffV(j) + finit;
                else
                    diffV(j) = diffV(j) + aux(j+k);
                end
            end
        end
        diffV(j) = diffV(j)/Mask;
    end
    aux = diffV;
end

Mask = 131; %131
order = 2;
diffH = zeros(1,length(SH));
aux = SH;
for i=1:order
    init = aux(1);
    finit = aux(length(aux));
    for j=1:length(diffH)
        for k=(0:Mask-1)-(Mask-1)/2
            if j + k<=0
                diffH(j) = diffH(j) + init;
            else
                if j+k>length(aux)
                    diffH(j) = diffH(j) + finit;
                else
                    diffH(j) = diffH(j) + aux(j+k);
                end
            end
        end
        diffH(j) = diffH(j)/Mask;
    end
    aux = diffH;
end

Mask = 21; %131
order = 2;
aux = diffH;
for i=1:order
    init = aux(1);
    finit = aux(length(aux));
    for j=1:length(diffH)
        for k=(0:Mask-1)-(Mask-1)/2
            if j + k<=0
                diffH(j) = diffH(j) + init;
            else
                if j+k>length(aux)
                    diffH(j) = diffH(j) + finit;
                else
                    diffH(j) = diffH(j) + aux(j+k);
                end
            end
        end
        diffH(j) = diffH(j)/Mask;
    end
    aux = diffH;
end
    
%derivo
order = 2;
aux = diffV;
for i=1:order
    init = aux(1);
    for j=1:length(diffV)
        if j-1==0
            diffV(j) = 0;
        else
            diffV(j) = aux(j) - aux(j-1);
        end
    end
    aux = diffV;
end

order = 2;
aux = diffH;
for i=1:order
    init = aux(1);
    for j=1:length(diffH)
        if j-1==0
            diffH(j) = 0;
        else
            diffH(j) = aux(j) - aux(j-1);
        end
    end
    aux = diffH;
end



%Hallo los picos de los perfiles que probablemente representan el borde de
%alguna de las modeas en los bordes de la región de importancia

[pcks,locs] = findpeaks((-1)*diffV);

porcentaje = 0.1;

umbral = porcentaje*max(pcks);

j0 = 1;
j1 = length(pcks);
aux = pcks(1);

while (aux<umbral)
    j0 = j0 + 1;
    aux = pcks(j0);
end
aux = pcks(length(pcks));
while (aux<umbral)
    j1 = j1-1;
    aux = pcks(j1);
end

j0 = locs(j0);
j1 = locs(j1);

[pcks,locs] = findpeaks((-1)*diffH);

umbral = porcentaje*max(pcks);

i0 = 1;
i1 = length(pcks);
aux = pcks(1);

while (aux<umbral)
    i0 = i0 + 1;
    aux = pcks(i0);
end
aux = pcks(length(pcks));
while (aux<umbral)
    i1 = i1-1;
    aux = pcks(i1);
end

i0 = locs(i0);
i1 = locs(i1);

%recorto la imagne y dejo en blanco la región de la figura original en la
%que se enconntraban las monedas

IR = I(j0:j1,i0:i1);

I(j0:j1,i0:i1) = 255;

[f,c] = size(IR);

figure(4);
subplot(1,2,1);imshow(I);subplot(1,2,2);imshow(IR);

%Hago un histograma acumulado de la imagen I

hist = histog(I);

hist(256) = hist(256) - f*c;

for i=1:length(hist)
    if i>1
        hist(i) = hist(i) + hist(i-1);
    end
end
%normalizo el histograma
hist = hist/hist(256);
figure(5);plot(1:length(hist),hist);

%Ahora hallo que niveles de grises representan al menos a lo mucho 80% de
%la imagen, para en la imagen recortada hacer que los pixeles con esos
%niveles de gris los transforme en blanco (255)

umbral = 0.2;
i=1;
while hist(i) <=umbral
    i = i + 1;
end

i = i - 1;

%Ahora hallo el histograma acumulado de la imagen recortada

hist = histog(IR);

for k=1:length(hist)
    if k>1
        hist(k) = hist(k) + hist(k-1);
    end
end
hist = hist/hist(256);

%hallo los niveles de gris que representan al menos 99.5% de la imagen

umbral = 0.005;
j=1;

while hist(j)<=umbral
    j = j + 1;
end
j = j - 2;

%utilizo el nivel de gris máximo y el nivel de gris minimo hallados
%anteriormente para mejorar el contraste de la imagen

m = 255/(i-j);

T = zeros(256,1);
for k=1:256
    T(k) = uint8((k-1-j)*m);
end

for i=1:f
    for j=1:c
        IR(i,j) = T(IR(i,j)+1);
    end
end

%ploteo la nueva imagen y hallo sus nuevos perfiles horizontal y vertical

figure(2);
subplot(2,2,1);imshow(IR);title(['Image, W=' num2str(c) 'H=' num2str(f)]);
subplot(2,2,2);imhist(IR);title('General Histogram');
SH = sum(IR);
SV = sum(IR');
subplot(2,2,3);plot(0:c-1,SH);xlabel('Horizontal Axe');ylabel('Cumulative sum of grays');title('Horizomtal profile');grid on;
subplot(2,2,4);plot(0:f-1,SV);xlabel('Vertical Axe');ylabel('Cumulative sum of grays');title('Vertical profile');grid on;

%alpico un filtro pasa bajas en ambos perfiles

Mask = 11; %131
order = 2;
diffV = zeros(1,length(SV));
aux = SV;
for i=1:order
    init = aux(1);
    finit = aux(length(aux));
    for j=1:length(diffV)
        for k=(0:Mask-1)-(Mask-1)/2
            if j + k<=0
                diffV(j) = diffV(j) + init;
            else
                if j+k>length(aux)
                    diffV(j) = diffV(j) + finit;
                else
                    diffV(j) = diffV(j) + aux(j+k);
                end
            end
        end
        diffV(j) = diffV(j)/Mask;
    end
    aux = diffV;
end

Mask = 11; %131
order = 2;
diffH = zeros(1,length(SH));
aux = SH;
for i=1:order
    init = aux(1);
    finit = aux(length(aux));
    for j=1:length(diffH)
        for k=(0:Mask-1)-(Mask-1)/2
            if j + k<=0
                diffH(j) = diffH(j) + init;
            else
                if j+k>length(aux)
                    diffH(j) = diffH(j) + finit;
                else
                    diffH(j) = diffH(j) + aux(j+k);
                end
            end
        end
        diffH(j) = diffH(j)/Mask;
    end
    aux = diffH;
end

figure(6);subplot(1,2,1);plot(1:length(diffH),diffH);
subplot(1,2,2);plot(1:length(diffV),diffV);

%defino el umbral para diferenciar el fondo de las monedas en el perfil
%horizontal utilizando los picos
umbral = min([diffH(1), diffH(length(diffH))]);

[pcks,locs] = findpeaks((-1)*diffH);

k=0;

for i=1:length(locs)
    auxm = locs(i);
    auxM = locs(i);
    while diffH(auxm)<umbral
        auxm = auxm - 1;
    end
    while diffH(auxM)<umbral
        auxM = auxM + 1;
    end
    if auxm ~= auxM
        if k>=1
            if HP(k)~=auxm
                k = k + 1;
                HP(k,1) = auxm;
                HP(k,2) = auxM;
            end
        else
            k = k + 1;
            HP(k,1) = auxm;
            HP(k,2) = auxM;
        end
    end
end

%Realizo el mismo proceso con el perfil vertical

umbral = min([diffV(1), diffV(length(diffV))]);

[pcks,locs] = findpeaks((-1)*diffV);

k=0;

for i=1:length(locs)
    auxm = locs(i);
    auxM = locs(i);
    while diffV(auxm)<umbral
        auxm = auxm - 1;
    end
    while diffV(auxM)<umbral
        auxM = auxM + 1;
    end
    if auxm ~= auxM
        if k>=1
            if VP(k)~=auxm
                k = k + 1;
                VP(k,1) = auxm;
                VP(k,2) = auxM;
            end
        else
            k = k + 1;
            VP(k,1) = auxm;
            VP(k,2) = auxM;
        end
    end
end
[fV,cV] = size(VP);
[fH,cH] = size(HP);
umbral = 200;
porcentaje = 0.05;
l=0;
%Analizamos las regiones halladas para saber si realmente hay monedas en
for i=1:fV
    for j=1:fH
        %auxI = filtro_mediana(IR(VP(i,1):VP(i,2),HP(j,1):HP(j,2)));
        %auxI = filtro_mediana(auxI);
        auxI = IR(VP(i,1):VP(i,2),HP(j,1):HP(j,2));
        hist = histog(auxI);
        Temp = 0;
        for k=1:(umbral+1)
            Temp = Temp + hist(k);
        end
        if Temp/(sum(hist))>porcentaje
            l = l + 1;
            P(l,1) = VP(i,1);
            P(l,2) = VP(i,2);
            P(l,3) = HP(j,1);
            P(l,4) = HP(j,2);
        end
    end
end
%filtro pasa altas en cada una de las regiones detectadas
[fP,cP] = size(P);
for i=1:fP
    auxI = acentuar(IR(P(i,1):P(i,2),P(i,3):P(i,4)));
    figure(3);
    subplot(ceil(fP/2),2,i);
    imshow(auxI);
end
%funciones utilizadas
function y = filtro_mediana(I)
    Mask = 7;
    [f,c] = size(I);
    for i=1:f
        for j=1:c
            for k=(1:Mask)-(Mask+1)/2
                a = k+(Mask+1)/2;
                for l=(1:Mask)-(Mask+1)/2
                    b = l + (Mask+1)/2;
                    if i + k<=0 | i + k>f
                        aux(b+Mask*(a-1))=255;
                    else
                        if j + l <=0 | j + l>c
                            aux(b+Mask*(a-1))=255;
                        else
                            aux(b+Mask*(a-1))=I(i+k,j+l);
                        end
                    end
                end
            end
            y(i,j) = uint8(median(aux));
        end
    end
end
function y = histog(I)
    y = zeros(256,1);
    [H_I,W_I] = size(I);
    for i=1:H_I
        for j=1:W_I
            y(uint16(I(i,j))+1) = y(uint16(I(i,j))+1) + 1;
        end
    end
end
function y = acentuar(I)
    n = 15;
    MaskS = (1/n^2)*ones(n,n);
    aux = uint8(conv2(I,MaskS,'same'));
    y = I-aux;
end