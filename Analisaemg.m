function [dados1F,resultados_rms,resultados_freq_mediana]=analisaemg_aula(dados1);
pkg load signal
dados1=detrend(dados1);

%% Filtrar os dados
freq = 2000;
[b,a] = butter(3,[20/freq/2 450/freq/2]);

dados1F(:,1)=filtfilt(b,a,dados1(:,1));

%% Corta cada contra��o dos dados 1

plot(dados1F)

[x,y,b]=ginput(20);
coord=[x y];       
inicioC1=round(coord(1,1));
inicioC2=round(coord(3,1));
inicioC3=round(coord(5,1));
inicioC4=round(coord(7,1));
inicioC5=round(coord(9,1));
inicioC6=round(coord(11,1));
inicioC7=round(coord(13,1));
inicioC8=round(coord(15,1));
inicioC9=round(coord(17,1));
inicioC10=round(coord(19,1));

fimC1=round(coord(2,1));
fimC2=round(coord(4,1));
fimC3=round(coord(6,1));
fimC4=round(coord(8,1));
fimC5=round(coord(10,1));
fimC6=round(coord(12,1));
fimC7=round(coord(14,1));
fimC8=round(coord(16,1));
fimC9=round(coord(18,1));
fimC10=round(coord(20,1));

rms(1,1)=sqrt(mean(dados1F(inicioC1:fimC1,1).^2))./max(abs(dados1F))*100;
rms(1,2)=sqrt(mean(dados1F(inicioC2:fimC2,1).^2))./max(abs(dados1F))*100;
rms(1,3)=sqrt(mean(dados1F(inicioC3:fimC3,1).^2))./max(abs(dados1F))*100;
rms(1,4)=sqrt(mean(dados1F(inicioC4:fimC4,1).^2))./max(abs(dados1F))*100;
rms(1,5)=sqrt(mean(dados1F(inicioC5:fimC5,1).^2))./max(abs(dados1F))*100;
rms(1,6)=sqrt(mean(dados1F(inicioC6:fimC6,1).^2))./max(abs(dados1F))*100;
rms(1,7)=sqrt(mean(dados1F(inicioC7:fimC7,1).^2))./max(abs(dados1F))*100;
rms(1,8)=sqrt(mean(dados1F(inicioC8:fimC8,1).^2))./max(abs(dados1F))*100;
rms(1,9)=sqrt(mean(dados1F(inicioC9:fimC9,1).^2))./max(abs(dados1F))*100;
rms(1,10)=sqrt(mean(dados1F(inicioC10:fimC10,1).^2))./max(abs(dados1F))*100;


%% Calcula a frequencia mediana
nlin=size(dados1F(inicioC1:fimC1,1),1);

%% Calcula a frequ�ncia mediana do sinal
Y=fft(dados1F(inicioC1:fimC1,1),nlin);
PSY=Y.*conj(Y)/nlin;
f = freq*(0:(nlin/2)-1)/nlin;

forca_total=trapz(f(1:nlin/2)',PSY(1:nlin/2));
n=3;
area_50=0;

while area_50<=forca_total*0.5;
    
    area_50(n-2)=trapz(f(2:n)',PSY(2:n));
    n=n+1;
end
Fm(1,1)=f(n-1);


nlin=size(dados1F(inicioC2:fimC2,1),1);

%% Calcula a frequ�ncia mediana do sinal
Y=fft(dados1F(inicioC2:fimC2,1),nlin);
PSY=Y.*conj(Y)/nlin;
f = freq*(0:(nlin/2)-1)/nlin;

forca_total=trapz(f(1:nlin/2)',PSY(1:nlin/2));
n=3;
area_50=0;

while area_50<=forca_total*0.5;
    
    area_50(n-2)=trapz(f(2:n)',PSY(2:n));
    n=n+1;
end
Fm(1,2)=f(n-1);

nlin=size(dados1F(inicioC3:fimC3,1),1);

%% Calcula a frequ�ncia mediana do sinal
Y=fft(dados1F(inicioC3:fimC3,1),nlin);
PSY=Y.*conj(Y)/nlin;
f = freq*(0:(nlin/2)-1)/nlin;

forca_total=trapz(f(1:nlin/2)',PSY(1:nlin/2));
n=3;
area_50=0;

while area_50<=forca_total*0.5;
    
    area_50(n-2)=trapz(f(2:n)',PSY(2:n));
    n=n+1;
end
Fm(1,3)=f(n-1);


nlin=size(dados1F(inicioC4:fimC4,1),1);

%% Calcula a frequ�ncia mediana do sinal
Y=fft(dados1F(inicioC4:fimC4,1),nlin);
PSY=Y.*conj(Y)/nlin;
f = freq*(0:(nlin/2)-1)/nlin;

forca_total=trapz(f(1:nlin/2)',PSY(1:nlin/2));
n=3;
area_50=0;

while area_50<=forca_total*0.5;
    
    area_50(n-2)=trapz(f(2:n)',PSY(2:n));
    n=n+1;
end
Fm(1,4)=f(n-1);


nlin=size(dados1F(inicioC5:fimC5,1),1);

%% Calcula a frequ�ncia mediana do sinal
Y=fft(dados1F(inicioC5:fimC5,1),nlin);
PSY=Y.*conj(Y)/nlin;
f = freq*(0:(nlin/2)-1)/nlin;

forca_total=trapz(f(1:nlin/2)',PSY(1:nlin/2));
n=3;
area_50=0;

while area_50<=forca_total*0.5;
    
    area_50(n-2)=trapz(f(2:n)',PSY(2:n));
    n=n+1;
end
Fm(1,5)=f(n-1);

nlin=size(dados1F(inicioC6:fimC6,1),1);

%% Calcula a frequ�ncia mediana do sinal
Y=fft(dados1F(inicioC6:fimC6,1),nlin);
PSY=Y.*conj(Y)/nlin;
f = freq*(0:(nlin/2)-1)/nlin;

forca_total=trapz(f(1:nlin/2)',PSY(1:nlin/2));
n=3;
area_50=0;

while area_50<=forca_total*0.5;
    
    area_50(n-2)=trapz(f(2:n)',PSY(2:n));
    n=n+1;
end
Fm(1,6)=f(n-1);


nlin=size(dados1F(inicioC7:fimC7,1),1);

%% Calcula a frequ�ncia mediana do sinal
Y=fft(dados1F(inicioC7:fimC7,1),nlin);
PSY=Y.*conj(Y)/nlin;
f = freq*(0:(nlin/2)-1)/nlin;

forca_total=trapz(f(1:nlin/2)',PSY(1:nlin/2));
n=3;
area_50=0;

while area_50<=forca_total*0.5;
    
    area_50(n-2)=trapz(f(2:n)',PSY(2:n));
    n=n+1;
end
Fm(1,7)=f(n-1);


nlin=size(dados1F(inicioC8:fimC8,1),1);

%% Calcula a frequ�ncia mediana do sinal
Y=fft(dados1F(inicioC8:fimC8,1),nlin);
PSY=Y.*conj(Y)/nlin;
f = freq*(0:(nlin/2)-1)/nlin;

forca_total=trapz(f(1:nlin/2)',PSY(1:nlin/2));
n=3;
area_50=0;

while area_50<=forca_total*0.5;
    
    area_50(n-2)=trapz(f(2:n)',PSY(2:n));
    n=n+1;
end
Fm(1,8)=f(n-1);

nlin=size(dados1F(inicioC9:fimC9,1),1);

%% Calcula a frequ�ncia mediana do sinal
Y=fft(dados1F(inicioC9:fimC9,1),nlin);
PSY=Y.*conj(Y)/nlin;
f = freq*(0:(nlin/2)-1)/nlin;

forca_total=trapz(f(1:nlin/2)',PSY(1:nlin/2));
n=3;
area_50=0;

while area_50<=forca_total*0.5;
    
    area_50(n-2)=trapz(f(2:n)',PSY(2:n));
    n=n+1;
end
Fm(1,9)=f(n-1);


nlin=size(dados1F(inicioC10:fimC10,1),1);

%% Calcula a frequ�ncia mediana do sinal
Y=fft(dados1F(inicioC10:fimC10,1),nlin);
PSY=Y.*conj(Y)/nlin;
f = freq*(0:(nlin/2)-1)/nlin;

forca_total=trapz(f(1:nlin/2)',PSY(1:nlin/2));
n=3;
area_50=0;

while area_50<=forca_total*0.5;
    
    area_50(n-2)=trapz(f(2:n)',PSY(2:n));
    n=n+1;
end
Fm(1,10)=f(n-1);

resultados_rms=rms;
resultados_freq_mediana=Fm;












