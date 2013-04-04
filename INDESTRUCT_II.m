diary on
format long
fprintf('\n\t\t\t ‘ ; .........................................................');
fprintf('\n\t\t\t ‘ ; .........................................................');
fprintf('\n\t\t\t ‘ ; .........................................................');
fprintf('\n\t\t\t ‘;|                                                         |');
fprintf('\n\t\t\t ‘;|          UNIVERSIDAD AUTÓNOMA DE NUEVO LEÓN             |');
fprintf('\n\t\t\t ‘;|                                                         |');
fprintf('\n\t\t\t ‘;|              FACULTAD DE INGENIERÍA CIVIL               |');
fprintf('\n\t\t\t ‘;|                 ANALISIS ESTRUCTURAL II                 |');
fprintf('\n\t\t\t ‘;|                                                         |');
fprintf('\n\t\t\t ‘;|                 PROGRAMA PARA RESOLVER:                 |');
fprintf('\n\t\t\t ‘;|                                                         |');
fprintf('\n\t\t\t ‘;|                MARCOS PLANOS CONSIDERANDO               |');
fprintf('\n\t\t\t ‘;|            DEFORMACIONES AXIALES Y CORTANTES            |');
fprintf('\n\t\t\t ‘;|                                                         |');
fprintf('\n\t\t\t ‘ ; .........................................................');
fprintf('\n\t\t\t ‘;|                     ¡INDESTRUCT!                        |');
fprintf('\n\t\t\t ‘ ; .........................................................');
fprintf('\n\t\t\t ‘;|               ING. EDWIN M. R. MARTÍNEZ                 |');
fprintf('\n\t\t\t ‘;|                                                         |');
fprintf('\n\t\t\t ‘;|                       2013-03-13                        |');
fprintf('\n\t\t\t ‘;|                                                         |');
fprintf('\n\t\t\t  ‘;..........................................................');
fprintf('\n\t\t\t ‘ ; .........................................................');
fprintf('\n\t\t\t ‘ ; .........................................................');
fprintf('\n\n\n');
DV =menu('ESCOJA UNA DE LAS SIGUIENTES OPCIONES DE ANÁLISIS:','TEORÍA DE PRIMER ORDEN','TEORÍA DE SEGUNDO ORDEN');
if DV==1
    clc
    clear all
    fprintf('\n\t\t\t ‘;|                     ¡INDESTRUCT I!                        |');
    fprintf('\n\n\t\t PROGRAMA PARA CALCULAR MARCOS PLANOS CON NUDOS RÍGIDOS Y CORTANTE');
    % -------- ARREGLO CG --------
    % EN ESTA ETAPA SE DESARROLLA LA MATRIZ DE COORDENADAS GENERALIZADAS DONDE
    % SE ESPECIFICAN LOS GRADOS DE LIBERTAD EXISTENTES DEPENDIENDO DEL TIPO DE
    % APOYO
    fprintf('\n\t\t\t COORDENADAS     Y|                                 ');
    fprintf('\n\t\t\t        GLOBALES  |                                 ');
    fprintf('\n\t\t\t                  |                                 ');
    fprintf('\n\t\t\t                  |_ _ _ _ _                        ');
    fprintf('\n\t\t\t                              X                     ');
    fprintf('\n\t\t\t                                                    ');
    %%
    nod = input('\n\n\n INGRESE EL NÚMERO DE NUDOS:');
    nnr = input(' NÚMERO DE NUDOS RESTRINGIDOS:');
    CG = ones(nod,3 );
    for i=1: nnr
        nudo_res = input('\n NÚMERO DEL NUDO RESTRINGIDO:');
        x1 = input('\n EXISTEN DESPLAZAMIENTO EN X,S o N:','s');
        if x1=='N'
            CG(nudo_res,1)=0;
        else,end
        y1 = input(' EXISTEN DESPLAZAMIENTO EN Y,S o N:','s');
        if y1=='N'
            CG(nudo_res,2)=0;
        else,end
        r1=input(' EXISTE ROTACIÓN, S o N:','s');
        if r1=='N'
            CG(nudo_res,3)=0;
        else,end
   end
   ngl=0;
   for i=1:nod
       for j=1:3
           if CG(i,j)~=0
               ngl=ngl+1;
               CG(i,j)=ngl;
         else,end
      end
   end
   CG
   ngl
   
   %--------------ARREGLO VC ----------------------------------------------
   mbr = input('\n INGRESE EL NUMERO DE MIEMBROS:');
   for i=1:mbr
       fprintf('\n MIEMBRO %d:',i);
       NUDO_INICIAL(i)=input('\n INGRESE EL NUDO INICIAL:');
       NUDO_FINAL(i)=input(' INGRESE EL NUDO FINAL:');
   end
   NUDO_INICIAL
   NUDO_FINAL
   for i=1:mbr
       for j=1:3
           VC(i,j)=CG(NUDO_INICIAL(i),j);
           VC(i,j+3)=CG(NUDO_FINAL(i),j);
       end
   end
   VC
   %--------------- ARREGLO L, SENO, COSENO -------------------------------
   fprintf('\n COORDENADAS EN LOS NUDOS:\n');
   for i=1:nod
       fprintf('\n NUDO %d:',i);
       x(i)=input('\n COORDENADA [m]:x:');
       y(i)=input(' COORDENADA [m]:y:');
   end
   for i=1:mbr
       Dx(i)=x(NUDO_FINAL(i))-x(NUDO_INICIAL(i));
       Dy(i)=y(NUDO_FINAL(i))-y(NUDO_INICIAL(i));
       L(i)=(Dx(i)^2+Dy(i)^2)^0.5;
       SENO(i)=Dy(i)/L(i);
       COSENO(i)=Dx(i)/L(i);
   end
   L
   SENO
   COSENO   
   ZZZ=menu('¿LOS DATOS INTRODUCIDOS SON CORRECTOS?','SÍ','NO');
   fprintf('\n\n\n'); 
   if ZZZ==2   
       fprintf('\n\n\t\t FIN DEL PROGRAMA, PRESIONE Y VOLVER A INICIAR');
       pause
       open('edwINDESTRUCT.m')
       break
   else
   end
   %----------------- K3-------------------------------------------------
   fprintf('\n CARACTERÍSTICAS DE LOS MIEMBROS:');
   ELAS=input('\n MÓDULO DE ELASTICIDAD [T/m^2]:');
   POISSON=input('\n MÓDULO DE POISSON (0.2):');
   GCOR=ELAS/(2*(1+POISSON));
   M_CORTANTE=GCOR
   T2_3=zeros(6,6);
   fprintf('\n\n UNICAMENTE EN LA SECCIÓN RECTANGULAR PRISMÁTICA SE PUEDEN MODELAR EXTREMOS INFINITAMENTE RÍGIDOS');
   for i=1:mbr
       if i==1
           fprintf('\n MIEMBRO %d:',i);
           fprintf('\n\n ESCOJA UNA DE LAS SIGUIENTES OPCIONES:');
           fprintf('\n\n\t\t (1) MIEMBRO PRISMÁTICO DE SECCIÓN RECTANGULAR');
           fprintf('\n\n\t\t (2) MIEMBRO CUADRADO DE SECCIÓN VARIABLE');
           fprintf('\n\n\t\t (3) MIEMBRO RECTANGULAR DE SECCIÓN VARIABLE');
           fprintf('\n\n\t\t (4) MIEMBRO CIRCULAR DE SECCIÓN PRISMÁTICA CONSTANTE');
           fprintf('\n\n\t\t (5) MIEMBRO CIRCULAR DE SECCIÓN VARIABLE');
           DD(i) = input('\n\n OPCION : ');
           fprintf('\n\n\n'); 
           % ELEMENTO PRISMÁTICO DE SECCIÓN CONSTANTE DE FORMA RECTANGULAR
          if DD(i)==1
              B(i)=input('\n BASE DEL ELEMENTO [m]:');
              H(i)=input(' ALTURA DEL ELEMENTO [m]:');
              % ELEMENTOS UNICAMENTE NECESARIOS PARA EL CONTEO DURANTE EL USO
              % DE ELEMENTOS DE CUALQUIERA DE LOS DISTINTOS TIPOS POSIBLES DE
              % SECCIONES
              Bi(i)=B(i);
              Hi(i)=H(i);
              Bf(i)=Bi(i);
              Hf(i)=Hi(i);
              Di(i)=B(i);
              Df(i)=Di(i);
              FF44(i)=0;
              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              BET(i)=input('\n VALOR DEL FACTOR DE FORMA (RECTANGULAR=1.2; CIRCULAR=10/9):');
              C1(i)=input('\n\n INGRESE EL VALOR DE C1 [m] (LONGITUD DEL SEGMENTO DE RIGIDEZ INFINITA EN EL NODO INICIAL):');
              C2(i)=input('\n\n INGRESE EL VALOR DE C2 [m] (LONGITUD DEL SEGMENTO DE RIGIDEZ INFINITA EN EL NODO FINAL):');
              F(i)=L(i)-C1(i)-C2(i)
              % CONSTANTES AUXILIARES PARA FUTUROS DESARROLLOS EN 3D
              AREA(i)=B(i)*H(i);
              INERx(i)=(B(i)*H(i)^3)/12;
              INER(i)=INERx(i);
              INERy(i)=(H(i)*B(i)^3)/12;
              FIx(i)=(12*ELAS*INERx(i)*BET(i))/(GCOR*AREA(i)*L(i)^2);
              FIy(i)=(12*ELAS*INERy(i)*BET(i))/(GCOR*AREA(i)*L(i)^2);
              % GENERACIÓN DE LAS VARIABLES
              FI(i)=(1/4)*FIx(i);
              C=((4*ELAS*INER(i))/L(i))*((1+FI(i))/(1+4*FI(i)));
              CP=C;
              A=((2*ELAS*INER(i))/L(i))*((1-2*FI(i))/(1+4*FI(i)));
              B=((6*ELAS*INER(i))/L(i)^2)*(1/(1+4*FI(i)));
              BP=B;
              T=(12*ELAS*INER(i))/L(i)^3*(1/(1+4*FI(i)));
              R=ELAS*AREA(i)/L(i);
              % MATRIZ DE RIGIDECES EN COORDENADAS LOCALES
              K2=[R 0 0 -R 0 0;
              0 T (B+C1(i)*T) 0 -T (B+C2(i)*T);
              0 (B+C1(i)*T) (C+2*C1(i)*B+C1(i)^2*T) 0 -(B+C1(i)*T) (A+C1(i)*B+C2(i)*B+C1(i)*C2(i)*T);
              -R 0 0 R 0 0;
              0 -T -(B+C1(i)*T) 0 T -(B+C2(i)*T);
              0 (B+C2(i)*T) (A+C1(i)*B+C2(i)*B+C1(i)*C2(i)*T) 0 -(B+C2(i)*T) (C+2*C2(i)*B+C2(i)^2*T)];
              % MATRIZ DE TRANSPORTE A COORDENADAS GLOBALES
              T2_3=zeros(6,6);
              T2_3(1,1)=COSENO(i);T2_3(1,2)=SENO(i);
              T2_3(2,1)=-SENO(i);T2_3(2,2)=COSENO(i);
              T2_3(3,3)=1;
              T2_3(4,4)=COSENO(i);T2_3(4,5)=SENO(i);
              T2_3(5,4)=-SENO(i);T2_3(5,5)=COSENO(i);
              T2_3(6,6)=1;
              % ELEMENTO CUADRADO SÓLIDO DE SECCIÓN VARIABLE
          elseif DD(i)==2
              Bi(i)=input('\n BASE INICIAL [m]:');
              Hi(i)=Bi(i);
              Bf(i)=input('\n BASE FINAL [m]:');
             % ELEMENTOS UNICAMENTE NECESARIOS PARA EL CONTEO DURANTE EL USO
             % DE ELEMENTOS DE CUALQUIERA DE LOS DISTINTOS TIPOS POSIBLES DE
             % SECCIONES
             Hf(i)=Bf(i);
             B(i)=Bi(i);
             H(i)=Hi(i);
             Di(i)=Bi(i);
             Df(i)=Di(i);
             BET(i)=1.2;
             INER(i)=(B(i)*H(i)^3)/12;
             FI(i)=1;
             AREA(i)=B(i)*H(i);
             BET(i)=1.2;
             C1(i)=0;
             C2(i)=0;
             F(i)=L(i);
             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             % ELEMENTOS DE LA MATRIZ DE FLEXIBILIDADES
             f11=(L(i)/(ELAS*Hi(i)*Hi(i)))*(Hi(i)/(Hf(i)-Hi(i)))*(1-Hi(i)/Hf(i));
             f22=(4*(L(i))^3/(ELAS*(Hi(i))^4))*(Hi(i)/Hf(i))^3+(6*L(i)/(5*GCOR*(Hi(i))^2))*(Hi(i)/(Hf(i)-Hi(i)))*(1-Hi(i)/Hf(i));
             f33=f22;
             f26=(4*(L(i))^2/(ELAS*(Hi(i))^4))*((Hi(i)/Hf(i))^3+(1/2)*(Hi(i)/(Hf(i)-Hi(i)))^2*(1+(Hi(i)/Hf(i))^2-2*(Hi(i)/Hf(i))));
             f35=f26;
             f66=(4*L(i)/(ELAS*(Hi(i))^4))*(Hi(i)/(Hf(i)-Hi(i)))*(1-(Hi(i)/Hf(i))^3);
             f55=f66;
             f44=(64*L(i)/(27*GCOR*(Hi(i))^4))*(Hi(i)/(Hf(i)-Hi(i)))*(1-(Hi(i)/Hf(i))^3);
             FF44(i)=f44;
             % ELEMENTOS DE LA MATRIZ DE RIGIDECES
             raz=1/f11;
             Detx=f22*f66-f26*f26;
             r11x=f22/Detx;
             r12x=(f26*L(i)-f22)/Detx;
             r22x=(f66*(L(i))^2-2*f26*L(i)+f22)/Detx;
             raax=(r11x+r22x+2*r12x)/(L(i))^2;
             rabx=(r11x+r12x)/L(i);
             rbax=(r22x+r12x)/L(i);
             % MATRIZ DE RIGIDECES EN COORDENADAS LOCALES
             K11=[raz 0 0 ; 0 raax rabx ; 0 rabx r11x];
             K12=[-raz 0 0 ; 0 -raax rbax ; 0 -rabx r12x];
             K21=K12';
             K22=[raz 0 0 ; 0 raax -rbax ; 0 -rbax r22x];
             K2=[K11 K12 ;K21 K22];
             % MATRIZ DE TRANSPORTE A COORDENADAS GLOBALES
             T2_3=zeros(6,6);
             T2_3(1,1)=COSENO(i);T2_3(1,2)=SENO(i);
             T2_3(2,1)=-SENO(i);T2_3(2,2)=COSENO(i);
             T2_3(3,3)=1;
             T2_3(4,4)=COSENO(i);T2_3(4,5)=SENO(i);
             T2_3(5,4)=-SENO(i);T2_3(5,5)=COSENO(i);
             T2_3(6,6)=1;
             % ELEMENTO RECTANGULAR SÓLIDO DE SECCIÓN VARIABLE
          elseif DD(i)==3
              Bi(i)=input('\n BASE DE TODO EL ELEMENTO [m]:');
              Hi(i)=input(' ALTURA INICIAL DEL ELEMENTO [m]:');
              Bf(i)=Bi(i);
              Hf(i)=input(' ALTURA FINAL DEL ELEMENTO [m]:');
              Bb=Bi(i);
              % ELEMENTOS UNICAMENTE NECESARIOS PARA EL CONTEO DURANTE EL USO
              % DE ELEMENTOS DE CUALQUIERA DE LOS DISTINTOS TIPOS POSIBLES DE
              % SECCIONES
              B(i)=Bi(i);
              H(i)=Hi(i);
              Di(i)=Bi(i);
              Df(i)=Di(i);
              INER(i)=(B(i)*H(i)^3)/12;
              FI(i)=1;
              AREA(i)=B(i)*H(i);
              BET(i)=1.2;
              C1(i)=0;
              C2(i)=0;
              F(i)=L(i);
              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              BET(i)=input('\n FACTOR DE FORMA (RECTANGULAR=1.2; CIRCULAR=10/9):');
              fprintf('\n\t\t\t ‘ ; ......................................................................................................................... ');
              fprintf('\n\t\t\t ‘ ; ......................................................................................................................... ');
              fprintf('\n\n\t\t     LA FUNCIÓN F(x) INDICA LA VARIACIÓN CON LA LONGITUD DE LA FORMA DEL ELEMENTO,');
              fprintf('\n\t\t\t                                                                    ');
              fprintf('\n\n\t\t     SI LA FUNCIÓN ES LINEAL: Y=F(x)=');
              fprintf('\n\t\t\t                                      (abs(Hi-Hf)/LL)*x                                       ');
              fprintf('\n\t\t\t ‘ ; ......................................................................................................................... ');
              fprintf('\n\t\t\t ‘ ; ......................................................................................................................... ');
              fprintf('\n\t\t\t ‘ ;                                                                                                                           ');
              fprintf('\n\t\t\t   ESTANDO EL EJE X SOBRE EL     Y|                                                                                            ');
              fprintf('\n\t\t\t    CENTROIDE DEL MIEMBRO         |                                                                                            ');
              fprintf('\n\t\t\t                                  |                                                                                            ');
              fprintf('\n\t\t\t                                  |_ _ _ _ _                                                                                   ');
              fprintf('\n\t\t\t                                           X                                                                                   ');
              fprintf('\n\t\t\t                                                                                                                               ');
              fprintf('\n\t\t\t ‘ ; ......................................................................................................................... ');
              fprintf('\n\t\t\t ‘ ; ......................................................................................................................... ');
              fprintf('\n\t\t\t                                                                                                                              ');
              fprintf('\n\t\t\t                                                                                                                              ');
              EE(i)=menu('ESCOJA UNA DE LAS SIGUIENTES OPCIONES:','FUNCIÓN LINEAL','OTRO TIPO DE VARIACIÓN');
              if EE(i)==1
                  % ELEMENTOS DE LA MATRIZ DE FLEXIBILIDADES
                  f11=(L(i)/(ELAS*Bb*Hi(i)))*(Hi(i)/(Hf(i)-Hi(i)))*log(Hf(i)/Hi(i));
                  f22=(6*(L(i))^3/(ELAS*Bb*(Hi(i))^3))*((Hi(i)/Hf(i))^2-(Hi(i)/(Hf(i)-Hi(i)))^3*(Hf(i)/Hi(i)-Hi(i)/Hf(i)-2*log(Hf(i)/Hi(i))))+(6*L(i)/(5*GCOR*Bb*Hi(i)))*(Hi(i)/(Hf(i)-Hi(i)))*log(Hf(i)/Hi(i));
                  f33=(12*(L(i))^3/(ELAS*(Bb)^3*Hi(i)))*(3/2+(1/2)*(Hf(i)/Hi(i))^2-2*(Hf(i)/Hi(i))+log(Hf(i)/Hi(i)))+(6*L(i)/(5*GCOR*Bb*Hi(i)))*(Hi(i)/(Hf(i)-Hi(i)))*log(Hf(i)/Hi(i));
                  f26=(6*L(i)/(ELAS*Bb*(Hi(i))^3))*(Hi(i)/(Hf(i)-Hi(i)))^2*(1+(Hi(i)/Hf(i))^2-2*(Hi(i)/Hf(i)));
                  f35=(12*(L(i))^2/(ELAS*(Bb)^3*Hi(i)))*(Hi(i)/(Hf(i)-Hi(i)))^2*(Hf(i)/Hi(i)-log(Hf(i)/Hi(i))-1);
                  f66=(6*L(i)/(ELAS*Bb*(Hi(i))^3))*(Hi(i)/(Hf(i)-Hi(i)))*(1-(Hi(i)/Hf(i))^2);
                  f55=(12*L(i)/(ELAS*(Bb)^3*Hi(i)))*(Hi(i)/(Hf(i)-Hi(i)))*log(Hf(i)/Hi(i));
                  % DEFINICIÓN DE UNA FUNCIÓN DE x
                  Hi=Hi(i);
                  Hf=Hf(i);
                  LL=L(i);
                  if Hi<Hf
                      dpc= @(x)(1./(((Bb.^3*(-Hi./2-((Hi-Hf)./(2*LL))*x))./3)*(1-(192*Bb./(pi.^5*(-Hi./2-((Hi-Hf)./(2*LL))*x)))*(tanh(pi*(-Hi./2-((Hi-Hf)./(2*LL))*x)./(2*Bb)))))); 
                  elseif Hi>Hf
                      dpc= @(x)(1./(((Bb.^3*(-Hi./2+((Hi-Hf)./(2*LL))*x))./3)*(1-(192*Bb./(pi.^5*(-Hi./2+((Hi-Hf)./(2*LL))*x)))*(tanh(pi*(-Hi./2+((Hi-Hf)./(2*LL))*x)./(2*Bb))))));
                  end
                  f44=quad(@(x)dpc, 0, LL)
                  FF44(i)=f44;
                  % ELEMENTOS DE LA MATRIZ DE RIGIDECES
                  raz=1/f11;
                  Detx=f22*f66-f26*f26;
                  r11x=f22/Detx;
                  r12x=(f26*L(i)-f22)/Detx;
                  r22x=(f66*(L(i))^2-2*f26*L(i)+f22)/Detx;
                  raax=(r11x+r22x+2*r12x)/(L(i))^2;
                  rabx=(r11x+r12x)/L(i);
                  rbax=(r22x+r12x)/L(i);
                  % MATRIZ DE RIGIDECES EN COORDENADAS LOCALES
                  K11=[raz 0 0 ; 0 raax rabx ; 0 rabx r11x];
                  K12=[-raz 0 0 ; 0 -raax rbax ; 0 -rabx r12x];
                  K21=K12';
                  K22=[raz 0 0 ; 0 raax -rbax ; 0 -rbax r22x];
                  K2=[K11 K12 ;K21 K22];
                  % MATRIZ DE TRANSPORTE A COORDENADAS GLOBALES
                  T2_3=zeros(6,6);
                  T2_3(1,1)=COSENO(i);T2_3(1,2)=SENO(i);
                  T2_3(2,1)=-SENO(i);T2_3(2,2)=COSENO(i);
                  T2_3(3,3)=1;
                  T2_3(4,4)=COSENO(i);T2_3(4,5)=SENO(i);
                  T2_3(5,4)=-SENO(i);T2_3(5,5)=COSENO(i);
                  T2_3(6,6)=1;
              else
                  % ELEMENTOS DE LA MATRIZ DE FLEXIBILIDADES
                  f11=(L(i)/(ELAS*Bb*Hi(i)))*(Hi(i)/(Hf(i)-Hi(i)))*log(Hf(i)/Hi(i));
                  f22=(6*(L(i))^3/(ELAS*Bb*(Hi(i))^3))*((Hi(i)/Hf(i))^2-(Hi(i)/(Hf(i)-Hi(i)))^3*(Hf(i)/Hi(i)-Hi(i)/Hf(i)-2*log(Hf(i)/Hi(i))))+(6*L(i)/(5*GCOR*Bb*Hi(i)))*(Hi(i)/(Hf(i)-Hi(i)))*log(Hf(i)/Hi(i));
                  f33=(12*(L(i))^3/(ELAS*(Bb)^3*Hi(i)))*(3/2+(1/2)*(Hf(i)/Hi(i))^2-2*(Hf(i)/Hi(i))+log(Hf(i)/Hi(i)))+(6*L(i)/(5*GCOR*Bb*Hi(i)))*(Hi(i)/(Hf(i)-Hi(i)))*log(Hf(i)/Hi(i));
                  f26=(6*L(i)/(ELAS*Bb*(Hi(i))^3))*(Hi(i)/(Hf(i)-Hi(i)))^2*(1+(Hi(i)/Hf(i))^2-2*(Hi(i)/Hf(i)));
                  f35=(12*(L(i))^2/(ELAS*(Bb)^3*Hi(i)))*(Hi(i)/(Hf(i)-Hi(i)))^2*(Hf(i)/Hi(i)-log(Hf(i)/Hi(i))-1);
                  f66=(6*L(i)/(ELAS*Bb*(Hi(i))^3))*(Hi(i)/(Hf(i)-Hi(i)))*(1-(Hi(i)/Hf(i))^2);
                  f55=(12*L(i)/(ELAS*(Bb)^3*Hi(i)))*(Hi(i)/(Hf(i)-Hi(i)))*log(Hf(i)/Hi(i));
                  % DEFINICIÓN DE UNA FUNCIÓN DE x
                  Hi=Hi(i);
                  Hf=Hf(i);
                  LL=L(i);
                  f44=input(' INTRODUZCA LA f44=');
                  FF44(i)=f44;
                  % ELEMENTOS DE LA MATRIZ DE RIGIDECES
                  raz=1/f11;
                  Detx=f22*f66-f26*f26;
                  r11x=f22/Detx;
                  r12x=(f26*L(i)-f22)/Detx;
                  r22x=(f66*(L(i))^2-2*f26*L(i)+f22)/Detx;
                  raax=(r11x+r22x+2*r12x)/(L(i))^2;
                  rabx=(r11x+r12x)/L(i);
                  rbax=(r22x+r12x)/L(i);
                  % MATRIZ DE RIGIDECES EN COORDENADAS LOCALES
                  K11=[raz 0 0 ; 0 raax rabx ; 0 rabx r11x];
                  K12=[-raz 0 0 ; 0 -raax rbax ; 0 -rabx r12x];
                  K21=K12';
                  K22=[raz 0 0 ; 0 raax -rbax ; 0 -rbax r22x];
                  K2=[K11 K12 ;K21 K22];
                  % MATRIZ DE TRANSPORTE A COORDENADAS GLOBALES
                  T2_3=zeros(6,6);
                  T2_3(1,1)=COSENO(i);T2_3(1,2)=SENO(i);
                  T2_3(2,1)=-SENO(i);T2_3(2,2)=COSENO(i);
                  T2_3(3,3)=1;
                  T2_3(4,4)=COSENO(i);T2_3(4,5)=SENO(i);
                  T2_3(5,4)=-SENO(i);T2_3(5,5)=COSENO(i);
                  T2_3(6,6)=1;
              end
              % ELEMENTO CIRCULAR SÓLIDO DE SECCIÓN CONSTATE
          elseif DD(i)==4
              Di(i)=input('\n DIÁMETRO DEL ELEMENTO [m]:');
              % ELEMENTOS UNICAMENTE NECESARIOS PARA EL CONTEO DURANTE EL USO
              % DE ELEMENTOS DE CUALQUIERA DE LOS DISTINTOS TIPOS POSIBLES DE
              % SECCIONES
              Df(i)=Di(i);
              B(i)=Di(i);
              H(i)=Di(i);
              Bi(i)=B(i);
              Hi(i)=H(i);
              Bf(i)=Bi(i);
              Hf(i)=Hi(i);
              FF44(i)=0;
              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              BET(i)=input('\n VALOR DEL FACTOR DE FORMA (RECTANGULAR=1.2; CIRCULAR=10/9):');
              AREA(i)=pi*Di(i)*Di(i)/4;
              INER(i)=pi*(Di(i)/2)^4/4;
              C1(i)=input('\n\n INGRESE EL VALOR DE C1 [m] (LONGITUD DEL SEGMENTO DE RIGIDEZ INFINITA EN EL NODO INICIAL):');
              C2(i)=input('\n\n INGRESE EL VALOR DE C2 [m] (LONGITUD DEL SEGMENTO DE RIGIDEZ INFINITA EN EL NODO FINAL):');
              F(i)=L(i)-C1(i)-C2(i)
              FI(i)=(3*ELAS*INER(i)*BET(i))/(GCOR*AREA(i)*L(i)^2);
              C=((4*ELAS*INER(i))/L(i))*((1+FI(i))/(1+4*FI(i)));
              CP=C;
              A=((2*ELAS*INER(i))/L(i))*((1-2*FI(i))/(1+4*FI(i)));
              B=((6*ELAS*INER(i))/L(i)^2)*(1/(1+4*FI(i)));
              BP=B;
              T=(12*ELAS*INER(i))/L(i)^3*(1/(1+4*FI(i)));
              R=ELAS*AREA(i)/L(i);
              K2=[R 0 0 -R 0 0;
              0 T (B+C1(i)*T) 0 -T (B+C2(i)*T);
              0 (B+C1(i)*T) (C+2*C1(i)*B+C1(i)^2*T) 0 -(B+C1(i)*T) (A+C1(i)*B+C2(i)*B+C1(i)*C2(i)*T);
              -R 0 0 R 0 0;
              0 -T -(B+C1(i)*T) 0 T -(B+C2(i)*T);
              0 (B+C2(i)*T) (A+C1(i)*B+C2(i)*B+C1(i)*C2(i)*T) 0 -(B+C2(i)*T) (C+2*C2(i)*B+C2(i)^2*T)];
              % MATRIZ DE TRANSPORTE A COORDENADAS GLOBALES
              T2_3=zeros(6,6);
              T2_3(1,1)=COSENO(i);T2_3(1,2)=SENO(i);
              T2_3(2,1)=-SENO(i);T2_3(2,2)=COSENO(i);
              T2_3(3,3)=1;
              T2_3(4,4)=COSENO(i);T2_3(4,5)=SENO(i);
              T2_3(5,4)=-SENO(i);T2_3(5,5)=COSENO(i);
              T2_3(6,6)=1;
              % ELEMENTO CIRCULAR SÓLIDO DE SECCIÓN VARIABLE
          elseif DD(i)==5
              Di(i)=input('\n DIÁMETRO INICIAL DEL ELEMENTO [m]:');
              Df(i)=input('\n DIÁMETRO FINAL DEL ELEMENTO [m]:');
              BET(i)=input('\n VALOR DEL FACTOR DE FORMA (RECTANGULAR=1.2; CIRCULAR=10/9):');
              % ELEMENTOS UNICAMENTE NECESARIOS PARA EL CONTEO DURANTE EL USO
              % DE ELEMENTOS DE CUALQUIERA DE LOS DISTINTOS TIPOS POSIBLES DE
              % SECCIONES
              B(i)=Di(i);
              H(i)=Df(i);
              Bi(i)=B(i);
              Hi(i)=H(i);
              Bf(i)=Bi(i);
              Hf(i)=Hi(i);
              BET(i)=10/9;
              C1(i)=0;
              C2(i)=0;
              F(i)=L(i);
              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              AREA(i)=pi*Di(i)*Di(i)/4;
              INER(i)=pi*(Di(i)/2)^4/4;
              FIx(i)=(12*ELAS*INER(i)*BET(i))/(GCOR*AREA(i)*L(i)^2)
              FIy(i)=FIx(i);
              FI(i)=(1/4)*FIx(i);
              % ELEMENTOS DE LA MATRIZ DE FLEXIBILIDADES
              f11=(4*L(i)/(pi*ELAS*(Di(i))^2))*(Di(i)/(Df(i)-Di(i)))*(1-Di(i)/Df(i));
              f22=(64*(L(i))^3/(3*pi*ELAS*(Di(i))^4))*(Di(i)/Df(i))^3+(40*L(i)/(9*pi*GCOR*(Di(i))^2))*(Di(i)/(Df(i)-Di(i)))*(1-(Di(i)/Df(i)));
              f33=f22;
              f26=(64*(L(i))^2/(3*pi*ELAS*(Di(i))^4))*((Di(i)/Df(i))^3+(1/2)*(Di(i)/(Df(i)-Di(i)))^2*(1+(Di(i)/Df(i))^2-2*(Di(i)/Df(i))));
              f35=f26;
              f66=(64*L(i)/(3*pi*ELAS*(Di(i))^4))*(Di(i)/(Df(i)-Di(i)))*(1-(Di(i)/Df(i))^3);
              f55=f66;
              f44=(32*L(i)/(3*pi*GCOR*Di(i)^4))*(Di(i)/(Df(i)-Di(i)))*(1-(Di(i)/Df(i))^3);
              FF44(i)=f44;
              % ELEMENTOS DE LA MATRIZ DE RIGIDECES
              raz=1/f11;
              Detx=f22*f66-f26*f26;
              r11x=f22/Detx;
              r12x=(f26*L(i)-f22)/Detx;
              r22x=(f66*(L(i))^2-2*f26*L(i)+f22)/Detx;
              raax=(r11x+r22x+2*r12x)/(L(i))^2;
              rabx=(r11x+r12x)/L(i);
              rbax=(r22x+r12x)/L(i);
              % MATRIZ DE RIGIDECES EN COORDENADAS LOCALES
              K11=[raz 0 0 ; 0 raax rabx ; 0 rabx r11x];
              K12=[-raz 0 0 ; 0 -raax rbax ; 0 -rabx r12x];
              K21=K12';
              K22=[raz 0 0 ; 0 raax -rbax ; 0 -rbax r22x];
              K2=[K11 K12 ;K21 K22];
              % MATRIZ DE TRANSPORTE A COORDENADAS GLOBALES
              T2_3=zeros(6,6);
              T2_3(1,1)=COSENO(i);T2_3(1,2)=SENO(i);
              T2_3(2,1)=-SENO(i);T2_3(2,2)=COSENO(i);
              T2_3(3,3)=1;
              T2_3(4,4)=COSENO(i);T2_3(4,5)=SENO(i);
              T2_3(5,4)=-SENO(i);T2_3(5,5)=COSENO(i);
              T2_3(6,6)=1;
          end
       elseif i==mbr
           fprintf('\n MIEMBRO %d:',i);
           fprintf('\n\n ESCOJA UNA DE LAS SIGUIENTES OPCIONES:');
           fprintf('\n\n\t\t (1) MIEMBRO PRISMÁTICO DE SECCIÓN RECTANGULAR');
           fprintf('\n\n\t\t (2) MIEMBRO CUADRADO DE SECCIÓN VARIABLE');
           fprintf('\n\n\t\t (3) MIEMBRO RECTANGULAR DE SECCIÓN VARIABLE');
           fprintf('\n\n\t\t (4) MIEMBRO CIRCULAR DE SECCIÓN PRISMÁTICA CONSTANTE');
           fprintf('\n\n\t\t (5) MIEMBRO CIRCULAR DE SECCIÓN VARIABLE');
           DD(i) = input('\n\n OPCION : ');
           fprintf('\n\n\n'); 
           % ELEMENTO PRISMÁTICO DE SECCIÓN CONSTANTE DE FORMA RECTANGULAR
           if DD(i)==1
               B(i)=input('\n BASE DEL ELEMENTO [m]:');
               H(i)=input(' ALTURA DEL ELEMENTO [m]:');
               % ELEMENTOS UNICAMENTE NECESARIOS PARA EL CONTEO DURANTE EL USO
               % DE ELEMENTOS DE CUALQUIERA DE LOS DISTINTOS TIPOS POSIBLES DE
               % SECCIONES
               Bi(i)=B(i);
               Hi(i)=H(i);
               Bf(i)=Bi(i);
               Hf(i)=Hi(i);
               Di(i)=B(i);
               Df(i)=Di(i);
               FF44(i)=0;
               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               BET(i)=input('\n VALOR DEL FACTOR DE FORMA (RECTANGULAR=1.2; CIRCULAR=10/9):');
               C1(i)=input('\n\n INGRESE EL VALOR DE C1 [m] (LONGITUD DEL SEGMENTO DE RIGIDEZ INFINITA EN EL NODO INICIAL):');
               C2(i)=input('\n\n INGRESE EL VALOR DE C2 [m] (LONGITUD DEL SEGMENTO DE RIGIDEZ INFINITA EN EL NODO FINAL):');
               F(i)=L(i)-C1(i)-C2(i)
               % CONSTANTES AUXILIARES PARA FUTUROS DESARROLLOS EN 3D
               AREA(i)=B(i)*H(i);
               INERx(i)=(B(i)*H(i)^3)/12;
               INER(i)=INERx(i);
               INERy(i)=(H(i)*B(i)^3)/12;
               FIx(i)=(12*ELAS*INERx(i)*BET(i))/(GCOR*AREA(i)*L(i)^2);
               FIy(i)=(12*ELAS*INERy(i)*BET(i))/(GCOR*AREA(i)*L(i)^2);
               % GENERACIÓN DE LAS VARIABLES
               FI(i)=(1/4)*FIx(i);
               C=((4*ELAS*INER(i))/L(i))*((1+FI(i))/(1+4*FI(i)));
               CP=C;
               A=((2*ELAS*INER(i))/L(i))*((1-2*FI(i))/(1+4*FI(i)));
               B=((6*ELAS*INER(i))/L(i)^2)*(1/(1+4*FI(i)));
               BP=B;
               T=(12*ELAS*INER(i))/L(i)^3*(1/(1+4*FI(i)));
               R=ELAS*AREA(i)/L(i);
               % MATRIZ DE RIGIDECES EN COORDENADAS LOCALES
               K2=[R 0 0 -R 0 0;
               0 T (B+C1(i)*T) 0 -T (B+C2(i)*T);
               0 (B+C1(i)*T) (C+2*C1(i)*B+C1(i)^2*T) 0 -(B+C1(i)*T) (A+C1(i)*B+C2(i)*B+C1(i)*C2(i)*T);
               -R 0 0 R 0 0;
               0 -T -(B+C1(i)*T) 0 T -(B+C2(i)*T);
               0 (B+C2(i)*T) (A+C1(i)*B+C2(i)*B+C1(i)*C2(i)*T) 0 -(B+C2(i)*T) (C+2*C2(i)*B+C2(i)^2*T)];
               % MATRIZ DE TRANSPORTE A COORDENADAS GLOBALES
               T2_3=zeros(6,6);
               T2_3(1,1)=COSENO(i);T2_3(1,2)=SENO(i);
               T2_3(2,1)=-SENO(i);T2_3(2,2)=COSENO(i);
               T2_3(3,3)=1;
               T2_3(4,4)=COSENO(i);T2_3(4,5)=SENO(i);
               T2_3(5,4)=-SENO(i);T2_3(5,5)=COSENO(i);
               T2_3(6,6)=1;
               % ELEMENTO CUADRADO SÓLIDO DE SECCIÓN VARIABLE
           elseif DD(i)==2
               Bi(i)=input('\n BASE INICIAL [m]:');
               Hi(i)=Bi(i);
               Bf(i)=input('\n BASE FINAL [m]:');
              % ELEMENTOS UNICAMENTE NECESARIOS PARA EL CONTEO DURANTE EL USO
              % DE ELEMENTOS DE CUALQUIERA DE LOS DISTINTOS TIPOS POSIBLES DE
              % SECCIONES
              Hf(i)=Bf(i);
              B(i)=Bi(i);
              H(i)=Hi(i);
              Di(i)=Bi(i);
              Df(i)=Di(i);
              BET(i)=1.2;
              INER(i)=(B(i)*H(i)^3)/12;
              FI(i)=1;
              AREA(i)=B(i)*H(i);
              C1(i)=0;
              C2(i)=0;
              F(i)=L(i);
              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              % ELEMENTOS DE LA MATRIZ DE FLEXIBILIDADES
              f11=(L(i)/(ELAS*Hi(i)*Hi(i)))*(Hi(i)/(Hf(i)-Hi(i)))*(1-Hi(i)/Hf(i));
              f22=(4*(L(i))^3/(ELAS*(Hi(i))^4))*(Hi(i)/Hf(i))^3+(6*L(i)/(5*GCOR*(Hi(i))^2))*(Hi(i)/(Hf(i)-Hi(i)))*(1-Hi(i)/Hf(i));
              f33=f22;
              f26=(4*(L(i))^2/(ELAS*(Hi(i))^4))*((Hi(i)/Hf(i))^3+(1/2)*(Hi(i)/(Hf(i)-Hi(i)))^2*(1+(Hi(i)/Hf(i))^2-2*(Hi(i)/Hf(i))));
              f35=f26;
              f66=(4*L(i)/(ELAS*(Hi(i))^4))*(Hi(i)/(Hf(i)-Hi(i)))*(1-(Hi(i)/Hf(i))^3);
              f55=f66;
              f44=(64*L(i)/(27*GCOR*(Hi(i))^4))*(Hi(i)/(Hf(i)-Hi(i)))*(1-(Hi(i)/Hf(i))^3);
              FF44(i)=f44;
              % ELEMENTOS DE LA MATRIZ DE RIGIDECES
              raz=1/f11;
              Detx=f22*f66-f26*f26;
              r11x=f22/Detx;
              r12x=(f26*L(i)-f22)/Detx;
              r22x=(f66*(L(i))^2-2*f26*L(i)+f22)/Detx;
              raax=(r11x+r22x+2*r12x)/(L(i))^2;
              rabx=(r11x+r12x)/L(i);
              rbax=(r22x+r12x)/L(i);
              % MATRIZ DE RIGIDECES EN COORDENADAS LOCALES
              K11=[raz 0 0 ; 0 raax rabx ; 0 rabx r11x];
              K12=[-raz 0 0 ; 0 -raax rbax ; 0 -rabx r12x];
              K21=K12';
              K22=[raz 0 0 ; 0 raax -rbax ; 0 -rbax r22x];
              K2=[K11 K12 ;K21 K22];
              % MATRIZ DE TRANSPORTE A COORDENADAS GLOBALES
              T2_3=zeros(6,6);
              T2_3(1,1)=COSENO(i);T2_3(1,2)=SENO(i);
              T2_3(2,1)=-SENO(i);T2_3(2,2)=COSENO(i);
              T2_3(3,3)=1;
              T2_3(4,4)=COSENO(i);T2_3(4,5)=SENO(i);
              T2_3(5,4)=-SENO(i);T2_3(5,5)=COSENO(i);
              T2_3(6,6)=1;
              % ELEMENTO RECTANGULAR SÓLIDO DE SECCIÓN VARIABLE
           elseif DD(i)==3
               Bi(i)=input('\n BASE DE TODO EL ELEMENTO [m]:');
               Hi(i)=input(' ALTURA INICIAL DEL ELEMENTO [m]:');
               Bf(i)=Bi(i);
               Hf(i)=input(' ALTURA FINAL DEL ELEMENTO [m]:');
               Bb=Bi(i);
               % ELEMENTOS UNICAMENTE NECESARIOS PARA EL CONTEO DURANTE EL USO
               % DE ELEMENTOS DE CUALQUIERA DE LOS DISTINTOS TIPOS POSIBLES DE
               % SECCIONES
               B(i)=Bi(i);
               H(i)=Hi(i);
               Di(i)=Bi(i);
               Df(i)=Di(i);
               INER(i)=(B(i)*H(i)^3)/12;
               FI(i)=1;
               AREA(i)=B(i)*H(i);
               BET(i)=1.2;
               C1(i)=0;
               C2(i)=0;
               F(i)=L(i);
               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               BET(i)=input('\n FACTOR DE FORMA (RECTANGULAR=1.2; CIRCULAR=10/9):');
               fprintf('\n\t\t\t ‘ ; ......................................................................................................................... ');
               fprintf('\n\t\t\t ‘ ; ......................................................................................................................... ');
               fprintf('\n\n\t\t     LA FUNCIÓN F(x) INDICA LA VARIACIÓN CON LA LONGITUD DE LA FORMA DEL ELEMENTO,');
               fprintf('\n\t\t\t                                                                    ');
               fprintf('\n\n\t\t     SI LA FUNCIÓN ES LINEAL: Y=F(x)=');
               fprintf('\n\t\t\t                                      (abs(Hi-Hf)/LL)*x                                       ');
               fprintf('\n\t\t\t ‘ ; ......................................................................................................................... ');
               fprintf('\n\t\t\t ‘ ; ......................................................................................................................... ');
               fprintf('\n\t\t\t ‘ ;                                                                                                                           ');
               fprintf('\n\t\t\t   ESTANDO EL EJE X SOBRE EL     Y|                                                                                            ');
               fprintf('\n\t\t\t    CENTROIDE DEL MIEMBRO         |                                                                                            ');
               fprintf('\n\t\t\t                                  |                                                                                            ');
               fprintf('\n\t\t\t                                  |_ _ _ _ _                                                                                   ');
               fprintf('\n\t\t\t                                           X                                                                                   ');
               fprintf('\n\t\t\t                                                                                                                               ');
               fprintf('\n\t\t\t ‘ ; ......................................................................................................................... ');
               fprintf('\n\t\t\t ‘ ; ......................................................................................................................... ');
               fprintf('\n\t\t\t                                                                                                                              ');
               fprintf('\n\t\t\t                                                                                                                              ');
               EE(i)=menu('ESCOJA UNA DE LAS SIGUIENTES OPCIONES:','FUNCIÓN LINEAL','OTRO TIPO DE VARIACIÓN');
               if EE(i)==1
                   % ELEMENTOS DE LA MATRIZ DE FLEXIBILIDADES
                   f11=(L(i)/(ELAS*Bb*Hi(i)))*(Hi(i)/(Hf(i)-Hi(i)))*log(Hf(i)/Hi(i));
                   f22=(6*(L(i))^3/(ELAS*Bb*(Hi(i))^3))*((Hi(i)/Hf(i))^2-(Hi(i)/(Hf(i)-Hi(i)))^3*(Hf(i)/Hi(i)-Hi(i)/Hf(i)-2*log(Hf(i)/Hi(i))))+(6*L(i)/(5*GCOR*Bb*Hi(i)))*(Hi(i)/(Hf(i)-Hi(i)))*log(Hf(i)/Hi(i));
                   f33=(12*(L(i))^3/(ELAS*(Bb)^3*Hi(i)))*(3/2+(1/2)*(Hf(i)/Hi(i))^2-2*(Hf(i)/Hi(i))+log(Hf(i)/Hi(i)))+(6*L(i)/(5*GCOR*Bb*Hi(i)))*(Hi(i)/(Hf(i)-Hi(i)))*log(Hf(i)/Hi(i));
                   f26=(6*L(i)/(ELAS*Bb*(Hi(i))^3))*(Hi(i)/(Hf(i)-Hi(i)))^2*(1+(Hi(i)/Hf(i))^2-2*(Hi(i)/Hf(i)));
                   f35=(12*(L(i))^2/(ELAS*(Bb)^3*Hi(i)))*(Hi(i)/(Hf(i)-Hi(i)))^2*(Hf(i)/Hi(i)-log(Hf(i)/Hi(i))-1);
                   f66=(6*L(i)/(ELAS*Bb*(Hi(i))^3))*(Hi(i)/(Hf(i)-Hi(i)))*(1-(Hi(i)/Hf(i))^2);
                   f55=(12*L(i)/(ELAS*(Bb)^3*Hi(i)))*(Hi(i)/(Hf(i)-Hi(i)))*log(Hf(i)/Hi(i));
                   % DEFINICIÓN DE UNA FUNCIÓN DE x
                   Hi=Hi(i);
                   Hf=Hf(i);
                   LL=L(i);
                   if Hi<Hf
                      dpc= @(x)(1./(((Bb.^3*(-Hi./2-((Hi-Hf)./(2*LL))*x))./3)*(1-(192*Bb./(pi.^5*(-Hi./2-((Hi-Hf)./(2*LL))*x)))*(tanh(pi*(-Hi./2-((Hi-Hf)./(2*LL))*x)./(2*Bb)))))); 
                  elseif Hi>Hf
                      dpc= @(x)(1./(((Bb.^3*(-Hi./2+((Hi-Hf)./(2*LL))*x))./3)*(1-(192*Bb./(pi.^5*(-Hi./2+((Hi-Hf)./(2*LL))*x)))*(tanh(pi*(-Hi./2+((Hi-Hf)./(2*LL))*x)./(2*Bb))))));
                  end
                   f44=quad(@(x)dpc, 0, LL)
                   FF44(i)=f44;
                   % ELEMENTOS DE LA MATRIZ DE RIGIDECES
                   raz=1/f11;
                   Detx=f22*f66-f26*f26;
                   r11x=f22/Detx;
                   r12x=(f26*L(i)-f22)/Detx;
                   r22x=(f66*(L(i))^2-2*f26*L(i)+f22)/Detx;
                   raax=(r11x+r22x+2*r12x)/(L(i))^2;
                   rabx=(r11x+r12x)/L(i);
                   rbax=(r22x+r12x)/L(i);
                   % MATRIZ DE RIGIDECES EN COORDENADAS LOCALES
                   K11=[raz 0 0 ; 0 raax rabx ; 0 rabx r11x];
                   K12=[-raz 0 0 ; 0 -raax rbax ; 0 -rabx r12x];
                   K21=K12';
                   K22=[raz 0 0 ; 0 raax -rbax ; 0 -rbax r22x];
                   K2=[K11 K12 ;K21 K22];
                   % MATRIZ DE TRANSPORTE A COORDENADAS GLOBALES
                   T2_3=zeros(6,6);
                   T2_3(1,1)=COSENO(i);T2_3(1,2)=SENO(i);
                   T2_3(2,1)=-SENO(i);T2_3(2,2)=COSENO(i);
                   T2_3(3,3)=1;
                   T2_3(4,4)=COSENO(i);T2_3(4,5)=SENO(i);
                   T2_3(5,4)=-SENO(i);T2_3(5,5)=COSENO(i);
                   T2_3(6,6)=1;
               else
                   % ELEMENTOS DE LA MATRIZ DE FLEXIBILIDADES
                   f11=(L(i)/(ELAS*Bb*Hi(i)))*(Hi(i)/(Hf(i)-Hi(i)))*log(Hf(i)/Hi(i));
                   f22=(6*(L(i))^3/(ELAS*Bb*(Hi(i))^3))*((Hi(i)/Hf(i))^2-(Hi(i)/(Hf(i)-Hi(i)))^3*(Hf(i)/Hi(i)-Hi(i)/Hf(i)-2*log(Hf(i)/Hi(i))))+(6*L(i)/(5*GCOR*Bb*Hi(i)))*(Hi(i)/(Hf(i)-Hi(i)))*log(Hf(i)/Hi(i));
                   f33=(12*(L(i))^3/(ELAS*(Bb)^3*Hi(i)))*(3/2+(1/2)*(Hf(i)/Hi(i))^2-2*(Hf(i)/Hi(i))+log(Hf(i)/Hi(i)))+(6*L(i)/(5*GCOR*Bb*Hi(i)))*(Hi(i)/(Hf(i)-Hi(i)))*log(Hf(i)/Hi(i));
                   f26=(6*L(i)/(ELAS*Bb*(Hi(i))^3))*(Hi(i)/(Hf(i)-Hi(i)))^2*(1+(Hi(i)/Hf(i))^2-2*(Hi(i)/Hf(i)));
                   f35=(12*(L(i))^2/(ELAS*(Bb)^3*Hi(i)))*(Hi(i)/(Hf(i)-Hi(i)))^2*(Hf(i)/Hi(i)-log(Hf(i)/Hi(i))-1);
                   f66=(6*L(i)/(ELAS*Bb*(Hi(i))^3))*(Hi(i)/(Hf(i)-Hi(i)))*(1-(Hi(i)/Hf(i))^2);
                   f55=(12*L(i)/(ELAS*(Bb)^3*Hi(i)))*(Hi(i)/(Hf(i)-Hi(i)))*log(Hf(i)/Hi(i));
                   % DEFINICIÓN DE UNA FUNCIÓN DE x
                   Hi=Hi(i);
                   Hf=Hf(i);
                   LL=L(i);
                   f44=input(' INTRODUZCA LA f44=');
                   FF44(i)=f44;
                   % ELEMENTOS DE LA MATRIZ DE RIGIDECES
                   raz=1/f11;
                   Detx=f22*f66-f26*f26;
                   r11x=f22/Detx;
                   r12x=(f26*L(i)-f22)/Detx;
                   r22x=(f66*(L(i))^2-2*f26*L(i)+f22)/Detx;
                   raax=(r11x+r22x+2*r12x)/(L(i))^2;
                   rabx=(r11x+r12x)/L(i);
                   rbax=(r22x+r12x)/L(i);
                   % MATRIZ DE RIGIDECES EN COORDENADAS LOCALES
                   K11=[raz 0 0 ; 0 raax rabx ; 0 rabx r11x];
                   K12=[-raz 0 0 ; 0 -raax rbax ; 0 -rabx r12x];
                   K21=K12';
                   K22=[raz 0 0 ; 0 raax -rbax ; 0 -rbax r22x];
                   K2=[K11 K12 ;K21 K22];
                   % MATRIZ DE TRANSPORTE A COORDENADAS GLOBALES
                   T2_3=zeros(6,6);
                   T2_3(1,1)=COSENO(i);T2_3(1,2)=SENO(i);
                   T2_3(2,1)=-SENO(i);T2_3(2,2)=COSENO(i);
                   T2_3(3,3)=1;
                   T2_3(4,4)=COSENO(i);T2_3(4,5)=SENO(i);
                   T2_3(5,4)=-SENO(i);T2_3(5,5)=COSENO(i);
                   T2_3(6,6)=1;
               end
              % ELEMENTO CIRCULAR SÓLIDO DE SECCIÓN CONSTATE
           elseif DD(i)==4
               Di(i)=input('\n DIÁMETRO DEL ELEMENTO [m]:');
               % ELEMENTOS UNICAMENTE NECESARIOS PARA EL CONTEO DURANTE EL USO
               % DE ELEMENTOS DE CUALQUIERA DE LOS DISTINTOS TIPOS POSIBLES DE
               % SECCIONES
               Df(i)=Di(i);
               B(i)=Di(i);
               H(i)=Di(i);
               Bi(i)=B(i);
               Hi(i)=H(i);
               Bf(i)=Bi(i);
               Hf(i)=Hi(i);
               FF44(i)=0;
               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               BET(i)=input('\n VALOR DEL FACTOR DE FORMA (RECTANGULAR=1.2; CIRCULAR=10/9):');
               AREA(i)=pi*Di(i)*Di(i)/4;
               INER(i)=pi*(Di(i)/2)^4/4;
               C1(i)=input('\n\n INGRESE EL VALOR DE C1 [m] (LONGITUD DEL SEGMENTO DE RIGIDEZ INFINITA EN EL NODO INICIAL):');
               C2(i)=input('\n\n INGRESE EL VALOR DE C2 [m] (LONGITUD DEL SEGMENTO DE RIGIDEZ INFINITA EN EL NODO FINAL):');
               F(i)=L(i)-C1(i)-C2(i)
               FI(i)=(3*ELAS*INER(i)*BET(i))/(GCOR*AREA(i)*L(i)^2);
               C=((4*ELAS*INER(i))/L(i))*((1+FI(i))/(1+4*FI(i)));
               CP=C;
               A=((2*ELAS*INER(i))/L(i))*((1-2*FI(i))/(1+4*FI(i)));
               B=((6*ELAS*INER(i))/L(i)^2)*(1/(1+4*FI(i)));
               BP=B;
               T=(12*ELAS*INER(i))/L(i)^3*(1/(1+4*FI(i)));
               R=ELAS*AREA(i)/L(i);
               K2=[R 0 0 -R 0 0;
               0 T (B+C1(i)*T) 0 -T (B+C2(i)*T);
               0 (B+C1(i)*T) (C+2*C1(i)*B+C1(i)^2*T) 0 -(B+C1(i)*T) (A+C1(i)*B+C2(i)*B+C1(i)*C2(i)*T);
               -R 0 0 R 0 0;
               0 -T -(B+C1(i)*T) 0 T -(B+C2(i)*T);
               0 (B+C2(i)*T) (A+C1(i)*B+C2(i)*B+C1(i)*C2(i)*T) 0 -(B+C2(i)*T) (C+2*C2(i)*B+C2(i)^2*T)];
               % MATRIZ DE TRANSPORTE A COORDENADAS GLOBALES
               T2_3=zeros(6,6);
               T2_3(1,1)=COSENO(i);T2_3(1,2)=SENO(i);
               T2_3(2,1)=-SENO(i);T2_3(2,2)=COSENO(i);
               T2_3(3,3)=1;
               T2_3(4,4)=COSENO(i);T2_3(4,5)=SENO(i);
               T2_3(5,4)=-SENO(i);T2_3(5,5)=COSENO(i);
               T2_3(6,6)=1;
               % ELEMENTO CIRCULAR SÓLIDO DE SECCIÓN VARIABLE
           elseif DD(i)==5
               Di(i)=input('\n DIÁMETRO INICIAL DEL ELEMENTO [m]:');
               Df(i)=input('\n DIÁMETRO FINAL DEL ELEMENTO [m]:');
               BET(i)=input('\n VALOR DEL FACTOR DE FORMA (RECTANGULAR=1.2; CIRCULAR=10/9):');
               % ELEMENTOS UNICAMENTE NECESARIOS PARA EL CONTEO DURANTE EL USO
               % DE ELEMENTOS DE CUALQUIERA DE LOS DISTINTOS TIPOS POSIBLES DE
               % SECCIONES
               B(i)=Di(i);
               H(i)=Df(i);
               Bi(i)=B(i);
               Hi(i)=H(i);
               Bf(i)=Bi(i);
               Hf(i)=Hi(i);
               C1(i)=0;
               C2(i)=0;
               F(i)=L(i);
               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               AREA(i)=pi*Di(i)*Di(i)/4;
               INER(i)=pi*(Di(i)/2)^4/4;
               FIx(i)=(12*ELAS*INER(i)*BET(i))/(GCOR*AREA(i)*L(i)^2)
               FIy(i)=FIx(i);
               FI(i)=(1/4)*FIx(i);
               % ELEMENTOS DE LA MATRIZ DE FLEXIBILIDADES
               f11=(4*L(i)/(pi*ELAS*(Di(i))^2))*(Di(i)/(Df(i)-Di(i)))*(1-Di(i)/Df(i));
               f22=(64*(L(i))^3/(3*pi*ELAS*(Di(i))^4))*(Di(i)/Df(i))^3+(40*L(i)/(9*pi*GCOR*(Di(i))^2))*(Di(i)/(Df(i)-Di(i)))*(1-(Di(i)/Df(i)));
               f33=f22;
               f26=(64*(L(i))^2/(3*pi*ELAS*(Di(i))^4))*((Di(i)/Df(i))^3+(1/2)*(Di(i)/(Df(i)-Di(i)))^2*(1+(Di(i)/Df(i))^2-2*(Di(i)/Df(i))));
               f35=f26;
               f66=(64*L(i)/(3*pi*ELAS*(Di(i))^4))*(Di(i)/(Df(i)-Di(i)))*(1-(Di(i)/Df(i))^3);
               f55=f66;
               f44=(32*L(i)/(3*pi*GCOR*Di(i)^4))*(Di(i)/(Df(i)-Di(i)))*(1-(Di(i)/Df(i))^3);
               FF44(i)=f44;
               % ELEMENTOS DE LA MATRIZ DE RIGIDECES
               raz=1/f11;
               Detx=f22*f66-f26*f26;
               r11x=f22/Detx;
               r12x=(f26*L(i)-f22)/Detx;
               r22x=(f66*(L(i))^2-2*f26*L(i)+f22)/Detx;
               raax=(r11x+r22x+2*r12x)/(L(i))^2;
               rabx=(r11x+r12x)/L(i);
               rbax=(r22x+r12x)/L(i);
               % MATRIZ DE RIGIDECES EN COORDENADAS LOCALES
               K11=[raz 0 0 ; 0 raax rabx ; 0 rabx r11x];
               K12=[-raz 0 0 ; 0 -raax rbax ; 0 -rabx r12x];
               K21=K12';
               K22=[raz 0 0 ; 0 raax -rbax ; 0 -rbax r22x];
               K2=[K11 K12 ;K21 K22];
               % MATRIZ DE TRANSPORTE A COORDENADAS GLOBALES
               T2_3=zeros(6,6);
               T2_3(1,1)=COSENO(i);T2_3(1,2)=SENO(i);
               T2_3(2,1)=-SENO(i);T2_3(2,2)=COSENO(i);
               T2_3(3,3)=1;
               T2_3(4,4)=COSENO(i);T2_3(4,5)=SENO(i);
               T2_3(5,4)=-SENO(i);T2_3(5,5)=COSENO(i);
               T2_3(6,6)=1;
           end
       else
           fprintf('\n MIEMBRO %d:',i);
           fprintf('\n\n ESCOJA UNA DE LAS SIGUIENTES OPCIONES:');
           fprintf('\n\n\t\t (1) EL MIEMBRO TIENE LAS MISMAS CARACTERÍSTICAS QUE EL MIEMBRO ANTERIOR');
           fprintf('\n\n\t\t (2) EL MIEMBRO TIENE CARACTERÍSTICAS DISTINTAS AL MIEMBRO ANTERIOR');
           NUB=input('\n\n OPCION : ');
           if NUB==1
               DD(i)=DD(i-1);        
               H(i)=H(i-1);
               B(i)=1;
               Bi(i)=Bi(i-1);
               Hi(i)=Hi(i-1);
               Bf(i)=Bf(i-1);
               Hf(i)=Hf(i-1);
               Di(i)=Di(i-1);
               Df(i)=Df(i-1);
               BET(i)=BET(i-1);
               C1(i)=C1(i-1);
               C2(i)=C2(i-1);
               F(i)=F(i-1);
               FI(i)=FI(i-1);
               INER(i)=INER(i-1);
               AREA(i)=AREA(i-1);
               FF44(i)=FF44(i-1);
               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               if DD(i)==1               
                   C=((4*ELAS*INER(i))/L(i))*((1+FI(i))/(1+4*FI(i)));
                   CP=C;
                   A=((2*ELAS*INER(i))/L(i))*((1-2*FI(i))/(1+4*FI(i)));
                   B=((6*ELAS*INER(i))/L(i)^2)*(1/(1+4*FI(i)));
                   BP=B;
                   T=(12*ELAS*INER(i))/L(i)^3*(1/(1+4*FI(i)));
                   R=ELAS*AREA(i)/L(i);
                   K2=[R 0 0 -R 0 0;
                   0 T (B+C1(i)*T) 0 -T (B+C2(i)*T);
                   0 (B+C1(i)*T) (C+2*C1(i)*B+C1(i)^2*T) 0 -(B+C1(i)*T) (A+C1(i)*B+C2(i)*B+C1(i)*C2(i)*T);
                   -R 0 0 R 0 0;
                   0 -T -(B+C1(i)*T) 0 T -(B+C2(i)*T);
                   0 (B+C2(i)*T) (A+C1(i)*B+C2(i)*B+C1(i)*C2(i)*T) 0 -(B+C2(i)*T) (C+2*C2(i)*B+C2(i)^2*T)];
                   % MATRIZ DE TRANSPORTE A COORDENADAS GLOBALES
                   T2_3=zeros(6,6);
                   T2_3(1,1)=COSENO(i);T2_3(1,2)=SENO(i);
                   T2_3(2,1)=-SENO(i);T2_3(2,2)=COSENO(i);
                   T2_3(3,3)=1;
                   T2_3(4,4)=COSENO(i);T2_3(4,5)=SENO(i);
                   T2_3(5,4)=-SENO(i);T2_3(5,5)=COSENO(i);
                   T2_3(6,6)=1;
                   % ELEMENTO CUADRADO SÓLIDO DE SECCIÓN VARIABLE
               elseif DD(i)==2
                   % ELEMENTOS DE LA MATRIZ DE FLEXIBILIDADES
                   f11=(L(i)/(ELAS*Hi(i)*Hi(i)))*(Hi(i)/(Hf(i)-Hi(i)))*(1-Hi(i)/Hf(i));
                   f22=(4*(L(i))^3/(ELAS*(Hi(i))^4))*(Hi(i)/Hf(i))^3+(6*L(i)/(5*GCOR*(Hi(i))^2))*(Hi(i)/(Hf(i)-Hi(i)))*(1-Hi(i)/Hf(i));
                   f33=f22;
                   f26=(4*(L(i))^2/(ELAS*(Hi(i))^4))*((Hi(i)/Hf(i))^3+(1/2)*(Hi(i)/(Hf(i)-Hi(i)))^2*(1+(Hi(i)/Hf(i))^2-2*(Hi(i)/Hf(i))));
                   f35=f26;
                   f66=(4*L(i)/(ELAS*(Hi(i))^4))*(Hi(i)/(Hf(i)-Hi(i)))*(1-(Hi(i)/Hf(i))^3);
                   f55=f66;
                   f44=(64*L(i)/(27*GCOR*(Hi(i))^4))*(Hi(i)/(Hf(i)-Hi(i)))*(1-(Hi(i)/Hf(i))^3);
                   % ELEMENTOS DE LA MATRIZ DE RIGIDECES
                   raz=1/f11;
                   Detx=f22*f66-f26*f26;
                   r11x=f22/Detx;
                   r12x=(f26*L(i)-f22)/Detx;
                   r22x=(f66*(L(i))^2-2*f26*L(i)+f22)/Detx;
                   raax=(r11x+r22x+2*r12x)/(L(i))^2;
                   rabx=(r11x+r12x)/L(i);
                   rbax=(r22x+r12x)/L(i);
                   % MATRIZ DE RIGIDECES EN COORDENADAS LOCALES
                   K11=[raz 0 0 ; 0 raax rabx ; 0 rabx r11x];
                   K12=[-raz 0 0 ; 0 -raax rbax ; 0 -rabx r12x];
                   K21=K12';
                   K22=[raz 0 0 ; 0 raax -rbax ; 0 -rbax r22x];
                   K2=[K11 K12 ;K21 K22];
                   % MATRIZ DE TRANSPORTE A COORDENADAS GLOBALES
                   T2_3=zeros(6,6);
                   T2_3(1,1)=COSENO(i);T2_3(1,2)=SENO(i);
                   T2_3(2,1)=-SENO(i);T2_3(2,2)=COSENO(i);
                   T2_3(3,3)=1;
                   T2_3(4,4)=COSENO(i);T2_3(4,5)=SENO(i);
                   T2_3(5,4)=-SENO(i);T2_3(5,5)=COSENO(i);
                   T2_3(6,6)=1;
                   % ELEMENTO RECTANGULAR SÓLIDO DE SECCIÓN VARIABLE
               elseif DD(i)==3
                   Bb=Bi(i);
                   % ELEMENTOS DE LA MATRIZ DE FLEXIBILIDADES
                   f11=(L(i)/(ELAS*Bb*Hi(i)))*(Hi(i)/(Hf(i)-Hi(i)))*log(Hf(i)/Hi(i));
                   f22=(6*(L(i))^3/(ELAS*Bb*(Hi(i))^3))*((Hi(i)/Hf(i))^2-(Hi(i)/(Hf(i)-Hi(i)))^3*(Hf(i)/Hi(i)-Hi(i)/Hf(i)-2*log(Hf(i)/Hi(i))))+(6*L(i)/(5*GCOR*Bb*Hi(i)))*(Hi(i)/(Hf(i)-Hi(i)))*log(Hf(i)/Hi(i));
                   f33=(12*(L(i))^3/(ELAS*(Bb)^3*Hi(i)))*(3/2+(1/2)*(Hf(i)/Hi(i))^2-2*(Hf(i)/Hi(i))+log(Hf(i)/Hi(i)))+(6*L(i)/(5*GCOR*Bb*Hi(i)))*(Hi(i)/(Hf(i)-Hi(i)))*log(Hf(i)/Hi(i));
                   f26=(6*L(i)/(ELAS*Bb*(Hi(i))^3))*(Hi(i)/(Hf(i)-Hi(i)))^2*(1+(Hi(i)/Hf(i))^2-2*(Hi(i)/Hf(i)));
                   f35=(12*(L(i))^2/(ELAS*(Bb)^3*Hi(i)))*(Hi(i)/(Hf(i)-Hi(i)))^2*(Hf(i)/Hi(i)-log(Hf(i)/Hi(i))-1);
                   f66=(6*L(i)/(ELAS*Bb*(Hi(i))^3))*(Hi(i)/(Hf(i)-Hi(i)))*(1-(Hi(i)/Hf(i))^2);
                   f55=(12*L(i)/(ELAS*(Bb)^3*Hi(i)))*(Hi(i)/(Hf(i)-Hi(i)))*log(Hf(i)/Hi(i));
                   f44=FF44(i);
                   % ELEMENTOS DE LA MATRIZ DE RIGIDECES
                   raz=1/f11;
                   Detx=f22*f66-f26*f26;
                   r11x=f22/Detx;
                   r12x=(f26*L(i)-f22)/Detx;
                   r22x=(f66*(L(i))^2-2*f26*L(i)+f22)/Detx;
                   raax=(r11x+r22x+2*r12x)/(L(i))^2;
                   rabx=(r11x+r12x)/L(i);
                   rbax=(r22x+r12x)/L(i);
                   % MATRIZ DE RIGIDECES EN COORDENADAS LOCALES
                   K11=[raz 0 0 ; 0 raax rabx ; 0 rabx r11x];
                   K12=[-raz 0 0 ; 0 -raax rbax ; 0 -rabx r12x];
                   K21=K12';
                   K22=[raz 0 0 ; 0 raax -rbax ; 0 -rbax r22x];
                   K2=[K11 K12 ;K21 K22];
                   % MATRIZ DE TRANSPORTE A COORDENADAS GLOBALES
                   T2_3=zeros(6,6);
                   T2_3(1,1)=COSENO(i);T2_3(1,2)=SENO(i);
                   T2_3(2,1)=-SENO(i);T2_3(2,2)=COSENO(i);
                   T2_3(3,3)=1;
                   T2_3(4,4)=COSENO(i);T2_3(4,5)=SENO(i);
                   T2_3(5,4)=-SENO(i);T2_3(5,5)=COSENO(i);
                   T2_3(6,6)=1;
                   % ELEMENTO CIRCULAR SÓLIDO DE SECCIÓN CONSTATE
               elseif DD(i)==4
                   C=((4*ELAS*INER(i))/L(i))*((1+FI(i))/(1+4*FI(i)));
                   CP=C;
                   A=((2*ELAS*INER(i))/L(i))*((1-2*FI(i))/(1+4*FI(i)));
                   B=((6*ELAS*INER(i))/L(i)^2)*(1/(1+4*FI(i)));
                   BP=B;
                   T=(12*ELAS*INER(i))/L(i)^3*(1/(1+4*FI(i)));
                   R=ELAS*AREA(i)/L(i);
                   K2=[R 0 0 -R 0 0;
                   0 T (B+C1(i)*T) 0 -T (B+C2(i)*T);
                   0 (B+C1(i)*T) (C+2*C1(i)*B+C1(i)^2*T) 0 -(B+C1(i)*T) (A+C1(i)*B+C2(i)*B+C1(i)*C2(i)*T);
                   -R 0 0 R 0 0;
                   0 -T -(B+C1(i)*T) 0 T -(B+C2(i)*T);
                   0 (B+C2(i)*T) (A+C1(i)*B+C2(i)*B+C1(i)*C2(i)*T) 0 -(B+C2(i)*T) (C+2*C2(i)*B+C2(i)^2*T)];
                   % MATRIZ DE TRANSPORTE A COORDENADAS GLOBALES
                   T2_3=zeros(6,6);
                   T2_3(1,1)=COSENO(i);T2_3(1,2)=SENO(i);
                   T2_3(2,1)=-SENO(i);T2_3(2,2)=COSENO(i);
                   T2_3(3,3)=1;
                   T2_3(4,4)=COSENO(i);T2_3(4,5)=SENO(i);
                   T2_3(5,4)=-SENO(i);T2_3(5,5)=COSENO(i);
                   T2_3(6,6)=1;
                   % ELEMENTO CIRCULAR SÓLIDO DE SECCIÓN VARIABLE
               elseif DD(i)==5
                   % ELEMENTOS DE LA MATRIZ DE FLEXIBILIDADES
                   f11=(4*L(i)/(pi*ELAS*(Di(i))^2))*(Di(i)/(Df(i)-Di(i)))*(1-Di(i)/Df(i));
                   f22=(64*(L(i))^3/(3*pi*ELAS*(Di(i))^4))*(Di(i)/Df(i))^3+(40*L(i)/(9*pi*GCOR*(Di(i))^2))*(Di(i)/(Df(i)-Di(i)))*(1-(Di(i)/Df(i)));
                   f33=f22;
                   f26=(64*(L(i))^2/(3*pi*ELAS*(Di(i))^4))*((Di(i)/Df(i))^3+(1/2)*(Di(i)/(Df(i)-Di(i)))^2*(1+(Di(i)/Df(i))^2-2*(Di(i)/Df(i))));
                   f35=f26;
                   f66=(64*L(i)/(3*pi*ELAS*(Di(i))^4))*(Di(i)/(Df(i)-Di(i)))*(1-(Di(i)/Df(i))^3);
                   f55=f66;
                   f44=(32*L(i)/(3*pi*GCOR*Di(i)^4))*(Di(i)/(Df(i)-Di(i)))*(1-(Di(i)/Df(i))^3);
                   % ELEMENTOS DE LA MATRIZ DE RIGIDECES
                   raz=1/f11;
                   Detx=f22*f66-f26*f26;
                   r11x=f22/Detx;
                   r12x=(f26*L(i)-f22)/Detx;
                   r22x=(f66*(L(i))^2-2*f26*L(i)+f22)/Detx;
                   raax=(r11x+r22x+2*r12x)/(L(i))^2;
                   rabx=(r11x+r12x)/L(i);
                   rbax=(r22x+r12x)/L(i);
                   % MATRIZ DE RIGIDECES EN COORDENADAS LOCALES
                   K11=[raz 0 0 ; 0 raax rabx ; 0 rabx r11x];
                   K12=[-raz 0 0 ; 0 -raax rbax ; 0 -rabx r12x];
                   K21=K12';
                   K22=[raz 0 0 ; 0 raax -rbax ; 0 -rbax r22x];
                   K2=[K11 K12 ;K21 K22];
                   % MATRIZ DE TRANSPORTE A COORDENADAS GLOBALES
                   T2_3=zeros(6,6);
                   T2_3(1,1)=COSENO(i);T2_3(1,2)=SENO(i);
                   T2_3(2,1)=-SENO(i);T2_3(2,2)=COSENO(i);
                   T2_3(3,3)=1;
                   T2_3(4,4)=COSENO(i);T2_3(4,5)=SENO(i);
                   T2_3(5,4)=-SENO(i);T2_3(5,5)=COSENO(i);
                   T2_3(6,6)=1;
               end
           elseif NUB==2
               fprintf('\n MIEMBRO %d:',i);
               fprintf('\n\n ESCOJA UNA DE LAS SIGUIENTES OPCIONES:');
               fprintf('\n\n\t\t (1) MIEMBRO PRISMÁTICO DE SECCIÓN RECTANGULAR');
               fprintf('\n\n\t\t (2) MIEMBRO CUADRADO DE SECCIÓN VARIABLE');
               fprintf('\n\n\t\t (3) MIEMBRO RECTANGULAR DE SECCIÓN VARIABLE');
               fprintf('\n\n\t\t (4) MIEMBRO CIRCULAR DE SECCIÓN PRISMÁTICA CONSTANTE');
               fprintf('\n\n\t\t (5) MIEMBRO CIRCULAR DE SECCIÓN VARIABLE');
               DD(i) = input('\n\n OPCION : ');
               fprintf('\n\n\n'); 
               % ELEMENTO PRISMÁTICO DE SECCIÓN CONSTANTE DE FORMA RECTANGULAR
               if DD(i)==1
                   B(i)=input('\n BASE DEL ELEMENTO [m]:');
                   H(i)=input(' ALTURA DEL ELEMENTO [m]:');
                   % ELEMENTOS UNICAMENTE NECESARIOS PARA EL CONTEO DURANTE EL USO
                   % DE ELEMENTOS DE CUALQUIERA DE LOS DISTINTOS TIPOS POSIBLES DE
                   % SECCIONES
                   Bi(i)=B(i);
                   Hi(i)=H(i);
                   Bf(i)=Bi(i);
                   Hf(i)=Hi(i);
                   Di(i)=B(i);
                   Df(i)=Di(i);
                   FF44(i)=0;
                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                   BET(i)=input('\n VALOR DEL FACTOR DE FORMA (RECTANGULAR=1.2; CIRCULAR=10/9):');
                   C1(i)=input('\n\n INGRESE EL VALOR DE C1 [m] (LONGITUD DEL SEGMENTO DE RIGIDEZ INFINITA EN EL NODO INICIAL):');
                   C2(i)=input('\n\n INGRESE EL VALOR DE C2 [m] (LONGITUD DEL SEGMENTO DE RIGIDEZ INFINITA EN EL NODO FINAL):');
                   F(i)=L(i)-C1(i)-C2(i)
                   % CONSTANTES AUXILIARES PARA FUTUROS DESARROLLOS EN 3D
                   AREA(i)=B(i)*H(i);
                   INERx(i)=(B(i)*H(i)^3)/12;
                   INER(i)=INERx(i);
                   INERy(i)=(H(i)*B(i)^3)/12;
                   FIx(i)=(12*ELAS*INERx(i)*BET(i))/(GCOR*AREA(i)*L(i)^2);
                   FIy(i)=(12*ELAS*INERy(i)*BET(i))/(GCOR*AREA(i)*L(i)^2);
                   % GENERACIÓN DE LAS VARIABLES
                   FI(i)=(1/4)*FIx(i);
                   C=((4*ELAS*INER(i))/L(i))*((1+FI(i))/(1+4*FI(i)));
                   CP=C;
                   A=((2*ELAS*INER(i))/L(i))*((1-2*FI(i))/(1+4*FI(i)));
                   B=((6*ELAS*INER(i))/L(i)^2)*(1/(1+4*FI(i)));
                   BP=B;
                   T=(12*ELAS*INER(i))/L(i)^3*(1/(1+4*FI(i)));
                   R=ELAS*AREA(i)/L(i);
                   % MATRIZ DE RIGIDECES EN COORDENADAS LOCALES
                   K2=[R 0 0 -R 0 0;
                   0 T (B+C1(i)*T) 0 -T (B+C2(i)*T);
                   0 (B+C1(i)*T) (C+2*C1(i)*B+C1(i)^2*T) 0 -(B+C1(i)*T) (A+C1(i)*B+C2(i)*B+C1(i)*C2(i)*T);
                   -R 0 0 R 0 0;
                   0 -T -(B+C1(i)*T) 0 T -(B+C2(i)*T);
                   0 (B+C2(i)*T) (A+C1(i)*B+C2(i)*B+C1(i)*C2(i)*T) 0 -(B+C2(i)*T) (C+2*C2(i)*B+C2(i)^2*T)];
                   % MATRIZ DE TRANSPORTE A COORDENADAS GLOBALES
                   T2_3=zeros(6,6);
                   T2_3(1,1)=COSENO(i);T2_3(1,2)=SENO(i);
                   T2_3(2,1)=-SENO(i);T2_3(2,2)=COSENO(i);
                   T2_3(3,3)=1;
                   T2_3(4,4)=COSENO(i);T2_3(4,5)=SENO(i);
                   T2_3(5,4)=-SENO(i);T2_3(5,5)=COSENO(i);
                   T2_3(6,6)=1;
                   % ELEMENTO CUADRADO SÓLIDO DE SECCIÓN VARIABLE
               elseif DD(i)==2
                   Bi(i)=input('\n BASE INICIAL [m]:');
                   Hi(i)=Bi(i);
                   Bf(i)=input('\n BASE FINAL [m]:');
                   % ELEMENTOS UNICAMENTE NECESARIOS PARA EL CONTEO DURANTE EL USO
                   % DE ELEMENTOS DE CUALQUIERA DE LOS DISTINTOS TIPOS POSIBLES DE
                   % SECCIONES
                   Hf(i)=Bf(i);
                   B(i)=Bi(i);
                   H(i)=Hi(i);
                   Di(i)=Bi(i);
                   Df(i)=Di(i);
                   BET(i)=1.2;
                   INER(i)=(B(i)*H(i)^3)/12;
                   FI(i)=1;
                   AREA(i)=B(i)*H(i);
                   C1(i)=C1(i-1);
                   C2(i)=C2(i-1);
                   F(i)=F(i-1);
                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                   % ELEMENTOS DE LA MATRIZ DE FLEXIBILIDADES
                   f11=(L(i)/(ELAS*Hi(i)*Hi(i)))*(Hi(i)/(Hf(i)-Hi(i)))*(1-Hi(i)/Hf(i));
                   f22=(4*(L(i))^3/(ELAS*(Hi(i))^4))*(Hi(i)/Hf(i))^3+(6*L(i)/(5*GCOR*(Hi(i))^2))*(Hi(i)/(Hf(i)-Hi(i)))*(1-Hi(i)/Hf(i));
                   f33=f22;
                   f26=(4*(L(i))^2/(ELAS*(Hi(i))^4))*((Hi(i)/Hf(i))^3+(1/2)*(Hi(i)/(Hf(i)-Hi(i)))^2*(1+(Hi(i)/Hf(i))^2-2*(Hi(i)/Hf(i))));
                   f35=f26;
                   f66=(4*L(i)/(ELAS*(Hi(i))^4))*(Hi(i)/(Hf(i)-Hi(i)))*(1-(Hi(i)/Hf(i))^3);
                   f55=f66;
                   f44=(64*L(i)/(27*GCOR*(Hi(i))^4))*(Hi(i)/(Hf(i)-Hi(i)))*(1-(Hi(i)/Hf(i))^3);
                   FF44(i)=f44;
                   % ELEMENTOS DE LA MATRIZ DE RIGIDECES
                   raz=1/f11;
                   Detx=f22*f66-f26*f26;
                   r11x=f22/Detx;
                   r12x=(f26*L(i)-f22)/Detx;
                   r22x=(f66*(L(i))^2-2*f26*L(i)+f22)/Detx;
                   raax=(r11x+r22x+2*r12x)/(L(i))^2;
                   rabx=(r11x+r12x)/L(i);
                   rbax=(r22x+r12x)/L(i);
                   % MATRIZ DE RIGIDECES EN COORDENADAS LOCALES
                   K11=[raz 0 0 ; 0 raax rabx ; 0 rabx r11x];
                   K12=[-raz 0 0 ; 0 -raax rbax ; 0 -rabx r12x];
                   K21=K12';
                   K22=[raz 0 0 ; 0 raax -rbax ; 0 -rbax r22x];
                   K2=[K11 K12 ;K21 K22];
                   % MATRIZ DE TRANSPORTE A COORDENADAS GLOBALES
                   T2_3=zeros(6,6);
                   T2_3(1,1)=COSENO(i);T2_3(1,2)=SENO(i);
                   T2_3(2,1)=-SENO(i);T2_3(2,2)=COSENO(i);
                   T2_3(3,3)=1;
                   T2_3(4,4)=COSENO(i);T2_3(4,5)=SENO(i);
                   T2_3(5,4)=-SENO(i);T2_3(5,5)=COSENO(i);
                   T2_3(6,6)=1;
                   % ELEMENTO RECTANGULAR SÓLIDO DE SECCIÓN VARIABLE
               elseif DD(i)==3
                   Bi(i)=input('\n BASE DE TODO EL ELEMENTO [m]:');
                   Hi(i)=input(' ALTURA INICIAL DEL ELEMENTO [m]:');
                   Bf(i)=Bi(i);
                   Hf(i)=input(' ALTURA FINAL DEL ELEMENTO [m]:');
                   Bb=Bi(i);
                   % ELEMENTOS UNICAMENTE NECESARIOS PARA EL CONTEO DURANTE EL USO
                   % DE ELEMENTOS DE CUALQUIERA DE LOS DISTINTOS TIPOS POSIBLES DE
                   % SECCIONES
                   B(i)=Bi(i);
                   H(i)=Hi(i);
                   Di(i)=Bi(i);
                   Df(i)=Di(i);
                   INER(i)=(B(i)*H(i)^3)/12;
                   FI(i)=1;
                   AREA(i)=B(i)*H(i);
                   BET(i)=1.2;
                   C1(i)=0;
                   C2(i)=0;
                   F(i)=L(i);
                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                   BET(i)=input('\n FACTOR DE FORMA (RECTANGULAR=1.2; CIRCULAR=10/9):');
                   fprintf('\n\t\t\t ‘ ; ......................................................................................................................... ');
                   fprintf('\n\t\t\t ‘ ; ......................................................................................................................... ');
                   fprintf('\n\n\t\t     LA FUNCIÓN F(x) INDICA LA VARIACIÓN CON LA LONGITUD DE LA FORMA DEL ELEMENTO,');
                   fprintf('\n\t\t\t                                                                    ');
                   fprintf('\n\n\t\t     SI LA FUNCIÓN ES LINEAL: Y=F(x)=');
                   fprintf('\n\t\t\t                                      (abs(Hi-Hf)/LL)*x                                       ');
                   fprintf('\n\t\t\t ‘ ; ......................................................................................................................... ');
                   fprintf('\n\t\t\t ‘ ; ......................................................................................................................... ');
                   fprintf('\n\t\t\t ‘ ;                                                                                                                           ');
                   fprintf('\n\t\t\t   ESTANDO EL EJE X SOBRE EL     Y|                                                                                            ');
                   fprintf('\n\t\t\t    CENTROIDE DEL MIEMBRO         |                                                                                            ');
                   fprintf('\n\t\t\t                                  |                                                                                            ');
                   fprintf('\n\t\t\t                                  |_ _ _ _ _                                                                                   ');
                   fprintf('\n\t\t\t                                           X                                                                                   ');
                   fprintf('\n\t\t\t                                                                                                                               ');
                   fprintf('\n\t\t\t ‘ ; ......................................................................................................................... ');
                   fprintf('\n\t\t\t ‘ ; ......................................................................................................................... ');
                   fprintf('\n\t\t\t                                                                                                                              ');
                   fprintf('\n\t\t\t                                                                                                                              ');
                   EE(i)=menu('ESCOJA UNA DE LAS SIGUIENTES OPCIONES:','FUNCIÓN LINEAL','OTRO TIPO DE VARIACIÓN');
                   if EE(i)==1
                       % ELEMENTOS DE LA MATRIZ DE FLEXIBILIDADES
                       f11=(L(i)/(ELAS*Bb*Hi(i)))*(Hi(i)/(Hf(i)-Hi(i)))*log(Hf(i)/Hi(i));
                       f22=(6*(L(i))^3/(ELAS*Bb*(Hi(i))^3))*((Hi(i)/Hf(i))^2-(Hi(i)/(Hf(i)-Hi(i)))^3*(Hf(i)/Hi(i)-Hi(i)/Hf(i)-2*log(Hf(i)/Hi(i))))+(6*L(i)/(5*GCOR*Bb*Hi(i)))*(Hi(i)/(Hf(i)-Hi(i)))*log(Hf(i)/Hi(i));
                       f33=(12*(L(i))^3/(ELAS*(Bb)^3*Hi(i)))*(3/2+(1/2)*(Hf(i)/Hi(i))^2-2*(Hf(i)/Hi(i))+log(Hf(i)/Hi(i)))+(6*L(i)/(5*GCOR*Bb*Hi(i)))*(Hi(i)/(Hf(i)-Hi(i)))*log(Hf(i)/Hi(i));
                       f26=(6*L(i)/(ELAS*Bb*(Hi(i))^3))*(Hi(i)/(Hf(i)-Hi(i)))^2*(1+(Hi(i)/Hf(i))^2-2*(Hi(i)/Hf(i)));
                       f35=(12*(L(i))^2/(ELAS*(Bb)^3*Hi(i)))*(Hi(i)/(Hf(i)-Hi(i)))^2*(Hf(i)/Hi(i)-log(Hf(i)/Hi(i))-1);
                       f66=(6*L(i)/(ELAS*Bb*(Hi(i))^3))*(Hi(i)/(Hf(i)-Hi(i)))*(1-(Hi(i)/Hf(i))^2);
                       f55=(12*L(i)/(ELAS*(Bb)^3*Hi(i)))*(Hi(i)/(Hf(i)-Hi(i)))*log(Hf(i)/Hi(i));
                       % DEFINICIÓN DE UNA FUNCIÓN DE x
                       Hi=Hi(i);
                       Hf=Hf(i);
                       LL=L(i);
                       if Hi<Hf
                           dpc= @(x)(1./(((Bb.^3*(-Hi./2-((Hi-Hf)./(2*LL))*x))./3)*(1-(192*Bb./(pi.^5*(-Hi./2-((Hi-Hf)./(2*LL))*x)))*(tanh(pi*(-Hi./2-((Hi-Hf)./(2*LL))*x)./(2*Bb)))))); 
                        elseif Hi>Hf
                            dpc= @(x)(1./(((Bb.^3*(-Hi./2+((Hi-Hf)./(2*LL))*x))./3)*(1-(192*Bb./(pi.^5*(-Hi./2+((Hi-Hf)./(2*LL))*x)))*(tanh(pi*(-Hi./2+((Hi-Hf)./(2*LL))*x)./(2*Bb))))));
                       end
                       f44=quad(@(x)dpc, 0, LL)
                       FF44(i)=f44;
                       % ELEMENTOS DE LA MATRIZ DE RIGIDECES
                       raz=1/f11;
                       Detx=f22*f66-f26*f26;
                       r11x=f22/Detx;
                       r12x=(f26*L(i)-f22)/Detx;
                       r22x=(f66*(L(i))^2-2*f26*L(i)+f22)/Detx;
                       raax=(r11x+r22x+2*r12x)/(L(i))^2;
                       rabx=(r11x+r12x)/L(i);
                       rbax=(r22x+r12x)/L(i);
                       % MATRIZ DE RIGIDECES EN COORDENADAS LOCALES
                       K11=[raz 0 0 ; 0 raax rabx ; 0 rabx r11x];
                       K12=[-raz 0 0 ; 0 -raax rbax ; 0 -rabx r12x];
                       K21=K12';
                       K22=[raz 0 0 ; 0 raax -rbax ; 0 -rbax r22x];
                       K2=[K11 K12 ;K21 K22];
                       % MATRIZ DE TRANSPORTE A COORDENADAS GLOBALES
                       T2_3=zeros(6,6);
                       T2_3(1,1)=COSENO(i);T2_3(1,2)=SENO(i);
                       T2_3(2,1)=-SENO(i);T2_3(2,2)=COSENO(i);
                       T2_3(3,3)=1;
                       T2_3(4,4)=COSENO(i);T2_3(4,5)=SENO(i);
                       T2_3(5,4)=-SENO(i);T2_3(5,5)=COSENO(i);
                       T2_3(6,6)=1;
                   else
                       % ELEMENTOS DE LA MATRIZ DE FLEXIBILIDADES
                       f11=(L(i)/(ELAS*Bb*Hi(i)))*(Hi(i)/(Hf(i)-Hi(i)))*log(Hf(i)/Hi(i));
                       f22=(6*(L(i))^3/(ELAS*Bb*(Hi(i))^3))*((Hi(i)/Hf(i))^2-(Hi(i)/(Hf(i)-Hi(i)))^3*(Hf(i)/Hi(i)-Hi(i)/Hf(i)-2*log(Hf(i)/Hi(i))))+(6*L(i)/(5*GCOR*Bb*Hi(i)))*(Hi(i)/(Hf(i)-Hi(i)))*log(Hf(i)/Hi(i));
                       f33=(12*(L(i))^3/(ELAS*(Bb)^3*Hi(i)))*(3/2+(1/2)*(Hf(i)/Hi(i))^2-2*(Hf(i)/Hi(i))+log(Hf(i)/Hi(i)))+(6*L(i)/(5*GCOR*Bb*Hi(i)))*(Hi(i)/(Hf(i)-Hi(i)))*log(Hf(i)/Hi(i));
                       f26=(6*L(i)/(ELAS*Bb*(Hi(i))^3))*(Hi(i)/(Hf(i)-Hi(i)))^2*(1+(Hi(i)/Hf(i))^2-2*(Hi(i)/Hf(i)));
                       f35=(12*(L(i))^2/(ELAS*(Bb)^3*Hi(i)))*(Hi(i)/(Hf(i)-Hi(i)))^2*(Hf(i)/Hi(i)-log(Hf(i)/Hi(i))-1);
                       f66=(6*L(i)/(ELAS*Bb*(Hi(i))^3))*(Hi(i)/(Hf(i)-Hi(i)))*(1-(Hi(i)/Hf(i))^2);
                       f55=(12*L(i)/(ELAS*(Bb)^3*Hi(i)))*(Hi(i)/(Hf(i)-Hi(i)))*log(Hf(i)/Hi(i));
                       % DEFINICIÓN DE UNA FUNCIÓN DE x
                       Hi=Hi(i);
                       Hf=Hf(i);
                       LL=L(i);
                       f44=input(' INTRODUZCA LA f44=');
                       FF44(i)=f44;
                       % ELEMENTOS DE LA MATRIZ DE RIGIDECES
                       raz=1/f11;
                       Detx=f22*f66-f26*f26;
                       r11x=f22/Detx;
                       r12x=(f26*L(i)-f22)/Detx;
                       r22x=(f66*(L(i))^2-2*f26*L(i)+f22)/Detx;
                       raax=(r11x+r22x+2*r12x)/(L(i))^2;
                       rabx=(r11x+r12x)/L(i);
                       rbax=(r22x+r12x)/L(i);
                       % MATRIZ DE RIGIDECES EN COORDENADAS LOCALES
                       K11=[raz 0 0 ; 0 raax rabx ; 0 rabx r11x];
                       K12=[-raz 0 0 ; 0 -raax rbax ; 0 -rabx r12x];
                       K21=K12';
                       K22=[raz 0 0 ; 0 raax -rbax ; 0 -rbax r22x];
                       K2=[K11 K12 ;K21 K22];
                       % MATRIZ DE TRANSPORTE A COORDENADAS GLOBALES
                       T2_3=zeros(6,6);
                       T2_3(1,1)=COSENO(i);T2_3(1,2)=SENO(i);
                       T2_3(2,1)=-SENO(i);T2_3(2,2)=COSENO(i);
                       T2_3(3,3)=1;
                       T2_3(4,4)=COSENO(i);T2_3(4,5)=SENO(i);
                       T2_3(5,4)=-SENO(i);T2_3(5,5)=COSENO(i);
                       T2_3(6,6)=1;
                   end
                  % ELEMENTO CIRCULAR SÓLIDO DE SECCIÓN CONSTATE
               elseif DD(i)==4
                   Di(i)=input('\n DIÁMETRO DEL ELEMENTO [m]:');
                   % ELEMENTOS UNICAMENTE NECESARIOS PARA EL CONTEO DURANTE EL USO
                   % DE ELEMENTOS DE CUALQUIERA DE LOS DISTINTOS TIPOS POSIBLES DE
                   % SECCIONES
                   Df(i)=Di(i);
                   B(i)=Di(i);
                   H(i)=Di(i);
                   Bi(i)=B(i);
                   Hi(i)=H(i);
                   Bf(i)=Bi(i);
                   Hf(i)=Hi(i);
                   FF44(i)=0;
                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                   BET(i)=input('\n VALOR DEL FACTOR DE FORMA (RECTANGULAR=1.2; CIRCULAR=10/9):');
                   AREA(i)=pi*Di(i)*Di(i)/4;
                   INER(i)=pi*(Di(i)/2)^4/4;
                   C1(i)=input('\n\n INGRESE EL VALOR DE C1 [m] (LONGITUD DEL SEGMENTO DE RIGIDEZ INFINITA EN EL NODO INICIAL):');
                   C2(i)=input('\n\n INGRESE EL VALOR DE C2 [m] (LONGITUD DEL SEGMENTO DE RIGIDEZ INFINITA EN EL NODO FINAL):');
                   F(i)=L(i)-C1(i)-C2(i)
                   FI(i)=(3*ELAS*INER(i)*BET(i))/(GCOR*AREA(i)*L(i)^2);
                   C=((4*ELAS*INER(i))/L(i))*((1+FI(i))/(1+4*FI(i)));
                   CP=C;
                   A=((2*ELAS*INER(i))/L(i))*((1-2*FI(i))/(1+4*FI(i)));
                   B=((6*ELAS*INER(i))/L(i)^2)*(1/(1+4*FI(i)));
                   BP=B;
                   T=(12*ELAS*INER(i))/L(i)^3*(1/(1+4*FI(i)));
                   R=ELAS*AREA(i)/L(i);
                   K2=[R 0 0 -R 0 0;
                   0 T (B+C1(i)*T) 0 -T (B+C2(i)*T);
                   0 (B+C1(i)*T) (C+2*C1(i)*B+C1(i)^2*T) 0 -(B+C1(i)*T) (A+C1(i)*B+C2(i)*B+C1(i)*C2(i)*T);
                   -R 0 0 R 0 0;
                   0 -T -(B+C1(i)*T) 0 T -(B+C2(i)*T);
                   0 (B+C2(i)*T) (A+C1(i)*B+C2(i)*B+C1(i)*C2(i)*T) 0 -(B+C2(i)*T) (C+2*C2(i)*B+C2(i)^2*T)];
                   % MATRIZ DE TRANSPORTE A COORDENADAS GLOBALES
                   T2_3=zeros(6,6);
                   T2_3(1,1)=COSENO(i);T2_3(1,2)=SENO(i);
                   T2_3(2,1)=-SENO(i);T2_3(2,2)=COSENO(i);
                   T2_3(3,3)=1;
                   T2_3(4,4)=COSENO(i);T2_3(4,5)=SENO(i);
                   T2_3(5,4)=-SENO(i);T2_3(5,5)=COSENO(i);
                   T2_3(6,6)=1;
                   % ELEMENTO CIRCULAR SÓLIDO DE SECCIÓN VARIABLE
               elseif DD(i)==5
                   Di(i)=input('\n DIÁMETRO INICIAL DEL ELEMENTO [m]:');
                   Df(i)=input('\n DIÁMETRO FINAL DEL ELEMENTO [m]:');
                   BET(i)=input('\n VALOR DEL FACTOR DE FORMA (RECTANGULAR=1.2; CIRCULAR=10/9):');
                   % ELEMENTOS UNICAMENTE NECESARIOS PARA EL CONTEO DURANTE EL USO
                   % DE ELEMENTOS DE CUALQUIERA DE LOS DISTINTOS TIPOS POSIBLES DE
                   % SECCIONES
                   B(i)=Di(i);
                   H(i)=Df(i);
                   Bi(i)=B(i);
                   Hi(i)=H(i);
                   Bf(i)=Bi(i);
                   Hf(i)=Hi(i);
                   C1(i)=0;
                   C2(i)=0;
                   F(i)=L(i);
                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                   AREA(i)=pi*Di(i)*Di(i)/4;
                   INER(i)=pi*(Di(i)/2)^4/4;
                   FIx(i)=(12*ELAS*INER(i)*BET(i))/(GCOR*AREA(i)*L(i)^2)
                   FIy(i)=FIx(i);
                   FI(i)=(1/4)*FIx(i);
                   % ELEMENTOS DE LA MATRIZ DE FLEXIBILIDADES
                   f11=(4*L(i)/(pi*ELAS*(Di(i))^2))*(Di(i)/(Df(i)-Di(i)))*(1-Di(i)/Df(i));
                   f22=(64*(L(i))^3/(3*pi*ELAS*(Di(i))^4))*(Di(i)/Df(i))^3+(40*L(i)/(9*pi*GCOR*(Di(i))^2))*(Di(i)/(Df(i)-Di(i)))*(1-(Di(i)/Df(i)));
                   f33=f22;
                   f26=(64*(L(i))^2/(3*pi*ELAS*(Di(i))^4))*((Di(i)/Df(i))^3+(1/2)*(Di(i)/(Df(i)-Di(i)))^2*(1+(Di(i)/Df(i))^2-2*(Di(i)/Df(i))));
                   f35=f26;
                   f66=(64*L(i)/(3*pi*ELAS*(Di(i))^4))*(Di(i)/(Df(i)-Di(i)))*(1-(Di(i)/Df(i))^3);
                   f55=f66;
                   f44=(32*L(i)/(3*pi*GCOR*Di(i)^4))*(Di(i)/(Df(i)-Di(i)))*(1-(Di(i)/Df(i))^3);
                   FF44(i)=f44;
                   % ELEMENTOS DE LA MATRIZ DE RIGIDECES
                   raz=1/f11;
                   Detx=f22*f66-f26*f26;
                   r11x=f22/Detx;
                   r12x=(f26*L(i)-f22)/Detx;
                   r22x=(f66*(L(i))^2-2*f26*L(i)+f22)/Detx;
                   raax=(r11x+r22x+2*r12x)/(L(i))^2;
                   rabx=(r11x+r12x)/L(i);
                   rbax=(r22x+r12x)/L(i);
                   % MATRIZ DE RIGIDECES EN COORDENADAS LOCALES
                   K11=[raz 0 0 ; 0 raax rabx ; 0 rabx r11x];
                   K12=[-raz 0 0 ; 0 -raax rbax ; 0 -rabx r12x];
                   K21=K12';
                   K22=[raz 0 0 ; 0 raax -rbax ; 0 -rbax r22x];
                   K2=[K11 K12 ;K21 K22];
                   % MATRIZ DE TRANSPORTE A COORDENADAS GLOBALES
                   T2_3=zeros(6,6);
                   T2_3(1,1)=COSENO(i);T2_3(1,2)=SENO(i);
                   T2_3(2,1)=-SENO(i);T2_3(2,2)=COSENO(i);
                   T2_3(3,3)=1;
                   T2_3(4,4)=COSENO(i);T2_3(4,5)=SENO(i);
                   T2_3(5,4)=-SENO(i);T2_3(5,5)=COSENO(i);
                   T2_3(6,6)=1;
               end
           end   
       end
       
      %------------------ ENSAMBLAJE DE LA MATRIZ--------------------------
      %--------------------------------------------------------------------
      K3=T2_3'*K2*T2_3;
      MAT_AUX=K3;
      for j=1:6
         if VC(i,j)==0
            MAT_AUX(j,:)=0;
            MAT_AUX(:,j)=0;
         end
      end
      K=zeros(ngl,ngl);
      for cont1=1:6
         for cont2=1:6
            if MAT_AUX(cont1,cont2)~=0
               UNO=VC(i,cont1);
               DOS=VC(i,cont2);
               TRES=MAT_AUX(cont1,cont2);
               K(UNO,DOS)=K(UNO,DOS)+TRES;
            end
         end
      end
      if i==1
         AUX=K;
      elseif i~=1
         AUX=AUX+K;
      else,end
   end
   fprintf('\n MATRIZ DE RIGIDEZ DEL SISTEMA:');
   MAT_RIGIDEZ=AUX
   %-------------- ARREGLO Q_TOTAL ( VECTOR DE CARGAS ) -------------------
   Q=zeros(ngl,1);
   V1=menu('¿EXISTEN CARGAS EN JUNTAS?','SÍ','NO');
   if V1==1
       njc=input('\n INGRESE EL NÚMERO DE JUNTAS CARGADAS:');
       for i=1:njc
           NC=input('\n NUDO CARGADO:');
           Q1(1)=input(' FUERZA HORIZONTAL=');
           Q1(2)=input(' FUERZA VERTICAL=');
           Q1(3)=input(' MOMENTO=');
           VCJ(i,:)=CG(NC,:);
           for m=1:3
               n=VCJ(i,m)
               if Q1(m)~=0
                   Q(n)=Q1(m);
               end
           end
       end
       Q_CJ=Q
   else
       Q_CJ=Q
   end
   Q=zeros(ngl,1);
   Q2_ALMAC=zeros(mbr,6);
   V3=menu('¿EXISTEN CARGAS EN LOS MIEMBROS?','SÍ','NO');
   if V3==1
       nmc=input('\n INGRESE EL NÚMERO DE MIEMBROS CARGADOS:');
       for ll=1:nmc
           clear Q2
           fprintf('\n\t EL PROGRAMA UNICAMENTE CONSIDERA CAMBIOS POR CORTANTE EN LAS CARGAS');
           fprintf('\n\t PUNTUALES LO CUAL NO ES INCORRECTO DEBIDO A QUE EN ESTE ESTADO DE CARGA');
           fprintf('\n\t ES DONDE EXISTE MAYOR IMPORTANCIA A LAS DEFORMACIONES DE CORTANTE');
           fprintf('\n\t ');
           fprintf('\n\t ');
           fprintf('\n\t ');
           MC=input('\n INGRESE EN ORDEN NUMÉRICO CADA MIEMBRO CARGADO:');
           fprintf('\n\t ');
           fprintf('\n\t ESCOJA UNA DE LAS SIGUIENTES OPCIONES DE CARGAS:');
           fprintf('\n\t ');
           fprintf('\n\n (1) CARGA UNIFORMEMENTE DISTRIBUIDA EN TODO EL MIEMBRO');
           fprintf('\n (2) CARGA TRIANGULAR CON ALTURA MÁXIMA EN EL CENTRO');
           fprintf('\n (3) CARGA TRIANGULAR CON CARGA CRECIENTE DESDE EL NUDO INICIAL');
           fprintf('\n (4) CARGA TRAPEZOIDAL');
           fprintf('\n (5) CARGA PUNTUAL EN EL CENTRO DEL MIEMBRO');
           fprintf('\n (6) CARGA PUNTUAL A DIFERENTE DISTANCIA');
           fprintf('\n (7) NINGUNA DE LAS ANTERIORES');
           V4=input('\n\n\n OPCIÓN:');
           % CARGA UNIFORMEMENTE DISTRIBUIDA EN TODO EL MIEMBRO
           if V4==1
               CAR=input('\n CARGA (T/m):');
               Q2(1)=0;
               Q2(2)=CAR*L(MC)/2;
               Q2(3)=CAR*L(MC)^2/12;
               Q2(4)=0;
               Q2(5)=Q2(2);
               Q2(6)=-Q2(3);
               % CARGA TRIANGULARMENTE DISTRIBUIDA CON ALTURA MÁXIMA AL CENTRO
               % LONGITUDINAL DEL MIEMBRO
           elseif V4==2
               CAR=input('\n CARGA MÁXIMA (T/m):');
               Q2(1)=0;
               Q2(2)=CAR*L(MC)/4;
               Q2(3)=5*CAR*L(MC)^2/96;
               Q2(4)=0;
               Q2(5)=Q2(2);
               Q2(6)=-Q2(3); 
               % CARGA TRIANGULAR CRECIENTE DESDE EL NUDO INICIAL
           elseif V4==3   
               CAR=input('\n CARGA MÁXIMA (T/m):');
               Q2(1)=0;
               Q2(2)=3*CAR*L(MC)/20;
               Q2(3)=CAR*L(MC)^2/30;
               Q2(4)=0;
               Q2(5)=7*CAR*L(MC)/20;
               Q2(6)=-CAR*L(MC)^2/20;
               % CARGA TRAPEZOIDALMENTE DISTRIBUIDA
           elseif V4==4
               CAR=input('\n CARGA MÁXIMA (T/m):');
               a=input('\n DISTANCIA a (m):');
               Q2(1)=0;
               Q2(2)=(CAR*L(MC)/2)*(1-(a/L(MC)));
               Q2(3)=(CAR*L(MC)^2/12)*(1-2*(a/L(MC))^2+(a/L(MC))^3);
               Q2(4)=0;
               Q2(5)=Q2(2);
               Q2(6)=-Q2(3); 
               % CARGA PUNTUAL EN EL CENTRO DEL MIEMBRO
           elseif V4==5
               CAR=input('\n CARGA (T):');
               Q2(1)=0;
               Q2(2)=CAR/2;
               Q2(3)=CAR*L(MC)/8;
               Q2(4)=0;
               Q2(5)=Q2(2);
               Q2(6)=-Q2(3);  
               % CARGA PUNTUAL EN CUALQUIER PUNTO DEL EJE LONGITUDINAL DEL
               % MIEMBRO
           elseif V4==6
               CAR=input('\n CARGA (T):');
               a=input('\n DISTANCIA a (m):');
               b=input('\n DISTANCIA b (m):');
               Q2(1)=0;
               Q2(2)=(CAR*b^2/L(MC)^2)*(3-(2*b)/(a+b));
               Q2(3)=(CAR*a*b^2/L(MC)^2);
               Q2(4)=0;
               Q2(5)=(CAR*a^2/L(MC)^2)*(3-(2*a)/(a+b));
               Q2(6)=-(CAR*a^2*b/L(MC)^2);
               % CARGA POR INGRESAR; RECOMENDADO PARA MIEMBROS
               % INCLINADOS U OTRAS CONDICIONES DE CARGA
           elseif V4==7
               Q2(1)=input('\n FUERZA AXIAL DEL NUDO INICIAL (T):');
               Q2(2)=input('\n FUERZA CORTANTE DEL NUDO INICIAL (T):');
               Q2(3)=input('\n MOMENTO DEL NUDO INICIAL (T*m):');
               Q2(4)=input('\n FUERZA AXIAL DEL NUDO FINAL (T):');
               Q2(5)=input('\n FUERZA CORTANTE DEL NUDO FINAL (T):');
               Q2(6)=input('\n MOMENTO DEL NUDO FINAL (T*m):');
           end
           Q2=Q2'
           for mm=1:6
               Q2_ALMAC(MC,mm)=Q2(mm)';
           end
           % GENERACIÓN DE LAS CARGAS EN COORDENADAS GLOBALES
           T2_3=zeros(6,6);
           T2_3(1,1)=COSENO(MC);T2_3(1,2)=SENO(MC);
           T2_3(2,1)=-SENO(MC);T2_3(2,2)=COSENO(MC);
           T2_3(3,3)=1;
           T2_3(4,4)=COSENO(MC);T2_3(4,5)=SENO(MC);
           T2_3(5,4)=-SENO(MC);T2_3(5,5)=COSENO(MC);
           T2_3(6,6)=1;
           clear Q3
           Q3=-T2_3'*Q2;
           for g=1:6
               h=VC(MC,g);
               if h~=0
                   Q(h)=Q3(g)+Q(h);
               end
           end
       end
       Q_CM=Q
       Q3
   else
       Q_CM=Q
   end
   fprintf('\n VECTOR DE CARGAS:');
   Q_TOTAL=Q_CJ+Q_CM   
   %------------ Q_DES ----------------------------------------------------
   fprintf('\n DESPLAZAMIENTOS GENERALIZADOS:');
   Q_des=MAT_RIGIDEZ\Q_TOTAL
   %------------- ARREGLO P1 Y P_FINAL ------------------------------------
   for i=1:mbr
       fprintf('\n MIEMBRO %d:',i);
       fprintf('\n DEFORMACIONES EN MIEMBRO P1:');
       clear P1 P2 P Q2_AUX
       for j=1:6
           if VC(i,j)~=0
               P1(j)=Q_des(VC(i,j));
           else % ELIMINACIÓN DE LOS GRADOS SIN LIBERTAD EN LA ESTRUCTURA
               P1(j)=0;
           end
       end
       P1 % VECTOR DE DESPLAZAMIENTOS d
       if DD(i)==1
           %
           % MATRIZ DE RIGIDECES EN COORDENADAS LOCALES
           % PARA CADA CLASE DE ELEMENTO
           %
           % PARA LA GENERACIÓN DE VARIABLES SE 
           % INCLUYEN LOS SIGUIENTES CÁLCULOS
           C=((4*ELAS*INER(i))/L(i))*((1+FI(i))/(1+4*FI(i)));
           CP=C;
           A=((2*ELAS*INER(i))/L(i))*((1-2*FI(i))/(1+4*FI(i)));
           B=((6*ELAS*INER(i))/L(i)^2)*(1/(1+4*FI(i)));
           BP=B;
           T=(12*ELAS*INER(i))/L(i)^3*(1/(1+4*FI(i)));
           R=ELAS*AREA(i)/L(i);
           K2=[R 0 0 -R 0 0;
           0 T (B+C1(i)*T) 0 -T (B+C2(i)*T);
           0 (B+C1(i)*T) (C+2*C1(i)*B+C1(i)^2*T) 0 -(B+C1(i)*T) (A+C1(i)*B+C2(i)*B+C1(i)*C2(i)*T);
           -R 0 0 R 0 0;
           0 -T -(B+C1(i)*T) 0 T -(B+C2(i)*T);
           0 (B+C2(i)*T) (A+C1(i)*B+C2(i)*B+C1(i)*C2(i)*T) 0 -(B+C2(i)*T) (C+2*C2(i)*B+C2(i)^2*T)];
           % MATRIZ DE TRANSPORTE A COORDENADAS GLOBALES
           T2_3=zeros(6,6);
           T2_3(1,1)=COSENO(i);T2_3(1,2)=SENO(i);
           T2_3(2,1)=-SENO(i);T2_3(2,2)=COSENO(i);
           T2_3(3,3)=1;
           T2_3(4,4)=COSENO(i);T2_3(4,5)=SENO(i);
           T2_3(5,4)=-SENO(i);T2_3(5,5)=COSENO(i);
           T2_3(6,6)=1;
           % ELEMENTO CUADRADO SÓLIDO DE SECCIÓN VARIABLE
       elseif DD(i)==2
           % ELEMENTOS DE LA MATRIZ DE FLEXIBILIDADES
           f11=(L(i)/(ELAS*Hi(i)*Hi(i)))*(Hi(i)/(Hf(i)-Hi(i)))*(1-Hi(i)/Hf(i));
           f22=(4*(L(i))^3/(ELAS*(Hi(i))^4))*(Hi(i)/Hf(i))^3+(6*L(i)/(5*GCOR*(Hi(i))^2))*(Hi(i)/(Hf(i)-Hi(i)))*(1-Hi(i)/Hf(i));
           f33=f22;
           f26=(4*(L(i))^2/(ELAS*(Hi(i))^4))*((Hi(i)/Hf(i))^3+(1/2)*(Hi(i)/(Hf(i)-Hi(i)))^2*(1+(Hi(i)/Hf(i))^2-2*(Hi(i)/Hf(i))));
           f35=f26;
           f66=(4*L(i)/(ELAS*(Hi(i))^4))*(Hi(i)/(Hf(i)-Hi(i)))*(1-(Hi(i)/Hf(i))^3);
           f55=f66;
           f44=(64*L(i)/(27*GCOR*(Hi(i))^4))*(Hi(i)/(Hf(i)-Hi(i)))*(1-(Hi(i)/Hf(i))^3);
           % ELEMENTOS DE LA MATRIZ DE RIGIDECES
           raz=1/f11;
           Detx=f22*f66-f26*f26;
           r11x=f22/Detx;
           r12x=(f26*L(i)-f22)/Detx;
           r22x=(f66*(L(i))^2-2*f26*L(i)+f22)/Detx;
           raax=(r11x+r22x+2*r12x)/(L(i))^2;
           rabx=(r11x+r12x)/L(i);
           rbax=(r22x+r12x)/L(i);
           % MATRIZ DE RIGIDECES EN COORDENADAS LOCALES
           K11=[raz 0 0 ; 0 raax rabx ; 0 rabx r11x];
           K12=[-raz 0 0 ; 0 -raax rbax ; 0 -rabx r12x];
           K21=K12';
           K22=[raz 0 0 ; 0 raax -rbax ; 0 -rbax r22x];
           K2=[K11 K12 ;K21 K22];
           % MATRIZ DE TRANSPORTE A COORDENADAS GLOBALES
           T2_3=zeros(6,6);
           T2_3(1,1)=COSENO(i);T2_3(1,2)=SENO(i);
           T2_3(2,1)=-SENO(i);T2_3(2,2)=COSENO(i);
           T2_3(3,3)=1;
           T2_3(4,4)=COSENO(i);T2_3(4,5)=SENO(i);
           T2_3(5,4)=-SENO(i);T2_3(5,5)=COSENO(i);
           T2_3(6,6)=1;
           % ELEMENTO RECTANGULAR SÓLIDO DE SECCIÓN VARIABLE
       elseif DD(i)==3
           % ELEMENTOS DE LA MATRIZ DE FLEXIBILIDADES
           f11=(L(i)/(ELAS*Bb*Hi(i)))*(Hi(i)/(Hf(i)-Hi(i)))*log(Hf(i)/Hi(i));
           f22=(6*(L(i))^3/(ELAS*Bb*(Hi(i))^3))*((Hi(i)/Hf(i))^2-(Hi(i)/(Hf(i)-Hi(i)))^3*(Hf(i)/Hi(i)-Hi(i)/Hf(i)-2*log(Hf(i)/Hi(i))))+(6*L(i)/(5*GCOR*Bb*Hi(i)))*(Hi(i)/(Hf(i)-Hi(i)))*log(Hf(i)/Hi(i));
           f33=(12*(L(i))^3/(ELAS*(Bb)^3*Hi(i)))*(3/2+(1/2)*(Hf(i)/Hi(i))^2-2*(Hf(i)/Hi(i))+log(Hf(i)/Hi(i)))+(6*L(i)/(5*GCOR*Bb*Hi(i)))*(Hi(i)/(Hf(i)-Hi(i)))*log(Hf(i)/Hi(i));
           f26=(6*L(i)/(ELAS*Bb*(Hi(i))^3))*(Hi(i)/(Hf(i)-Hi(i)))^2*(1+(Hi(i)/Hf(i))^2-2*(Hi(i)/Hf(i)));
           f35=(12*(L(i))^2/(ELAS*(Bb)^3*Hi(i)))*(Hi(i)/(Hf(i)-Hi(i)))^2*(Hf(i)/Hi(i)-log(Hf(i)/Hi(i))-1);
           f66=(6*L(i)/(ELAS*Bb*(Hi(i))^3))*(Hi(i)/(Hf(i)-Hi(i)))*(1-(Hi(i)/Hf(i))^2);
           f55=(12*L(i)/(ELAS*(Bb)^3*Hi(i)))*(Hi(i)/(Hf(i)-Hi(i)))*log(Hf(i)/Hi(i));
           f44=FF44(i);
           % ELEMENTOS DE LA MATRIZ DE RIGIDECES
           raz=1/f11;
           Detx=f22*f66-f26*f26;
           r11x=f22/Detx;
           r12x=(f26*L(i)-f22)/Detx;
           r22x=(f66*(L(i))^2-2*f26*L(i)+f22)/Detx;
           raax=(r11x+r22x+2*r12x)/(L(i))^2;
           rabx=(r11x+r12x)/L(i);
           rbax=(r22x+r12x)/L(i);
           % MATRIZ DE RIGIDECES EN COORDENADAS LOCALES
           K11=[raz 0 0 ; 0 raax rabx ; 0 rabx r11x];
           K12=[-raz 0 0 ; 0 -raax rbax ; 0 -rabx r12x];
           K21=K12';
           K22=[raz 0 0 ; 0 raax -rbax ; 0 -rbax r22x];
           K2=[K11 K12 ;K21 K22];
           % MATRIZ DE TRANSPORTE A COORDENADAS GLOBALES
           T2_3=zeros(6,6);
           T2_3(1,1)=COSENO(i);T2_3(1,2)=SENO(i);
           T2_3(2,1)=-SENO(i);T2_3(2,2)=COSENO(i);
           T2_3(3,3)=1;
           T2_3(4,4)=COSENO(i);T2_3(4,5)=SENO(i);
           T2_3(5,4)=-SENO(i);T2_3(5,5)=COSENO(i);
           T2_3(6,6)=1;
           % ELEMENTO CIRCULAR SÓLIDO DE SECCIÓN CONSTATE
       elseif DD(i)==4
           C=((4*ELAS*INER(i))/L(i))*((1+FI(i))/(1+4*FI(i)));
           CP=C;
           A=((2*ELAS*INER(i))/L(i))*((1-2*FI(i))/(1+4*FI(i)));
           B=((6*ELAS*INER(i))/L(i)^2)*(1/(1+4*FI(i)));
           BP=B;
           T=(12*ELAS*INER(i))/L(i)^3*(1/(1+4*FI(i)));
           R=ELAS*AREA(i)/L(i);
           K2=[R 0 0 -R 0 0;
           0 T (B+C1(i)*T) 0 -T (B+C2(i)*T);
           0 (B+C1(i)*T) (C+2*C1(i)*B+C1(i)^2*T) 0 -(B+C1(i)*T) (A+C1(i)*B+C2(i)*B+C1(i)*C2(i)*T);
           -R 0 0 R 0 0;
           0 -T -(B+C1(i)*T) 0 T -(B+C2(i)*T);
           0 (B+C2(i)*T) (A+C1(i)*B+C2(i)*B+C1(i)*C2(i)*T) 0 -(B+C2(i)*T) (C+2*C2(i)*B+C2(i)^2*T)];
           % MATRIZ DE TRANSPORTE A COORDENADAS GLOBALES
           T2_3=zeros(6,6);
           T2_3(1,1)=COSENO(i);T2_3(1,2)=SENO(i);
           T2_3(2,1)=-SENO(i);T2_3(2,2)=COSENO(i);
           T2_3(3,3)=1;
           T2_3(4,4)=COSENO(i);T2_3(4,5)=SENO(i);
           T2_3(5,4)=-SENO(i);T2_3(5,5)=COSENO(i);
           T2_3(6,6)=1;
           % ELEMENTO CIRCULAR SÓLIDO DE SECCIÓN VARIABLE
       elseif DD(i)==5
           % ELEMENTOS DE LA MATRIZ DE FLEXIBILIDADES
           f11=(4*L(i)/(pi*ELAS*(Di(i))^2))*(Di(i)/(Df(i)-Di(i)))*(1-Di(i)/Df(i));
           f22=(64*(L(i))^3/(3*pi*ELAS*(Di(i))^4))*(Di(i)/Df(i))^3+(40*L(i)/(9*pi*GCOR*(Di(i))^2))*(Di(i)/(Df(i)-Di(i)))*(1-(Di(i)/Df(i)));
           f33=f22;
           f26=(64*(L(i))^2/(3*pi*ELAS*(Di(i))^4))*((Di(i)/Df(i))^3+(1/2)*(Di(i)/(Df(i)-Di(i)))^2*(1+(Di(i)/Df(i))^2-2*(Di(i)/Df(i))));
           f35=f26;
           f66=(64*L(i)/(3*pi*ELAS*(Di(i))^4))*(Di(i)/(Df(i)-Di(i)))*(1-(Di(i)/Df(i))^3);
           f55=f66;
           f44=(32*L(i)/(3*pi*GCOR*Di(i)^4))*(Di(i)/(Df(i)-Di(i)))*(1-(Di(i)/Df(i))^3);
           % ELEMENTOS DE LA MATRIZ DE RIGIDECES
           raz=1/f11;
           Detx=f22*f66-f26*f26;
           r11x=f22/Detx;
           r12x=(f26*L(i)-f22)/Detx;
           r22x=(f66*(L(i))^2-2*f26*L(i)+f22)/Detx;
           raax=(r11x+r22x+2*r12x)/(L(i))^2;
           rabx=(r11x+r12x)/L(i);
           rbax=(r22x+r12x)/L(i);
           % MATRIZ DE RIGIDECES EN COORDENADAS LOCALES
           K11=[raz 0 0 ; 0 raax rabx ; 0 rabx r11x];
           K12=[-raz 0 0 ; 0 -raax rbax ; 0 -rabx r12x];
           K21=K12';
           K22=[raz 0 0 ; 0 raax -rbax ; 0 -rbax r22x];
           K2=[K11 K12 ;K21 K22];
           % MATRIZ DE TRANSPORTE A COORDENADAS GLOBALES
           T2_3=zeros(6,6);
           T2_3(1,1)=COSENO(i);T2_3(1,2)=SENO(i);
           T2_3(2,1)=-SENO(i);T2_3(2,2)=COSENO(i);
           T2_3(3,3)=1;
           T2_3(4,4)=COSENO(i);T2_3(4,5)=SENO(i);
           T2_3(5,4)=-SENO(i);T2_3(5,5)=COSENO(i);
           T2_3(6,6)=1;
       end
       % LO ANTERIOR SE REALIZA PARA LA GENERACIÓN NUEVAMENTE DE LA MATRIZ
       % DE RIGIDEZ
       fprintf('\n RESULTADOS FINALES DEL MIEMBRO %d:',i);
       fprintf('\n MATRIZ DE RIGIDEZ EN COORDENADAS LOCALES:\n');
       K2 % K EN COORDENADAS LOCALES
       fprintf('\n MATRIZ DE RIGIDEZ EN COORDENADAS GLOBALES:\n');
       K3=T2_3'*K2*T2_3 % K EN COORDENADAS GLOBALES
       fprintf('\n VECTOR DE CARGAS EN COORDENADAS GLOBALES:\n');
       P2=K3*P1';% K*d=q CARGAS EN COORDENADAS GLOBALES
       fprintf('\n ACCIONES EN COORDENADAS LOCALES (PROBLEMA DE LA ESTRUCTURA LIBERADA):\n');
       P=T2_3*P2% R*q=q' (P.COMPLEMENTARIO) % CARGA EN COORDENADAS LOCALES 
       fprintf('\n ACCIONES DE EMPOTRAMIENTO EN COORDENADAS LOCALES (PROBLEMA DE LA ESTRUCTURA RESTRINGIDA):\n');
       for j=1:6
           Q2_AUX(j)=Q2_ALMAC(i,j); % VALORES CORRESPONDIENTES A q DE EMPOTRAMIENTO 
         %                          EN COORDENADAS GLOBALES
       end
       Q2_AUX=Q2_AUX' % q_emp
       fprintf('\n ACCIONES FINALES DEL MIEMBRO EN COORDENADAS LOCALES ( P.P.+P.C. ):');
       P_FINAL=Q2_AUX+P % ELEMENTOS MECÁNICOS FINALES: q_emp+q'
   end
elseif DV==2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% DIVISIÓN DE LOS PROGRAMAS APLICADOS %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    clc
    clear all
    fprintf('\n\t\t\t ‘;|                     ¡INDESTRUCT II!                       |');
    fprintf('\n\n\t\t');
    fprintf('\n\n\t\t PROGRAMA PARA CALCULAR MARCOS PLANOS CON NUDOS RÍGIDOS Y CORTANTE');
    fprintf('\n\n\t\t CONSIDERANDO LA TEORÍA DE SEGUNDO ORDEN');
    % -------- ARREGLO CG --------
    % EN ESTA ETAPA SE DESARROLLA LA MATRIZ DE COORDENADAS GENERALIZADAS DONDE
    % SE ESPECIFICAN LOS GRADOS DE LIBERTAD EXISTENTES DEPENDIENDO DEL TIPO DE
    % APOYO
    fprintf('\n\t\t\t COORDENADAS     Y|                                 ');
    fprintf('\n\t\t\t        GLOBALES  |                                 ');
    fprintf('\n\t\t\t                  |                                 ');
    fprintf('\n\t\t\t                  |_ _ _ _ _                        ');
    fprintf('\n\t\t\t                              X                     ');
    fprintf('\n\t\t\t                                                    ');
    %%
    nod = input('\n\n\n INGRESE EL NÚMERO DE NUDOS:');
    nnr = input(' NÚMERO DE NUDOS RESTRINGIDOS:');
    CG = ones(nod,3 );
    for i=1: nnr
        nudo_res = input('\n NÚMERO DEL NUDO RESTRINGIDO:');
        x1 = input('\n EXISTEN DESPLAZAMIENTO EN X,S o N:','s');
        if x1=='N'
            CG(nudo_res,1)=0;
        else,end
        y1 = input(' EXISTEN DESPLAZAMIENTO EN Y,S o N:','s');
        if y1=='N'
            CG(nudo_res,2)=0;
        else,end
        r1=input(' EXISTE ROTACIÓN, S o N:','s');
        if r1=='N'
            CG(nudo_res,3)=0;
        else,end
    end
    ngl=0;
    for i=1:nod
        for j=1:3
            if CG(i,j)~=0
                ngl=ngl+1;
                CG(i,j)=ngl;
            else,end
        end
    end
    CG
    ngl
    %--------------ARREGLO VC ----------------------------------------------
    mbr = input('\n INGRESE EL NUMERO DE MIEMBROS:');
    for i=1:mbr
        fprintf('\n MIEMBRO %d:',i);
        NUDO_INICIAL(i)=input('\n INGRESE EL NUDO INICIAL:');
        NUDO_FINAL(i)=input(' INGRESE EL NUDO FINAL:');
    end
    NUDO_INICIAL
    NUDO_FINAL
    for i=1:mbr
        for j=1:3
            VC(i,j)=CG(NUDO_INICIAL(i),j);
            VC(i,j+3)=CG(NUDO_FINAL(i),j);
        end
    end
    VC
    %--------------- ARREGLO L, SENO, COSENO -------------------------------
    fprintf('\n COORDENADAS EN LOS NUDOS:\n');
    for i=1:nod
        fprintf('\n NUDO %d:',i);
        x(i)=input('\n COORDENADA [m]:x:');
        y(i)=input(' COORDENADA [m]:y:');
    end
    for i=1:mbr
        Dx(i)=x(NUDO_FINAL(i))-x(NUDO_INICIAL(i));
        Dy(i)=y(NUDO_FINAL(i))-y(NUDO_INICIAL(i));
        L(i)=(Dx(i)^2+Dy(i)^2)^0.5;
        SENO(i)=Dy(i)/L(i);
        COSENO(i)=Dx(i)/L(i);
    end
    L
    SENO
    COSENO
    ZZZ=menu('¿LOS DATOS INTRODUCIDOS SON CORRECTOS?','SÍ','NO');
    fprintf('\n\n\n'); 
    if ZZZ==2     
        fprintf('\n\n\t\t FIN DEL PROGRAMA, PRESIONE Y VOLVER A INICIAR');
        pause
        open('edwINDESTRUCT.m')
        break
    else
    end   
   %----------------- K3-------------------------------------------------
   fprintf('\n CARACTERÍSTICAS DE LOS MIEMBROS:');
   ELAS=input('\n MÓDULO DE ELASTICIDAD [T/m^2]:');
   POISSON=input('\n MÓDULO DE POISSON (0.2):');
   GCOR=ELAS/(2*(1+POISSON));
   M_CORTANTE=GCOR
   T2_3=zeros(6,6);
   fprintf('\n\n UNICAMENTE EN LA SECCIÓN RECTANGULAR PRISMÁTICA SE PUEDEN MODELAR EXTREMOS INFINITAMENTE RÍGIDOS');
   for i=1:mbr
       if i==1
           fprintf('\n MIEMBRO %d:',i);
           fprintf('\n\n ESCOJA UNA DE LAS SIGUIENTES OPCIONES:');
           fprintf('\n\n\t\t (1) MIEMBRO PRISMÁTICO DE SECCIÓN RECTANGULAR');
           fprintf('\n\n\t\t (2) MIEMBRO CUADRADO DE SECCIÓN VARIABLE');
           fprintf('\n\n\t\t (3) MIEMBRO RECTANGULAR DE SECCIÓN VARIABLE');
           fprintf('\n\n\t\t (4) MIEMBRO CIRCULAR DE SECCIÓN PRISMÁTICA CONSTANTE');
           fprintf('\n\n\t\t (5) MIEMBRO CIRCULAR DE SECCIÓN VARIABLE');
           DD(i) = input('\n\n OPCION : ');
           fprintf('\n\n\n'); 
           % ELEMENTO PRISMÁTICO DE SECCIÓN CONSTANTE DE FORMA RECTANGULAR
          if DD(i)==1
              B(i)=input('\n BASE DEL ELEMENTO [m]:');
              H(i)=input(' ALTURA DEL ELEMENTO [m]:');
              % ELEMENTOS UNICAMENTE NECESARIOS PARA EL CONTEO DURANTE EL USO
              % DE ELEMENTOS DE CUALQUIERA DE LOS DISTINTOS TIPOS POSIBLES DE
              % SECCIONES
              Bi(i)=B(i);
              Hi(i)=H(i);
              Bf(i)=Bi(i);
              Hf(i)=Hi(i);
              Di(i)=B(i);
              Df(i)=Di(i);
              FF44(i)=0;
              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              BET(i)=input('\n VALOR DEL FACTOR DE FORMA (RECTANGULAR=1.2; CIRCULAR=10/9):');
              C1(i)=input('\n\n INGRESE EL VALOR DE C1 [m] (LONGITUD DEL SEGMENTO DE RIGIDEZ INFINITA EN EL NODO INICIAL):');
              C2(i)=input('\n\n INGRESE EL VALOR DE C2 [m] (LONGITUD DEL SEGMENTO DE RIGIDEZ INFINITA EN EL NODO FINAL):');
              F(i)=L(i)-C1(i)-C2(i)
              % CONSTANTES AUXILIARES PARA FUTUROS DESARROLLOS EN 3D
              AREA(i)=B(i)*H(i);
              INERx(i)=(B(i)*H(i)^3)/12;
              INER(i)=INERx(i);
              INERy(i)=(H(i)*B(i)^3)/12;
              FIx(i)=(12*ELAS*INERx(i)*BET(i))/(GCOR*AREA(i)*L(i)^2);
              FIy(i)=(12*ELAS*INERy(i)*BET(i))/(GCOR*AREA(i)*L(i)^2);
              % GENERACIÓN DE LAS VARIABLES
              FI(i)=(1/4)*FIx(i);
              C=((4*ELAS*INER(i))/L(i))*((1+FI(i))/(1+4*FI(i)));
              CP=C;
              A=((2*ELAS*INER(i))/L(i))*((1-2*FI(i))/(1+4*FI(i)));
              B=((6*ELAS*INER(i))/L(i)^2)*(1/(1+4*FI(i)));
              BP=B;
              T=(12*ELAS*INER(i))/L(i)^3*(1/(1+4*FI(i)));
              R=ELAS*AREA(i)/L(i);
              % MATRIZ DE RIGIDECES EN COORDENADAS LOCALES
              K2=[R 0 0 -R 0 0;
              0 T (B+C1(i)*T) 0 -T (B+C2(i)*T);
              0 (B+C1(i)*T) (C+2*C1(i)*B+C1(i)^2*T) 0 -(B+C1(i)*T) (A+C1(i)*B+C2(i)*B+C1(i)*C2(i)*T);
              -R 0 0 R 0 0;
              0 -T -(B+C1(i)*T) 0 T -(B+C2(i)*T);
              0 (B+C2(i)*T) (A+C1(i)*B+C2(i)*B+C1(i)*C2(i)*T) 0 -(B+C2(i)*T) (C+2*C2(i)*B+C2(i)^2*T)];
              % MATRIZ DE TRANSPORTE A COORDENADAS GLOBALES
              T2_3=zeros(6,6);
              T2_3(1,1)=COSENO(i);T2_3(1,2)=SENO(i);
              T2_3(2,1)=-SENO(i);T2_3(2,2)=COSENO(i);
              T2_3(3,3)=1;
              T2_3(4,4)=COSENO(i);T2_3(4,5)=SENO(i);
              T2_3(5,4)=-SENO(i);T2_3(5,5)=COSENO(i);
              T2_3(6,6)=1;
              % ELEMENTO CUADRADO SÓLIDO DE SECCIÓN VARIABLE
          elseif DD(i)==2
              Bi(i)=input('\n BASE INICIAL [m]:');
              Hi(i)=Bi(i);
              Bf(i)=input('\n BASE FINAL [m]:');
             % ELEMENTOS UNICAMENTE NECESARIOS PARA EL CONTEO DURANTE EL USO
             % DE ELEMENTOS DE CUALQUIERA DE LOS DISTINTOS TIPOS POSIBLES DE
             % SECCIONES
             Hf(i)=Bf(i);
             B(i)=Bi(i);
             H(i)=Hi(i);
             Di(i)=Bi(i);
             Df(i)=Di(i);
             BET(i)=1.2;
             INER(i)=(B(i)*H(i)^3)/12;
             FI(i)=1;
             AREA(i)=B(i)*H(i);
             BET(i)=1.2;
             C1(i)=0;
             C2(i)=0;
             F(i)=L(i);
             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             % ELEMENTOS DE LA MATRIZ DE FLEXIBILIDADES
             f11=(L(i)/(ELAS*Hi(i)*Hi(i)))*(Hi(i)/(Hf(i)-Hi(i)))*(1-Hi(i)/Hf(i));
             f22=(4*(L(i))^3/(ELAS*(Hi(i))^4))*(Hi(i)/Hf(i))^3+(6*L(i)/(5*GCOR*(Hi(i))^2))*(Hi(i)/(Hf(i)-Hi(i)))*(1-Hi(i)/Hf(i));
             f33=f22;
             f26=(4*(L(i))^2/(ELAS*(Hi(i))^4))*((Hi(i)/Hf(i))^3+(1/2)*(Hi(i)/(Hf(i)-Hi(i)))^2*(1+(Hi(i)/Hf(i))^2-2*(Hi(i)/Hf(i))));
             f35=f26;
             f66=(4*L(i)/(ELAS*(Hi(i))^4))*(Hi(i)/(Hf(i)-Hi(i)))*(1-(Hi(i)/Hf(i))^3);
             f55=f66;
             f44=(64*L(i)/(27*GCOR*(Hi(i))^4))*(Hi(i)/(Hf(i)-Hi(i)))*(1-(Hi(i)/Hf(i))^3);
             FF44(i)=f44;
             % ELEMENTOS DE LA MATRIZ DE RIGIDECES
             raz=1/f11;
             Detx=f22*f66-f26*f26;
             r11x=f22/Detx;
             r12x=(f26*L(i)-f22)/Detx;
             r22x=(f66*(L(i))^2-2*f26*L(i)+f22)/Detx;
             raax=(r11x+r22x+2*r12x)/(L(i))^2;
             rabx=(r11x+r12x)/L(i);
             rbax=(r22x+r12x)/L(i);
             % MATRIZ DE RIGIDECES EN COORDENADAS LOCALES
             K11=[raz 0 0 ; 0 raax rabx ; 0 rabx r11x];
             K12=[-raz 0 0 ; 0 -raax rbax ; 0 -rabx r12x];
             K21=K12';
             K22=[raz 0 0 ; 0 raax -rbax ; 0 -rbax r22x];
             K2=[K11 K12 ;K21 K22];
             % MATRIZ DE TRANSPORTE A COORDENADAS GLOBALES
             T2_3=zeros(6,6);
             T2_3(1,1)=COSENO(i);T2_3(1,2)=SENO(i);
             T2_3(2,1)=-SENO(i);T2_3(2,2)=COSENO(i);
             T2_3(3,3)=1;
             T2_3(4,4)=COSENO(i);T2_3(4,5)=SENO(i);
             T2_3(5,4)=-SENO(i);T2_3(5,5)=COSENO(i);
             T2_3(6,6)=1;
             % ELEMENTO RECTANGULAR SÓLIDO DE SECCIÓN VARIABLE
          elseif DD(i)==3
              Bi(i)=input('\n BASE DE TODO EL ELEMENTO [m]:');
              Hi(i)=input(' ALTURA INICIAL DEL ELEMENTO [m]:');
              Bf(i)=Bi(i);
              Hf(i)=input(' ALTURA FINAL DEL ELEMENTO [m]:');
              Bb=Bi(i);
              % ELEMENTOS UNICAMENTE NECESARIOS PARA EL CONTEO DURANTE EL USO
              % DE ELEMENTOS DE CUALQUIERA DE LOS DISTINTOS TIPOS POSIBLES DE
              % SECCIONES
              B(i)=Bi(i);
              H(i)=Hi(i);
              Di(i)=Bi(i);
              Df(i)=Di(i);
              INER(i)=(B(i)*H(i)^3)/12;
              FI(i)=1;
              AREA(i)=B(i)*H(i);
              BET(i)=1.2;
              C1(i)=0;
              C2(i)=0;
              F(i)=L(i);
              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              BET(i)=input('\n FACTOR DE FORMA (RECTANGULAR=1.2; CIRCULAR=10/9):');
              fprintf('\n\t\t\t ‘ ; ......................................................................................................................... ');
              fprintf('\n\t\t\t ‘ ; ......................................................................................................................... ');
              fprintf('\n\n\t\t     LA FUNCIÓN F(x) INDICA LA VARIACIÓN CON LA LONGITUD DE LA FORMA DEL ELEMENTO,');
              fprintf('\n\t\t\t                                                                    ');
              fprintf('\n\n\t\t     SI LA FUNCIÓN ES LINEAL: Y=F(x)=');
              fprintf('\n\t\t\t                                      (abs(Hi-Hf)/LL)*x                                       ');
              fprintf('\n\t\t\t ‘ ; ......................................................................................................................... ');
              fprintf('\n\t\t\t ‘ ; ......................................................................................................................... ');
              fprintf('\n\t\t\t ‘ ;                                                                                                                           ');
              fprintf('\n\t\t\t   ESTANDO EL EJE X SOBRE EL     Y|                                                                                            ');
              fprintf('\n\t\t\t    CENTROIDE DEL MIEMBRO         |                                                                                            ');
              fprintf('\n\t\t\t                                  |                                                                                            ');
              fprintf('\n\t\t\t                                  |_ _ _ _ _                                                                                   ');
              fprintf('\n\t\t\t                                           X                                                                                   ');
              fprintf('\n\t\t\t                                                                                                                               ');
              fprintf('\n\t\t\t ‘ ; ......................................................................................................................... ');
              fprintf('\n\t\t\t ‘ ; ......................................................................................................................... ');
              fprintf('\n\t\t\t                                                                                                                              ');
              fprintf('\n\t\t\t                                                                                                                              ');
              EE(i)=menu('ESCOJA UNA DE LAS SIGUIENTES OPCIONES:','FUNCIÓN LINEAL','OTRO TIPO DE VARIACIÓN');
              if EE(i)==1
                  % ELEMENTOS DE LA MATRIZ DE FLEXIBILIDADES
                  f11=(L(i)/(ELAS*Bb*Hi(i)))*(Hi(i)/(Hf(i)-Hi(i)))*log(Hf(i)/Hi(i));
                  f22=(6*(L(i))^3/(ELAS*Bb*(Hi(i))^3))*((Hi(i)/Hf(i))^2-(Hi(i)/(Hf(i)-Hi(i)))^3*(Hf(i)/Hi(i)-Hi(i)/Hf(i)-2*log(Hf(i)/Hi(i))))+(6*L(i)/(5*GCOR*Bb*Hi(i)))*(Hi(i)/(Hf(i)-Hi(i)))*log(Hf(i)/Hi(i));
                  f33=(12*(L(i))^3/(ELAS*(Bb)^3*Hi(i)))*(3/2+(1/2)*(Hf(i)/Hi(i))^2-2*(Hf(i)/Hi(i))+log(Hf(i)/Hi(i)))+(6*L(i)/(5*GCOR*Bb*Hi(i)))*(Hi(i)/(Hf(i)-Hi(i)))*log(Hf(i)/Hi(i));
                  f26=(6*L(i)/(ELAS*Bb*(Hi(i))^3))*(Hi(i)/(Hf(i)-Hi(i)))^2*(1+(Hi(i)/Hf(i))^2-2*(Hi(i)/Hf(i)));
                  f35=(12*(L(i))^2/(ELAS*(Bb)^3*Hi(i)))*(Hi(i)/(Hf(i)-Hi(i)))^2*(Hf(i)/Hi(i)-log(Hf(i)/Hi(i))-1);
                  f66=(6*L(i)/(ELAS*Bb*(Hi(i))^3))*(Hi(i)/(Hf(i)-Hi(i)))*(1-(Hi(i)/Hf(i))^2);
                  f55=(12*L(i)/(ELAS*(Bb)^3*Hi(i)))*(Hi(i)/(Hf(i)-Hi(i)))*log(Hf(i)/Hi(i));
                  % DEFINICIÓN DE UNA FUNCIÓN DE x
                  Hi=Hi(i);
                  Hf=Hf(i);
                  LL=L(i);
                  if Hi<Hf
                      dpc= @(x)(1./(((Bb.^3*(-Hi./2-((Hi-Hf)./(2*LL))*x))./3)*(1-(192*Bb./(pi.^5*(-Hi./2-((Hi-Hf)./(2*LL))*x)))*(tanh(pi*(-Hi./2-((Hi-Hf)./(2*LL))*x)./(2*Bb)))))); 
                  elseif Hi>Hf
                      dpc= @(x)(1./(((Bb.^3*(-Hi./2+((Hi-Hf)./(2*LL))*x))./3)*(1-(192*Bb./(pi.^5*(-Hi./2+((Hi-Hf)./(2*LL))*x)))*(tanh(pi*(-Hi./2+((Hi-Hf)./(2*LL))*x)./(2*Bb))))));
                  end
                  f44=quad(@(x)dpc, 0, LL)
                  FF44(i)=f44;
                  % ELEMENTOS DE LA MATRIZ DE RIGIDECES
                  raz=1/f11;
                  Detx=f22*f66-f26*f26;
                  r11x=f22/Detx;
                  r12x=(f26*L(i)-f22)/Detx;
                  r22x=(f66*(L(i))^2-2*f26*L(i)+f22)/Detx;
                  raax=(r11x+r22x+2*r12x)/(L(i))^2;
                  rabx=(r11x+r12x)/L(i);
                  rbax=(r22x+r12x)/L(i);
                  % MATRIZ DE RIGIDECES EN COORDENADAS LOCALES
                  K11=[raz 0 0 ; 0 raax rabx ; 0 rabx r11x];
                  K12=[-raz 0 0 ; 0 -raax rbax ; 0 -rabx r12x];
                  K21=K12';
                  K22=[raz 0 0 ; 0 raax -rbax ; 0 -rbax r22x];
                  K2=[K11 K12 ;K21 K22];
                  % MATRIZ DE TRANSPORTE A COORDENADAS GLOBALES
                  T2_3=zeros(6,6);
                  T2_3(1,1)=COSENO(i);T2_3(1,2)=SENO(i);
                  T2_3(2,1)=-SENO(i);T2_3(2,2)=COSENO(i);
                  T2_3(3,3)=1;
                  T2_3(4,4)=COSENO(i);T2_3(4,5)=SENO(i);
                  T2_3(5,4)=-SENO(i);T2_3(5,5)=COSENO(i);
                  T2_3(6,6)=1;
              else
                  % ELEMENTOS DE LA MATRIZ DE FLEXIBILIDADES
                  f11=(L(i)/(ELAS*Bb*Hi(i)))*(Hi(i)/(Hf(i)-Hi(i)))*log(Hf(i)/Hi(i));
                  f22=(6*(L(i))^3/(ELAS*Bb*(Hi(i))^3))*((Hi(i)/Hf(i))^2-(Hi(i)/(Hf(i)-Hi(i)))^3*(Hf(i)/Hi(i)-Hi(i)/Hf(i)-2*log(Hf(i)/Hi(i))))+(6*L(i)/(5*GCOR*Bb*Hi(i)))*(Hi(i)/(Hf(i)-Hi(i)))*log(Hf(i)/Hi(i));
                  f33=(12*(L(i))^3/(ELAS*(Bb)^3*Hi(i)))*(3/2+(1/2)*(Hf(i)/Hi(i))^2-2*(Hf(i)/Hi(i))+log(Hf(i)/Hi(i)))+(6*L(i)/(5*GCOR*Bb*Hi(i)))*(Hi(i)/(Hf(i)-Hi(i)))*log(Hf(i)/Hi(i));
                  f26=(6*L(i)/(ELAS*Bb*(Hi(i))^3))*(Hi(i)/(Hf(i)-Hi(i)))^2*(1+(Hi(i)/Hf(i))^2-2*(Hi(i)/Hf(i)));
                  f35=(12*(L(i))^2/(ELAS*(Bb)^3*Hi(i)))*(Hi(i)/(Hf(i)-Hi(i)))^2*(Hf(i)/Hi(i)-log(Hf(i)/Hi(i))-1);
                  f66=(6*L(i)/(ELAS*Bb*(Hi(i))^3))*(Hi(i)/(Hf(i)-Hi(i)))*(1-(Hi(i)/Hf(i))^2);
                  f55=(12*L(i)/(ELAS*(Bb)^3*Hi(i)))*(Hi(i)/(Hf(i)-Hi(i)))*log(Hf(i)/Hi(i));
                  % DEFINICIÓN DE UNA FUNCIÓN DE x
                  Hi=Hi(i);
                  Hf=Hf(i);
                  LL=L(i);
                  f44=input(' INTRODUZCA LA f44=');
                  FF44(i)=f44;
                  % ELEMENTOS DE LA MATRIZ DE RIGIDECES
                  raz=1/f11;
                  Detx=f22*f66-f26*f26;
                  r11x=f22/Detx;
                  r12x=(f26*L(i)-f22)/Detx;
                  r22x=(f66*(L(i))^2-2*f26*L(i)+f22)/Detx;
                  raax=(r11x+r22x+2*r12x)/(L(i))^2;
                  rabx=(r11x+r12x)/L(i);
                  rbax=(r22x+r12x)/L(i);
                  % MATRIZ DE RIGIDECES EN COORDENADAS LOCALES
                  K11=[raz 0 0 ; 0 raax rabx ; 0 rabx r11x];
                  K12=[-raz 0 0 ; 0 -raax rbax ; 0 -rabx r12x];
                  K21=K12';
                  K22=[raz 0 0 ; 0 raax -rbax ; 0 -rbax r22x];
                  K2=[K11 K12 ;K21 K22];
                  % MATRIZ DE TRANSPORTE A COORDENADAS GLOBALES
                  T2_3=zeros(6,6);
                  T2_3(1,1)=COSENO(i);T2_3(1,2)=SENO(i);
                  T2_3(2,1)=-SENO(i);T2_3(2,2)=COSENO(i);
                  T2_3(3,3)=1;
                  T2_3(4,4)=COSENO(i);T2_3(4,5)=SENO(i);
                  T2_3(5,4)=-SENO(i);T2_3(5,5)=COSENO(i);
                  T2_3(6,6)=1;
              end
              % ELEMENTO CIRCULAR SÓLIDO DE SECCIÓN CONSTATE
          elseif DD(i)==4
              Di(i)=input('\n DIÁMETRO DEL ELEMENTO [m]:');
              % ELEMENTOS UNICAMENTE NECESARIOS PARA EL CONTEO DURANTE EL USO
              % DE ELEMENTOS DE CUALQUIERA DE LOS DISTINTOS TIPOS POSIBLES DE
              % SECCIONES
              Df(i)=Di(i);
              B(i)=Di(i);
              H(i)=Di(i);
              Bi(i)=B(i);
              Hi(i)=H(i);
              Bf(i)=Bi(i);
              Hf(i)=Hi(i);
              FF44(i)=0;
              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              BET(i)=input('\n VALOR DEL FACTOR DE FORMA (RECTANGULAR=1.2; CIRCULAR=10/9):');
              AREA(i)=pi*Di(i)*Di(i)/4;
              INER(i)=pi*(Di(i)/2)^4/4;
              C1(i)=input('\n\n INGRESE EL VALOR DE C1 [m] (LONGITUD DEL SEGMENTO DE RIGIDEZ INFINITA EN EL NODO INICIAL):');
              C2(i)=input('\n\n INGRESE EL VALOR DE C2 [m] (LONGITUD DEL SEGMENTO DE RIGIDEZ INFINITA EN EL NODO FINAL):');
              F(i)=L(i)-C1(i)-C2(i)
              FI(i)=(3*ELAS*INER(i)*BET(i))/(GCOR*AREA(i)*L(i)^2);
              C=((4*ELAS*INER(i))/L(i))*((1+FI(i))/(1+4*FI(i)));
              CP=C;
              A=((2*ELAS*INER(i))/L(i))*((1-2*FI(i))/(1+4*FI(i)));
              B=((6*ELAS*INER(i))/L(i)^2)*(1/(1+4*FI(i)));
              BP=B;
              T=(12*ELAS*INER(i))/L(i)^3*(1/(1+4*FI(i)));
              R=ELAS*AREA(i)/L(i);
              K2=[R 0 0 -R 0 0;
              0 T (B+C1(i)*T) 0 -T (B+C2(i)*T);
              0 (B+C1(i)*T) (C+2*C1(i)*B+C1(i)^2*T) 0 -(B+C1(i)*T) (A+C1(i)*B+C2(i)*B+C1(i)*C2(i)*T);
              -R 0 0 R 0 0;
              0 -T -(B+C1(i)*T) 0 T -(B+C2(i)*T);
              0 (B+C2(i)*T) (A+C1(i)*B+C2(i)*B+C1(i)*C2(i)*T) 0 -(B+C2(i)*T) (C+2*C2(i)*B+C2(i)^2*T)];
              % MATRIZ DE TRANSPORTE A COORDENADAS GLOBALES
              T2_3=zeros(6,6);
              T2_3(1,1)=COSENO(i);T2_3(1,2)=SENO(i);
              T2_3(2,1)=-SENO(i);T2_3(2,2)=COSENO(i);
              T2_3(3,3)=1;
              T2_3(4,4)=COSENO(i);T2_3(4,5)=SENO(i);
              T2_3(5,4)=-SENO(i);T2_3(5,5)=COSENO(i);
              T2_3(6,6)=1;
              % ELEMENTO CIRCULAR SÓLIDO DE SECCIÓN VARIABLE
          elseif DD(i)==5
              Di(i)=input('\n DIÁMETRO INICIAL DEL ELEMENTO [m]:');
              Df(i)=input('\n DIÁMETRO FINAL DEL ELEMENTO [m]:');
              BET(i)=input('\n VALOR DEL FACTOR DE FORMA (RECTANGULAR=1.2; CIRCULAR=10/9):');
              % ELEMENTOS UNICAMENTE NECESARIOS PARA EL CONTEO DURANTE EL USO
              % DE ELEMENTOS DE CUALQUIERA DE LOS DISTINTOS TIPOS POSIBLES DE
              % SECCIONES
              B(i)=Di(i);
              H(i)=Df(i);
              Bi(i)=B(i);
              Hi(i)=H(i);
              Bf(i)=Bi(i);
              Hf(i)=Hi(i);
              BET(i)=10/9;
              C1(i)=0;
              C2(i)=0;
              F(i)=L(i);
              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              AREA(i)=pi*Di(i)*Di(i)/4;
              INER(i)=pi*(Di(i)/2)^4/4;
              FIx(i)=(12*ELAS*INER(i)*BET(i))/(GCOR*AREA(i)*L(i)^2)
              FIy(i)=FIx(i);
              FI(i)=(1/4)*FIx(i);
              % ELEMENTOS DE LA MATRIZ DE FLEXIBILIDADES
              f11=(4*L(i)/(pi*ELAS*(Di(i))^2))*(Di(i)/(Df(i)-Di(i)))*(1-Di(i)/Df(i));
              f22=(64*(L(i))^3/(3*pi*ELAS*(Di(i))^4))*(Di(i)/Df(i))^3+(40*L(i)/(9*pi*GCOR*(Di(i))^2))*(Di(i)/(Df(i)-Di(i)))*(1-(Di(i)/Df(i)));
              f33=f22;
              f26=(64*(L(i))^2/(3*pi*ELAS*(Di(i))^4))*((Di(i)/Df(i))^3+(1/2)*(Di(i)/(Df(i)-Di(i)))^2*(1+(Di(i)/Df(i))^2-2*(Di(i)/Df(i))));
              f35=f26;
              f66=(64*L(i)/(3*pi*ELAS*(Di(i))^4))*(Di(i)/(Df(i)-Di(i)))*(1-(Di(i)/Df(i))^3);
              f55=f66;
              f44=(32*L(i)/(3*pi*GCOR*Di(i)^4))*(Di(i)/(Df(i)-Di(i)))*(1-(Di(i)/Df(i))^3);
              FF44(i)=f44;
              % ELEMENTOS DE LA MATRIZ DE RIGIDECES
              raz=1/f11;
              Detx=f22*f66-f26*f26;
              r11x=f22/Detx;
              r12x=(f26*L(i)-f22)/Detx;
              r22x=(f66*(L(i))^2-2*f26*L(i)+f22)/Detx;
              raax=(r11x+r22x+2*r12x)/(L(i))^2;
              rabx=(r11x+r12x)/L(i);
              rbax=(r22x+r12x)/L(i);
              % MATRIZ DE RIGIDECES EN COORDENADAS LOCALES
              K11=[raz 0 0 ; 0 raax rabx ; 0 rabx r11x];
              K12=[-raz 0 0 ; 0 -raax rbax ; 0 -rabx r12x];
              K21=K12';
              K22=[raz 0 0 ; 0 raax -rbax ; 0 -rbax r22x];
              K2=[K11 K12 ;K21 K22];
              % MATRIZ DE TRANSPORTE A COORDENADAS GLOBALES
              T2_3=zeros(6,6);
              T2_3(1,1)=COSENO(i);T2_3(1,2)=SENO(i);
              T2_3(2,1)=-SENO(i);T2_3(2,2)=COSENO(i);
              T2_3(3,3)=1;
              T2_3(4,4)=COSENO(i);T2_3(4,5)=SENO(i);
              T2_3(5,4)=-SENO(i);T2_3(5,5)=COSENO(i);
              T2_3(6,6)=1;
          end
       elseif i==mbr
           fprintf('\n MIEMBRO %d:',i);
           fprintf('\n\n ESCOJA UNA DE LAS SIGUIENTES OPCIONES:');
           fprintf('\n\n\t\t (1) MIEMBRO PRISMÁTICO DE SECCIÓN RECTANGULAR');
           fprintf('\n\n\t\t (2) MIEMBRO CUADRADO DE SECCIÓN VARIABLE');
           fprintf('\n\n\t\t (3) MIEMBRO RECTANGULAR DE SECCIÓN VARIABLE');
           fprintf('\n\n\t\t (4) MIEMBRO CIRCULAR DE SECCIÓN PRISMÁTICA CONSTANTE');
           fprintf('\n\n\t\t (5) MIEMBRO CIRCULAR DE SECCIÓN VARIABLE');
           DD(i) = input('\n\n OPCION : ');
           fprintf('\n\n\n'); 
           % ELEMENTO PRISMÁTICO DE SECCIÓN CONSTANTE DE FORMA RECTANGULAR
           if DD(i)==1
               B(i)=input('\n BASE DEL ELEMENTO [m]:');
               H(i)=input(' ALTURA DEL ELEMENTO [m]:');
               % ELEMENTOS UNICAMENTE NECESARIOS PARA EL CONTEO DURANTE EL USO
               % DE ELEMENTOS DE CUALQUIERA DE LOS DISTINTOS TIPOS POSIBLES DE
               % SECCIONES
               Bi(i)=B(i);
               Hi(i)=H(i);
               Bf(i)=Bi(i);
               Hf(i)=Hi(i);
               Di(i)=B(i);
               Df(i)=Di(i);
               FF44(i)=0;
               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               BET(i)=input('\n VALOR DEL FACTOR DE FORMA (RECTANGULAR=1.2; CIRCULAR=10/9):');
               C1(i)=input('\n\n INGRESE EL VALOR DE C1 [m] (LONGITUD DEL SEGMENTO DE RIGIDEZ INFINITA EN EL NODO INICIAL):');
               C2(i)=input('\n\n INGRESE EL VALOR DE C2 [m] (LONGITUD DEL SEGMENTO DE RIGIDEZ INFINITA EN EL NODO FINAL):');
               F(i)=L(i)-C1(i)-C2(i)
               % CONSTANTES AUXILIARES PARA FUTUROS DESARROLLOS EN 3D
               AREA(i)=B(i)*H(i);
               INERx(i)=(B(i)*H(i)^3)/12;
               INER(i)=INERx(i);
               INERy(i)=(H(i)*B(i)^3)/12;
               FIx(i)=(12*ELAS*INERx(i)*BET(i))/(GCOR*AREA(i)*L(i)^2);
               FIy(i)=(12*ELAS*INERy(i)*BET(i))/(GCOR*AREA(i)*L(i)^2);
               % GENERACIÓN DE LAS VARIABLES
               FI(i)=(1/4)*FIx(i);
               C=((4*ELAS*INER(i))/L(i))*((1+FI(i))/(1+4*FI(i)));
               CP=C;
               A=((2*ELAS*INER(i))/L(i))*((1-2*FI(i))/(1+4*FI(i)));
               B=((6*ELAS*INER(i))/L(i)^2)*(1/(1+4*FI(i)));
               BP=B;
               T=(12*ELAS*INER(i))/L(i)^3*(1/(1+4*FI(i)));
               R=ELAS*AREA(i)/L(i);
               % MATRIZ DE RIGIDECES EN COORDENADAS LOCALES
               K2=[R 0 0 -R 0 0;
               0 T (B+C1(i)*T) 0 -T (B+C2(i)*T);
               0 (B+C1(i)*T) (C+2*C1(i)*B+C1(i)^2*T) 0 -(B+C1(i)*T) (A+C1(i)*B+C2(i)*B+C1(i)*C2(i)*T);
               -R 0 0 R 0 0;
               0 -T -(B+C1(i)*T) 0 T -(B+C2(i)*T);
               0 (B+C2(i)*T) (A+C1(i)*B+C2(i)*B+C1(i)*C2(i)*T) 0 -(B+C2(i)*T) (C+2*C2(i)*B+C2(i)^2*T)];
               % MATRIZ DE TRANSPORTE A COORDENADAS GLOBALES
               T2_3=zeros(6,6);
               T2_3(1,1)=COSENO(i);T2_3(1,2)=SENO(i);
               T2_3(2,1)=-SENO(i);T2_3(2,2)=COSENO(i);
               T2_3(3,3)=1;
               T2_3(4,4)=COSENO(i);T2_3(4,5)=SENO(i);
               T2_3(5,4)=-SENO(i);T2_3(5,5)=COSENO(i);
               T2_3(6,6)=1;
               % ELEMENTO CUADRADO SÓLIDO DE SECCIÓN VARIABLE
           elseif DD(i)==2
               Bi(i)=input('\n BASE INICIAL [m]:');
               Hi(i)=Bi(i);
               Bf(i)=input('\n BASE FINAL [m]:');
              % ELEMENTOS UNICAMENTE NECESARIOS PARA EL CONTEO DURANTE EL USO
              % DE ELEMENTOS DE CUALQUIERA DE LOS DISTINTOS TIPOS POSIBLES DE
              % SECCIONES
              Hf(i)=Bf(i);
              B(i)=Bi(i);
              H(i)=Hi(i);
              Di(i)=Bi(i);
              Df(i)=Di(i);
              BET(i)=1.2;
              INER(i)=(B(i)*H(i)^3)/12;
              FI(i)=1;
              AREA(i)=B(i)*H(i);
              C1(i)=0;
              C2(i)=0;
              F(i)=L(i);
              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              % ELEMENTOS DE LA MATRIZ DE FLEXIBILIDADES
              f11=(L(i)/(ELAS*Hi(i)*Hi(i)))*(Hi(i)/(Hf(i)-Hi(i)))*(1-Hi(i)/Hf(i));
              f22=(4*(L(i))^3/(ELAS*(Hi(i))^4))*(Hi(i)/Hf(i))^3+(6*L(i)/(5*GCOR*(Hi(i))^2))*(Hi(i)/(Hf(i)-Hi(i)))*(1-Hi(i)/Hf(i));
              f33=f22;
              f26=(4*(L(i))^2/(ELAS*(Hi(i))^4))*((Hi(i)/Hf(i))^3+(1/2)*(Hi(i)/(Hf(i)-Hi(i)))^2*(1+(Hi(i)/Hf(i))^2-2*(Hi(i)/Hf(i))));
              f35=f26;
              f66=(4*L(i)/(ELAS*(Hi(i))^4))*(Hi(i)/(Hf(i)-Hi(i)))*(1-(Hi(i)/Hf(i))^3);
              f55=f66;
              f44=(64*L(i)/(27*GCOR*(Hi(i))^4))*(Hi(i)/(Hf(i)-Hi(i)))*(1-(Hi(i)/Hf(i))^3);
              FF44(i)=f44;
              % ELEMENTOS DE LA MATRIZ DE RIGIDECES
              raz=1/f11;
              Detx=f22*f66-f26*f26;
              r11x=f22/Detx;
              r12x=(f26*L(i)-f22)/Detx;
              r22x=(f66*(L(i))^2-2*f26*L(i)+f22)/Detx;
              raax=(r11x+r22x+2*r12x)/(L(i))^2;
              rabx=(r11x+r12x)/L(i);
              rbax=(r22x+r12x)/L(i);
              % MATRIZ DE RIGIDECES EN COORDENADAS LOCALES
              K11=[raz 0 0 ; 0 raax rabx ; 0 rabx r11x];
              K12=[-raz 0 0 ; 0 -raax rbax ; 0 -rabx r12x];
              K21=K12';
              K22=[raz 0 0 ; 0 raax -rbax ; 0 -rbax r22x];
              K2=[K11 K12 ;K21 K22];
              % MATRIZ DE TRANSPORTE A COORDENADAS GLOBALES
              T2_3=zeros(6,6);
              T2_3(1,1)=COSENO(i);T2_3(1,2)=SENO(i);
              T2_3(2,1)=-SENO(i);T2_3(2,2)=COSENO(i);
              T2_3(3,3)=1;
              T2_3(4,4)=COSENO(i);T2_3(4,5)=SENO(i);
              T2_3(5,4)=-SENO(i);T2_3(5,5)=COSENO(i);
              T2_3(6,6)=1;
              % ELEMENTO RECTANGULAR SÓLIDO DE SECCIÓN VARIABLE
           elseif DD(i)==3
               Bi(i)=input('\n BASE DE TODO EL ELEMENTO [m]:');
               Hi(i)=input(' ALTURA INICIAL DEL ELEMENTO [m]:');
               Bf(i)=Bi(i);
               Hf(i)=input(' ALTURA FINAL DEL ELEMENTO [m]:');
               Bb=Bi(i);
               % ELEMENTOS UNICAMENTE NECESARIOS PARA EL CONTEO DURANTE EL USO
               % DE ELEMENTOS DE CUALQUIERA DE LOS DISTINTOS TIPOS POSIBLES DE
               % SECCIONES
               B(i)=Bi(i);
               H(i)=Hi(i);
               Di(i)=Bi(i);
               Df(i)=Di(i);
               INER(i)=(B(i)*H(i)^3)/12;
               FI(i)=1;
               AREA(i)=B(i)*H(i);
               BET(i)=1.2;
               C1(i)=0;
               C2(i)=0;
               F(i)=L(i);
               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               BET(i)=input('\n FACTOR DE FORMA (RECTANGULAR=1.2; CIRCULAR=10/9):');
               fprintf('\n\t\t\t ‘ ; ......................................................................................................................... ');
               fprintf('\n\t\t\t ‘ ; ......................................................................................................................... ');
               fprintf('\n\n\t\t     LA FUNCIÓN F(x) INDICA LA VARIACIÓN CON LA LONGITUD DE LA FORMA DEL ELEMENTO,');
               fprintf('\n\t\t\t                                                                    ');
               fprintf('\n\n\t\t     SI LA FUNCIÓN ES LINEAL: Y=F(x)=');
               fprintf('\n\t\t\t                                      (abs(Hi-Hf)/LL)*x                                       ');
               fprintf('\n\t\t\t ‘ ; ......................................................................................................................... ');
               fprintf('\n\t\t\t ‘ ; ......................................................................................................................... ');
               fprintf('\n\t\t\t ‘ ;                                                                                                                           ');
               fprintf('\n\t\t\t   ESTANDO EL EJE X SOBRE EL     Y|                                                                                            ');
               fprintf('\n\t\t\t    CENTROIDE DEL MIEMBRO         |                                                                                            ');
               fprintf('\n\t\t\t                                  |                                                                                            ');
               fprintf('\n\t\t\t                                  |_ _ _ _ _                                                                                   ');
               fprintf('\n\t\t\t                                           X                                                                                   ');
               fprintf('\n\t\t\t                                                                                                                               ');
               fprintf('\n\t\t\t ‘ ; ......................................................................................................................... ');
               fprintf('\n\t\t\t ‘ ; ......................................................................................................................... ');
               fprintf('\n\t\t\t                                                                                                                              ');
               fprintf('\n\t\t\t                                                                                                                              ');
               EE(i)=menu('ESCOJA UNA DE LAS SIGUIENTES OPCIONES:','FUNCIÓN LINEAL','OTRO TIPO DE VARIACIÓN');
               if EE(i)==1
                   % ELEMENTOS DE LA MATRIZ DE FLEXIBILIDADES
                   f11=(L(i)/(ELAS*Bb*Hi(i)))*(Hi(i)/(Hf(i)-Hi(i)))*log(Hf(i)/Hi(i));
                   f22=(6*(L(i))^3/(ELAS*Bb*(Hi(i))^3))*((Hi(i)/Hf(i))^2-(Hi(i)/(Hf(i)-Hi(i)))^3*(Hf(i)/Hi(i)-Hi(i)/Hf(i)-2*log(Hf(i)/Hi(i))))+(6*L(i)/(5*GCOR*Bb*Hi(i)))*(Hi(i)/(Hf(i)-Hi(i)))*log(Hf(i)/Hi(i));
                   f33=(12*(L(i))^3/(ELAS*(Bb)^3*Hi(i)))*(3/2+(1/2)*(Hf(i)/Hi(i))^2-2*(Hf(i)/Hi(i))+log(Hf(i)/Hi(i)))+(6*L(i)/(5*GCOR*Bb*Hi(i)))*(Hi(i)/(Hf(i)-Hi(i)))*log(Hf(i)/Hi(i));
                   f26=(6*L(i)/(ELAS*Bb*(Hi(i))^3))*(Hi(i)/(Hf(i)-Hi(i)))^2*(1+(Hi(i)/Hf(i))^2-2*(Hi(i)/Hf(i)));
                   f35=(12*(L(i))^2/(ELAS*(Bb)^3*Hi(i)))*(Hi(i)/(Hf(i)-Hi(i)))^2*(Hf(i)/Hi(i)-log(Hf(i)/Hi(i))-1);
                   f66=(6*L(i)/(ELAS*Bb*(Hi(i))^3))*(Hi(i)/(Hf(i)-Hi(i)))*(1-(Hi(i)/Hf(i))^2);
                   f55=(12*L(i)/(ELAS*(Bb)^3*Hi(i)))*(Hi(i)/(Hf(i)-Hi(i)))*log(Hf(i)/Hi(i));
                   % DEFINICIÓN DE UNA FUNCIÓN DE x
                   Hi=Hi(i);
                   Hf=Hf(i);
                   LL=L(i);
                   if Hi<Hf
                       dpc= @(x)(1./(((Bb.^3*(-Hi./2-((Hi-Hf)./(2*LL))*x))./3)*(1-(192*Bb./(pi.^5*(-Hi./2-((Hi-Hf)./(2*LL))*x)))*(tanh(pi*(-Hi./2-((Hi-Hf)./(2*LL))*x)./(2*Bb)))))); 
                   elseif Hi>Hf
                       dpc= @(x)(1./(((Bb.^3*(-Hi./2+((Hi-Hf)./(2*LL))*x))./3)*(1-(192*Bb./(pi.^5*(-Hi./2+((Hi-Hf)./(2*LL))*x)))*(tanh(pi*(-Hi./2+((Hi-Hf)./(2*LL))*x)./(2*Bb))))));
                   end
                   f44=quad(@(x)dpc, 0, LL)
                   FF44(i)=f44;
                   % ELEMENTOS DE LA MATRIZ DE RIGIDECES
                   raz=1/f11;
                   Detx=f22*f66-f26*f26;
                   r11x=f22/Detx;
                   r12x=(f26*L(i)-f22)/Detx;
                   r22x=(f66*(L(i))^2-2*f26*L(i)+f22)/Detx;
                   raax=(r11x+r22x+2*r12x)/(L(i))^2;
                   rabx=(r11x+r12x)/L(i);
                   rbax=(r22x+r12x)/L(i);
                   % MATRIZ DE RIGIDECES EN COORDENADAS LOCALES
                   K11=[raz 0 0 ; 0 raax rabx ; 0 rabx r11x];
                   K12=[-raz 0 0 ; 0 -raax rbax ; 0 -rabx r12x];
                   K21=K12';
                   K22=[raz 0 0 ; 0 raax -rbax ; 0 -rbax r22x];
                   K2=[K11 K12 ;K21 K22];
                   % MATRIZ DE TRANSPORTE A COORDENADAS GLOBALES
                   T2_3=zeros(6,6);
                   T2_3(1,1)=COSENO(i);T2_3(1,2)=SENO(i);
                   T2_3(2,1)=-SENO(i);T2_3(2,2)=COSENO(i);
                   T2_3(3,3)=1;
                   T2_3(4,4)=COSENO(i);T2_3(4,5)=SENO(i);
                   T2_3(5,4)=-SENO(i);T2_3(5,5)=COSENO(i);
                   T2_3(6,6)=1;
               else
                   % ELEMENTOS DE LA MATRIZ DE FLEXIBILIDADES
                   f11=(L(i)/(ELAS*Bb*Hi(i)))*(Hi(i)/(Hf(i)-Hi(i)))*log(Hf(i)/Hi(i));
                   f22=(6*(L(i))^3/(ELAS*Bb*(Hi(i))^3))*((Hi(i)/Hf(i))^2-(Hi(i)/(Hf(i)-Hi(i)))^3*(Hf(i)/Hi(i)-Hi(i)/Hf(i)-2*log(Hf(i)/Hi(i))))+(6*L(i)/(5*GCOR*Bb*Hi(i)))*(Hi(i)/(Hf(i)-Hi(i)))*log(Hf(i)/Hi(i));
                   f33=(12*(L(i))^3/(ELAS*(Bb)^3*Hi(i)))*(3/2+(1/2)*(Hf(i)/Hi(i))^2-2*(Hf(i)/Hi(i))+log(Hf(i)/Hi(i)))+(6*L(i)/(5*GCOR*Bb*Hi(i)))*(Hi(i)/(Hf(i)-Hi(i)))*log(Hf(i)/Hi(i));
                   f26=(6*L(i)/(ELAS*Bb*(Hi(i))^3))*(Hi(i)/(Hf(i)-Hi(i)))^2*(1+(Hi(i)/Hf(i))^2-2*(Hi(i)/Hf(i)));
                   f35=(12*(L(i))^2/(ELAS*(Bb)^3*Hi(i)))*(Hi(i)/(Hf(i)-Hi(i)))^2*(Hf(i)/Hi(i)-log(Hf(i)/Hi(i))-1);
                   f66=(6*L(i)/(ELAS*Bb*(Hi(i))^3))*(Hi(i)/(Hf(i)-Hi(i)))*(1-(Hi(i)/Hf(i))^2);
                   f55=(12*L(i)/(ELAS*(Bb)^3*Hi(i)))*(Hi(i)/(Hf(i)-Hi(i)))*log(Hf(i)/Hi(i));
                   % DEFINICIÓN DE UNA FUNCIÓN DE x
                   Hi=Hi(i);
                   Hf=Hf(i);
                   LL=L(i);
                   f44=input(' INTRODUZCA LA f44=');
                   FF44(i)=f44;
                   % ELEMENTOS DE LA MATRIZ DE RIGIDECES
                   raz=1/f11;
                   Detx=f22*f66-f26*f26;
                   r11x=f22/Detx;
                   r12x=(f26*L(i)-f22)/Detx;
                   r22x=(f66*(L(i))^2-2*f26*L(i)+f22)/Detx;
                   raax=(r11x+r22x+2*r12x)/(L(i))^2;
                   rabx=(r11x+r12x)/L(i);
                   rbax=(r22x+r12x)/L(i);
                   % MATRIZ DE RIGIDECES EN COORDENADAS LOCALES
                   K11=[raz 0 0 ; 0 raax rabx ; 0 rabx r11x];
                   K12=[-raz 0 0 ; 0 -raax rbax ; 0 -rabx r12x];
                   K21=K12';
                   K22=[raz 0 0 ; 0 raax -rbax ; 0 -rbax r22x];
                   K2=[K11 K12 ;K21 K22];
                   % MATRIZ DE TRANSPORTE A COORDENADAS GLOBALES
                   T2_3=zeros(6,6);
                   T2_3(1,1)=COSENO(i);T2_3(1,2)=SENO(i);
                   T2_3(2,1)=-SENO(i);T2_3(2,2)=COSENO(i);
                   T2_3(3,3)=1;
                   T2_3(4,4)=COSENO(i);T2_3(4,5)=SENO(i);
                   T2_3(5,4)=-SENO(i);T2_3(5,5)=COSENO(i);
                   T2_3(6,6)=1;
               end
              % ELEMENTO CIRCULAR SÓLIDO DE SECCIÓN CONSTATE
           elseif DD(i)==4
               Di(i)=input('\n DIÁMETRO DEL ELEMENTO [m]:');
               % ELEMENTOS UNICAMENTE NECESARIOS PARA EL CONTEO DURANTE EL USO
               % DE ELEMENTOS DE CUALQUIERA DE LOS DISTINTOS TIPOS POSIBLES DE
               % SECCIONES
               Df(i)=Di(i);
               B(i)=Di(i);
               H(i)=Di(i);
               Bi(i)=B(i);
               Hi(i)=H(i);
               Bf(i)=Bi(i);
               Hf(i)=Hi(i);
               FF44(i)=0;
               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               BET(i)=input('\n VALOR DEL FACTOR DE FORMA (RECTANGULAR=1.2; CIRCULAR=10/9):');
               AREA(i)=pi*Di(i)*Di(i)/4;
               INER(i)=pi*(Di(i)/2)^4/4;
               C1(i)=input('\n\n INGRESE EL VALOR DE C1 [m] (LONGITUD DEL SEGMENTO DE RIGIDEZ INFINITA EN EL NODO INICIAL):');
               C2(i)=input('\n\n INGRESE EL VALOR DE C2 [m] (LONGITUD DEL SEGMENTO DE RIGIDEZ INFINITA EN EL NODO FINAL):');
               F(i)=L(i)-C1(i)-C2(i)
               FI(i)=(3*ELAS*INER(i)*BET(i))/(GCOR*AREA(i)*L(i)^2);
               C=((4*ELAS*INER(i))/L(i))*((1+FI(i))/(1+4*FI(i)));
               CP=C;
               A=((2*ELAS*INER(i))/L(i))*((1-2*FI(i))/(1+4*FI(i)));
               B=((6*ELAS*INER(i))/L(i)^2)*(1/(1+4*FI(i)));
               BP=B;
               T=(12*ELAS*INER(i))/L(i)^3*(1/(1+4*FI(i)));
               R=ELAS*AREA(i)/L(i);
               K2=[R 0 0 -R 0 0;
               0 T (B+C1(i)*T) 0 -T (B+C2(i)*T);
               0 (B+C1(i)*T) (C+2*C1(i)*B+C1(i)^2*T) 0 -(B+C1(i)*T) (A+C1(i)*B+C2(i)*B+C1(i)*C2(i)*T);
               -R 0 0 R 0 0;
               0 -T -(B+C1(i)*T) 0 T -(B+C2(i)*T);
               0 (B+C2(i)*T) (A+C1(i)*B+C2(i)*B+C1(i)*C2(i)*T) 0 -(B+C2(i)*T) (C+2*C2(i)*B+C2(i)^2*T)];
              % MATRIZ DE TRANSPORTE A COORDENADAS GLOBALES
               T2_3=zeros(6,6);
               T2_3(1,1)=COSENO(i);T2_3(1,2)=SENO(i);
               T2_3(2,1)=-SENO(i);T2_3(2,2)=COSENO(i);
               T2_3(3,3)=1;
               T2_3(4,4)=COSENO(i);T2_3(4,5)=SENO(i);
               T2_3(5,4)=-SENO(i);T2_3(5,5)=COSENO(i);
               T2_3(6,6)=1;
               % ELEMENTO CIRCULAR SÓLIDO DE SECCIÓN VARIABLE
           elseif DD(i)==5
               Di(i)=input('\n DIÁMETRO INICIAL DEL ELEMENTO [m]:');
               Df(i)=input('\n DIÁMETRO FINAL DEL ELEMENTO [m]:');
               BET(i)=input('\n VALOR DEL FACTOR DE FORMA (RECTANGULAR=1.2; CIRCULAR=10/9):');
               % ELEMENTOS UNICAMENTE NECESARIOS PARA EL CONTEO DURANTE EL USO
               % DE ELEMENTOS DE CUALQUIERA DE LOS DISTINTOS TIPOS POSIBLES DE
               % SECCIONES
               B(i)=Di(i);
               H(i)=Df(i);
               Bi(i)=B(i);
               Hi(i)=H(i);
               Bf(i)=Bi(i);
               Hf(i)=Hi(i);
               C1(i)=0;
               C2(i)=0;
               F(i)=L(i);
               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               AREA(i)=pi*Di(i)*Di(i)/4;
               INER(i)=pi*(Di(i)/2)^4/4;
               FIx(i)=(12*ELAS*INER(i)*BET(i))/(GCOR*AREA(i)*L(i)^2)
               FIy(i)=FIx(i);
               FI(i)=(1/4)*FIx(i);
               % ELEMENTOS DE LA MATRIZ DE FLEXIBILIDADES
               f11=(4*L(i)/(pi*ELAS*(Di(i))^2))*(Di(i)/(Df(i)-Di(i)))*(1-Di(i)/Df(i));
               f22=(64*(L(i))^3/(3*pi*ELAS*(Di(i))^4))*(Di(i)/Df(i))^3+(40*L(i)/(9*pi*GCOR*(Di(i))^2))*(Di(i)/(Df(i)-Di(i)))*(1-(Di(i)/Df(i)));
               f33=f22;
               f26=(64*(L(i))^2/(3*pi*ELAS*(Di(i))^4))*((Di(i)/Df(i))^3+(1/2)*(Di(i)/(Df(i)-Di(i)))^2*(1+(Di(i)/Df(i))^2-2*(Di(i)/Df(i))));
               f35=f26;
               f66=(64*L(i)/(3*pi*ELAS*(Di(i))^4))*(Di(i)/(Df(i)-Di(i)))*(1-(Di(i)/Df(i))^3);
               f55=f66;
               f44=(32*L(i)/(3*pi*GCOR*Di(i)^4))*(Di(i)/(Df(i)-Di(i)))*(1-(Di(i)/Df(i))^3);
               FF44(i)=f44;
               % ELEMENTOS DE LA MATRIZ DE RIGIDECES
               raz=1/f11;
               Detx=f22*f66-f26*f26;
               r11x=f22/Detx;
               r12x=(f26*L(i)-f22)/Detx;
               r22x=(f66*(L(i))^2-2*f26*L(i)+f22)/Detx;
               raax=(r11x+r22x+2*r12x)/(L(i))^2;
               rabx=(r11x+r12x)/L(i);
               rbax=(r22x+r12x)/L(i);
               % MATRIZ DE TRANSPORTE A COORDENADAS GLOBALES
               K11=[raz 0 0 ; 0 raax rabx ; 0 rabx r11x];
               K12=[-raz 0 0 ; 0 -raax rbax ; 0 -rabx r12x];
               K21=K12';
               K22=[raz 0 0 ; 0 raax -rbax ; 0 -rbax r22x];
               K2=[K11 K12 ;K21 K22];
               % MATRIZ DE TRANSPORTE A COORDENADAS LOCALES
               T2_3=zeros(6,6);
               T2_3(1,1)=COSENO(i);T2_3(1,2)=SENO(i);
               T2_3(2,1)=-SENO(i);T2_3(2,2)=COSENO(i);
               T2_3(3,3)=1;
               T2_3(4,4)=COSENO(i);T2_3(4,5)=SENO(i);
               T2_3(5,4)=-SENO(i);T2_3(5,5)=COSENO(i);
               T2_3(6,6)=1;
           end
       else
           fprintf('\n MIEMBRO %d:',i);
           fprintf('\n\n ESCOJA UNA DE LAS SIGUIENTES OPCIONES:');
           fprintf('\n\n\t\t (1) EL MIEMBRO TIENE LAS MISMAS CARACTERÍSTICAS QUE EL MIEMBRO ANTERIOR');
           fprintf('\n\n\t\t (2) EL MIEMBRO TIENE CARACTERÍSTICAS DISTINTAS AL MIEMBRO ANTERIOR');
           NUB=input('\n\n OPCION : ');
           if NUB==1
               DD(i)=DD(i-1);        
               H(i)=H(i-1);
               B(i)=1;
               Bi(i)=Bi(i-1);
               Hi(i)=Hi(i-1);
               Bf(i)=Bf(i-1);
               Hf(i)=Hf(i-1);
               Di(i)=Di(i-1);
               Df(i)=Df(i-1);
               BET(i)=BET(i-1);
               C1(i)=C1(i-1);
               C2(i)=C2(i-1);
               F(i)=F(i-1);
               FI(i)=FI(i-1);
               INER(i)=INER(i-1);
               AREA(i)=AREA(i-1);
               FF44(i)=FF44(i-1);
               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               if DD(i)==1  
                   C=((4*ELAS*INER(i))/L(i))*((1+FI(i))/(1+4*FI(i)));
                   CP=C;
                   A=((2*ELAS*INER(i))/L(i))*((1-2*FI(i))/(1+4*FI(i)));
                   B=((6*ELAS*INER(i))/L(i)^2)*(1/(1+4*FI(i)));
                   BP=B;
                   T=(12*ELAS*INER(i))/L(i)^3*(1/(1+4*FI(i)));
                   R=ELAS*AREA(i)/L(i);
                   K2=[R 0 0 -R 0 0;
                   0 T (B+C1(i)*T) 0 -T (B+C2(i)*T);
                   0 (B+C1(i)*T) (C+2*C1(i)*B+C1(i)^2*T) 0 -(B+C1(i)*T) (A+C1(i)*B+C2(i)*B+C1(i)*C2(i)*T);
                   -R 0 0 R 0 0;
                   0 -T -(B+C1(i)*T) 0 T -(B+C2(i)*T);
                   0 (B+C2(i)*T) (A+C1(i)*B+C2(i)*B+C1(i)*C2(i)*T) 0 -(B+C2(i)*T) (C+2*C2(i)*B+C2(i)^2*T)];
                   % MATRIZ DE TRANSPORTE A COORDENADAS GLOBALES
                   T2_3=zeros(6,6);
                   T2_3(1,1)=COSENO(i);T2_3(1,2)=SENO(i);
                   T2_3(2,1)=-SENO(i);T2_3(2,2)=COSENO(i);
                   T2_3(3,3)=1;
                   T2_3(4,4)=COSENO(i);T2_3(4,5)=SENO(i);
                   T2_3(5,4)=-SENO(i);T2_3(5,5)=COSENO(i);
                   T2_3(6,6)=1;
                   % ELEMENTO CUADRADO SÓLIDO DE SECCIÓN VARIABLE
               elseif DD(i)==2
                   % ELEMENTOS DE LA MATRIZ DE FLEXIBILIDADES
                   f11=(L(i)/(ELAS*Hi(i)*Hi(i)))*(Hi(i)/(Hf(i)-Hi(i)))*(1-Hi(i)/Hf(i));
                   f22=(4*(L(i))^3/(ELAS*(Hi(i))^4))*(Hi(i)/Hf(i))^3+(6*L(i)/(5*GCOR*(Hi(i))^2))*(Hi(i)/(Hf(i)-Hi(i)))*(1-Hi(i)/Hf(i));
                   f33=f22;
                   f26=(4*(L(i))^2/(ELAS*(Hi(i))^4))*((Hi(i)/Hf(i))^3+(1/2)*(Hi(i)/(Hf(i)-Hi(i)))^2*(1+(Hi(i)/Hf(i))^2-2*(Hi(i)/Hf(i))));
                   f35=f26;
                   f66=(4*L(i)/(ELAS*(Hi(i))^4))*(Hi(i)/(Hf(i)-Hi(i)))*(1-(Hi(i)/Hf(i))^3);
                   f55=f66;
                   f44=(64*L(i)/(27*GCOR*(Hi(i))^4))*(Hi(i)/(Hf(i)-Hi(i)))*(1-(Hi(i)/Hf(i))^3);
                   % ELEMENTOS DE LA MATRIZ DE RIGIDECES
                   raz=1/f11;
                   Detx=f22*f66-f26*f26;
                   r11x=f22/Detx;
                   r12x=(f26*L(i)-f22)/Detx;
                   r22x=(f66*(L(i))^2-2*f26*L(i)+f22)/Detx;
                   raax=(r11x+r22x+2*r12x)/(L(i))^2;
                   rabx=(r11x+r12x)/L(i);
                   rbax=(r22x+r12x)/L(i);
                   % MATRIZ DE RIGIDECES EN COORDENADAS LOCALES
                   K11=[raz 0 0 ; 0 raax rabx ; 0 rabx r11x];
                   K12=[-raz 0 0 ; 0 -raax rbax ; 0 -rabx r12x];
                   K21=K12';
                   K22=[raz 0 0 ; 0 raax -rbax ; 0 -rbax r22x];
                   K2=[K11 K12 ;K21 K22];
                   % MATRIZ DE TRANSPORTE A COORDENADAS GLOBALES
                   T2_3=zeros(6,6);
                   T2_3(1,1)=COSENO(i);T2_3(1,2)=SENO(i);
                   T2_3(2,1)=-SENO(i);T2_3(2,2)=COSENO(i);
                   T2_3(3,3)=1;
                   T2_3(4,4)=COSENO(i);T2_3(4,5)=SENO(i);
                   T2_3(5,4)=-SENO(i);T2_3(5,5)=COSENO(i);
                   T2_3(6,6)=1;
                   % ELEMENTO RECTANGULAR SÓLIDO DE SECCIÓN VARIABLE
               elseif DD(i)==3
                   Bb=Bi(i);
                   % ELEMENTOS DE LA MATRIZ DE FLEXIBILIDADES
                   f11=(L(i)/(ELAS*Bb*Hi(i)))*(Hi(i)/(Hf(i)-Hi(i)))*log(Hf(i)/Hi(i));
                   f22=(6*(L(i))^3/(ELAS*Bb*(Hi(i))^3))*((Hi(i)/Hf(i))^2-(Hi(i)/(Hf(i)-Hi(i)))^3*(Hf(i)/Hi(i)-Hi(i)/Hf(i)-2*log(Hf(i)/Hi(i))))+(6*L(i)/(5*GCOR*Bb*Hi(i)))*(Hi(i)/(Hf(i)-Hi(i)))*log(Hf(i)/Hi(i));
                   f33=(12*(L(i))^3/(ELAS*(Bb)^3*Hi(i)))*(3/2+(1/2)*(Hf(i)/Hi(i))^2-2*(Hf(i)/Hi(i))+log(Hf(i)/Hi(i)))+(6*L(i)/(5*GCOR*Bb*Hi(i)))*(Hi(i)/(Hf(i)-Hi(i)))*log(Hf(i)/Hi(i));
                   f26=(6*L(i)/(ELAS*Bb*(Hi(i))^3))*(Hi(i)/(Hf(i)-Hi(i)))^2*(1+(Hi(i)/Hf(i))^2-2*(Hi(i)/Hf(i)));
                   f35=(12*(L(i))^2/(ELAS*(Bb)^3*Hi(i)))*(Hi(i)/(Hf(i)-Hi(i)))^2*(Hf(i)/Hi(i)-log(Hf(i)/Hi(i))-1);
                   f66=(6*L(i)/(ELAS*Bb*(Hi(i))^3))*(Hi(i)/(Hf(i)-Hi(i)))*(1-(Hi(i)/Hf(i))^2);
                   f55=(12*L(i)/(ELAS*(Bb)^3*Hi(i)))*(Hi(i)/(Hf(i)-Hi(i)))*log(Hf(i)/Hi(i));
                   f44=FF44(i);
                   % ELEMENTOS DE LA MATRIZ DE RIGIDECES
                   raz=1/f11;
                   Detx=f22*f66-f26*f26;
                   r11x=f22/Detx;
                   r12x=(f26*L(i)-f22)/Detx;
                   r22x=(f66*(L(i))^2-2*f26*L(i)+f22)/Detx;
                   raax=(r11x+r22x+2*r12x)/(L(i))^2;
                   rabx=(r11x+r12x)/L(i);
                   rbax=(r22x+r12x)/L(i);
                   % MATRIZ DE RIGIDECES EN COORDENADAS LOCALES
                   K11=[raz 0 0 ; 0 raax rabx ; 0 rabx r11x];
                   K12=[-raz 0 0 ; 0 -raax rbax ; 0 -rabx r12x];
                   K21=K12';
                   K22=[raz 0 0 ; 0 raax -rbax ; 0 -rbax r22x];
                   K2=[K11 K12 ;K21 K22];
                   % MATRIZ DE TRANSPORTE A COORDENADAS GLOBALES
                   T2_3=zeros(6,6);
                   T2_3(1,1)=COSENO(i);T2_3(1,2)=SENO(i);
                   T2_3(2,1)=-SENO(i);T2_3(2,2)=COSENO(i);
                   T2_3(3,3)=1;
                   T2_3(4,4)=COSENO(i);T2_3(4,5)=SENO(i);
                   T2_3(5,4)=-SENO(i);T2_3(5,5)=COSENO(i);
                   T2_3(6,6)=1;
                   % ELEMENTO CIRCULAR SÓLIDO DE SECCIÓN CONSTATE
               elseif DD(i)==4
                   C=((4*ELAS*INER(i))/L(i))*((1+FI(i))/(1+4*FI(i)));
                   CP=C;
                   A=((2*ELAS*INER(i))/L(i))*((1-2*FI(i))/(1+4*FI(i)));
                   B=((6*ELAS*INER(i))/L(i)^2)*(1/(1+4*FI(i)));
                   BP=B;
                   T=(12*ELAS*INER(i))/L(i)^3*(1/(1+4*FI(i)));
                   R=ELAS*AREA(i)/L(i);
                   K2=[R 0 0 -R 0 0;
                   0 T (B+C1(i)*T) 0 -T (B+C2(i)*T);
                   0 (B+C1(i)*T) (C+2*C1(i)*B+C1(i)^2*T) 0 -(B+C1(i)*T) (A+C1(i)*B+C2(i)*B+C1(i)*C2(i)*T);
                   -R 0 0 R 0 0;
                   0 -T -(B+C1(i)*T) 0 T -(B+C2(i)*T);
                   0 (B+C2(i)*T) (A+C1(i)*B+C2(i)*B+C1(i)*C2(i)*T) 0 -(B+C2(i)*T) (C+2*C2(i)*B+C2(i)^2*T)];
                   % MATRIZ DE TRANSPORTE A COORDENADAS GLOBALES
                   T2_3=zeros(6,6);
                   T2_3(1,1)=COSENO(i);T2_3(1,2)=SENO(i);
                   T2_3(2,1)=-SENO(i);T2_3(2,2)=COSENO(i);
                   T2_3(3,3)=1;
                   T2_3(4,4)=COSENO(i);T2_3(4,5)=SENO(i);
                   T2_3(5,4)=-SENO(i);T2_3(5,5)=COSENO(i);
                   T2_3(6,6)=1;
                   % ELEMENTO CIRCULAR SÓLIDO DE SECCIÓN VARIABLE
               elseif DD(i)==5
                   % ELEMENTOS DE LA MATRIZ DE FLEXIBILIDADES
                   f11=(4*L(i)/(pi*ELAS*(Di(i))^2))*(Di(i)/(Df(i)-Di(i)))*(1-Di(i)/Df(i));
                   f22=(64*(L(i))^3/(3*pi*ELAS*(Di(i))^4))*(Di(i)/Df(i))^3+(40*L(i)/(9*pi*GCOR*(Di(i))^2))*(Di(i)/(Df(i)-Di(i)))*(1-(Di(i)/Df(i)));
                   f33=f22;
                   f26=(64*(L(i))^2/(3*pi*ELAS*(Di(i))^4))*((Di(i)/Df(i))^3+(1/2)*(Di(i)/(Df(i)-Di(i)))^2*(1+(Di(i)/Df(i))^2-2*(Di(i)/Df(i))));
                   f35=f26;
                   f66=(64*L(i)/(3*pi*ELAS*(Di(i))^4))*(Di(i)/(Df(i)-Di(i)))*(1-(Di(i)/Df(i))^3);
                   f55=f66;
                   f44=(32*L(i)/(3*pi*GCOR*Di(i)^4))*(Di(i)/(Df(i)-Di(i)))*(1-(Di(i)/Df(i))^3);
                   % ELEMENTOS DE LA MATRIZ DE RIGIDECES
                   raz=1/f11;
                   Detx=f22*f66-f26*f26;
                   r11x=f22/Detx;
                   r12x=(f26*L(i)-f22)/Detx;
                   r22x=(f66*(L(i))^2-2*f26*L(i)+f22)/Detx;
                   raax=(r11x+r22x+2*r12x)/(L(i))^2;
                   rabx=(r11x+r12x)/L(i);
                   rbax=(r22x+r12x)/L(i);
                   % MATRIZ DE RIGIDECES EN COORDENADAS LOCALES
                   K11=[raz 0 0 ; 0 raax rabx ; 0 rabx r11x];
                   K12=[-raz 0 0 ; 0 -raax rbax ; 0 -rabx r12x];
                   K21=K12';
                   K22=[raz 0 0 ; 0 raax -rbax ; 0 -rbax r22x];
                   K2=[K11 K12 ;K21 K22];
                   % MATRIZ DE TRANSPORTE A COORDENADAS GLOBALES
                   T2_3=zeros(6,6);
                   T2_3(1,1)=COSENO(i);T2_3(1,2)=SENO(i);
                   T2_3(2,1)=-SENO(i);T2_3(2,2)=COSENO(i);
                   T2_3(3,3)=1;
                   T2_3(4,4)=COSENO(i);T2_3(4,5)=SENO(i);
                   T2_3(5,4)=-SENO(i);T2_3(5,5)=COSENO(i);
                   T2_3(6,6)=1;
               end
           elseif NUB==2
               fprintf('\n MIEMBRO %d:',i);
               fprintf('\n\n ESCOJA UNA DE LAS SIGUIENTES OPCIONES:');
               fprintf('\n\n\t\t (1) MIEMBRO PRISMÁTICO DE SECCIÓN RECTANGULAR');
               fprintf('\n\n\t\t (2) MIEMBRO CUADRADO DE SECCIÓN VARIABLE');
               fprintf('\n\n\t\t (3) MIEMBRO RECTANGULAR DE SECCIÓN VARIABLE');
               fprintf('\n\n\t\t (4) MIEMBRO CIRCULAR DE SECCIÓN PRISMÁTICA CONSTANTE');
               fprintf('\n\n\t\t (5) MIEMBRO CIRCULAR DE SECCIÓN VARIABLE');
               DD(i) = input('\n\n OPCION : ');
               fprintf('\n\n\n'); 
               % ELEMENTO PRISMÁTICO DE SECCIÓN CONSTANTE DE FORMA RECTANGULAR
              if DD(i)==1
                  B(i)=input('\n BASE DEL ELEMENTO [m]:');
                  H(i)=input(' ALTURA DEL ELEMENTO [m]:');
                  % ELEMENTOS UNICAMENTE NECESARIOS PARA EL CONTEO DURANTE EL USO
                  % DE ELEMENTOS DE CUALQUIERA DE LOS DISTINTOS TIPOS POSIBLES DE
                  % SECCIONES
                  Bi(i)=B(i);
                  Hi(i)=H(i);
                  Bf(i)=Bi(i);
                  Hf(i)=Hi(i);
                  Di(i)=B(i);
                  Df(i)=Di(i);
                  FF44(i)=0;
                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                  BET(i)=input('\n VALOR DEL FACTOR DE FORMA (RECTANGULAR=1.2; CIRCULAR=10/9):');
                  C1(i)=input('\n\n INGRESE EL VALOR DE C1 [m] (LONGITUD DEL SEGMENTO DE RIGIDEZ INFINITA EN EL NODO INICIAL):');
                  C2(i)=input('\n\n INGRESE EL VALOR DE C2 [m] (LONGITUD DEL SEGMENTO DE RIGIDEZ INFINITA EN EL NODO FINAL):');
                  F(i)=L(i)-C1(i)-C2(i)
                  % CONSTANTES AUXILIARES PARA FUTUROS DESARROLLOS EN 3D
                  AREA(i)=B(i)*H(i);
                  INERx(i)=(B(i)*H(i)^3)/12;
                  INER(i)=INERx(i);
                  INERy(i)=(H(i)*B(i)^3)/12;
                  FIx(i)=(12*ELAS*INERx(i)*BET(i))/(GCOR*AREA(i)*L(i)^2);
                  FIy(i)=(12*ELAS*INERy(i)*BET(i))/(GCOR*AREA(i)*L(i)^2);
                  % GENERACIÓN DE LAS VARIABLES
                  FI(i)=(1/4)*FIx(i);
                  C=((4*ELAS*INER(i))/L(i))*((1+FI(i))/(1+4*FI(i)));
                  CP=C;
                  A=((2*ELAS*INER(i))/L(i))*((1-2*FI(i))/(1+4*FI(i)));
                  B=((6*ELAS*INER(i))/L(i)^2)*(1/(1+4*FI(i)));
                  BP=B;
                  T=(12*ELAS*INER(i))/L(i)^3*(1/(1+4*FI(i)));
                  R=ELAS*AREA(i)/L(i);
                  % MATRIZ DE RIGIDECES EN COORDENADAS LOCALES
                  K2=[R 0 0 -R 0 0;
                  0 T (B+C1(i)*T) 0 -T (B+C2(i)*T);
                  0 (B+C1(i)*T) (C+2*C1(i)*B+C1(i)^2*T) 0 -(B+C1(i)*T) (A+C1(i)*B+C2(i)*B+C1(i)*C2(i)*T);
                  -R 0 0 R 0 0;
                  0 -T -(B+C1(i)*T) 0 T -(B+C2(i)*T);
                  0 (B+C2(i)*T) (A+C1(i)*B+C2(i)*B+C1(i)*C2(i)*T) 0 -(B+C2(i)*T) (C+2*C2(i)*B+C2(i)^2*T)];
                  % MATRIZ DE TRANSPORTE A COORDENADAS GLOBALES
                  T2_3=zeros(6,6);
                  T2_3(1,1)=COSENO(i);T2_3(1,2)=SENO(i);
                  T2_3(2,1)=-SENO(i);T2_3(2,2)=COSENO(i);
                  T2_3(3,3)=1;
                  T2_3(4,4)=COSENO(i);T2_3(4,5)=SENO(i);
                  T2_3(5,4)=-SENO(i);T2_3(5,5)=COSENO(i);
                  T2_3(6,6)=1;
                  % ELEMENTO CUADRADO SÓLIDO DE SECCIÓN VARIABLE
              elseif DD(i)==2
                  Bi(i)=input('\n BASE INICIAL [m]:');
                  Hi(i)=Bi(i);
                  Bf(i)=input('\n BASE FINAL [m]:');
                  % ELEMENTOS UNICAMENTE NECESARIOS PARA EL CONTEO DURANTE EL USO
                  % DE ELEMENTOS DE CUALQUIERA DE LOS DISTINTOS TIPOS POSIBLES DE
                  % SECCIONES
                  Hf(i)=Bf(i);
                  B(i)=Bi(i);
                  H(i)=Hi(i);
                  Di(i)=Bi(i);
                  Df(i)=Di(i);
                  BET(i)=1.2;
                  INER(i)=(B(i)*H(i)^3)/12;
                  FI(i)=1;
                  AREA(i)=B(i)*H(i);
                  C1(i)=C1(i-1);
                  C2(i)=C2(i-1);
                  F(i)=F(i-1);
                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                  % ELEMENTOS DE LA MATRIZ DE FLEXIBILIDADES
                  f11=(L(i)/(ELAS*Hi(i)*Hi(i)))*(Hi(i)/(Hf(i)-Hi(i)))*(1-Hi(i)/Hf(i));
                  f22=(4*(L(i))^3/(ELAS*(Hi(i))^4))*(Hi(i)/Hf(i))^3+(6*L(i)/(5*GCOR*(Hi(i))^2))*(Hi(i)/(Hf(i)-Hi(i)))*(1-Hi(i)/Hf(i));
                  f33=f22;
                  f26=(4*(L(i))^2/(ELAS*(Hi(i))^4))*((Hi(i)/Hf(i))^3+(1/2)*(Hi(i)/(Hf(i)-Hi(i)))^2*(1+(Hi(i)/Hf(i))^2-2*(Hi(i)/Hf(i))));
                  f35=f26;
                  f66=(4*L(i)/(ELAS*(Hi(i))^4))*(Hi(i)/(Hf(i)-Hi(i)))*(1-(Hi(i)/Hf(i))^3);
                  f55=f66;
                  f44=(64*L(i)/(27*GCOR*(Hi(i))^4))*(Hi(i)/(Hf(i)-Hi(i)))*(1-(Hi(i)/Hf(i))^3);
                  FF44(i)=0;
                  % ELEMENTOS DE LA MATRIZ DE RIGIDECES
                  raz=1/f11;
                  Detx=f22*f66-f26*f26;
                  r11x=f22/Detx;
                  r12x=(f26*L(i)-f22)/Detx;
                  r22x=(f66*(L(i))^2-2*f26*L(i)+f22)/Detx;
                  raax=(r11x+r22x+2*r12x)/(L(i))^2;
                  rabx=(r11x+r12x)/L(i);
                  rbax=(r22x+r12x)/L(i);
                  % MATRIZ DE RIGIDECES EN COORDENADAS LOCALES
                  K11=[raz 0 0 ; 0 raax rabx ; 0 rabx r11x];
                  K12=[-raz 0 0 ; 0 -raax rbax ; 0 -rabx r12x];
                  K21=K12';
                  K22=[raz 0 0 ; 0 raax -rbax ; 0 -rbax r22x];
                  K2=[K11 K12 ;K21 K22];
                  % MATRIZ DE TRANSPORTE A COORDENADAS GLOBALES
                  T2_3=zeros(6,6);
                  T2_3(1,1)=COSENO(i);T2_3(1,2)=SENO(i);
                  T2_3(2,1)=-SENO(i);T2_3(2,2)=COSENO(i);
                  T2_3(3,3)=1;
                  T2_3(4,4)=COSENO(i);T2_3(4,5)=SENO(i);
                  T2_3(5,4)=-SENO(i);T2_3(5,5)=COSENO(i);
                  T2_3(6,6)=1;
                  % ELEMENTO RECTANGULAR SÓLIDO DE SECCIÓN VARIABLE
              elseif DD(i)==3
                  Bi(i)=input('\n BASE DE TODO EL ELEMENTO [m]:');
                  Hi(i)=input(' ALTURA INICIAL DEL ELEMENTO [m]:');
                  Bf(i)=Bi(i);
                  Hf(i)=input(' ALTURA FINAL DEL ELEMENTO [m]:');
                  Bb=Bi(i);
                  % ELEMENTOS UNICAMENTE NECESARIOS PARA EL CONTEO DURANTE EL USO
                  % DE ELEMENTOS DE CUALQUIERA DE LOS DISTINTOS TIPOS POSIBLES DE
                  % SECCIONES
                  B(i)=Bi(i);
                  H(i)=Hi(i);
                  Di(i)=Bi(i);
                  Df(i)=Di(i);
                  INER(i)=(B(i)*H(i)^3)/12;
                  FI(i)=1;
                  AREA(i)=B(i)*H(i);
                  BET(i)=1.2;
                  C1(i)=0;
                  C2(i)=0;
                  F(i)=L(i);
                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                  BET(i)=input('\n FACTOR DE FORMA (RECTANGULAR=1.2; CIRCULAR=10/9):');
                  fprintf('\n\t\t\t ‘ ; ......................................................................................................................... ');
                  fprintf('\n\t\t\t ‘ ; ......................................................................................................................... ');
                  fprintf('\n\n\t\t     LA FUNCIÓN F(x) INDICA LA VARIACIÓN CON LA LONGITUD DE LA FORMA DEL ELEMENTO,');
                  fprintf('\n\t\t\t                                                                    ');
                  fprintf('\n\n\t\t     SI LA FUNCIÓN ES LINEAL: Y=F(x)=');
                  fprintf('\n\t\t\t                                      (abs(Hi-Hf)/LL)*x                                       ');
                  fprintf('\n\t\t\t ‘ ; ......................................................................................................................... ');
                  fprintf('\n\t\t\t ‘ ; ......................................................................................................................... ');
                  fprintf('\n\t\t\t ‘ ;                                                                                                                           ');
                  fprintf('\n\t\t\t   ESTANDO EL EJE X SOBRE EL     Y|                                                                                            ');
                  fprintf('\n\t\t\t    CENTROIDE DEL MIEMBRO         |                                                                                            ');
                  fprintf('\n\t\t\t                                  |                                                                                            ');
                  fprintf('\n\t\t\t                                  |_ _ _ _ _                                                                                   ');
                  fprintf('\n\t\t\t                                           X                                                                                   ');
                  fprintf('\n\t\t\t                                                                                                                               ');
                  fprintf('\n\t\t\t ‘ ; ......................................................................................................................... ');
                  fprintf('\n\t\t\t ‘ ; ......................................................................................................................... ');
                  fprintf('\n\t\t\t                                                                                                                              ');
                  fprintf('\n\t\t\t                                                                                                                              ');
                  EE(i)=menu('ESCOJA UNA DE LAS SIGUIENTES OPCIONES:','FUNCIÓN LINEAL','OTRO TIPO DE VARIACIÓN');
                  if EE(i)==1
                      % ELEMENTOS DE LA MATRIZ DE FLEXIBILIDADES
                      f11=(L(i)/(ELAS*Bb*Hi(i)))*(Hi(i)/(Hf(i)-Hi(i)))*log(Hf(i)/Hi(i));
                      f22=(6*(L(i))^3/(ELAS*Bb*(Hi(i))^3))*((Hi(i)/Hf(i))^2-(Hi(i)/(Hf(i)-Hi(i)))^3*(Hf(i)/Hi(i)-Hi(i)/Hf(i)-2*log(Hf(i)/Hi(i))))+(6*L(i)/(5*GCOR*Bb*Hi(i)))*(Hi(i)/(Hf(i)-Hi(i)))*log(Hf(i)/Hi(i));
                      f33=(12*(L(i))^3/(ELAS*(Bb)^3*Hi(i)))*(3/2+(1/2)*(Hf(i)/Hi(i))^2-2*(Hf(i)/Hi(i))+log(Hf(i)/Hi(i)))+(6*L(i)/(5*GCOR*Bb*Hi(i)))*(Hi(i)/(Hf(i)-Hi(i)))*log(Hf(i)/Hi(i));
                      f26=(6*L(i)/(ELAS*Bb*(Hi(i))^3))*(Hi(i)/(Hf(i)-Hi(i)))^2*(1+(Hi(i)/Hf(i))^2-2*(Hi(i)/Hf(i)));
                      f35=(12*(L(i))^2/(ELAS*(Bb)^3*Hi(i)))*(Hi(i)/(Hf(i)-Hi(i)))^2*(Hf(i)/Hi(i)-log(Hf(i)/Hi(i))-1);
                      f66=(6*L(i)/(ELAS*Bb*(Hi(i))^3))*(Hi(i)/(Hf(i)-Hi(i)))*(1-(Hi(i)/Hf(i))^2);
                      f55=(12*L(i)/(ELAS*(Bb)^3*Hi(i)))*(Hi(i)/(Hf(i)-Hi(i)))*log(Hf(i)/Hi(i));
                      % DEFINICIÓN DE UNA FUNCIÓN DE x
                      Hi=Hi(i);
                      Hf=Hf(i);
                      LL=L(i);
                      if Hi< Hf
                          dpc= @(x)(1./(((Bb.^3*(Hi./2-(abs(Hi-Hf)./(2*LL))*x))./3)*(1-(192*Bb./(pi.^5*(Hi./2-(abs(Hi-Hf)./(2*LL))*x)))*(tanh(pi*(Hi./2-(abs(Hi-Hf)./(2*LL))*x)./(2*Bb)))))); 
                      elseif Hi > Hf
                          dpc= @(x)(1./(((Bb.^3*(Hi./2+(abs(Hi-Hf)./(2*LL))*x))./3)*(1-(192*Bb./(pi.^5*(Hi./2+(abs(Hi-Hf)./(2*LL))*x)))*(tanh(pi*(Hi./2+(abs(Hi-Hf)./(2*LL))*x)./(2*Bb))))));
                      end
                      f44=quad(@(x)dpc, 0, LL)
                      FF44(i)=f44;
                      % ELEMENTOS DE LA MATRIZ DE RIGIDECES
                      raz=1/f11;
                      Detx=f22*f66-f26*f26;
                      r11x=f22/Detx;
                      r12x=(f26*L(i)-f22)/Detx;
                      r22x=(f66*(L(i))^2-2*f26*L(i)+f22)/Detx;
                      raax=(r11x+r22x+2*r12x)/(L(i))^2;
                      rabx=(r11x+r12x)/L(i);
                      rbax=(r22x+r12x)/L(i);
                      % MATRIZ DE RIGIDECES EN COORDENADAS LOCALES
                      K11=[raz 0 0 ; 0 raax rabx ; 0 rabx r11x];
                      K12=[-raz 0 0 ; 0 -raax rbax ; 0 -rabx r12x];
                      K21=K12';
                      K22=[raz 0 0 ; 0 raax -rbax ; 0 -rbax r22x];
                      K2=[K11 K12 ;K21 K22];
                      % MATRIZ DE TRANSPORTE A COORDENADAS GLOBALES
                      T2_3=zeros(6,6);
                      T2_3(1,1)=COSENO(i);T2_3(1,2)=SENO(i);
                      T2_3(2,1)=-SENO(i);T2_3(2,2)=COSENO(i);
                      T2_3(3,3)=1;
                      T2_3(4,4)=COSENO(i);T2_3(4,5)=SENO(i);
                      T2_3(5,4)=-SENO(i);T2_3(5,5)=COSENO(i);
                      T2_3(6,6)=1;
                  else
                      % ELEMENTOS DE LA MATRIZ DE FLEXIBILIDADES
                      f11=(L(i)/(ELAS*Bb*Hi(i)))*(Hi(i)/(Hf(i)-Hi(i)))*log(Hf(i)/Hi(i));
                      f22=(6*(L(i))^3/(ELAS*Bb*(Hi(i))^3))*((Hi(i)/Hf(i))^2-(Hi(i)/(Hf(i)-Hi(i)))^3*(Hf(i)/Hi(i)-Hi(i)/Hf(i)-2*log(Hf(i)/Hi(i))))+(6*L(i)/(5*GCOR*Bb*Hi(i)))*(Hi(i)/(Hf(i)-Hi(i)))*log(Hf(i)/Hi(i));
                      f33=(12*(L(i))^3/(ELAS*(Bb)^3*Hi(i)))*(3/2+(1/2)*(Hf(i)/Hi(i))^2-2*(Hf(i)/Hi(i))+log(Hf(i)/Hi(i)))+(6*L(i)/(5*GCOR*Bb*Hi(i)))*(Hi(i)/(Hf(i)-Hi(i)))*log(Hf(i)/Hi(i));
                      f26=(6*L(i)/(ELAS*Bb*(Hi(i))^3))*(Hi(i)/(Hf(i)-Hi(i)))^2*(1+(Hi(i)/Hf(i))^2-2*(Hi(i)/Hf(i)));
                      f35=(12*(L(i))^2/(ELAS*(Bb)^3*Hi(i)))*(Hi(i)/(Hf(i)-Hi(i)))^2*(Hf(i)/Hi(i)-log(Hf(i)/Hi(i))-1);
                      f66=(6*L(i)/(ELAS*Bb*(Hi(i))^3))*(Hi(i)/(Hf(i)-Hi(i)))*(1-(Hi(i)/Hf(i))^2);
                      f55=(12*L(i)/(ELAS*(Bb)^3*Hi(i)))*(Hi(i)/(Hf(i)-Hi(i)))*log(Hf(i)/Hi(i));
                      % DEFINICIÓN DE UNA FUNCIÓN DE x
                      Hi=Hi(i);
                      Hf=Hf(i);
                      LL=L(i);
                      f44=input(' INTRODUZCA LA f44=');
                      FF44(i)=f44;
                      % ELEMENTOS DE LA MATRIZ DE RIGIDECES
                      raz=1/f11;
                      Detx=f22*f66-f26*f26;
                      r11x=f22/Detx;
                      r12x=(f26*L(i)-f22)/Detx;
                      r22x=(f66*(L(i))^2-2*f26*L(i)+f22)/Detx;
                      raax=(r11x+r22x+2*r12x)/(L(i))^2;
                      rabx=(r11x+r12x)/L(i);
                      rbax=(r22x+r12x)/L(i);
                      % MATRIZ DE RIGIDECES EN COORDENADAS LOCALES
                      K11=[raz 0 0 ; 0 raax rabx ; 0 rabx r11x];
                      K12=[-raz 0 0 ; 0 -raax rbax ; 0 -rabx r12x];
                      K21=K12';
                      K22=[raz 0 0 ; 0 raax -rbax ; 0 -rbax r22x];
                      K2=[K11 K12 ;K21 K22];
                      % MATRIZ DE TRANSPORTE A COORDENADAS GLOBALES
                      T2_3=zeros(6,6);
                      T2_3(1,1)=COSENO(i);T2_3(1,2)=SENO(i);
                      T2_3(2,1)=-SENO(i);T2_3(2,2)=COSENO(i);
                      T2_3(3,3)=1;
                      T2_3(4,4)=COSENO(i);T2_3(4,5)=SENO(i);
                      T2_3(5,4)=-SENO(i);T2_3(5,5)=COSENO(i);
                      T2_3(6,6)=1;
                  end
                  % ELEMENTO CIRCULAR SÓLIDO DE SECCIÓN CONSTATE
              elseif DD(i)==4
                  Di(i)=input('\n DIÁMETRO DEL ELEMENTO [m]:');
                  % ELEMENTOS UNICAMENTE NECESARIOS PARA EL CONTEO DURANTE EL USO
                  % DE ELEMENTOS DE CUALQUIERA DE LOS DISTINTOS TIPOS POSIBLES DE
                  % SECCIONES
                  Df(i)=Di(i);
                  B(i)=Di(i);
                  H(i)=Di(i);
                  Bi(i)=B(i);
                  Hi(i)=H(i);
                  Bf(i)=Bi(i);
                  Hf(i)=Hi(i);
                  FF44(i)=0;
                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                  BET(i)=input('\n VALOR DEL FACTOR DE FORMA (RECTANGULAR=1.2; CIRCULAR=10/9):');
                  AREA(i)=pi*Di(i)*Di(i)/4;
                  INER(i)=pi*(Di(i)/2)^4/4;
                  C1(i)=input('\n\n INGRESE EL VALOR DE C1 [m] (LONGITUD DEL SEGMENTO DE RIGIDEZ INFINITA EN EL NODO INICIAL):');
                  C2(i)=input('\n\n INGRESE EL VALOR DE C2 [m] (LONGITUD DEL SEGMENTO DE RIGIDEZ INFINITA EN EL NODO FINAL):');
                  F(i)=L(i)-C1(i)-C2(i)
                  FI(i)=(3*ELAS*INER(i)*BET(i))/(GCOR*AREA(i)*L(i)^2);
                  C=((4*ELAS*INER(i))/L(i))*((1+FI(i))/(1+4*FI(i)));
                  CP=C;
                  A=((2*ELAS*INER(i))/L(i))*((1-2*FI(i))/(1+4*FI(i)));
                  B=((6*ELAS*INER(i))/L(i)^2)*(1/(1+4*FI(i)));
                  BP=B;
                  T=(12*ELAS*INER(i))/L(i)^3*(1/(1+4*FI(i)));
                  R=ELAS*AREA(i)/L(i);
                  K2=[R 0 0 -R 0 0;
                  0 T (B+C1(i)*T) 0 -T (B+C2(i)*T);
                  0 (B+C1(i)*T) (C+2*C1(i)*B+C1(i)^2*T) 0 -(B+C1(i)*T) (A+C1(i)*B+C2(i)*B+C1(i)*C2(i)*T);
                  -R 0 0 R 0 0;
                  0 -T -(B+C1(i)*T) 0 T -(B+C2(i)*T);
                  0 (B+C2(i)*T) (A+C1(i)*B+C2(i)*B+C1(i)*C2(i)*T) 0 -(B+C2(i)*T) (C+2*C2(i)*B+C2(i)^2*T)];
                  % MATRIZ DE TRANSPORTE A COORDENADAS GLOBALES
                  T2_3=zeros(6,6);
                  T2_3(1,1)=COSENO(i);T2_3(1,2)=SENO(i);
                  T2_3(2,1)=-SENO(i);T2_3(2,2)=COSENO(i);
                  T2_3(3,3)=1;
                  T2_3(4,4)=COSENO(i);T2_3(4,5)=SENO(i);
                  T2_3(5,4)=-SENO(i);T2_3(5,5)=COSENO(i);
                  T2_3(6,6)=1;
                  % ELEMENTO CIRCULAR SÓLIDO DE SECCIÓN VARIABLE
              elseif DD(i)==5
                  Di(i)=input('\n DIÁMETRO INICIAL DEL ELEMENTO [m]:');
                  Df(i)=input('\n DIÁMETRO FINAL DEL ELEMENTO [m]:');
                  BET(i)=input('\n VALOR DEL FACTOR DE FORMA (RECTANGULAR=1.2; CIRCULAR=10/9):');
                  % ELEMENTOS UNICAMENTE NECESARIOS PARA EL CONTEO DURANTE EL USO
                  % DE ELEMENTOS DE CUALQUIERA DE LOS DISTINTOS TIPOS POSIBLES DE
                  % SECCIONES
                  B(i)=Di(i);
                  H(i)=Df(i);
                  Bi(i)=B(i);
                  Hi(i)=H(i);
                  Bf(i)=Bi(i);
                  Hf(i)=Hi(i);
                  C1(i)=0;
                  C2(i)=0;
                  F(i)=L(i);
                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                  AREA(i)=pi*Di(i)*Di(i)/4;
                  INER(i)=pi*(Di(i)/2)^4/4;
                  FIx(i)=(12*ELAS*INER(i)*BET(i))/(GCOR*AREA(i)*L(i)^2)
                  FIy(i)=FIx(i);
                  FI(i)=(1/4)*FIx(i);
                  % ELEMENTOS DE LA MATRIZ DE FLEXIBILIDADES
                  f11=(4*L(i)/(pi*ELAS*(Di(i))^2))*(Di(i)/(Df(i)-Di(i)))*(1-Di(i)/Df(i));
                  f22=(64*(L(i))^3/(3*pi*ELAS*(Di(i))^4))*(Di(i)/Df(i))^3+(40*L(i)/(9*pi*GCOR*(Di(i))^2))*(Di(i)/(Df(i)-Di(i)))*(1-(Di(i)/Df(i)));
                  f33=f22;
                  f26=(64*(L(i))^2/(3*pi*ELAS*(Di(i))^4))*((Di(i)/Df(i))^3+(1/2)*(Di(i)/(Df(i)-Di(i)))^2*(1+(Di(i)/Df(i))^2-2*(Di(i)/Df(i))));
                  f35=f26;
                  f66=(64*L(i)/(3*pi*ELAS*(Di(i))^4))*(Di(i)/(Df(i)-Di(i)))*(1-(Di(i)/Df(i))^3);
                  f55=f66;
                  f44=(32*L(i)/(3*pi*GCOR*Di(i)^4))*(Di(i)/(Df(i)-Di(i)))*(1-(Di(i)/Df(i))^3);
                  FF44(i)=f44;
                  % ELEMENTOS DE LA MATRIZ DE RIGIDECES
                  raz=1/f11;
                  Detx=f22*f66-f26*f26;
                  r11x=f22/Detx;
                  r12x=(f26*L(i)-f22)/Detx;
                  r22x=(f66*(L(i))^2-2*f26*L(i)+f22)/Detx;
                  raax=(r11x+r22x+2*r12x)/(L(i))^2;
                  rabx=(r11x+r12x)/L(i);
                  rbax=(r22x+r12x)/L(i);
                  % MATRIZ DE RIGIDECES EN COORDENADAS LOCALES
                  K11=[raz 0 0 ; 0 raax rabx ; 0 rabx r11x];
                  K12=[-raz 0 0 ; 0 -raax rbax ; 0 -rabx r12x];
                  K21=K12';
                  K22=[raz 0 0 ; 0 raax -rbax ; 0 -rbax r22x];
                  K2=[K11 K12 ;K21 K22];
                  % MATRIZ DE TRANSPORTE A COORDENADAS GLOBALES
                  T2_3=zeros(6,6);
                  T2_3(1,1)=COSENO(i);T2_3(1,2)=SENO(i);
                  T2_3(2,1)=-SENO(i);T2_3(2,2)=COSENO(i);
                  T2_3(3,3)=1;
                  T2_3(4,4)=COSENO(i);T2_3(4,5)=SENO(i);
                  T2_3(5,4)=-SENO(i);T2_3(5,5)=COSENO(i);
                  T2_3(6,6)=1;
              end
           end   
       end       
      %------------------ ENSAMBLAJE DE LA MATRIZ--------------------------
      %--------------------------------------------------------------------      %
      K3=T2_3'*K2*T2_3; % MATRIZ DE RIGIDEZ EN COORDENADAS LOCALES
      MAT_AUX=K3;
      for j=1:6
         if VC(i,j)==0
            MAT_AUX(j,:)=0;
            MAT_AUX(:,j)=0;
         end
      end
      K=zeros(ngl,ngl);
      for cont1=1:6
         for cont2=1:6
            if MAT_AUX(cont1,cont2)~=0
               UNO=VC(i,cont1);
               DOS=VC(i,cont2);
               TRES=MAT_AUX(cont1,cont2);
               K(UNO,DOS)=K(UNO,DOS)+TRES;
            end
         end
      end
      if i==1
         AUX=K;
      elseif i~=1
         AUX=AUX+K;
      else,end
   end
   fprintf('\n MATRIZ DE RIGIDEZ DE PRIMER ORDEN DEL SISTEMA:');
   MAT_RIGIDEZ=AUX   
   %-------------- ARREGLO Q_TOTAL ( VECTOR DE CARGAS ) -------------------
   Q=zeros(ngl,1);
   V1=menu('¿EXISTEN CARGAS EN JUNTAS?','SÍ','NO');
   if V1==1
      njc=input('\n INGRESE EL NÚMERO DE JUNTAS CARGADAS:');
      for i=1:njc
         NC=input('\n NUDO CARGADO:');
         Q1(1)=input(' FUERZA HORIZONTAL=');
         Q1(2)=input(' FUERZA VERTICAL=');
         Q1(3)=input(' MOMENTO=');
         VCJ(i,:)=CG(NC,:);
         for m=1:3
            n=VCJ(i,m)
            if Q1(m)~=0
               Q(n)=Q1(m);
            end
         end
      end
      Q_CJ=Q
   else
      Q_CJ=Q
   end
   Q=zeros(ngl,1);
   Q2_ALMAC=zeros(mbr,6);
   V3=menu('¿EXISTEN CARGAS EN LOS MIEMBROS?','SÍ','NO');
   if V3==1
      nmc=input('\n INGRESE EL NÚMERO DE MIEMBROS CARGADOS:');
      for ll=1:nmc
         clear Q2
         fprintf('\n\t EL PROGRAMA UNICAMENTE CONSIDERA CAMBIOS POR CORTANTE EN LAS CARGAS');
         fprintf('\n\t PUNTUALES LO CUAL NO ES INCORRECTO DEBIDO A QUE EN ESTE ESTADO DE CARGA');
         fprintf('\n\t ES DONDE EXISTE MAYOR IMPORTANCIA A LAS DEFORMACIONES DE CORTANTE');
         fprintf('\n\t ');
         fprintf('\n\t ');
         fprintf('\n\t ');
         MC=input('\n INGRESE EN ORDEN NUMÉRICO CADA MIEMBRO CARGADO:');
         fprintf('\n\t ');
         fprintf('\n\t ESCOJA UNA DE LAS SIGUIENTES OPCIONES DE CARGAS:');
         fprintf('\n\t ');
         fprintf('\n\n (1) CARGA UNIFORMEMENTE DISTRIBUIDA EN TODO EL MIEMBRO');
         fprintf('\n (2) CARGA TRIANGULAR CON ALTURA MÁXIMA EN EL CENTRO');
         fprintf('\n (3) CARGA TRIANGULAR CON CARGA CRECIENTE DESDE EL NUDO INICIAL');
         fprintf('\n (4) CARGA TRAPEZOIDAL');
         fprintf('\n (5) CARGA PUNTUAL EN EL CENTRO DEL MIEMBRO');
         fprintf('\n (6) CARGA PUNTUAL A DIFERENTE DISTANCIA');
         fprintf('\n (7) NINGUNA DE LAS ANTERIORES');
         V4=input('\n\n\n OPCIÓN:');
         % CARGA UNIFORMEMENTE DISTRIBUIDA EN TODO EL MIEMBRO
         if V4==1
            CAR=input('\n CARGA (T/m):');
            Q2(1)=0;
            Q2(2)=CAR*L(MC)/2;
            Q2(3)=CAR*L(MC)^2/12;
            Q2(4)=0;
            Q2(5)=Q2(2);
            Q2(6)=-Q2(3);
            % CARGA TRIANGULARMENTE DISTRIBUIDA CON ALTURA MÁXIMA AL CENTRO
            % LONGITUDINAL DEL MIEMBRO
         elseif V4==2
            CAR=input('\n CARGA MÁXIMA (T/m):');
            Q2(1)=0;
            Q2(2)=CAR*L(MC)/4;
            Q2(3)=5*CAR*L(MC)^2/96;
            Q2(4)=0;
            Q2(5)=Q2(2);
            Q2(6)=-Q2(3); 
            % CARGA TRIANGULAR CRECIENTE DESDE EL NUDO INICIAL
         elseif V4==3
            CAR=input('\n CARGA MÁXIMA (T/m):');
            Q2(1)=0;
            Q2(2)=3*CAR*L(MC)/20;
            Q2(3)=CAR*L(MC)^2/30;
            Q2(4)=0;
            Q2(5)=7*CAR*L(MC)/20;
            Q2(6)=-CAR*L(MC)^2/20;
            % CARGA TRAPEZOIDALMENTE DISTRIBUIDA
         elseif V4==4
            CAR=input('\n CARGA MÁXIMA (T/m):');
            a=input('\n DISTANCIA a (m):');
            Q2(1)=0;
            Q2(2)=(CAR*L(MC)/2)*(1-(a/L(MC)));
            Q2(3)=(CAR*L(MC)^2/12)*(1-2*(a/L(MC))^2+(a/L(MC))^3);
            Q2(4)=0;
            Q2(5)=Q2(2);
            Q2(6)=-Q2(3); 
            % CARGA PUNTUAL EN EL CENTRO DEL MIEMBRO
         elseif V4==5
             CAR=input('\n CARGA (T):');
             Q2(1)=0;
             Q2(2)=CAR/2;
             Q2(3)=CAR*L(MC)/8;
             Q2(4)=0;
             Q2(5)=Q2(2);
             Q2(6)=-Q2(3);  
            % CARGA PUNTUAL EN CUALQUIER PUNTO DEL EJE LONGITUDINAL DEL
            % MIEMBRO
         elseif V4==6
             CAR=input('\n CARGA (T):');
             a=input('\n DISTANCIA a (m):');
             b=input('\n DISTANCIA b (m):');
             Q2(1)=0;
             Q2(2)=(CAR*b^2/L(MC)^2)*(3-(2*b)/(a+b));
             Q2(3)=(CAR*a*b^2/L(MC)^2);
             Q2(4)=0;
             Q2(5)=(CAR*a^2/L(MC)^2)*(3-(2*a)/(a+b));
             Q2(6)=-(CAR*a^2*b/L(MC)^2);           
            % CARGA POR INGRESAR; RECOMENDADO PARA MIEMBROS
            % INCLINADOS U OTRAS CONDICIONES DE CARGA
         elseif V4==7
            Q2(1)=input('\n FUERZA AXIAL DEL NUDO INICIAL (T):');
            Q2(2)=input('\n FUERZA CORTANTE DEL NUDO INICIAL (T):');
            Q2(3)=input('\n MOMENTO DEL NUDO INICIAL (T*m):');
            Q2(4)=input('\n FUERZA AXIAL DEL NUDO FINAL (T):');
            Q2(5)=input('\n FUERZA CORTANTE DEL NUDO FINAL (T):');
            Q2(6)=input('\n MOMENTO DEL NUDO FINAL (T*m):');
         end
         Q2=Q2'
         for mm=1:6
            Q2_ALMAC(MC,mm)=Q2(mm)';
         end
         % GENERACIÓN DE LAS CARGAS EN COORDENADAS GLOBALES
         T2_3=zeros(6,6);
         T2_3(1,1)=COSENO(MC);T2_3(1,2)=SENO(MC);
         T2_3(2,1)=-SENO(MC);T2_3(2,2)=COSENO(MC);
         T2_3(3,3)=1;
         T2_3(4,4)=COSENO(MC);T2_3(4,5)=SENO(MC);
         T2_3(5,4)=-SENO(MC);T2_3(5,5)=COSENO(MC);
         T2_3(6,6)=1;
         clear Q3
         Q3=-T2_3'*Q2;
         for g=1:6
            h=VC(MC,g);
            if h~=0
               Q(h)=Q3(g)+Q(h);
            end
         end
      end
      Q_CM=Q
      Q3
   else
      Q_CM=Q
   end
   fprintf('\n VECTOR DE CARGAS EN COORDENADAS GLOBALES:');
   Q_TOTAL=Q_CJ+Q_CM
   
   %------------ Q_DES ----------------------------------------------------
   fprintf('\n DESPLAZAMIENTOS GENERALIZADOS DE PRIMER ORDEN EN COORDENADAS GLOBALES:');
   Q_des=MAT_RIGIDEZ\Q_TOTAL
   
   %------------- ARREGLO DE DESPLAZAMIENTOS Y P_FINAL ------------------------------------
   for i=1:mbr
      fprintf('\n MIEMBRO %d:',i);
      fprintf('\n DEFORMACIONES EN MIEMBRO P1:');
      clear P1 P2 P Q2_AUX
      for j=1:6
         if VC(i,j)~=0
            P1(j)=Q_des(VC(i,j));
         else % ELIMINACIÓN DE LOS GRADOS SIN LIBERTAD EN LA ESTRUCTURA
            P1(j)=0;
         end
      end
      P1 % VECTOR DE DESPLAZAMIENTOS d
      if DD(i)==1
          %
          % MATRIZ DE RIGIDECES EN COORDENADAS LOCALES
          % PARA CADA CLASE DE ELEMENTO
          %
          % PARA LA GENERACIÓN DE VARIABLES SE 
          % INCLUYEN LOS SIGUIENTES CÁLCULOS
          C=((4*ELAS*INER(i))/L(i))*((1+FI(i))/(1+4*FI(i)));
          CP=C;
          A=((2*ELAS*INER(i))/L(i))*((1-2*FI(i))/(1+4*FI(i)));
          B=((6*ELAS*INER(i))/L(i)^2)*(1/(1+4*FI(i)));
          BP=B;
          T=(12*ELAS*INER(i))/L(i)^3*(1/(1+4*FI(i)));
          R=ELAS*AREA(i)/L(i);
          K2=[R 0 0 -R 0 0;
          0 T (B+C1(i)*T) 0 -T (B+C2(i)*T);
          0 (B+C1(i)*T) (C+2*C1(i)*B+C1(i)^2*T) 0 -(B+C1(i)*T) (A+C1(i)*B+C2(i)*B+C1(i)*C2(i)*T);
          -R 0 0 R 0 0;
          0 -T -(B+C1(i)*T) 0 T -(B+C2(i)*T);
          0 (B+C2(i)*T) (A+C1(i)*B+C2(i)*B+C1(i)*C2(i)*T) 0 -(B+C2(i)*T) (C+2*C2(i)*B+C2(i)^2*T)];
          % MATRIZ DE TRANSPORTE A COORDENADAS GLOBALES
          T2_3=zeros(6,6);
          T2_3(1,1)=COSENO(i);T2_3(1,2)=SENO(i);
          T2_3(2,1)=-SENO(i);T2_3(2,2)=COSENO(i);
          T2_3(3,3)=1;
          T2_3(4,4)=COSENO(i);T2_3(4,5)=SENO(i);
          T2_3(5,4)=-SENO(i);T2_3(5,5)=COSENO(i);
          T2_3(6,6)=1;
          % ELEMENTO CUADRADO SÓLIDO DE SECCIÓN VARIABLE
      elseif DD(i)==2
          % ELEMENTOS DE LA MATRIZ DE FLEXIBILIDADES
          f11=(L(i)/(ELAS*Hi(i)*Hi(i)))*(Hi(i)/(Hf(i)-Hi(i)))*(1-Hi(i)/Hf(i));
          f22=(4*(L(i))^3/(ELAS*(Hi(i))^4))*(Hi(i)/Hf(i))^3+(6*L(i)/(5*GCOR*(Hi(i))^2))*(Hi(i)/(Hf(i)-Hi(i)))*(1-Hi(i)/Hf(i));
          f33=f22;
          f26=(4*(L(i))^2/(ELAS*(Hi(i))^4))*((Hi(i)/Hf(i))^3+(1/2)*(Hi(i)/(Hf(i)-Hi(i)))^2*(1+(Hi(i)/Hf(i))^2-2*(Hi(i)/Hf(i))));
          f35=f26;
          f66=(4*L(i)/(ELAS*(Hi(i))^4))*(Hi(i)/(Hf(i)-Hi(i)))*(1-(Hi(i)/Hf(i))^3);
          f55=f66;
          f44=(64*L(i)/(27*GCOR*(Hi(i))^4))*(Hi(i)/(Hf(i)-Hi(i)))*(1-(Hi(i)/Hf(i))^3);
          % ELEMENTOS DE LA MATRIZ DE RIGIDECES
          raz=1/f11;
          Detx=f22*f66-f26*f26;
          r11x=f22/Detx;
          r12x=(f26*L(i)-f22)/Detx;
          r22x=(f66*(L(i))^2-2*f26*L(i)+f22)/Detx;
          raax=(r11x+r22x+2*r12x)/(L(i))^2;
          rabx=(r11x+r12x)/L(i);
          rbax=(r22x+r12x)/L(i);
          % MATRIZ DE RIGIDECES EN COORDENADAS LOCALES
          K11=[raz 0 0 ; 0 raax rabx ; 0 rabx r11x];
          K12=[-raz 0 0 ; 0 -raax rbax ; 0 -rabx r12x];
          K21=K12';
          K22=[raz 0 0 ; 0 raax -rbax ; 0 -rbax r22x];
          K2=[K11 K12 ;K21 K22];
          % MATRIZ DE TRANSPORTE A COORDENADAS GLOBALES
          T2_3=zeros(6,6);
          T2_3(1,1)=COSENO(i);T2_3(1,2)=SENO(i);
          T2_3(2,1)=-SENO(i);T2_3(2,2)=COSENO(i);
          T2_3(3,3)=1;
          T2_3(4,4)=COSENO(i);T2_3(4,5)=SENO(i);
          T2_3(5,4)=-SENO(i);T2_3(5,5)=COSENO(i);
          T2_3(6,6)=1;
          % ELEMENTO RECTANGULAR SÓLIDO DE SECCIÓN VARIABLE
      elseif DD(i)==3
          % ELEMENTOS DE LA MATRIZ DE FLEXIBILIDADES
          f11=(L(i)/(ELAS*Bb*Hi(i)))*(Hi(i)/(Hf(i)-Hi(i)))*log(Hf(i)/Hi(i));
          f22=(6*(L(i))^3/(ELAS*Bb*(Hi(i))^3))*((Hi(i)/Hf(i))^2-(Hi(i)/(Hf(i)-Hi(i)))^3*(Hf(i)/Hi(i)-Hi(i)/Hf(i)-2*log(Hf(i)/Hi(i))))+(6*L(i)/(5*GCOR*Bb*Hi(i)))*(Hi(i)/(Hf(i)-Hi(i)))*log(Hf(i)/Hi(i));
          f33=(12*(L(i))^3/(ELAS*(Bb)^3*Hi(i)))*(3/2+(1/2)*(Hf(i)/Hi(i))^2-2*(Hf(i)/Hi(i))+log(Hf(i)/Hi(i)))+(6*L(i)/(5*GCOR*Bb*Hi(i)))*(Hi(i)/(Hf(i)-Hi(i)))*log(Hf(i)/Hi(i));
          f26=(6*L(i)/(ELAS*Bb*(Hi(i))^3))*(Hi(i)/(Hf(i)-Hi(i)))^2*(1+(Hi(i)/Hf(i))^2-2*(Hi(i)/Hf(i)));
          f35=(12*(L(i))^2/(ELAS*(Bb)^3*Hi(i)))*(Hi(i)/(Hf(i)-Hi(i)))^2*(Hf(i)/Hi(i)-log(Hf(i)/Hi(i))-1);
          f66=(6*L(i)/(ELAS*Bb*(Hi(i))^3))*(Hi(i)/(Hf(i)-Hi(i)))*(1-(Hi(i)/Hf(i))^2);
          f55=(12*L(i)/(ELAS*(Bb)^3*Hi(i)))*(Hi(i)/(Hf(i)-Hi(i)))*log(Hf(i)/Hi(i));
          f44=FF44(i);
          % ELEMENTOS DE LA MATRIZ DE RIGIDECES
          raz=1/f11;
          Detx=f22*f66-f26*f26;
          r11x=f22/Detx;
          r12x=(f26*L(i)-f22)/Detx;
          r22x=(f66*(L(i))^2-2*f26*L(i)+f22)/Detx;
          raax=(r11x+r22x+2*r12x)/(L(i))^2;
          rabx=(r11x+r12x)/L(i);
          rbax=(r22x+r12x)/L(i);
          % MATRIZ DE RIGIDECES EN COORDENADAS LOCALES
          K11=[raz 0 0 ; 0 raax rabx ; 0 rabx r11x];
          K12=[-raz 0 0 ; 0 -raax rbax ; 0 -rabx r12x];
          K21=K12';
          K22=[raz 0 0 ; 0 raax -rbax ; 0 -rbax r22x];
          K2=[K11 K12 ;K21 K22];
          % MATRIZ DE TRANSPORTE A COORDENADAS GLOBALES
          T2_3=zeros(6,6);
          T2_3(1,1)=COSENO(i);T2_3(1,2)=SENO(i);
          T2_3(2,1)=-SENO(i);T2_3(2,2)=COSENO(i);
          T2_3(3,3)=1;
          T2_3(4,4)=COSENO(i);T2_3(4,5)=SENO(i);
          T2_3(5,4)=-SENO(i);T2_3(5,5)=COSENO(i);
          T2_3(6,6)=1;
          % ELEMENTO CIRCULAR SÓLIDO DE SECCIÓN CONSTATE
      elseif DD(i)==4
          C=((4*ELAS*INER(i))/L(i))*((1+FI(i))/(1+4*FI(i)));
          CP=C;
          A=((2*ELAS*INER(i))/L(i))*((1-2*FI(i))/(1+4*FI(i)));
          B=((6*ELAS*INER(i))/L(i)^2)*(1/(1+4*FI(i)));
          BP=B;
          T=(12*ELAS*INER(i))/L(i)^3*(1/(1+4*FI(i)));
          R=ELAS*AREA(i)/L(i);
          K2=[R 0 0 -R 0 0;
          0 T (B+C1(i)*T) 0 -T (B+C2(i)*T);
          0 (B+C1(i)*T) (C+2*C1(i)*B+C1(i)^2*T) 0 -(B+C1(i)*T) (A+C1(i)*B+C2(i)*B+C1(i)*C2(i)*T);
          -R 0 0 R 0 0;
          0 -T -(B+C1(i)*T) 0 T -(B+C2(i)*T);
          0 (B+C2(i)*T) (A+C1(i)*B+C2(i)*B+C1(i)*C2(i)*T) 0 -(B+C2(i)*T) (C+2*C2(i)*B+C2(i)^2*T)];
          % MATRIZ DE TRANSPORTE A COORDENADAS GLOBALES
          T2_3=zeros(6,6);
          T2_3(1,1)=COSENO(i);T2_3(1,2)=SENO(i);
          T2_3(2,1)=-SENO(i);T2_3(2,2)=COSENO(i);
          T2_3(3,3)=1;
          T2_3(4,4)=COSENO(i);T2_3(4,5)=SENO(i);
          T2_3(5,4)=-SENO(i);T2_3(5,5)=COSENO(i);
          T2_3(6,6)=1;
          % ELEMENTO CIRCULAR SÓLIDO DE SECCIÓN VARIABLE
      elseif DD(i)==5
          % ELEMENTOS DE LA MATRIZ DE FLEXIBILIDADES
          f11=(4*L(i)/(pi*ELAS*(Di(i))^2))*(Di(i)/(Df(i)-Di(i)))*(1-Di(i)/Df(i));
          f22=(64*(L(i))^3/(3*pi*ELAS*(Di(i))^4))*(Di(i)/Df(i))^3+(40*L(i)/(9*pi*GCOR*(Di(i))^2))*(Di(i)/(Df(i)-Di(i)))*(1-(Di(i)/Df(i)));
          f33=f22;
          f26=(64*(L(i))^2/(3*pi*ELAS*(Di(i))^4))*((Di(i)/Df(i))^3+(1/2)*(Di(i)/(Df(i)-Di(i)))^2*(1+(Di(i)/Df(i))^2-2*(Di(i)/Df(i))));
          f35=f26;
          f66=(64*L(i)/(3*pi*ELAS*(Di(i))^4))*(Di(i)/(Df(i)-Di(i)))*(1-(Di(i)/Df(i))^3);
          f55=f66;
          f44=(32*L(i)/(3*pi*GCOR*Di(i)^4))*(Di(i)/(Df(i)-Di(i)))*(1-(Di(i)/Df(i))^3);
          % ELEMENTOS DE LA MATRIZ DE RIGIDECES
          raz=1/f11;
          Detx=f22*f66-f26*f26;
          r11x=f22/Detx;
          r12x=(f26*L(i)-f22)/Detx;
          r22x=(f66*(L(i))^2-2*f26*L(i)+f22)/Detx;
          raax=(r11x+r22x+2*r12x)/(L(i))^2;
          rabx=(r11x+r12x)/L(i);
          rbax=(r22x+r12x)/L(i);
          % MATRIZ DE RIGIDECES EN COORDENADAS LOCALES
          K11=[raz 0 0 ; 0 raax rabx ; 0 rabx r11x];
          K12=[-raz 0 0 ; 0 -raax rbax ; 0 -rabx r12x];
          K21=K12';
          K22=[raz 0 0 ; 0 raax -rbax ; 0 -rbax r22x];
          K2=[K11 K12 ;K21 K22];
          % MATRIZ DE TRANSPORTE A COORDENADAS GLOBALES
          T2_3=zeros(6,6);
          T2_3(1,1)=COSENO(i);T2_3(1,2)=SENO(i);
          T2_3(2,1)=-SENO(i);T2_3(2,2)=COSENO(i);
          T2_3(3,3)=1;
          T2_3(4,4)=COSENO(i);T2_3(4,5)=SENO(i);
          T2_3(5,4)=-SENO(i);T2_3(5,5)=COSENO(i);
          T2_3(6,6)=1;
      end
      % LO ANTERIOR SE REALIZA PARA LA GENERACIÓN NUEVAMENTE DE LA MATRIZ
      % DE RIGIDEZ
      fprintf('\n RESULTADOS FINALES DE PRIMER ORDEN EN EL MIEMBRO %d:',i);
      K2; % K EN COORDENADAS LOCALES
      K3=T2_3'*K2*T2_3; % K EN COORDENADAS GLOBALES
      P2=K3*P1';% K*d=q CARGAS EN COORDENADAS GLOBALES
      P=T2_3*P2; % R*q=q' (P.COMPLEMENTARIO) % CARGA EN COORDENADAS LOCALES 
      for j=1:6
         Q2_AUX(j)=Q2_ALMAC(i,j); % VALORES CORRESPONDIENTES A q DE EMPOTRAMIENTO 
         %                          EN COORDENADAS GLOBALES
      end
      Q2_AUX=Q2_AUX'; % q_emp
      fprintf('\n ACCIONES FINALES DEL MIEMBRO EN COORDENADAS LOCALES ( P.P.+P.C. ):');
      P_FINAL=Q2_AUX+P % ELEMENTOS MECÁNICOS FINALES: q_emp+q'   
   %   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % A PARTIR DE ESTE MOMENTO SE REALIZA EL SEGUNDO ANÁLISIS DE LA
   % ESTRUCTURA
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %
   %-------------- CÁLCULO DE LOS COEFICIENTES DE RIGIDEZ DE --------------
   %------------------------ LA MATRIZ GEOMETRICA -------------------------
   F1(i)=P_FINAL(1); % NORMAL 1
   F2(i)=P_FINAL(2); % CORTANTE 1
   F3(i)=P_FINAL(3); % MOMENTO 1
   F4(i)=P_FINAL(4); % NORMAL 2
   F5(i)=P_FINAL(5); % CORTANTE 2
   F6(i)=P_FINAL(6);% MOMENTO 2
   %---------------- SE CALCULAN LOS COEFICIENTES DE RIGIDEZ --------------
   %
   % Guide to STABILITY DESIGN CRITERIA for METAL STRUCTURES, Fifth Edition
   %Theodore V. Galambos
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %
   % Análisis de estructuras con métodos matriciales, Arturo Tena Colunga
   %
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   rg1x=(6*F4(i))/(5*L(i));
   rg2x=F4(i)/10;
   rg3x=(2*F4(i)*L(i))/15;
   rg4x=(F4(i)*L(i))/30;
   rv1=F2(i)/L(i);
   %-------------- SE GENERA LA MATRIZ DE RIGIDEZ GEOMÉTRICA---------------
   KG=[0 rv1 0 0 -rv1 0 ; 
       rv1 rg1x rg2x -rv1 -rg1x rg2x ;
       0 rg2x rg3x 0 -rg2x -rg4x ;
       0 -rv1 0 0 rv1 0 ;
       -rv1 -rg1x -rg2x rv1 rg1x -rg2x ;
       0 rg2x -rg4x 0 -rg2x rg3x];
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %------------------ ENSAMBLAJE DE LA MATRIZ--------------------------
      %------------ CONSIDERANDO EFECTOS DE SEGUNDO ORDEN------------------
      if DD(i)==1
          %
          % MATRIZ DE RIGIDECES EN COORDENADAS LOCALES
          % PARA CADA CLASE DE ELEMENTO
          % PARA LA GENERACIÓN DE VARIABLES SE 
          % INCLUYEN LOS SIGUIENTES CÁLCULOS
          C=((4*ELAS*INER(i))/L(i))*((1+FI(i))/(1+4*FI(i)));
          CP=C;
          A=((2*ELAS*INER(i))/L(i))*((1-2*FI(i))/(1+4*FI(i)));
          B=((6*ELAS*INER(i))/L(i)^2)*(1/(1+4*FI(i)));
          BP=B;
          T=(12*ELAS*INER(i))/L(i)^3*(1/(1+4*FI(i)));
          R=ELAS*AREA(i)/L(i);
          KE=[R 0 0 -R 0 0;
          0 T (B+C1(i)*T) 0 -T (B+C2(i)*T);
          0 (B+C1(i)*T) (C+2*C1(i)*B+C1(i)^2*T) 0 -(B+C1(i)*T) (A+C1(i)*B+C2(i)*B+C1(i)*C2(i)*T);
          -R 0 0 R 0 0;
          0 -T -(B+C1(i)*T) 0 T -(B+C2(i)*T);
          0 (B+C2(i)*T) (A+C1(i)*B+C2(i)*B+C1(i)*C2(i)*T) 0 -(B+C2(i)*T) (C+2*C2(i)*B+C2(i)^2*T)];
          % MATRIZ DE TRANSPORTE A COORDENADAS GLOBALES
          T2_3=zeros(6,6);
          T2_3(1,1)=COSENO(i);T2_3(1,2)=SENO(i);
          T2_3(2,1)=-SENO(i);T2_3(2,2)=COSENO(i);
          T2_3(3,3)=1;
          T2_3(4,4)=COSENO(i);T2_3(4,5)=SENO(i);
          T2_3(5,4)=-SENO(i);T2_3(5,5)=COSENO(i);
          T2_3(6,6)=1;
          % ELEMENTO CUADRADO SÓLIDO DE SECCIÓN VARIABLE
      elseif DD(i)==2
          % ELEMENTOS DE LA MATRIZ DE FLEXIBILIDADES
          f11=(L(i)/(ELAS*Hi(i)*Hi(i)))*(Hi(i)/(Hf(i)-Hi(i)))*(1-Hi(i)/Hf(i));
          f22=(4*(L(i))^3/(ELAS*(Hi(i))^4))*(Hi(i)/Hf(i))^3+(6*L(i)/(5*GCOR*(Hi(i))^2))*(Hi(i)/(Hf(i)-Hi(i)))*(1-Hi(i)/Hf(i));
          f33=f22;
          f26=(4*(L(i))^2/(ELAS*(Hi(i))^4))*((Hi(i)/Hf(i))^3+(1/2)*(Hi(i)/(Hf(i)-Hi(i)))^2*(1+(Hi(i)/Hf(i))^2-2*(Hi(i)/Hf(i))));
          f35=f26;
          f66=(4*L(i)/(ELAS*(Hi(i))^4))*(Hi(i)/(Hf(i)-Hi(i)))*(1-(Hi(i)/Hf(i))^3);
          f55=f66;
          f44=(64*L(i)/(27*GCOR*(Hi(i))^4))*(Hi(i)/(Hf(i)-Hi(i)))*(1-(Hi(i)/Hf(i))^3);
          % ELEMENTOS DE LA MATRIZ DE RIGIDECES
          raz=1/f11;
          Detx=f22*f66-f26*f26;
          r11x=f22/Detx;
          r12x=(f26*L(i)-f22)/Detx;
          r22x=(f66*(L(i))^2-2*f26*L(i)+f22)/Detx;
          raax=(r11x+r22x+2*r12x)/(L(i))^2;
          rabx=(r11x+r12x)/L(i);
          rbax=(r22x+r12x)/L(i);
          % MATRIZ DE RIGIDECES EN COORDENADAS LOCALES
          K11=[raz 0 0 ; 0 raax rabx ; 0 rabx r11x];
          K12=[-raz 0 0 ; 0 -raax rbax ; 0 -rabx r12x];
          K21=K12';
          K22=[raz 0 0 ; 0 raax -rbax ; 0 -rbax r22x];
          KE=[K11 K12 ;K21 K22];
          % MATRIZ DE TRANSPORTE A COORDENADAS GLOBALES
          T2_3=zeros(6,6);
          T2_3(1,1)=COSENO(i);T2_3(1,2)=SENO(i);
          T2_3(2,1)=-SENO(i);T2_3(2,2)=COSENO(i);
          T2_3(3,3)=1;
          T2_3(4,4)=COSENO(i);T2_3(4,5)=SENO(i);
          T2_3(5,4)=-SENO(i);T2_3(5,5)=COSENO(i);
          T2_3(6,6)=1;
          % ELEMENTO RECTANGULAR SÓLIDO DE SECCIÓN VARIABLE
      elseif DD(i)==3
          % ELEMENTOS DE LA MATRIZ DE FLEXIBILIDADES
          f11=(L(i)/(ELAS*Bb*Hi(i)))*(Hi(i)/(Hf(i)-Hi(i)))*log(Hf(i)/Hi(i));
          f22=(6*(L(i))^3/(ELAS*Bb*(Hi(i))^3))*((Hi(i)/Hf(i))^2-(Hi(i)/(Hf(i)-Hi(i)))^3*(Hf(i)/Hi(i)-Hi(i)/Hf(i)-2*log(Hf(i)/Hi(i))))+(6*L(i)/(5*GCOR*Bb*Hi(i)))*(Hi(i)/(Hf(i)-Hi(i)))*log(Hf(i)/Hi(i));
          f33=(12*(L(i))^3/(ELAS*(Bb)^3*Hi(i)))*(3/2+(1/2)*(Hf(i)/Hi(i))^2-2*(Hf(i)/Hi(i))+log(Hf(i)/Hi(i)))+(6*L(i)/(5*GCOR*Bb*Hi(i)))*(Hi(i)/(Hf(i)-Hi(i)))*log(Hf(i)/Hi(i));
          f26=(6*L(i)/(ELAS*Bb*(Hi(i))^3))*(Hi(i)/(Hf(i)-Hi(i)))^2*(1+(Hi(i)/Hf(i))^2-2*(Hi(i)/Hf(i)));
          f35=(12*(L(i))^2/(ELAS*(Bb)^3*Hi(i)))*(Hi(i)/(Hf(i)-Hi(i)))^2*(Hf(i)/Hi(i)-log(Hf(i)/Hi(i))-1);
          f66=(6*L(i)/(ELAS*Bb*(Hi(i))^3))*(Hi(i)/(Hf(i)-Hi(i)))*(1-(Hi(i)/Hf(i))^2);
          f55=(12*L(i)/(ELAS*(Bb)^3*Hi(i)))*(Hi(i)/(Hf(i)-Hi(i)))*log(Hf(i)/Hi(i));
          f44=FF44(i);
          % ELEMENTOS DE LA MATRIZ DE RIGIDECES
          raz=1/f11;
          Detx=f22*f66-f26*f26;
          r11x=f22/Detx;
          r12x=(f26*L(i)-f22)/Detx;
          r22x=(f66*(L(i))^2-2*f26*L(i)+f22)/Detx;
          raax=(r11x+r22x+2*r12x)/(L(i))^2;
          rabx=(r11x+r12x)/L(i);
          rbax=(r22x+r12x)/L(i);
          % MATRIZ DE RIGIDECES EN COORDENADAS LOCALES
          K11=[raz 0 0 ; 0 raax rabx ; 0 rabx r11x];
          K12=[-raz 0 0 ; 0 -raax rbax ; 0 -rabx r12x];
          K21=K12';
          K22=[raz 0 0 ; 0 raax -rbax ; 0 -rbax r22x];
          KE=[K11 K12 ;K21 K22];
          % MATRIZ DE TRANSPORTE A COORDENADAS GLOBALES
          T2_3=zeros(6,6);
          T2_3(1,1)=COSENO(i);T2_3(1,2)=SENO(i);
          T2_3(2,1)=-SENO(i);T2_3(2,2)=COSENO(i);
          T2_3(3,3)=1;
          T2_3(4,4)=COSENO(i);T2_3(4,5)=SENO(i);
          T2_3(5,4)=-SENO(i);T2_3(5,5)=COSENO(i);
          T2_3(6,6)=1;
          % ELEMENTO CIRCULAR SÓLIDO DE SECCIÓN CONSTATE
      elseif DD(i)==4
          C=((4*ELAS*INER(i))/L(i))*((1+FI(i))/(1+4*FI(i)));
          CP=C;
          A=((2*ELAS*INER(i))/L(i))*((1-2*FI(i))/(1+4*FI(i)));
          B=((6*ELAS*INER(i))/L(i)^2)*(1/(1+4*FI(i)));
          BP=B;
          T=(12*ELAS*INER(i))/L(i)^3*(1/(1+4*FI(i)));
          R=ELAS*AREA(i)/L(i);
          KE=[R 0 0 -R 0 0;
          0 T (B+C1(i)*T) 0 -T (B+C2(i)*T);
          0 (B+C1(i)*T) (C+2*C1(i)*B+C1(i)^2*T) 0 -(B+C1(i)*T) (A+C1(i)*B+C2(i)*B+C1(i)*C2(i)*T);
          -R 0 0 R 0 0;
          0 -T -(B+C1(i)*T) 0 T -(B+C2(i)*T);
          0 (B+C2(i)*T) (A+C1(i)*B+C2(i)*B+C1(i)*C2(i)*T) 0 -(B+C2(i)*T) (C+2*C2(i)*B+C2(i)^2*T)];
          % MATRIZ DE TRANSPORTE A COORDENADAS GLOBALES
          T2_3=zeros(6,6);
          T2_3(1,1)=COSENO(i);T2_3(1,2)=SENO(i);
          T2_3(2,1)=-SENO(i);T2_3(2,2)=COSENO(i);
          T2_3(3,3)=1;
          T2_3(4,4)=COSENO(i);T2_3(4,5)=SENO(i);
          T2_3(5,4)=-SENO(i);T2_3(5,5)=COSENO(i);
          T2_3(6,6)=1;
          % ELEMENTO CIRCULAR SÓLIDO DE SECCIÓN VARIABLE
      elseif DD(i)==5
          % ELEMENTOS DE LA MATRIZ DE FLEXIBILIDADES
          f11=(4*L(i)/(pi*ELAS*(Di(i))^2))*(Di(i)/(Df(i)-Di(i)))*(1-Di(i)/Df(i));
          f22=(64*(L(i))^3/(3*pi*ELAS*(Di(i))^4))*(Di(i)/Df(i))^3+(40*L(i)/(9*pi*GCOR*(Di(i))^2))*(Di(i)/(Df(i)-Di(i)))*(1-(Di(i)/Df(i)));
          f33=f22;
          f26=(64*(L(i))^2/(3*pi*ELAS*(Di(i))^4))*((Di(i)/Df(i))^3+(1/2)*(Di(i)/(Df(i)-Di(i)))^2*(1+(Di(i)/Df(i))^2-2*(Di(i)/Df(i))));
          f35=f26;
          f66=(64*L(i)/(3*pi*ELAS*(Di(i))^4))*(Di(i)/(Df(i)-Di(i)))*(1-(Di(i)/Df(i))^3);
          f55=f66;
          f44=(32*L(i)/(3*pi*GCOR*Di(i)^4))*(Di(i)/(Df(i)-Di(i)))*(1-(Di(i)/Df(i))^3);
          % ELEMENTOS DE LA MATRIZ DE RIGIDECES
          raz=1/f11;
          Detx=f22*f66-f26*f26;
          r11x=f22/Detx;
          r12x=(f26*L(i)-f22)/Detx;
          r22x=(f66*(L(i))^2-2*f26*L(i)+f22)/Detx;
          raax=(r11x+r22x+2*r12x)/(L(i))^2;
          rabx=(r11x+r12x)/L(i);
          rbax=(r22x+r12x)/L(i);
          % MATRIZ DE RIGIDECES EN COORDENADAS LOCALES
          K11=[raz 0 0 ; 0 raax rabx ; 0 rabx r11x];
          K12=[-raz 0 0 ; 0 -raax rbax ; 0 -rabx r12x];
          K21=K12';
          K22=[raz 0 0 ; 0 raax -rbax ; 0 -rbax r22x];
          KE=[K11 K12 ;K21 K22];
          % MATRIZ DE TRANSPORTE A COORDENADAS GLOBALES
          T2_3=zeros(6,6);
          T2_3(1,1)=COSENO(i);T2_3(1,2)=SENO(i);
          T2_3(2,1)=-SENO(i);T2_3(2,2)=COSENO(i);
          T2_3(3,3)=1;
          T2_3(4,4)=COSENO(i);T2_3(4,5)=SENO(i);
          T2_3(5,4)=-SENO(i);T2_3(5,5)=COSENO(i);
          T2_3(6,6)=1;
      end
      K2=KE+KG; % MATRIZ DE RIGIDEZ EN COORDENADAS LOCALES
      %--------------------------------------------------------------------      %
      K3=T2_3'*K2*T2_3; % MATRIZ DE RIGIDEZ EN COORDENADAS GLOBALES
      MAT_AUX=K3;
      for j=1:6
         if VC(i,j)==0
            MAT_AUX(j,:)=0;
            MAT_AUX(:,j)=0;
         end
      end
      K=zeros(ngl,ngl);
      for cont1=1:6
         for cont2=1:6
            if MAT_AUX(cont1,cont2)~=0
               UNO=VC(i,cont1);
               DOS=VC(i,cont2);
               TRES=MAT_AUX(cont1,cont2);
               K(UNO,DOS)=K(UNO,DOS)+TRES;
            end
         end
      end
      if i==1
         AUX=K;
      elseif i~=1
         AUX=AUX+K;
      else
      end
   end
   fprintf('\n MATRIZ DE RIGIDEZ DE SEGUNDO ORDEN DEL SISTEMA:');
   MAT_RIGIDEZ=AUX   
   
   fprintf('\n DESPLAZAMIENTOS GENERALIZADOS DE SEGUNDO ORDEN EN COORDENADAS GLOBALES:');
   Q_des=MAT_RIGIDEZ\Q_TOTAL
   
   %------------- ARREGLO DE DESPLAZAMIENTOS Y P_FINAL ------------------------------------
   for i=1:mbr
      fprintf('\n MIEMBRO %d:',i);
      fprintf('\n DEFORMACIONES EN MIEMBRO P1:');
      clear P1 P2 P Q2_AUX
      for j=1:6
         if VC(i,j)~=0
            P1(j)=Q_des(VC(i,j));
         else % ELIMINACIÓN DE LOS GRADOS SIN LIBERTAD EN LA ESTRUCTURA
            P1(j)=0;
         end
      end
      fprintf('\n DESPLAZAMIENTOS DEL MIEMBRO EN COORDENADAS LOCALES:');
      P1 % VECTOR DE DESPLAZAMIENTOS d
   %
   %-------------- CÁLCULO DE LOS COEFICIENTES DE RIGIDEZ DE --------------
   %------------------------ LA MATRIZ GEOMETRICA -------------------------
   rg1x=(6*F4(i))/(5*L(i));
   rg2x=F4(i)/10;
   rg3x=(2*F4(i)*L(i))/15;
   rg4x=(F4(i)*L(i))/30;
   rv1=F2(i)/L(i);
   %-------------- SE GENERA LA MATRIZ DE RIGIDEZ GEOMÉTRICA---------------
   KG=[0 rv1 0 0 -rv1 0 ; 
       rv1 rg1x rg2x -rv1 -rg1x rg2x ;
       0 rg2x rg3x 0 -rg2x -rg4x ;
       0 -rv1 0 0 rv1 0 ;
       -rv1 -rg1x -rg2x rv1 rg1x -rg2x ;
       0 rg2x -rg4x 0 -rg2x rg3x];
      if DD(i)==1
          %
          % MATRIZ DE RIGIDECES EN COORDENADAS LOCALES
          % PARA CADA CLASE DE ELEMENTO
          %
          % PARA LA GENERACIÓN DE VARIABLES SE 
          % INCLUYEN LOS SIGUIENTES CÁLCULOS
          C=((4*ELAS*INER(i))/L(i))*((1+FI(i))/(1+4*FI(i)));
          CP=C;
          A=((2*ELAS*INER(i))/L(i))*((1-2*FI(i))/(1+4*FI(i)));
          B=((6*ELAS*INER(i))/L(i)^2)*(1/(1+4*FI(i)));
          BP=B;
          T=(12*ELAS*INER(i))/L(i)^3*(1/(1+4*FI(i)));
          R=ELAS*AREA(i)/L(i);
          KE=[R 0 0 -R 0 0;
          0 T (B+C1(i)*T) 0 -T (B+C2(i)*T);
          0 (B+C1(i)*T) (C+2*C1(i)*B+C1(i)^2*T) 0 -(B+C1(i)*T) (A+C1(i)*B+C2(i)*B+C1(i)*C2(i)*T);
          -R 0 0 R 0 0;
          0 -T -(B+C1(i)*T) 0 T -(B+C2(i)*T);
          0 (B+C2(i)*T) (A+C1(i)*B+C2(i)*B+C1(i)*C2(i)*T) 0 -(B+C2(i)*T) (C+2*C2(i)*B+C2(i)^2*T)];
          % MATRIZ DE TRANSPORTE A COORDENADAS GLOBALES
          T2_3=zeros(6,6);
          T2_3(1,1)=COSENO(i);T2_3(1,2)=SENO(i);
          T2_3(2,1)=-SENO(i);T2_3(2,2)=COSENO(i);
          T2_3(3,3)=1;
          T2_3(4,4)=COSENO(i);T2_3(4,5)=SENO(i);
          T2_3(5,4)=-SENO(i);T2_3(5,5)=COSENO(i);
          T2_3(6,6)=1;
          % ELEMENTO CUADRADO SÓLIDO DE SECCIÓN VARIABLE
      elseif DD(i)==2
          % ELEMENTOS DE LA MATRIZ DE FLEXIBILIDADES
          f11=(L(i)/(ELAS*Hi(i)*Hi(i)))*(Hi(i)/(Hf(i)-Hi(i)))*(1-Hi(i)/Hf(i));
          f22=(4*(L(i))^3/(ELAS*(Hi(i))^4))*(Hi(i)/Hf(i))^3+(6*L(i)/(5*GCOR*(Hi(i))^2))*(Hi(i)/(Hf(i)-Hi(i)))*(1-Hi(i)/Hf(i));
          f33=f22;
          f26=(4*(L(i))^2/(ELAS*(Hi(i))^4))*((Hi(i)/Hf(i))^3+(1/2)*(Hi(i)/(Hf(i)-Hi(i)))^2*(1+(Hi(i)/Hf(i))^2-2*(Hi(i)/Hf(i))));
          f35=f26;
          f66=(4*L(i)/(ELAS*(Hi(i))^4))*(Hi(i)/(Hf(i)-Hi(i)))*(1-(Hi(i)/Hf(i))^3);
          f55=f66;
          f44=(64*L(i)/(27*GCOR*(Hi(i))^4))*(Hi(i)/(Hf(i)-Hi(i)))*(1-(Hi(i)/Hf(i))^3);
          % ELEMENTOS DE LA MATRIZ DE RIGIDECES
          raz=1/f11;
          Detx=f22*f66-f26*f26;
          r11x=f22/Detx;
          r12x=(f26*L(i)-f22)/Detx;
          r22x=(f66*(L(i))^2-2*f26*L(i)+f22)/Detx;
          raax=(r11x+r22x+2*r12x)/(L(i))^2;
          rabx=(r11x+r12x)/L(i);
          rbax=(r22x+r12x)/L(i);
          % MATRIZ DE RIGIDECES EN COORDENADAS LOCALES
          K11=[raz 0 0 ; 0 raax rabx ; 0 rabx r11x];
          K12=[-raz 0 0 ; 0 -raax rbax ; 0 -rabx r12x];
          K21=K12';
          K22=[raz 0 0 ; 0 raax -rbax ; 0 -rbax r22x];
          KE=[K11 K12 ;K21 K22];
          % MATRIZ DE TRANSPORTE A COORDENADAS GLOBALES
          T2_3=zeros(6,6);
          T2_3(1,1)=COSENO(i);T2_3(1,2)=SENO(i);
          T2_3(2,1)=-SENO(i);T2_3(2,2)=COSENO(i);
          T2_3(3,3)=1;
          T2_3(4,4)=COSENO(i);T2_3(4,5)=SENO(i);
          T2_3(5,4)=-SENO(i);T2_3(5,5)=COSENO(i);
          T2_3(6,6)=1;
          % ELEMENTO RECTANGULAR SÓLIDO DE SECCIÓN VARIABLE
      elseif DD(i)==3
          % ELEMENTOS DE LA MATRIZ DE FLEXIBILIDADES
          f11=(L(i)/(ELAS*Bb*Hi(i)))*(Hi(i)/(Hf(i)-Hi(i)))*log(Hf(i)/Hi(i));
          f22=(6*(L(i))^3/(ELAS*Bb*(Hi(i))^3))*((Hi(i)/Hf(i))^2-(Hi(i)/(Hf(i)-Hi(i)))^3*(Hf(i)/Hi(i)-Hi(i)/Hf(i)-2*log(Hf(i)/Hi(i))))+(6*L(i)/(5*GCOR*Bb*Hi(i)))*(Hi(i)/(Hf(i)-Hi(i)))*log(Hf(i)/Hi(i));
          f33=(12*(L(i))^3/(ELAS*(Bb)^3*Hi(i)))*(3/2+(1/2)*(Hf(i)/Hi(i))^2-2*(Hf(i)/Hi(i))+log(Hf(i)/Hi(i)))+(6*L(i)/(5*GCOR*Bb*Hi(i)))*(Hi(i)/(Hf(i)-Hi(i)))*log(Hf(i)/Hi(i));
          f26=(6*L(i)/(ELAS*Bb*(Hi(i))^3))*(Hi(i)/(Hf(i)-Hi(i)))^2*(1+(Hi(i)/Hf(i))^2-2*(Hi(i)/Hf(i)));
          f35=(12*(L(i))^2/(ELAS*(Bb)^3*Hi(i)))*(Hi(i)/(Hf(i)-Hi(i)))^2*(Hf(i)/Hi(i)-log(Hf(i)/Hi(i))-1);
          f66=(6*L(i)/(ELAS*Bb*(Hi(i))^3))*(Hi(i)/(Hf(i)-Hi(i)))*(1-(Hi(i)/Hf(i))^2);
          f55=(12*L(i)/(ELAS*(Bb)^3*Hi(i)))*(Hi(i)/(Hf(i)-Hi(i)))*log(Hf(i)/Hi(i));
          f44=FF44(i);
          % ELEMENTOS DE LA MATRIZ DE RIGIDECES
          raz=1/f11;
          Detx=f22*f66-f26*f26;
          r11x=f22/Detx;
          r12x=(f26*L(i)-f22)/Detx;
          r22x=(f66*(L(i))^2-2*f26*L(i)+f22)/Detx;
          raax=(r11x+r22x+2*r12x)/(L(i))^2;
          rabx=(r11x+r12x)/L(i);
          rbax=(r22x+r12x)/L(i);
          % MATRIZ DE RIGIDECES EN COORDENADAS LOCALES
          K11=[raz 0 0 ; 0 raax rabx ; 0 rabx r11x];
          K12=[-raz 0 0 ; 0 -raax rbax ; 0 -rabx r12x];
          K21=K12';
          K22=[raz 0 0 ; 0 raax -rbax ; 0 -rbax r22x];
          KE=[K11 K12 ;K21 K22];
          % MATRIZ DE TRANSPORTE A COORDENADAS GLOBALES
          T2_3=zeros(6,6);
          T2_3(1,1)=COSENO(i);T2_3(1,2)=SENO(i);
          T2_3(2,1)=-SENO(i);T2_3(2,2)=COSENO(i);
          T2_3(3,3)=1;
          T2_3(4,4)=COSENO(i);T2_3(4,5)=SENO(i);
          T2_3(5,4)=-SENO(i);T2_3(5,5)=COSENO(i);
          T2_3(6,6)=1;
          % ELEMENTO CIRCULAR SÓLIDO DE SECCIÓN CONSTATE
      elseif DD(i)==4
          C=((4*ELAS*INER(i))/L(i))*((1+FI(i))/(1+4*FI(i)));
          CP=C;
          A=((2*ELAS*INER(i))/L(i))*((1-2*FI(i))/(1+4*FI(i)));
          B=((6*ELAS*INER(i))/L(i)^2)*(1/(1+4*FI(i)));
          BP=B;
          T=(12*ELAS*INER(i))/L(i)^3*(1/(1+4*FI(i)));
          R=ELAS*AREA(i)/L(i);
          KE=[R 0 0 -R 0 0;
          0 T (B+C1(i)*T) 0 -T (B+C2(i)*T);
          0 (B+C1(i)*T) (C+2*C1(i)*B+C1(i)^2*T) 0 -(B+C1(i)*T) (A+C1(i)*B+C2(i)*B+C1(i)*C2(i)*T);
          -R 0 0 R 0 0;
          0 -T -(B+C1(i)*T) 0 T -(B+C2(i)*T);
          0 (B+C2(i)*T) (A+C1(i)*B+C2(i)*B+C1(i)*C2(i)*T) 0 -(B+C2(i)*T) (C+2*C2(i)*B+C2(i)^2*T)];
          % MATRIZ DE TRANSPORTE A COORDENADAS GLOBALES
          T2_3=zeros(6,6);
          T2_3(1,1)=COSENO(i);T2_3(1,2)=SENO(i);
          T2_3(2,1)=-SENO(i);T2_3(2,2)=COSENO(i);
          T2_3(3,3)=1;
          T2_3(4,4)=COSENO(i);T2_3(4,5)=SENO(i);
          T2_3(5,4)=-SENO(i);T2_3(5,5)=COSENO(i);
          T2_3(6,6)=1;
          % ELEMENTO CIRCULAR SÓLIDO DE SECCIÓN VARIABLE
      elseif DD(i)==5
          % ELEMENTOS DE LA MATRIZ DE FLEXIBILIDADES
          f11=(4*L(i)/(pi*ELAS*(Di(i))^2))*(Di(i)/(Df(i)-Di(i)))*(1-Di(i)/Df(i));
          f22=(64*(L(i))^3/(3*pi*ELAS*(Di(i))^4))*(Di(i)/Df(i))^3+(40*L(i)/(9*pi*GCOR*(Di(i))^2))*(Di(i)/(Df(i)-Di(i)))*(1-(Di(i)/Df(i)));
          f33=f22;
          f26=(64*(L(i))^2/(3*pi*ELAS*(Di(i))^4))*((Di(i)/Df(i))^3+(1/2)*(Di(i)/(Df(i)-Di(i)))^2*(1+(Di(i)/Df(i))^2-2*(Di(i)/Df(i))));
          f35=f26;
          f66=(64*L(i)/(3*pi*ELAS*(Di(i))^4))*(Di(i)/(Df(i)-Di(i)))*(1-(Di(i)/Df(i))^3);
          f55=f66;
          f44=(32*L(i)/(3*pi*GCOR*Di(i)^4))*(Di(i)/(Df(i)-Di(i)))*(1-(Di(i)/Df(i))^3);
          % ELEMENTOS DE LA MATRIZ DE RIGIDECES
          raz=1/f11;
          Detx=f22*f66-f26*f26;
          r11x=f22/Detx;
          r12x=(f26*L(i)-f22)/Detx;
          r22x=(f66*(L(i))^2-2*f26*L(i)+f22)/Detx;
          raax=(r11x+r22x+2*r12x)/(L(i))^2;
          rabx=(r11x+r12x)/L(i);
          rbax=(r22x+r12x)/L(i);
          % MATRIZ DE RIGIDECES EN COORDENADAS LOCALES
          K11=[raz 0 0 ; 0 raax rabx ; 0 rabx r11x];
          K12=[-raz 0 0 ; 0 -raax rbax ; 0 -rabx r12x];
          K21=K12';
          K22=[raz 0 0 ; 0 raax -rbax ; 0 -rbax r22x];
          KE=[K11 K12 ;K21 K22];
          % MATRIZ DE TRANSPORTE A COORDENADAS GLOBALES
          T2_3=zeros(6,6);
          T2_3(1,1)=COSENO(i);T2_3(1,2)=SENO(i);
          T2_3(2,1)=-SENO(i);T2_3(2,2)=COSENO(i);
          T2_3(3,3)=1;
          T2_3(4,4)=COSENO(i);T2_3(4,5)=SENO(i);
          T2_3(5,4)=-SENO(i);T2_3(5,5)=COSENO(i);
          T2_3(6,6)=1;
      end
      % LO ANTERIOR SE REALIZA PARA LA GENERACIÓN NUEVAMENTE DE LA MATRIZ
      % DE RIGIDEZ
      fprintf('\n RESULTADOS FINALES DE SEGUNDO ORDEN EN EL MIEMBRO %d:',i);
      fprintf('\n MATRIZ DE RIGIDEZ EN COORDENADAS LOCALES:\n');
      K2=KE+KG % K EN COORDENADAS LOCALES
      fprintf('\n MATRIZ DE RIGIDEZ EN COORDENADAS GLOBALES:\n');
      K3=T2_3'*K2*T2_3 % K EN COORDENADAS GLOBALES
      fprintf('\n VECTOR DE CARGAS EN COORDENADAS GLOBALES:\n');
      P2=K3*P1';% K*d=q CARGAS EN COORDENADAS GLOBALES
      fprintf('\n ACCIONES EN COORDENADAS LOCALES (PROBLEMA DE LA ESTRUCTURA LIBERADA):\n');
      P=T2_3*P2% R*q=q' (P.COMPLEMENTARIO) % CARGA EN COORDENADAS LOCALES 
      fprintf('\n ACCIONES DE EMPOTRAMIENTO EN COORDENADAS LOCALES (PROBLEMA DE LA ESTRUCTURA RESTRINGIDA):\n');
      for j=1:6
         Q2_AUX(j)=Q2_ALMAC(i,j); % VALORES CORRESPONDIENTES A q DE EMPOTRAMIENTO 
         %                          EN COORDENADAS GLOBALES
      end
      Q2_AUX=Q2_AUX' % q_emp
      fprintf('\n ACCIONES FINALES DEL MIEMBRO EN COORDENADAS LOCALES ( P.P.+P.C. ):');
      P_FINAL=Q2_AUX+P % ELEMENTOS MECÁNICOS FINALES: q_emp+q'
   end
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% FIN DEL PROGRAMA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
diary off
format short
fprintf('\n NO OLVIDE CONSULTAR EL DIARIO DEL CÓDIGO');
pause
open('diary')