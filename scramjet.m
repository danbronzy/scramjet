function SCRAMJET(gamma1, T1, P1)
global gamma;
close all
gamma = gamma1;
AMach=zeros(3,20);
ATt=zeros(3,20);
AT=zeros(3,20);
APt=zeros(3,20);
AP=zeros(3,20);
for M1=8:10

    %engine parameters
    turn1 = 7.82;%degrees
    turn2 =8.8;%degrees
    turn3 = turn1 + turn2;
    thetas = [turn1,turn2,turn3]*pi/180;

    Ms = zeros(4,1);
    Ts = zeros(4,1);
    Ps = zeros(4,1);
    Pts = zeros(4,1);

    %input conditions
    Ms(1) = M1;
    Ps(1) = P1;
    Ts(1) = T1;
    Tt = Ts(1)*(1 + (gamma-1)/2*Ms(1)^2);
    Pts(1) = Ps(1) * (1 + (gamma - 1)/2*Ms(1)^2)^(gamma/(gamma - 1));


%after oblique 1
    for station = [2,3,4]
        beta = findBeta(Ms(station - 1), thetas(station - 1));
        Ms(station) = findNextM(Ms(station - 1), beta, thetas(station - 1));
        Ps(station) = Ps(station-1)*findPRatio(Ms(station - 1), beta);
        Ts(station) = Ts(station - 1)*findTRatio(Ms(station - 1), beta);
        Pts(station) = Pts(station - 1)*findPtRatio(Ms(station - 1), beta);
    end

%Final Conditions of inlet
    Pt4 = Pts(4);
    M4 = Ms(4);
    P4 = Ps(4);
    T4 = Ts(4);
    Tt4 = Tt;


%%Combustor
    Heat=[100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500,1600,1700,1800,1900,2000];
    step=100;
    CombM=[];
    CombTt=[];
    CombPt=[];
    CombP=[];
    CombT=[];
    %Loop that runs through the different values for Q input
    for Qin = 100:step:2000
        Cp=1.005;
        Tt5=Qin/Cp+Tt4;
        Pt5=sqrt(Tt4).*Pt4./sqrt(Tt5);


%MFP Mach Number Derivation
        syms Ma
        Mach=((1+gamma*M4^2)/(1+gamma*Ma^2))^2*(Ma/M4)^2*((1+((gamma-1)/2)*Ma^2)/(1+((gamma-1)/2)*M4^2))==Tt5/Tt4;
        solvething=solve(Mach,Ma);
        ginger=vpa(solvething,4);
        M5=ginger(4);

        P5=(1+gamma*M4^2)*P4/(1+gamma*M5^2);
        T5=((1+gamma*M4^2)/(1+gamma*M5^2))^2*(M5/M4)^2*T4;
 %Appending calculated values to a matrix
        CombM(Qin/step)=M5;
        CombTt(Qin/step)=Tt5;
        CombPt(Qin/step)=Pt5;
        CombP(Qin/step)=P5;
        CombT(Qin/step)=T5;
    end
    %Combining different Mach number values into one graph
    AMach(M1-7,:) = CombM;
    ATt(M1-7,:) = CombTt;
    AT(M1-7,:) = CombT;
    APt(M1-7,:) = CombPt;
    AP(M1-7,:) = CombP;

end
%AMach
figure
plot(Heat,AMach(1,:),'LineWidth',3)
hold on
plot(Heat,AMach(2,:),'LineWidth',3)
plot(Heat,AMach(3,:),'LineWidth',3)
title( 'Mach vs. Heat Addition' )
xlabel( 'Q(kJ)' )
ylabel( 'Mach Number' )
legend( 'M=8' , 'M=9' , 'M=10')

figure
plot(Heat,ATt(1,:),'LineWidth',3)
hold on
plot(Heat,ATt(2,:),'LineWidth',3)
plot(Heat,ATt(3,:),'LineWidth',3)
title( 'Total Temperature vs. Heat Addition' )
xlabel( 'Q(kJ)' )
ylabel( 'Total Temperature (K)' )
legend( 'M=8' , 'M=9' , 'M=10')

figure
plot(Heat,APt(1,:),'LineWidth',3)
hold on
plot(Heat,APt(2,:),'LineWidth',3)
plot(Heat,APt(3,:),'LineWidth',3)
title( 'Total Pressure vs. Heat Addition' )
xlabel( 'Q(kJ)' )
ylabel( 'Total Pressure (Pa)' )
legend( 'M=8' , 'M=9' , 'M=10')

figure
plot(Heat,AP(1,:),'LineWidth',3)
hold on
plot(Heat,AP(2,:),'LineWidth',3)
plot(Heat,AP(3,:),'LineWidth',3)
title( 'Static Pressure vs. Heat' )
xlabel( 'Q(kJ)' )
ylabel( 'Pressure (Pa)' )
legend( 'M=8' , 'M=9' , 'M=10')

figure
plot(Heat,AT(1,:),'LineWidth',3)
hold on
plot(Heat,AT(2,:),'LineWidth',3)
plot(Heat,AT(3,:),'LineWidth',3)
title( 'Temperature vs. Heat' )
xlabel( 'Q(kJ)' )
ylabel( 'Temperature (K)' )
legend( 'M=8' , 'M=9' , 'M=10')
end
function beta = findBeta(M, theta)
    global gamma;
    syms b;
    eqn = tan(theta) == 2*cot(b)*((M^2*sin(b)^2 - 1)/(M^2*(gamma + cos(2*b)) + 2));
    beta = vpasolve(eqn, b,pi/180*[0,90]);
end

function nextM = findNextM(M, beta, theta)
    global gamma;
    mb2 = M^2*sin(beta)^2;
    nextM = 1/sin(beta - theta)*sqrt((1 + (gamma - 1)/2*mb2)/...
        (gamma*mb2 - (gamma - 1)/2));
end
function pRatio = findPRatio(M, beta)
    global gamma;
    mb2 = M^2*sin(beta)^2;
    pRatio = 1 + 2*gamma/(gamma+1)*(mb2 - 1);
end
function TRatio = findTRatio(M, beta)
    global gamma;
    mb2 = M^2*sin(beta)^2;
    TRatio = ((2*gamma*mb2 - (gamma-1))*((gamma-1)*mb2 + 2))...
        /((gamma+1)^2*mb2);
end
function PtRatio = findPtRatio(M, beta)
    global gamma;
    mb2 = M^2*sin(beta)^2;
    PtRatio = ((gamma+1)*mb2/((gamma-1)*mb2 + 2))^(gamma/(gamma-1))...
        *((gamma+1)/(2*gamma*mb2-(gamma-1)))^(1/(gamma-1));
end
