function var=odeEquation(t,y,a)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
var = zeros(12,1);
FSH	=	y	(	1	)	;
LH	=	y	(	2	)	;
FSHp	=	y	(	3	)	;
LHp	=	y	(	4	)	;
phi	=	y	(	5	)	;
omega	=	y	(	6	)	;
lamda	=	y	(	7	)	;
S	=	y	(	8	)	;
Ty	=	y	(	9	)	;
T	=	y	(	10	)	;
E2	=	y	(	11	)	;
P4	=	y	(	12	)	;


CFE	=	a	(	1	)	;
CFI	=	a	(	2	)	;
CFP	=	a	(	3	)	;
CLE	=	a	(	4	)	;
CLP	=	a	(	5	)	;
CLT	=	a	(	6	)	;
deltaF	=	a	(	7	)	;
deltaL	=	a	(	8	)	;
kF	=	a	(	9	)	;
KFI	=	a	(	10	)	;
KiLP	=	a	(	11	)	;
KLT	=	a	(	12	)	;
KmL	=	a	(	13	)	;
kL	=	a	(	14	)	;
n	=	a	(	15	)	;
V	=	a	(	16	)	;
v0L	=	a	(	17	)	;
v1L	=	a	(	18	)	;
vF	=	a	(	19	)	;
cphiF	=	a	(	20	)	;
deltas	=	a	(	21	)	;
f0	=	a	(	22	)	;
f1	=	a	(	23	)	;
f2	=	a	(	24	)	;
h1	=	a	(	25	)	;
h2	=	a	(	26	)	;
hp	=	a	(	27	)	;
hs	=	a	(	28	)	;
l	=	a	(	29	)	;
m	=	a	(	30	)	;
s	=	a	(	31	)	;
w	=	a	(	32	)	;
CTF2	=	a	(	33	)	;
deltaE	=	a	(	34	)	;
deltaP	=	a	(	35	)	;
deltaT	=	a	(	36	)	;
e0	=	a	(	37	)	;
eta	=	a	(	38	)	;
h3	=	a	(	39	)	;
k1	=	a	(	40	)	;
k2	=	a	(	41	)	;
k3	=	a	(	42	)	;
p	=	a	(	43	)	;
t0	=	a	(	44	)	;
t1	=	a	(	45	)	;
t2	=	a	(	46	)	;
tao1	=	a	(	47	)	;
tao2	=	a	(	48	)	;
tao3	=	a	(	49	)	;
tg1	=	a	(	50	)	;
tg2	=	a	(	51	)	;
psi	=	a	(	52	)	;
G1=1;
G2=1;

var(3)=vF/(1+CFI*S*lamda/(KFI+S*lamda))-kF*(1+CFP*P4)/(1+CFE*E2^2)*FSHp;
var(1)=1/V*kF*(1+CFP*P4)/(1+CFE*E2^2)*FSHp-deltaF*FSH;
var(4)=(v0L*T/(KLT+T)+v1L*(E2.^n)/(KmL.^n+E2.^n))/(1+P4/(KiLP*(1+CLT*T)))-kL*(1+CLP*P4)/(1+CLE*E2)*LHp;
var(2)=1/V*kL*(1+CLP*P4)/(1+CLE*E2)*LHp-deltaL*LH;
var(5)=f0*T/273.67+(f1*FSH^2/((h1/(1+0.19878*T/273.67))^2+FSH^2)-f2*LH^2/((h2/(1+cphiF*FSH))^2+LH^2))*phi;
var(6)=(f2*LH^2/((h2/(1+cphiF*FSH))^2+LH^2))*phi-w*S*omega;
var(7)=w*S*omega-l*(1-S)*lamda;
var(8)=s*(LH.^m)/((hs.^m)+(LH.^m))*(1-S)-deltas*S;
D=k1*LH^2+k2*LH+k3;
F1=(LH^2)/D;
F2=LH/D;
var(10)=t0-deltaT*T+(t1*G1*(F1+CTF2*F2)+t2*G1*G2*F1)*(phi+tao1*omega+tao2*S*lamda+tao3*(1-(phi+omega+lamda)/psi));
var(9)=tg1*G1*G2*F1-tg2*FSH/(h3+FSH)*Ty;
var(11)=e0-deltaE*E2+tg2*FSH/(h3+FSH)*Ty*(phi+eta*lamda*S);
var(12)=-deltaP*P4+p*LH/(LH+hp)*lamda*S;

end

