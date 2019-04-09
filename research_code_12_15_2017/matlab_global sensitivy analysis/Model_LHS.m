%% The results should be compared to the PRCC results section in
%% Supplementary Material D and Table D.1 for different N (specified by
%% "runs" in the script below
clear all;
close all;

%% Sample size N
runs = 100;

%% LHS MATRIX  %%
parameter;
bootmax = 50;
tp=length(time_points);
prcc_matrix = zeros(bootmax,52);

% old_dir = 'pcos_norm';  % current directory

id = 'MATLAB:nearlySingularMatrix';
warning('off',id);

%start_time = tic;
for boots=41:bootmax

CFE_LHS	=lhsnorm(	0.0022729	,	0.000045458	,runs);
CFI_LHS	=lhsnorm(	1.9488	,	0.038976	,runs);
CFP_LHS	=lhsnorm(	60.428	,	1.20856	,runs);
CLE_LHS	=lhsnorm(	0.0010404	,	0.000020808	,runs);
CLP_LHS	=lhsnorm(	0.0099415	,	0.00019883	,runs);
CLT_LHS	=lhsnorm(	0.0095942	,	0.000191884	,runs);
deltaF_LHS	=lhsnorm(	8.21	,	0.1642	,runs);
deltal_LHS	=lhsnorm(	14	,	0.28	,runs);
kF_LHS	=lhsnorm(	2.5412	,	0.050824	,runs);
KFI_LHS	=lhsnorm(	107.01	,	2.1402	,runs);
KiLP_LHS	=lhsnorm(	0.34952	,	0.0069904	,runs);
KLT_LHS	=lhsnorm(	420	,	8.4	,runs);
KmL_LHS	=lhsnorm(	183.56	,	3.6712	,runs);
kL_LHS	=lhsnorm(	0.74567	,	0.0149134	,runs);
n_LHS	=lhsnorm(	8	,	0.16	,runs);
V_LHS	=lhsnorm(	2.5	,	0.05	,runs);
v0L_LHS	=lhsnorm(	1051.7	,	21.034	,runs);
v1L_LHS	=lhsnorm(	34838	,	696.76	,runs);
vF_LHS	=lhsnorm(	3236.6	,	64.732	,runs);
cphiF_LHS	=lhsnorm(	0.01127	,	0.0002254	,runs);
ss_LHS	=lhsnorm(	0.74702	,	0.0149404	,runs);
f0_LHS	=lhsnorm(	0.0025112	,	0.000050224	,runs);
f1_LHS	=lhsnorm(	4.3764	,	0.087528	,runs);
f2_LHS	=lhsnorm(	27.812	,	0.55624	,runs);
h1_LHS	=lhsnorm(	590.32	,	11.8064	,runs);
h2_LHS	=lhsnorm(	1815.3	,	36.306	,runs);
hp_LHS	=lhsnorm(	20.764	,	0.41528	,runs);
hs_LHS	=lhsnorm(	12.329	,	0.24658	,runs);
l_LHS	=lhsnorm(	0.49017	,	0.0098034	,runs);
m_LHS	=lhsnorm(	4	,	0.08	,runs);
s_LHS	=lhsnorm(	2.378	,	0.04756	,runs);
w_LHS	=lhsnorm(	0.23173	,	0.0046346	,runs);
CTF2_LHS	=lhsnorm(	123.8136	,	2.476272	,runs);
deltaE_LHS	=lhsnorm(	1.1	,	0.022	,runs);
deltaP_LHS	=lhsnorm(	0.5	,	0.01	,runs);
deltaT_LHS	=lhsnorm(	5.5	,	0.11	,runs);
e0_LHS	=lhsnorm(	44.512	,	0.89024	,runs);
eta_LHS	=lhsnorm(	1.1087	,	0.022174	,runs);
h3_LHS	=lhsnorm(	17.796	,	0.35592	,runs);
k1_LHS	=lhsnorm(	1.09	,	0.0218	,runs);
k2_LHS	=lhsnorm(	22.28645	,	0.445729	,runs);
k3_LHS	=lhsnorm(	113.9188	,	2.278376	,runs);
p_LHS	=lhsnorm(	0.3734	,	0.007468	,runs);
t0_LHS	=lhsnorm(	741.68	,	14.8336	,runs);
t1_LHS	=lhsnorm(	0.57088	,	0.0114176	,runs);
t2_LHS	=lhsnorm(	1.3481	,	0.026962	,runs);
tao1_LHS	=lhsnorm(	5.3989	,	0.107978	,runs);
tao2_LHS	=lhsnorm(	0	,	0	,runs);
tao3_LHS	=lhsnorm(	430.91	,	8.6182	,runs);
tg1_LHS	=lhsnorm(	6.6548	,	0.133096	,runs);
tg2_LHS	=lhsnorm(	186.27	,	3.7254	,runs);
psi_LHS	=lhsnorm(	2004.3	,	40.086	,runs);

%% LHS MATRIX and PARAMETER LABELS
dummy_LHS=lhsdesign(runs,1)*9 + 1;
LHSmatrix(:,	1	)=	CFE_LHS	;
LHSmatrix(:,	2	)=	CFI_LHS	;
LHSmatrix(:,	3	)=	CFP_LHS	;
LHSmatrix(:,	4	)=	CLE_LHS	;
LHSmatrix(:,	5	)=	CLP_LHS	;
LHSmatrix(:,	6	)=	CLT_LHS	;
LHSmatrix(:,	7	)=	deltaF_LHS	;
LHSmatrix(:,	8	)=	deltal_LHS	;
LHSmatrix(:,	9	)=	kF_LHS	;
LHSmatrix(:,	10	)=	KFI_LHS	;
LHSmatrix(:,	11	)=	KiLP_LHS	;
LHSmatrix(:,	12	)=	KLT_LHS	;
LHSmatrix(:,	13	)=	KmL_LHS	;
LHSmatrix(:,	14	)=	kL_LHS	;
LHSmatrix(:,	15	)=	n_LHS	;
LHSmatrix(:,	16	)=	V_LHS	;
LHSmatrix(:,	17	)=	v0L_LHS	;
LHSmatrix(:,	18	)=	v1L_LHS	;
LHSmatrix(:,	19	)=	vF_LHS	;
LHSmatrix(:,	20	)=	cphiF_LHS	;
LHSmatrix(:,	21	)=	ss_LHS	;
LHSmatrix(:,	22	)=	f0_LHS	;
LHSmatrix(:,	23	)=	f1_LHS	;
LHSmatrix(:,	24	)=	f2_LHS	;
LHSmatrix(:,	25	)=	h1_LHS	;
LHSmatrix(:,	26	)=	h2_LHS	;
LHSmatrix(:,	27	)=	hp_LHS	;
LHSmatrix(:,	28	)=	hs_LHS	;
LHSmatrix(:,	29	)=	l_LHS	;
LHSmatrix(:,	30	)=	m_LHS	;
LHSmatrix(:,	31	)=	s_LHS	;
LHSmatrix(:,	32	)=	w_LHS	;
LHSmatrix(:,	33	)=	CTF2_LHS	;
LHSmatrix(:,	34	)=	deltaE_LHS	;
LHSmatrix(:,	35	)=	deltaP_LHS	;
LHSmatrix(:,	36	)=	deltaT_LHS	;
LHSmatrix(:,	37	)=	e0_LHS	;
LHSmatrix(:,	38	)=	eta_LHS	;
LHSmatrix(:,	39	)=	h3_LHS	;
LHSmatrix(:,	40	)=	k1_LHS	;
LHSmatrix(:,	41	)=	k2_LHS	;
LHSmatrix(:,	42	)=	k3_LHS	;
LHSmatrix(:,	43	)=	p_LHS	;
LHSmatrix(:,	44	)=	t0_LHS	;
LHSmatrix(:,	45	)=	t1_LHS	;
LHSmatrix(:,	46	)=	t2_LHS	;
LHSmatrix(:,	47	)=	tao1_LHS	;
%LHSmatrix(:,	48	)=	tao2_LHS	;
LHSmatrix(:,	48	)=	tao3_LHS	;
LHSmatrix(:,	49	)=	tg1_LHS	;
LHSmatrix(:,	50	)=	tg2_LHS	;
LHSmatrix(:,	51	)=	psi_LHS	;
LHSmatrix(:,	52	)=dummy_LHS;

FSH_lhs = zeros(tp,runs);
LH_lhs = FSH_lhs;
T_lhs = FSH_lhs;
E2_lhs = FSH_lhs;
P4_lhs = FSH_lhs;


for x=1:runs %Run solution x times choosing different values
%     f=@ODE_LHS;
    %x
    %LHSmatrix(x,:)
    dlmwrite('progress.txt',['boot:',num2str(boots),'; run:',num2str(x)],'delimiter','');

    options = odeset('RelTol',1e-3,'AbsTol',1e-8);

    pars = abs(LHSmatrix(x,:));
    [~,y1]=ode23s(@(t,y) ODE_LHS(t,y,pars),[0 6*31],y0,options);
    init=y1(end,:);
   out= ode23s(@(t,y) ODE_LHS(t,y,pars),tspan,init,options);

   A = deval(out,time_points)';  % [time y]
%A=[y];
     %% Save the outputs at ALL time points [tspan]
%FSH_lhs(:,x)=	Anew(:,	2	);LH_lhs(:,x)=	Anew(:,3	);FSHp_lhs(:,x)=	Anew(:,4	);LHp_lhs(:,x)=	Anew(:,5	);phi_lhs(:,x)=	Anew(:,6	);
%omega_lhs(:,x)=Anew(:,7	);lamda_lhs(:,x)=	Anew(:,8	);S_lhs(:,x)=	Anew(:,9	);Ty_lhs(:,x)=Anew(:,10	);
%T_lhs(:,x)=	Anew(:,11	);E2_lhs(:,x)=Anew(:,12	);P4_lhs(:,x)=Anew(:,13	);
    %% Save only the outputs at the time points of interest [time_points]:
    %% MORE EFFICIENT
    FSH_lhs(:,x)=	A(:,	1	);
    LH_lhs(:,x)=	A(:,	2	);
    % FSHp_lhs(:,x)=	A(index_points,	4	);LHp_lhs(:,x)=	A(index_points,	5	);
    % phi_lhs(:,x)=	A(index_points,	6	);omega_lhs(:,x)=	A(index_points,	7	);lamda_lhs(:,x)=	A(index_points,	8	);S_lhs(:,x)=	A(index_points,9	);
    % Ty_lhs(:,x)=	A(index_points,	10	);
    T_lhs(:,x)=	A(:,10	);
    E2_lhs(:,x)=	A(:,	11	);
    P4_lhs(:,x)=	A(:,	12	);

end
%end_time = [num2str(toc/60),' minutes elapsed.'];
%dlmwrite('run_time.txt',end_time,'delimiter','');
%% Save the workspace


% CALCULATE PRCC
[prcc,~,~]=PRCC(LHSmatrix,FSH_lhs,1:length(time_points),PRCC_var,0);

%T=prcc.*((runs-2-3)./(1-prcc.^2)).^(1/2); T_matrix(boots,:)=T;
%CC_PLOT(LHSmatrix,Q_lhs,1:length(time_points),'lin',PRCC_var,y_var_label)
%RCC_PLOT(LHSmatrix,Q_lhs,1:length(time_points),'lin',PRCC_var,y_var_label)
%PRCC_PLOT(LHSmatrix,Q_lhs,1:length(time_points),PRCC_var,y_var_label)

% for i=1:tp
% prcc_matrix((boots-1)*tp+(1:tp),:)=prcc;
prcc_matrix = prcc;
%%prcc_matrix(i,:)=prcc(i,:);
% end
%filename = ['Pmatrix', int2str(boots), '.csv'];
%dlmwrite(filename,prcc_matrix,'\t');

%prcc_matrix(boots,:)=prcc;
current_string = ['b',num2str(boots)];
mkdir(current_string)
cd(current_string)

dlmwrite(['FSH.dat'],FSH_lhs,'delimiter','\t');
dlmwrite(['LH.dat'],LH_lhs,'delimiter','\t');
dlmwrite(['T.dat'],T_lhs,'delimiter','\t');
dlmwrite(['E2.dat'],E2_lhs,'delimiter','\t');
dlmwrite(['P4.dat'],P4_lhs,'delimiter','\t');
dlmwrite(['LHSmatrix.dat'],LHSmatrix,'delimiter','\t');
dlmwrite('Pmatrix.txt',prcc_matrix,'delimiter','\t'); % rows: time points (63); columns: parameters (52)
cd('../')
end

% save Model_LHS.mat;
