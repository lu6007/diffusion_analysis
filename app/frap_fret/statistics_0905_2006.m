% Use multiple comparison statistics to compare the diffusion coefficients
% between Lyn, Kras-Src biosensors and the CD treatment.
function statistics_09_05_2006()
%data 07/13/2006
cyto = [0.7821; 1.1625; 0.7529; 0.8838; 0.6405; 0.9361]; % cyto_2064
t_cyto = get_tag(cyto,'cyto');

%07/13/2006
lyn = [0.0631; 0.0598; 0.0827; 0.0962; 0.0678; 0.0649; ... %mem14
    0.1532; 0.1450; 0.1616; 0.1773; 0.2629; 0.3287;... %mem17
    0.2276; 0.1289; 0.0959; ... %mem4
    0.0698; 0.0904; 0.0985; 0.1261; 0.1140; 0.1288; ... %mem5
    0.0868; 0.0969; 0.1032; 0.1293; 0.1120; 0.1298]; %mem16
%08/29/2006
lyn = [lyn; 0.1206; 0.0816; 0.0604; 0.0851;... %24_lyn31
    0.0658; 0.0398;0.1068;0.0453;0.0776;... %24_lyn41
    0.0627; 0.0653; 0.1798; 0.1637;... %24_lyn51
    0.0809;0.0931;0.1126]; %24_lyn61
t_lyn = get_tag(lyn,'lyn');


%12/14/2006
kras = [0.3170; 0.3656; 0.0896; 0.3974; 0.1424; ... % 24_kras11
    0.1077; 0.1103; 0.1730; 0.1477; 0.1828; 0.1438;... % 24_kras21
    0.1599; 0.1222; 0.1515; 0.1139; 0.1452; 0.1149]; % 24_kras31
t_kras = get_tag(kras,'kras');

lyn_bcd = [0.1751; 0.1701; 0.1420; 0.1295;...
    0.1183; 0.3204; 0.2642; 0.2074; 0.2521; ...
    0.1602; 0.1231; 0.0928; 0.0977; 0.2300; ...
    0.1597; 0.1511; 0.1615; 0.1691; 0.1945; 0.1147];
t_lyn_bcd = get_tag(lyn_bcd, 'lyn bcd');

kras_bcd = [0.2376; 0.2880; 0.3408; 0.2261;...
    0.1715; 0.2166; 0.2333; 0.2441; 0.1212; 0.1520; ...
    0.1420; 0.1261; 0.1690; 0.1995; 0.236; 0.2129; ...
    0.1457; 0.1708; 0.1420; 0.2120; 0.2581; 0.2617];
t_kras_bcd = get_tag(kras_bcd,'kras bcd');

lyn_cytod = [0.0438; 0.0327; 0.0434; 0.0333; 0.0642; 0.0541; ...
    0.0661; 0.0781;0.1020; 0.1472; 0.0679;...
    0.0478; 0.0691; 0.0475; 0.0445; 0.0740; ...
    0.0946; 0.0429; 0.0474; 0.0642; 0.0553; 0.0768;...
    0.0467; 0.0662; 0.2357; 0.1407; 0.089; ...
    0.0604; 0.0689; 0.0668; 0.0717; 0.1109];
t_lyn_cytod = get_tag(lyn_cytod, 'lyn cytoD');

lyn_nacodozole = [ 0.1108; 0.1401; 0.1599; 0.1097; 0.137; ...
    0.0848; 0.1540; 0.0775; 0.0946];
t_lyn_nacodozole = get_tag(lyn_nacodozole, 'lyn naco');

d_data = [lyn; kras; lyn_bcd; kras_bcd]; %lyn_cytod; lyn_nacodozole];%; cyto];
t_tag = char(t_lyn, t_kras, t_lyn_bcd, t_kras_bcd);%, ...
    %t_lyn_cytod, t_lyn_nacodozole); %, t_cyto);
[p,a,stat] = anova1(d_data, t_tag);
[c,m,h,nms]=multcompare(stat, 'ctype', 'bonferroni', 'alpha', 0.05);
p
c
m
h
nms

function vec_tag = get_tag(vec, tag_name);
n = size(vec,1);
vec_tag = tag_name;
for i = 1:n-1,
    vec_tag = [tag_name; vec_tag];
end;



