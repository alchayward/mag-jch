%%% Generate Graphs %%%


p = struct(); %parameters structers
p.lattice_dim = [6 6];
p.model = 'BH';
p.maxParticlesPerSite = 4;
p.nParticles = 4;
p.nAtoms=1;
p.alpha = [4 36];
p.onsiteStrength = [0,0,10,0];
%p.onsiteStrength = [1,sqrt(2),-5,-10];
p.hoppingStrength = 1;
p.ddStrength = 0;
p.twist = [0 0];
p.A = [-1 0];
p.sysPath = './systems/';
p.atomLevels = 3;

    
h=HubbardLibrary;
s2=h.Lattice(p);
s2 = s2.Initilize();

ham = s2.MakeHamiltonian();
[V, D] = eigs(ham,8,'SR');
s2.V = V;
s2.D = D;

save('./saves/pfaffbh.mat','s2');

p.model = 'GJCH';
p.onsiteStrength = [sqrt(2),1,-5,-10];

s1=h.Lattice(p);
s1 = s1.Initilize();

ham = s1.MakeHamiltonian();
[V, D] = eigs(ham,8,'SR');
s1.V = V;
s1.D = D;
save('./saves/pfaffgjch.mat','s1');
