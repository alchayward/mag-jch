function p = Pfaffians()
    p.new_jch_system = @new_jch_system;
    p.make_jch_ham_params = @make_jch_ham_params;
    
end
function h = make_jch_ham_params(s,kappa,detuning,beta_2)
h = s.MakeHamiltonian(kappa,[1,beta_2,detuning,2*detuning]);
end

function s = new_jch_system(lattice_dim,nParticles,detuning,beta_2)
p = struct(); %parameters structers
p.lattice_dim = lattice_dim;
p.model = 'GJCH';
p.maxParticlesPerSite = nParticles;
p.nParticles = nParticles;
p.nAtoms=1;
p.alpha = [nParticles prod(lattice_dim)];
p.onsiteStrength = [1,beta_2,detuning,2*detuning];
%p.onsiteStrength = [1,sqrt(2),-5,-10];
p.hoppingStrength = 1;
p.ddStrength = 0;
p.twist = [0 0];
p.A = [-1 0];
p.sysPath = './systems/';
p.atomLevels = 3;

h=HubbardLibrary;
s=h.Lattice(p);
s=s.Initilize();
end


