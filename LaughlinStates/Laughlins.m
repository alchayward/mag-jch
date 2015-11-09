function p = Pfaffians()
    p.new_jch_system = @new_jch_system;
    p.make_jch_ham_params = @make_jch_ham_params;
    
end
function h = make_jch_ham_params(s,kappa,detuning)
h = s.MakeHamiltonian(kappa,[1,detuning]);
end

function s = new_jch_system(lattice_dim,nParticles,detuning)
p = struct(); %parameters structers
p.lattice_dim = lattice_dim;
p.model = 'GJCH';
p.maxParticlesPerSite = nParticles;
p.nParticles = nParticles;
p.nAtoms=1;
p.alpha = [nParticles 2*prod(lattice_dim)];
p.onsiteStrength = [1,detuning];
p.hoppingStrength = 1;
p.ddStrength = 0;
p.twist = [0 0];
p.A = [-1 0];
p.sysPath = './systems/';
p.atomLevels = 2;

h=HubbardLibrary;
s=h.Lattice(p);
s=s.Initilize();
end


