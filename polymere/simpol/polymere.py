#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Moteur tridimensionnel, commun aux deux voies.

Une bille vaut `taille_bin` paires de bases. Les forces permanentes sont
celles d'usage pour la chromatine gros-grain : liaison harmonique le long
de la chaine, rigidite angulaire, repulsion molle a courte portee qui
autorise un pas de temps long, et confinement spherique a la densite
nucleaire.

S'y ajoutent, selon la voie :
  A  des liens transitoires entre les deux jambes de chaque facteur
     d'extrusion, mis a jour a chaque bloc depuis la trajectoire 1D ;
  B  une attraction entre billes partageant un type de liant, d'intensite
     proportionnelle au poids du type.

Les deux peuvent etre actives ensemble, ce qui donne une troisieme voie a
comparer aux deux autres.

POINTS DELICATS, APPRIS A L'USAGE
  La masse vaut 100 amu : a 1 amu, la vitesse thermique est telle que
  l'integrateur diverge au pas de 70 fs.
  Le confinement porte sur sqrt(x^2+y^2+z^2+1e-6) : sans le terme constant
  la derivee est indefinie pour une bille posee a l'origine.
  Les liens d'extrusion sont montes en raideur progressivement : depuis une
  configuration aleatoire, deux billes distantes de cinquante rangs
  ramenees d'emblee au contact font diverger l'integrateur.
"""
import numpy as np
import openmm
import openmm.unit as unit


def marche_depart(n, rng, pas=1.0):
    """Marche aleatoire grossierement auto-evitante."""
    p = np.zeros((n, 3))
    p[0] = rng.normal(0, 0.3, 3)
    for i in range(1, n):
        q = p[i - 1]
        for _ in range(60):
            v = rng.normal(0, 1, 3)
            v = v / np.linalg.norm(v) * pas
            q = p[i - 1] + v
            if i < 3 or np.min(np.linalg.norm(p[:i - 1] - q, axis=1)) > 0.7:
                break
        p[i] = q
    return p


class Polymere:
    def __init__(self, n, densite=0.2, paires_lef=None, liants=None,
                 poids=None, eps_liant=0.9, portee_liant=1.6,
                 plateforme=None, graine=0, pas_fs=70.0):
        """`paires_lef` est l'union de tous les couples de billes qu'un
        facteur d'extrusion reliera au cours de la trajectoire. OpenMM
        interdit de changer les billes d'une liaison deja posee : on les
        pose donc toutes une fois pour toutes, a raideur nulle, et on ne
        bascule que la raideur. C'est la methode de polychrom."""
        self.n = int(n)
        self.rng = np.random.default_rng(graine)
        self.rayon = (3 * self.n / (4 * np.pi * densite)) ** (1 / 3)

        syst = openmm.System()
        for _ in range(self.n):
            syst.addParticle(100.0 * unit.amu)

        liens = openmm.HarmonicBondForce()
        for i in range(self.n - 1):
            liens.addBond(i, i + 1, 1.0, 2.0 / (0.05 ** 2))
        self.map_lef = {}
        for i, j in (paires_lef or []):
            cle = (int(min(i, j)), int(max(i, j)))
            if cle not in self.map_lef and cle[0] != cle[1]:
                self.map_lef[cle] = liens.getNumBonds()
                liens.addBond(cle[0], cle[1], 1.0, 0.0)
        syst.addForce(liens)
        self.liens = liens
        self.k_lef = 2.0 / (0.12 ** 2)
        self.actifs = set()

        angles = openmm.CustomAngleForce("ka*(1 - cos(theta - 3.141592))")
        angles.addGlobalParameter("ka", 5.0)
        for i in range(self.n - 2):
            angles.addAngle(i, i + 1, i + 2)
        syst.addForce(angles)

        rep = openmm.CustomNonbondedForce("step(r0 - r)*e*(1 - r/r0)^2")
        rep.addGlobalParameter("e", 7.5)
        rep.addGlobalParameter("r0", 1.0)
        rep.setNonbondedMethod(openmm.CustomNonbondedForce.CutoffNonPeriodic)
        rep.setCutoffDistance(1.0)
        for _ in range(self.n):
            rep.addParticle([])
        syst.addForce(rep)

        # voie B : attraction entre billes partageant un type de liant
        self.n_types = 0
        if liants is not None and liants.size and poids is not None:
            s = np.asarray(liants, dtype=float)
            w = np.asarray(poids, dtype=float)
            w = w / (w.max() + 1e-12)
            self.n_types = s.shape[1]
            termes = " + ".join(f"c{c}1*c{c}2*{w[c]:.6f}"
                                for c in range(self.n_types))
            att = openmm.CustomNonbondedForce(
                f"-eps*({termes})*step(rc - r)*(1 - (r/rc)^2)^2")
            att.addGlobalParameter("eps", float(eps_liant))
            att.addGlobalParameter("rc", float(portee_liant))
            for c in range(self.n_types):
                att.addPerParticleParameter(f"c{c}")
            att.setNonbondedMethod(openmm.CustomNonbondedForce.CutoffNonPeriodic)
            att.setCutoffDistance(float(portee_liant))
            for i in range(self.n):
                att.addParticle([float(x) for x in s[i]])
            syst.addForce(att)

        conf = openmm.CustomExternalForce(
            "step(r - R)*kc*(r - R)^2; r = sqrt(x*x + y*y + z*z + 1e-6)")
        conf.addGlobalParameter("kc", 25.0)
        conf.addGlobalParameter("R", self.rayon)
        for i in range(self.n):
            conf.addParticle(i, [])
        syst.addForce(conf)

        self.integ = openmm.LangevinMiddleIntegrator(
            300 * unit.kelvin, 1.0 / unit.picosecond, pas_fs * unit.femtosecond)
        self.integ.setRandomNumberSeed(int(graine) % 2147483647)
        if plateforme is None:
            noms = [openmm.Platform.getPlatform(i).getName()
                    for i in range(openmm.Platform.getNumPlatforms())]
            plateforme = next(p for p in ("CUDA", "OpenCL", "CPU", "Reference")
                              if p in noms)
        self.plateforme = plateforme
        self.graine = int(graine) % 2147483647
        plat = openmm.Platform.getPlatformByName(plateforme)
        if plateforme == "CPU":
            # Une tache du job array ne recoit qu'un coeur : un seul fil
            # colle a l'allocation Slurm, et l'ordre de sommation des
            # forces ne depend plus du nombre de fils disponibles. Sans
            # cela, deux appels a la meme ligne de commande divergeaient.
            self.ctx = openmm.Context(syst, self.integ, plat,
                                      {"Threads": "1"})
        else:
            self.ctx = openmm.Context(syst, self.integ, plat)

    # ------------------------------------------------------------------
    def demarrer(self, positions=None, liens_initiaux=None):
        p = marche_depart(self.n, self.rng) if positions is None else positions
        self.ctx.setPositions(p)
        # graine explicite : sans elle, le tirage des vitesses initiales
        # est aleatoire a chaque appel et la trajectoire n'est pas rejouable
        self.ctx.setVelocitiesToTemperature(300 * unit.kelvin, self.graine)
        openmm.LocalEnergyMinimizer.minimize(self.ctx, maxIterations=800)
        if liens_initiaux:
            for frac in (0.02, 0.1, 0.3, 0.6, 1.0):
                self.poser_liens(liens_initiaux, frac)
                self.integ.step(400)
        self.integ.step(3000)

    def poser_liens(self, couples, frac=1.0):
        """Bascule la raideur des liens d'extrusion actifs a cet instant."""
        if not self.map_lef:
            return
        vise = set()
        for i, j in couples:
            cle = (int(min(i, j)), int(max(i, j)))
            if cle in self.map_lef and cle[0] != cle[1]:
                vise.add(cle)
        for cle in self.actifs - vise:
            b = self.map_lef[cle]
            self.liens.setBondParameters(b, cle[0], cle[1], 1.0, 0.0)
        for cle in vise:
            b = self.map_lef[cle]
            self.liens.setBondParameters(b, cle[0], cle[1], 1.0,
                                         frac * self.k_lef)
        self.liens.updateParametersInContext(self.ctx)
        self.actifs = vise

    def bloc(self, n_pas, couples=None):
        if couples is not None and self.map_lef:
            self.poser_liens(couples)
        self.integ.step(int(n_pas))
        q = self.ctx.getState(getPositions=True).getPositions(asNumpy=True)
        return np.asarray(q.value_in_unit(unit.nanometer))


# ======================================================================
def simuler(n, n_blocs, pas_par_bloc, traj_lef=None, liants=None, poids=None,
            graine=0, **kw):
    """Une trajectoire complete. Renvoie (n_blocs, n, 3)."""
    paires = None
    if traj_lef is not None:
        paires = {(int(min(a, b)), int(max(a, b)))
                  for inst in traj_lef for a, b in inst}
    pol = Polymere(n, paires_lef=paires, liants=liants, poids=poids,
                   graine=graine, **kw)
    init = [tuple(x) for x in traj_lef[0]] if traj_lef is not None else None
    pol.demarrer(liens_initiaux=init)
    sorties = []
    for b in range(n_blocs):
        couples = ([tuple(x) for x in traj_lef[b % len(traj_lef)]]
                   if traj_lef is not None else None)
        sorties.append(pol.bloc(pas_par_bloc, couples))
    return np.array(sorties), pol.plateforme
