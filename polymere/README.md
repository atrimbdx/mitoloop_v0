# polymere/ - orthogonal check of AlphaGenome 2 by polymer simulation

`simpol/` is the simulation package of the TYRP1 project (path A: loop
extrusion stopped by oriented CTCF sites, then an OpenMM coarse-grained
polymer; path B: binders fitted on the contact map). It is copied unchanged,
except `simpol/config.py`, which also reads loci described in a JSON file
named by `SIMPOL_LOCI`. See `simpol/README.md` for the model.

Only path A is used here. Its single input is the position and orientation of
the CTCF sites: it never sees the Hi-C map nor the rE2G elements, so a contact
it predicts, or loses after a deletion, is independent of the data that
called the loop. Path B is fitted on the map and would be circular.

| file | role |
|---|---|
| `preparer_locus.py preparer` | window, anchors and hub of one gene (from `ag2_balayage.py locus` cibles.tsv and the tissue's loops), measured matrix (choice of the parameter set only), oriented CTCF sites, and the deletion conditions: SANSELEMENT (sites inside the elements on anchor), SANSANCRE<kb> (one per distal loop summit: sites within 10 kb of it); a deletion that removes no site is null by construction |
| `preparer_locus.py plan` | plan of the deletion runs, with the parameter set that best matches the measured virtual 4C (simpol comparer.py) |
| `trajectoires.sbatch` | one trajectory per job-array task |
| `contact_deletion.py` | promoter-anchor contact frequency per condition, log2 change, bootstrap interval, null from splits of the wild-type replicates |

Reading with AlphaGenome 2: a deletion that lowers expression (AG2) while
CTCF-only extrusion keeps the contact points to an activity element; a
deletion that also loses the contact in path A points to an architectural one.
The CTCF sites come from ENCODE, not from the patients' cells: a proxy.
