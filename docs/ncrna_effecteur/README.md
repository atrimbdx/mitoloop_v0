# Chantier « non codant effecteur » et phénotypes mitochondriaux

Ce chantier a démarré dans une autre session, le 19 novembre 2026. Il se poursuit dans une discussion à part.

## La question

Quels RNU, miRNA, lncRNA ou autres transcrits non traduits peuvent donner un phénotype mitochondrial ?
Comment les chercher chez les patients non résolus après WGS ?

## Fichiers

- `brief_ncrna.md` : la question, le contexte, ce qui existe déjà dans reannot (`petits_arn`) et les
  contraintes.
- `strategie_ncrna.md` : le rapport, en version 2 après contre-revue. Il contient le catalogue, la stratégie,
  les 13 décisions D0 à D12 et les références.
- `catalogue_ncrna_mito.tsv` : 329 entrées.
  - Statuts : `confirme` (141), `conteste` (4), `refute` (1), `non_verifie` (12),
    `non_verifie_priorite_basse` (26) et `contexte_non_verifie` (145, gènes protéiques et expansions).
- `introns_U12_genes_mito_hg38.tsv` : gènes mitochondriaux à intron mineur (U12). Les symboles HGNC sont à
  mettre à jour, par exemple GARS devenu GARS1.

## Ce qui reste à faire

1. **Verdict des 12 gènes restés sans verdict** : CHASERR, SNHG14, CARMN, Maenli, MDRL, MALAT1, MIR23A,
   MIR29B1, MIR195, MIR15B, MIR141, MIR30A.
2. **Contre-revue inachevée**, interrompue par la limite d'usage :
   - Les 7 corrections de fond de la v2 n'ont pas été contre-vérifiées. Elles portent sur le barème des cas
     types (§ e.5), le PP4 de l'ADNmt, MT-TM et les régions récessives de RNU2-2.
   - Deux sections n'ont pas été relues : « étapes » (§ e.6 à e.9) et « synthèse, risques, décisions ».
   - Les 33 remarques mineures ne sont pas appliquées.
3. **Points « à relire »**, listés à la fin de la section (h) du rapport.
4. **Décisions de l'utilisateur**, D0 à D12. Les plus urgentes : D0 (CRAM et conteneurs sur la VM), D1, D2
   et D12.

## Contraintes (rappel)

- Aucun entraînement ni ajustement de seuil sur des sorties d'AlphaGenome, SpliceAI, SpliceAI2 ou
  AlphaFold 3.
- Aucune donnée patient vers un service externe.
- OMIM jamais redistribué.
- Fichiers bruts d'éditeurs (données supplémentaires d'articles, dans `reannot/donnees/`) et tables qui en dérivent : usage interne seulement ; avant toute diffusion d'un dépôt hors de l'équipe, vérifier leur licence ou les retirer du dépôt et de son historique (règle détaillée : `reannot/donnees/LISEZMOI_editeurs.md`).
- Pour les travaux en agents multiples : effort moyen, et peu d'agents, à cause des limites d'usage.
