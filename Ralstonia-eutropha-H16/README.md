*Ralstonia eutropha H16* RBA model
================
Michael Jahn,
2020-04-28

***

This resource balance analysis (RBA) model was generated using the [RBApy package from Bulovic et al.](https://doi.org/10.1016/j.ymben.2019.06.001),
using a manually curated metabolic reconstruction originally developed by [Park et al., 2011](http://bmcsystbiol.biomedcentral.com/articles/10.1186/1752-0509-5-101). The original model was refactored and many reactions were curated and harmonized to community standards as e.g. available through the [Bigg model database](bigg.ucsd.edu/). The genome scale model that serves as base for the RBA model is maintained on [github.com/m-jahn/](https://github.com/m-jahn/genome-scale-models).

Sources for parameter estimation are:

 - [Goelzer at al., 2015](https://doi.org/10.1016/j.ymben.2015.10.003): set of parameters for transcription, translation, and other fundamental cellular functions.
 - [Park et al., 2011](http://bmcsystbiol.biomedcentral.com/articles/10.1186/1752-0509-5-101) for the cell composition of *R. eutropha* (original biomass equation)


## Generation of the XML model

The *Ralstonia eutropha H16* RBA model was generated in a step-wise manner,
based on information from a hand-curated model.


### Step 1: Generation of a default model

For initial creation of the RBA model, the following input files were used. The fasta file headers have to be formatted according to this pattern: `rba|id|name|set_name|stoichiometry`.

 - `data/sbml.xml`: metabolic netwok from the hand-curated model. Some minor modifications were made to make it fit before RBA model construction (see `generate_model.py`). For example, empty gene associations were replaced with `UNKNOWN` to prevent RBA reading those reactions as spontaneous.
 - `data/ribosome.fasta`: composition of ribosome (rRNA and proteins). rRNA sequences were retrieved from NCBI, ribosomal protein sequence and subunit structure from uniprot. 68 ribosomal proteins were included. Ribosomal stoichiometry was assumed to be one subunit each for 30S and 50S subunits, except for L7/L12 that forms 2-4 dimers per ribosome (`n >= 4`). *TODO: update stochiometry for elongation and initiation factors.*
 - `data/chaperones.fasta`: composition of molecular machines involved in protein folding. These include the protein subunits `groES`, `groEL`, `dnaJ`, `dnaK`, `grpE`, `trigger factor tig`, `htpG`, and `hfq` with their respective number od subunits obtained from uniprot.
 - `data/trnas.fasta`: composition of tRNAs. tRNA sequences were retrieved from NCBI, and only one representative tRNA (the first) per amino acid was used (20 out of 61 tRNAs).
 - `params.in`: pipeline parameters.

The initial RBApy model generation was performed via the following customized function. The function automatically applies custom changes to the SBML model, imports the `params.in` file with basic settings, and generates the model files.

```
python3 generate_model.py
```

### Step 2: Generation of a model improved through helper files

This first run downloaded `data/uniprot.csv`, generated helper files `data/*.tsv`, and generated a first XML model. The following modifications were then made to helper files.

 - `location_map.tsv`: the uniprot compartment names were mapped to generic names (such as `Cytoplasm`)
 - `metabolites.tsv`: this table contains target metabolite concentrations for basic biomolecules and cofactors (for example `ATP`, `dATP`, `thiamin`, and so on). Molecules that are already contained in macrocomponents were set to zero. One exception are some nucleotide species that are present in more detail compared to macrocomponent reactions.
 - `macrocomponents.tsv`: This table was supplemented with all major cell components (such as lipopolysaccharides, peptidoglycan, cofactors, and so on) and corresponding target concentrations from the original model's biomass equation. Target fluxes to other metabolites can be defined here as well, such as flux to PHB.
 - `medium.tsv`: added new minimal medium.
 - `subunits.tsv`: updated some protein stoichiometries (ATP synthase).
 - `unknown_proteins.tsv`: (empty) the SBML model contains no genes that could not be mapped

Two additional macroprocesses were added, **transcription of DNA to mRNA**, and **DNA replication**. 
The RBApy source code was adapted to automatically load protein associations from `transcription.fasta` and `replication.fasta` (see [github fork extra processes](https://github.com/m-jahn/RBApy/tree/extra_processes)). Adding machinery for transcription (TSC) was based on the following assumptions:

 - We need **a capacity constraint**, effectively the transcription rate per RNA polymerase unit that is added to parameters. On bionumbers, several values in the range of 20-80 nt per second per unit RNA polymerase are given. We take the most recent estimate for *E. coli*, 62 nt/s or 223200 nt/h (source: Epshtein V, Nudler E. Cooperation between RNA polymerase molecules in transcription elongation. Science. 2003 May 2 300(5620):801-5.). PubMed ID12730602.
- We need **proteins that catalyze transcription** (RNA polymerase complex). The canonical RNA polymerase in bacteria consists of 2 alpha, 1 beta, 1 beta prime, and 1 omega subunit. When a sigma factor is associated with the core, the holoenzyme is formed, which can initiate transcription. Four proteins were added as important sigma factors, rpoD1, D2, rpoS, and rpoN, each with partial contribution of 0.25 (source: uniprot.org, subunit structure).

Adding machinery for replication (REP) was based on the following assumptions:

- DNA replication is mainly catalyzed by DNA-dependent DNA polymerase III. We need to add a rate for this enzyme. Bionumbers gives estimates from 600-1000 nt/s, so on average 800 nt/s or 2880000 nt/h.
- DNA polymerase III contains a core (composed of alpha, epsilon and theta chains) that associates with a tau subunit. This core dimerizes to form the POLIII complex. PolIII associates with the gamma complex (composed of gamma, delta, delta prime, psi and chi chains) and with the beta chain to form the complete DNA polymerase III complex.
Additionally DNA helicase (hexamer), DNA (primase monomer), and DNA beta clamp (dimer) was added. The theta subunit (holE) is missing (`'Q0KD75': 2, 'Q0K8W5': 2, 'Q0K937': 3, 'Q0K709': 1, 'Q0KBB9': 1, 'Q0K7F6': 1, 'Q0K9E9': 6, 'Q0K866': 1, 'Q0KFR7': 2`). Source: uniprot.org, subunit structure.

### Step 3: Generate and solve full model

In the last step, the full model can be generated by running one more iteration of:
```
python3 generate_model.py
```

The model can then be solved by running:

```
python3 solve_rba_model.py .
```

### Step 4: Calibration

The model can be calibrated with proteomics data, fluxomics data, or estimations of the `k_app` parameter. `k_app` is a collective term that links fluxes to enzyme abundances. It's theoretical maximum is `k_cat`, the turnover number of an enzyme, but it also depends on the metabolite concentration (therefore saturation), and reaction equilibirium.

**k_app**

`k_app` values e.g. obtained from BRENDA can optionally be included using the following command in `generate_model.py`:    `reutropha.set_enzyme_efficiencies('data/enzyme_efficiency.tsv')`

**Protein fraction per compartment**
<!--
Protein fraction per compartment was calculated using the 'RBA estim' functions. Two tables, `calibration/experiment.csv` and `calibration/protein.csv` were prepaeed according to the RBApy manual. The estimation procedure is performed using the following script in the `rba/estim` folder (adapt file path):

```
cd path/to/RBApy/rba/estim
python3 prot_per_compartment.py
```
-->

Protein fraction per compartment was calculated as described in the R notebook [Ralstonia model constraints](https://github.com/m-jahn/R-notebooks).
