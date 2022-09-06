"""Semi-automatic generation of Ralstonia eutropha RBA model"""

# python 2 compatibility
from __future__ import absolute_import, division, print_function

# package imports
import re
import rba
import cobra


# MAIN FUNCTION --------------------------------------------------------
#
# model creation using files in data/:
def main():
    
    # make some inital modifications to sbml required for RBA
    #import_sbml_model("../../genome-scale-models/Ralstonia_eutropha/sbml/RehMBEL1391_sbml_L3V1_maint.xml")
    
    # inital run of model generation creates helper files
    reutropha = rba.RbaModel.from_data('params.in')
    
    # set a growth medium
    reutropha.set_medium('data/medium.tsv')
    
    # add replication and transcription machinery
    update_processes(reutropha)
    
    # set k_app default efficiencies
    set_default_efficiencies(reutropha)
    
    # set k_app for selected reactions
    reutropha.set_enzyme_efficiencies('calibration/kapp_consensus.csv')
    
    # set maintenance demand
    set_maintenance(reutropha)
    
    # set total protein constraints for cell compartments
    set_compartment_params(reutropha)
    
    # export to files
    reutropha.write()



# the following function imports the genome scale model from remote
# location and implements some changes important for making RBA model
def import_sbml_model(model_path):
    
    # import model from external dir
    model = cobra.io.read_sbml_model(model_path)
    
    # replace all empty gene associations with UNKNOWN in order to
    # avoid spontaneous reactions (see RBApy manual)
    # An exception are gases that diffuse more or less freely.
    for r in model.reactions:
        if r.gene_reaction_rule == '':
            if r.id in (['H2Ot', 'CO2t', 'O2t']):
                r.gene_reaction_rule = ''
            elif re.search('t[0-9]?r?$|abc$|pts$', r.id):
                r.gene_reaction_rule = 'UNKNOWN_TRANSPORTER'
            else:
                r.gene_reaction_rule = 'UNKNOWN'
    
    # remove original, superfluous maintenance reaction
    model.remove_reactions(['Maintenance'])
    
    # export model as sbml.xml
    cobra.io.write_sbml_model(model, 'data/sbml.xml')


# customize some process parameters
def update_processes(model):
       
    # Adjust target (steady-state) concentrations for DNA and RNA from 
    # RBApy default to the ones from R.e. biomass equation (Park et al., 2011)
    fn = model.parameters.functions.get_by_id('mrna_concentration')
    fn.parameters.get_by_id('CONSTANT').value = 0.06
    
    fn = model.parameters.functions.get_by_id('dna_concentration')
    fn.parameters.get_by_id('CONSTANT').value = 0.031
    
    # Remove test processes from default model that are not neccessary
    # for pr in ['test_process_0', 'test_process_1', 'test_process_2']:
        # pr = model.processes.processes.get_by_id(pr)
        # model.processes.processes.remove(pr)


# set k_app default efficiencies analogous to E.coli
def set_default_efficiencies(model):
    
    # different options for default enzyme kcat:
    # 12.5 * 3600 = 45000/h; Bulovic et al., 2019, RBApy
    # 65 * 3600 = 234000/h; Lloyd et al., 2018, CobraME
    # 172 * 3600 = 619200/h; Salvy et al., 2020, ETFL
    # or median of k_app parameter estimation
    
    fn = model.parameters.functions.get_by_id('default_efficiency')
    fn.parameters.get_by_id('CONSTANT').value = 5770
    
    fn = model.parameters.functions.get_by_id('default_transporter_efficiency')
    fn.parameters.get_by_id('CONSTANT').value = 5770


# set maintenance ATP consumption. Maintenance can consist of a growth-
# dependent (constant) and independent term (coef)
def set_maintenance(model):
    
    # In the original publication for the R. eutropha genome scale model,
    # the growth- and non growth-associated maintenance in minimal medium were
    # determined as 15.30 mmol gDCW-1 and 3.00 mmol gDCW-1 h-1, respectively. 
    # A constant NGAM leads to higher yield at higher growth rate.
    
    fn = model.parameters.functions.get_by_id('maintenance_atp')
    fn.parameters.get_by_id('LINEAR_CONSTANT').value = 3
    fn.parameters.get_by_id('LINEAR_COEF').value = 150
    fn.parameters.get_by_id('X_MIN').value = 0
    fn.parameters.get_by_id('X_MAX').value = 1
    
    # Set flux to PHB depending on NH4+ concentration
    # first construct target from input reaction ID and flux boundary
    boundary_id = 'R_PHBt_flux_boundary'
    phb_target = rba.xml.targets.TargetReaction('R_PHBt')
    phb_target.value = boundary_id
    # add target to model
    mp = model.targets.target_groups.get_by_id('metabolite_production')
    mp.reaction_fluxes.append(phb_target)
    # add linear relationship to model
    par_phb = {'LINEAR_COEF': -1.0, 'LINEAR_CONSTANT': 3.25, 'X_MIN':0, 'X_MAX':3.25, 'Y_MIN':0, 'Y_MAX':3.25}
    phb = rba.xml.Function(boundary_id, 'linear', par_phb)
    phb.variable = 'M_nh4_e'
    model.parameters.functions.append(phb)


# set total protein constraints for cell compartments
def set_compartment_params(model):
    
    # The values used for compartments were determined using MS-based proteomics.
    # For details, see R notebook 'Ralstonia eutropha model constraints'
    # at https://github.com/m-jahn/R-notebooks
    # ------------------------------------------------------------------
    #
    # remove old parameter constraints from model
    pars = ([
        'Cytoplasm_density',
        'fraction_protein_Cytoplasm', 
        'fraction_protein_Cell_membrane', 
        'fraction_non_enzymatic_protein_Cytoplasm', 
        'fraction_non_enzymatic_protein_Cell_membrane'])
    
    for par in pars:
        fn = model.parameters.functions.get_by_id(par)
        model.parameters.functions.remove(fn)
    
    # set amino_acid_concentration to constant, growth rate independent value.
    # The default value of is 4.8972 mmol/gDCW. Assuming an average molecular weight
    # per amino acid of of 110 mg / mmol, total protein mass = 4.8972*110/1000 = 0.5387 g/gDCW.
    # This is a bit lower than the reported estimates for R. eutropha of 0.68 g/gDCW, or 
    # cyanobacteria with 0.65 g/gDCW. A constant protein pool of 0.68g/gDCW would
    # translate to an average amino acid concentration = 6.1818 mmol/gDCW
    fn = model.parameters.functions.get_by_id('amino_acid_concentration')
    fn.parameters.get_by_id('LINEAR_CONSTANT').value = 6.1818
    fn.parameters.get_by_id('LINEAR_COEF').value = 0
    fn.parameters.get_by_id('X_MIN').value = 0
    fn.parameters.get_by_id('X_MAX').value = 1
    
    # set cytoplasm density to a fraction of total amino acid 
    # concentration analogously to Cell_membrane density
    Cytoplasm_density = rba.xml.Aggregate('Cytoplasm_density', type_='multiplication')
    Cytoplasm_density.function_references.append(rba.xml.FunctionReference('amino_acid_concentration'))
    Cytoplasm_density.function_references.append(rba.xml.FunctionReference('fraction_protein_Cytoplasm'))
    model.parameters.aggregates.append(Cytoplasm_density)
    
    # set secreted proteins to zero
    fn = model.parameters.functions.get_by_id('fraction_protein_Secreted')
    fn.parameters.get_by_id('CONSTANT').value = 0
    
    # set fraction of cytoplasmic proteins
    par_frac_cp = {'LINEAR_COEF': 0.1060, 'LINEAR_CONSTANT': 0.8684, 'X_MIN':0, 'X_MAX':2, 'Y_MIN':0, 'Y_MAX':1}
    model.parameters.functions.append(
        rba.xml.Function('fraction_protein_Cytoplasm', 'linear', par_frac_cp)
    )
    
    # set fraction of membrane proteins
    par_frac_mp = {'LINEAR_COEF': -0.1060, 'LINEAR_CONSTANT': 0.1316, 'X_MIN':0, 'X_MAX':2, 'Y_MIN':0, 'Y_MAX':1}
    model.parameters.functions.append(
        rba.xml.Function('fraction_protein_Cell_membrane', 'linear', par_frac_mp)
    )
    
    # set non-enzymatic fraction of protein for Cytoplasm
    par_ne_cp = {'LINEAR_COEF': -0.5657, 'LINEAR_CONSTANT': 0.5374, 'X_MIN':0, 'X_MAX':2, 'Y_MIN':0, 'Y_MAX':1}
    model.parameters.functions.append(
        rba.xml.Function('fraction_non_enzymatic_protein_Cytoplasm', 'linear', par_ne_cp)
    )
    
    # set non-enzymatic fraction of protein for Cell membrane
    par_ne_mp = {'LINEAR_COEF': -0.2144, 'LINEAR_CONSTANT': 0.8415, 'X_MIN':0, 'X_MAX':2, 'Y_MIN':0, 'Y_MAX':1}
    model.parameters.functions.append(
        rba.xml.Function('fraction_non_enzymatic_protein_Cell_membrane', 'linear', par_ne_mp)
    )


if __name__ == "__main__":
    main()
