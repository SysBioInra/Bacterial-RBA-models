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
    # inital run of model generation creates helper files
    vnat_rba = rba.RbaModel.from_data('params.in')

    # set a growth medium
    vnat_rba.set_medium('data/medium.tsv')
    
    # add replication and transcription machinery
    update_processes(vnat_rba)
    
    # set k_app default efficiencies
    set_default_efficiencies(vnat_rba)
    
    # set k_app for selected reactions
    vnat_rba.set_enzyme_efficiencies('calibration_data/kapp_estimate_reformatted.csv')

    # set maintenance demand
    set_maintenance(vnat_rba)
    
    # set total protein constraints for cell compartments
    # set_compartment_params(vnat_rba)
    # set_compartment_params_simple(vnat_rba)
    set_compartment_params(vnat_rba)

    # export to files
    # vnat_rba.solve()
    vnat_rba.write()

    make_composition_file(vnat_rba)

# customize some process parameters
def update_processes(model):
       
    # Adjust target (steady-state) concentrations for DNA and RNA from 
    # RBApy default to the ones from R.e. biomass equation (Park et al., 2011)
    fn = model.parameters.functions.get_by_id('mrna_concentration')
    fn.parameters.get_by_id('CONSTANT').value = 0.06
    
    fn = model.parameters.functions.get_by_id('dna_concentration')
    fn.parameters.get_by_id('CONSTANT').value = 0.031
    
    # Remove test processes from default model that are not neccessary
    for pr in ['test_process_0', 'test_process_1', 'test_process_2']:
        pr = model.processes.processes.get_by_id(pr)
        model.processes.processes.remove(pr)

# set k_app default efficiencies analogous to E.coli
def set_default_efficiencies(model):
    
    # different options for default enzyme kcat:
    # 12.5 * 3600 = 45000/h; Bulovic et al., 2019, RBApy
    # 65 * 3600 = 234000/h; Lloyd et al., 2018, CobraME
    # 172 * 3600 = 619200/h; Salvy et al., 2020, ETFL
    # or median of k_app parameter estimation
    
    fn = model.parameters.functions.get_by_id('default_efficiency')
    fn.parameters.get_by_id('CONSTANT').value = 619200
    
    fn = model.parameters.functions.get_by_id('default_transporter_efficiency')
    fn.parameters.get_by_id('CONSTANT').value = 619200

# set maintenance ATP consumption. Maintenance can consist of a growth-
# dependent (constant) and independent term (coef)
def set_maintenance(model):
    # ngam mmol gDCW-1
    # gam mmol gDCW-1 h-1
    fn = model.parameters.functions.get_by_id('maintenance_atp')
    fn.parameters.get_by_id('LINEAR_CONSTANT').value = 50
    fn.parameters.get_by_id('LINEAR_COEF').value = 0
    fn.parameters.get_by_id('X_MIN').value = 0
    fn.parameters.get_by_id('X_MAX').value = 1

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
    fn.parameters.get_by_id('LINEAR_CONSTANT').value = 4.2727
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
    par_frac_cp = {'LINEAR_COEF': 0.0321, 'LINEAR_CONSTANT': 0.8064, 'X_MIN':0, 'X_MAX':2, 'Y_MIN':0, 'Y_MAX':1}
    model.parameters.functions.append(
        rba.xml.Function('fraction_protein_Cytoplasm', 'linear', par_frac_cp)
    )
    
    # set fraction of membrane proteins
    par_frac_mp = {'LINEAR_COEF': -0.0321, 'LINEAR_CONSTANT': 0.1935, 'X_MIN':0, 'X_MAX':2, 'Y_MIN':0, 'Y_MAX':1}    
    model.parameters.functions.append(
        rba.xml.Function('fraction_protein_Cell_membrane', 'linear', par_frac_mp)
    )
    
    # set non-enzymatic fraction of protein for Cytoplasm
    par_ne_cp = {'LINEAR_COEF': -0.0517, 'LINEAR_CONSTANT': 0.3066, 'X_MIN':0, 'X_MAX':2, 'Y_MIN':0, 'Y_MAX':1}
    model.parameters.functions.append(
        rba.xml.Function('fraction_non_enzymatic_protein_Cytoplasm', 'linear', par_ne_cp)
    )
    
    # set non-enzymatic fraction of protein for Cell membrane
    par_ne_mp = {'LINEAR_COEF': -0.1120, 'LINEAR_CONSTANT': 0.7761, 'X_MIN':0, 'X_MAX':2, 'Y_MIN':0, 'Y_MAX':1}
    model.parameters.functions.append(
        rba.xml.Function('fraction_non_enzymatic_protein_Cell_membrane', 'linear', par_ne_mp)
    )

def set_compartment_params_simple(model):
    
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
    fn.parameters.get_by_id('LINEAR_CONSTANT').value = 4.2727
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
    par_frac_cp = {'LINEAR_COEF': 0.0331, 'LINEAR_CONSTANT': 0.8496, 'X_MIN':0, 'X_MAX':2, 'Y_MIN':0, 'Y_MAX':1}
    model.parameters.functions.append(
        rba.xml.Function('fraction_protein_Cytoplasm', 'linear', par_frac_cp)
    )
    
    # set fraction of membrane proteins
    par_frac_mp = {'LINEAR_COEF': -0.0331, 'LINEAR_CONSTANT': 0.1503, 'X_MIN':0, 'X_MAX':2, 'Y_MIN':0, 'Y_MAX':1}
    model.parameters.functions.append(
        rba.xml.Function('fraction_protein_Cell_membrane', 'linear', par_frac_mp)
    )
    
    # set non-enzymatic fraction of protein for Cytoplasm
    par_ne_cp = {'LINEAR_COEF': -0.0666, 'LINEAR_CONSTANT': 0.2530, 'X_MIN':0, 'X_MAX':2, 'Y_MIN':0, 'Y_MAX':1}
    model.parameters.functions.append(
        rba.xml.Function('fraction_non_enzymatic_protein_Cytoplasm', 'linear', par_ne_cp)
    )
    
    # set non-enzymatic fraction of protein for Cell membrane
    par_ne_mp = {'LINEAR_COEF': -0.1635, 'LINEAR_CONSTANT': 0.6766, 'X_MIN':0, 'X_MAX':2, 'Y_MIN':0, 'Y_MAX':1}
    model.parameters.functions.append(
        rba.xml.Function('fraction_non_enzymatic_protein_Cell_membrane', 'linear', par_ne_mp)
    )

def make_composition_file(model):
    out_file_handle = open('model/compositions_enzymes.tsv', 'w')

    for e in model.enzymes.enzymes:
        out_file_handle.write(e.id)
        for r in e.machinery_composition.reactants:
            out_file_handle.write("\t"+str(int(r.stoichiometry))+"x"+r.species)
        out_file_handle.write("\n")

    out_file_handle = open('model/compositions_machineries.tsv', 'w')

    for p in model.processes.processes:
        out_file_handle.write(p.id)
        for r in p.machinery.machinery_composition.reactants:
            if r.species.startswith("PN"):
                out_file_handle.write("\t"+str(int(r.stoichiometry))+"x"+r.species)
        out_file_handle.write("\n")



if __name__ == "__main__":
    main()
