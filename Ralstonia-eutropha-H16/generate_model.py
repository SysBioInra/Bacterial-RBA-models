"""Semi-automatic generation of Ralstonia eutropha RBA model"""

# python 2 compatibility
from __future__ import absolute_import, division, print_function

# package imports
import rba
import cobra


# MAIN FUNCTION --------------------------------------------------------
#
# model creation using files in data/:
def main():
    
    # make some inital modifications to sbml required for RBA
    #import_sbml_model("../../genome-scale-models/Ralstonia_eutropha/sbml/RehMBEL1391_sbml_L3V1.xml")
    
    # inital run of model generation creates helper files
    reutropha = rba.RbaModel.from_data('params.in')
    
    # set a growth medium
    reutropha.set_medium('data/medium.tsv')
    
    # set k_app default efficiencies
    set_default_efficiencies(reutropha)
    
    # set k_app for selected reactions (from GECKO, optional)
    #reutropha.set_enzyme_efficiencies('data/enzyme_efficiency.tsv')
    
    # set total protein constraints for cell compartments
    set_compartment_params(reutropha)
    
    # export to files
    reutropha.write()


# the following function imports the genome scale model from remote
# location and implements some changes important for making RBA model
# import model
def import_sbml_model(model_path):
    
    model = cobra.io.read_sbml_model(model_path)
    
    # add (dummy) tRNA loading reaction for Asparagin 
    # in the original model, tRNA-Asn is made from tRNA-Asp directly
    # Reactions PRUK and RBPC (Rubisco enzyme) are added
    model.add_metabolites(
    cobra.Metabolite(
        id = 'trnaasn_c',
        name = 'tRNA(Asp)',
        compartment = 'c',
        charge = 0,
        formula = 'C10H17O10PR2'))
    
    
    ASNTRS = cobra.Reaction('ASNTRS')
    ASNTRS.name = 'Asparaginyl-tRNA synthetase'
    ASNTRS_string = 'asn__L_c + atp_c + trnaasn_c --> amp_c + asntrna_c + h_c + ppi_c'
    model.add_reactions([ASNTRS])
    model.reactions.ASNTRS.build_reaction_from_string(ASNTRS_string)
    model.reactions.ASNTRS.gene_reaction_rule = 'H16_A0453'
    
    # replace all empty gene associations with UNKNOWN in order to
    # avoid spontaneous reactions (see RBApy manual)
    for r in model.reactions:
        if r.gene_reaction_rule == '':
            r.gene_reaction_rule = 'UNKNOWN'
    
    # manual curation of reactions that take part in artificial cycles
    model.reactions.FDH.bounds = (0.0, 1000.0)
    model.remove_reactions(['FRUpts2', 'ACt2r', 'ACACt2', 'MDH2'])
    
    # export model as sbml.xml
    cobra.io.write_sbml_model(model, 'data/sbml.xml')


# set k_app default efficiencies analogous to E.coli
def set_default_efficiencies(model):
    
    # different options for default enzyme kcat:
    # 12.5 * 3600 = 45000/h; Bulovic et al., 2019, RBApy
    # 65 * 3600 = 234000/h; Lloyd et al., 2018, CobraME
    # 172 * 3600 = 619200/h; Salvy et al., 2020, ETFL
    
    # systematically testing different default values revealed
    # that transporter efficiency has little effect on µ compared to
    # non transporter enzymes
    fn = model.parameters.functions.get_by_id('default_efficiency')
    fn.parameters.get_by_id('CONSTANT').value = 45000
    
    fn = model.parameters.functions.get_by_id('default_transporter_efficiency')
    fn.parameters.get_by_id('CONSTANT').value = 360000


# set total protein constraints for cell compartments
def set_compartment_params(model):
    
    # THE FOLLOWING VALUES WERE EXPERIMENTALLY DETERMINED USING MS
    # For details, see R notebook 'Ralstonia eutropha model constraints'
    # at https://github.com/m-jahn/R-notebooks
    # ------------------------------------------------------------------
    #
    # remove old parameter constraints from model
    pars = ([
        'fraction_protein_Cytoplasm', 
        'fraction_protein_Cell_membrane', 
        'fraction_non_enzymatic_protein_Cytoplasm', 
        'fraction_non_enzymatic_protein_Cell_membrane'])
    
    for par in pars:
        fn = model.parameters.functions.get_by_id(par)
        model.parameters.functions.remove(fn)
    
    # set secreted proteins to zero
    fn = model.parameters.functions.get_by_id('fraction_protein_Secreted')
    fn.parameters.get_by_id('CONSTANT').value = 0
    
    # set fraction of cytoplasmic proteins
    par_frac_cp = {'LINEAR_COEF': 0.178, 'LINEAR_CONSTANT': 0.848, 'X_MIN':0, 'X_MAX':2, 'Y_MIN':0, 'Y_MAX':1}
    model.parameters.functions.append(
        rba.xml.Function('fraction_protein_Cytoplasm', 'linear', par_frac_cp)
    )
    
    # set fraction of membrane proteins
    par_frac_mp = {'LINEAR_COEF': -0.178, 'LINEAR_CONSTANT': 0.152, 'X_MIN':0, 'X_MAX':2, 'Y_MIN':0, 'Y_MAX':1}
    model.parameters.functions.append(
        rba.xml.Function('fraction_protein_Cell_membrane', 'linear', par_frac_mp)
    )
    
    # set non-enzymatic fraction of protein for Cytoplasm
    par_ne_cp = {'LINEAR_COEF': -0.613, 'LINEAR_CONSTANT': 0.610, 'X_MIN':0, 'X_MAX':2, 'Y_MIN':0, 'Y_MAX':1}
    model.parameters.functions.append(
        rba.xml.Function('fraction_non_enzymatic_protein_Cytoplasm', 'linear', par_ne_cp)
    )
    
    # set non-enzymatic fraction of protein for Cell membrane
    par_ne_mp = {'LINEAR_COEF': -0.204, 'LINEAR_CONSTANT': 0.908, 'X_MIN':0, 'X_MAX':2, 'Y_MIN':0, 'Y_MAX':1}
    model.parameters.functions.append(
        rba.xml.Function('fraction_non_enzymatic_protein_Cell_membrane', 'linear', par_ne_mp)
    )


if __name__ == "__main__":
    main()
