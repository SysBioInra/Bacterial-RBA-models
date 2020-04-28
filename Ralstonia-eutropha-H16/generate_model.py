"""Semi-automatic generation of Ralstonia eutropha RBA model"""

# python 2 compatibility
from __future__ import absolute_import, division, print_function

# package imports
import rba
import cobra


# the following function imports the genome scale model from remote
# location and implements some changes important for making RBA model
# import model
def modify_model():
    
    model_path = "../../genome-scale-models/Ralstonia_eutropha/sbml/RehMBEL1391_sbml_L3V1.xml"
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
    model.add_reactions([ASNTRS])
    model.reactions.ASNTRS.build_reaction_from_string('asn__L_c + atp_c + trnaasn_c --> amp_c + asntrna_c + h_c + ppi_c')
    model.reactions.ASNTRS.gene_reaction_rule = 'H16_A0453'
    
    # export model as sbml.xml
    cobra.io.write_sbml_model(model, 'data/sbml.xml')


def export_minimal_medium():
    
    # Ralstonia eutropha minimal medium with the
    # following core components, not including carbon
    # or nitrogen source
    minimal_medium = {
        'EX_mg2_e': 10.0,
        'EX_pi_e': 1000.0,
        'EX_cobalt2_e': 10.0,
        'EX_cl_e': 10.0,
        'EX_k_e': 10.0,
        'EX_fe3_e': 10.0,
        'EX_so4_e': 10.0,
        'EX_na_e': 10.0,
        'EX_o2_e': 18.5,
        'EX_mobd_e': 10.0,
        'EX_h2o_e': 1000.0,
        'EX_h_e': 100.0
        }
    return(minimal_medium)


# MAIN FUNCTION --------------------------------------------------------
#
# input for initial model creation are the following files in data/:
# - genome scale model in SBML format
# - ribosome.fasta
# - chaperones.fasta
# - trnas.fasta
def main():
    
    # make some inital modifications to sbml required for RBA
    #modify_model()
    
    # inital run of model generation creates helper files
    reutropha = rba.RbaModel.from_data('params.in')
    
    # set a growth medium
    reutropha.set_medium('data/medium.tsv')
    
    # export to files
    reutropha.write()



if __name__ == "__main__":
    main()
