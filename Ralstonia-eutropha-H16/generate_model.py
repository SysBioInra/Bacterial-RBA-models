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
    import_sbml_model("../../genome-scale-models/Ralstonia_eutropha/sbml/RehMBEL1391_sbml_L3V1.xml")
    
    # inital run of model generation creates helper files
    reutropha = rba.RbaModel.from_data('params.in')
    
    # set a growth medium
    reutropha.set_medium('data/medium.tsv')
    
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
    model.add_reactions([ASNTRS])
    model.reactions.ASNTRS.build_reaction_from_string('asn__L_c + atp_c + trnaasn_c --> amp_c + asntrna_c + h_c + ppi_c')
    model.reactions.ASNTRS.gene_reaction_rule = 'H16_A0453'
    
    # replace all empty gene associations with UNKNOWN in order to
    # avoid spontaneus reactions (see RBApy manual)
    for r in model.reactions:
        if r.gene_reaction_rule == '':
            r.gene_reaction_rule = 'UNKNOWN'
    
    # export model as sbml.xml
    cobra.io.write_sbml_model(model, 'data/sbml.xml')


if __name__ == "__main__":
    main()
