"""Module generating subtilis model."""

# python 2 compatibility
from __future__ import absolute_import, division, print_function

# imports
import rba


def main():
    subtilis = rba.RbaModel.from_data('params.in')
    subtilis.set_medium('data/curated_medium.tsv')
    subtilis.set_enzyme_efficiencies('data/catalytic_activity_medium_2.csv')
    add_flagella_constraint(subtilis)
    subtilis.write()


def add_flagella_constraint(subtilis):
    subtilis.targets.target_groups.append(flagella_activation())
    for fn in flagella_activation_functions():
        subtilis.parameters.functions.append(fn)
    subtilis.parameters.aggregates.append(flagella_activation_aggregate())


def flagella_activation():
    target_group = rba.xml.TargetGroup('flagella_activation')
    target = rba.xml.TargetReaction('Th')
    target.value = 'flagella_proton_flux'
    target_group.reaction_fluxes.append(target)
    return target_group


def flagella_activation_functions():
    return [
        rba.xml.Function('flagella_speed', 'constant', {'CONSTANT': 5.81}),
        rba.xml.Function('flagella_h_consumption', 'constant',
                         {'CONSTANT': 0.9415}),
        rba.xml.Function('number_flagella', 'linear',
                         {'LINEAR_COEF': 4.5197, 'LINEAR_CONSTANT': 3.7991,
                          'X_MIN': 0.25, 'X_MAX': 1.6,
                          'Y_MIN': float('-Inf'), 'Y_MAX': float('Inf')})
        ]


def flagella_activation_aggregate():
    aggregate = rba.xml.Aggregate('flagella_proton_flux', 'multiplication')
    aggregate.function_references.append(
        rba.xml.FunctionReference('number_flagella')
    )
    aggregate.function_references.append(
        rba.xml.FunctionReference('flagella_speed')
    )
    aggregate.function_references.append(
        rba.xml.FunctionReference('flagella_h_consumption')
    )
    return aggregate


if __name__ == "__main__":
    main()
