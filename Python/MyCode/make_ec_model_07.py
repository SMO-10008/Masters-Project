# This script creates an enzymatically constrained version of E. coli genome-scale metabolic model iML1515
# based on the GECKO method.
# It does not take any inputs, if one wants to change the upper bound for enzyme amount, it is marked with TODO
# and can be changed directly in the code.
# Shantala Marie Ødegården, 2021.

# import
from cobra.io import read_sbml_model, write_sbml_model


# read files with enzymatic data, store the data as objects
def read_data():
    import numpy as np

    # classes to save the information from data from files:

    # Kcat class
    class Kcat:
        def __init__(self, react_id, kcat_value):
            self.react_id = react_id
            self.kcat_value = kcat_value

    # Enzyme class
    class Enzyme:
        def __init__(self, name, mw, bigg_genes, alt_names):
            self.name = name
            self.mw = mw
            self.bigg_genes = bigg_genes
            self.alt_names = alt_names

    # read kcat file
    f = open('kcat.txt')
    all_lines_kcat = f.readlines()
    f.close()
    all_lines_array_kcat = np.array(all_lines_kcat)

    # list of Kcat objects
    all_kcats = []

    # read kcat file, line by line
    for i in range(1, len(all_lines_array_kcat)):
        current_line = all_lines_array_kcat[i]
        current_line_list = current_line.split("\t")
        current_react_id = current_line_list[0]
        current_kcat_value = float(current_line_list[7])
        # add to list of Kcat objects
        all_kcats.append(Kcat(current_react_id, current_kcat_value))

    # read enzyme file
    f = open('enzymes_01.txt', 'r')
    all_lines_enzyme = f.readlines()
    f.close()
    all_lines_array_enzyme = np.array(all_lines_enzyme)

    # list of Enzyme objects
    all_enzymes = []

    # read enzyme file, line by line
    for i in range(1, len(all_lines_array_enzyme)):
        current_line = all_lines_array_enzyme[i]
        current_line_list = current_line.split("\t")

        # split line into their cell values - enzyme name, molecular weight, and genes
        current_enzyme = current_line_list[0]
        current_mw = current_line_list[1]
        current_genes = current_line_list[2]

        # get molecular weight:
        # handle multiple MW. will choose the highest value
        current_mw_list_string = current_mw.split('//')

        # if there is no mw, set to nan - will be set to the average of all mws later
        for mw in current_mw_list_string:
            if len(mw) == 0:
                current_mw_list_string = ['nan']

        # list of mws of an enzyme:
        current_mw_list_float = []

        # convert from string to float
        for mw in current_mw_list_string:
            current_mw_list_float.append(float(mw))

        # if multiple mw, chose biggest mw
        if len(current_mw_list_float) > 1:
            current_mw_float = max(current_mw_list_float)
        else:
            current_mw_float = current_mw_list_float[0]

        # Da = g/mol, so convert kDa to Da (and g/mol) :
        current_mw_float = current_mw_float / 1000

        # get genes:
        current_genes = current_genes[:-1]
        current_genes_list = current_genes.split('//')
        processed_genes_list = []

        # set gene to '---' if no genes are listed for the enzyme
        for gene in current_genes_list:
            if len(gene) == 0:
                processed_genes_list.append('---')
            else:
                # remove leading and trailing spaces, and add to processed list
                processed_genes_list.append(gene.strip())

        # remove gene names that are not BiGG identifiers, as these will be used to match enzymes to reactions
        # all BiGG identifiers in the model starts with 'b', except s0001, but this one is not among the enzyme genes
        # the 'b' is also always followed by a number
        final_genes_list = []
        # other gene names go in separate list
        current_alt_names = []
        for gene in processed_genes_list:
            # print(gene, gene.startswith('b'))
            gene_sub = gene[1:]
            if gene.startswith('b') and gene_sub[0].isdigit():
                final_genes_list.append(gene)
            else:
                current_alt_names.append(gene)

        # add to list of Enzyme objects
        all_enzymes.append(Enzyme(current_enzyme, current_mw_float, final_genes_list, current_alt_names))

    # find average, minimum and maximum mw
    # note: this is for all enzymes from which data was collected
    # which may be more enzymes than just those actually included in the model
    counter = 0
    mw_sum = 0
    mw_min = 1000
    mw_max = 0
    for e in all_enzymes:
        if not np.isnan(e.mw):
            mw_sum += e.mw
            counter += 1
            # check if max or min mw
            if e.mw > mw_max:
                mw_max = e.mw
            if e.mw < mw_min:
                mw_min = e.mw
    mw_avg = mw_sum / counter
    # print(mw_avg, mw_min, mw_max)

    # set all enzymes without mw to the average of about 1.2e-18 grams (78 kDa)
    # not ideal, as mws range from about 1.1e-19 grams (7kDa) to 1.6e-17 grams (977 kDa)
    # but they need to have some mw
    for e in all_enzymes:
        if np.isnan(e.mw):
            e.mw = mw_avg
    # print('number of enzyme objects: ', len(all_enzymes))
    # print('number of kcat objects: ', len(all_kcats))
    # print('read data complete')
    return all_enzymes, all_kcats


# split reversible reactions into two separate reactions
def split_reversible_reactions(model):
    from cobra import Reaction

    # split reversible reactions into a couple going in each direction
    for rxn in model.reactions:
        # if reaction is reversible:
        # this should exclude export reactions
        if rxn.lower_bound < 0 and not rxn.id.startswith('EX'):
            # create an opposite reaction, id: [reaction id]_inv, name: [reaction name] inverse
            rxn_inv = Reaction(rxn.id + '_inv')
            rxn_inv.name = rxn.name + ' inverse'
            rxn_inv.lower_bound = 0
            rxn_inv.upper_bound = -rxn.lower_bound
            for m in rxn.metabolites:
                rxn_inv.add_metabolites({m: -rxn.get_coefficient(m)})
            model.add_reactions([rxn_inv])
            rxn.lower_bound = 0
            rxn_inv.gene_reaction_rule = rxn.gene_reaction_rule

    print('split reversible reactions complete')


# incorporate enzymes into the model
def add_enzymes(model, all_enzymes, all_kcats):
    from cobra import Reaction, Metabolite

    # a pseudo-metabolite that represents the enzyme pool
    # all the enzymes draw mass from this pool
    # set compartment to e and use import reaction - this seems to work
    # as opposed to setting compartment to c and have no import - I guess the metabolite must come from somewhere
    # may be a better way to do this - look into later
    enzyme_pool_met = Metabolite('enzyme_pool', name='Enzyme pool pseudo-metabolite', compartment='e',
                                 formula=None, charge=None)

    # enzyme pool import reaction
    # seems to fix infeasible solution when enzyme metabolite compartments are set to 'c'
    # can be used to restrict the total enzyme flux
    # importing enzymes does not represent what really happens in the cell - look for different ways later
    EX_pool_enzyme = Reaction('EX_pool_enzyme')
    EX_pool_enzyme.lower_bound = 0

    # set upper bound to constrain model with enzyme amount
    # TODO : set ub
    # alternative 1: use p*f*s
    # source: 70.1 % of dry weight is protein
    p = 0.701
    f = 0.4461
    s = 0.51
    ub = round(p * f * s, 2)
    # alternative 2: set ub directly
    # ub = 0.7

    EX_pool_enzyme.upper_bound = ub
    EX_pool_enzyme.add_metabolites({enzyme_pool_met: 1})
    model.add_reactions([EX_pool_enzyme])

    # counters
    num_enzymes_added = 0
    enzymes_in_model = []
    enzymatic_reactions = []

    # add enzymes to reactions
    for rxn in model.reactions:
        for g in rxn.genes:
            for e in all_enzymes:
                for g_enz in e.bigg_genes:
                    if g.id == g_enz:
                        # find reaction kcat
                        for k in all_kcats:
                            if k.react_id == rxn.id:
                                if rxn.id not in enzymatic_reactions:
                                    enzymatic_reactions.append(rxn.id)
                                # add enzymes as metabolites, in both directions if reaction is reversible
                                # 5 decimals on stoichiometric coefficient
                                # put formula and charge to None - would mess up mass balance if added
                                # had to specify compartment to be able to save model, just set it to 'c' for all
                                enzyme_id = 'enzyme_' + g_enz

                                # check if metabolite is already added to the model
                                enzyme_test_bool = 0
                                for met in model.metabolites:
                                    if met.id == enzyme_id:
                                        # if this happens, enzyme metabolite is already added
                                        enzyme_test_bool = 1

                                # if already added, get existing metabolite
                                if enzyme_test_bool == 1:
                                    enzyme_metabolite = model.metabolites.get_by_id(enzyme_id)

                                else:
                                    # if not already added, make new metabolite
                                    enzyme_metabolite = Metabolite(enzyme_id, name=e.name, compartment='c',
                                                                   formula=None, charge=None)
                                    num_enzymes_added += 1
                                    enzymes_in_model.append(e)

                                # assuming that kcat is the same for both directions
                                # kcat at this point is applied by reaction
                                # if different enzymes catalyze same reaction
                                # same kcat is applied, but is likely different in real life

                                coef = -round(1 / k.kcat_value, 5)
                                rxn.add_metabolites({enzyme_metabolite: coef})

                                # add a reaction where enzyme draws mass from enzyme pool
                                # add new pool reaction if it doesn't already exist
                                enz_rxn_id = 'POOL_' + enzyme_id
                                rxn_test_bool = 0
                                for rxn_test in model.reactions:
                                    if rxn_test.id == enz_rxn_id:
                                        rxn_test_bool = 1
                                if rxn_test_bool == 0:
                                    enzyme_reaction = Reaction(enz_rxn_id)
                                    enzyme_reaction.name = e.name + ' from enzyme pool'
                                    enzyme_reaction.lower_bound = 0
                                    enzyme_reaction.upper_bound = 1000
                                    enzyme_reaction.add_metabolites(
                                        {enzyme_pool_met: -e.mw, enzyme_metabolite: 1})
                                    model.add_reactions([enzyme_reaction])

    # print('number of enzymes added: ', num_enzymes_added)
    # print(len(enzymes_in_model))
    all_enzymes_mass = 0
    model_enzymes_mass = 0
    for e in all_enzymes:
        all_enzymes_mass += e.mw
    for e in enzymes_in_model:
        model_enzymes_mass += e.mw
    # print(all_enzymes_mass, model_enzymes_mass)

    print('add enzymes complete')


# make arm reactions so that isozyme reactions may happen separately, but still keep original upper bound
def make_arm_reactions(model):
    from typing import List, Any
    from cobra import Reaction, Metabolite
    import re
    # check reactions with multiple added enzymes
    # make "arm reactions" for isozymes

    # class for keeping track of partial rules of reactions - used for creating new reactions for isozymes
    # this ensures that the part of reaction.gene_reaction_rule is preserved, prevents data loss from the processing
    class SubRule:
        def __init__(self, original_rule, processed_rule):
            self.original_rule = original_rule
            self.processed_rule = processed_rule

    # for isozyme sub-reactions created, the original reactions must be removed, they will be added to this list
    reactions_to_remove: List[Any] = []

    for rxn in model.reactions:
        # list of SubRule objects
        sub_rule_list = []

        # this list contains the bigg names of the genes of the added enzymes
        enzyme_gene_list = []
        for met in rxn.metabolites:
            if met.id.startswith('enzyme') and not met.id.endswith('pool'):
                enzyme_id = met.id
                enzyme_gene = enzyme_id.replace('enzyme_', '')
                enzyme_gene_list.append(enzyme_gene)

        # now to find out if any of the enzyme genes are isozymes
        # only look into reactions with more than one enzyme in it
        if len(enzyme_gene_list) > 1:
            # print('reaction with more than one enzyme: ', rxn)
            rule_original = rxn.gene_reaction_rule
            rule_list_original = rule_original.split(' or ')
            rule_list = rule_original.split(' or ')
            # print('BEFORE: ', rule_list)

            # make the rules into a list of lists, containing genes as:
            # [solo_enzyme1, [complex_enzymeA, complex_enzymeB]]
            for i in range(len(rule_list)):
                sub_rule = rule_list[i]
                complex_list = sub_rule.split(' and ')
                # if there are any 'and's
                # if len(complex_list) > 1:
                complex_list_processed = []
                for j in range(len(complex_list)):
                    complex_unprocessed = complex_list[j]
                    complex_processed = re.sub('[^A-Za-z0-9]+', '', complex_unprocessed)
                    complex_list_processed.append(complex_processed)
                rule_list[i] = complex_list_processed

            # print('AFTER: ', rule_list)

            # only proceed if there is one or more 'or's in the rule
            if len(rule_list) > 1:

                # make a list of lists based on the rule list, but will only contain the added enzymes
                # will not add genes that are not added enzymes
                enzyme_rule_list = []

                # print('original rule ', rule_original)
                # print('rule list', rule_list)
                # print('enzyme gene list', enzyme_gene_list)

                for i in range(len(rule_list)):
                    sub_rule = rule_list[i]
                    sub_rule_enzymes = []
                    for j in range(len(sub_rule)):
                        # print('sub rule gene ', sub_rule[j])
                        sub_rule_gene = sub_rule[j]
                        if sub_rule_gene in enzyme_gene_list:
                            # print('is in enzyme list: ', sub_rule_gene)
                            gene_to_add = sub_rule_gene
                            sub_rule_enzymes.append(gene_to_add)
                    enzyme_rule_list.append(sub_rule_enzymes)
                    # print('enzyme rule list: ', enzyme_rule_list)

                # print('o: ', rule_list_original, len(rule_list_original))
                # print('e: ', enzyme_rule_list, len(enzyme_rule_list))

                # make subrule objects and add to list
                for i in range(len(rule_list_original)):
                    sub_rule_obj = SubRule(rule_list_original[i], enzyme_rule_list[i])
                    sub_rule_list.append(sub_rule_obj)

                # sort processed subrules and delete those that turned out that have no enzymes added from the rule
                # as these won't need their own isozyme reaction
                for sub_rule in sub_rule_list:
                    p = sub_rule.processed_rule
                    p.sort()
                    if len(p) == 0:
                        sub_rule_list.remove(sub_rule)

                # remove duplicate subrules that may have been formed from sets of enzymes that are common to complexes
                # but that lack additional enzymes that have not been added but are part of the model

                # identify duplicates to remove
                unique_rules = []
                rules_to_remove = []
                for i in range(len(sub_rule_list)):
                    sub_rule = sub_rule_list[i]
                    p = sub_rule.processed_rule
                    if p not in unique_rules:
                        unique_rules.append(p)
                    else:
                        # print('duplicate found')
                        rules_to_remove.append(sub_rule)

                # remove the identified duplicate rules
                for i in range(len(rules_to_remove)):
                    sub_rule_list.remove(rules_to_remove[i])

                # sub_rule_list_processed should now indicate which enzymes go together in which sub-reactions
                # some might not need their own reaction still - possible that only one subunit that is part of
                # multiple complexes was added as enzyme. so will only proceed if more than one sub-rule remains
                # this doesn't seem to be the case, but it is something that could have happened with different data

                if len(sub_rule_list) > 1:
                    # now it's time to create and add the sub-reactions for the isozymes

                    # create an intermediate pseudo-metabolite - should be one for each reaction
                    ipm_id = rxn.id + '_INTERMEDIATE'
                    ipm_name = rxn.id + 'Intermediate pseudo-metabolite '
                    ipm_met = Metabolite(ipm_id, name=ipm_name, compartment='c', formula=None, charge=None)

                    # make reactions:  intermediate pseudomet + enzyme(s) --> prod
                    for i in range(len(sub_rule_list)):
                        sub_rule = sub_rule_list[i]
                        sub_reaction_id = rxn.id + '_' + 'ARM_2_' + str(i)
                        sub_reaction = Reaction(sub_reaction_id)
                        sub_reaction.name = rxn.id + ' ' + 'Arm reaction 2_' + str(i)
                        sub_reaction.lower_bound = rxn.lower_bound
                        sub_reaction.upper_bound = rxn.upper_bound
                        sub_reaction.gene_reaction_rule = sub_rule.original_rule

                        # add original products:
                        for product in rxn.products:
                            sub_reaction.add_metabolites({product: rxn.get_coefficient(product)})

                        # add intermediate as reactant
                        sub_reaction.add_metabolites({ipm_met: -1})

                        # add enzyme metabolites
                        for g in sub_rule.processed_rule:
                            for e in all_enzymes:
                                for g_enz in e.bigg_genes:
                                    if g == g_enz:
                                        # find reaction kcat
                                        for k in all_kcats:
                                            if k.react_id == rxn.id:
                                                enzyme_id = 'enzyme_' + g_enz

                                                # check if metabolite is already added to the model
                                                enzyme_test_bool: int = 0
                                                for met in model.metabolites:
                                                    if met.id == enzyme_id:
                                                        # if this happens, enzyme metabolite is already added
                                                        enzyme_test_bool = 1

                                                # if already added, get existing metabolite
                                                if enzyme_test_bool == 1:
                                                    enzyme_metabolite = model.metabolites.get_by_id(enzyme_id)

                                                else:
                                                    # if not already added, make new metabolite
                                                    enzyme_metabolite = Metabolite(enzyme_id, name=e.name,
                                                                                   compartment='c',
                                                                                   formula=None, charge=None)

                                                # assuming that kcat is the same for both directions
                                                # kcat at this point is applied by reaction
                                                # if different enzymes catalyze same reaction
                                                # same kcat is applied, but is likely different in real life

                                                coef = -round(1 / k.kcat_value, 5)
                                                sub_reaction.add_metabolites({enzyme_metabolite: coef})

                        # print('sub 2: ', sub_reaction)
                        # add sub-reaction to model
                        model.add_reactions([sub_reaction])

                    # make part 1 of reaction: intermediate pseudo-met + enz --> products
                    sub_reaction_id = rxn.id + '_ARM_1'
                    sub_reaction = Reaction(sub_reaction_id)
                    sub_reaction.name = rxn.id + ' Arm reaction 1'
                    sub_reaction.lower_bound = rxn.lower_bound
                    sub_reaction.upper_bound = rxn.upper_bound
                    sub_reaction.gene_reaction_rule = rxn.gene_reaction_rule

                    # add original reactants to reaction
                    for reactant in rxn.reactants:
                        if not reactant.id.startswith('enzyme'):
                            sub_reaction.add_metabolites({reactant: rxn.get_coefficient(reactant)})
                    # add intermediate pseudo-metabolite to reaction
                    sub_reaction.add_metabolites({ipm_met: 1})

                    # print('sub 1: ', sub_reaction)
                    # add sub reaction 1 to model
                    model.add_reactions([sub_reaction])

                    # add to reactions to remove
                    reactions_to_remove.append(rxn)

    # print('number of reactions before removing: ', len(model.reactions))
    # print('number of reactions to remove:', len(reactions_to_remove))
    # delete original reaction
    model.remove_reactions(reactions_to_remove)
    # print('number of reactions after removing: ', len(model.reactions))
    print('make arm reactions complete')


def add_bofs(model):
    from cobra import Model, Reaction, Metabolite
    import numpy as np
    # get biomass functions

    # read file with biomass functions, add lines to array
    f = open('Biomasses_edited.txt')
    all_lines = f.readlines()
    f.close()

    all_lines_array = np.array(all_lines)

    # lists for storing the data of each line
    metabolite_list = []
    # class_list = []
    # mw_list = []
    wt_list = []
    core_list = []
    mean_list = []
    carbstarv_list = []
    nitstarv_list = []

    # read line by line
    for i in range(len(all_lines_array)):
        current_line = all_lines_array[i]
        current_line_list = current_line.split('\t')
        # print(len(current_line_list))
        # print(current_line_list)
        if i == 0:
            metabolite_list = current_line_list[1:-2]
        # elif i == 1:
        #     class_list = current_line_list[1:]
        # elif i == 2:
        #     for j in range(1, len(current_line_list)):
        #         mw_list.append(float(current_line_list[j]))
        elif i == 3:
            for j in range(1, len(current_line_list) - 2):
                wt_list.append(current_line_list[j])
        elif i == 4:
            for j in range(1, len(current_line_list) - 2):
                core_list.append(current_line_list[j])
        elif i == 5:
            for j in range(1, len(current_line_list) - 2):
                mean_list.append(current_line_list[j])
        elif i == 6:
            for j in range(1, len(current_line_list) - 2):
                carbstarv_list.append(current_line_list[j])
        elif i == 7:
            for j in range(1, len(current_line_list) - 2):
                nitstarv_list.append(current_line_list[j])

    # make new biomass reactions UL (unlimited), CL (carbon limited), NL (nitrogen limited)
    biomass_ul = Reaction('BIOMASS_Ec_iML1515_UL_75p37M')
    biomass_ul.name = 'E. coli biomass objective function (iML1515) - Unlimited - with 75.37 GAM estimate'
    biomass_ul.lower_bound = 0
    biomass_ul.upper_bound = 1000

    biomass_cl = Reaction('BIOMASS_Ec_iML1515_CL_75p37M')
    biomass_cl.name = 'E. coli biomass objective function (iML1515) - Carbon Limited - with 75.37 GAM estimate'
    biomass_cl.lower_bound = 0
    biomass_cl.upper_bound = 1000

    biomass_nl = Reaction('BIOMASS_Ec_iML1515_NL_75p37M')
    biomass_nl.name = 'E. coli biomass objective function (iML1515) - Nitrogen Limited - with 75.37 GAM estimate'
    biomass_nl.lower_bound = 0
    biomass_nl.upper_bound = 1000

    # negative for reactants, positive for products
    coef_sign = -1

    # add metabolites to reactions
    for i in range(len(metabolite_list)):
        # if arrow, skip adding metabolites and change coefficient sign, as metabolites added after will be products
        if metabolite_list[i] == '->':
            coef_sign = 1
        else:
            ul_coef = float(mean_list[i])
            cl_coef = float(carbstarv_list[i])
            nl_coef = float(nitstarv_list[i])
            met = model.metabolites.get_by_id(metabolite_list[i])
            # print(metabolite_list[i], met)
            biomass_ul.add_metabolites({met: ul_coef * coef_sign})
            biomass_cl.add_metabolites({met: cl_coef * coef_sign})
            biomass_nl.add_metabolites({met: nl_coef * coef_sign})

    biomass_reactions = [biomass_ul, biomass_cl, biomass_nl]

    # add reactions to model
    model.add_reactions(biomass_reactions)
    print('add bofs complete')


def tune_kcats(model):
    from cobra import Reaction, Metabolite
    ub = model.reactions.get_by_id('EX_pool_enzyme').upper_bound
    # tune kcats for enzymes that take up a high percentage of enzyme mass
    # want to get growth to be about 0.6 at glucose uptake 8, o2 uptake 12
    target_growth = 0.6

    # save orginal bounds
    original_glc_ub = model.reactions.get_by_id('EX_glc__D_e').upper_bound
    original_glc_lb = model.reactions.get_by_id('EX_glc__D_e').lower_bound
    original_o2_ub = model.reactions.get_by_id('EX_o2_e').upper_bound
    original_o2_lb = model.reactions.get_by_id('EX_o2_e').lower_bound

    model.objective = model.reactions.get_by_id('BIOMASS_Ec_iML1515_UL_75p37M')
    model.reactions.get_by_id('EX_glc__D_e').upper_bound = -8
    model.reactions.get_by_id('EX_glc__D_e').lower_bound = -8
    model.reactions.get_by_id('EX_o2_e').upper_bound = -12
    model.reactions.get_by_id('EX_o2_e').lower_bound = -12
    sol = model.optimize()
    # print('original solution at glucose upt 8 and O2 upt 12: ', sol)

    while sol.objective_value < target_growth:
        # find enzyme that takes up the most of the enzyme mass
        biggest_mass_contribution = 0
        bmc_enzyme: Metabolite
        biggest_p = 0
        p_sum = 0
        for rxn in model.reactions:
            if rxn.id.startswith('POOL'):
                if rxn.flux > 0:
                    coef = -rxn.get_coefficient('enzyme_pool')
                    # p is mass contribution percentage
                    p = rxn.flux * coef / ub
                    p_sum += p
                    if p > biggest_mass_contribution:
                        biggest_mass_contribution = p
                        bmc_enzyme = model.metabolites.get_by_id(rxn.id.lstrip('POOL_'))
                        biggest_p = p

        # print('p_sum: %s, enzyme: %s, percentage: %s' % (p_sum, bmc_enzyme.id, biggest_p))
        # print('bmc enz id:', bmc_enzyme.id)
        # tune kcat (effectively coefficient) of that enzyme in that reaction
        # find reaction
        bmc_reaction: Reaction
        bmc = 0
        for rxn in bmc_enzyme.reactions:
            if not rxn.id.startswith('POOL'):
                coef = rxn.get_coefficient(bmc_enzyme.id)
                mc = rxn.flux * (-coef)
                if mc > bmc:
                    bmc = mc
                    bmc_reaction = rxn
        # print(bmc_reaction.id)
        # print(bmc_enzyme.id)
        # print('old rxn: ', bmc_reaction.reaction)
        # remove existing amount and add new amount
        old_coef = bmc_reaction.get_coefficient(bmc_enzyme.id)
        # "new_coef = old_coef * 0.1" resulted in some strange rounding errors, but the two next lines fixes the issue
        num_decimals = str(old_coef)[::-1].find('.')
        new_coef = round(old_coef * 0.1, num_decimals + 1)
        bmc_reaction.add_metabolites({bmc_enzyme: -old_coef})
        bmc_reaction.add_metabolites({bmc_enzyme: new_coef})
        # print(old_coef)
        # print(new_coef)
        # print('new rxn: ', bmc_reaction.reaction)
        # print("{0},{1},{2},{3}".format(bmc_reaction.id, bmc_enzyme.id, old_coef, new_coef))

        sol = model.optimize()
        # print('new solution: ', sol.objective_value)

    # reset glucose and o2 bounds
    model.reactions.get_by_id('EX_glc__D_e').upper_bound = original_glc_ub
    model.reactions.get_by_id('EX_glc__D_e').lower_bound = original_glc_lb
    model.reactions.get_by_id('EX_o2_e').upper_bound = original_o2_ub
    model.reactions.get_by_id('EX_o2_e').lower_bound = original_o2_lb

    print('tune kcats complete')


if __name__ == "__main__":
    all_enzymes, all_kcats = read_data()
    model = read_sbml_model('iML1515.xml')
    split_reversible_reactions(model)
    # print('#rxns and #mets after rev split: ', len(model.reactions), len(model.metabolites))
    add_enzymes(model, all_enzymes, all_kcats)
    # print('#rxns and #mets after adding enzymes: ', len(model.reactions), len(model.metabolites))
    make_arm_reactions(model)
    # print('#rxns and #mets after arm reactions: ', len(model.reactions), len(model.metabolites))
    add_bofs(model)
    # print('#rxns and #mets after adding BOFs: ', len(model.reactions), len(model.metabolites))
    tune_kcats(model)
    write_sbml_model(model, 'rounded_ec_iML1515.xml')
