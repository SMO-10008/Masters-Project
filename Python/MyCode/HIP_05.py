# This script creates BOFs for new uptake rates, using a set of uptake rates and corresponding BOF coefficients
# based on the HIP method.
# Shantala Marie Ødegården, 2021.

# imports
from cobra.io import read_sbml_model, write_sbml_model
from cobra import Model, Reaction, Metabolite
import numpy as np
import math


# class for getting mw data
class MW:
    def __init__(self, met_id, mw, original_coef, scaled_coef):
        self.met_id = met_id
        self.mw = mw
        self.original_coef = original_coef
        self.scaled_coef = scaled_coef


# read model
# model = read_sbml_model('iML1515_withBOFs.xml')
model = read_sbml_model('rounded_ec_iML1515.xml')

# set solver
model.solver = 'gurobi'

# get biomass reactions
biomass_ul = model.reactions.get_by_id('BIOMASS_Ec_iML1515_UL_75p37M')
biomass_nl = model.reactions.get_by_id('BIOMASS_Ec_iML1515_NL_75p37M')
biomass_cl = model.reactions.get_by_id('BIOMASS_Ec_iML1515_CL_75p37M')

biomass_all = [biomass_ul, biomass_nl, biomass_cl]


# artificial in vivo uptake rates
c_uptake_ul = 18
c_uptake_nl = 13.5
c_uptake_cl = 1.5
n_uptake_ul = 8.5
n_uptake_nl = 1.5
n_uptake_cl = 0.68

c_uptakes = [c_uptake_ul, c_uptake_nl, c_uptake_cl]
n_uptakes = [n_uptake_ul, n_uptake_nl, n_uptake_cl]

# matrix of c and n uptake for UL, NL, CL
old_uptake_vals = np.array([c_uptakes, n_uptakes])

# make matrix with zeros, will consist of biomass coefficients for
# row 0: UL, row 1: NL, row 2: CL
# coefficients will all be positive here, doesn't matter if they are reactants or products
old_biomasses = np.zeros((len(biomass_ul.metabolites), len(biomass_all)))

n_old_vals = len(biomass_all)

# add biomass coefficients to matrix
i = 0
for met in biomass_ul.metabolites:
    for j in range(0, n_old_vals):
        biomass_rxn = biomass_all[j]
        old_biomasses[i, j] = abs(biomass_rxn.get_coefficient(met.id))
    i += 1

# UL = A, NL = B, CL = C. new = D Find relative vectors: AB and AC
# these vectors are used later to find a linear combination of AB and AC from A to the new environment (C,N)
# resulting in vector weightings s, t such that AD = s*AB + t*AC
new_origin = old_uptake_vals[:, 0]
vectors = old_uptake_vals[:, 1:]
relative_vectors = np.zeros((len(vectors), len(vectors[0])))
# len(vectors[1]) should give number of "columns"
for i in range(0, len(vectors[0])):
    relative_vectors[:, i] = vectors[:, i] - new_origin


for i in range(0, 200, 5):
    i_adj = 0.1 * i
    for j in range(0, 200, 5):
        j_adj = 0.1 * j
        print(i_adj, j_adj)
        new_uptake_vals = np.array((i_adj, j_adj))
        # distances from each point (UL c-uptake, UL n-uptake) (NL c, NL n) (CL c, CL n) to the new point (new c, new n
        distances = np.zeros(n_old_vals)
        # k = 0: UL, k = 1: NL, k = 2: CL
        for k in range(0, n_old_vals):
            # eucledian distance: sqrt((x1-x2)^2 + (y1-y2)^2)
            distances[k] = math.sqrt((old_uptake_vals[0, k] - new_uptake_vals[0])**2
                                     + (old_uptake_vals[1, k] - new_uptake_vals[1])**2)

        # vector from new origin (UL) point in C,N uptake plane to the new env point (C,N)
        relative_vector_new = new_uptake_vals - new_origin

        # vector weights: s and t in AD = s*AB + t*AC
        vector_weights = np.linalg.solve(relative_vectors, relative_vector_new)

        new_origin_biomass = old_biomasses[:, 0]
        relative_biomasses = np.zeros((len(biomass_ul.metabolites), len(biomass_all)))
        for k in range(0, len(relative_biomasses[0])):
            relative_biomasses[:, k] = old_biomasses[:, k] - new_origin_biomass

        new_biomass = new_origin_biomass + np.matmul(relative_biomasses[:, 1:], vector_weights)

        # make new biomass reaction
        biomass_rxn_new = Reaction('BIOMASS_Ec_iML1515_HIP_C' + str(i_adj).replace('.', '-') +
                                   '_N' + str(j_adj).replace('.', '-'))
        biomass_rxn_new.name = 'E. coli biomass objective function (iML1515) - HIP: (C,N) = ' \
                               '(' + str(i_adj).replace('.', '-') + ', ' \
                               + str(j_adj).replace('.', '-') + ')'
        biomass_rxn_new.lower_bound = 0
        biomass_rxn_new.upper_bound = 1000

        # if coefficient negative: use 1% of UL
        k = 0
        for met in biomass_ul.metabolites:
            # get biomass coefficient in UL
            ul_coef = biomass_ul.get_coefficient(met.id)
            sign = 0
            # check if UL coefficient is pos or neg, indicating if metabolite is product or reactant
            if ul_coef > 0:
                sign = 1
            elif ul_coef < 0:
                sign = -1
            # get coefficient from new biomass
            new_coef = new_biomass[k]
            if new_coef >= 0:
                # new coefficient is positive and will be used
                biomass_rxn_new.add_metabolites({met: new_biomass[k]*sign})
            else:
                # new coefficient is negative, so use 1% of UL coefficient instead
                biomass_rxn_new.add_metabolites({met: ul_coef*0.01})
            k += 1
        # add metabolites with new biomass to

        model.add_reactions([biomass_rxn_new])


# read file with biomass functions, add lines to array
f = open('Biomasses_edited.txt')
all_lines = f.readlines()
f.close()

all_lines_array = np.array(all_lines)

# lists for storing the data of each line
metabolite_list = []
class_list = []
mw_list = []
wt_list = []
core_list = []
mean_list =[]
carbstarv_list = []
nitstarv_list = []

# read line by line
for i in range(len(all_lines_array)):
    current_line = all_lines_array[i]
    current_line_list = current_line.split('\t')
    if i == 0:
        metabolite_list = current_line_list[1:-2]
    elif i == 1:
        class_list = current_line_list[1:]
    elif i == 2:
        for j in range(1, len(current_line_list)-2):
            if current_line_list[j] == '':
                mw_list.append(current_line_list[j])
            else:
                mw_list.append(float(current_line_list[j]))
    elif i == 3:
        for j in range(1, len(current_line_list)-2):
            wt_list.append(current_line_list[j])
    elif i == 4:
        for j in range(1, len(current_line_list)-2):
            core_list.append(current_line_list[j])
    elif i == 5:
        for j in range(1, len(current_line_list)-2):
            mean_list.append(current_line_list[j])
    elif i == 6:
        for j in range(1, len(current_line_list)-2):
            carbstarv_list.append(current_line_list[j])
    elif i == 7:
        for j in range(1, len(current_line_list)-2):
            nitstarv_list.append(current_line_list[j])

# make list of objects containing metabolite id and molecular weight
mw_obj_list = []

# make objects, set original and scaled coef to 0 for now
for i in range(len(metabolite_list)):
    if not metabolite_list[i] == '->':
        mw_obj_list.append(MW(metabolite_list[i], mw_list[i], 0, 0))


# normalize new biomass reactions
for biomass_hip in model.reactions:
    if biomass_hip.id.startswith('BIOMASS_Ec_iML1515_HIP_C'):
        print('ID: ', biomass_hip.id)
        # SCALE TO 1 g/(h * CDW)
        mass_before = 0

        # find original mass and set original coefficient in mw object
        for met in biomass_hip.metabolites:
            for mw_obj in mw_obj_list:
                if met.id == mw_obj.met_id:
                    met_coef = biomass_hip.get_coefficient(met.id)
                    mass_before += met_coef * mw_obj.mw
                    mw_obj.original_coef = met_coef

        print('mass before: ', mass_before)

        # get factor to scale each coefficient
        normalization_factor = -1000 / mass_before

        # set scaled coef for each object
        for mw_obj in mw_obj_list:
            mw_obj.scaled_coef = mw_obj.original_coef * normalization_factor

        # check that the new mass is -1000 because 1000 mg should "disappear" and turn into biomass
        sum_mass_scaled = 0
        for mw_obj in mw_obj_list:
            sum_mass_scaled += mw_obj.mw * mw_obj.scaled_coef

        mass_after = round(sum_mass_scaled * 100) / 100
        print('mass after: ', mass_after)

        # print of old reaction
        # for met in biomass_hip.metabolites:
        #     print(met.id, biomass_hip.get_coefficient(met.id))

        # set new coefficients in biomass reaction
        for met in biomass_hip.metabolites:
            for mw_obj in mw_obj_list:
                if met.id == mw_obj.met_id:
                    # remove old amount
                    biomass_hip.add_metabolites({met: -mw_obj.original_coef})
                    # add new amount
                    biomass_hip.add_metabolites({met: mw_obj.scaled_coef})

        # print new reaction
        # for met in biomass_hip.metabolites:
        #     print(met.id, met.formula, biomass_hip.get_coefficient(met.id))

write_sbml_model(model, "rounded_ec_iML1515_HIP_05.xml")
