from gsmodutils.test.utils import ModelTestSelector
import json
import os
import numpy as np
from collections import defaultdict


@ModelTestSelector(models="*", designs="*", conditions="*")
def test_atpase_validity(model, project, log):
    """ 
    Test to check if the model is capable of ATP production without media.
    Should be ran on all designs and condtions.
    """
    for reaction in model.exchanges: #constrain all transporters to allow no uptake
        reaction.lower_bound = 0

    model.reactions.get_by_id('ATPASE-RXN').bounds = (1,1)
    solution = model.solver.optimize()
    log.assertion(solution == "infeasible", "Valid model conditions", "Model capable of ATP production from nothing")
    
    
def test_acetate(model, project, log):
    """ Test that acetate production hasn't changed significantly """
    s = model.optimize()
    condition = s.fluxes["EX_ACET"] < 2.0 or s.fluxes["EX_ACET"] > 2.8
    log.warning(condition, "Acetate production outside expected range, value - {}".format(s.fluxes["EX_ACET"]))
    
    
@ModelTestSelector(designs=["fluxm_flux_minimisation", "fluxm_flux_minimisation_h2_restricted"])
def test_minimisation_reactions(model, project, log):
    """ Checks scrumpy reactions are contained in the flux minimisation model"""
    #define list of enzyme-catalysed reactions based on PGDB from ScrumPy (i.e. exclude spontaneous and diffusion reactions)
    react_list = os.path.join(project.project_path, 'reacs_for_min.json')
    with open(react_list,'r') as f:
        scrumpy_min_reactions = json.load(f)
    #list of cobra model reaction ids
    cobra_reactions = [x.id for x in model.reactions]
    log.assertion(
        len(set(scrumpy_min_reactions).intersection(cobra_reactions)) == len(scrumpy_min_reactions),
        "ALL SCRUMPY ENZYME-CATALYSED REACTIONS ARE IN COBRA MODEL",
        "SOME SCRUMPY ENZYME-CATALYSED REACTIONS ARE NOT IN COBRA MODEL")


def co_uptake_spec(model, product_rxns, co_uptake_range=np.arange(15,50,1), carbon_source='EX_CARBON-MONOXIDE'):
    products = defaultdict(list)
    #fba scan
    csource = model.reactions.get_by_id(carbon_source)
    for co_uptake in co_uptake_range:
        csource.bounds = (-co_uptake,-co_uptake)
        sol = model.optimize()
        v = sol.fluxes
        for rxn in product_rxns:
            products[rxn] += v[rxn]
    return products
    

def strictly_increasing(l):
    return all(x<y for x, y in zip(l, l[1:]))

def strictly_decreasing(l):
    return all(x>y for x, y in zip(l, l[1:]))

def non_increasing(l):
    return all(x>=y for x, y in zip(l, l[1:]))

def non_decreasing(l):
    return all(x<=y for x, y in zip(l, l[1:]))


@ModelTestSelector(designs=["fluxm_flux_minimisation"])
def test_minimisation_fluxm_products(model, project, log):
    """
    This tests the product spectrum in hydrogen non-limited conditions
    """
    test_product_rxns = ["EX_ACET", 'EX_ETOH', 'EX_D-LACTATE', 'EX_BUTANEDIOL', 'EX_HYDROGEN-MOLECULE', 'EX_CARBON-DIOXIDE']
    products = co_uptake_spec(model, test_product_rxns)
    
    log.assertion(
        strictly_increasing(products['EX_CARBON-DIOXIDE']),
        "Carbon dioxide production is not monotonically increasing",
        "Carbon dioxide production is monotonically increases",
        desc="Carbon dioxide prodction should increase"
    )
    
    # NO lactate of butanediol should be produced
    log.assertion(
        len([x for x in products['EX_D-LACTATE'] if x != 0]) == 0,
        "No lactate is produced",
        "lactate is produced",
        desc="Test for lactate productions"
    )
    
    log.assertion(
        len([x for x in products['EX_BUTANEDIOL'] if x != 0]) == 0,
        "No butanediol is produced",
        "butanediol is produced",
        desc="Test for butanediol productions"
    )
    
    # Hydrogen should increase monitonically
    log.assertion(
        strictly_increasing(products['EX_HYDROGEN-MOLECULE']),
        "hydrogen is produced",
        "No hydrogen is produced",
        desc="Test for hydrogen production"
    )
    
    
@ModelTestSelector(designs=["fluxm_flux_minimisation_h2_restricted"])
def test_minimisation_fluxm_h2_restricted_products(model, project, log):
    """
    This tests the product spectrum in hydrogen conditions
    """
    test_product_rxns = ["EX_ACET", 'EX_ETOH', 'EX_D-LACTATE', 'EX_BUTANEDIOL', 'EX_HYDROGEN-MOLECULE', 'EX_CARBON-DIOXIDE']
    products = co_uptake_spec(model, test_product_rxns)
    
    log.assertion(
        len([x for x in products['EX_HYDROGEN-MOLECULE'] if x != 0]) == 0,
        "No hydrogen is produced",
        "hydrogen is produced",
        desc="Test for hydrogen production"
    )
    
    log.assertion(
        strictly_increasing(products['EX_CARBON-DIOXIDE']),
        "Carbon dioxide production is not monotonically increasing",
        "Carbon dioxide production is monotonically increases",
        desc="Carbon dioxide prodction should increase"
    )
    
    # Lactate production should increase
    # Check that the acetate and lactate production are not different in hydrogen limited and hydorgren producing setups
    no_h2_model = project.load_design("fluxm_flux_minimisation")
    no_h2_products = co_uptake_spec(no_h2_model, test_product_rxns)
    
    log.assertion(
        tuple(products["EX_ACET"]) == tuple(no_h2_products["EX_ACET"]),
        "Acetate production remains unchanged",
        "Acetate production changes",
        desc="check to see if there is a difference in acetate production in non-hydrogen limited conditions"
    )
    
    log.assertion(
        tuple(products["EX_ETOH"]) == tuple(no_h2_products["EX_ETOH"]),
        "Ethanol production remains unchanged",
        "Ethanol production changes",
        desc="check to see if there is a difference in acetate production in non-hydrogen limited conditions"
    )
    

@ModelTestSelector(models="*", designs="*", conditions="*")
def test_essential_reactions(model, project, log):
    
    singles_path = os.path.join(project.project_path, 'essential_singles.json')
    with open(singles_path) as sp:
        essential_reactions = json.load(sp)
        
    for rxn in essential_reactions:
        try:
            reaction = model.reactions.get_by_id(rxn)
        except KeyError:
            log.error("Reaction {} is listed as essential but not found in the model".format(rxn))
            continue
        
        bnds = reaction.bounds
        reaction.bounds = (0,0)
        
        reaction2 = None
        try:
            reaction2 = model.reactions.get_by_id(rxn +"_reverse")
            r2bnds = reaction2.bounds
            reaction2.bounds = (0,0)
        except KeyError:
            pass
        
        sol = model.optimize()
        reaction.bounds = bnds
        
        if reaction2 is not None:
            reaction2.bounds = r2bnds
        
        log.assertion(
            sol.status == "infeasible" or sol.objective_value < 1e-3,
            "Essential reaction {} remains essential".format(rxn),
            "Reaction lited as essential {} is not essential - solution = {}".format(rxn, sol.objective_value)
        )
        
        