from gsmodutils.test.utils import ModelTestSelector
import json
import os


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
