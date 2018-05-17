from gsmodutils.test.utils import ModelTestSelector


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