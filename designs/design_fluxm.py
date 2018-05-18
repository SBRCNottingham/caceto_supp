import cobra
import os
import json


def gsmdesign_flux_minimisation(m, project):
    """ 
    Minimisation of flux accross reactions as objective function
    minimize absolute sum of enzyme-catalysed network flux (Row 1 of Table 3 in manuscript)
    """
    biomass = 0.028
    ngam = 2.154
    
    # biomass and atp maintenance constraints
    m.reactions.get_by_id('EX_BIOMASS').bounds = (biomass,biomass)
    m.reactions.get_by_id('ATPASE-RXN').bounds = (ngam,ngam)
    
    cobra.manipulation.convert_to_irreversible(m) # model to irreversible version for minimisation objective

    #define list of enzyme-catalysed reactions based on PGDB from ScrumPy (i.e. exclude spontaneous and diffusion reactions)
    react_list = os.path.join(project.project_path, 'reacs_for_min.json')
    
    with open(react_list,'r') as f:
        scrumpy_min_reactions = json.load(f)

    #list of cobra model reaction ids
    cobra_reactions = [x.id for x in m.reactions]

    #list of reverse enzyme-catalysed reactions for irreversible model objective function
    reverse_scrumpy_min_reactions = [x+'_reverse' for x in scrumpy_min_reactions if x+'_reverse' in cobra_reactions]

    #add reversed reaction ids to objective function list and define list of coresponding reactions objects
    min_reaction_ids = scrumpy_min_reactions+reverse_scrumpy_min_reactions
    min_reactions = [m.reactions.get_by_id(x) for x in min_reaction_ids]

    #define objective function
    m.objective.direction = 'min'
    m.objective = min_reactions

    #constrain biomass transporters (feature of ScrumPy models) to 0
    for reaction in m.reactions: 
        if '_bm_tx' in reaction.id:
            reaction.upper_bound = 0
            reaction.lower_bound = 0

    #enable production of biomass byproduct
    m.reactions.adenosyl_homocysteine_bm_tx.bounds = (-1000,0) 
    return m

gsmdesign_flux_minimisation.name = "flux minimisation"


def gsmdesign_flux_minimisation_h2_restricted(model, project):
    """ flux minimisation with h2 restriction """
    model.reactions.get_by_id('EX_HYDROGEN-MOLECULE').bounds = (0,0)
    return model

gsmdesign_flux_minimisation_h2_restricted.name = "flux minimisation h2 rest"
gsmdesign_flux_minimisation_h2_restricted.parent = "fluxm_flux_minimisation"