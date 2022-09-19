import pandas as pd
import reactome2py.content as rc

'''
This module gets the current Reactome release pathway hierarchy for use in the pathway scatter plot.

'''

def _traverse_pathway(children: list,  top_pathway: str, paths_found: set, pathways_df: pd.DataFrame) -> pd.DataFrame:
    for child in children:
        # skip non-pathways
        if child['type'] != 'Pathway':
            continue
        # keep track of newly discovered pathways
        if child['stId'] in paths_found:
            continue
        else:
            paths_found.add(child['stId'])
        # append newly discovered pathway to hierarchy dataframe
        pathways_df.loc[-1] = [child['name'], child['stId'], top_pathway]
        pathways_df.index = pathways_df.index + 1
        if 'children' in child:
            pathways_df = _traverse_pathway(child['children'], top_pathway, paths_found, pathways_df)    
    return pathways_df

def get_reactome_pathways(species: str = "Homo sapiens") -> pd.DataFrame:
    # Species implementation in reactome2py dev branch
    hierarchy = rc.event_species()
    pathways_df = pd.DataFrame(columns = ['Pathway','Stable ID','Top Level Pathway'])
    paths_found = set()

    for level in hierarchy:
        if level['type'] == 'TopLevelPathway':
            top_pathway = level['name']
            paths_found.add(level['stId'])
            pathways_df.loc[-1] = [level['name'], level['stId'], top_pathway]
            pathways_df.index = pathways_df.index + 1
            pathways_df = _traverse_pathway(level['children'], top_pathway, paths_found, pathways_df) 

        else:
            pass # reactome2py.content.event_species should only contain top-level paths
    
    pathways_df = pathways_df.reset_index(drop=True)
    return pathways_df
