import pandas as pd
from reactome2py import analysis
'''
From a list of genes, return a Pandas DataFrame with pathway enrichment results
'''
def pathway_enrichment(gene_query: list, projection: bool = True) -> pd.DataFrame:
    query = ",".join(gene_query)
    response = analysis.identifiers(ids=query, projection=projection)
    token = response['summary']['token']
    results = analysis.pathway2df(token)
    convert = {
        "#Entities found": "int64",
        "#Entities total": "int64",
        "Entities ratio": "float64",
        "Entities pValue": "float64",
        "Entities FDR": "float64",
        "#Reactions found": "int64",
        "#Reactions total": "int64",
        "Reactions ratio": "float64"
    }
    results = results.astype(convert)
    return results