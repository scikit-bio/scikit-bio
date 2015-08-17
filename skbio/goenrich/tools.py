"""
helper functions for `goenrich`
"""
import pandas as pd

def generate_background(annotations, df, go_id, entry_id):
    """ generate the backgound from pandas datasets

    >>> O = ontology(...)
    >>> annotations = goenrich.read.gene2go(...)
    >>> background = generate_background(annotations, df, 'GO_ID', 'GeneID')
    >>> propagate(O, background, ...)

    :param annotations: pd.DataFrame containing the annotations
    :param df: pd.DataFrame containing the background genes
    :param go_id: GO id column name
    :param entry_id: Gene id column name
    :returns: dictionary with the background annotations
    """
    return {k : set(v) for k,v in
            pd.merge(annotations, df[[entry_id]]).groupby(go_id)[entry_id]}
