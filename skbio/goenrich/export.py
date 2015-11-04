r"""
goenrich export (:mod:`skbio.goenrich.export`)
==============================================

.. currentmodule:: skbio.goenrich.export

This module provides utility functions to deal with files and I/O in
general.asdf

Functions
---------

.. autosummary::
    :toctree:

    to_frame
    to_graphviz
"""
import networkx as nx
import numpy as np
import pandas as pd

import skbio.goenrich

def to_frame(nodes, **kwargs):
    """ export node attributes and key-values pairs to pd.DataFrame

    >>> to_frame(nodes, term = terms, pvalues = ps, ...)

    Parameters
    ----------
    nodes :  list of dict
        list of dictionaries with node attributes
    **kwargs :
        additional columns

    Returns
    -------
    df : pd.DataFrame
    """
    names, namespaces = zip(*[(n['name'], n['namespace']) for n in nodes])
    kwargs.update({'name' : names, 'namespace' : namespaces})
    return pd.DataFrame(kwargs)

def to_graphviz(G, gvfile, graph_label='', **kwargs):
    """ export graph of signifcant findings to dot file.
    A png can be generated from the commandline using graphviz

    >>> import subprocess
    >>> subprocess.call(['dot', '-Tpng', 'filpath.dot', '>', 'filepath.png'])

    Parameters
    ----------
    G : nx.Graph
        the graph to be exported
    gvfile : filehandle, filepath
        output path
    graph_label : str
        For an empty label pass graph_label=''.
    """
    for n in G:
        node = G.node[n]
        attr = {}
        attr['shape'] = 'record'
        if not np.isnan(node.get('q', float('NaN'))):
            attr['color'] = 'red' if node['significant'] else 'black'
            attr['label'] = "{name}\\n{x} / {n} genes\\nq = {q:E}".format(name=node['name'],
                    q=node['q'], x=node['x'], n=node['n'])
        else:
            attr['color'] = 'black'
            attr['label'] = """{name}""".format(name=node['name'])
        G.node[n] = attr

    A = nx.to_agraph(G)
    A.graph_attr['label'] = graph_label
    A.graph_attr['labelloc'] = 't'
    
    if hasattr(gvfile, 'write'):
        A.write(gvfile)
    else:
        with open(gvfile, 'w') as f:
            A.write(f)
