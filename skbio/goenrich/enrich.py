import networkx as nx
import numpy as np
from scipy.stats import hypergeom
from statsmodels.stats.multitest import fdrcorrection

import goenrich.export

def analyze(O, query, background_attribute, **kwargs):
    """ run enrichment analysis for query

    >>> O = goenrich.obo.ontology('db/go-basic.obo')
    >>> gene2go = goenrich.read.gene2go('db/gene2go.gz')
    >>> values = {k: set(v) for k,v in gene2go.groupby('GO_ID')['GeneID']}
    >>> goenrich.enrich.propagate(O, values, 'gene2go')
    >>> df = goenrich.enrich.analyze(O, query, ...)

    :param O: Ontology graph after backgroud was set
    :param query: array like of ids
    :returns: pandas.DataFrame with results
    """
    options = {
            'show' : 'top20'
    }
    options.update(kwargs)
    _query = set(query)
    terms, nodes = zip(*O.nodes(data=True))
    M = len({x for n in nodes for x in n[background_attribute]}) # all ids used
    N = len(_query)
    ps, xs, ns = calculate_pvalues(nodes, _query, background_attribute,
            M, **options)
    qs, rejs = multiple_testing_correction(ps, **options)
    df = goenrich.export.to_frame(nodes, term=terms, q=qs, rejected=rejs,
            p=ps, x=xs, n=ns, M=M, N=N)
    if 'gvfile' in options:
        show = options['show']
        if show.startswith('top'):
            top = int(show.replace('top', ''))
            sig = df.sort('q').head(top)['term']
        else:
            raise NotImplementedError(show)
        G = induced_subgraph(O, sig)
        for term, node, q, x, n, rej in zip(terms, nodes, qs, xs, ns, rejs):
            if term in G:
                G.node[term].update({'name' : node['name'], 'x' : x,
                    'q' : q, 'n' : n, 'significant' : rej})
        G.reverse(copy=False)
        goenrich.export.to_graphviz(G, **options)
    return df
    
def propagate(O, values, attribute):
    """ Propagate values trough the hierarchy
    
    >>> O = goenrich.obo.ontology('db/go-basic.obo')
    >>> gene2go = goenrich.read.gene2go('db/gene2go.gz')
    >>> values = {k: set(v) for k,v in gene2go.groupby('GO_ID')['GeneID']}
    >>> goenrich.enrich.propagate(O, values, 'gene2go')

    Uses topological sorting of the vertices. Since degrees are
    usually low performance is almost linear time.

    :param O: ontology graph
    :param values: mapping of nodes to set of ids
    :param attribute: name of the attribute
    """
    for n in nx.topological_sort(O):
        current = O.node[n].setdefault(attribute, set())
        current.update(values.get(n, set()))
        for p in O[n]:
            O.node[p].setdefault(attribute, set()).update(current)

def induced_subgraph(O, terms):
    """  Extracts a subgraph from O including the provided terms
    and all higher hierarchy

    >>> df = goenrich.enrich.analyze(O, ...)
    >>> G = goenrich.induced_subgraph(O, df[df.rejected]['terms'])

    :param O: ontology graph
    :param terms: a list of terms to extract
    """
    roots = O.graph['roots'].values()
    nodes = set()
    for term in terms:
        if term in nodes:
            continue
        for root in roots:
            for path in nx.all_simple_paths(O, term, root):
                nodes.update(path)
    return O.subgraph(nodes)

def calculate_pvalues(nodes, query, background_attribute, M,
        min_category_size=3, max_category_size=500,
        max_category_depth=5, **kwargs):
    """ calculate pvalues for all categories in the graph
    
    :param G: ontology graph after background was set
    :param query: set of identifiers
    :param background_attribute: node attribute assoc. with the background set
    :param min_category_size: categories smaller than this number are ignored
    :param max_category_size: categories larger than this number are ignored
    :returns: pvalues, x, n
    """
    N = len(query)
    vals = []
    for node in nodes:
        background = node[background_attribute]
        n = len(background)
        hits = query.intersection(background)
        x = len(hits)
        if ((node['depth'] > max_category_depth)
            or (n <= min_category_size)
            or (n > max_category_size)):
            vals.append((float('NaN'), x, n))
        else:
            vals.append((hypergeom.sf(x, M, n, N), x, n))
    return zip(*vals)

def multiple_testing_correction(ps, alpha=0.05,
        method='benjamini-hochberg', **kwargs):
    """ correct pvalues for multiple testing and add corrected `q` value
    
    :param ps: list of pvalues
    :param alpha: significance level default : 0.05
    :param method: multiple testing correction method [bonferroni|benjamini-hochberg]
    :returns (q, rej): two lists of q-values and rejected nodes
    """
    _p = np.array(ps)
    q = _p.copy()
    rej = _p.copy()
    mask = ~np.isnan(_p)
    p = _p[mask]
    if method == 'bonferroni':
        q[mask] = p / len(p)
        rej[mask] = q[mask] < alpha
    elif method == 'benjamini-hochberg':
        _rej, _q = fdrcorrection(p, alpha)
        rej[mask] = _rej
        q[mask] = _q
    else:
        raise ValueError(method)
    return q, rej
