"""
parsers for different go-annotation formats
"""
import pandas as pd
from skbio.io import create_format

_EXPERIMENTAL_EVIDENCE = ('EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'IEP')

goa = create_format('goa')
_GENE_ASSOCIATION_COLUMNS = ('db', 'db_object_id', 'db_object_symbol',
                            'qualifier', 'go_id', 'db_reference',
                            'evidence_code', 'with_from', 'aspect',
                            'db_object_name', 'db_object_synonym',
                            'db_object_type', 'taxon', 'date', 'assigned_by',
                            'annotation_extension', 'gene_product_form_id')
@goa.reader(pd.DataFrame, monkey_patch=False)
def _goa_to_pd_dataframe(filename, experimental=True, **kwds):
    """ read go-annotation file
    
    :param entry_id: protein or gene identifier column
    :param category_id: GO term column
    :param background: background annotation set
    """
    defaults = {'comment' : '!',
            'names': _GENE_ASSOCIATION_COLUMNS}

    if experimental and 'usecols' in kwds:
        kwds['usecols'] += ('evidence_code', )

    defaults.update(kwds)
    result = pd.read_table(filename, **defaults)

    if experimental:
        retain_mask = result.evidence_code.isin(_EXPERIMENTAL_EVIDENCE)
        result.drop(result.index[~retain_mask], inplace=True)

    return result


gene2go = create_format('gene2go', encoding='binary')
_GENE2GO_COLUMNS = ('tax_id', 'GeneID', 'GO_ID', 'Evidence', 'Qualifier', 'GO_term', 'PubMed', 'Category')
@gene2go.reader(pd.DataFrame, monkey_patch=False)
def gene2go(filename, experimental=False, tax_id=9606, **kwds):
    """ read go-annotation file
        
    :param entry_id: protein or gene identifier column
    :param category_id: GO term column
    :param background: background annotation set
    """
    defaults = { 'comment' : '#',
            'names' : _GENE2GO_COLUMNS }
    defaults.update(kwds)
    result = pd.read_table(filename, **defaults)
    
    retain_mask = result.tax_id == tax_id
    result.drop(result.index[~retain_mask], inplace=True)

    if experimental:
        retain_mask = result.Evidence.isin(_EXPERIMENTAL_EVIDENCE)
        result.drop(result.index[~retain_mask], inplace=True)

    return result
