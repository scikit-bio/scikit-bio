r"""
Ordination results format (:mod:`skbio.io.ordination`)
======================================================

.. currentmodule:: skbio.io.ordination

The ordination results file format (``ordination``) stores the results of an
ordination method in a human-readable, text-based format. The format supports
storing the results of various ordination methods available in scikit-bio,
including (but not necessarily limited to) PCoA, CA, RDA, and CCA.

Format Support
--------------
**Has Sniffer: Yes**

+----------+----------+------------------------------------------------------+
|**Reader**|**Writer**|                   **Object Class**                   |
+----------+----------+------------------------------------------------------+
|Yes       |Yes       |:mod:`skbio.stats.ordination.OrdinationResults`       |
+----------+----------+------------------------------------------------------+

Format Specification
--------------------
The format is text-based, consisting of six attributes that describe the
ordination results:

- ``Eigvals``: 1-D
- ``Proportion explained``: 1-D
- ``Species``: 2-D
- ``Site``: 2-D
- ``Biplot``: 2-D
- ``Site constraints``: 2-D

The attributes in the file *must* be in this order.

Each attribute is defined in its own section of the file, where sections are
separated by a blank (or whitespace-only) line. Each attribute begins with a
header line, which contains the attribute's name (as listed above), followed by
a tab character, followed by one or more tab-separated dimensions (integers)
that describe the shape of the attribute's data.

The attribute's data follows its header line, and is stored in tab-separated
format. ``Species``, ``Site``, and ``Site constraints`` store species and site
IDs, respectively, as the first column, followed by the 2-D data array.

An example of this file format might look like::

    Eigvals<tab>9
    0.366135830393<tab>0.186887643052<tab>0.0788466514249<tab>0.082287840501<tab>0.0351348475787<tab>0.0233265839374<tab>0.0099048981912<tab>0.00122461669234<tab>0.000417454724117

    Proportion explained<tab>9
    0.466910908242<tab>0.238326522327<tab>0.100548371868<tab>0.104936712434<tab>0.0448053488135<tab>0.0297469834643<tab>0.0126311183626<tab>0.00156167969536<tab>0.000532354793515

    Species<tab>9<tab>9
    Species0<tab>0.110350890177<tab>0.282399990052<tab>-0.203028976154<tab>-0.00192462284409<tab>-0.082232863384<tab>0.0857314258364<tab>-0.0122038907184<tab>-0.04251987936660.00466719926338
    Species1<tab>0.141359038961<tab>0.303495645402<tab>0.395441211576<tab>-0.14126625534<tab>-0.0268859204718<tab>0.143253061936<tab>0.0430260301697<tab>0.0476377655759<tab>-0.00228172378295
    Species2<tab>-1.01552204222<tab>0.0958317865043<tab>-0.198262718034<tab>-0.1048010300670.130025239749<tab>0.0244045261332<tab>0.0464701211285<tab>0.0269279200532<tab>0.0350103177576
    Species3<tab>-1.03620650502<tab>0.109624974112<tab>0.220984718362<tab>0.223640072997<tab>-0.243745876054<tab>-0.0259064859794<tab>-0.0534088909011<tab>-0.03156111959930.0256448427552
    Species4<tab>1.05371722248<tab>0.537178749104<tab>-0.438075060322<tab>0.223480553581<tab>-0.323948461806<tab>0.124644870822<tab>-0.119275907223<tab>0.0416254660329<tab>-0.0381955235096
    Species5<tab>0.998558655<tab>0.573960582723<tab>0.679918103399<tab>-0.3899633807170.299077945999<tab>0.328451006171<tab>0.21215881857<tab>-0.0829871883001<tab>-0.0439653996462
    Species6<tab>0.255245719879<tab>-0.178168259149<tab>-0.204127155429<tab>0.433397565801<tab>0.0707099230629<tab>-0.18817306522<tab>0.126908756045<tab>0.0044937289123<tab>-0.0122511718244
    Species7<tab>0.146555872394<tab>-0.857362497037<tab>-0.0152499051659<tab>0.0527604990862<tab>0.354475793915<tab>-0.0416813697787<tab>-0.199011239586<tab>-0.00213723187073<tab>-0.00782946141667
    Species8<tab>0.413705102117<tab>-0.707948964322<tab>0.21569736034<tab>-0.690314241725-0.148431001217<tab>-0.334251868558<tab>-0.00628707445028<tab>-0.00364123416731<tab>-0.0122722164511

    Site<tab>10<tab>9
    Site0<tab>0.710587311248<tab>-3.08166800613<tab>0.219651379947<tab>-1.24528801163<tab>-1.07293546227<tab>-0.506241907472<tab>0.244126652455<tab>-3.63164833508<tab>1.16311896657
    Site1<tab>0.584771352278<tab>-3.00669301091<tab>-0.947448656768<tab>2.69965142856<tab>2.13682885838<tab>0.813520011254<tab>0.471530333182<tab>0.908423015086<tab>-1.34724387844
    Site2<tab>0.762734278287<tab>-3.15258603503<tab>2.13924426714<tab>-3.1162748358<tab>-2.30660936925<tab>-0.698929858809<tab>-1.39062619586<tab>4.84117591747<tab>0.562102984837
    Site3<tab>1.11230735331<tab>1.07150585141<tab>-1.87527740873<tab>0.666370241998<tab>-1.10153224699<tab>1.43517552491<tab>-1.10619960297<tab>0.0137029328454<tab>-0.0371803939101
    Site4<tab>-0.979116769996<tab>-0.0603144289026<tab>-0.696277367656<tab>-0.6126467033080.983006619615<tab>0.315662442163<tab>0.574110232297<tab>0.328630035672<tab>0.868027697443
    Site5<tab>1.04322560423<tab>0.459426970165<tab>-0.639802790578<tab>0.287156643872<tab>-0.573935423877<tab>-1.44980634943<tab>1.70166994063<tab>0.306164261447<tab>-0.442115969758
    Site6<tab>-0.954490118162<tab>-0.0847021660539<tab>0.132509124411<tab>-0.42143341064<tab>-0.111552348931<tab>-0.394242454835<tab>-0.673963982894<tab>-0.379018566362<tab>-1.7472502885
    Site7<tab>0.947268764751<tab>-0.108370567311<tab>0.526107182587<tab>-0.00565282365567<tab>1.26272400228<tab>-1.06565692165<tab>-1.46326596729<tab>-0.154459216567<tab>0.778139732463
    Site8<tab>-1.14808173207<tab>0.490449274267<tab>0.478353666755<tab>1.17015870919<tab>-1.00599224074<tab>0.0735071441404<tab>0.0860462673715<tab>0.0417647558417<tab>0.935819560428
    Site9<tab>1.03291557934<tab>1.0350490304<tab>2.74691777314<tab>-1.28083971649<tab>0.363002636972<tab>1.98647950015<tab>1.05356145232<tab>-0.24813142226<tab>-0.463165215106

    Biplot<tab>3<tab>3
    -0.169746767979<tab>0.63069090084<tab>0.760769036049
    -0.994016563505<tab>0.0609533148724<tab>-0.0449369418179
    0.184352565909<tab>-0.974867543612<tab>0.0309865007541

    Site constraints<tab>10<tab>9
    Site0<tab>0.692138797603<tab>-3.08053663489<tab>-0.328747278055<tab>-1.24528801163<tab>-1.07293546227<tab>-0.506241907472<tab>0.244126652455<tab>-3.63164833508<tab>1.16311896657
    Site1<tab>0.664559513865<tab>-3.06214571808<tab>0.230249303805<tab>2.69965142856<tab>2.13682885838<tab>0.813520011254<tab>0.471530333182<tab>0.908423015086<tab>-1.34724387844
    Site2<tab>0.636980230127<tab>-3.04375480127<tab>0.789245885666<tab>-3.1162748358<tab>-2.30660936925<tab>-0.698929858809<tab>-1.39062619586<tab>4.84117591747<tab>0.562102984837
    Site3<tab>1.10887578995<tab>0.500396915484<tab>-1.55606822404<tab>0.666370241998<tab>-1.10153224699<tab>1.43517552491<tab>-1.10619960297<tab>0.0137029328454<tab>-0.0371803939101
    Site4<tab>-0.970016224052<tab>0.0654867737684<tab>-1.1206070781<tab>-0.612646703308<tab>0.983006619615<tab>0.315662442163<tab>0.574110232297<tab>0.328630035672<tab>0.868027697443
    Site5<tab>1.05371722248<tab>0.537178749104<tab>-0.438075060322<tab>0.287156643872<tab>-0.573935423877<tab>-1.44980634943<tab>1.70166994063<tab>0.306164261447<tab>-0.442115969758
    Site6<tab>-1.02517479153<tab>0.102268607388<tab>-0.00261391438256<tab>-0.42143341064<tab>-0.111552348931<tab>-0.394242454835<tab>-0.673963982894<tab>-0.379018566362<tab>-1.7472502885
    Site7<tab>0.998558655<tab>0.573960582723<tab>0.679918103399<tab>-0.00565282365567<tab>1.26272400228<tab>-1.06565692165<tab>-1.46326596729<tab>-0.154459216567<tab>0.778139732463
    Site8<tab>-1.080333359<tab>0.139050441007<tab>1.11537924934<tab>1.17015870919<tab>-1.00599224074<tab>0.0735071441404<tab>0.0860462673715<tab>0.0417647558417<tab>0.935819560428
    Site9<tab>0.943400087524<tab>0.610742416342<tab>1.79791126712<tab>-1.28083971649<tab>0.363002636972<tab>1.98647950015<tab>1.05356145232<tab>-0.24813142226<tab>-0.463165215106


If a given result attribute is not present (e.g. ``Biplot``), it should still
be defined and declare its dimensions as 0. For example::

    Biplot<tab>0<tab>0

All attributes are optional except for ``Eigvals``.

Examples
--------
Assume we have the following tab-delimited text file storing the
ordination results in ``ordination`` format::

    Eigvals<tab>9
    0.366135830393<tab>0.186887643052<tab>0.0788466514249<tab>0.082287840501<tab>0.0351348475787<tab>0.0233265839374<tab>0.0099048981912<tab>0.00122461669234<tab>0.000417454724117

    Proportion explained<tab>9
    0.466910908242<tab>0.238326522327<tab>0.100548371868<tab>0.104936712434<tab>0.0448053488135<tab>0.0297469834643<tab>0.0126311183626<tab>0.00156167969536<tab>0.000532354793515

    Species<tab>9<tab>9
    Species0<tab>0.110350890177<tab>0.282399990052<tab>-0.203028976154<tab>-0.00192462284409<tab>-0.082232863384<tab>0.0857314258364<tab>-0.0122038907184<tab>-0.04251987936660.00466719926338
    Species1<tab>0.141359038961<tab>0.303495645402<tab>0.395441211576<tab>-0.14126625534<tab>-0.0268859204718<tab>0.143253061936<tab>0.0430260301697<tab>0.0476377655759<tab>-0.00228172378295
    Species2<tab>-1.01552204222<tab>0.0958317865043<tab>-0.198262718034<tab>-0.1048010300670.130025239749<tab>0.0244045261332<tab>0.0464701211285<tab>0.0269279200532<tab>0.0350103177576
    Species3<tab>-1.03620650502<tab>0.109624974112<tab>0.220984718362<tab>0.223640072997<tab>-0.243745876054<tab>-0.0259064859794<tab>-0.0534088909011<tab>-0.03156111959930.0256448427552
    Species4<tab>1.05371722248<tab>0.537178749104<tab>-0.438075060322<tab>0.223480553581<tab>-0.323948461806<tab>0.124644870822<tab>-0.119275907223<tab>0.0416254660329<tab>-0.0381955235096
    Species5<tab>0.998558655<tab>0.573960582723<tab>0.679918103399<tab>-0.3899633807170.299077945999<tab>0.328451006171<tab>0.21215881857<tab>-0.0829871883001<tab>-0.0439653996462
    Species6<tab>0.255245719879<tab>-0.178168259149<tab>-0.204127155429<tab>0.433397565801<tab>0.0707099230629<tab>-0.18817306522<tab>0.126908756045<tab>0.0044937289123<tab>-0.0122511718244
    Species7<tab>0.146555872394<tab>-0.857362497037<tab>-0.0152499051659<tab>0.0527604990862<tab>0.354475793915<tab>-0.0416813697787<tab>-0.199011239586<tab>-0.00213723187073<tab>-0.00782946141667
    Species8<tab>0.413705102117<tab>-0.707948964322<tab>0.21569736034<tab>-0.690314241725-0.148431001217<tab>-0.334251868558<tab>-0.00628707445028<tab>-0.00364123416731<tab>-0.0122722164511

    Site<tab>10<tab>9
    Site0<tab>0.710587311248<tab>-3.08166800613<tab>0.219651379947<tab>-1.24528801163<tab>-1.07293546227<tab>-0.506241907472<tab>0.244126652455<tab>-3.63164833508<tab>1.16311896657
    Site1<tab>0.584771352278<tab>-3.00669301091<tab>-0.947448656768<tab>2.69965142856<tab>2.13682885838<tab>0.813520011254<tab>0.471530333182<tab>0.908423015086<tab>-1.34724387844
    Site2<tab>0.762734278287<tab>-3.15258603503<tab>2.13924426714<tab>-3.1162748358<tab>-2.30660936925<tab>-0.698929858809<tab>-1.39062619586<tab>4.84117591747<tab>0.562102984837
    Site3<tab>1.11230735331<tab>1.07150585141<tab>-1.87527740873<tab>0.666370241998<tab>-1.10153224699<tab>1.43517552491<tab>-1.10619960297<tab>0.0137029328454<tab>-0.0371803939101
    Site4<tab>-0.979116769996<tab>-0.0603144289026<tab>-0.696277367656<tab>-0.6126467033080.983006619615<tab>0.315662442163<tab>0.574110232297<tab>0.328630035672<tab>0.868027697443
    Site5<tab>1.04322560423<tab>0.459426970165<tab>-0.639802790578<tab>0.287156643872<tab>-0.573935423877<tab>-1.44980634943<tab>1.70166994063<tab>0.306164261447<tab>-0.442115969758
    Site6<tab>-0.954490118162<tab>-0.0847021660539<tab>0.132509124411<tab>-0.42143341064<tab>-0.111552348931<tab>-0.394242454835<tab>-0.673963982894<tab>-0.379018566362<tab>-1.7472502885
    Site7<tab>0.947268764751<tab>-0.108370567311<tab>0.526107182587<tab>-0.00565282365567<tab>1.26272400228<tab>-1.06565692165<tab>-1.46326596729<tab>-0.154459216567<tab>0.778139732463
    Site8<tab>-1.14808173207<tab>0.490449274267<tab>0.478353666755<tab>1.17015870919<tab>-1.00599224074<tab>0.0735071441404<tab>0.0860462673715<tab>0.0417647558417<tab>0.935819560428
    Site9<tab>1.03291557934<tab>1.0350490304<tab>2.74691777314<tab>-1.28083971649<tab>0.363002636972<tab>1.98647950015<tab>1.05356145232<tab>-0.24813142226<tab>-0.463165215106

    Biplot<tab>0<tab>0

    Site constraints<tab>0<tab>0

Load the ordination results from the file:

>>> from StringIO import StringIO
>>> from skbio.stats.ordination import OrdinationResults
>>> or_f = StringIO(
...  "Eigvals\t9\n"
...  "0.366135830393\t0.186887643052\t0.0788466514249\t0.082287840501\t0.0351348475787\t0.0233265839374\t0.0099048981912\t0.00122461669234\t0.000417454724117\n"
...  "\n"
...  "Proportion explained\t9\n"
...  "0.466910908242\t0.238326522327\t0.100548371868\t0.104936712434\t0.0448053488135\t0.0297469834643\t0.0126311183626\t0.00156167969536\t0.000532354793515\n"
...  "\n"
...  "Species\t9\t9\n"
...  "Species0\t0.110350890177\t0.282399990052\t-0.203028976154\t-0.00192462284409\t-0.082232863384\t0.0857314258364\t-0.0122038907184\t-0.0425198793666\t0.00466719926338\n"
...  "Species1\t0.141359038961\t0.303495645402\t0.395441211576\t-0.14126625534\t-0.0268859204718\t0.143253061936\t0.0430260301697\t0.0476377655759\t-0.00228172378295\n"
...  "Species2\t-1.01552204222\t0.0958317865043\t-0.198262718034\t-0.104801030067\t0.130025239749\t0.0244045261332\t0.0464701211285\t0.0269279200532\t0.0350103177576\n"
...  "Species3\t-1.03620650502\t0.109624974112\t0.220984718362\t0.223640072997\t-0.243745876054\t-0.0259064859794\t-0.0534088909011\t-0.0315611195993\t0.0256448427552\n"
...  "Species4\t1.05371722248\t0.537178749104\t-0.438075060322\t0.223480553581\t-0.323948461806\t0.124644870822\t-0.119275907223\t0.0416254660329\t-0.0381955235096\n"
...  "Species5\t0.998558655\t0.573960582723\t0.679918103399\t-0.389963380717\t0.299077945999\t0.328451006171\t0.21215881857\t-0.0829871883001\t-0.0439653996462\n"
...  "Species6\t0.255245719879\t-0.178168259149\t-0.204127155429\t0.433397565801\t0.0707099230629\t-0.18817306522\t0.126908756045\t0.0044937289123\t-0.0122511718244\n"
...  "Species7\t0.146555872394\t-0.857362497037\t-0.0152499051659\t0.0527604990862\t0.354475793915\t-0.0416813697787\t-0.199011239586\t-0.00213723187073\t-0.00782946141667\n"
...  "Species8\t0.413705102117\t-0.707948964322\t0.21569736034\t-0.690314241725\t-0.148431001217\t-0.334251868558\t-0.00628707445028\t-0.00364123416731\t-0.0122722164511\n"
...  "\n"
...  "Site\t10\t9\n"
...  "Site0\t0.710587311248\t-3.08166800613\t0.219651379947\t-1.24528801163\t-1.07293546227\t-0.506241907472\t0.244126652455\t-3.63164833508\t1.16311896657\n"
...  "Site1\t0.584771352278\t-3.00669301091\t-0.947448656768\t2.69965142856\t2.13682885838\t0.813520011254\t0.471530333182\t0.908423015086\t-1.34724387844\n"
...  "Site2\t0.762734278287\t-3.15258603503\t2.13924426714\t-3.1162748358\t-2.30660936925\t-0.698929858809\t-1.39062619586\t4.84117591747\t0.562102984837\n"
...  "Site3\t1.11230735331\t1.07150585141\t-1.87527740873\t0.666370241998\t-1.10153224699\t1.43517552491\t-1.10619960297\t0.0137029328454\t-0.0371803939101\n"
...  "Site4\t-0.979116769996\t-0.0603144289026\t-0.696277367656\t-0.612646703308\t0.983006619615\t0.315662442163\t0.574110232297\t0.328630035672\t0.868027697443\n"
...  "Site5\t1.04322560423\t0.459426970165\t-0.639802790578\t0.287156643872\t-0.573935423877\t-1.44980634943\t1.70166994063\t0.306164261447\t-0.442115969758\n"
...  "Site6\t-0.954490118162\t-0.0847021660539\t0.132509124411\t-0.42143341064\t-0.111552348931\t-0.394242454835\t-0.673963982894\t-0.379018566362\t-1.7472502885\n"
...  "Site7\t0.947268764751\t-0.108370567311\t0.526107182587\t-0.00565282365567\t1.26272400228\t-1.06565692165\t-1.46326596729\t-0.154459216567\t0.778139732463\n"
...  "Site8\t-1.14808173207\t0.490449274267\t0.478353666755\t1.17015870919\t-1.00599224074\t0.0735071441404\t0.0860462673715\t0.0417647558417\t0.935819560428\n"
...  "Site9\t1.03291557934\t1.0350490304\t2.74691777314\t-1.28083971649\t0.363002636972\t1.98647950015\t1.05356145232\t-0.24813142226\t-0.463165215106\n"
...  "\n"
...  "Biplot\t0\t0\n"
...  "\n"
...  "Site constraints\t0\t0\n")
>>> ord_res = OrdinationResults.read(or_f)

"""
# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function
from future.builtins import zip

import numpy as np

from skbio.stats.ordination import OrdinationResults
from skbio.io import (register_reader, register_writer, register_sniffer,
                      OrdinationFormatError)


@register_sniffer('ordination')
def _ordination_sniffer(fh):
    # Smells an ordination file if *all* of the following lines are present
    # *from the beginning* of the file:
    #   - eigvals header (minimally parsed)
    #   - another line (contents ignored)
    #   - a whitespace-only line
    #   - proportion explained header (minimally parsed)
    try:
        _parse_header(fh, 'Eigvals', 1)
        next_line = next(fh, None)

        if next_line is not None:
            _check_empty_line(fh)
            _parse_header(fh, 'Proportion explained', 1)
            return True, {}
    except OrdinationFormatError:
        pass

    return False, {}


@register_reader('ordination', OrdinationResults)
def _ordination_to_ordination_results(fh):
    eigvals = _parse_vector_section(fh, 'Eigvals')
    if eigvals is None:
        raise OrdinationFormatError("At least one eigval must be present.")
    _check_empty_line(fh)

    prop_expl = _parse_vector_section(fh, 'Proportion explained')
    _check_length_against_eigvals(prop_expl, eigvals,
                                  'proportion explained values')
    _check_empty_line(fh)

    species, species_ids = _parse_array_section(fh, 'Species')
    _check_length_against_eigvals(species, eigvals,
                                  'coordinates per species')
    _check_empty_line(fh)

    site, site_ids = _parse_array_section(fh, 'Site')
    _check_length_against_eigvals(site, eigvals,
                                  'coordinates per site')
    _check_empty_line(fh)

    # biplot does not have ids to parse (the other arrays do)
    biplot, _ = _parse_array_section(fh, 'Biplot', has_ids=False)
    _check_empty_line(fh)

    cons, cons_ids = _parse_array_section(fh, 'Site constraints')

    if cons_ids is not None and site_ids is not None:
        if cons_ids != site_ids:
            raise OrdinationFormatError(
                "Site constraints ids and site ids must be equal: %s != %s" %
                (cons_ids, site_ids))

    return OrdinationResults(
        eigvals=eigvals, species=species, site=site, biplot=biplot,
        site_constraints=cons, proportion_explained=prop_expl,
        species_ids=species_ids, site_ids=site_ids)


def _parse_header(fh, header_id, num_dimensions):
    line = next(fh, None)
    if line is None:
        raise OrdinationFormatError(
            "Reached end of file while looking for %s header." % header_id)

    header = line.strip().split('\t')
    # +1 for the header ID
    if len(header) != num_dimensions + 1 or header[0] != header_id:
        raise OrdinationFormatError("%s header not found." % header_id)
    return header


def _check_empty_line(fh):
    """Check that the next line in `fh` is empty or whitespace-only."""
    line = next(fh, None)
    if line is None:
        raise OrdinationFormatError(
            "Reached end of file while looking for blank line separating "
            "sections.")

    if line.strip():
        raise OrdinationFormatError("Expected an empty line.")


def _check_length_against_eigvals(data, eigvals, label):
    if data is not None:
        num_vals = data.shape[-1]
        num_eigvals = eigvals.shape[-1]

        if num_vals != num_eigvals:
            raise OrdinationFormatError(
                "There should be as many %s as eigvals: %d != %d" %
                (label, num_vals, num_eigvals))


def _parse_vector_section(fh, header_id):
    header = _parse_header(fh, header_id, 1)

    # Parse how many values we are waiting for
    num_vals = int(header[1])
    if num_vals == 0:
        # The ordination method didn't generate the vector, so set it to None
        vals = None
    else:
        # Parse the line with the vector values
        line = next(fh, None)
        if line is None:
            raise OrdinationFormatError(
                "Reached end of file while looking for line containing values "
                "for %s section." % header_id)
        vals = np.asarray(line.strip().split('\t'), dtype=np.float64)
        if len(vals) != num_vals:
            raise OrdinationFormatError(
                "Expected %d values in %s section, but found %d." %
                (num_vals, header_id, len(vals)))
    return vals


def _parse_array_section(fh, header_id, has_ids=True):
    """Parse an array section of `fh` identified by `header_id`."""
    # Parse the array header
    header = _parse_header(fh, header_id, 2)

    # Parse the dimensions of the array
    rows = int(header[1])
    cols = int(header[2])

    ids = None
    if rows == 0 and cols == 0:
        # The ordination method didn't generate the array data for 'header', so
        # set it to None
        data = None
    elif rows == 0 or cols == 0:
        # Both dimensions should be 0 or none of them are zero
        raise OrdinationFormatError("One dimension of %s is 0: %d x %d" %
                                    (header_id, rows, cols))
    else:
        # Parse the data
        data = np.empty((rows, cols), dtype=np.float64)

        if has_ids:
            ids = []

        for i in range(rows):
            # Parse the next row of data
            line = next(fh, None)
            if line is None:
                raise OrdinationFormatError(
                    "Reached end of file while looking for row %d in %s "
                    "section." % (i + 1, header_id))
            vals = line.strip().split('\t')

            if has_ids:
                ids.append(vals[0])
                vals = vals[1:]

            if len(vals) != cols:
                raise OrdinationFormatError(
                    "Expected %d values, but found %d in row %d." %
                    (cols, len(vals), i + 1))
            data[i, :] = np.asarray(vals, dtype=np.float64)
    return data, ids


@register_writer('ordination', OrdinationResults)
def _ordination_results_to_ordination(obj, fh):
    _write_vector_section(fh, 'Eigvals', obj.eigvals)
    _write_vector_section(fh, 'Proportion explained', obj.proportion_explained)
    _write_array_section(fh, 'Species', obj.species, obj.species_ids)
    _write_array_section(fh, 'Site', obj.site, obj.site_ids)
    _write_array_section(fh, 'Biplot', obj.biplot)
    _write_array_section(fh, 'Site constraints', obj.site_constraints,
                         obj.site_ids, include_section_separator=False)


def _write_vector_section(fh, header_id, vector):
    if vector is None:
        shape = 0
    else:
        shape = vector.shape[0]
    fh.write("%s\t%d\n" % (header_id, shape))

    if vector is not None:
        fh.write(_format_vector(vector))
    fh.write("\n")


def _write_array_section(fh, header_id, data, ids=None,
                         include_section_separator=True):
    # write section header
    if data is None:
        shape = (0, 0)
    else:
        shape = data.shape
    fh.write("%s\t%d\t%d\n" % (header_id, shape[0], shape[1]))

    # write section data
    if data is not None:
        if ids is None:
            for vals in data:
                fh.write(_format_vector(vals))
        else:
            for id_, vals in zip(ids, data):
                fh.write(_format_vector(vals, id_))

    if include_section_separator:
        fh.write("\n")


def _format_vector(vector, id_=None):
    formatted_vector = '\t'.join(np.asarray(vector, dtype=np.str))

    if id_ is None:
        return "%s\n" % formatted_vector
    else:
        return "%s\t%s\n" % (id_, formatted_vector)
