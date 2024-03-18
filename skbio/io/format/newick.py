r"""Newick format (:mod:`skbio.io.format.newick`)
=============================================

.. currentmodule:: skbio.io.format.newick

Newick format (``newick``) stores spanning-trees with weighted edges and node
names in a minimal file format [1]_. This is useful for representing
phylogenetic trees and taxonomies. Newick was created as an informal
specification on June 26, 1986 [2]_.

Format Support
--------------
**Has Sniffer: Yes**

+------+------+---------------------------------------------------------------+
|Reader|Writer|                          Object Class                         |
+======+======+===============================================================+
|Yes   |Yes   |:mod:`skbio.tree.TreeNode`                                     |
+------+------+---------------------------------------------------------------+

Format Specification
--------------------
A Newick file represents a tree using the following grammar. See below for an
explanation of the format in plain English.

Formal Grammar
^^^^^^^^^^^^^^
.. code-block:: none

          NEWICK ==> NODE ;
            NODE ==> FORMATTING SUBTREE FORMATTING NODE_INFO FORMATTING
         SUBTREE ==> ( CHILDREN ) | null
       NODE_INFO ==> LABEL | LENGTH | LABEL FORMATTING LENGTH | null
      FORMATTING ==> [ COMMENT_CHARS ] | whitespace | null
        CHILDREN ==> NODE | CHILDREN , NODE
           LABEL ==> ' ALL_CHARS ' | SAFE_CHARS
          LENGTH ==> : FORMATTING NUMBER
   COMMENT_CHARS ==> any
       ALL_CHARS ==> any
      SAFE_CHARS ==> any except: ,;:()[] and whitespace
          NUMBER ==> a decimal or integer

.. note:: The ``_`` character inside of SAFE_CHARS will be converted to a
   blank space in ``skbio.tree.TreeNode`` and vice versa.

   ``'`` is considered the escape character. To escape ``'`` use a
   preceding ``'``.

   The implementation of newick in scikit-bio allows nested comments. To
   escape ``[`` or ``]`` from within COMMENT_CHARS, use a preceding ``'``.

Explanation
^^^^^^^^^^^
The Newick format defines a tree by creating a minimal representation of nodes
and their relationships to each other.

Basic Symbols
~~~~~~~~~~~~~
There are several symbols which define nodes, the first of which is the
semi-colon (``;``). The semi-colon creates a root node to its left. Recall that
there can only be one root in a tree.

The next symbol is the comma (``,``), which creates a node to its right.
However, these two alone are not enough. For example imagine the following
string: ``, , , ;``. It is evident that there is a root, but the other 3 nodes,
defined by commas, have no relationship. For this reason, it is not a valid
Newick string to have more than one node at the root level.

To provide these relationships, there is another structure:
paired parenthesis (``( )``). These are inserted at the location of an existing
node and give it the ability to have children. Placing ``( )`` in a node's
location will create a child inside the parenthesis on the left-most
inner edge.

Application of Rules
~~~~~~~~~~~~~~~~~~~~
Adding a comma within the parenthesis will create two children: ``( , )``
(also known as a bifurcating node). Notice that only one comma is needed
because the parenthesis have already created a child. Adding more commas will
create more children who are siblings to each other. For example, writing
``( , , , )`` will create a multifurcating node with 4 child nodes who are
siblings to each other.

The notation for a root can be used to create a complete tree. The ``;`` will
create a root node where parenthesis can be placed: ``( );``. Adding commas
will create more children: ``( , );``. These rules can be applied recursively
ad. infinitum: ``(( , ), ( , ));``.

Adding Node Information
~~~~~~~~~~~~~~~~~~~~~~~
Information about a node can be added to improve the clarity and meaning of a
tree. Each node may have a label and/or a length (to the parent). Newick always
places the node information at the right-most edge of a node's position.

Starting with labels, ``(( , ), ( , ));`` would become
``((D, E)B, (F, G)C)A;``. There is a named root ``A`` and the root's children
(from left to right) are ``B`` and ``C``. ``B`` has the children ``D`` and
``E``, and ``C`` has the children ``F`` and ``G``.

Length represents the distance (or weight of the edge) that connects a node to
its parent. This must be a decimal or integer. As an example, suppose ``D`` is
rather estranged from ``B``, and ``E`` is very close. That can be written as:
``((D:10, E:0.5)B, (F, G)C)A;``. Notice that the colon (``:``) separates the
label from the length. If the length is provided but the label is omitted, a
colon must still precede the length (``(:0.25,:0.5):0.0;``). Without this, the
length would be interpreted as a label (which happens to be a number).

.. note:: Internally scikit-bio will cast a length to ``float`` which
   technically means that even exponent strings (``1e-3``) are supported)

Advanced Label and Length Rules
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
More characters can be used to create more descriptive labels. When creating a
label there are some rules that must be considered due to limitations in the
Newick format. The following characters are not allowed within a standard
label: parenthesis, commas, square-brackets, colon, semi-colon, and whitespace.
These characters are also disallowed from occurring within a length, which has
a much stricter format: decimal or integer. Many of these characters are
symbols which define the structure of a Newick tree and are thus disallowed for
obvious reasons. The symbols not yet mentioned are square-brackets (``[ ]``)
and whitespace (space, tab, and newline).

What if these characters are needed within a label? In the simple case of
spaces, an underscore (``_``) will be translated as a space on read and vice
versa on write.

What if a literal underscore or any of the others mentioned are needed?
A label can be escaped (meaning that its contents are understood as regular
text) using single-quotes (``'``). When a label is surrounded by single-quotes,
any character is permissible. If a single-quote is needed inside of an escaped
label or anywhere else, it can be escaped with another single-quote.
For example, ``A_1`` is written ``'A_1'`` and ``'A'_1`` would be ``'''A''_1'``.

Inline Comments
~~~~~~~~~~~~~~~
Square-brackets define a comment, which are the least commonly used part of
the Newick format. Comments are not included in the generated objects and exist
only as human readable text ignored by the parser. The implementation in
scikit-bio allows for nested comments (``[comment [nested]]``). Unpaired
square-brackets can be escaped with a single-quote preceding the bracket when
inside an existing comment. (This is identical to escaping a single-quote).
The single-quote has the highest operator precedence, so there is no need to
worry about starting a comment from within a properly escaped label.

Whitespace
~~~~~~~~~~
Whitespace is not allowed within any un-escaped label or in any length, but it
is permitted anywhere else.

Caveats
~~~~~~~
Newick cannot always provide a unique representation of any tree, in other
words, the same tree can be written multiple ways. For example: ``(A, B);`` is
isomorphic to ``(B, A);``. The implementation in scikit-bio maintains the given
sibling order in its object representations.

Newick has no representation of an unrooted tree. Some biological packages make
the assumption that when a trifurcated root exists in an otherwise bifurcated
tree that the tree must be unrooted. In scikit-bio, ``skbio.tree.TreeNode``
will always be rooted at the ``newick`` root (``;``).

Format Parameters
-----------------
The only supported format parameter is `convert_underscores`. This is `True` by
default. When `False`, underscores found in unescaped labels will not be
converted to spaces. This is useful when reading the output of an external
program in which the underscores were not escaped. This parameter only affects
`read` operations. It does not exist for `write` operations; they will always
properly escape underscores.

Examples
--------
This is a simple Newick string.

>>> from io import StringIO
>>> from skbio import read
>>> from skbio.tree import TreeNode
>>> f = StringIO("((D, E)B, (F, G)C)A;")
>>> tree = read(f, format="newick", into=TreeNode)
>>> f.close()
>>> print(tree.ascii_art())
                    /-D
          /B-------|
         |          \-E
-A-------|
         |          /-F
          \C-------|
                    \-G

This is a complex Newick string.

>>> f = StringIO("[example](a:0.1, 'b_b''':0.2, (c:0.3, d_d:0.4)e:0.5)f:0.0;")
>>> tree = read(f, format="newick", into=TreeNode)
>>> f.close()
>>> print(tree.ascii_art())
          /-a
         |
-f-------|--b_b'
         |
         |          /-c
          \e-------|
                    \-d d

Notice that the node originally labeled ``d_d`` became ``d d``. Additionally
``'b_b'''`` became ``b_b'``. Note that the underscore was preserved in `b_b'`.

References
----------
.. [1] http://evolution.genetics.washington.edu/phylip/newick_doc.html
.. [2] http://evolution.genetics.washington.edu/phylip/newicktree.html


"""  # noqa: D205, D415

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from skbio.io import create_format, NewickFormatError
from skbio.tree import TreeNode

newick = create_format("newick")


@newick.sniffer()
def _newick_sniffer(fh):
    # Strategy:
    #   The following conditions preclude a file from being newick:
    #       * It is an empty file.
    #       * There is whitespace inside of a label (handled by tokenizer)
    #       * : is followed by anything that is an operator
    #       * ( is not preceded immediately by , or another (
    #       * The parens are unbalanced when ; is found.
    #   If 100 tokens (or less if EOF occurs earlier) then it is probably
    #   newick, or at least we can't prove it isn't.
    operators = set(",;:()")
    empty = True
    last_token = ","
    indent = 0
    try:
        # 100 tokens ought to be enough for anybody.
        for token, _ in zip(_tokenize_newick(fh), range(100)):
            if token not in operators:
                pass
            elif token == "," and last_token != ":" and indent > 0:
                pass
            elif token == ":" and last_token != ":":
                pass
            elif token == ";" and last_token != ":" and indent == 0:
                pass
            elif token == ")" and last_token != ":":
                indent -= 1
            elif token == "(" and (last_token == "(" or last_token == ","):
                indent += 1
            else:
                raise NewickFormatError()

            last_token = token
            empty = False

    except NewickFormatError:
        return False, {}
    return not empty, {}


@newick.reader(TreeNode)
def _newick_to_tree_node(fh, convert_underscores=True):
    tree_stack = []
    current_depth = 0
    last_token = ""
    next_is_distance = False
    root = TreeNode()
    tree_stack.append((root, current_depth))
    for token in _tokenize_newick(fh, convert_underscores=convert_underscores):
        # Check for a label
        if last_token not in "(,):":
            if not next_is_distance:
                tree_stack[-1][0].name = last_token if last_token else None
            else:
                next_is_distance = False
        # Check for a distance
        if token == ":":
            next_is_distance = True
        elif last_token == ":":
            try:
                tree_stack[-1][0].length = float(token)
            except ValueError:
                raise NewickFormatError(
                    "Could not read length as numeric type" ": %s." % token
                )

        elif token == "(":
            current_depth += 1
            tree_stack.append((TreeNode(), current_depth))
        elif token == ",":
            tree_stack.append((TreeNode(), current_depth))
        elif token == ")":
            if len(tree_stack) < 2:
                raise NewickFormatError(
                    "Could not parse file as newick." " Parenthesis are unbalanced."
                )
            children = []
            # Pop all nodes at this depth as they belong to the remaining
            # node on the top of the stack as children.
            while current_depth == tree_stack[-1][1]:
                node, _ = tree_stack.pop()
                children.insert(0, node)
            parent = tree_stack[-1][0]
            if parent.children:
                raise NewickFormatError(
                    "Could not parse file as newick." " Contains unnested children."
                )
            # This is much faster than TreeNode.extend
            for child in children:
                child.parent = parent
            parent.children = children
            current_depth -= 1
        elif token == ";":
            if len(tree_stack) == 1:
                return root
            break

        last_token = token

    raise NewickFormatError(
        "Could not parse file as newick."
        " `(Parenthesis)`, `'single-quotes'`,"
        " `[comments]` may be unbalanced, or tree may be"
        " missing its root."
    )


@newick.writer(TreeNode)
def _tree_node_to_newick(obj, fh):
    operators = set(",:_;()[]")
    current_depth = 0
    nodes_left = [(obj, 0)]
    while len(nodes_left) > 0:
        entry = nodes_left.pop()
        node, node_depth = entry
        if node.children and node_depth >= current_depth:
            fh.write("(")
            nodes_left.append(entry)
            nodes_left += ((child, node_depth + 1) for child in reversed(node.children))
            current_depth = node_depth + 1
        else:
            if node_depth < current_depth:
                fh.write(")")
                current_depth -= 1

            # Note we don't check for None because there is no way to represent
            # an empty string as a label in Newick. Therefore, both None and ''
            # are considered to be the absence of a label.
            label = node._node_label()
            if label:
                escaped = "%s" % label.replace("'", "''")
                if any(t in operators for t in label):
                    fh.write("'")
                    fh.write(escaped)
                    fh.write("'")
                else:
                    fh.write(escaped.replace(" ", "_"))
            if node.length is not None:
                fh.write(":")
                fh.write("%s" % node.length)
            if nodes_left and nodes_left[-1][1] == current_depth:
                fh.write(",")

    fh.write(";\n")


def _tokenize_newick(fh, convert_underscores=True):
    structure_tokens = set("(),;:")
    not_escaped = True
    label_start = False
    last_non_ws_char = ""
    last_char = ""
    comment_depth = 0
    metadata_buffer = []
    # Strategy:
    # We will iterate by character.
    # Comments in newick are defined as:
    # [This is a comment]
    # Nested comments are allowed.
    #
    # The following characters indicate structure:
    #      ( ) , ; :
    #
    # Whitespace is never allowed in a newick label, so an exception will be
    # thrown.
    #
    # We use ' to indicate a literal string. It has the highest precedence of
    # any operator.
    for line in fh:
        for character in line:
            # We will start by handling the comment case.
            # This code branch will probably never execute in practice.
            # Using a comment_depth we can handle nested comments.
            # Additionally if we are inside an escaped literal string, then
            # we don't want to consider it a comment.
            if character == "[" and not_escaped:
                # Sometimes we might not want to nest a comment, so we will use
                # our escape character. This is not explicitly mentioned in
                # any format specification, but seems like what a reasonable
                # person might do.
                if last_non_ws_char != "'" or comment_depth == 0:
                    # Once again, only advance our depth if [ has not been
                    # escaped inside our comment.
                    comment_depth += 1
            if comment_depth > 0:
                # Same as above, but in reverse
                if character == "]" and last_non_ws_char != "'":
                    comment_depth -= 1
                last_non_ws_char = character
                continue
            # We are not in a comment block if we are below here.

            # If we are inside of an escaped string literal, then ( ) , ; are
            # meaningless to the structure.
            # Otherwise, we are ready to submit our metadata token.
            if not_escaped and character in structure_tokens:
                label_start = False
                metadata = "".join(metadata_buffer)
                # If the following condition is True, then we must have just
                # closed a literal. We know this because last_non_ws_char is
                # either None or the last non-whitespace character.
                # last_non_ws_char is None when we have just escaped an escape
                # and at the first iteration.
                if last_non_ws_char == "'" or not convert_underscores:
                    # Make no modifications.
                    yield metadata
                elif metadata:
                    # Underscores are considered to be spaces when not in an
                    # escaped literal string.
                    yield metadata.replace("_", " ")
                # Clear our buffer for the next metadata token and yield our
                # current structure token.
                metadata_buffer = []
                yield character
            # We will now handle escaped string literals.
            # They are inconvenient because any character inside of them is
            # valid, especially whitespace.
            # We also need to allow ' to be escaped by '. e.g. '' -> '
            elif character == "'":
                not_escaped = not not_escaped
                label_start = True
                if last_non_ws_char == "'":
                    # We are escaping our escape, so it should be added to our
                    # metadata_buffer which will represent some future token.
                    metadata_buffer.append(character)
                    # We do not want a running chain of overcounts, so we need
                    # to clear the last character and continue iteration from
                    # the top. Without this, the following would happen:
                    # ''' ' -> '' <open literal>
                    # What we want is:
                    # ''' ' -> '<open literal> <close literal>
                    last_non_ws_char = ""
                    last_char = ""
                    continue

            elif not character.isspace() or not not_escaped:
                if label_start and last_char.isspace() and not_escaped:
                    raise NewickFormatError(
                        "Newick files cannot have"
                        " unescaped whitespace in their"
                        " labels."
                    )
                metadata_buffer.append(character)
                label_start = True

            # This is equivalent to an `else` however it prevents coverage from
            # mis-identifying the `continue` as uncalled because cpython will
            # optimize it to a jump that is slightly different from the normal
            # jump it would have done anyways.
            elif True:
                # Skip the last statement
                last_char = character
                continue

            last_char = character
            # This line is skipped in the following cases:
            #    * comment_depth > 0, i.e. we are in a comment.
            #    * We have just processed the sequence '' and we don't want
            #      the sequence ''' to result in ''.
            #    * We have encountered whitespace that is not properly escaped.
            last_non_ws_char = character
