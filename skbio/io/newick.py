from skbio.io import register_reader, register_writer, register_sniffer
from skbio.tree import TreeNode
# Newick format, an informal standard agreed to in 1986
def _tokenize_newick(fh):
    """Produces valid tokens in the set of ['(', ')', ';', ',', <metadata>]."""
    not_escaped = True
    last_character = None
    comment_depth = 0
    metadata_buffer = []
    # Strategy:
    # We will iterate by character.
    # Comments in newick are defined as:
    # [This is a comment]
    # Nested comments are allowed.
    #
    # The following characters indicate structure:
    #      ( ) , ; : '
    #
    # Whitespace is never allowed in a newick label, so we will ignore all
    # instances of whitespace unless properly escaped by ''.
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
                if last_character != "'" or comment_depth == 0:
                    # Once again, only advance our depth if [ has not been
                    # escaped inside our comment.
                    comment_depth += 1
            if comment_depth > 0:
                # Same as above, but in reverse
                if character == "]" and last_character != "'":
                    comment_depth -= 1
                last_character = character
                continue
            # We are not in a comment block if we are below here.

            # If we are inside of an escaped string literal, then ( ) , ; are
            # meaningless to the structure.
            # Otherwise, we are ready to submit our metadata token.
            if not_escaped and (character == "(" or
                                character == ")" or
                                character == ":" or
                                character == "," or
                                character == ";"):
                metadata = ''.join(metadata_buffer)
                # If the following condition is True, then we must have just
                # closed a literal. We know this because last_character is
                # either None or the last non-whitespace character.
                # last_character is None when we have just escaped an escape
                # and at the first iteration.
                if last_character == "'":
                    # Make no modifications.
                    yield metadata
                elif metadata:
                    # Underscores are considered to be spaces when not in a
                    # escaped literal string.
                    # At first glance it may appear that this could clobber
                    # escaped literals, however a label must be entirely quoted
                    # if it is a literal. This means concatenation of literal
                    # and standard label buffers is illegal.
                    yield metadata.replace('_', ' ')
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
                if last_character == "'":
                    # We are escaping our escape, so it should be added to our
                    # metadata_buffer which will represent some future token.
                    metadata_buffer.append(character)
                    # We do not want a running chain of overcounts, so we need
                    # to clear the last character and continue iteration from
                    # the top. Without this, the following would happen:
                    # ''' ' -> '' <open literal>
                    # What we want is:
                    # ''' ' -> '<open literal> <close literal>
                    last_character = None
                    continue

            elif not character.isspace() or not not_escaped:
                metadata_buffer.append(character)

            # This is equivalent to an `else` however it prevents coverage from
            # mis-identifying the `continue` as uncalled because cpython will
            # optimize it to a jump that is slightly different from the normal
            # jump it would have done anyways.
            elif True:
                continue

            # This line is skipped in the following cases:
            #    * comment_depth > 0, i.e. we are in a comment.
            #    * We have just processed the sequence '' and we don't want
            #      the sequence ''' to result in ''.
            #    * We have encountered whitespace that is not properly escaped.
            last_character = character


@register_reader('newick', TreeNode)
def _newick_to_tree_node(fh):
    tree_stack = []
    current_depth = 0
    last_token = ''
    next_is_distance = False
    root = TreeNode()
    tree_stack.append((root, current_depth))
    for token in _tokenize_newick(fh):
        # Check for a label
        if last_token not in '(,):':
            if not next_is_distance:
                tree_stack[-1][0].name = last_token if last_token else None
            else:
                next_is_distance = False
        # Check for a distance
        if token == ':':
            next_is_distance = True
        elif last_token == ':':
            tree_stack[-1][0].length = float(token)

        elif token == '(':
            current_depth += 1
            tree_stack.append((TreeNode(), current_depth))
        elif token == ',':
            tree_stack.append((TreeNode(), current_depth))
        elif token == ')':
            children = []
            # Pop all nodes at this depth as they belong to the remaining
            # node on the top of the stack as children.
            while current_depth == tree_stack[-1][1]:
                node, _ = tree_stack.pop()
                children.append(node)
            parent = tree_stack[-1][0]
            # This is much faster than TreeNode.extend
            for child in children:
                child.parent = parent
            tree_stack[-1][0].children = children
            current_depth -= 1
        elif token == ';':
            break

        last_token = token

    if len(tree_stack) == 1:
        return root

    raise Exception()


@register_writer("newick", TreeNode)
def _tree_node_to_newick(obj, fh):
    current_node = obj
    tree_stack = []


@register_sniffer("newick")
def _newick_sniffer(fh):
    pass
