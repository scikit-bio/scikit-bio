# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import unittest

from skbio.util._docstring import (
    _note_into_doc,
    _note_into_doc_param,
)


class TestDocstring(unittest.TestCase):
    def setUp(self):
        self.docstring = """An example function.

        Here is an example function with a typical docstring.

        Parameters
        ----------
        param1 : int
            The first parameter.
        param2 : str, optional
            The second parameter.

        Returns
        -------
        float
            The return value.

        """

    def test_note_into_doc(self):
        note = ".. note:: This is a note.\n"

        # normal case
        obs = _note_into_doc(note, self.docstring)
        exp = """An example function.

        .. note:: This is a note.

        Here is an example function with a typical docstring.

        """
        self.assertIn(exp, obs)

        # no docstring
        obs = _note_into_doc(note, "")
        self.assertEqual(obs, note)

        # single-line docstring
        obs = _note_into_doc(note, "Hello!")
        self.assertEqual(obs, f"Hello!\n\n{note}")
        obs = _note_into_doc(note, "Hello!\n")
        self.assertEqual(obs, f"Hello!\n\n{note}")

        # single-paragraph docstring
        obs = _note_into_doc(note, "Hello!\nHello!\nHello!\n")
        self.assertEqual(obs, f"Hello!\nHello!\nHello!\n\n{note}")

        # single-paragraph docstring with indentation
        doc = """Hello! I am the
            header paragraph.
            """
        obs = _note_into_doc(note, doc)
        exp = """Hello! I am the
            header paragraph.

            .. note:: This is a note.
            """
        self.assertEqual(obs, exp)

    def test_note_into_doc_param(self):
        note = ".. note:: This is a note."

        # first parameter
        obs = _note_into_doc_param(note, self.docstring, "param1")
        exp = """
        param1 : int
            The first parameter.

            .. note:: This is a note.

        param2 : str, optional
        """
        self.assertIn(exp, obs)

        # last parameter
        obs = _note_into_doc_param(note, self.docstring, "param2")
        exp = """
        param2 : str, optional
            The second parameter.

            .. note:: This is a note.

        Returns
        """
        self.assertIn(exp, obs)

        # missing parameter
        with self.assertRaises(ValueError):
            _note_into_doc_param(note, self.docstring, "param3")

        # no more section
        doc = self.docstring[:self.docstring.find("Returns")]
        obs = _note_into_doc_param(note, doc, "param2")
        exp = """
        param2 : str, optional
            The second parameter.

            .. note:: This is a note.

        """
        self.assertTrue(obs.endswith(exp))


if __name__ == '__main__':
    unittest.main()
