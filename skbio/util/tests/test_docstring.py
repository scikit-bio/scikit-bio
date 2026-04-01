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
    _array_api_compat_section,
    _insert_into_notes_section,
    _deprecation_note,
    _renaming_note,
)


class TestNoteIntoDoc(unittest.TestCase):
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

    def test_normal_case(self):
        note = ".. note:: This is a note.\n"
        obs = _note_into_doc(note, self.docstring)
        exp = """An example function.

        .. note:: This is a note.

        Here is an example function with a typical docstring.

        """
        self.assertIn(exp, obs)

    def test_no_docstring(self):
        note = ".. note:: This is a note.\n"
        obs = _note_into_doc(note, "")
        self.assertEqual(obs, note)

    def test_single_line_docstring(self):
        note = ".. note:: This is a note.\n"
        obs = _note_into_doc(note, "Hello!")
        self.assertEqual(obs, f"Hello!\n\n{note}")

    def test_single_line_with_trailing_newline(self):
        note = ".. note:: This is a note.\n"
        obs = _note_into_doc(note, "Hello!\n")
        self.assertEqual(obs, f"Hello!\n\n{note}")

    def test_single_paragraph(self):
        note = ".. note:: This is a note.\n"
        obs = _note_into_doc(note, "Hello!\nHello!\nHello!\n")
        self.assertEqual(obs, f"Hello!\nHello!\nHello!\n\n{note}")

    def test_single_paragraph_with_indentation(self):
        note = ".. note:: This is a note.\n"
        doc = """Hello! I am the
            header paragraph.
            """
        obs = _note_into_doc(note, doc)
        exp = """Hello! I am the
            header paragraph.

            .. note:: This is a note.
            """
        self.assertEqual(obs, exp)

    def test_note_trailing_whitespace_stripped(self):
        """Note with extra trailing whitespace should be normalized."""
        note = ".. note:: Trailing.   \n\n"
        obs = _note_into_doc(note, "")
        self.assertEqual(obs, ".. note:: Trailing.\n")


class TestNoteIntoDocParam(unittest.TestCase):
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

    def test_first_parameter(self):
        note = ".. note:: This is a note."
        obs = _note_into_doc_param(note, self.docstring, "param1")
        exp = """
        param1 : int
            The first parameter.

            .. note:: This is a note.

        param2 : str, optional
        """
        self.assertIn(exp, obs)

    def test_last_parameter(self):
        note = ".. note:: This is a note."
        obs = _note_into_doc_param(note, self.docstring, "param2")
        exp = """
        param2 : str, optional
            The second parameter.

            .. note:: This is a note.

        Returns
        """
        self.assertIn(exp, obs)

    def test_missing_parameter(self):
        note = ".. note:: This is a note."
        with self.assertRaises(ValueError):
            _note_into_doc_param(note, self.docstring, "param3")

    def test_no_section_after_param(self):
        note = ".. note:: This is a note."
        doc = self.docstring[:self.docstring.find("Returns")]
        obs = _note_into_doc_param(note, doc, "param2")
        exp = """
        param2 : str, optional
            The second parameter.

            .. note:: This is a note.

        """
        self.assertTrue(obs.endswith(exp))


class TestArrayApiCompatSection(unittest.TestCase):
    def test_single_backend_default_devices(self):
        result = _array_api_compat_section(["numpy"])
        self.assertIn("array API standard", result)
        self.assertIn("NumPy", result)
        self.assertIn("CPU", result)
        self.assertIn("\u2713", result)

    def test_multiple_backends_default_devices(self):
        result = _array_api_compat_section(["numpy", "torch", "jax"])
        self.assertIn("NumPy", result)
        self.assertIn("PyTorch", result)
        self.assertIn("JAX", result)
        self.assertIn("CPU", result)
        self.assertIn("GPU", result)
        # NumPy doesn't support GPU
        self.assertIn("n/a", result)

    def test_explicit_devices(self):
        result = _array_api_compat_section(["numpy", "torch"], devices=["cpu"])
        self.assertIn("CPU", result)
        # GPU column should not appear
        self.assertNotIn("GPU", result)

    def test_unknown_backend_fallback(self):
        """Unknown backend should capitalize its name and default to cpu."""
        result = _array_api_compat_section(["mylib"])
        self.assertIn("Mylib", result)
        self.assertIn("CPU", result)
        self.assertIn("\u2713", result)

    def test_cupy_backend(self):
        result = _array_api_compat_section(["cupy"])
        self.assertIn("CuPy", result)
        self.assertIn("GPU", result)
        self.assertIn("\u2713", result)

    def test_dask_backend(self):
        result = _array_api_compat_section(["dask"])
        self.assertIn("Dask", result)
        self.assertIn("CPU", result)

    def test_all_known_backends(self):
        result = _array_api_compat_section(
            ["numpy", "cupy", "torch", "jax", "dask"])
        for name in ("NumPy", "CuPy", "PyTorch", "JAX", "Dask"):
            self.assertIn(name, result)

    def test_table_structure(self):
        """Table should have proper RST grid table delimiters."""
        result = _array_api_compat_section(["numpy"])
        lines = result.strip().splitlines()
        # Table lines start with + or |
        table_lines = [l for l in lines if l.startswith("+") or l.startswith("|")]
        self.assertTrue(len(table_lines) >= 5)
        # Header separator uses =
        header_seps = [l for l in table_lines if "=" in l]
        self.assertTrue(len(header_seps) >= 1)

    def test_explicit_devices_filters_supported(self):
        """When devices are explicit, only those devices the backend supports appear."""
        result = _array_api_compat_section(["numpy"], devices=["cpu", "gpu"])
        # NumPy only supports cpu, so gpu should show n/a
        self.assertIn("n/a", result)


class TestInsertIntoNotesSection(unittest.TestCase):
    def test_existing_notes_section(self):
        doc = """My function.

        Parameters
        ----------
        x : int
            A value.

        Notes
        -----
        Some existing notes.

        """
        note = "This is new.\n"
        result = _insert_into_notes_section(note, doc)
        # The new note should appear inside the Notes section
        self.assertIn("This is new.", result)
        # The existing notes should still be there
        self.assertIn("Some existing notes.", result)
        # The new note should come before the existing notes
        new_pos = result.index("This is new.")
        old_pos = result.index("Some existing notes.")
        self.assertLess(new_pos, old_pos)

    def test_no_notes_section(self):
        doc = """My function.

        Parameters
        ----------
        x : int
            A value.

        """
        note = "This is new.\n"
        result = _insert_into_notes_section(note, doc)
        # A Notes section should be created
        self.assertIn("Notes", result)
        self.assertIn("-----", result)
        self.assertIn("This is new.", result)

    def test_single_line_doc_no_notes(self):
        doc = "Short docstring."
        note = "Added note."
        result = _insert_into_notes_section(note, doc)
        self.assertIn("Notes", result)
        self.assertIn("-----", result)
        self.assertIn("Added note.", result)

    def test_note_trailing_whitespace_stripped(self):
        doc = """My function.

        Notes
        -----
        Existing.

        """
        note = "New note.   \n\n"
        result = _insert_into_notes_section(note, doc)
        self.assertIn("New note.", result)
        # Trailing whitespace from note should not appear
        self.assertNotIn("New note.   ", result)

    def test_multiline_note(self):
        doc = """My function.

        Notes
        -----
        Existing.

        """
        note = "Line one.\nLine two.\n"
        result = _insert_into_notes_section(note, doc)
        self.assertIn("Line one.", result)
        self.assertIn("Line two.", result)

    def test_notes_with_blank_lines_after_header(self):
        """Blank lines between Notes header and content should be handled."""
        doc = """My function.

        Notes
        -----

        Existing content here.

        """
        note = "Inserted."
        result = _insert_into_notes_section(note, doc)
        inserted_pos = result.index("Inserted.")
        existing_pos = result.index("Existing content here.")
        self.assertLess(inserted_pos, existing_pos)


class TestDeprecationNote(unittest.TestCase):
    def test_version_only(self):
        result = _deprecation_note(ver="0.5.0")
        self.assertEqual(result, ".. deprecated:: 0.5.0")

    def test_version_and_message(self):
        result = _deprecation_note(ver="0.5.0", msg="Use bar instead.")
        self.assertEqual(result, ".. deprecated:: 0.5.0 Use bar instead.")

    def test_no_version(self):
        result = _deprecation_note()
        self.assertEqual(result, ".. deprecated:: None")

    def test_no_message(self):
        result = _deprecation_note(ver="1.0", msg=None)
        self.assertEqual(result, ".. deprecated:: 1.0")

    def test_empty_message(self):
        result = _deprecation_note(ver="1.0", msg="")
        self.assertEqual(result, ".. deprecated:: 1.0")


class TestRenamingNote(unittest.TestCase):
    def test_no_version(self):
        result = _renaming_note("old_func")
        self.assertEqual(result, "Alias: ``old_func``")

    def test_with_version(self):
        result = _renaming_note("old_func", ver="0.6.0")
        self.assertIn("versionchanged:: 0.6.0", result)
        self.assertIn("Renamed from ``old_func``", result)
        self.assertIn("kept as an alias", result)
        self.assertNotIn("deprecated", result)
        self.assertTrue(result.endswith("."))

    def test_with_version_and_warn(self):
        result = _renaming_note("old_func", ver="0.6.0", warn=True)
        self.assertIn("versionchanged:: 0.6.0", result)
        self.assertIn("Renamed from ``old_func``", result)
        self.assertIn("but is deprecated", result)
        self.assertTrue(result.endswith("."))

    def test_warn_without_version(self):
        """warn is irrelevant without ver — should return alias form."""
        result = _renaming_note("old_func", warn=True)
        self.assertEqual(result, "Alias: ``old_func``")

    def test_empty_version_string(self):
        result = _renaming_note("old_func", ver="")
        self.assertEqual(result, "Alias: ``old_func``")


if __name__ == '__main__':
    unittest.main()