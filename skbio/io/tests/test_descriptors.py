# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import unittest
from io import StringIO
from skbio.io.descriptors import Read, Write, _docstring_vars
from skbio.io.registry import create_format, io_registry
from skbio.io._exception import UnrecognizedFormatError


class TestRead(unittest.TestCase):
    """Tests for the Read descriptor."""
    
    def setUp(self):
        """Set up test fixtures."""
        # Create a mock format for testing
        self.test_format = create_format('test_format')
        
        # Create a test class with the Read descriptor
        class TestClass:
            read = Read()
            
            def __init__(self, data):
                self.data = data
        
        self.TestClass = TestClass
        
        # Register a reader for the test format
        @self.test_format.reader(TestClass)
        def _test_reader(fh, cls=None):
            if cls is None:
                cls = TestClass
            content = fh.read()
            return cls(content)
        
    def tearDown(self):
        """Clean up test format."""
        io_registry.remove_format('test_format')
    
    def test_read_descriptor_creates_method(self):
        """Test that accessing read creates a _read_method."""
        self.assertNotIn('_read_method', self.TestClass.__dict__)
        method = self.TestClass.read
        self.assertIn('_read_method', self.TestClass.__dict__)
        self.assertIsNotNone(method)
    
    def test_read_method_caching(self):
        """Test that _read_method is cached on the class."""
        method1 = self.TestClass.read
        method2 = self.TestClass.read
        self.assertIs(method1, method2)
    
    def test_read_returns_correct_class(self):
        """Test that read returns an instance of the calling class."""
        fh = StringIO("test data")
        obj = self.TestClass.read(fh, format='test_format')
        self.assertIsInstance(obj, self.TestClass)
        self.assertEqual(obj.data, "test data")
    
    def test_read_with_subclass(self):
        """Test that subclasses get their own read method."""
        class SubClass(self.TestClass):
            pass
        
        fh = StringIO("subclass data")
        obj = SubClass.read(fh, format='test_format')
        self.assertIsInstance(obj, SubClass)
        self.assertEqual(obj.data, "subclass data")
    
    def test_read_subclass_returns_subclass_instance(self):
        """Test that reading with a subclass returns subclass instance."""
        class SubClass(self.TestClass):
            def custom_method(self):
                return "custom"
        
        fh = StringIO("data")
        obj = SubClass.read(fh, format='test_format')
        self.assertIsInstance(obj, SubClass)
        self.assertTrue(hasattr(obj, 'custom_method'))
        self.assertEqual(obj.custom_method(), "custom")
    
    def test_read_method_has_docstring(self):
        """Test that the generated read method has a docstring."""
        method = self.TestClass.read
        self.assertIsNotNone(method.__doc__)
        self.assertIn('TestClass', method.__doc__)

    def test_tabular_msa_read_docstring_mentions_constructor(self):
        """Test that TabularMSA read docstring mentions constructor."""
        test_format = create_format('tabular_msa_docstring_test_format')

        class TabularMSA:
            read = Read()

        @test_format.reader(TabularMSA)
        def _reader(fh, cls=None):
            if cls is None:
                cls = TabularMSA
            return cls()

        try:
            doc = TabularMSA.read.__doc__
            self.assertIn('constructor', doc)
            self.assertIn('GrammaredSequence', doc)
        finally:
            io_registry.remove_format('tabular_msa_docstring_test_format')
    
    def test_read_passes_kwargs(self):
        """Test that kwargs are passed through to the reader."""
        @self.test_format.reader(self.TestClass, override=True)
        def _test_reader_with_kwargs(fh, cls=None, custom_arg=None):
            if cls is None:
                cls = self.TestClass
            return cls(f"{fh.read()}-{custom_arg}")
        
        fh = StringIO("test")
        obj = self.TestClass.read(fh, format='test_format', custom_arg='value')
        self.assertEqual(obj.data, "test-value")
    
    def test_read_with_format_inference(self):
        """Test reading without explicitly providing format."""
        @self.test_format.sniffer()
        def _test_sniffer(fh):
            content = fh.read(4)
            return content == "test", {}
        
        fh = StringIO("test data")
        obj = self.TestClass.read(fh)
        self.assertIsInstance(obj, self.TestClass)


class TestWrite(unittest.TestCase):
    """Tests for the Write descriptor."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.test_format = create_format('test_write_format')
        
        class TestClass:
            default_write_format = 'test_write_format'
            write = Write()
            
            def __init__(self, data):
                self.data = data
        
        self.TestClass = TestClass
        
        @self.test_format.writer(TestClass)
        def _test_writer(obj, fh):
            fh.write(obj.data)
        
    def tearDown(self):
        """Clean up test format."""
        io_registry.remove_format('test_write_format')
    
    def test_write_descriptor_on_class(self):
        """Test accessing write on the class itself."""
        method = self.TestClass.write
        self.assertIn('_write_method', self.TestClass.__dict__)
    
    def test_write_descriptor_on_instance(self):
        """Test accessing write on an instance."""
        obj = self.TestClass("test data")
        method = obj.write
        self.assertIn('_write_method', obj.__dict__)
    
    def test_write_method_caching_on_class(self):
        """Test that write method is cached on the class."""
        method1 = self.TestClass.write
        method2 = self.TestClass.write
        self.assertIs(method1, method2)
    
    def test_write_method_caching_on_instance(self):
        """Test that write method is cached on instances."""
        obj = self.TestClass("data")
        method1 = obj.write
        method2 = obj.write
        self.assertIs(method1, method2)
    
    def test_write_works(self):
        """Test that write actually writes data."""
        obj = self.TestClass("test data")
        fh = StringIO()
        obj.write(fh, format='test_write_format')
        self.assertEqual(fh.getvalue(), "test data")
    
    def test_write_with_subclass(self):
        """Test that subclasses can write."""
        class SubClass(self.TestClass):
            pass
        
        obj = SubClass("subclass data")
        fh = StringIO()
        obj.write(fh, format='test_write_format')
        self.assertEqual(fh.getvalue(), "subclass data")
    
    def test_write_uses_default_format(self):
        """Test that write uses default_write_format when format not provided."""
        obj = self.TestClass("data")
        fh = StringIO()
        result = obj.write(fh)  # No format specified
        self.assertEqual(fh.getvalue(), "data")
    
    def test_write_raises_without_default_format(self):
        """Test that write raises error when no default format and none provided."""
        class NoDefaultClass:
            write = Write()
            
            def __init__(self, data):
                self.data = data
        
        obj = NoDefaultClass("data")
        fh = StringIO()
        with self.assertRaises(ValueError) as cm:
            obj.write(fh)
        self.assertIn("no default write format", str(cm.exception))
    
    def test_write_method_has_docstring(self):
        """Test that generated write method has a docstring."""
        obj = self.TestClass("data")
        method = obj.write
        self.assertIsNotNone(method.__doc__)
        self.assertIn('TestClass', method.__doc__)
    
    def test_write_returns_file_path(self):
        """Test that write returns the file path/handle passed in."""
        obj = self.TestClass("data")
        fh = StringIO()
        result = obj.write(fh, format='test_write_format')
        self.assertIs(result, fh)
    
    def test_write_passes_kwargs(self):
        """Test that kwargs are passed through to the writer."""
        @self.test_format.writer(self.TestClass, override=True)
        def _test_writer_with_kwargs(obj, fh, prefix=""):
            fh.write(f"{prefix}{obj.data}")
        
        obj = self.TestClass("data")
        fh = StringIO()
        obj.write(fh, format='test_write_format', prefix="PREFIX:")
        self.assertEqual(fh.getvalue(), "PREFIX:data")


class TestDocstringVars(unittest.TestCase):
    """Tests for the _docstring_vars helper function."""
    
    def setUp(self):
        """Set up test format and class."""
        self.test_format = create_format('docstring_test_format')
        
        class TestClass:
            default_write_format = 'docstring_test_format'
        
        self.TestClass = TestClass
        
        @self.test_format.reader(TestClass)
        def _reader(fh, cls=None):
            pass
        
        @self.test_format.writer(TestClass)
        def _writer(obj, fh):
            pass
    
    def tearDown(self):
        """Clean up."""
        io_registry.remove_format('docstring_test_format')
    
    def test_docstring_vars_for_read(self):
        """Test _docstring_vars returns correct values for read."""
        name, supported_fmts, default, see = _docstring_vars(self.TestClass, 'read')
        self.assertEqual(name, 'TestClass')
        self.assertIn('docstring_test_format', supported_fmts)
        self.assertIn('skbio.io.format.docstring_test_format', see)
    
    def test_docstring_vars_for_write(self):
        """Test _docstring_vars returns correct values for write."""
        name, supported_fmts, default, see = _docstring_vars(self.TestClass, 'write')
        self.assertEqual(name, 'TestClass')
        self.assertIn('docstring_test_format', supported_fmts)
        self.assertEqual(default, 'docstring_test_format')
    
    def test_docstring_vars_invalid_func(self):
        """Test _docstring_vars raises error for invalid func."""
        with self.assertRaises(ValueError) as cm:
            _docstring_vars(self.TestClass, 'invalid')
        self.assertIn("'read' or 'write'", str(cm.exception))
    
    def test_docstring_vars_no_formats(self):
        """Test _docstring_vars when class has no registered formats."""
        class NoFormatClass:
            pass
        
        name, supported_fmts, default, see = _docstring_vars(NoFormatClass, 'read')
        self.assertEqual(name, 'NoFormatClass')
        self.assertEqual(supported_fmts, "")
        self.assertEqual(default, "None")


class TestInheritanceIntegration(unittest.TestCase):
    """Integration tests for descriptor behavior with inheritance."""
    
    def setUp(self):
        """Set up test format and class hierarchy."""
        self.test_format = create_format('inheritance_test_format')
        
        class ParentClass:
            default_write_format = 'inheritance_test_format'
            read = Read()
            write = Write()
            
            def __init__(self, data):
                self.data = data
        
        class ChildClass(ParentClass):
            pass
        
        class GrandchildClass(ChildClass):
            pass
        
        self.ParentClass = ParentClass
        self.ChildClass = ChildClass
        self.GrandchildClass = GrandchildClass
        
        @self.test_format.reader(ParentClass)
        def _reader(fh, cls=None):
            if cls is None:
                cls = ParentClass
            return cls(fh.read())
        
        @self.test_format.writer(ParentClass)
        def _writer(obj, fh):
            fh.write(f"[{obj.__class__.__name__}:{obj.data}]")

        @self.test_format.sniffer()
        def _sniffer(fh):
            content = fh.read(10)
            return content.startswith("TESTDATA:"), {}
    
    def tearDown(self):
        """Clean up."""
        io_registry.remove_format('inheritance_test_format')
    
    def test_child_inherits_read(self):
        """Test that child class inherits parent's read capability."""
        fh = StringIO("child data")
        obj = self.ChildClass.read(fh, format='inheritance_test_format')
        self.assertIsInstance(obj, self.ChildClass)
        self.assertEqual(obj.data, "child data")
    
    def test_grandchild_inherits_read(self):
        """Test that grandchild class inherits parent's read capability."""
        fh = StringIO("grandchild data")
        obj = self.GrandchildClass.read(fh, format='inheritance_test_format')
        self.assertIsInstance(obj, self.GrandchildClass)
        self.assertEqual(obj.data, "grandchild data")
    
    def test_child_inherits_write(self):
        """Test that child class inherits parent's write capability."""
        obj = self.ChildClass("child data")
        fh = StringIO()
        obj.write(fh, format='inheritance_test_format')
        self.assertEqual(fh.getvalue(), "[ChildClass:child data]")
    
    def test_grandchild_inherits_write(self):
        """Test that grandchild class inherits parent's write capability."""
        obj = self.GrandchildClass("grandchild data")
        fh = StringIO()
        obj.write(fh, format='inheritance_test_format')
        self.assertEqual(fh.getvalue(), "[GrandchildClass:grandchild data]")
    
    def test_each_class_has_own_cached_method(self):
        """Test that each class in hierarchy caches its own methods."""
        parent_read = self.ParentClass.read
        child_read = self.ChildClass.read
        grandchild_read = self.GrandchildClass.read
        
        # They should all be different cached methods
        self.assertIsNot(parent_read, child_read)
        self.assertIsNot(child_read, grandchild_read)
        self.assertIsNot(parent_read, grandchild_read)

    def test_list_read_formats_parent(self):
        """Test that parent class lists its read formats."""
        formats = io_registry.list_read_formats(self.ParentClass)
        self.assertIn('inheritance_test_format', formats)
    
    def test_list_read_formats_child(self):
        """Test that child class inherits parent's read formats."""
        formats = io_registry.list_read_formats(self.ChildClass)
        self.assertIn('inheritance_test_format', formats)
    
    def test_list_read_formats_grandchild(self):
        """Test that grandchild class inherits parent's read formats."""
        formats = io_registry.list_read_formats(self.GrandchildClass)
        self.assertIn('inheritance_test_format', formats)
    
    def test_list_write_formats_parent(self):
        """Test that parent class lists its write formats."""
        formats = io_registry.list_write_formats(self.ParentClass)
        self.assertIn('inheritance_test_format', formats)
    
    def test_list_write_formats_child(self):
        """Test that child class inherits parent's write formats."""
        formats = io_registry.list_write_formats(self.ChildClass)
        self.assertIn('inheritance_test_format', formats)
    
    def test_list_write_formats_grandchild(self):
        """Test that grandchild class inherits parent's write formats."""
        formats = io_registry.list_write_formats(self.GrandchildClass)
        self.assertIn('inheritance_test_format', formats)
    
    def test_sniff_with_parent_into(self):
        """Test that sniffing works when into is parent class."""
        fh = StringIO("TESTDATA:some content")
        format_name, kwargs = io_registry.sniff(fh, into=self.ParentClass)
        self.assertEqual(format_name, 'inheritance_test_format')
    
    def test_sniff_with_child_into(self):
        """Test that sniffing works when into is child class."""
        fh = StringIO("TESTDATA:some content")
        format_name, kwargs = io_registry.sniff(fh, into=self.ChildClass)
        self.assertEqual(format_name, 'inheritance_test_format')
    
    def test_sniff_with_grandchild_into(self):
        """Test that sniffing works when into is grandchild class."""
        fh = StringIO("TESTDATA:some content")
        format_name, kwargs = io_registry.sniff(fh, into=self.GrandchildClass)
        self.assertEqual(format_name, 'inheritance_test_format')
    
    def test_read_with_format_inference_child(self):
        """Test that child class can read with format inference."""
        fh = StringIO("TESTDATA:child data")
        obj = self.ChildClass.read(fh)  # No format specified
        self.assertIsInstance(obj, self.ChildClass)
        self.assertEqual(obj.data, "TESTDATA:child data")
    
    def test_read_with_format_inference_grandchild(self):
        """Test that grandchild class can read with format inference."""
        fh = StringIO("TESTDATA:grandchild data")
        obj = self.GrandchildClass.read(fh)  # No format specified
        self.assertIsInstance(obj, self.GrandchildClass)
        self.assertEqual(obj.data, "TESTDATA:grandchild data")
    
    def test_format_listing_no_duplicates(self):
        """Test that format listing doesn't return duplicates for inherited formats."""
        formats = io_registry.list_read_formats(self.ChildClass)
        # Count occurrences of the format
        count = formats.count('inheritance_test_format')
        self.assertEqual(count, 1, "Format should appear only once in the list")
    
    def test_multiple_inheritance_format_listing(self):
        """Test format listing with multiple inheritance."""
        # Create another format and parent class
        other_format = create_format('other_test_format')
        
        class OtherParent:
            read = Read()
            write = Write()
            
            def __init__(self, data):
                self.data = data
        
        @other_format.reader(OtherParent)
        def _other_reader(fh, cls=None):
            if cls is None:
                cls = OtherParent
            return cls(fh.read())
        
        # Create a class with multiple inheritance
        class MultiChild(self.ParentClass, OtherParent):
            pass
        
        try:
            formats = io_registry.list_read_formats(MultiChild)
            # Should inherit formats from both parents
            self.assertIn('inheritance_test_format', formats)
            self.assertIn('other_test_format', formats)
        finally:
            io_registry.remove_format('other_test_format')
    
    def test_sniff_filters_by_child_class(self):
        """Test that sniff only considers formats appropriate for the child class."""
        # Create a format that's only for a different class
        unrelated_format = create_format('unrelated_format')
        
        class UnrelatedClass:
            read = Read()
            
            def __init__(self, data):
                self.data = data
        
        @unrelated_format.reader(UnrelatedClass)
        def _unrelated_reader(fh, cls=None):
            if cls is None:
                cls = UnrelatedClass
            return cls(fh.read())
        
        @unrelated_format.sniffer()
        def _unrelated_sniffer(fh):
            # This sniffer would match, but shouldn't be considered
            return True, {}
        
        try:
            # Even though unrelated_sniffer always returns True,
            # it shouldn't match because ChildClass can't use it
            fh = StringIO("some data")
            with self.assertRaises(UnrecognizedFormatError):
                io_registry.sniff(fh, into=self.ChildClass)
        finally:
            io_registry.remove_format('unrelated_format')


class TestInheritanceEdgeCases(unittest.TestCase):
    """Test edge cases in inheritance behavior."""
    
    def setUp(self):
        """Set up test formats."""
        self.format1 = create_format('format1_edge')
        self.format2 = create_format('format2_edge')
        
        class GrandParent:
            read = Read()
            write = Write()
            
            def __init__(self, data):
                self.data = data
        
        class Parent(GrandParent):
            pass
        
        class Child(Parent):
            pass
        
        self.GrandParent = GrandParent
        self.Parent = Parent
        self.Child = Child
        
        # Register reader only for GrandParent in format1
        @self.format1.reader(GrandParent)
        def _format1_reader(fh, cls=None):
            if cls is None:
                cls = GrandParent
            return cls("format1:" + fh.read())
        
        # Register reader only for Parent in format2
        @self.format2.reader(Parent)
        def _format2_reader(fh, cls=None):
            if cls is None:
                cls = Parent
            return cls("format2:" + fh.read())
    
    def tearDown(self):
        """Clean up."""
        io_registry.remove_format('format1_edge')
        io_registry.remove_format('format2_edge')
    
    def test_child_inherits_from_grandparent(self):
        """Test that Child can use GrandParent's reader."""
        fh = StringIO("data")
        obj = self.Child.read(fh, format='format1_edge')
        self.assertIsInstance(obj, self.Child)
        self.assertEqual(obj.data, "format1:data")
    
    def test_child_inherits_from_parent(self):
        """Test that Child can use Parent's reader."""
        fh = StringIO("data")
        obj = self.Child.read(fh, format='format2_edge')
        self.assertIsInstance(obj, self.Child)
        self.assertEqual(obj.data, "format2:data")
    
    def test_child_lists_both_formats(self):
        """Test that Child lists formats from both Parent and GrandParent."""
        formats = io_registry.list_read_formats(self.Child)
        self.assertIn('format1_edge', formats)
        self.assertIn('format2_edge', formats)
    
    def test_parent_has_both_formats(self):
        """Test that Parent has both its own and inherited formats."""
        formats = io_registry.list_read_formats(self.Parent)
        self.assertIn('format1_edge', formats)
        self.assertIn('format2_edge', formats)
    
    def test_grandparent_has_only_own_format(self):
        """Test that GrandParent only has its own format."""
        formats = io_registry.list_read_formats(self.GrandParent)
        self.assertIn('format1_edge', formats)
        self.assertNotIn('format2_edge', formats)


if __name__ == '__main__':
    unittest.main()
