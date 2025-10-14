# test_descriptors.py
import unittest
from io import StringIO
from unittest.mock import Mock, patch, MagicMock
from skbio.io.descriptors import Read, Write, _docstring_vars
from skbio.io.registry import create_format, io_registry


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


if __name__ == '__main__':
    unittest.main()