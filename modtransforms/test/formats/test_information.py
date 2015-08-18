import types
import unittest
from modtransforms.formats import information

class TestFileTypes(unittest.TestCase):
    """Test case for modtransforms.formats.information.

    """
    def shortDescription(self):
        return None

    def test_types(self):
        """Assert file_types is type dict with string keys and function values.

        """
        self.assertIsInstance(information.file_types, dict)
        self.assertIsInstance(information.file_types.keys(), list)
        self.assertIsInstance(information.file_types.values(), list)
        for key in information.file_types.keys():
            self.assertIsInstance(key, str)
        for value in information.file_types.values():
            self.assertIsInstance(value, types.FunctionType)
    
