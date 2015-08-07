import os
import types
import unittest
from modtransforms.mod import stream

class TestGenMod(unittest.TestCase):
    """Unit tests for modtransforms.mod.stream.gen_mod.

    gen_mod returns a MOD file python generator.

    """
    @classmethod
    def setUpClass(self):
        """Initialize test with specific mod file subset.

        """
        self.mod_file = os.path.join(os.path.dirname(__file__),
                                     "../test_data/CAST_to_BL6_chr19.mod")
        self.generator = stream.gen_mod(self.mod_file)
        self.data = self.generator.next()
        
    @classmethod
    def tearDownClass(self):
        pass

    def shortDescription(self):
        return None
    
    def test_return_type(self):
        """Assert gen_mod creates python generator type.

        """
        self.assertIsInstance(self.generator, types.GeneratorType)

    def test_ignore_comments(self):
        """Assert MOD file comments are not included in data.

        """
        self.assertTrue('#' not in self.data)

    def test_data(self):
        """Assert data is first line of chr19 from CAST MOD file.

        """
        self.assertTrue(self.data, "s       chr19   3072851 G/A")

