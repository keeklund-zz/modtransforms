import os
import pickle
import unittest
import pysam
from modtransforms.formats.bamsam import update_header

class TestUpdateHeader(unittest.TestCase):
    """Test class for modtransforms.formats.bamsam.update_header.

    """
    @classmethod
    def setUpClass(self):
        self.file_dir = os.path.dirname(__file__)
        self.to_sizes_file = os.path.join(self.file_dir,
                                     "../test_data/chrom.sizes.mm9.txt")
        self.from_align = os.path.join(self.file_dir,
                                       "../test_data/cast_ATAC_chr19.bam")
        self.cast = pysam.AlignmentFile(self.from_align, 'rb')
        self.header = update_header(self.cast.header, self.to_sizes_file)
        self.known_value = pickle.load(
            open(os.path.join(self.file_dir, "../test_data/mm9_header.dat"), 'r'))
                                                         

    @classmethod
    def tearDownClass(self):
        pass

    def shortDescription(self):
        return None

    def test_return_type(self):
        """Assert return type is dict.

        """
        self.assertIsInstance(self.header, dict)
    
    def test_assert_raises(self):
        """Assert AssertionError on input types.

        """
        self.assertRaises(AssertionError,
                          update_header,
                          [],
                          self.to_sizes_file)
        self.assertRaises(AssertionError,
                          update_header,
                          self.cast.header,
                          "/not/real/file")

    def test_data_validity(self):
        """Assert data equals known values.

        """
        self.assertDictEqual(self.header, self.known_value)
