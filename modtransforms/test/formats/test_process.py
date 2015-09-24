import os
import random
import unittest
import pysam
from modtransforms.formats.process import build_modification_index
from modtransforms.formats.process import get_positions_and_deltas
from modtransforms.mod.transform import build_transform, find_delta
from modtransforms.utils.logger import build_logger

class SetUpClass(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.alignment_file = os.path.join(os.path.dirname(__file__),
                                           "../test_data/cast_ATAC_chr19.bam")
        self.input_ = pysam.AlignmentFile(self.alignment_file, 'rb')
        self.lines = random.sample([i for i in self.input_], 100)
        self.mod_file = os.path.join(
            os.path.dirname(__file__), "../test_data/CAST_to_BL6_chr19.mod")
        self.logger = build_logger()
        self.curr_chrom = 'chr19'
        self.chrom_mods = build_transform(self.mod_file, self.logger)
        self.positions, self.deltas = get_positions_and_deltas(self.chrom_mods,
                                                               self.curr_chrom,
                                                               self.logger)
        self.output = (self.positions, self.deltas)

    @classmethod
    def tearDownClass(self):
        pass

    def shortDescription(self):
        return None


# this will probably move
class TestUpdateCigar(SetUpClass):
    """Test class for modtransforms.formats.process.update_cigar.

    """
    def test_honk(self):
        pass

    
# this will probably move
# add assertRaises
class TestBuildModificationIndex(SetUpClass):
    """Test class for modtransforms.formats.process.build_modificantion_index.

    """
    def test_return_types(self):
        """Assert return type is list of tuples.

        """
        for line in self.lines:
            start_delta = find_delta(self.positions,
                                     self.deltas,
                                     int(line.reference_start))
            mod_index = build_modification_index(self.positions,
                                                 self.deltas,
                                                 line,
                                                 start_delta)
            self.assertIsInstance(mod_index, list)
            for mi in mod_index:
                self.assertIsInstance(mi, tuple)
                self.assertIsInstance(mi[0], int)
                self.assertIsInstance(mi[1], int)
                    
            
# this will probably move 
class TestGetPositionsAndDeltas(SetUpClass):
    """Test class for modtransforms.formats.process.get_positions_and_deltas.

    """
    def test_return_types(self):
        """Assert return type is tupe type: (tuple, tuple).

        """
        self.assertIsInstance(self.output, tuple)
        self.assertIsInstance(self.output[0], tuple)
        self.assertIsInstance(self.output[1], tuple)
        self.assertIsInstance(self.output[0][0], int)
        self.assertIsInstance(self.output[0][1], int)

    def test_assert_raises(self):
        self.assertRaises(AssertionError,
                          get_positions_and_deltas,
                          [],
                          self.curr_chrom,
                          self.logger)
        self.assertRaises(AssertionError,
                          get_positions_and_deltas,
                          self.chrom_mods,
                          [],
                          self.logger)
        self.assertRaises(AssertionError,
                          get_positions_and_deltas,
                          self.chrom_mods,
                          self.curr_chrom,
                          [])

    
class TestBam(SetUpClass):
    """Test class for modtransforms.formats.process.bam.

    """
    def test_return_type(self):
        pass

    def test_assert_raises(self):
        pass

    
