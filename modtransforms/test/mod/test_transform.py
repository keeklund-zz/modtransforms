import os 
import pickle
import random
import unittest
from modtransforms.mod import transform 
from modtransforms.utils import logger

class SetUpClass(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.mod_file = os.path.join(os.path.dirname(__file__),
                                     "../test_data/CAST_to_BL6_chr19.mod")
        self.logger = logger.build_logger()
        self.reverse = False
        self.transform = transform.build_transform(self.mod_file,
                                                   self.logger,
                                                   self.reverse)
        self.positions, self.deltas = zip(*self.transform.get('chr19'))
        self.position = 0
        self.mod_file = os.path.join(os.path.dirname(__file__),
                                     "../test_data/CAST_to_BL6_chr19.mod")
        self.logger = logger.build_logger()
        self.reverseF = False
        self.reverseT = True
        self.for_trans = transform.build_transform(self.mod_file,
                                                   self.logger,
                                                   self.reverseF)
        self.rev_trans = transform.build_transform(self.mod_file,
                                                   self.logger,
                                                   self.reverseT)
        self.for_file = os.path.join(os.path.dirname(__file__),
                                     "../test_data/CAST_to_BL6_chr19.dat")
        self.rev_file = os.path.join(os.path.dirname(__file__),
                                     "../test_data/CAST_to_BL6_chr19_rev.dat")
        self.for_data = pickle.load(open(self.for_file, 'r'))
        self.rev_data = pickle.load(open(self.rev_file, 'r'))
        
    @classmethod
    def tearDownClass(self):
        self.rev_trans = None
        self.for_trans = None
        
    def shortDescription(self):
        return None
    

class TestDoNothing(SetUpClass):
    """Test class for modtransforms.mod.transforms.do_nothing function.

    """
    def test_return_type(self):
        """Assert return type is int.

        """
        self.assertIsInstance(transform.do_nothing("", 0), int)

    def test_assert_raises(self):
        """Assert AssertionError on input types.

        """
        self.assertRaises(AssertionError, transform.do_nothing, 0, 0)
        self.assertRaises(AssertionError, transform.do_nothing, "", "")

    def test_output_equals(self):
        """Assert output equals specific values.

        """
        self.assertEqual(transform.do_nothing("", 0), 0)
        self.assertEqual(transform.do_nothing("ATCG", 1), 1)


class TestAddition(SetUpClass):
    """Test class for modtransforms.mod.transforms.addition function.

    """
    def test_return_type(self):
        """Assert return type is int.

        """
        delta = transform.addition("", 0)
        self.assertIsInstance(delta, int)

    def test_assert_raises(self):
        """Assert AssertionError on input types.

        """
        self.assertRaises(AssertionError, transform.addition, 0, 0)
        self.assertRaises(AssertionError, transform.addition, "", "")

    def test_output_equals(self):
        """Assert output equals specific values.

        """
        self.assertEqual(transform.addition("", 0), 0)
        self.assertEqual(transform.addition("ATCG", 1), 5)


class TestSubtraction(SetUpClass):
    """Test class for modtransforms.mod.transforms.subtraction function.

    """
    def test_return_type(self):
        """Assert return type is int.

        """
        delta = transform.subtraction("", 0)
        self.assertIsInstance(delta, int)

    def test_assert_raises(self):
        """Assert AssertionError on input types.

        """
        self.assertRaises(AssertionError, transform.subtraction, 0, 0)
        self.assertRaises(AssertionError, transform.subtraction, "", "")

    def test_output_equals(self):
        """Assert output equals specific values.

        """
        self.assertEqual(transform.subtraction("", 0), 0)
        self.assertEqual(transform.subtraction("ATCG", 1), -3)
        

class TestErrorHandler(SetUpClass):
    """Test class for modtransforms.mod.transforms.error_handler function.

    """
    def test_return_type(self):
        """Assert return type is int.

        """
        self.assertRaises(SystemExit, transform.error_handler, "", 0)

    def test_assert_raises(self):
        """Assert AssertionError on input types.

        """
        self.assertRaises(AssertionError, transform.error_handler, 0, 0)
        self.assertRaises(AssertionError, transform.error_handler, "", "")
        

class TestBuildTransform(SetUpClass):
    """Test class for modtransforms.mod.transform.build_transform function.

    """
    def test_assert_raises(self):
        """Assert AssertionError on incorrect input types.

        """
        self.assertRaises(AssertionError,
                          transform.build_transform,
                          "/tmp/clearlyfakefile.mod",
                          self.logger,
                          self.reverseF)
        self.assertRaises(AssertionError,
                          transform.build_transform,
                          self.mod_file,
                          "",
                          self.reverseF)
        self.assertRaises(AssertionError,
                          transform.build_transform,
                          self.mod_file,
                          self.logger,
                          "")

    def test_build_transform_return_types(self):
        """Assert return type dict with key type str and value type list of ints.

        """
        self.assertIsInstance(self.for_trans, dict)
        self.assertIsInstance(self.rev_trans, dict)
        for key in self.for_trans.keys():
            self.assertIsInstance(key, str)
        for key in self.rev_trans.keys():
            self.assertIsInstance(key, str)
        for value in self.for_trans.values():
            self.assertIsInstance(value, list)
        for value in self.rev_trans.values():
            self.assertIsInstance(value, list)
        for value in self.for_trans.values():
            for tup in value:
                self.assertIsInstance(tup, tuple)
                self.assertIsInstance(tup[0], int)
                self.assertIsInstance(tup[1], int)
        for value in self.rev_trans.values():
            for tup in value:
                self.assertIsInstance(tup, tuple)
                self.assertIsInstance(tup[0], int)
                self.assertIsInstance(tup[1], int)
                
    def test_build_transform_validity(self):
        """Assert data equals known values.

        """
        self.assertEqual(self.for_trans.keys(), ['chr19',])
        self.assertEqual(self.rev_trans.keys(), ['chr19',])

        # spot check 100 values - whole list takes too long
        for rand in random.sample(range(len(self.for_trans.get('chr19'))), 100):
            self.assertEqual(self.for_trans.get('chr19')[rand],
                             self.for_data.get('chr19')[rand])
            self.assertEqual(self.rev_trans.get('chr19')[rand],
                             self.rev_data.get('chr19')[rand])


class TestFindDelta(SetUpClass):
    """Test class for modtransforms.mod.transform.find_delta.

    """
    def test_assert_raises(self):
        """Assert AssertionError on incorrect input types.

        """
        self.assertRaises(AssertionError,
                          transform.find_delta,
                          "",
                          self.deltas,
                          self.position)
        self.assertRaises(AssertionError,
                          transform.find_delta,
                          self.positions,
                          "",
                          self.position)
        self.assertRaises(AssertionError,
                          transform.find_delta,
                          self.positions,
                          self.deltas,
                          "")
        
    def test_return_type(self):
        """Assert return type is int.

        """
        position = random.randint(0, len(self.transform.get('chr19')))
        self.assertIsInstance(
            transform.find_delta(self.positions, self.deltas, position),
            int)
                              
    def test_data_validity(self):
        """Assert data equals known value.

        """
        self.assertEqual(
            transform.find_delta(self.positions, self.deltas, 32624884),
            -15328)
        self.assertEqual(
            transform.find_delta(self.positions, self.deltas, 32624885),
            -15328)

