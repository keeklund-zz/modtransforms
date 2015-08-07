import logging
import unittest
from modtransforms.utils import logger 

class TestBuildLogger(unittest.TestCase):
    """Test all functions in modtransforms.utils.logger.

    """
    @classmethod
    def setUpClass(self):
        self.logger = logger.build_logger()

    @classmethod
    def tearDownClass(self):
        pass

    def shortDescription(self):
        return None
    
    def test_class(self):
        """Assert correct logger class.

        """
        self.assertEqual(self.logger.getLoggerClass(), logging.Logger)

