import io
import unittest
from unittest.mock import Mock, patch, PropertyMock
from src.config import read_config

class TestConfig(unittest.TestCase):

    def test_read_config(self):
        config = io.StringIO(
            "foo='bar'\n"
            "bar=42.2\n"
            "#Hello this is a comment\n"
            "\n"
            "baz=True\n"
        )
        expected = {
            "foo": "bar",
            "bar": 42.2,
            "baz": True,
        }
        self.assertEquals(read_config(config), expected)


###################
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestConfig))
    return suite

if __name__ == '__main__':
    unittest.main()
