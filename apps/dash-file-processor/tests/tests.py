import unittest
from aktools.conf import *
from aktools.tools import *


class NamesTestCase(unittest.TestCase):

    # Test connected_hostnames() on short and long files
    def test_connected_hostnames_sf(self):
        result = connected_hostnames(
            "data/input_test_case_1.txt", 1607880434801, 1607880438820, "Steeve"
        )
        self.assertEqual(result, {"Hanny": 1, "Hannibal": 2})

    def test_connected_hostnames_lf(self):
        result = connected_hostnames(
            "data/input-file.txt", 1565647204351, 1565733598341, "Dristen"
        )
        self.assertEqual(
            result,
            {
                "Aadison": 1,
                "Wilkens": 1,
                "Kahlina": 1,
                "Alei": 1,
                "Zhanasia": 1,
                "Jamor": 1,
                "Joy": 1,
            },
        )

    # Test connected_to() on short and long files
    def test_connect_to_sf(self):
        result = connected_to(
            "data/input_test_case_1.txt", 1607880434801, 1607880438820, "Steeve"
        )
        self.assertEqual(result, {"Hannibal": 1})

    def test_connect_to_lf(self):
        result = connected_to(
            "data/input-file.txt", 1565647204351, 1565733598341, "Jadon"
        )
        self.assertEqual(
            result,
            {
                "Ahmya": 1,
                "Kayleann": 1,
                "Shainah": 1,
                "Aniyah": 1,
                "Eveleigh": 1,
                "Caris": 1,
                "Rahniya": 1,
                "Remiel": 1,
            },
        )

    # Test received_from() on short and long files
    def test_received_from_sf(self):
        result = received_from(
            "data/input_test_case_1.txt", 1607880434801, 1607880438820, "Steeve"
        )
        self.assertEqual(result, {"Hannibal": 1, "Hanny": 1})

    def test_received_from_lf(self):
        result = received_from(
            "data/input-file.txt", 1565647204351, 1565733598341, "Dristen"
        )
        self.assertEqual(
            result,
            {
                "Joy": 1,
                "Jamor": 1,
                "Zhanasia": 1,
                "Alei": 1,
                "Kahlina": 1,
                "Wilkens": 1,
                "Aadison": 1,
            },
        )

    # Test generated_conn
    def test_generated_conn(self):
        result = generated_conn(
            "data/input_test_case_1.txt", 1607880434801, 1607880438820
        )
        self.assertEqual(result, {"Hannibal": 3, "Steeve": 2, "Hanny": 1})
