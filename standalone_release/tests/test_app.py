import unittest
from app import attributes, strip_prefix

class ParserTests(unittest.TestCase):
    def test_gff_attributes(self):
        parsed = attributes("ID=gene:Solyc01g000010.4;Alias=ACS7;Note=ethylene enzyme")
        self.assertEqual(parsed["Alias"], "ACS7")
        self.assertEqual(strip_prefix(parsed["ID"]), "Solyc01g000010.4")

if __name__ == "__main__":
    unittest.main()
