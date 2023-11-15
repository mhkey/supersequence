import unittest
from supersequencer import super_sequence

class TestSupersequence(unittest.TestCase):
    def test_basic(self):
        seq = [
            [20006, 20471, 211158, 20462, 20461, 21371, 21122, 21141, 21271, 21151],
            [2000442, 21271],
            [20006, 20471, 211158, 20462, 20461, 21371, 21122, 21141, 21271],
            [20006, 2000442, 20009, 20471, 211158, 20462, 20461, 21371, 21122, 21141, 21271],
            [20005, 2000442, 20009, 20471, 211158, 20462, 20461, 21371, 21122, 21141, 21271]

        ]
        s = super_sequence(seq)
        self.assertEqual(s, [20006, 20005, 2000442, 20009, 20471, 211158, 20462, 20461, 21371, 21122, 21141, 21271, 21151])