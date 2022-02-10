import unittest
import numpy as np
from evaluate import error_profile

class ErrorProfileTest(unittest.TestCase):

    def test_matches(self):

        ref = ['A','C','T','A','G','C','T','A']
        ali = ['|','|','|','|','|','|','|','|']
        que = ['A','C','T','A','G','C','T','A']
        arr = np.array([ref, ali, que])
        
        signatures = error_profile(arr)
        self.assertEqual(signatures['AC>CT'], 1)
        self.assertEqual(signatures['CT>TA'], 2)
        self.assertEqual(signatures['TA>AG'], 1)
        self.assertEqual(signatures['AG>GC'], 1)
        self.assertEqual(signatures['GC>CT'], 1)
        t = 0
        for v in signatures.values():
            if v == 0:
                t += 1
        self.assertEqual(t, len(signatures) - 5)

    def test_mismatches(self):

        ref = ['A','C','T','A','G','C','T','A']
        ali = ['|','.','.','|','|','.','|','|']
        que = ['A','G','C','A','G','T','T','A']
        arr = np.array([ref, ali, que])
        
        signatures = error_profile(arr)
        self.assertEqual(signatures['AG>CC'], 1)
        self.assertEqual(signatures['GC>TA'], 1)
        self.assertEqual(signatures['CA>AG'], 1)
        self.assertEqual(signatures['AG>GT'], 1)
        self.assertEqual(signatures['GT>CT'], 1)
        self.assertEqual(signatures['TT>TA'], 1)
        t = 0
        for v in signatures.values():
            if v == 0:
                t += 1
        self.assertEqual(t, len(signatures) - 6)

    def test_deletions(self):

        ref = ['A','C','T','A','G','C','T','A']
        ali = ['|',' ','|','|',' ',' ','|','|']
        que = ['A','-','T','A','-','-','T','A']
        arr = np.array([ref, ali, que])
        
        signatures = error_profile(arr)
        self.assertEqual(signatures['A->CT'], 2)
        self.assertEqual(signatures['AT>TA'], 2)
        self.assertEqual(signatures['TA>AT'], 1)
        self.assertEqual(signatures['A->GT'], 1)
        t = 0
        for v in signatures.values():
            if v == 0:
                t += 1
        self.assertEqual(t, len(signatures) - 4)

    def test_insertions(self):

        ref = ['A','-','T','A','-','-','T','A']
        ali = ['|',' ','|','|',' ',' ','|','|']
        que = ['A','C','T','A','G','C','T','A']
        arr = np.array([ref, ali, que])
        
        signatures = error_profile(arr)
        self.assertEqual(signatures['AC>-T'], 1)
        self.assertEqual(signatures['CT>TA'], 2)
        self.assertEqual(signatures['TA>AG'], 1)
        self.assertEqual(signatures['AG>-C'], 1)
        self.assertEqual(signatures['GC>-T'], 1)

        t = 0
        for v in signatures.values():
            if v == 0:
                t += 1
        self.assertEqual(t, len(signatures) - 5)

    def test_mix(self):

        ref = ['A','-','T','A','-','G','T','A']
        ali = ['|',' ',' ','|',' ','.',' ','|']
        que = ['A','C','-','A','G','C','-','A']
        arr = np.array([ref, ali, que])
        
        signatures = error_profile(arr)
        self.assertEqual(signatures['AC>-A'], 1)
        self.assertEqual(signatures['C->TA'], 2)
        self.assertEqual(signatures['CA>AG'], 1)
        self.assertEqual(signatures['AG>-C'], 1)
        self.assertEqual(signatures['GC>GA'], 1)
        
        t = 0
        for v in signatures.values():
            if v == 0:
                t += 1
        self.assertEqual(t, len(signatures) - 5)

    def test_wrong_input(self):

        ref = ['-','-','T','A','-','G','T','A']
        ali = [' ',' ',' ','|',' ','.',' ','|']
        que = ['C','C','-','A','G','C','-','A']
        arr = np.array([ref, ali, que])

        with self.assertRaises(ValueError):
            error_profile(arr)
        
        ref = ['C','-','T','A','-','G','T','-']
        ali = ['|',' ',' ','|',' ','.',' ',' ']
        que = ['C','C','-','A','G','C','-','A']
        arr = np.array([ref, ali, que])

        with self.assertRaises(ValueError):
            error_profile(arr)
        
        ref = ['-','-','T','A','-','G','T','-']
        ali = [' ',' ',' ','|',' ','.',' ',' ']
        que = ['C','C','-','A','G','C','-','A']
        arr = np.array([ref, ali, que])

        with self.assertRaises(ValueError):
            error_profile(arr)