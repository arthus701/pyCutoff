import unittest
from pathlib import Path

from pyCutoff.singletj import singletj

path = str(Path(__file__).parent) + '/'


class TestLines(unittest.TestCase):
    def test_lines(self):
        with open(path + 'TAPE8.TXT', 'r') as fh:
            it = 0
            for line in fh.readlines():
                values = line.split()
                gdlat = float(values[0])
                gdlon = float(values[2])
                rigidity = float(values[5])

                res = singletj(gdlat, gdlon, rigidity)

                self.assertTrue(round(res[1], 2) == float(values[1]))
                self.assertTrue(1 - res[7] == int(values[9]))
                if res[7] == 1:
                    self.assertTrue(round(res[4], 2) == float(values[6]))
                    self.assertTrue(round(res[5], 2) == float(values[7]))

                it += 1
                if 4 <= it:
                    break


if __name__ == '__main__':
    unittest.main()
