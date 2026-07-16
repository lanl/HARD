import tempfile
import unittest
from pathlib import Path

from verify_lib import find_last_output


class FindLastOutputTest(unittest.TestCase):
    def test_combines_csv_files_without_cat(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for rank in (0, 1):
                for step in (1, 2):
                    path = root / f"output-2D-{rank}-0000{step}.csv"
                    path.write_text(f"rank{rank}\n")

            output, combined = find_last_output(dir=directory)

            self.assertTrue(combined)
            self.assertEqual(Path(output).read_text(), "rank0\nrank1\n")


if __name__ == "__main__":
    unittest.main()
