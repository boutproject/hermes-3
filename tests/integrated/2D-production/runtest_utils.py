import hashlib
from pathlib import Path
import shutil
import urllib.request
import zipfile

url = "https://zenodo.org/records/22113649/files/test-2D-production-2026-08-26.zip"

expected_hash = "f0e02d0599c56a3a95ad28702e50b94ab12bfe3d3b25f6d6a0800851e177cd0c"
expected_filenames = [
    "BOUT.restart.0.nc",
    "BOUT.restart.1.nc",
    "BOUT.restart.2.nc",
    "BOUT.restart.3.nc",
    "BOUT.restart.4.nc",
    "BOUT.restart.5.nc",
    "BOUT.restart.6.nc",
    "BOUT.restart.7.nc",
    "BOUT.restart.8.nc",
    "BOUT.restart.9.nc",
    "grid_test2_allpump.nc",
]


class DataDownloader:
    def __init__(self, directory=Path(__file__).parent):
        self.directory = directory
        self.data = directory / "data"
        self.zipfile_path = self.directory / "test-2D-production.zip"

    def download_data(self):
        ## Download files
        tmp_path = self.zipfile_path.with_name(self.zipfile_path.name + ".tmp")

        with urllib.request.urlopen(url, timeout=60) as response:
            if response.status != 200:
                raise RuntimeError(
                    f"2D-Production test: download failed - HTTP {response.status}"
                )

            # Copy bits of the file from response to a temp file
            # This ensures no partial files are left if the download fails
            with open(tmp_path, "wb") as out_file:
                shutil.copyfileobj(response, out_file)

        # Rename temp file with the correct name
        tmp_path.replace(self.zipfile_path)

    def extractZip(self):
        with zipfile.ZipFile(self.zipfile_path, "r") as zf:
            zip_contents = set(zf.namelist())
            try:
                # Extract only expected grids
                for filename in expected_filenames:
                    if filename in zip_contents:
                        dest = self.data if filename.startswith("BOUT.restart") else self.directory
                        zf.extract(filename, path=dest)
            except Exception as e:
                print("2D-Production test: extracting test grids failed:", e)
                raise

    def checkHash(self, doRaise=True):
        # Check hash
        try:
            with open(self.zipfile_path, "rb") as f:
                file_hash = hashlib.sha256(f.read()).hexdigest()
        except FileNotFoundError:
            if doRaise:
                raise
            return False

        if doRaise:
            if file_hash != expected_hash:
                raise RuntimeError(
                    "2D-Production test: downloaded zip file hash does not match expected value"
                )
        return file_hash == expected_hash
