
from pywfn import reader

class MoldelReader(reader.Reader):
    def __init__(self, path:str):
        super().__init__(path)