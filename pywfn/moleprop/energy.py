from pywfn.base import Mole


class Calculator:
    def __init__(self,mole:Mole) -> None:
        self.mole=mole

    def VIP(self,molP:Mole): #First vertical ionization potential
        return molP.energy-self.mole.energy