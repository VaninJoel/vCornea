
from cc3d import CompuCellSetup
        


from vCornea_TENTATIVE0Steppables import ConstraintInitializerSteppable

CompuCellSetup.register_steppable(steppable=ConstraintInitializerSteppable(frequency=1))




from vCornea_TENTATIVE0Steppables import GrowthSteppable

CompuCellSetup.register_steppable(steppable=GrowthSteppable(frequency=1))




from vCornea_TENTATIVE0Steppables import MitosisSteppable

CompuCellSetup.register_steppable(steppable=MitosisSteppable(frequency=1))




from vCornea_TENTATIVE0Steppables import DeathSteppable

CompuCellSetup.register_steppable(steppable=DeathSteppable(frequency=1))



        
from vCornea_TENTATIVE0Steppables import DifferentiationSteppable

CompuCellSetup.register_steppable(steppable=DifferentiationSteppable(frequency=1))




from vCornea_TENTATIVE0Steppables import PlotSteppable

CompuCellSetup.register_steppable(steppable=PlotSteppable(frequency=1))




from vCornea_TENTATIVE0Steppables import SecretionSteppable

CompuCellSetup.register_steppable(steppable=SecretionSteppable(frequency=1))



CompuCellSetup.run()
