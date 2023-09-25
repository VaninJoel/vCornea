#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# Simulation Version 2 Paper 
# 
# ----- ADDITIONS -----
#       * Making the simulation more physically accurate
#           * With a constant diffusion and decay for EGF
#           * Correction for boundary conditions for EGF and movement bias       
#           * Adding a wall to the simulation to avoid cells to crawl on the edge of the simulation
#           * Changing cell volume and superficial constraints to be more realistic
#           * Update of differentiation rule based on the area of contact calculation and not only contact
#           * Fixing the calculation of growth rate to the correct form 
#           * Adding a calculated cell death rate to be more realistic to the corneal epithelium turnover time  
#           * Adding doubling limit for cell growth currently using lower bound value 8 hours for Corneal Epithelium
#           * Setting Basal cells to divide on the horizontal axis since that is the worst case scenario for maintaing the simulation 
#           * Cell death only starts after 120 mcs or 12 hours, to allow better equilibration time for simulation initialization   
#         
#       * Development of new capabilities  
#           * Adding calculation of area of contact between cells
#           * Controling EGF_coefdiffusion for cell type and EGF_decay values to be explored in the simulation
#           * Adition of time of simulation control 
#           * Completely different method for data collection and plotting. 
#               Better run time by limiting the use of realtime plotting on CC3D player, maintaining data in ram for faster and correct CVS file format creation.
#               * This allow for better data analysis and plotting later on
#               New plots for:       
#               * Thickness - of wing cells both in the periphery and limbus
#               * Growth - to follow how much each of the growth constraints are contributing to the total growth of the cells
#       
#       TODO:   
#           !New tests should be done for Basal vertical axis division and true random division
#           !Trasform as much of the defined constant parameters into random normal distributions based on their upper and lower limits
#           !Cleanup code and remove unused variables and functions
#           !Add more comments to the code           
#                  
#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

from cc3d.cpp.PlayerPython import * 
from cc3d import CompuCellSetup
from cc3d.core.PySteppables import *
import numpy as np
import random
import csv
from pathlib import Path
from Parameters import *


# GLOBAL PARAMETERS
# RANDOM_SEED = random.seed(12)

    #---- Time Scales ---
HOURtoMCS   = 10
DAYtoMCS    = 240
WEEKtoMCS   = 1680
MONTHtoMCS  = 7300
YEARtoMCS   = 87600

    #---- Spatial Scales ---
UMtoVOXEL   = 2

    #---- Plots Variables ---
DEATHCOUNT  = 0

class ConstraintInitializerSteppable(SteppableBasePy):

# Static Parameters for Simulation Initiation
    current_script_directory = Path(__file__).parent

#   CELLS PARMETERS
#---STEM---
    InitSTEM_LambdaSurface      = InitSTEM_LambdaSurface
    InitSTEM_TargetSurface      = InitSTEM_TargetSurface # 25=pi*r^2 | 5x4 sides

    InitSTEM_LambdaVolume       = InitSTEM_LambdaVolume
    InitSTEM_TargetVolume       = InitSTEM_TargetVolume 

    DensitySTEM_HalfMaxValue    = DensitySTEM_HalfMaxValue
    EGF_STEM_HalfMaxValue       = EGF_STEM_HalfMaxValue

    InitSTEM_LambdaChemo        = InitSTEM_LambdaChemo  #10000  
    
    # Growth Scalars    
    EGF_GrowthScalar_STEM       = EGF_GrowthScalar_STEM #  
    DensityGrowthScalar_STEM    = DensityGrowthScalar_STEM # 

#---BASAL---
    InitBASAL_LambdaSurface     = InitBASAL_LambdaSurface
    InitBASAL_TargetSurface     = InitBASAL_TargetSurface # 25=pi*r^2 = 18 | 5x4 sides =20
    
    InitBASAL_LambdaVolume      = InitBASAL_LambdaVolume
    InitBASAL_TargetVolume      = InitBASAL_TargetVolume

    InitBASAL_LambdaChemo       = InitBASAL_LambdaChemo #125
    InitBASAL_Division          = InitBASAL_Division

    DensityBASAL_HalfMaxValue   = DensityBASAL_HalfMaxValue # 15.3 is very close to the correct value
    EGF_BASAL_HalfMaxValue      = EGF_BASAL_HalfMaxValue

    # Growth Scalars 
    EGF_GrowthScalar_BASAL      = EGF_GrowthScalar_BASAL  #  
    DensityGrowthScalar_BASAL   = DensityGrowthScalar_BASAL   # 

#---WING---    
    InitWING_LambdaSurface      = InitWING_LambdaSurface  
    InitWING_TargetSurface      = InitWING_TargetSurface

    InitWING_LambdaVolume       = InitWING_LambdaVolume
    InitWING_TargetVolume       = InitWING_TargetVolume

    InitWING_EGFLambdaChemo     = InitWING_EGFLambdaChemo

#---SUPERFICIAL---
    InitSUPER_LambdaSurface     = InitSUPER_LambdaSurface
    InitSUPER_TargetSurface     = InitSUPER_TargetSurface

    InitSUPER_LambdaVolume      = InitSUPER_LambdaVolume
    InitSUPER_TargetVolume      = InitSUPER_TargetVolume

    EGF_SUPERDiffCoef           = EGF_SUPERDiffCoef    

    SloughProbability           = SloughProbability

    # Death Scalars
    DeathTimeScalar             = DeathTimeScalar
    DeathVolumeScalar           = DeathVolumeScalar
    SloughScalar                = SloughScalar 
# FIELDS
    MovementBias                = object()
    MovementBiasScreteAmount    = MovementBiasScreteAmount    
    MovementBiasUptake          = MovementBiasUptake

    EGF_Field                   = object()

    EGF_ScreteAmount            = EGF_ScreteAmount   

    EGF_FieldUptakeBASAL        = EGF_FieldUptakeBASAL
    EGF_FieldUptakeSTEM         = EGF_FieldUptakeSTEM
    EGF_FieldUptakeSuper        = EGF_FieldUptakeSuper
    EGF_GlobalDecay             = EGF_GlobalDecay
    

# WOUND
    InjuryType                  = InjuryType
    IsInjury                    = IsInjury
    InjuryTime                  = InjuryTime
    
    #---INJURY AREA---
    InjuryX_Center              = InjuryX_Center 
    InjuryY_Center              = InjuryY_Center
    InjuryRadius                = InjuryRadius

# DEBUGGING
    GrowthControl               = GrowthControl
    MitosisControl              = MitosisControl
    DeathControl                = DeathControl
    DifferentiationControl      = DifferentiationControl

#   PLOTS
    CC3D_PLOT                   = True
    CellCount                   = CellCount
    XHypTracker                 = XHypTracker 
    SloughTracker               = SloughTracker
    PressureTracker             = PressureTracker
    PressurePlot                = PressurePlot
    VolumeTracker               = VolumeTracker
    EGF_SeenByCell              = EGF_SeenByCell
    CenterBiasPlot              = CenterBiasPlot
    CenterBias                  = CenterBias 
    DivisionTracker             = DivisionTracker
    GrowthPlot                  = GrowthPlot
    ThicknessPlot               = ThicknessPlot
    VolumeSurfaceDetail         = VolumeSurfaceDetail
    MitosisPlot                 = MitosisPlot

# TIME OF SIMULATION
    SimTime                     = SimTime

    def __init__(self,frequency=1):
        SteppableBasePy.__init__(self,frequency)       

        #  Checking for which if anny plot should be displayed
        if self.PressureTracker:
            self.track_cell_level_scalar_attribute(field_name='Pressure', attribute_name='Pressure')
        if self.VolumeTracker:
            self.track_cell_level_scalar_attribute(field_name='Volume', attribute_name='Volume')
        if self.EGF_SeenByCell:
            self.track_cell_level_scalar_attribute(field_name='EGF_Seen', attribute_name='EGF')
        if self.CenterBias:
            self.track_cell_level_scalar_attribute(field_name='CenterBias', attribute_name='CenterBias')
        if self.DivisionTracker:
            self.track_cell_level_scalar_attribute(field_name='DivisionCount', attribute_name='DivisionCount')

        # self.MovementBias =self.get_field_secretor("BASALMVBIAS")
        # self.TearField = self.get_field_secretor("TEAR")

    def start(self):         
       
        # FIELDS
        self.MovementBias =self.get_field_secretor("BASALMVBIAS")
        self.EGF_Field = self.get_field_secretor("EGF")
        
        # STROMA
        for cell in self.cell_list_by_type(self.STROMA):

            cell.targetVolume = 5000
            cell.lambdaVolume = 60.0          
          
        # MEMBANE    
        for cell in self.cell_list_by_type(self.MEMB):
            
            cell.targetVolume = 1
            cell.lambdaVolume = 1000.0
            cell.lambdaSurface = 1.0
            cell.targetSurface = 100.0

        # LIMB    
        for cell in self.cell_list_by_type(self.LIMB):
            
            cell.targetVolume = 1
            cell.lambdaVolume = 1000.0
            cell.lambdaSurface = 1.0
            cell.targetSurface = 100.0

        # KERATO    
        # for cell in self.cell_list_by_type(self.KERATO):

        #     cell.targetVolume = 35
        #     cell.lambdaVolume = 2.0
            # cell.lambdaSurface = 2.0
            # cell.targetSurface = 35.0

        # STEM
        for cell in self.cell_list_by_type(self.STEM):

            cell.targetVolume = self.InitSTEM_TargetVolume
            cell.lambdaVolume = self.InitSTEM_LambdaVolume
            cell.lambdaSurface = self.InitSTEM_LambdaSurface
            cell.targetSurface = self.InitSTEM_TargetSurface

            cell.dict["LambdaChemo"] = self.InitSTEM_LambdaChemo
            ChemotaxisData = self.chemotaxisPlugin.addChemotaxisData(cell, "BASALMVBIAS")
            ChemotaxisData.setLambda(cell.dict["LambdaChemo"])           
           
        # BASAL          
        for cell in self.cell_list_by_type(self.BASAL):

            cell.targetVolume = self.InitBASAL_TargetVolume
            cell.lambdaVolume = self.InitBASAL_LambdaVolume
            cell.lambdaSurface = self.InitBASAL_LambdaSurface
            cell.targetSurface = self.InitBASAL_TargetSurface

            # cell.dict['DivisionCount'] = random.randint(0, self.InitBASAL_Division)
            cell.dict['DivisionCount'] = self.InitBASAL_Division

            cell.dict["LambdaChemo"] = self.InitBASAL_LambdaChemo
            ChemotaxisData = self.chemotaxisPlugin.addChemotaxisData(cell, "BASALMVBIAS")            
            ChemotaxisData.setLambda(cell.dict["LambdaChemo"])

        # WING  
        for cell in self.cell_list_by_type(self.WING):

            cell.targetVolume = self.InitWING_TargetVolume
            cell.lambdaVolume = self.InitWING_LambdaVolume
            cell.lambdaSurface = self.InitWING_LambdaSurface
            cell.targetSurface = self.InitWING_TargetSurface

            cell.dict['EGF_LambdaChemo'] = self.InitWING_EGFLambdaChemo
            ChemotaxisData = self.chemotaxisPlugin.addChemotaxisData(cell, "EGF")            
            ChemotaxisData.setLambda(cell.dict["EGF_LambdaChemo"])

            
        # SUPERFICIAL
        for cell in self.cell_list_by_type(self.SUPER):
            cell.targetVolume = self.InitSUPER_TargetVolume
            cell.lambdaVolume = self.InitSUPER_LambdaVolume

            cell.lambdaSurface = self.InitSUPER_LambdaSurface
            cell.targetSurface = self.InitSUPER_TargetSurface
            cell.dict['Slough'] = self.SloughProbability

        self.get_xml_element("EGF_Coef_SUPER").cdata = self.EGF_SUPERDiffCoef
        self.get_xml_element("EGF_GlobalDecay").cdata = self.EGF_GlobalDecay  

        # Leaving it here for DEBUGGING 
        # for cell in self.cell_list_by_type(self.SUPER):
        #     self.delete_cell(cell)      

        # for x in range(0,197,5):
        #     self.cell_field[x:x+5, 55:60, 0] = self.new_cell(self.SUPER)

    def step(self, mcs):
        
        # for cell in self.cell_list_by_type(self.WING):
        #     self.delete_cell(cell)
        # for cell in self.cell_list_by_type(self.SUPER):
        #     self.delete_cell(cell)
        # for cell in self.cell_list_by_type(self.WALL):
        # #    self.TearField.secreteOutsideCellAtBoundaryOnContactWith(cell, 1, [self.MEDIUM])

        #     if cell.yCOM == 79:
        #         self.TearField.secreteOutsideCellAtBoundaryOnContactWith(cell, 10, [self.MEDIUM])

        # ---- CELL PARAMETERS UPDATE ----        
        for cell in self.cell_list_by_type(self.BASAL, self.STEM, self.WING, self.SUPER):           
            cell.dict['Pressure']   = abs(cell.pressure)
            cell.dict['Volume']     = cell.volume
            cell.dict['EGF']        = self.EGF_Field.amountSeenByCell(cell)

            # UPTAKE OF FIELD
                # BASAL
            if cell.type == self.BASAL:              
                cell.dict['Bias_Uptake'] = self.MovementBias.uptakeInsideCellTotalCount(cell, self.MovementBiasUptake, self.MovementBiasUptake).tot_amount            
                
                cell.dict['EGF_Uptake'] =+ abs(self.EGF_Field.uptakeInsideCellTotalCount(cell,self.EGF_FieldUptakeBASAL, self.EGF_FieldUptakeBASAL).tot_amount)            
                self.EGF_Field.uptakeInsideCellAtBoundary(cell, self.EGF_FieldUptakeBASAL, self.EGF_FieldUptakeBASAL)

                # print(cell.dict['EGF_Uptake'], "BASAL", cell.id)
             
                # STEM 
            elif cell.type == self.STEM:
                cell.dict['EGF_Uptake'] =+ abs(self.EGF_Field.uptakeInsideCellTotalCount(cell,self.EGF_FieldUptakeSTEM, self.EGF_FieldUptakeSTEM/2).tot_amount)            

                self.EGF_Field.uptakeInsideCellAtBoundary(cell, self.EGF_FieldUptakeSTEM, self.EGF_FieldUptakeSTEM)               
               
                # print(cell.dict['EGF_Uptake'], "STEM", cell.id)

            elif cell.type == self.SUPER: # why did I put this here?
                cell.dict['EGF_Uptake'] = self.EGF_Field.uptakeInsideCellTotalCount(cell,self.EGF_FieldUptakeSuper, self.EGF_FieldUptakeSuper/10).tot_amount            

                self.EGF_Field.uptakeInsideCellAtBoundary(cell, self.EGF_FieldUptakeSuper, self.EGF_FieldUptakeSuper) 
         
        # print("--------------------")    
        # for cell in self.cell_list_by_type(self.MEMB):

        #     self.MovementBias.secreteOutsideCellAtBoundaryOnContactWith(cell, self.MovementBiasValue, [self.MEDIUM, self.WING, self.SUPER])

class GrowthSteppable(SteppableBasePy):

    def __init__(self,frequency=1):

        SteppableBasePy.__init__(self, frequency)

        self.GrowthControl = ConstraintInitializerSteppable.GrowthControl        

        self.DensityBASAL_HalfMaxValue = ConstraintInitializerSteppable.DensityBASAL_HalfMaxValue
        self.DensitySTEM_HalfMaxValue = ConstraintInitializerSteppable.DensitySTEM_HalfMaxValue

        self.EGF_BASAL_HalfMaxValue = ConstraintInitializerSteppable.EGF_BASAL_HalfMaxValue
        self.EGF_STEM_HalfMaxValue = ConstraintInitializerSteppable.EGF_STEM_HalfMaxValue
        
        self.EGF_GrowthScalar_BASAL = ConstraintInitializerSteppable.EGF_GrowthScalar_BASAL
        self.EGF_GrowthScalar_STEM = ConstraintInitializerSteppable.EGF_GrowthScalar_STEM

        self.DensityGrowthScalar_BASAL = ConstraintInitializerSteppable.DensityGrowthScalar_BASAL
        self.DensityGrowthScalar_STEM = ConstraintInitializerSteppable.DensityGrowthScalar_STEM

        self.BASAL_initial_volume = ConstraintInitializerSteppable.InitBASAL_TargetVolume
        self.STEM_initial_volume = ConstraintInitializerSteppable.InitSTEM_TargetVolume

        self.BASAL_doubling = (1/(HOURtoMCS*8)*self.BASAL_initial_volume) # lower bound
        self.STEM_doubling = (1/(HOURtoMCS*8)*self.STEM_initial_volume)
  
    def step(self, mcs):              
        
        if self.GrowthControl:
            # ---- BASAL ----
            for cell in self.cell_list_by_type(self.BASAL):
            
                # GROWTH THROUGH EGF
                cell.dict['EGF_Growth'] = ((cell.dict['EGF']**4/(self.EGF_BASAL_HalfMaxValue**4 + cell.dict['EGF']**4))) # Hill Promoter EGF                
            
                # GROWTH CONTACT INHIBITION  
                cell.dict['DensityGrowth'] = ((self.DensityBASAL_HalfMaxValue**4/(self.DensityBASAL_HalfMaxValue**4 + cell.dict['Pressure']**4))) # Hill Inhibitor Pressure        

                # TOTAL GROWTH
                cell.dict['TotalGrowth'] = (self.BASAL_doubling * (cell.dict['DensityGrowth'] * cell.dict['EGF_Growth'])**2 /
                                             1**2 + (cell.dict['DensityGrowth'] * cell.dict['EGF_Growth'])**2)                  
                
                cell.targetVolume += cell.dict['TotalGrowth']
           
            # ---- STEM ----
            for cell in self.cell_list_by_type(self.STEM):
                
                # GROWTH THROUGH TEAR 
                cell.dict['EGF_Growth'] = ((cell.dict['EGF']**4/(self.EGF_STEM_HalfMaxValue**4 + cell.dict['EGF']**4))) # Hill Promoter EGF                

                # GROWTH CONTACT INHIBITION 
                cell.dict['DensityGrowth'] = ((self.DensitySTEM_HalfMaxValue**4/(self.DensitySTEM_HalfMaxValue**4 + cell.dict['Pressure']**4))) # Hill Inhibitor Pressure                        
                
                # TOTAL GROWTH
                cell.dict['TotalGrowth'] = (self.STEM_doubling * (cell.dict['DensityGrowth'] * cell.dict['EGF_Growth'])**2 /
                                             1**2 + (cell.dict['DensityGrowth'] * cell.dict['EGF_Growth'])**2)
                
                cell.targetVolume += cell.dict['TotalGrowth']

class MitosisSteppable(MitosisSteppableBase):
    
    STEM_to_divide = 0
    BASAL_to_divide = 0
    cells_divided_count_history = {'Basal' : [],
                                   'Stem' : [] }  # to store the counts at each MCS

    def __init__(self,frequency=1):
        MitosisSteppableBase.__init__(self,frequency)
        
        self.MitosisControl = ConstraintInitializerSteppable.MitosisControl
        # self.STEM_to_divide = 0
        # self.BASAL_to_divide = 0

    def step(self, mcs):

        if self.MitosisControl:

            self.STEM_to_divide = 0
            self.BASAL_to_divide = 0
            cells_to_divide=[]

            #---BASAL---            
            for cell in self.cell_list_by_type(self.BASAL):
                if cell.volume>50:
                    cells_to_divide.append(cell)
            
            #---STEM---
            for cell in self.cell_list_by_type(self.STEM):
                if cell.volume>50:
                    cells_to_divide.append(cell)

            #---WING---
            # for cell in self.cell_list_by_type(self.WING):
            #     if cell.volume>25:
            #         cell.targetVolume = 25
                    # cells_to_divide.append(cell)
                        


            # Division orientation 
            for cell in cells_to_divide:
                NEIGHBOR_DICT = self.get_cell_neighbor_data_list(cell).neighbor_count_by_type()                        
                    
                # STEM 
                if cell.type == self.STEM:
                    self.STEM_to_divide += 1
                    self.set_parent_child_position_flag(1)
                    if (self.MEMB in NEIGHBOR_DICT.keys()):
                        self.divide_cell_orientation_vector_based(cell,1,0,0) # Orientation for Stem division (Vertical)
                    else:
                        self.divide_cell_random_orientation(cell)
               
                # BASAL cell constraints in division    
                elif cell.type == self.BASAL:
                    self.BASAL_to_divide += 1
                    if cell.dict['DivisionCount'] > 0:
                        cell.dict['DivisionCount'] -= 1                
                        self.divide_cell_random_orientation(cell)
                        # self.divide_cell_orientation_vector_based(cell,1,0,0)
                    else:
                        cell.type = self.WING # Differentiation Due to Last Division
                        cell.lambdaSurface = ConstraintInitializerSteppable.InitWING_LambdaSurface
                        cell.targetSurface = ConstraintInitializerSteppable.InitWING_TargetSurface                     
                        # self.divide_cell_random_orientation(cell)  
                # WING
                # else: 
                #     self.divide_cell_random_orientation(cell)
            
            # cells_to_divide = []        
                # Other valid options
                # self.divide_cell_orientation_vector_based(cell,1,1,0)
                # self.divide_cell_along_major_axis(cell)
                # self.divide_cell_along_minor_axis(cell)

            # Append counts to history list
            self.cells_divided_count_history['Basal'].append(self.BASAL_to_divide)
            self.cells_divided_count_history['Stem'].append(self.STEM_to_divide)

    def get_divided_cells_counts(self):
        return self.cells_divided_count_history

    def clear_divided_cells_counts(self):       
        for key in self.cells_divided_count_history:
           self.cells_divided_count_history[key] = []
        

    def update_attributes(self):        
        self.parent_cell.targetVolume /= 2.0            
        self.clone_parent_2_child()            

class DeathSteppable(SteppableBasePy):

    def __init__(self, frequency=1):

        SteppableBasePy.__init__(self, frequency)

        # CONSTANT CELL DEATH RATE CONTROLLER
        self.DeathControl = ConstraintInitializerSteppable.DeathControl

        # WOUND EVENT 
        self.injuryType = ConstraintInitializerSteppable.InjuryType
        self.injury     = ConstraintInitializerSteppable.IsInjury     
        self.injuryTime = ConstraintInitializerSteppable.InjuryTime

        # INJURY AREA
        self.x_center = ConstraintInitializerSteppable.InjuryX_Center 
        self.y_center = ConstraintInitializerSteppable.InjuryY_Center
        self.radius   = ConstraintInitializerSteppable.InjuryRadius

        # DEATH SCALER 
        self.deathTimeScalar    = ConstraintInitializerSteppable.DeathTimeScalar
        self.deathVolumeScalar  = ConstraintInitializerSteppable.DeathVolumeScalar                
        self.sloughScalar       = ConstraintInitializerSteppable.SloughScalar

        self.SUPER_TargetVolume = ConstraintInitializerSteppable.InitSUPER_TargetVolume
        self.SUPER_TargetSurface = ConstraintInitializerSteppable.InitSUPER_TargetSurface
     
    def step(self, mcs):
        global DEATHCOUNT

        # --- Basal cell death rate ---
        if self.injury:

            # if (mcs > 10 and mcs < 300): # sustained injury         
            if (mcs == self.injuryTime): 

                # ---- PRODUCTION OF TEAR ----
                # self.get_xml_element('Tear_GlbDiff').cdata = 60 # GlobalDiffusion change

                cells_to_kill = set()
                for i in np.arange(-self.radius,  self.radius, 1.0):
                    y_position = self.y_center + i
                    if not (self.y_center + i < 0 or self.y_center + i> self.dim.y):
                        for j in np.arange(-np.sqrt(self.radius**2-i**2), np.sqrt(self.radius**2-i**2),1.0):
                            x_position = self.x_center + j
                            if not (self.x_center + j < 0 or self.x_center + j> self.dim.x):
                                cell = self.cell_field[x_position, y_position, 0]
                                if cell:
                                    cells_to_kill.add(cell.id)
                # print(cells_to_kill)
                # print("###############\n", "Cell that will die on injury", len(cells_to_kill))

                if self.injuryType == 'A':
                    # --- ABRASION ---    
                    for cellid in cells_to_kill:
                        cell = self.fetch_cell_by_id(cellid)                
                        self.delete_cell(cell)
                else:
                    # --- CHEMICAL ---    
                    for cellid in cells_to_kill:
                        cell = self.fetch_cell_by_id(cellid)                
                        cell.targetVolume = 0
                        cell.lambdaVolume = 100
                        
            # # ---- PRODUCTION OF TEAR ----           
            # if (mcs > (self.injuryTime + HOURtoMCS * 10)):
            #     self.get_xml_element('Tear_GlbDiff').cdata = 10                
         
        # --- Constant Cell Death Rate --- 
        if self.DeathControl: 
            if mcs > 120:   
                deathsum = 0            
                for cell in self.cell_list_by_type(self.SUPER):
                    NEIGHBOR_DICT = self.get_cell_neighbor_data_list(cell).neighbor_count_by_type()  # Checking Neighbors                        
                                    
                    #  Considering only incontact with MEDIUM
                    if (self.MEDIUM in NEIGHBOR_DICT.keys()) : # Only cells at the edge has basal death for now
                            
                        # cell.targetVolume -= (1/(WEEKtoMCS)) * 25 * (random.random()) # adding a stochastic death rate
                        # cell.targetVolume -= (1/(WEEKtoMCS)) * self.sloughScalar

                        # On average the whole process for cell turnover takes 7-14 days. 
                        # If we assume that the average layer cells are 4 cells underneath the superficial
                        # at the periphery and  it can take between 1.75 days/layer. If that is maintained
                        # and we consider a 8 layer cells in the limbus we still have the 1.75 rate for the 14 days. 
                        cell.targetVolume -= ((self.SUPER_TargetVolume/1.75)/(DAYtoMCS)) * self.sloughScalar
                        cell.targetSurface -= ((self.SUPER_TargetSurface/1.75)/(DAYtoMCS)) * 0.25  
                                            

                        if cell.volume < 15: # Minimum Cell Size Before Slough | if no slogh cell will disapear in 672 MCS ~3 days(2.8)
                            cell.dict['Slough'] = (1 - np.exp(1/(-(HOURtoMCS*self.deathTimeScalar)*(cell.volume*self.deathVolumeScalar))))
                        
                            if (random.random() < cell.dict['Slough']):
                                deathsum += 1
                                self.delete_cell(cell)
                DEATHCOUNT = deathsum

            for cell in self.cell_list_by_type(self.SUPER, self.WING, self.BASAL, self.STEM):
                if (cell.volume < 3):
                    self.delete_cell(cell)
  
class DifferentiationSteppable(SteppableBasePy):

    def __init__(self, frequency=1):      
        SteppableBasePy.__init__(self, frequency)  

        self.DifferentiationControl = ConstraintInitializerSteppable.DifferentiationControl    
      
    def start(self):
       pass
   
    def step(self, mcs): 
      
        if self.DifferentiationControl:

            if mcs > 120:
                # SUPER CELL
                # for cell in self.cell_list_by_type(self.SUPER):
                #     # PARAMETERS
                #     NEIGHBOR_DICT = self.get_cell_neighbor_data_list(cell).neighbor_count_by_type()
                            
                # WING CELL
                for cell in self.cell_list_by_type(self.WING):
                    # PARAMETERS
                    NEIGHBOR_DICT = self.get_cell_neighbor_data_list(cell).neighbor_count_by_type()
                                
                    # DIFFERENTIATION RULES
                    if (self.MEDIUM in NEIGHBOR_DICT.keys()) and (self.WING in NEIGHBOR_DICT.keys()):
                        cell.type = self.SUPER
                        self.initializeDifferentiatedCell(cell)

                    # Tear based Differention
                    # elif (cell.dict['Tear']> 5) and (self.SUPER in NEIGHBOR_DICT.keys()): # Amount seen by cell = cell.dict['Tear']
                    #     cell.type = self.SUPER
                    # elif (self.TEAR in NEIGHBOR_DICT.keys()) and (self.WING in NEIGHBOR_DICT.keys()): # Old version
                        # cell.type = self.SUPER
                        
                # BASAL CELL
                for cell in self.cell_list_by_type(self.BASAL):
                    # Fetch the boundary pixels of the current Basal cell
                    boundary_pixel_list = self.get_cell_boundary_pixel_list(cell)
                    
                    memb_contact_area = 0  # Counter for the area of contact with MEMB cells
                    
                    # Iterate through each boundary pixel to check for MEMB neighbors
                    for boundary_pixel_tracker_data in boundary_pixel_list:
                        boundary_pixel = boundary_pixel_tracker_data.pixel
                        for neighbor_pixel in self.get_pixel_neighbors(boundary_pixel):
                            if neighbor_pixel.type == self.MEMB:
                                memb_contact_area += 1  # Increment the counter for every MEMB neighbor pixel

                    # Differentiation Rules: The cell will differentiate if its contact area with MEMB is 2 or less
                    if memb_contact_area <= 5:
                        cell.type = self.WING
                        self.initializeDifferentiatedCell(cell)
                      
                                
                # STEM CELL
                for cell in self.cell_list_by_type(self.STEM):
                    # Fetch the boundary pixels of the current Basal cell
                    boundary_pixel_list = self.get_cell_boundary_pixel_list(cell)
                    
                    memb_contact_area = 0  # Counter for the area of contact with MEMB cells
                    
                    # Iterate through each boundary pixel to check for MEMB neighbors
                    for boundary_pixel_tracker_data in boundary_pixel_list:
                        boundary_pixel = boundary_pixel_tracker_data.pixel
                        for neighbor_pixel in self.get_pixel_neighbors(boundary_pixel):
                            if neighbor_pixel.type == self.LIMB:
                                memb_contact_area += 1  # Increment the counter for every MEMB neighbor pixel

                    # Differentiation Rules: The cell will differentiate if its contact area with MEMB is 2 or less
                    if memb_contact_area <= 5:
                        cell.type = self.BASAL
                        self.initializeDifferentiatedCell(cell)

    def get_pixel_neighbors(self, pixel):
        # Extracting x, y, z coordinates of the pixel
        x, y, z = pixel.x, pixel.y, pixel.z
        neighbors = []
        
        # Define possible relative coordinates for neighbors in a 2D grid
        # relative_coords = [(-1, -1, 0), (-1, 0, 0), (-1, 1, 0), 
        #                     (0, -1, 0),              (0, 1, 0), 
        #                     (1, -1, 0), (1, 0, 0), (1, 1, 0)]

        # Faster version with less neighbors to check (since it is basally constrained)
        relative_coords = [
                            (0, -1, 0),              (0, 1, 0), 
                            (1, -1, 0), (1, 0, 0), (1, 1, 0)]
        
        for dx, dy, dz in relative_coords:
            # Check if the neighboring cell type matches our desired cell type
            neighbor_pixel = self.cellField[x + dx, y + dy, z + dz]
            if neighbor_pixel:
                neighbors.append(neighbor_pixel)
        
        return neighbors                   

    def initializeDifferentiatedCell(self, cell):   
        
        # --- BASAL ---
        if cell.type == self.BASAL:
            cell.lambdaSurface = ConstraintInitializerSteppable.InitBASAL_LambdaSurface
            cell.targetSurface = ConstraintInitializerSteppable.InitBASAL_TargetSurface
            cell.dict["LambdaChemo"] = ConstraintInitializerSteppable.InitBASAL_LambdaChemo
            cell.dict['DivisionCount'] = ConstraintInitializerSteppable.InitBASAL_Division
            ChemotaxisData = self.chemotaxisPlugin.addChemotaxisData(cell, "BASALMVBIAS")            
            ChemotaxisData.setLambda(cell.dict["LambdaChemo"])

        # --- STEM ---
        elif cell.type == self.STEM:
            cell.lambdaSurface = ConstraintInitializerSteppable.InitSTEM_LambdaSurface
            cell.targetSurface = ConstraintInitializerSteppable.InitSTEM_TargetSurface
            cell.dict["LambdaChemo"] = ConstraintInitializerSteppable.InitSTEM_LambdaChemo
            ChemotaxisData = self.chemotaxisPlugin.addChemotaxisData(cell, "BASALMVBIAS")            
            ChemotaxisData.setLambda(cell.dict["LambdaChemo"])

        # --- WING ---
        elif cell.type == self.WING:
            cell.lambdaSurface = ConstraintInitializerSteppable.InitWING_LambdaSurface
            cell.targetSurface = ConstraintInitializerSteppable.InitWING_TargetSurface
            cell.lamdaVolume = ConstraintInitializerSteppable.InitWING_LambdaVolume
            cell.targetVolume = ConstraintInitializerSteppable.InitWING_TargetVolume
            cell.dict["LambdaChemo"] = 0

        # --- SUPER ---
        elif cell.type == self.SUPER:
            cell.lambdaSurface = ConstraintInitializerSteppable.InitSUPER_LambdaSurface
            cell.targetSurface = ConstraintInitializerSteppable.InitSUPER_TargetSurface
            cell.dict["LambdaChemo"] = 0

class SecretionSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)

        self.EGF_SecreteAmount          = ConstraintInitializerSteppable.EGF_ScreteAmount
        self.MovementBiasScreteAmount   = ConstraintInitializerSteppable.MovementBiasScreteAmount
    
    def start(self):
        self.MovementBias   = self.get_field_secretor("BASALMVBIAS")
        self.EGF_Field      = self.get_field_secretor("EGF")

    def step(self, mcs):

        for cell in self.cell_list_by_type(self.MEMB):
            self.MovementBias.secreteOutsideCellAtBoundaryOnContactWith(cell, self.MovementBiasScreteAmount, [self.MEDIUM, self.WING, self.SUPER])

        # for cell in self.cell_list_by_type(self.WALL):
            # Only WALL at the top edge of the simulation
            # if cell.yCOM == 68: # using hex latice             
            # if cell.yCOM > 78: 
            #     self.EGF_Field.secreteOutsideCellAtBoundaryOnContactWith(cell, self.EGF_SecreteAmount, [self.MEDIUM])

class PlotSteppable(SteppableBasePy):

    def __init__(self, frequency=10):      
        SteppableBasePy.__init__(self, frequency)

        self.CC3D_PLOT      = ConstraintInitializerSteppable.CC3D_PLOT 
        self.mitosis = MitosisSteppable()

        self.cellCount      = ConstraintInitializerSteppable.CellCount
        self.xhypTracker    = ConstraintInitializerSteppable.XHypTracker
        self.sloughTracker  = ConstraintInitializerSteppable.SloughTracker
        self.centerBiasPlot = ConstraintInitializerSteppable.CenterBiasPlot
        self.pressurePlot   = ConstraintInitializerSteppable.PressurePlot
        self.growthPlot     = ConstraintInitializerSteppable.GrowthPlot
        self.thicknessPlot  = ConstraintInitializerSteppable.ThicknessPlot
        self.VolumeSurfaceDetailPlot = ConstraintInitializerSteppable.VolumeSurfaceDetail
        self.MitosisPlot    = ConstraintInitializerSteppable.MitosisPlot

        self.EGF_GrowthScalar_STEM = ConstraintInitializerSteppable.EGF_GrowthScalar_STEM
        self.EGF_DensityScalar_STEM = ConstraintInitializerSteppable.DensityGrowthScalar_STEM

        self.EGF_GrowthScalar_BASAL = ConstraintInitializerSteppable.EGF_GrowthScalar_BASAL
        self.EGF_DensityScalar_BASAL = ConstraintInitializerSteppable.DensityGrowthScalar_BASAL

        self.BASAL_to_divide = MitosisSteppable.BASAL_to_divide
        self.STEM_to_divide = MitosisSteppable.STEM_to_divide

        self.current_script_directory = Path(__file__).parent

        self.SimTime = ConstraintInitializerSteppable.SimTime

        self.data_points =  {
            'cell_count': {
                "Superficial": [],
                "Wing": [],
                "Basal": [],
                "Stem": [],
                "Time": []
            },

            'slough_tracker': {
                "Wing": [],
                "Individual Volume": [],
                "Average Vol": [],
                "Slough Count": [],
                "Time": []
            },

            'center_bias': {
                "Wing": [],
                "Individual Volume": [],
                "Average Vol": [],
                "Slough Count": [],
                "Time": []
            },

            'pressure': {
                "Stem Average": [],
                "Basal Average": [],
                "Time": []
                # "Basal Individual": [],
                # "Wings": [],
                # "Average Vol": [],
                # "Slough Count": []
            },

            'growth': {
                "Stem Average Total": [],
                "Basal Average Total": [],
                "Stem Average EGF": [],
                "Basal Average EGF": [],
                "Stem Average Density": [],
                "Basal Average Density": [],
                "Stem Average TargetVolume": [],
                "Basal Average TargetVolume": [],
                "Stem Average EGFseen" : [],
                "Basal Average EGFseen" : [],
                "Stem Average Pressure" : [],
                "Basal Average Pressure" : [],
                "Time": []
            },

            'thickness': {
                "Limbus Avg Thickness" : [],
                "Limbus Left Thickness" : [],
                "Limbus Right Thickness" : [],
                "Periphery Avg Thickness" : [],
                "Periphery Left Thickness" : [],
                "Periphery Right Thickness" : [],
                "Time": []
            },

            'volume_surface': {
                'Superficial Volume': [],
                'Superficial Surface': [],
                'Wing Volume': [],
                'Wing Surface': [],
                'Basal Volume': [],
                'Basal Surface': [],                
                'Stem Volume': [],
                'Stem Surface': [],
                'Time': []
            },
            
            'mitosis': {
                'BASAL': [],
                'STEM': [],
                'Time': []
            }
                            }

    def start(self):

        if self.CC3D_PLOT:
        # ---- Cell Count Plot ----
            if self.cellCount:            
            
                self.plot_cell_count = self.add_new_plot_window(title='Cell Count by Type',
                                                                x_axis_title='Time (Hours)',
                                                                y_axis_title='Number of Cells',
                                                                x_scale_type='linear', y_scale_type='linear',
                                                                grid=True, config_options={'legend':True})

                self.plot_cell_count.add_plot("Superficial", style='Lines', color='cyan')
                self.plot_cell_count.add_plot("Wing", style='Lines', color='blue')
                self.plot_cell_count.add_plot("Basal", style='Lines', color='pink')
                self.plot_cell_count.add_plot("Stem", style='Lines', color='red')

            # ---- Slough Tracker Plot ----
            if self.sloughTracker:
                
                self.plot_slough_tracker = self.add_new_plot_window(title='Slough Tracker Tracker',
                                                                x_axis_title='Time (Hour)',
                                                                y_axis_title='Volume(Cyan and Red) / # Cells(Green)',
                                                                x_scale_type='linear', y_scale_type='linear',
                                                                grid=True, config_options={'legend':True})

                self.plot_slough_tracker.add_plot("Wing", style='Dots', color='blue')
                self.plot_slough_tracker.add_plot("Individual Volume", style='Dots', color='Cyan')
                self.plot_slough_tracker.add_plot("Average Vol", style='lines', color='red')
                self.plot_slough_tracker.add_plot("Slough Count", style='lines', color='green')
            
            # ---- Chemochine Center Movement Bias Plot ----
            if self.centerBiasPlot:
                
                self.plot_center_bias = self.add_new_plot_window(title='Chemochine Center Movement Bias',
                                                                x_axis_title='Time (Hour)',
                                                                y_axis_title='Strenth of the Signal',
                                                                x_scale_type='linear', y_scale_type='linear',
                                                                grid=True, config_options={'legend':True})

                self.plot_center_bias.add_plot("Wing", style='Dots', color='blue')
                self.plot_center_bias.add_plot("Individual Volume", style='Dots', color='Cyan')
                self.plot_center_bias.add_plot("Average Vol", style='lines', color='red')
                self.plot_center_bias.add_plot("Slough Count", style='lines', color='green')

            if self.pressurePlot:
                
                self.plot_pressure = self.add_new_plot_window(title='System Pressure',
                                                                x_axis_title='Time (Hour)',
                                                                y_axis_title='Pressure',
                                                                x_scale_type='linear', y_scale_type='linear',
                                                                grid=True, config_options={'legend':True})

                self.plot_pressure.add_plot("Stem Average", style='lines', color='red')
                self.plot_pressure.add_plot("Basal Average", style='lines', color='pink')
                # self.plot_pressure.add_plot("Basal Individual", style='lines', color='green')
                # self.plot_pressure.add_plot("Wings", style='Dots', color='Cyan')
                # self.plot_pressure.add_plot("Average Vol", style='lines', color='red')
                # self.plot_pressure.add_plot("Slough Count", style='lines', color='green')

            if self.growthPlot:
                
                self.plot_growth = self.add_new_plot_window(title='System Growth',
                                                                x_axis_title='Time (Hour)',
                                                                y_axis_title='x100',
                                                                x_scale_type='linear', y_scale_type='linear',
                                                                grid=True, config_options={'legend':True})

                self.plot_growth.add_plot("Stem Average Total", style='lines', color='red')
                self.plot_growth.add_plot("Basal Average Total", style='lines', color='pink')
                self.plot_growth.add_plot("Stem Average EGF", style='lines', color='grey')
                self.plot_growth.add_plot("Basal Average EGF", style='lines', color='Orange')
                self.plot_growth.add_plot("Stem Average Density", style='lines', color='Purple')
                self.plot_growth.add_plot("Basal Average Density", style='lines', color='Brown')
                self.plot_growth.add_plot("Stem Average TargetVolume", style='lines', color='Magenta')
                self.plot_growth.add_plot("Basal Average TargetVolume", style='lines', color='cyan')

                self.plot_growth.add_plot("Stem Average EGFseen", style='lines', color='green')
                self.plot_growth.add_plot("Basal Average EGFseen", style='lines', color='teal')
                self.plot_growth.add_plot("Stem Average Pressure", style='lines', color='Blue')
                self.plot_growth.add_plot("Basal Average Pressure", style='lines', color='lightblue')


                # self.plot_growth.add_plot("Wing Count", style='lines', color='Blue')

            if self.thicknessPlot:

                self.plot_thickness = self.add_new_plot_window(title='Wing Thickness',
                                                                x_axis_title='Time(Hours)',
                                                                y_axis_title='Average Y position',
                                                                x_scale_type='linear', y_scale_type='linear',
                                                                grid=True, config_options={'legend':True})
                
                self.plot_thickness.add_plot("Limbus Avg Thickness", style='lines', color='Blue')
                self.plot_thickness.add_plot("Limbus Left Thickness", style='lines', color='lightblue')
                self.plot_thickness.add_plot("Limbus Right Thickness", style='lines', color='darkblue')

                self.plot_thickness.add_plot("Periphery Avg Thickness", style='lines', color='red')
                self.plot_thickness.add_plot("Periphery Left Thickness", style='lines', color='pink')
                self.plot_thickness.add_plot("Periphery Right Thickness", style='lines', color='purple')

            if self.VolumeSurfaceDetailPlot:

                self.plot_volume_surface = self.add_new_plot_window(title='Volume Surface Detail',
                                                                x_axis_title='Time(Hours)',
                                                                y_axis_title='Volume Surface',
                                                                x_scale_type='linear', y_scale_type='linear',
                                                                grid=True, config_options={'legend':True})
                
                self.plot_volume_surface.add_plot('Superficial Volume', style='lines', color='cyan')
                self.plot_volume_surface.add_plot('Superficial Surface', style='lines', color='blue')
                self.plot_volume_surface.add_plot('Wing Volume', style='lines', color='pink')
                self.plot_volume_surface.add_plot('Wing Surface', style='lines', color='red')
                self.plot_volume_surface.add_plot('Basal Volume', style='lines', color='green')
                self.plot_volume_surface.add_plot('Basal Surface', style='lines', color='teal')
                self.plot_volume_surface.add_plot('Stem Volume', style='lines', color='purple')
                self.plot_volume_surface.add_plot('Stem Surface', style='lines', color='brown')

            if self.MitosisPlot:

                self.plot_mitosis = self.add_new_plot_window(title='Mitosis',
                                                                x_axis_title='Time(Hours)',
                                                                y_axis_title='Number of Cells',
                                                                x_scale_type='linear', y_scale_type='linear',
                                                                grid=True, config_options={'legend':True})
                
                self.plot_mitosis.add_plot('BASAL', style='lines', color='cyan')
                self.plot_mitosis.add_plot('STEM', style='lines', color='blue')


    def step(self, mcs):
        global DEATHCOUNT
        HOURtoMCS_factor = mcs / HOURtoMCS

        # ---- End of Simulation ----
        if mcs == (self.SimTime + 1) :

            if self.cellCount:                
                self.write_csv_for_category('cell_count', mcs)                
            
            if self.sloughTracker:
                self.write_csv_for_category('slough_tracker', mcs)                

            if self.pressurePlot:
                self.write_csv_for_category('pressure', mcs)                

            if self.growthPlot:
                self.write_csv_for_category('growth', mcs)                

            if self.thicknessPlot:
                self.write_csv_for_category('thickness', mcs)

            if self.VolumeSurfaceDetailPlot:
                self.write_csv_for_category('volume_surface', mcs)
            
            if self.MitosisPlot:
                self.write_csv_for_category('mitosis', mcs)
                
            CompuCellSetup.stop_simulation()

        # ---- Plots ----
        if self.cellCount and mcs % 10 == 0:

            # CSV
            self.data_points['cell_count']['Superficial'].append(len(self.cell_list_by_type(self.SUPER)))
            self.data_points['cell_count']['Wing'].append(len(self.cell_list_by_type(self.WING)))
            self.data_points['cell_count']['Basal'].append(len(self.cell_list_by_type(self.BASAL)))
            self.data_points['cell_count']['Stem'].append(len(self.cell_list_by_type(self.STEM)))            
            self.data_points['cell_count']['Time'].append(HOURtoMCS_factor)

            if self.CC3D_PLOT:
                # CC3D Plots
                # for cell_type, count in self.data_points['cell_count'].items():
                #     self.plot_data_point(self.plot_cell_count, cell_type, HOURtoMCS_factor, count)
                self.plot_cell_count.add_data_point("Superficial", HOURtoMCS_factor, len(self.cell_list_by_type(self.SUPER)))
                self.plot_cell_count.add_data_point("Wing", HOURtoMCS_factor, len(self.cell_list_by_type(self.WING)))
                self.plot_cell_count.add_data_point("Basal", HOURtoMCS_factor, len(self.cell_list_by_type(self.BASAL)))
                self.plot_cell_count.add_data_point("Stem", HOURtoMCS_factor, len(self.cell_list_by_type(self.STEM)))
            
        if self.sloughTracker and mcs%10 == 0:
            
            listVol = []
            
            for cell in self.cell_list_by_type(self.SUPER):
                listVol.append(cell.volume)

            SumVol = sum(listVol)

            listMinVol = [i for i in listVol if i <= 15]

            for i in listMinVol:                
                # CSV
                self.data_points['slough_tracker']['Individual Volume'].append(i)                
                
            AvrVol = SumVol/len(self.cell_list_by_type(self.SUPER))           
            # CSV
            self.data_points['slough_tracker']['Average Vol'].append(AvrVol)
            self.data_points['slough_tracker']['Slough Count'].append(DEATHCOUNT)
            self.data_points['slough_tracker']['Wing'].append(len(self.cell_list_by_type(self.WING)))
            self.data_points['slough_tracker']['Individual Volume'].append(HOURtoMCS_factor)
            elf.data_points['slough_tracker']['Time'].append(HOURtoMCS_factor)

            if self.CC3D_PLOT:
                # CC3D Plots
                # for cell_type, count in self.data_points['slough_tracker'].items():                    
                #     self.plot_data_point(self.plot_slough_tracker, cell_type, HOURtoMCS_factor, count)
                self.plot_slough_tracker.add_data_point("Wing", HOURtoMCS_factor, len(self.cell_list_by_type(self.WING)))
                self.plot_slough_tracker.add_data_point("Individual Volume", HOURtoMCS_factor, len(listMinVol))
                self.plot_slough_tracker.add_data_point("Average Vol", HOURtoMCS_factor, AvrVol)
                self.plot_slough_tracker.add_data_point("Slough Count", HOURtoMCS_factor, DEATHCOUNT)

        if self.pressurePlot and mcs % 100 == 0:
            BASAL_avg_pressure  = []
            STEM_avg_pressure   = []
            for cell in self.cell_list_by_type(self.BASAL,self.STEM):
                if cell.type == self.BASAL:
                    BASAL_avg_pressure.append(abs(cell.pressure))
                    # self.plot_pressure.add_data_point("Basal Individual", mcs/HOURtoMCS, abs(cell.pressure))

                elif cell.type == self.STEM:
                    STEM_avg_pressure.append(abs(cell.pressure))

            # CSV
            self.data_points['pressure']['Stem Average'].append(np.mean(STEM_avg_pressure))
            self.data_points['pressure']['Basal Average'].append(np.mean(BASAL_avg_pressure))
            self.data_points['pressure']['Time'].append(HOURtoMCS_factor)            

            if self.CC3D_PLOT:
                # CC3D Plots
                # for cell_type, count in self.data_points['pressure'].items():
                #     self.plot_data_point(self.plot_pressure, cell_type, HOURtoMCS_factor, count)
                self.plot_pressure.add_data_point("Stem Average", HOURtoMCS_factor, np.mean(STEM_avg_pressure))
                self.plot_pressure.add_data_point("Basal Average", HOURtoMCS_factor, np.mean(BASAL_avg_pressure))

        if self.growthPlot and mcs % 10 == 0:            

            BASAL_avg_EGF_growth  = []
            BASAL_avg_density_growth  = []
            BASAL_avg_total_growth  = []
            BASAL_avg_targetVolume  = []
            BASAL_avg_EGFseen  = []
            BASAL_avg_pressure  = []

            STEM_avg_EGF_growth  = []
            STEM_avg_density_growth  = []
            STEM_avg_total_growth  = []
            STEM_avg_targetVolume  = []
            STEM_avg_EGFseen  = []
            STEM_avg_pressure  = []
                                       

            for cell in self.cell_list_by_type(self.BASAL,self.STEM):
                if cell.type == self.BASAL:
                    BASAL_avg_EGF_growth.append(cell.dict['EGF_Growth'])
                    BASAL_avg_density_growth.append(cell.dict['DensityGrowth'])
                    BASAL_avg_total_growth.append(cell.dict['TotalGrowth'])
                    BASAL_avg_targetVolume.append(cell.targetVolume)
                    BASAL_avg_EGFseen.append(cell.dict['EGF'])
                    BASAL_avg_pressure.append(cell.dict['Pressure'])

                elif cell.type == self.STEM:
                    STEM_avg_EGF_growth.append(cell.dict['EGF_Growth'])
                    STEM_avg_density_growth.append(cell.dict['DensityGrowth'])
                    STEM_avg_total_growth.append(cell.dict['TotalGrowth'])
                    STEM_avg_targetVolume.append(cell.targetVolume)
                    STEM_avg_EGFseen.append(cell.dict['EGF'])
                    STEM_avg_pressure.append(cell.dict['Pressure'])

            # CSV
            self.data_points['growth']['Stem Average Total'].append(np.mean(STEM_avg_total_growth))
            self.data_points['growth']['Basal Average Total'].append(np.mean(BASAL_avg_total_growth))
            self.data_points['growth']['Stem Average EGF'].append(np.mean(STEM_avg_EGF_growth))
            self.data_points['growth']['Basal Average EGF'].append(np.mean(BASAL_avg_EGF_growth))
            self.data_points['growth']['Stem Average Density'].append(np.mean(STEM_avg_density_growth))
            self.data_points['growth']['Basal Average Density'].append(np.mean(BASAL_avg_density_growth))
            self.data_points['growth']['Stem Average TargetVolume'].append(np.mean(STEM_avg_targetVolume))
            self.data_points['growth']['Basal Average TargetVolume'].append(np.mean(BASAL_avg_targetVolume))
            self.data_points['growth']['Stem Average EGFseen'].append(np.mean(STEM_avg_EGFseen))
            self.data_points['growth']['Basal Average EGFseen'].append(np.mean(BASAL_avg_EGFseen))
            self.data_points['growth']['Stem Average Pressure'].append(np.mean(STEM_avg_pressure))
            self.data_points['growth']['Basal Average Pressure'].append(np.mean(BASAL_avg_pressure))
            self.data_points['growth']['Time'].append(HOURtoMCS_factor)

            if self.CC3D_PLOT:
                # for cell_type, count in self.data_points['growth'].items():
                #     self.plot_data_point(self.plot_growth, cell_type, HOURtoMCS_factor, count) 
                self.plot_growth.add_data_point("Stem Average Total", HOURtoMCS_factor, np.mean(STEM_avg_total_growth))
                self.plot_growth.add_data_point("Basal Average Total", HOURtoMCS_factor, np.mean(BASAL_avg_total_growth))
                self.plot_growth.add_data_point("Stem Average EGF", HOURtoMCS_factor, np.mean(STEM_avg_EGF_growth))
                self.plot_growth.add_data_point("Basal Average EGF", HOURtoMCS_factor, np.mean(BASAL_avg_EGF_growth))
                self.plot_growth.add_data_point("Stem Average Density", HOURtoMCS_factor, np.mean(STEM_avg_density_growth))
                self.plot_growth.add_data_point("Basal Average Density", HOURtoMCS_factor, np.mean(BASAL_avg_density_growth))
                self.plot_growth.add_data_point("Stem Average TargetVolume", HOURtoMCS_factor, np.mean(STEM_avg_targetVolume))
                self.plot_growth.add_data_point("Basal Average TargetVolume", HOURtoMCS_factor, np.mean(BASAL_avg_targetVolume))
                self.plot_growth.add_data_point("Stem Average EGFseen", HOURtoMCS_factor, np.mean(STEM_avg_EGFseen))
                self.plot_growth.add_data_point("Basal Average EGFseen", HOURtoMCS_factor, np.mean(BASAL_avg_EGFseen))
                self.plot_growth.add_data_point("Stem Average Pressure", HOURtoMCS_factor, np.mean(STEM_avg_pressure))
                self.plot_growth.add_data_point("Basal Average Pressure", HOURtoMCS_factor, np.mean(BASAL_avg_pressure))
  
        if self.thicknessPlot and mcs % 10 == 0:
            
            thickness_dict = self.compute_thickness_for_all_columns(self.WING)
            left_thicknesses = [thickness for x, thickness in thickness_dict.items() if x <= 110]
            right_thicknesses = [thickness for x, thickness in thickness_dict.items() if x > 110]

            if left_thicknesses:  # Prevent division by zero
                # CSV
                self.data_points['thickness']['Limbus Avg Thickness'].append(np.mean(left_thicknesses))
                self.data_points['thickness']['Limbus Left Thickness'].append(left_thicknesses[0])
                self.data_points['thickness']['Limbus Right Thickness'].append(left_thicknesses[-1])
            else:               
                # CSV
                self.data_points['thickness']['Periphery Avg Thickness'].append(0)                
                self.data_points['thickness']['Periphery Left Thickness'].append(0)
                self.data_points['thickness']['Periphery Right Thickness'].append(0)                

            if right_thicknesses:               
                # CSV
                self.data_points['thickness']['Periphery Avg Thickness'].append(np.mean(right_thicknesses))
                self.data_points['thickness']['Periphery Left Thickness'].append(right_thicknesses[0])
                self.data_points['thickness']['Periphery Right Thickness'].append(right_thicknesses[-1])
            else:
                # CSV
                self.data_points['thickness']['Periphery Avg Thickness'].append(0)
                self.data_points['thickness']['Periphery Left Thickness'].append(0)
                self.data_points['thickness']['Periphery Right Thickness'].append(0)

            self.data_points['thickness']['Time'].append(HOURtoMCS_factor)


            if self.CC3D_PLOT:
                # for cell_type, count in self.data_points['thickness'].items():
                #     self.plot_data_point(self.plot_thickness, cell_type, HOURtoMCS_factor, count)
                if left_thicknesses:
                    self.plot_thickness.add_data_point("Limbus Avg Thickness", HOURtoMCS_factor, np.mean(left_thicknesses))
                    self.plot_thickness.add_data_point("Limbus Left Thickness", HOURtoMCS_factor, left_thicknesses[0])
                    self.plot_thickness.add_data_point("Limbus Right Thickness", HOURtoMCS_factor, left_thicknesses[-1])
                if right_thicknesses:
                    self.plot_thickness.add_data_point("Periphery Left Thickness", HOURtoMCS_factor, right_thicknesses[0])
                    self.plot_thickness.add_data_point("Periphery Right Thickness", HOURtoMCS_factor, right_thicknesses[-1])
                    self.plot_thickness.add_data_point("Periphery Avg Thickness", HOURtoMCS_factor, np.mean(right_thicknesses))
        
        if self.VolumeSurfaceDetailPlot and mcs % 10 == 0:

            SUPER_vol = []
            SUPER_surf = []
            WING_vol = []
            WING_surf = []
            BASAL_vol = []
            BASAL_surf = []
            STEM_vol = []
            STEM_surf = []

            for cell in self.cell_list_by_type(self.SUPER,self.WING,self.BASAL,self.STEM):
                if cell.type == self.SUPER:
                    SUPER_vol.append(cell.volume)
                    SUPER_surf.append(cell.surface)
                elif cell.type == self.WING:
                    WING_vol.append(cell.volume)
                    WING_surf.append(cell.surface)
                elif cell.type == self.BASAL:
                    BASAL_vol.append(cell.volume)
                    BASAL_surf.append(cell.surface)
                elif cell.type == self.STEM:
                    STEM_vol.append(cell.volume)
                    STEM_surf.append(cell.surface)
            # CSV
            self.data_points['volume_surface']['Superficial Volume'].append(np.mean(SUPER_vol))
            self.data_points['volume_surface']['Superficial Surface'].append(np.mean(SUPER_surf))
            self.data_points['volume_surface']['Wing Volume'].append(np.mean(WING_vol))
            self.data_points['volume_surface']['Wing Surface'].append(np.mean(WING_surf))
            self.data_points['volume_surface']['Basal Volume'].append(np.mean(BASAL_vol))
            self.data_points['volume_surface']['Basal Surface'].append(np.mean(BASAL_surf))
            self.data_points['volume_surface']['Stem Volume'].append(np.mean(STEM_vol))
            self.data_points['volume_surface']['Stem Surface'].append(np.mean(STEM_surf))
            self.data_points['volume_surface']['Time'].append(HOURtoMCS_factor)

            if self.CC3D_PLOT:
                # for cell_type, count in self.data_points['volume_surface'].items():
                #     self.plot_data_point(self.plot_thickness, cell_type, HOURtoMCS_factor, count)
                self.plot_volume_surface.add_data_point('Superficial Volume', HOURtoMCS_factor, np.mean(SUPER_vol))
                self.plot_volume_surface.add_data_point('Superficial Surface', HOURtoMCS_factor, np.mean(SUPER_surf))
                self.plot_volume_surface.add_data_point('Wing Volume', HOURtoMCS_factor, np.mean(WING_vol))
                self.plot_volume_surface.add_data_point('Wing Surface', HOURtoMCS_factor, np.mean(WING_surf))
                self.plot_volume_surface.add_data_point('Basal Volume', HOURtoMCS_factor, np.mean(BASAL_vol))
                self.plot_volume_surface.add_data_point('Basal Surface', HOURtoMCS_factor, np.mean(BASAL_surf))
                self.plot_volume_surface.add_data_point('Stem Volume', HOURtoMCS_factor, np.mean(STEM_vol))
                self.plot_volume_surface.add_data_point('Stem Surface', HOURtoMCS_factor, np.mean(STEM_surf))
            
        if self.MitosisPlot and mcs % 10 == 0:
                
                history_count = self.mitosis.get_divided_cells_counts()
                
                # CSV
                self.data_points['mitosis']['BASAL'].append(sum(history_count['Basal']))
                self.data_points['mitosis']['STEM'].append(sum(history_count['Stem']))
                self.data_points['mitosis']['Time'].append(HOURtoMCS_factor)
    
                if self.CC3D_PLOT:
                    # for cell_type, count in self.data_points['mitosis'].items():
                    #     self.plot_data_point(self.plot_mitosis, cell_type, HOURtoMCS_factor, count)
                    self.plot_mitosis.add_data_point("BASAL", HOURtoMCS_factor, sum(history_count['Basal']))
                    self.plot_mitosis.add_data_point("STEM", HOURtoMCS_factor, sum(history_count['Stem']))
                
                self.mitosis.clear_divided_cells_counts()

    def compute_thickness_for_all_columns(self, cell_type):
        thickness_dict = {}
        x_coords = [round(cell.xCOM) for cell in self.cell_list_by_type(self.WING)]
        clustered_x_coords = self.cluster_x_coords(x_coords)

        for x_coord in clustered_x_coords:
            min_y = float('inf')
            max_y = float('-inf')

            for cell in self.cell_list_by_type(cell_type):
                com_x = round(cell.xCOM)
                com_y = cell.yCOM

                if abs(com_x - x_coord) <= 4:
                    if com_y < min_y:
                        min_y = com_y
                    if com_y > max_y:
                        max_y = com_y

            thickness_dict[x_coord] = max_y - min_y  # Changed to calculate thickness

        return thickness_dict

    def cluster_x_coords(self, x_coords, threshold=1):
        if not x_coords:
            return []
        
        sorted_coords = sorted(x_coords)
        clusters = [[sorted_coords[0]]]

        for x in sorted_coords[1:]:
            # If the current x-coordinate is close to the last x-coordinate in the current cluster
            if abs(x - clusters[-1][-1]) <= threshold:
                clusters[-1].append(x)
            else:
                clusters.append([x])

        # Return the average x-coordinate for each cluster
        return [int(sum(cluster) / len(cluster)) for cluster in clusters]
    
    def plot_data_point(self, plot_obj, name, x, y):
        """
        Adds a data point to the CC3D plot if CC3D_PLOT is True.
        """
        plot_obj.add_data_point(name, x, y)

    def write_csv_for_category(self, category, mcs):
        """
        Writes a CSV file for the given category from the self.data_points dictionary.
        The file will be named after the category and saved in the current_script_directory.
        """        
        data = self.data_points[category]
        file_name = f"{category}_{mcs}.csv"
        file_path = self.current_script_directory.joinpath(file_name)
        
        # Assuming all sub-dictionaries have the same length as 'Time'
        num_entries = len(data['Time'])
        
        with open(file_path, 'w', newline='') as csvfile:
            csvwriter = csv.writer(csvfile)
            
            # Write the header (column names)
            headers = list(data.keys())
            csvwriter.writerow(headers)
            
            # Write the rows of data
            for i in range(num_entries):
                row = [data[header][i] for header in headers]
                csvwriter.writerow(row)