import casie.gui.properties
import casie.controller.base

import casie.simulation.simple.cell


class Controller(casie.controller.base.Controller):
    
    def __init__(self, controller_id = 0, label = None, io_data = None):
        
        casie.controller.base.Controller.__init__(self, controller_id, "generator", "cell", 
                                                  label, io_data)

        self.geometry = None

    def activate(self):
        population = casie.simulation.simple.cell.CellPopulation()
        
        for eachConnector in self.input_data_connectors:
            controller_type, data = eachConnector.activate()
            if controller_type[0] == "generator" and controller_type[1] == "part":
                cells = self.generate_cells(data)
                population.add_cells(cells)
        return self.own_type, population 
    
    def generate_cells(self, geometry):
        points = geometry[0]
        cells = []
        
        for eachPoint in points:
            cells.append(casie.simulation.simple.cell.Cell(eachPoint))
            
        return cells
    
if casie.gui.properties._is_gui_session:
    
    class Panel(casie.controller.base.Panel):
    
        ''' Panel for holding the configuration options and canvas panels
        '''
        
        def __init__(self, parent):
            casie.controller.base.Panel.__init__(self, parent, "generator", "cell")
    
        def add_gui_controls(self):
            pass