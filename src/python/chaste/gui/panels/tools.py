import wx

import chaste.gui.properties
import chaste.gui.panels.base
import chaste.utility.recursion

class Panel(chaste.gui.panels.base.Panel):
    
    ''' Panel for holding controller specific configuration option panels
    '''
    
    def __init__(self, parent):
        
        ''' 
        Set up the panel. Do all control and sizing here.
        '''
        
        chaste.gui.panels.base.Panel.__init__(self, parent)
        
        self.panels = []
        self.panels.append(chaste.gui.panels.base.Panel(self))
        self.active_panel = self.panels[0]
        
        self.hbox = wx.BoxSizer()
        self.hbox.Add(self.panels[0], 1, wx.EXPAND)
        self.add_hidden_panels()
        self.SetSizer(self.hbox)
        self.hbox.Fit(self) 
        
    def add_hidden_panels(self):
        
        ''' 
        Add a hidden panel for each tool in the registry
        '''
        
        for eachTool in chaste.gui.properties._tools:
            class_name = eachTool[1] + "Panel"
            new_panel_class = chaste.utility.recursion.get_class(class_name) 
            new_panel = new_panel_class(self)
            new_panel.Hide()
            self.panels.append(new_panel)
            self.hbox.Add(new_panel, 1, wx.EXPAND)
            
    def get_panel(self, tool_name):
        
        ''' 
        Return the panel corresponding to the tool name
        '''
        
        query_panel = next((panel for panel in self.panels if (panel.name == tool_name)), None)
        return query_panel
    
    def get_active_panel(self):
        
        ''' 
        Return the current active panel
        '''     
        
        return self.active_panel
    
    def show_panel(self, panel):
        
        ''' 
        Show the specified panel
        '''  
        
        for eachPanel in self.panels:
            if eachPanel == panel:
                self.hide_panels()
                eachPanel.Show()
                self.active_panel = eachPanel
                self.Layout()
                break
        
    def hide_panels(self):
  
        ''' Hide all the panels
        '''        
        
        map(lambda x: x.Hide(), self.panels) 
            
    def reset(self):
        
        ''' Reset this panel
        '''
               
        self.hide_panels()
        self.panels[0].Show()