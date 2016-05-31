import logging
import wx

import casie.gui.properties
import casie.gui.panels.base

class Panel(casie.gui.panels.base.Panel):
    
    ''' 
    wx Panel for session management functionality
    '''
    
    def __init__(self, parent):
        
        ''' 
        Set up the panel, add the controls, bind the events
        '''
        
        casie.gui.panels.base.Panel.__init__(self, parent)
        
    def add_controls(self):
        
        ''' 
        Add wx controls to the panel
        '''  
        
        self.wd_label = wx.StaticText(parent = self, label="Working Directory: ")
        self.wd_text = wx.StaticText(parent = self, label = self.GetTopLevelParent().work_dir)
        self.change_wd_button = wx.Button(parent = self, label="Change")
        
        self.active_tool = wx.StaticText(parent = self, label="Choose Tool: ")
        tools = [x[0] for x in casie.gui.properties._tools]
        tools.insert(0, "None")
        self.tool_choices = wx.ComboBox(parent = self, choices = tools)
        self.tool_choices.SetStringSelection(tools[0]) 

    def size_controls(self):
        
        ''' Size the controls
        '''  
        
        centre_fmt = [0, wx.CENTER, 3]
        expand_fmt = [1, wx.EXPAND, 3]
        buttonSizer = wx.FlexGridSizer(3, cols=1, vgap=3, hgap=3)
        buttonSizer.Add(self.wd_label, *centre_fmt)  
        buttonSizer.Add(self.wd_text, *centre_fmt) 
        buttonSizer.Add(self.change_wd_button, *centre_fmt)    
        buttonSizer.Add(self.active_tool, *centre_fmt)  
        buttonSizer.Add(self.tool_choices, *centre_fmt)  
        
        hbox = wx.BoxSizer()
        hbox.AddSpacer(10)
        hbox.Add(buttonSizer, *expand_fmt)
        hbox.AddSpacer(10)
        
        vbox = wx.BoxSizer(wx.VERTICAL)
        vbox.AddSpacer(10)
        vbox.Add(hbox, *expand_fmt)
        vbox.AddSpacer(10)

        self.SetSizer(vbox)
        vbox.Fit(self) 
        
    def bind_events(self):
        
        ''' 
        Bind the events
        '''  
        
        self.change_wd_button.Bind(wx.EVT_BUTTON, self.on_change_wd)
        self.tool_choices.Bind(wx.EVT_COMBOBOX, self.on_change_tool)
        
    def on_change_wd(self, event):
        
        ''' 
        Open a dialogue for the user to select a new working directory.
        ''' 
             
        dir_dialogue = wx.DirDialog(self, style=wx.FD_OPEN)
        dir_dialogue.ShowModal()
        
        self.GetTopLevelParent().work_dir = dir_dialogue.GetPath()
        dir_dialogue.Destroy()
        
        self.wd_text.Label = self.GetTopLevelParent().work_dir
        self.refresh_controls()
        
        logging.info("Changed working directory to: " + str(self.GetTopLevelParent().work_dir))
        
    def on_change_tool(self, event):
        
        ''' 
        Change the active tool
        '''       
        
        new_name = self.tool_choices.GetValue()
        self.tool_choices.SetStringSelection(new_name) 
        self.GetTopLevelParent().change_tool(new_name)
        self.refresh_controls()
        
        logging.info("Changed active tool to: " + new_name)
            
    def refresh_controls(self):
        
        ''' 
        Refresh the text values of the panel controls
        '''  
        
        pass

    def reset(self):
        
        ''' 
        Refresh the panel controls
        '''  
        
        self.refresh_controls()