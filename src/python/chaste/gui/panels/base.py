import wx

class Panel(wx.Panel):
    
    ''' 
    A wx panel with some extra useful methods
    '''
    
    def __init__(self, parent, border = True):
        
        ''' 
        Set up the panel, add the controls, bind the events
        '''
        
        if border:
            wx.Panel.__init__(self, parent, style=wx.SIMPLE_BORDER)
        else:
            wx.Panel.__init__(self, parent)
            
        self.SetBackgroundColour((255, 255, 255, 1.0))
        
        self.add_controls()
        self.size_controls()
        self.bind_events()
        self.name = "None"
        
    def add_controls(self):
        
        ''' 
        Add wx controls to the panel
        '''  
        
        pass

    def size_controls(self):
        
        ''' 
        Size the controls
        '''  
        
        pass 
        
    def bind_events(self):
        
        ''' 
        Bind the events
        '''  
        
        pass

    def reset(self):
        
        ''' Default panel reset functionality
        '''  
        
        pass
    
    def get_file_name(self):
        
        ''' 
        Return a user specified file path
        ''' 
        
        file_dialogue = wx.FileDialog(self, style=wx.FD_OPEN)
        file_dialogue.ShowModal()
        file_name = file_dialogue.GetPath()
        file_dialogue.Destroy()
        return file_name