import wx

class MainMenuBar(wx.MenuBar):
    
    ''' 
    Main Menu Bar
    '''  
    
    def __init__(self):
        wx.MenuBar.__init__(self)
        
        self.file_menu = FileMenu()
        self.Append(self.file_menu, self.file_menu.label)

class FileMenu(wx.Menu):
    
    ''' 
    File Menu
    '''    
    
    def __init__(self):
        
        wx.Menu.__init__(self)
        self.label = "&File"
        self.exit = self.Append(wx.ID_EXIT,"E&xit"," Terminate the program")