import logging
import wx
import wx.lib.newevent
import casie.gui.panels.base

# create event type
wxLogEvent, EVT_WX_LOG_EVENT = wx.lib.newevent.NewEvent()

class wxLogHandler(logging.Handler):
    
    """
    A handler class which sends log strings to a wx object
    """
    
    def __init__(self, wxDest=None):
        
        """
        Initialize the handler
        @param wxDest: the destination object to post the event to 
        @type wxDest: wx.Window
        """
        logging.Handler.__init__(self)
        self.wxDest = wxDest
        self.level = logging.INFO

    def flush(self):
        
        """
        does nothing for this handler
        """

    def emit(self, record):
        
        """
        Emit a record.

        """
        try:
            msg = self.format(record)
            evt = wxLogEvent(message=msg,levelname=record.levelname)            
            wx.PostEvent(self.wxDest,evt)
        except (KeyboardInterrupt, SystemExit):
            raise
        except:
            self.handleError(record)

class Panel(casie.gui.panels.base.Panel):
    
    ''' A wx panel with some useful methods for getting TopLevelParent attributes.
    '''
    
    def __init__(self, parent):
        
        ''' Set up the panel, add the controls, bind the events
        '''
        
        casie.gui.panels.base.Panel.__init__(self, parent)
        self.setup_event_handler()
        
    def setup_event_handler(self):
        logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(levelname)-8s %(message)s', datefmt='%d-%m %H:%M')
        handler = wxLogHandler(self)
        handler.setLevel(logging.INFO)
        formatter = logging.Formatter('%(asctime)s %(levelname)-8s %(message)s', datefmt='%d-%m %H:%M')
        handler.setFormatter(formatter)
        logging.getLogger('').addHandler(handler)
        
    def add_controls(self):
        
        ''' Add wx controls to the panel
        '''  
        
        self.log_window = wx.TextCtrl(self, wx.ID_ANY, size=(1000,150), style = wx.TE_MULTILINE|wx.TE_READONLY|wx.HSCROLL)

    def size_controls(self):
        
        ''' Size the controls
        '''  
        
        sizer = wx.BoxSizer()
        sizer.Add(self.log_window, 1, wx.EXPAND, 5)
        sizer.Fit(self)
        
    def bind_events(self):
        
        ''' Bind the events
        '''  
        
        self.Bind(EVT_WX_LOG_EVENT, self.onLogEvent)
        
    def onLogEvent(self,event):
        
        '''
        Add event.message to text window
        '''
        
        msg = event.message.strip("\r")+"\n"
        self.log_window.AppendText(msg)
        event.Skip()