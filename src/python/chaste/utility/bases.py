class SimpleIOBase(object):
    
    """
    Class for ensuring consistency in the definition of simple tools with a single input and single output.
    Classes defining operations with low numbers of inputs can inherit from this one and overwrite the 
    update method.
    
    @param self.input the input to the tool
    @return self.output the tool output 
    """
    
    def __init__(self):
        self.input = None
        self.output = None
        
    def set_input(self, my_input):
        self.input = my_input
        
    def update(self):
        pass
    
    def get_output(self):
        return self.output()
    