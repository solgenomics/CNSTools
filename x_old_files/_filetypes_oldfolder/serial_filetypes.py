from abc import ABCMeta, abstractmethod
class Serial_Filetype(object):
    """Abstract parent class of all of the filetype functions. The __init__ function here should be called by each subclass."""
    __metaclass__ = ABCMeta

    Entry_class = None

    def __init__(self,file_name=None,file_object=None,lines=None):
        """Accepts a file object, file name, or list of line strings. Sets up an entry list (self.entries) and then calls the abstract add_lines() function with a list of lines to be parsed into entrys."""
        self.entries = []
        if lines:
            self.add_lines(lines)
        elif file_object==None and file_name==None:
            return None
        else:
            with (file_object if file_object!=None else open(file_name)) as file:
                self.add_lines(file.readlines())

    def add_entry(self,*args,**kwargs): 
        """Creates a new entry by passing along the arguements to instantiate the Entry_class then appends that new instance to self.entries. Entry_class must be defined in the subclass for the add_entry Serial_Filetype function to work, otherwise it will return None.""" 
        if self.Entry_class:
            new_entry = self.Entry_class(*args,**kwargs) 
            self.entries.append(new_entry)
            return new_entry
        else: 
            return None

    def save_file(self,save_name):
        with open(save_name,"w") as out:
            out.write("\n".join(self.get_lines()))

    @abstractmethod
    def add_lines(self,file): 
        """Abstract. Should take a list of line strings and parse them into instances of the Entry_class then append them to self.entries"""
        pass
    @abstractmethod
    def get_lines(self): 
        """Abstract. Should return the entry data formatted as a list of line strings"""
        pass