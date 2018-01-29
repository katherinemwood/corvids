__author__ = 'sean'
from Tkinter import *
import math, traceback
from tkFileDialog import askopenfilename, asksaveasfile

class FileLoader(Frame):
    def __init__(self, master, funcToPassTo):
        self.filename=""
        self.master = master
        self.entry = Entry(master=master, width=50)
        self.entry.delete(0,END)
        self.entry.insert(0,"Browse to your desired Path")
        self.functionToPassTo = funcToPassTo
        self.cbutton= Button(master, text="OK", command=self.finish)
        self.master.bind('<Return>', self.finish)
        self.cbutton.grid(row=1, column=16)
        self.bbutton= Button(self.master, text="Browse", command=self.browse_open_files)
        self.bbutton.grid(row=1, column=19)
        self.entry.grid(row=1, column = 8)

    def browse_open_files(self):

        Tk().withdraw()
        self.filename = askopenfilename(filetypes=(("RecreateData Files", "*.rd"),
                                           ("All files", "*.*")), initialdir="DATA/")
        self.entry.delete(0,END)
        self.entry.insert(0, self.filename)

    def finish(self,*event):
        if len(self.filename)>0:
            self.functionToPassTo(self.filename)
        self.master.destroy()

class FileSaver(Frame):
    def __init__(self, master, funcToPassTo):
        self.filename=""
        self.master = master
        self.entry = Entry(master=master, width=50)
        self.entry.delete(0,END)
        self.entry.insert(0,"Browse to your desired Path")
        self.functionToPassTo = funcToPassTo
        self.cbutton= Button(master, text="OK", command=self.finish)
        self.master.bind('<Return>', self.finish)
        self.cbutton.grid(row=1, column=16)
        self.bbutton= Button(self.master, text="Browse", command=self.browse_save_files)
        self.bbutton.grid(row=1, column=19)
        self.entry.grid(row=1, column = 8)

    def browse_save_files(self):

        Tk().withdraw()
        self.file = asksaveasfile(mode='w', defaultextension=".rd", initialdir="DATA/")
        if self.file:
            self.file.close()
            self.filename = self.file.name
            self.entry.delete(0,END)
            self.entry.insert(0, self.filename)

    def finish(self,*event):
        if len(self.filename)>0:
            self.functionToPassTo(self.filename)
        self.master.destroy()

def lambda_new_window(newWindowName, frameEnvirons, WindowTitle):
    '''
    Produce a new TkInter window with the title WindowTitle that exists in the namespace of frameEnvirons with the name
      frameEnvirons.newWindowName
    this is a function to be used inside a Lambda: The lambda has to exist on one line, so this makes it possible to
    neatly extend a lambda over multiple lines
    :param newWindowName:
    :param frameEnvirons:
    :param WindowTitle:
    :return:
    '''
    vars(frameEnvirons)[newWindowName] = Toplevel(frameEnvirons.master)
    vars(frameEnvirons)[newWindowName].minsize(width=400, height=50)
    vars(frameEnvirons)[newWindowName].title(WindowTitle)

def lambda_make_file_loader(frameEnvirons, frameName, fileLoaderName, filenameFunction):
    '''
    makes a FileLoader that handles browsing for and selecting a .RD file name
    another function to be used inside a lambda (see above)
    :param frameEnvirons: the namespace within which the file loader will be created
    :param frameName: the name of the file loader on the frameEnvirons object
    :param fileLoaderName: the internal name of the resulting file loader
    :param filenameFunction: the function to be run on the fileName set inside FileLoader()
    :return:
    '''
    vars(frameEnvirons)[fileLoaderName]=FileLoader(vars(frameEnvirons)[frameName], filenameFunction)

def lambda_make_file_saver(frameEnvirons, frameName, fileSaverName, filenameFunction):
    '''
    makes a FileSaver that handles browsing for and selecting a .RD file name
    another function to be used inside a lambda (see above)
    :param frameEnvirons: the namespace within which the file loader will be created
    :param frameName: the name of the file loader on the frameEnvirons object
    :param fileSaverName: the internal name of the resulting file loader
    :param filenameFunction: the function to be run on the fileName set inside FileSaver()
    :return:
    '''
    vars(frameEnvirons)[fileSaverName]=FileSaver(vars(frameEnvirons)[frameName], filenameFunction)


def fixed_callback(*args):
    print args

def bindings(widget, seq):
    return [x for x in widget.bind(seq).splitlines() if x.strip()]

def _funcid(binding):
    return binding.split()[1][3:]

def remove_binding(widget, seq, index=None, funcid=None):
    b = bindings(widget, seq)

    if not index is None:
        try:
            binding = b[index]
            widget.unbind(seq, _funcid(binding))
            b.remove(binding)
        except IndexError:
            print 'Binding #%d not defined.' % index
            return

    elif funcid:
        binding = None
        for x in b:
            if _funcid(x) == funcid:
                binding = x
                b.remove(binding)
                widget.unbind(seq, funcid)
                break
        if not binding:
            print 'Binding "%s" not defined.' % funcid
            return
    else:
        raise ValueError, 'No index or function id defined.'

    for x in b:
        widget.bind(seq, '+'+x, 1)

def scroll_lambda(frame):
    return lambda event, frame=frame: mouse_wheel(event,frame)

def mouse_wheel(event, frame):
    if(event.type == "MouseWheel"):
        # remove_binding(event.widget,"<Button-4>")
        # remove_binding(event.widget,"<Button-5>")
        if math.fabs(event.delta) == 120:
            frame.yview_scroll(-1*event.delta/120,"units")
        else:
            frame.yview_scroll(-1*event.delta,"units")
    elif(event.type == "4" or event.type == "5"):
        # remove_binding(event.widget,"<MouseWheel>")
        frame.yview_scroll(int(2*(int(event.num)-4.5)),"units")

def getTkinterLocation():
    """Returns the location of the Tkinter module."""
    import Tkinter
    if Tkinter.__file__.endswith('pyc'):
        return Tkinter.__file__[:-1]
    return Tkinter.__file__


def inTkinterMainloop():
    """Returns true if we're called in the context of Tkinter's
mainloop(), and false otherwise."""
    stack = traceback.extract_stack()
    tkinter_file = getTkinterLocation()
    for (file_name, lineno, function_name, text) in stack:
        if (file_name, function_name) == (tkinter_file, 'mainloop'):
            return 1
    return 0

class MenuItems(Frame):

    def __init__(self, master):
        Frame.__init__(self, master)
        self.menubar = None
        self.fileMenu = None
        self.parent = master
        self.initUI()
        self.postedMenu=None
        self.master.bind_all("<Button-1>", self.check_to_unpost)

    def check_to_unpost(self, event):
        if self.postedMenu and (not event.widget == self.postedMenu):
            self.postedMenu.unpost()

    def initUI(self):

        self.parent.title("CORVIDS Vers 1.0")

        self.menubar = Menu(self.parent)
        self.parent.config(menu=self.menubar)

        self.fileMenu = Menu(self.menubar, tearoff=False)
        self.menubar.add_cascade(label="File", menu=self.fileMenu)
        self.parent.update_idletasks()
        self.master.bind("<Alt-f>", lambda entry, menu=self.fileMenu: self.openMenu(menu))

    def addFunction(self, funcName, func):
        '''
        Adds a function 'func' to the MenuItems' namespace under the name 'funcName'.  Allows us to add new functions
        to the menu while keeping the MenuItem's code agnostic of application.
        :param funcName: The name of the new function as it will be called in the code
        :param func: The function to be added to the namespace of MenuItems
        :return:
        '''
        vars(self)[funcName]=func


    def openMenu(self, menu):
        if(self.postedMenu):
            self.postedMenu.unpost()
            if(self.postedMenu == menu):
                self.postedMenu = None
                return
        x =self.parent.winfo_x()
        y =self.parent.winfo_y()
        menu.post(x,y)
        self.postedMenu = menu

    def addMenu(self, menuName, text, hotKey=None):
        vars(self)[menuName]=Menu(self.menubar, tearoff=False)
        self.menubar.add_cascade(label=text, menu=vars(self)[menuName])
        if hotKey:
            self.master.bind(hotKey, lambda entry, menu=vars(self)[menuName]: self.openMenu(menu))

    def addTopLevelMenuItem(self, menu, text, function, *params, **kwargs):
        menu.add_command(label=text, command=function)
        if(kwargs.has_key("hotKey")):
            self.master.bind(kwargs["hotKey"],lambda event: function())

    def addQuit(self):
        self.fileMenu.add_command(label="Exit", command=self.onExit)
        self.master.bind("<Control-q>", lambda entry: self.onExit())

    def onExit(self):
        self.parent.quit()
        # Weird behavior with matplotlib.pyplot where after opening a window it seems to run a new mainloop
        #    and won't exit unless another root.quit() is called (multiple graphs will incriment the required quits).
        #    But no biggie since if there's only one mainloop it won't do anything to try to quit out of it again!
        self.parent.after(10, self.onExit)